"""check -> changeset -> dry-run -> apply. Local data changes only in apply()."""
import hashlib
import json
import os
import shutil
from datetime import datetime

from . import config, sources
from .db import audit, connect, tx

GENOME_COLS = [
    "accession_version", "genbank_accessions", "refseq_accessions", "taxon_id", "genome_name",
    "taxon_lineage_ids", "taxon_lineage_names", "genome_status", "genome_length", "segment",
    "strain", "subtype", "genotype", "host_name", "isolation_country", "geographic_location",
    "collection_date", "isolation_source", "publication", "submitter", "genome_quality",
    "sequence_hash", "bvbrc_genome_id", "pathoplexus_accession", "release_date", "update_date",
]
COMPARE = ["accession_version"] + config.TRACKED_FIELDS


# ------------------------------------------------------------------ watchlist
def add_watch(name: str, taxon: str | None = None, filters: dict | None = None) -> int:
    taxon = taxon or name
    tid = sources.ncbi_resolve_taxon_id(taxon)
    sp = sources.ncbi_species_scope(tid)
    if sp:  # watch the species so records tagged at the species node are not missed
        tid, taxon = sp
    with tx() as con:
        cur = con.execute(
            "INSERT INTO watchlist(name, taxon, taxon_id, filters_json) VALUES(?,?,?,?)",
            (name, taxon, tid, json.dumps(filters or {})))
        audit(con, "watch_add", f"{name} taxon_id={tid}")
        return cur.lastrowid


# ---------------------------------------------------------------------- merge
def merge(per_source: dict[str, dict[str, dict]]) -> dict[str, tuple[dict, dict]]:
    """NCBI defines the accession set; other sources enrich. -> {acc: (record, provenance)}"""
    merged = {}
    for acc in per_source["ncbi"]:
        rec, prov = {}, {}
        for src in config.SOURCE_PRIORITY:
            for k, v in per_source.get(src, {}).get(acc, {}).items():
                if isinstance(v, (list, tuple, set)):  # sources occasionally return arrays (e.g. segment)
                    v = "; ".join(str(x) for x in v)
                if k not in rec and v not in (None, ""):
                    rec[k], prov[k] = v, src
        merged[acc] = (rec, prov)
    return merged


def _fingerprint(kind: str, acc: str, rec: dict) -> str:
    keep = {k: rec.get(k) for k in COMPARE}
    return hashlib.sha1(json.dumps([kind, acc, keep], sort_keys=True, default=str).encode()).hexdigest()


def _diff(old: dict | None, new: dict) -> dict:
    d = {}
    for k in COMPARE:
        o = old.get(k) if old else None
        n = new.get(k)
        if str(o or "") != str(n or ""):
            d[k] = [o, n]
    return d


# ---------------------------------------------------------------------- check
def check(watch_name: str, progress=None) -> int | None:
    """Fetch sources, diff against local DB, store a changeset (no local data change).
    Returns changeset id or None when nothing changed."""
    p = progress or (lambda *a: None)
    p(2, "워치리스트 확인 중…")
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE name=?", (watch_name,)).fetchone()
        if not w:
            raise SystemExit(f"unknown watch: {watch_name}")
        notes0 = []
        sp = sources.ncbi_species_scope(w["taxon_id"])
        if sp:
            con.execute("UPDATE watchlist SET taxon=?, taxon_id=? WHERE id=?", (sp[1], sp[0], w["id"]))
            con.commit()
            audit(con, "watch_scope", f"{watch_name}: taxon widened {w['taxon']} ({w['taxon_id']}) -> {sp[1]} ({sp[0]})")
            notes0.append(f"taxon widened to species {sp[1]} ({sp[0]})")
            w = con.execute("SELECT * FROM watchlist WHERE id=?", (w["id"],)).fetchone()
        per_source, notes = sources.fetch_all(w["name"], w["taxon"], w["taxon_id"],
                                              json.loads(w["filters_json"] or "{}"), progress=p)
        notes = notes0 + notes
        p(75, "로컬 데이터와 비교 중…")
        if not per_source["ncbi"]:
            raise SystemExit("NCBI returned no records; refusing to compute removals")
        merged = merge(per_source)
        local = {r["accession"]: dict(r) for r in
                 con.execute("SELECT * FROM genome WHERE watch_id=?", (w["id"],))}
        rejected = {r[0] for r in con.execute(
            "SELECT i.fingerprint FROM changeset_item i JOIN changeset c ON c.id=i.changeset_id "
            "WHERE c.watch_id=? AND i.decision='rejected'", (w["id"],))}
        held = {r[0] for r in con.execute(
            "SELECT i.fingerprint FROM changeset_item i JOIN changeset c ON c.id=i.changeset_id "
            "WHERE c.watch_id=? AND i.decision='hold'", (w["id"],))}

        items = []
        for acc, (rec, prov) in merged.items():
            old = local.get(acc)
            if old is None:
                kind, diff = "new", {k: [None, rec.get(k)] for k in COMPARE if rec.get(k) is not None}
            else:
                diff = _diff(old, rec)
                kind = "changed" if diff else None
            if kind:
                items.append((acc, kind, diff, rec, prov))
        for acc, old in local.items():
            if acc not in merged and old.get("genome_status") != "Removed":
                items.append((acc, "removed", {"genome_status": [old.get("genome_status"), "Removed"]}, old, {}))

        fresh = []
        for acc, kind, diff, rec, prov in items:
            fp = _fingerprint(kind, acc, rec)
            if fp in rejected:
                continue
            fresh.append((acc, kind, diff, rec, prov, fp))

        p(88, "changeset 저장 중…")
        now = datetime.now().isoformat(timespec="seconds")
        for src, recs in per_source.items():
            for acc, rec in recs.items():
                con.execute("INSERT OR REPLACE INTO source_record(source,accession,watch_id,payload_json,fetched_at)"
                            " VALUES(?,?,?,?,?)", (src, acc, w["id"], json.dumps(rec), now))
        con.execute("UPDATE watchlist SET last_checked=? WHERE id=?", (now, w["id"]))
        if not fresh:
            audit(con, "check", f"{watch_name}: no changes; {'; '.join(notes)}")
            con.commit()
            return None

        con.execute("UPDATE changeset SET status='DISCARDED' WHERE watch_id=? AND status IN ('PENDING','READY','BLOCKED')",
                    (w["id"],))
        summary = {"new": 0, "changed": 0, "removed": 0, "notes": notes,
                   "sources": {s: len(r) for s, r in per_source.items()}}
        cur = con.execute("INSERT INTO changeset(watch_id, summary_json) VALUES(?,?)", (w["id"], "{}"))
        cs = cur.lastrowid
        for acc, kind, diff, rec, prov, fp in fresh:
            summary[kind] += 1
            con.execute(
                "INSERT INTO changeset_item(changeset_id,accession,kind,diff_json,new_record_json,provenance_json,fingerprint,decision)"
                " VALUES(?,?,?,?,?,?,?,?)",
                (cs, acc, kind, json.dumps(diff), json.dumps(rec, default=str), json.dumps(prov), fp,
                 "hold" if fp in held else "pending"))
        con.execute("UPDATE changeset SET summary_json=? WHERE id=?", (json.dumps(summary), cs))
        audit(con, "check", f"{watch_name}: changeset {cs} {summary}")
        con.commit()
        return cs
    finally:
        con.close()


# -------------------------------------------------------------------- dry-run
def _versioned(acc: str, rec: dict) -> str:
    v = rec.get("accession_version")
    return f"{acc}.{v}" if v else acc


def dryrun(changeset_id: int, sample: int | None = None, progress=None) -> dict:
    p = progress or (lambda *a: None)
    p(92, "Dry-run: 링크·디스크 검사 중…")
    sample = sample if sample is not None else config.DRYRUN_SAMPLE_SIZE
    con = connect()
    try:
        items = con.execute("SELECT * FROM changeset_item WHERE changeset_id=? AND kind IN ('new','changed')",
                            (changeset_id,)).fetchall()
        checks = []
        recs = {i["accession"]: json.loads(i["new_record_json"]) for i in items}
        total_bytes = sum((r.get("genome_length") or 0) for r in recs.values())

        free = shutil.disk_usage(config.DATA_DIR if config.DATA_DIR.exists() else ".").free
        need = total_bytes * 2 + config.MIN_FREE_BYTES
        checks.append({"name": "disk_space", "ok": free >= need,
                       "detail": f"free {free // 2**20} MiB, need ~{need // 2**20} MiB for {len(recs)} genomes ({total_bytes} bp)"})

        ok_dl = True
        picked = list(recs)[:sample]
        if picked:
            try:
                got = sources.download_ncbi_fasta([_versioned(a, recs[a]) for a in picked])
                for a in picked:
                    va = _versioned(a, recs[a])
                    seq = got.get(va)
                    exp = recs[a].get("genome_length")
                    good = bool(seq) and (exp is None or len(seq) == exp)
                    ok_dl &= good
                    checks.append({"name": f"sample {va}", "ok": good,
                                   "detail": (f"{len(seq)} bp, md5 {hashlib.md5(seq.encode()).hexdigest()}"
                                              + ("" if good else f" (expected {exp})")) if seq else "not returned by NCBI"})
            except Exception as e:
                ok_dl = False
                checks.append({"name": "download_endpoint", "ok": False, "detail": str(e)})
        elif not recs:
            checks.append({"name": "download_endpoint", "ok": True, "detail": "nothing to download"})

        ready = all(c["ok"] for c in checks)
        report = {"checked_at": datetime.now().isoformat(timespec="seconds"), "checks": checks}
        con.execute("UPDATE changeset SET status=?, dryrun_json=? WHERE id=?",
                    ("READY" if ready else "BLOCKED", json.dumps(report), changeset_id))
        audit(con, "dryrun", f"changeset {changeset_id}: {'READY' if ready else 'BLOCKED'}")
        con.commit()
        return report
    finally:
        con.close()


# ---------------------------------------------------------------------- apply
def decide(item_ids: list[int], decision: str):
    assert decision in ("approved", "rejected", "hold", "pending")
    with tx() as con:
        for i in item_ids:
            con.execute("UPDATE changeset_item SET decision=? WHERE id=? AND decision!='applied'", (decision, i))
        audit(con, "decide", f"{decision}: {item_ids}")


def _gc(seq: str) -> float:
    return round(100 * sum(seq.count(b) for b in "GCgc") / len(seq), 3) if seq else 0.0


def _read_staged(f, versioned: str):
    """Sequence from a staged FASTA written by apply(), or None if missing/mismatched."""
    try:
        lines = f.read_text().splitlines()
    except OSError:
        return None
    if not lines or not lines[0].startswith(">" + versioned + " ") and lines[0][1:] != versioned:
        return None
    return "".join(lines[1:])


def apply(changeset_id: int, progress=None) -> dict:
    """Apply approved items atomically: files first (with backups), then one DB transaction."""
    p = progress or (lambda *a: None)
    p(3, "적용 준비 중…")
    con = connect()
    backups: list[tuple[str, str | None]] = []
    try:
        cs = con.execute("SELECT * FROM changeset WHERE id=?", (changeset_id,)).fetchone()
        if cs["status"] != "READY":
            raise SystemExit(f"changeset {changeset_id} is {cs['status']}; run dry-run until READY")
        items = con.execute("SELECT * FROM changeset_item WHERE changeset_id=? AND decision='approved'",
                            (changeset_id,)).fetchall()
        if not items:
            return {"applied": 0}
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (cs["watch_id"],)).fetchone()
        stage = config.DATA_DIR / "staging" / str(changeset_id)
        snap = config.DATA_DIR / "snapshots" / f"cs{changeset_id}-{datetime.now():%Y%m%d-%H%M%S}"
        stage.mkdir(parents=True, exist_ok=True)

        recs = {i["accession"]: json.loads(i["new_record_json"]) for i in items if i["kind"] != "removed"}
        meta: dict[str, dict] = {}  # staged files only; sequences are never all held in memory
        failed: dict[str, str] = {}
        todo = []
        for a, r in recs.items():  # resume: reuse files staged by an interrupted run if they still verify
            v = _versioned(a, r)
            f = stage / f"{a}.fna"
            seq = _read_staged(f, v)
            if seq is not None and (not r.get("genome_length") or len(seq) == r["genome_length"]):
                meta[a] = {"md5": hashlib.md5(seq.encode()).hexdigest(), "len": len(seq), "gc": _gc(seq), "file": f}
            else:
                todo.append((a, v))
        for start in range(0, len(todo), 100):
            batch = todo[start:start + 100]
            p(5 + 65 * start // max(len(todo), 1), f"서열 다운로드 중… {start}/{len(todo)}")
            try:
                got = sources.download_ncbi_fasta([v for _, v in batch])
            except Exception as e:  # one bad batch must not lose the whole run
                for a, v in batch:
                    failed[a] = f"download error: {e}"
                continue
            for a, v in batch:
                s = got.get(v)
                exp = recs[a].get("genome_length")
                if not s:
                    failed[a] = "not returned by NCBI"
                    continue
                if exp and len(s) != exp:
                    failed[a] = f"length {len(s)} != expected {exp}"
                    continue
                f = stage / f"{a}.fna"
                f.write_text(f">{v} {recs[a].get('genome_name', '')}\n" +
                             "\n".join(s[i:i + 70] for i in range(0, len(s), 70)) + "\n")
                meta[a] = {"md5": hashlib.md5(s.encode()).hexdigest(), "len": len(s), "gc": _gc(s), "file": f}
        if failed and not meta and recs:
            raise SystemExit(f"all downloads failed, e.g. {next(iter(failed.items()))}")
        p(72, f"서열 검증 완료 ({len(failed)}건 실패), 파일·DB 기록 중…")
        taxon_dir = config.DATA_DIR / "genomes" / (w["name"].replace(" ", "_"))
        taxon_dir.mkdir(parents=True, exist_ok=True)
        snap.mkdir(parents=True, exist_ok=True)
        now = datetime.now().isoformat(timespec="seconds")
        con.execute("BEGIN")
        next_n: dict = {}
        done_ids = []
        for it in items:
            acc, rec = it["accession"], json.loads(it["new_record_json"])
            if acc in failed:
                continue
            rec = {k: ("; ".join(map(str, v)) if isinstance(v, list) else v) for k, v in rec.items()}
            old = con.execute("SELECT * FROM genome WHERE accession=?", (acc,)).fetchone()
            if old:
                con.execute("INSERT INTO genome_version(genome_id,snapshot_json,changeset_id) VALUES(?,?,?)",
                            (old["genome_id"], json.dumps(dict(old)), changeset_id))
            if it["kind"] == "removed":
                con.execute("UPDATE genome SET genome_status='Removed', date_modified=? WHERE accession=?", (now, acc))
                done_ids.append(it["id"])
                continue
            if old:
                gid = old["genome_id"]
            else:
                tid = rec.get("taxon_id") or w["taxon_id"] or 0
                if tid not in next_n:
                    nums = [int(r[0].split(".", 1)[1]) for r in con.execute(
                        "SELECT genome_id FROM genome WHERE taxon_id=? OR genome_id LIKE ?", (tid, f"{tid}.%"))
                        if r[0].split(".", 1)[1].isdigit()]
                    next_n[tid] = max(nums, default=0)
                next_n[tid] += 1
                gid = f"{tid}.{next_n[tid]}"
            m = meta[acc]
            md5 = rec["sequence_md5"] = m["md5"]
            path = taxon_dir / f"{gid}.fna"
            tmp = m["file"]
            if path.exists():
                bak = snap / path.name
                shutil.copy2(path, bak)
                backups.append((str(path), str(bak)))
            else:
                backups.append((str(path), None))
            os.replace(tmp, path)
            cols = [c for c in GENOME_COLS if c in rec] + ["sequence_md5"]
            vals = [rec[c] for c in cols]
            if old:
                con.execute(f"UPDATE genome SET {','.join(c + '=?' for c in cols)}, date_modified=? WHERE accession=?",
                            vals + [now, acc])
            else:
                con.execute(f"INSERT INTO genome(genome_id,watch_id,accession,{','.join(cols)},date_inserted,date_modified)"
                            f" VALUES(?,?,?,{','.join('?' * len(cols))},?,?)",
                            [gid, w["id"], acc] + vals + [now, now])
            con.execute("INSERT OR REPLACE INTO genome_sequence(accession,genome_id,version,length,gc_content,md5,file_path)"
                        " VALUES(?,?,?,?,?,?,?)", (acc, gid, rec.get("accession_version"), m["len"], m["gc"], md5, str(path)))
            for f, s in json.loads(it["provenance_json"] or "{}").items():
                con.execute("INSERT OR REPLACE INTO field_provenance(genome_id,field,source) VALUES(?,?,?)", (gid, f, s))
            done_ids.append(it["id"])
        for i in range(0, len(done_ids), 500):
            chunk = done_ids[i:i + 500]
            con.execute(f"UPDATE changeset_item SET decision='applied' WHERE id IN ({','.join('?' * len(chunk))})", chunk)
        left = con.execute("SELECT COUNT(*) FROM changeset_item WHERE changeset_id=? AND decision IN ('pending','hold','approved')",
                           (changeset_id,)).fetchone()[0]
        if not left:
            con.execute("UPDATE changeset SET status='APPLIED', applied_at=? WHERE id=?", (now, changeset_id))
        audit(con, "apply", f"changeset {changeset_id}: {len(done_ids)} applied, {len(failed)} failed"
              + (f" e.g. {list(failed.items())[:3]}" if failed else ""))
        con.commit()
        if not failed:
            shutil.rmtree(stage, ignore_errors=True)
        p(100, "완료")
        return {"applied": len(done_ids), "failed": failed, "remaining": left}
    except BaseException:
        con.rollback()
        for path, bak in backups:  # restore files
            if bak:
                shutil.copy2(bak, path)
            elif os.path.exists(path):
                os.remove(path)
        raise
    finally:
        con.close()
