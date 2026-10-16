"""check -> changeset -> dry-run -> apply. Local data changes only in apply()."""
import hashlib
import json
import os
import shutil
from pathlib import Path
from datetime import datetime

from . import config, sources
from .db import audit, connect, tx

GENOME_COLS = [
    "accession_version", "genbank_accessions", "refseq_accessions", "taxon_id", "genome_name",
    "taxon_lineage_ids", "taxon_lineage_names", "genome_status", "genome_length", "segment",
    "strain", "subtype", "genotype", "host_name", "isolation_country", "geographic_location",
    "collection_date", "isolation_source", "publication", "submitter", "genome_quality",
    "sequence_hash", "bvbrc_genome_id", "pathoplexus_accession", "release_date", "update_date",
    "sra_accession", "bioproject_accession", "biosample_accession", "assembly_method", "sequencing_platform",
    "clade", "lineage", "serotype", "outbreak",
]
COMPARE = ["accession_version"] + config.TRACKED_FIELDS


# ------------------------------------------------------------------ watchlist
def add_watch(name: str, taxon: str | None = None, filters: dict | None = None) -> int:
    tid, taxon = sources.ncbi_resolve_virus(taxon or name)  # ValueError if not a virus / no records
    sp = sources.ncbi_species_scope(tid)
    if sp:  # watch the species so records tagged at the species node are not missed
        tid, taxon = sp
    with tx() as con:
        cur = con.execute(
            "INSERT INTO watchlist(name, taxon, taxon_id, filters_json) VALUES(?,?,?,?)",
            (name, taxon, tid, json.dumps(filters or {})))
        audit(con, "watch_add", f"{name} taxon_id={tid}")
        return cur.lastrowid


def watch_impact(watch_id: int) -> dict:
    """What deleting a watch would remove (shown to the user before they confirm)."""
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("워치를 찾을 수 없습니다.")
        one = lambda sql: con.execute(sql, (watch_id,)).fetchone()[0]
        return {"name": w["name"], "genomes": one("SELECT COUNT(*) FROM genome WHERE watch_id=?"),
                "fasta": one("SELECT COUNT(*) FROM genome_sequence WHERE genome_id IN (SELECT genome_id FROM genome WHERE watch_id=?)"),
                "changesets": one("SELECT COUNT(*) FROM changeset WHERE watch_id=?"),
                "items": one("SELECT COUNT(*) FROM changeset_item WHERE changeset_id IN (SELECT id FROM changeset WHERE watch_id=?)")}
    finally:
        con.close()


def update_watch(watch_id: int, name: str | None = None, taxon: str | None = None, enabled: bool | None = None) -> list[str]:
    """Rename, re-point to another taxon (validated like add_watch) and/or (de)activate automatic checks.
    Returns human-readable notes; raises ValueError for invalid input."""
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("워치를 찾을 수 없습니다.")
        sets: dict = {}
        notes: list[str] = []
        if name is not None and name.strip() and name.strip() != w["name"]:
            new = name.strip()
            if con.execute("SELECT 1 FROM watchlist WHERE name=? AND id!=?", (new, watch_id)).fetchone():
                raise ValueError(f"'{new}' 이름의 워치가 이미 있습니다.")
            sets["name"] = new
            notes.append(f"이름 {w['name']} → {new}")
        if taxon is not None and taxon.strip() and taxon.strip().lower() not in {str(w["taxon"]).lower(), str(w["taxon_id"])}:
            tid, tname = sources.ncbi_resolve_virus(taxon.strip())  # ValueError if not a virus / no records
            sp = sources.ncbi_species_scope(tid)
            if sp:
                tid, tname = sp
            if tid != w["taxon_id"]:
                sets.update(taxon=tname, taxon_id=tid, pathoplexus_organism=None)
                notes.append(f"수집 분류군 {w['taxon']} ({w['taxon_id']}) → {tname} ({tid})")
                n = con.execute("SELECT COUNT(*) FROM genome WHERE watch_id=?", (watch_id,)).fetchone()[0]
                if n:
                    notes.append(f"이미 반영된 유전체 {n}건 중 새 분류군에 없는 것은 다음 확인에서 '삭제됨' 후보로 올라올 수 있습니다")
        if enabled is not None and int(bool(enabled)) != int(w["enabled"]):
            sets["enabled"] = int(bool(enabled))
            notes.append("자동 확인 " + ("켬" if enabled else "끔"))
        if not sets:
            return []
        con.execute(f"UPDATE watchlist SET {','.join(k + '=?' for k in sets)} WHERE id=?", [*sets.values(), watch_id])
        audit(con, "watch_update", f"{w['name']}: " + "; ".join(notes))
        con.commit()
        return notes
    finally:
        con.close()


def delete_watch(watch_id: int) -> dict:
    """Delete a watch together with everything it collected (genomes, FASTA files, changesets, history).
    Irreversible: the web UI asks the user to type the watch name first."""
    con = connect()
    files: list[str] = []
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("워치를 찾을 수 없습니다.")
        impact = watch_impact(watch_id)
        files = [r[0] for r in con.execute(
            "SELECT file_path FROM genome_sequence WHERE file_path IS NOT NULL AND genome_id IN (SELECT genome_id FROM genome WHERE watch_id=?)",
            (watch_id,))]
        con.execute("BEGIN")
        mine = "(SELECT genome_id FROM genome WHERE watch_id=?)"
        for table in ("field_provenance", "genome_version", "genome_sequence"):
            con.execute(f"DELETE FROM {table} WHERE genome_id IN {mine}", (watch_id,))
        con.execute("DELETE FROM changeset_item WHERE changeset_id IN (SELECT id FROM changeset WHERE watch_id=?)", (watch_id,))
        for table in ("genome", "changeset", "source_record"):
            con.execute(f"DELETE FROM {table} WHERE watch_id=?", (watch_id,))
        con.execute("DELETE FROM watchlist WHERE id=?", (watch_id,))
        audit(con, "watch_delete", f"{w['name']}: {impact['genomes']} genomes, {impact['fasta']} fasta, {impact['changesets']} changesets")
        con.commit()
    except BaseException:
        con.rollback()
        raise
    finally:
        con.close()
    root = (config.DATA_DIR / "genomes").resolve()
    removed = 0
    for f in files:  # only ever delete inside data/genomes
        try:
            p = config.data_path(f).resolve()
            if root in p.parents and p.is_file():
                p.unlink()
                removed += 1
        except OSError:
            pass
    for d in {config.data_path(f).parent for f in files}:
        try:
            d.rmdir()  # only succeeds when empty
        except OSError:
            pass
    return {**impact, "files_removed": removed}


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
        ppx = sources.pathoplexus_resolve(w["name"], w["taxon"], w["taxon_id"])
        if ppx and ppx[0] != w["pathoplexus_organism"]:
            con.execute("UPDATE watchlist SET pathoplexus_organism=? WHERE id=?", (ppx[0], w["id"]))
        p(75, "로컬 데이터와 비교 중…")
        if not per_source["ncbi"]:
            raise SystemExit(f"NCBI returned no records for '{w['taxon']}' (taxon {w['taxon_id']}); check that the watch is a virus taxon. "
                             "Refusing to compute removals")
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


# Only these changes can alter the sequence itself; everything else is metadata.
SEQ_FIELDS = {"accession_version", "genome_length", "sequence_hash"}


def _needs_sequence(con, it) -> bool:
    """New records always need their FASTA. A changed record only does when the sequence may have
    changed (version/length/hash) or the local file is missing; a metadata-only change keeps the file."""
    if it["kind"] == "removed":
        return False
    if it["kind"] == "new":
        return True
    if SEQ_FIELDS & set(json.loads(it["diff_json"] or "{}")):
        return True
    row = con.execute("SELECT file_path FROM genome_sequence WHERE accession=?", (it["accession"],)).fetchone()
    return not (row and row["file_path"] and config.data_path(row["file_path"]).exists())


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
        need_dl = [i for i in items if _needs_sequence(con, i)]
        meta_only = len(items) - len(need_dl)
        recs = {i["accession"]: json.loads(i["new_record_json"]) for i in need_dl}
        checks.append({"name": "sequence_plan", "ok": True,
                       "detail": f"서열 다운로드 {len(need_dl)}건 · 메타데이터만 변경(서열 다운로드 생략) {meta_only}건"})
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


def _write_fasta(path, header: str, seq: str) -> None:
    """UTF-8 and LF line endings on every OS (plain write_text would produce CRLF on Windows)."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(header + "\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n")


def _read_staged(f, versioned: str):
    """Sequence from a staged FASTA written by apply(), or None if missing/mismatched."""
    try:
        lines = f.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    if not lines or not lines[0].startswith(">" + versioned + " ") and lines[0][1:] != versioned:
        return None
    return "".join(lines[1:])


def apply(changeset_id: int, progress=None, fetch_sequences: bool = True) -> dict:
    """Apply approved items atomically: files first (with backups), then one DB transaction.

    fetch_sequences=False applies metadata only: new records are added without a FASTA (fetch later with
    fetch_missing_sequences) and changes that need a new sequence stay approved for a later run."""
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

        need = {i["accession"] for i in items if _needs_sequence(con, i)}
        recs = {i["accession"]: json.loads(i["new_record_json"]) for i in items if i["accession"] in need} if fetch_sequences else {}
        no_fasta = {i["accession"] for i in items if i["kind"] == "new"} if not fetch_sequences else set()
        deferred = {i["accession"] for i in items if i["kind"] == "changed" and i["accession"] in need} if not fetch_sequences else set()
        metaonly = {i["accession"] for i in items if i["kind"] == "changed" and i["accession"] not in need}
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
                _write_fasta(f, f">{v} {recs[a].get('genome_name', '')}", s)
                meta[a] = {"md5": hashlib.md5(s.encode()).hexdigest(), "len": len(s), "gc": _gc(s), "file": f}
        if failed and not meta and recs:
            raise SystemExit(f"all downloads failed, e.g. {next(iter(failed.items()))}")
        p(72, f"서열 검증 완료 ({len(failed)}건 실패), 파일·DB 기록 중…")
        taxon_dir = config.DATA_DIR / "genomes" / config.safe_name(w["name"])
        taxon_dir.mkdir(parents=True, exist_ok=True)
        snap.mkdir(parents=True, exist_ok=True)
        now = datetime.now().isoformat(timespec="seconds")
        con.execute("BEGIN")
        next_n: dict = {}
        done_ids = []
        for it in items:
            acc, rec = it["accession"], json.loads(it["new_record_json"])
            if acc in failed or acc in deferred:
                continue
            rec = {k: ("; ".join(map(str, v)) if isinstance(v, list) else v) for k, v in rec.items()}
            if it["kind"] == "changed":  # a value the sources no longer report is cleared, not kept forever
                for k in json.loads(it["diff_json"] or "{}"):
                    if k in GENOME_COLS and k not in rec:
                        rec[k] = None
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
            if acc in no_fasta:  # metadata first; the FASTA can be fetched later
                ncols = [c for c in GENOME_COLS if c in rec]
                con.execute(f"INSERT INTO genome(genome_id,watch_id,accession,{','.join(ncols)},date_inserted,date_modified)"
                            f" VALUES(?,?,?,{','.join('?' * len(ncols))},?,?)",
                            [gid, w["id"], acc] + [rec[c] for c in ncols] + [now, now])
                for f, src in json.loads(it["provenance_json"] or "{}").items():
                    con.execute("INSERT OR REPLACE INTO field_provenance(genome_id,field,source) VALUES(?,?,?)", (gid, f, src))
                done_ids.append(it["id"])
                continue
            if acc in metaonly and old:
                mcols = [c for c in GENOME_COLS if c in rec]
                con.execute(f"UPDATE genome SET {','.join(c + '=?' for c in mcols)}, date_modified=? WHERE accession=?",
                            [rec[c] for c in mcols] + [now, acc])
                for f, src in json.loads(it["provenance_json"] or "{}").items():
                    con.execute("INSERT OR REPLACE INTO field_provenance(genome_id,field,source) VALUES(?,?,?)", (gid, f, src))
                done_ids.append(it["id"])
                continue
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
                        " VALUES(?,?,?,?,?,?,?)", (acc, gid, rec.get("accession_version"), m["len"], m["gc"], md5, config.rel_data_path(path)))
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
        return {"applied": len(done_ids), "failed": failed, "remaining": left, "metadata_only": len(metaonly - set(failed)),
                "no_fasta": len(no_fasta - set(failed)), "deferred": len(deferred)}
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



def fetch_missing_sequences(watch_id: int, genome_ids: list[str] | None = None, progress=None) -> dict:
    """Download FASTA for genomes that were applied metadata-only. Committed per batch, so an
    interrupted run simply continues with what is still missing."""
    p = progress or (lambda *a: None)
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        rows = con.execute(
            "SELECT g.genome_id, g.accession, g.accession_version, g.genome_name, g.genome_length FROM genome g "
            "LEFT JOIN genome_sequence s ON s.accession=g.accession WHERE g.watch_id=? "
            "AND IFNULL(g.genome_status,'')!='Removed' AND (s.file_path IS NULL OR s.file_path='') ORDER BY g.genome_id",
            (watch_id,)).fetchall()
        if genome_ids:
            wanted = set(genome_ids)
            rows = [r for r in rows if r["genome_id"] in wanted]
        if not rows:
            return {"fetched": 0, "failed": {}, "total": 0}
        taxon_dir = config.DATA_DIR / "genomes" / config.safe_name(w["name"])
        taxon_dir.mkdir(parents=True, exist_ok=True)
        failed: dict[str, str] = {}
        fetched = 0
        for start in range(0, len(rows), 100):
            batch = rows[start:start + 100]
            p(5 + 90 * start // len(rows), f"서열 다운로드 중… {start}/{len(rows)}")
            ver = {r["accession"]: _versioned(r["accession"], {"accession_version": r["accession_version"]}) for r in batch}
            try:
                got = sources.download_ncbi_fasta(list(ver.values()))
            except Exception as e:
                for r in batch:
                    failed[r["accession"]] = f"download error: {e}"
                continue
            for r in batch:
                a = r["accession"]
                s = got.get(ver[a])
                if not s:
                    failed[a] = "not returned by NCBI"
                    continue
                if r["genome_length"] and len(s) != r["genome_length"]:
                    failed[a] = f"length {len(s)} != expected {r['genome_length']}"
                    continue
                path = taxon_dir / f"{r['genome_id']}.fna"
                tmp = path.with_suffix(".fna.tmp")
                _write_fasta(tmp, f">{ver[a]} {r['genome_name'] or ''}", s)
                os.replace(tmp, path)
                md5 = hashlib.md5(s.encode()).hexdigest()
                con.execute("INSERT OR REPLACE INTO genome_sequence(accession,genome_id,version,length,gc_content,md5,file_path)"
                            " VALUES(?,?,?,?,?,?,?)", (a, r["genome_id"], r["accession_version"], len(s), _gc(s), md5, config.rel_data_path(path)))
                con.execute("UPDATE genome SET sequence_md5=? WHERE accession=?", (md5, a))
                fetched += 1
            con.commit()
        audit(con, "fetch_sequences", f"{w['name']}: {fetched} fetched, {len(failed)} failed")
        con.commit()
        p(100, "완료")
        return {"fetched": fetched, "failed": failed, "total": len(rows)}
    finally:
        con.close()
