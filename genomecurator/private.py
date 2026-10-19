"""The private store: local, unpublished genomes that belong in a final analysis dataset but must never mix with the public curated set.

Why a store of its own (and not rows in `genome`):
* the curated set is defined by NCBI (INSDC) accessions - a private genome would be reported as "gone from NCBI" by the next check;
* the public page, the Superset view, TSV downloads and backups are all built from `genome`: keeping private data out of that table keeps it
  out of every one of them by construction.

Input looks like what the Nextstrain mpox workflow uses: a metadata TSV (columns named as in the export - `strain`, `date`, `country`,
`division`, `location`, `host`, `clade`, ... plus an id column `PPX_accession` or `accession`) and a FASTA whose ids are those ids. Importing
is dry-run first (`check`): structural problems are errors that stop the import (missing / duplicate / colliding ids, rows without a sequence,
sequences without a row, bad sequence characters); format quirks (partial dates, a missing region, a length that differs) are normalised
and reported as warnings. Sequences are stored under `data/private/<virus>/`. QC (sequence statistics and Nextclade) runs on them like on
public genomes, under the accession `private:<id>`.

Private rows leave the tool only through an export that asks for them explicitly (`export --with-private`) and into separate files
(`metadata-private.tsv`, `sequences-private.fasta`) - the form the workflow's own `private_metadata` / `private_sequences` expect."""
import csv
import io
import json
import re
from datetime import datetime
from pathlib import Path

from . import config, export, log, reconcile
from .db import audit, connect

_log = log.get("private")
ID_COLUMNS = ("PPX_accession", "accession")
REQUIRED = ("strain", "date")
ACGT = set("ACGTRYSWKMBDHVN-")        # IUPAC nucleotide codes (and gaps); anything else is not a nucleotide sequence
RESERVED = {"qc", "source_nextstrain", "source_private"}


def prefix(pid: str) -> str:
    return f"private:{pid}"


def private_dir(virus: str) -> Path:
    return config.DATA_DIR / "private" / export.slug(virus)


def _safe(pid: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", pid)


def _read_tsv(text: str) -> tuple[list[str], list[dict]]:
    rd = csv.DictReader(io.StringIO(text.lstrip("﻿")), delimiter="\t")
    cols = [c.strip() for c in (rd.fieldnames or [])]
    rows = []
    for r in rd:
        rows.append({(k or "").strip(): (v or "").strip() for k, v in r.items() if k is not None})
    return cols, rows


def _read_fasta(text: str) -> tuple[dict[str, str], list[str]]:
    seqs: dict[str, str] = {}
    dups: list[str] = []
    cur, buf = None, []
    for ln in text.splitlines():
        if ln.startswith(">"):
            if cur is not None:
                (dups.append(cur) if cur in seqs else None)
                seqs[cur] = "".join(buf).upper()
            cur, buf = ln[1:].split()[0] if ln[1:].strip() else "", []
        elif cur is not None:
            buf.append(ln.strip())
    if cur is not None:
        (dups.append(cur) if cur in seqs else None)
        seqs[cur] = "".join(buf).upper()
    return seqs, dups


def _public_ids(con, watch_id: int) -> set[str]:
    ids = set()
    for r in con.execute("SELECT accession, accession_version, pathoplexus_accession FROM genome WHERE watch_id=?", (watch_id,)):
        ids.add(r[0])
        if r[1]:
            ids.add(f"{r[0]}.{r[1]}")
        if r[2]:
            ids.add(str(r[2]).split(".")[0])
    return ids


def check(watch_id: int, metadata_text: str, fasta_text: str, id_column: str | None = None) -> dict:
    """Dry run: parse, validate and normalise. Nothing is written. -> {ok, errors, warnings, rows (normalised), n, new, update, id_column}"""
    con = connect()
    try:
        w = con.execute("SELECT name FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("바이러스를 찾을 수 없습니다.")
        public = _public_ids(con, watch_id)
        have = {r[0] for r in con.execute("SELECT id FROM private_genome WHERE watch_id=?", (watch_id,))}
    finally:
        con.close()
    errors, warnings = [], []
    cols, rows = _read_tsv(metadata_text)
    seqs, dup_seq = _read_fasta(fasta_text)
    idc = id_column or next((c for c in ID_COLUMNS if c in cols), None)
    if not rows:
        errors.append("메타데이터에 행이 없습니다.")
    if not idc or idc not in cols:
        errors.append(f"id 열이 없습니다. 메타데이터에 {' 또는 '.join(ID_COLUMNS)} 열이 있어야 합니다 (FASTA 헤더와 같은 값).")
    for need in REQUIRED:
        if need not in cols:
            errors.append(f"필수 열 `{need}`가 없습니다.")
    if errors:
        return {"ok": False, "errors": errors, "warnings": warnings, "rows": [], "n": 0, "new": 0, "update": 0, "id_column": idc}
    seen: dict[str, int] = {}
    out, bad_seq = [], 0
    for i, r in enumerate(rows, start=2):                         # line numbers as in the file (header = line 1)
        pid = r.get(idc, "")
        if not pid:
            errors.append(f"{i}행: id가 비어 있습니다.")
            continue
        if not re.fullmatch(r"[A-Za-z0-9._-]+", pid):
            errors.append(f"{i}행: id `{pid}`에 허용되지 않는 문자가 있습니다 (영문·숫자·`._-`만).")
            continue
        if pid in seen:
            errors.append(f"{i}행: id `{pid}`가 {seen[pid]}행과 겹칩니다.")
            continue
        seen[pid] = i
        if pid in public:
            errors.append(f"{i}행: id `{pid}`는 공개 데이터(NCBI/Pathoplexus)의 accession과 같습니다. 비공개 id는 겹치면 안 됩니다.")
            continue
        s = seqs.get(pid)
        if s is None:
            errors.append(f"{i}행: id `{pid}`에 해당하는 FASTA 서열이 없습니다.")
            continue
        if not s:
            errors.append(f"{i}행: `{pid}`의 서열이 비어 있습니다.")
            continue
        badc = sorted(set(s) - ACGT)
        if badc:
            bad_seq += 1
            errors.append(f"{i}행: `{pid}`의 서열에 뉴클레오타이드가 아닌 문자가 있습니다: {''.join(badc[:8])}")
            continue
        rec = dict(r)
        d = (r.get("date") or "").strip()
        if not r.get("strain"):
            rec["strain"] = pid
            warnings.append(f"{i}행: `strain`이 비어 id로 채웠습니다.")
        if d:
            if re.fullmatch(r"\d{4}-(\d{2}|XX)-(\d{2}|XX)", d):
                rec["date"] = d
            else:
                nd = export.std_date(d)
                if nd:
                    rec["date"] = nd
                    if nd != d:
                        warnings.append(f"{i}행: 날짜 `{d}` → `{nd}`(Nextstrain 형식)로 바꿨습니다.")
                else:
                    rec["date"] = ""
                    warnings.append(f"{i}행: 날짜 `{d}`를 해석하지 못해 비웠습니다.")
        else:
            warnings.append(f"{i}행: 날짜가 비어 있습니다 (시간 기반 분석에서 제외될 수 있습니다).")
        if not r.get("region") and r.get("country"):
            rg = export.region_of(r["country"])
            if rg:
                rec["region"] = rg
            else:
                warnings.append(f"{i}행: 국가 `{r['country']}`의 대륙(region)을 알 수 없어 비웠습니다.")
        elif not r.get("country"):
            warnings.append(f"{i}행: 국가가 비어 있습니다.")
        if r.get("length") and r["length"].isdigit() and int(r["length"]) != len(s):
            warnings.append(f"{i}행: `length` {r['length']}과 서열 길이 {len(s):,}가 달라 서열 길이로 바꿨습니다.")
        rec["length"] = str(len(s))
        for k in list(rec):
            if k in RESERVED:
                rec.pop(k)
        rec["_id"], rec["_seq"] = pid, s
        out.append(rec)
    for pid in dup_seq:
        errors.append(f"FASTA에 id `{pid}`가 두 번 이상 있습니다.")
    extra = sorted(set(seqs) - set(seen))
    if extra:
        errors.append(f"메타데이터에 없는 FASTA 서열 {len(extra)}개: {', '.join(extra[:5])}{' …' if len(extra) > 5 else ''}")
    return {"ok": not errors, "errors": errors, "warnings": warnings, "rows": out, "n": len(rows), "new": sum(1 for r in out if r["_id"] not in have),
            "update": sum(1 for r in out if r["_id"] in have), "id_column": idc, "bad_seq": bad_seq}


def import_texts(watch_id: int, metadata_text: str, fasta_text: str, source_name: str = "", id_column: str | None = None) -> dict:
    """Import after a clean check. Errors stop everything (all or nothing); an existing id is updated."""
    rep = check(watch_id, metadata_text, fasta_text, id_column)
    if not rep["ok"]:
        raise ValueError("가져올 수 없습니다:\n" + "\n".join(rep["errors"][:20]))
    con = connect()
    try:
        w = con.execute("SELECT name FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        d = private_dir(w["name"])
        d.mkdir(parents=True, exist_ok=True)
        now = datetime.now().isoformat(timespec="seconds")
        import hashlib
        for r in rep["rows"]:
            pid, seq = r["_id"], r["_seq"]
            meta = {k: v for k, v in r.items() if not k.startswith("_")}
            f = d / f"{_safe(pid)}.fna"
            f.write_text(f">{pid}\n" + "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n", encoding="utf-8", newline="\n")
            rel = f.relative_to(config.DATA_DIR).as_posix()
            md5 = hashlib.md5(seq.encode()).hexdigest()
            old = con.execute("SELECT imported FROM private_genome WHERE watch_id=? AND id=?", (watch_id, pid)).fetchone()
            con.execute("INSERT OR REPLACE INTO private_genome(watch_id,id,metadata_json,file_path,seq_md5,length,source_file,imported,updated) VALUES(?,?,?,?,?,?,?,?,?)",
                        (watch_id, pid, json.dumps(meta, ensure_ascii=False), rel, md5, len(seq), source_name, old["imported"] if old else now, now))
        audit(con, "private_import", f"{w['name']}: 신규 {rep['new']}건, 갱신 {rep['update']}건 ({source_name})")
        con.commit()
    finally:
        con.close()
    _log.info("private import watch=%s new=%d update=%d", watch_id, rep["new"], rep["update"])
    return {"new": rep["new"], "update": rep["update"], "warnings": rep["warnings"]}


def list_rows(watch_id: int) -> list[dict]:
    con = connect()
    try:
        out = []
        qcs: dict[str, dict] = {}
        for q in con.execute("SELECT accession, tool, status, clade, lineage FROM qc_result WHERE watch_id=? AND accession LIKE 'private:%'", (watch_id,)):
            qcs.setdefault(q["accession"], {})[q["tool"]] = dict(q)
        for r in con.execute("SELECT id, metadata_json, length, imported, updated, seq_md5 FROM private_genome WHERE watch_id=? ORDER BY id", (watch_id,)):
            m = json.loads(r["metadata_json"] or "{}")
            q = qcs.get(prefix(r["id"]), {})
            nc, st = q.get("nextclade") or {}, q.get("seqstats") or {}
            out.append({"id": r["id"], "strain": m.get("strain"), "date": m.get("date"), "country": m.get("country"), "clade": m.get("clade") or nc.get("clade") or "",
                        "length": r["length"], "imported": r["imported"], "updated": r["updated"], "qc_stats": st.get("status", ""), "qc_nextclade": nc.get("status", ""),
                        "nextclade_clade": nc.get("clade") or "", "clade_differs": bool(m.get("clade") and nc.get("clade") and m["clade"] != nc["clade"])})
        return out
    finally:
        con.close()


def delete(watch_id: int, pid: str) -> bool:
    con = connect()
    try:
        r = con.execute("SELECT file_path FROM private_genome WHERE watch_id=? AND id=?", (watch_id, pid)).fetchone()
        if not r:
            return False
        con.execute("DELETE FROM private_genome WHERE watch_id=? AND id=?", (watch_id, pid))
        con.execute("DELETE FROM qc_result WHERE accession=?", (prefix(pid),))
        audit(con, "private_delete", f"{pid}")
        con.commit()
    finally:
        con.close()
    try:
        config.data_path(r["file_path"]).unlink()
    except OSError:
        pass
    return True


def rows_for_export(watch_id: int) -> list[dict]:
    """[{id, meta, file_path, qc:{tool: {...}}}] for the export module."""
    con = connect()
    try:
        qcs: dict[str, dict] = {}
        for q in con.execute("SELECT accession, tool, status, clade, lineage, metrics_json FROM qc_result WHERE watch_id=? AND accession LIKE 'private:%'", (watch_id,)):
            qcs.setdefault(q["accession"], {})[q["tool"]] = {"status": q["status"], "clade": q["clade"], "lineage": q["lineage"], "m": json.loads(q["metrics_json"] or "{}")}
        return [{"id": r["id"], "meta": json.loads(r["metadata_json"] or "{}"), "file_path": r["file_path"], "qc": qcs.get(prefix(r["id"]), {})}
                for r in con.execute("SELECT id, metadata_json, file_path FROM private_genome WHERE watch_id=? ORDER BY id", (watch_id,))]
    finally:
        con.close()
