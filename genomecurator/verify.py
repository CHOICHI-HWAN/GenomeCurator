"""Integrity checks on the local curated data (no network). Each check returns
{name, status: ok|warn|fail, count, detail, examples}. Used by `genomecurator verify` and the 정합성 page."""
import hashlib
import re
from datetime import date

from . import config
from .db import connect


def _fasta_seq(path) -> str | None:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    return "".join(lines[1:]) if lines and lines[0].startswith(">") else None


def run(watch_id: int | None = None, check_files: bool = True) -> list[dict]:
    con = connect()
    out: list[dict] = []
    w = "" if watch_id is None else " AND g.watch_id=%d" % int(watch_id)

    def add(name, bad, total_detail="", fail=True, examples=None):
        n = len(bad) if not isinstance(bad, int) else bad
        out.append({"name": name, "status": "ok" if n == 0 else ("fail" if fail else "warn"), "count": n,
                    "detail": total_detail, "examples": (list(bad)[:8] if not isinstance(bad, int) else examples or [])})

    try:
        q = lambda sql: [r[0] for r in con.execute(sql)]
        total = con.execute("SELECT COUNT(*) FROM genome g WHERE 1=1" + w).fetchone()[0]
        out.append({"name": "유전체 수", "status": "ok", "count": total, "detail": "검사 대상", "examples": []})
        add("watch가 없는 유전체", q("SELECT g.accession FROM genome g LEFT JOIN watchlist w ON w.id=g.watch_id WHERE w.id IS NULL" + w))
        add("유전체 없는 서열 행", q("SELECT s.accession FROM genome_sequence s LEFT JOIN genome g ON g.genome_id=s.genome_id WHERE g.genome_id IS NULL"))
        add("accession_version 누락", q("SELECT g.accession FROM genome g WHERE g.accession_version IS NULL" + w), fail=False)
        add("genome_length 이상(없음/0 이하)", q("SELECT g.accession FROM genome g WHERE (g.genome_length IS NULL OR g.genome_length<=0) AND g.genome_status!='Removed'" + w), fail=False)
        year = date.today().year
        add("수집일이 미래", [a for a, d in con.execute("SELECT g.accession, g.collection_date FROM genome g WHERE g.collection_date GLOB '[12][0-9][0-9][0-9]*'" + w)
                             if int(d[:4]) > year], fail=False)
        add("수집일이 공개일보다 늦음", q("SELECT g.accession FROM genome g WHERE length(g.collection_date)>=10 AND length(g.release_date)>=10 "
                                    "AND substr(g.collection_date,1,10) > substr(g.release_date,1,10)" + w), fail=False)
        add("수집일 형식 오류", [a for a, d in con.execute("SELECT g.accession, g.collection_date FROM genome g WHERE g.collection_date IS NOT NULL AND g.collection_date!=''" + w)
                              if not re.match(r"^\d{4}(-\d{2}(-\d{2})?)?(/\d{4}(-\d{2}(-\d{2})?)?)?$", d)], fail=False)
        add("서열 행의 길이 ≠ 유전체 길이", q("SELECT g.accession FROM genome g JOIN genome_sequence s ON s.genome_id=g.genome_id "
                                          "WHERE s.length IS NOT NULL AND g.genome_length IS NOT NULL AND s.length!=g.genome_length" + w))
        if check_files:
            missing, wrong_len, wrong_md5, no_file = [], [], [], 0
            rows = con.execute("SELECT g.accession, s.file_path, s.length, s.md5 FROM genome g JOIN genome_sequence s ON s.genome_id=g.genome_id WHERE 1=1" + w).fetchall()
            for acc, fp, ln, md5 in rows:
                if not fp:
                    no_file += 1
                    continue
                seq = _fasta_seq(config.data_path(fp))
                if seq is None:
                    missing.append(acc)
                    continue
                if ln is not None and len(seq) != ln:
                    wrong_len.append(acc)
                if md5 and hashlib.md5(seq.encode()).hexdigest() != md5:
                    wrong_md5.append(acc)
            add("FASTA 파일 없음/읽기 실패", missing)
            add("FASTA 실제 길이 ≠ 기록 길이", wrong_len)
            add("FASTA md5 ≠ 기록 md5", wrong_md5)
            out.append({"name": "FASTA 검사 범위", "status": "ok", "count": len(rows) - no_file,
                        "detail": f"파일 {len(rows) - no_file}건 검사, 서열 미수집 {total - len(rows) + no_file}건 (Get FASTA로 받을 수 있음)", "examples": []})
        add("Removed 상태(NCBI에서 사라짐)", q("SELECT g.accession FROM genome g WHERE g.genome_status='Removed'" + w), fail=False)
        add("출처(field_provenance)가 없는 유전체", q("SELECT g.accession FROM genome g WHERE NOT EXISTS (SELECT 1 FROM field_provenance f WHERE f.genome_id=g.genome_id)" + w), fail=False)
    finally:
        con.close()
    return out
