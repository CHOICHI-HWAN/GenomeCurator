"""Sequence-level verification: do NCBI, BV-BRC and Pathoplexus hold the *same* sequence for an accession?

Metadata can disagree for many harmless reasons, but identical sequence bodies prove that two records are the
same entity. Only accessions whose conflict could stem from a different sequence (length / version / hash) are
checked by default; the user may also pick items. Sequences are downloaded in small batches, compared after
removing headers/whitespace and upper-casing (SHA-256), and discarded - only the verdict is stored.
NCBI is the reference (it defines the accession set)."""
import hashlib
import json
import os
from datetime import datetime
from urllib.parse import quote

from . import config, sources
from .db import connect

SEQ_FIELDS = {"genome_length", "accession_version", "sequence_hash"}   # conflicts that may be a different sequence
BATCH = 20
MAX_ACCESSIONS = int(os.environ.get("GC_SEQCHECK_MAX", "2000"))
BVBRC_SEQ = "https://www.bv-brc.org/api/genome_sequence/"


def _norm(seq: str) -> str:
    return "".join(seq.split()).upper()


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


_RC = str.maketrans("ACGTUMRWSYKVHDBN", "TGCAAKYWSRMBDHVN")


def _revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def compare(ref: str, other: str) -> dict:
    """ref/other are normalized sequences -> verdict with where they start/stop differing."""
    if ref == other:
        return {"status": "same", "length": len(other), "ref_length": len(ref), "prefix": len(ref), "suffix": 0, "detail": ""}
    if len(ref) == len(other) and ref == _revcomp(other):
        return {"status": "diff", "length": len(other), "ref_length": len(ref), "prefix": 0, "suffix": 0,
                "detail": "서열 자체는 같지만 방향이 반대(역상보)입니다"}
    n = min(len(ref), len(other))
    p = 0
    while p < n and ref[p] == other[p]:
        p += 1
    s = 0
    while s < n - p and ref[-1 - s] == other[-1 - s]:
        s += 1
    mid_r, mid_o = len(ref) - p - s, len(other) - p - s
    detail = f"앞 {p:,} bp·뒤 {s:,} bp 동일, 가운데 NCBI {mid_r:,} bp ↔ 이 소스 {mid_o:,} bp가 다름"
    if len(ref) != len(other):
        detail = f"길이 {len(ref):,} ↔ {len(other):,} bp. " + detail
    return {"status": "diff", "length": len(other), "ref_length": len(ref), "prefix": p, "suffix": s, "detail": detail}


def targets(watch_id: int, accessions: list[str] | None = None) -> list[str]:
    """Accessions worth a sequence check: the chosen ones, else those with a length/version/hash conflict."""
    con = connect()
    try:
        if accessions:
            return sorted(set(accessions))[:MAX_ACCESSIONS]
        ph = ",".join("?" * len(SEQ_FIELDS))
        rows = con.execute(f"SELECT DISTINCT accession FROM source_conflict WHERE watch_id=? AND kind='conflict' AND field IN ({ph}) "
                           "ORDER BY accession", (watch_id, *sorted(SEQ_FIELDS))).fetchall()
        return [r[0] for r in rows]
    finally:
        con.close()


def _fetch_bvbrc(c, genome_ids: dict[str, str]) -> dict[str, str]:
    """{accession: normalized sequence} for BV-BRC genome ids {accession: genome_id}."""
    out: dict[str, list[str]] = {}
    inv = {v: k for k, v in genome_ids.items()}
    ids = list(inv)
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        url = f"{BVBRC_SEQ}?in(genome_id,({','.join(quote(x) for x in chunk)}))&select(genome_id,sequence_id,length,sequence)&limit(1000)"
        r = sources._get(c, url, headers={"Accept": "application/json"})
        r.raise_for_status()
        for row in r.json():
            if row.get("sequence"):
                out.setdefault(inv.get(row["genome_id"], ""), []).append(row["sequence"])
    return {a: _norm("".join(v)) for a, v in out.items() if a}


def _fetch_ppx(c, slug: str, pp_acc: str | None, acc: str, all_revisions: bool = False) -> dict[str, str]:
    """{Pathoplexus accessionVersion: normalized sequence}: the latest revision (or all revisions of the INSDC accession)."""
    params = {"insdcAccessionBase": acc} if all_revisions or not pp_acc else {"accessionVersion": pp_acc}
    r = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/unalignedNucleotideSequences", params=params, timeout=300)
    r.raise_for_status()
    return {k: _norm(v) for k, v in sources.parse_fasta(r.text).items()}


# ------------------------------------------------------------------ dry-run (plan)
def _ncbi_lengths(c, accs_ver: list[str]) -> dict[str, int]:
    """{accession.version: length} from the light dataset_report (no sequence). Absent = no such record."""
    out: dict[str, int] = {}
    for i in range(0, len(accs_ver), 100):
        token = None
        while True:   # the report is paged (20 per page by default): without this most records look "missing"
            params = {"page_size": 1000, **({"page_token": token} if token else {})}
            j = sources._ncbi_get(c, f"{config.NCBI_BASE}/accession/{','.join(accs_ver[i:i + 100])}/dataset_report", **params).json()
            for rep in j.get("reports", []):
                out[rep["accession"]] = rep.get("length")
            token = j.get("next_page_token")
            if not token:
                break
    return out


def _bvbrc_lengths(c, genome_ids: list[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for i in range(0, len(genome_ids), 50):
        chunk = genome_ids[i:i + 50]
        url = f"{BVBRC_SEQ}?in(genome_id,({','.join(quote(x) for x in chunk)}))&select(genome_id,length)&limit(2000)"
        r = sources._get(c, url, headers={"Accept": "application/json"})
        r.raise_for_status()
        for row in r.json():
            out[row["genome_id"]] = out.get(row["genome_id"], 0) + (row.get("length") or 0)
    return out


def _ppx_lengths(c, slug: str, bases: list[str]) -> dict[str, int]:
    """{insdc base accession: length of the LATEST Pathoplexus revision}."""
    out: dict[str, int] = {}
    for i in range(0, len(bases), 100):
        params = [("insdcAccessionBase", b) for b in bases[i:i + 100]] + [
            ("versionStatus", "LATEST_VERSION"), ("fields", "accessionVersion,insdcAccessionBase,length"), ("limit", "100000")]
        r = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details", params=params, timeout=300)
        r.raise_for_status()
        for row in r.json().get("data", []):
            out[row["insdcAccessionBase"]] = row.get("length")
    return out


def build_plan(watch_id: int, accessions: list[str] | None = None, progress=None) -> dict:
    """Dry-run: confirm every record/link exists and learn its length WITHOUT downloading sequences, then decide per
    (accession, source): `download` (same length as NCBI: only the body can tell), `diff_len` (length differs: already
    proven different, nothing to download), `missing` (source has no such record), `error`, or `done` (verified before)."""
    p = progress or (lambda *a: None)
    accs = targets(watch_id, accessions)
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        if not w:
            raise ValueError("워치를 찾을 수 없습니다.")
        want = set(accs)
        recs: dict[str, dict[str, dict]] = {}
        for src, acc, pj in con.execute("SELECT source, accession, payload_json FROM source_record WHERE watch_id=?", (watch_id,)):
            if acc in want:
                recs.setdefault(acc, {})[src] = json.loads(pj)
        prev = {(a, s): (st, ln, rl) for a, s, st, ln, rl in con.execute(
            "SELECT accession, source, status, length, ref_length FROM sequence_check WHERE watch_id=?", (watch_id,))}
    finally:
        con.close()
    slug = w["pathoplexus_organism"]
    errs: dict[str, str] = {}
    with sources._client() as c:
        p(5, "NCBI 링크·길이 확인 중…")
        ncbi_ver = {a: recs.get(a, {}).get("ncbi", {}).get("accession_version") for a in accs}
        try:
            n_len = _ncbi_lengths(c, [f"{a}.{v}" for a, v in ncbi_ver.items() if v])
        except Exception as e:
            n_len, errs["ncbi"] = {}, f"NCBI 조회 실패: {str(e)[:120]}"
        p(35, "BV-BRC 링크·길이 확인 중…")
        bv_id = {a: recs[a]["bvbrc"].get("bvbrc_genome_id") for a in accs if recs.get(a, {}).get("bvbrc", {}).get("bvbrc_genome_id")}
        try:
            b_len = _bvbrc_lengths(c, sorted(set(bv_id.values()))) if bv_id else {}
        except Exception as e:
            b_len, errs["bvbrc"] = {}, f"BV-BRC 조회 실패: {str(e)[:120]}"
        p(65, "Pathoplexus 링크·길이 확인 중…")
        pp_has = [a for a in accs if recs.get(a, {}).get("pathoplexus")]
        try:
            p_len = _ppx_lengths(c, slug, pp_has) if (pp_has and slug) else {}
        except Exception as e:
            p_len, errs["pathoplexus"] = {}, f"Pathoplexus 조회 실패: {str(e)[:120]}"
    items: dict[str, dict] = {}
    for a in accs:
        v = ncbi_ver.get(a)
        ref_len = n_len.get(f"{a}.{v}") if v else None
        row = {"ncbi_version": v, "ncbi_length": ref_len, "sources": {}}
        if ref_len is None:
            row["sources"]["ncbi"] = {"state": "error" if "ncbi" in errs else "missing",
                                      "detail": errs.get("ncbi") or f"NCBI에 {a}.{v} 레코드가 없음"}
        items[a] = row
        if ref_len is None:
            continue
        cand = {"bvbrc": (a in bv_id, b_len.get(bv_id.get(a)) if a in bv_id else None),
                "pathoplexus": (a in pp_has and bool(slug), p_len.get(a))}
        for src, (present, ln) in cand.items():
            if not present:
                continue
            if src in errs:
                row["sources"][src] = {"state": "error", "detail": errs[src]}
            elif ln is None:
                row["sources"][src] = {"state": "missing", "detail": f"{ {'bvbrc': 'BV-BRC', 'pathoplexus': 'Pathoplexus'}[src] }에 이 서열이 없음"}
            elif ln != ref_len:
                row["sources"][src] = {"state": "diff_len", "length": ln, "detail": f"길이 {ref_len:,} ↔ {ln:,} bp: 길이가 달라 서열이 다름이 확정 (다운로드 생략)"}
            elif prev.get((a, src), ("",))[0] == "same" and prev[(a, src)][1] == ln and prev[(a, src)][2] == ref_len:
                row["sources"][src] = {"state": "done", "length": ln, "detail": "이전에 검증됨 (재검증 생략)"}
            else:
                row["sources"][src] = {"state": "download", "length": ln, "detail": ""}
    summ: dict = {"accessions": len(accs), "sources": {}, "requests": 0, "bytes": 0}
    need_ncbi = set()
    for a, row in items.items():
        for src, d in row["sources"].items():
            if src == "ncbi":
                summ["sources"].setdefault("ncbi", {}).setdefault(d["state"], 0)
                summ["sources"]["ncbi"][d["state"]] += 1
                continue
            summ["sources"].setdefault(src, {}).setdefault(d["state"], 0)
            summ["sources"][src][d["state"]] += 1
            if d["state"] == "download":
                summ["bytes"] += d["length"] * 2     # source sequence + NCBI reference of the same length
                need_ncbi.add(a)
    dl = {s: sum(1 for r in items.values() if r["sources"].get(s, {}).get("state") == "download") for s in ("bvbrc", "pathoplexus")}
    summ["requests"] = -(-len(need_ncbi) // 100) + -(-dl["bvbrc"] // BATCH) + dl["pathoplexus"]
    summ["seconds"] = round(summ["requests"] * max(config.REQUEST_DELAY_S, 0.3) + summ["bytes"] / 5e6)
    summ["avoided_bytes"] = sum((d.get("length") or 0) * 2 for r in items.values() for s, d in r["sources"].items()
                                if s != "ncbi" and d["state"] in ("diff_len", "done"))
    plan = {"watch_id": watch_id, "created": datetime.now().isoformat(timespec="seconds"), "summary": summ, "items": items}
    con = connect()
    try:
        con.execute("INSERT OR REPLACE INTO sequence_plan(watch_id, created, plan_json) VALUES(?,?,?)",
                    (watch_id, plan["created"], json.dumps(plan, ensure_ascii=False)))
        con.commit()
    finally:
        con.close()
    p(100, "계획 완료")
    return plan


def load_plan(watch_id: int) -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT plan_json FROM sequence_plan WHERE watch_id=?", (watch_id,)).fetchone()
        return json.loads(r[0]) if r else None
    finally:
        con.close()


def execute(watch_id: int, progress=None) -> dict:
    """Carry out the saved plan: record the verdicts that need no download, download only the `download` items."""
    p = progress or (lambda *a: None)
    plan = load_plan(watch_id)
    if not plan:
        raise ValueError("먼저 계획(dry-run)을 만들어 주세요.")
    items = plan["items"]
    now = datetime.now().isoformat(timespec="seconds")
    tally = {"same": 0, "diff": 0, "missing": 0, "error": 0, "skipped": 0}
    rows: list[tuple] = []

    def add(acc, src, status, length=None, sha=None, ref_len=None, prefix=None, suffix=None, detail=""):
        rows.append((watch_id, acc, src, status, length, sha, ref_len, prefix, suffix, detail, now))
        if src != "ncbi":
            tally[status] += 1

    todo: list[str] = []
    for a, row in items.items():
        for src, d in row["sources"].items():
            if src == "ncbi":
                add(a, "ncbi", d["state"], detail=d["detail"])
            elif d["state"] == "diff_len":
                add(a, src, "diff", d["length"], None, row["ncbi_length"], None, None, d["detail"])
            elif d["state"] in ("missing", "error"):
                add(a, src, d["state"], detail=d["detail"])
            elif d["state"] == "done":
                tally["skipped"] += 1
            elif d["state"] == "download":
                if a not in todo:
                    todo.append(a)
    con = connect()
    try:
        w = con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone()
        recs: dict[str, dict[str, dict]] = {}
        want = set(todo)
        for src, acc, pj in con.execute("SELECT source, accession, payload_json FROM source_record WHERE watch_id=?", (watch_id,)):
            if acc in want:
                recs.setdefault(acc, {})[src] = json.loads(pj)
    finally:
        con.close()
    slug = w["pathoplexus_organism"]
    with sources._client() as c:
        for i in range(0, len(todo), BATCH):
            batch = todo[i:i + BATCH]
            p(int(100 * i / max(len(todo), 1)), f"서열 내려받아 비교 중… {i}/{len(todo)}")
            ver = {a: items[a]["ncbi_version"] for a in batch}
            try:
                ncbi = {k.rsplit(".", 1)[0]: _norm(s) for k, s in sources.download_ncbi_fasta([f"{a}.{v}" for a, v in ver.items()]).items()}
                ncbi_err = ""
            except Exception as e:
                ncbi, ncbi_err = {}, f"NCBI 다운로드 실패: {str(e)[:120]}"
            bv_ids = {a: recs[a]["bvbrc"]["bvbrc_genome_id"] for a in batch if items[a]["sources"].get("bvbrc", {}).get("state") == "download"}
            try:
                bv, bv_err = (_fetch_bvbrc(c, bv_ids) if bv_ids else {}), ""
            except Exception as e:
                bv, bv_err = {}, f"BV-BRC 다운로드 실패: {str(e)[:120]}"
            for a in batch:
                ref = ncbi.get(a)
                if ref is None:
                    add(a, "ncbi", "error" if ncbi_err else "missing", detail=ncbi_err or "NCBI에서 서열을 받지 못함")
                    continue
                add(a, "ncbi", "ref", len(ref), _sha(ref), len(ref))
                if a in bv_ids:
                    if bv_err:
                        add(a, "bvbrc", "error", detail=bv_err)
                    elif a in bv:
                        v = compare(ref, bv[a])
                        add(a, "bvbrc", v["status"], v["length"], _sha(bv[a]), v["ref_length"], v["prefix"], v["suffix"], v["detail"])
                    else:
                        add(a, "bvbrc", "missing", detail="BV-BRC가 서열을 주지 않음")
                if items[a]["sources"].get("pathoplexus", {}).get("state") == "download":
                    try:
                        got = _fetch_ppx(c, slug, recs[a]["pathoplexus"].get("pathoplexus_accession"), a)
                        seq = next(iter(got.values()), None)
                        if seq is None:
                            add(a, "pathoplexus", "missing", detail="Pathoplexus가 서열을 주지 않음")
                        else:
                            v = compare(ref, seq)
                            if v["status"] == "diff":
                                for k, s2 in _fetch_ppx(c, slug, None, a, all_revisions=True).items():
                                    if s2 == ref:
                                        v["detail"] += f" · 수정본 {k}이(가) NCBI와 동일"
                                        break
                            add(a, "pathoplexus", v["status"], v["length"], _sha(seq), v["ref_length"], v["prefix"], v["suffix"], v["detail"])
                    except Exception as e:
                        add(a, "pathoplexus", "error", detail=f"Pathoplexus 다운로드 실패: {str(e)[:120]}")
    con = connect()
    try:
        con.executemany("INSERT OR REPLACE INTO sequence_check(watch_id,accession,source,status,length,sha256,ref_length,prefix_same,suffix_same,detail,checked_at)"
                        " VALUES(?,?,?,?,?,?,?,?,?,?,?)", rows)
        con.execute("DELETE FROM sequence_plan WHERE watch_id=?", (watch_id,))
        con.execute("INSERT INTO audit_log(action, detail) VALUES('seqcheck', ?)",
                    (f"{w['name']}: {len(items)} accessions ({len(todo)} downloaded), " + ", ".join(f"{k} {v}" for k, v in tally.items()),))
        con.commit()
    finally:
        con.close()
    p(100, "서열 검증 완료")
    return {"checked": len(items), "downloaded": len(todo), **tally}


def run(watch_id: int, accessions: list[str] | None = None, progress=None) -> dict:
    """Plan and execute in one go (scripts/tests); the web UI shows the plan for confirmation in between."""
    build_plan(watch_id, accessions, progress)
    return execute(watch_id, progress)


def results(watch_id: int, status: str = "", limit: int = 300) -> list[dict]:
    con = connect()
    try:
        q, a = "SELECT * FROM sequence_check WHERE watch_id=? AND source!='ncbi'", [watch_id]
        if status:
            q += " AND status=?"
            a.append(status)
        q += " ORDER BY CASE status WHEN 'diff' THEN 0 WHEN 'error' THEN 1 WHEN 'missing' THEN 2 ELSE 3 END, accession LIMIT ?"
        a.append(limit)
        return [dict(r) for r in con.execute(q, a)]
    finally:
        con.close()


def summary(watch_id: int) -> dict:
    con = connect()
    try:
        by = {r[0]: r[1] for r in con.execute("SELECT status, COUNT(*) FROM sequence_check WHERE watch_id=? AND source!='ncbi' GROUP BY status", (watch_id,))}
        n = con.execute("SELECT COUNT(DISTINCT accession) FROM sequence_check WHERE watch_id=?", (watch_id,)).fetchone()[0]
        return {"accessions": n, "same": by.get("same", 0), "diff": by.get("diff", 0), "missing": by.get("missing", 0),
                "error": by.get("error", 0), "targets": len(targets(watch_id))}
    finally:
        con.close()


def verdicts(watch_id: int) -> dict[str, dict]:
    """{accession: {'same': bool (every compared source identical), 'diff': [sources], 'sources': {src: status}}}"""
    con = connect()
    try:
        out: dict[str, dict] = {}
        for acc, src, st in con.execute("SELECT accession, source, status FROM sequence_check WHERE watch_id=? AND source!='ncbi'", (watch_id,)):
            d = out.setdefault(acc, {"sources": {}, "diff": [], "same": True})
            d["sources"][src] = st
            if st == "diff":
                d["diff"].append(src)
            if st != "same":
                d["same"] = False
        return out
    finally:
        con.close()
