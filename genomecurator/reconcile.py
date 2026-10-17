"""Cross-source reconciliation: do NCBI, BV-BRC and Pathoplexus agree about the same accession?

For every accession that two sources both report a value for, each field is classified as
  match     identical after trimming
  cosmetic  same fact, written differently (date precision, "USA: CA" vs "USA, CA", host with age/sex, ...)
  conflict  a real disagreement (different version, length, strain, date, ...)
Fields whose meaning differs by source (genome_name, genome_status, release/update dates) are not compared.
The merge keeps using source priority; this module only measures and records the disagreement."""
import json
import re
from datetime import datetime
from itertools import combinations

from .db import connect

SOURCES = ["ncbi", "bvbrc", "pathoplexus"]
LABEL = {"ncbi": "NCBI", "bvbrc": "BV-BRC", "pathoplexus": "Pathoplexus"}
# field -> comparison kind
FIELDS = {
    "accession_version": "num", "genome_length": "num", "collection_date": "date", "isolation_country": "country",
    "geographic_location": "place", "host_name": "host", "strain": "strain", "isolation_source": "contain",
    "submitter": "tokens", "sra_accession": "set", "biosample_accession": "text", "bioproject_accession": "text",
    "clade": "text", "lineage": "text", "sequence_hash": "text",
}
EXCLUDED = {"genome_name": "소스마다 이름 규칙이 다름", "genome_status": "Complete/Partial 정의가 소스마다 다름",
            "release_date": "각 포털의 공개일이라 서로 다른 게 정상", "update_date": "각 포털의 갱신일"}
COUNTRY_ALIAS = {"zaire": "democratic republic of the congo", "drc": "democratic republic of the congo",
                 "usa": "united states", "us": "united states", "u.s.a.": "united states", "uk": "united kingdom",
                 "england": "united kingdom", "korea, south": "south korea", "republic of korea": "south korea",
                 "viet nam": "vietnam", "russian federation": "russia", "ivory coast": "cote d'ivoire"}


def _s(v) -> str:
    if isinstance(v, (list, tuple, set)):
        v = "; ".join(str(x) for x in v)
    return str(v).strip()


def _words(v: str) -> str:
    return re.sub(r"[\s:;,/_\-]+", " ", v.casefold()).strip()


def _country(v: str) -> str:
    c = v.casefold().strip()
    return COUNTRY_ALIAS.get(c, c)


MONTHS = {m: i for i, m in enumerate(["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], 1)}


def _date_parts(v: str):
    start = v.split("/")[0].strip()
    m = re.match(r"^(?:(\d{1,2})-)?([A-Za-z]{3})[a-z]*-(\d{4})$", start)   # 'Mar-2026', '11-Mar-2026'
    if m and m.group(2).lower() in MONTHS:
        mo = f"{MONTHS[m.group(2).lower()]:02d}"
        return (m.group(3), mo, f"{int(m.group(1)):02d}" if m.group(1) else None, False)
    m = re.match(r"^(\d{4})(?:-(\d{2}))?(?:-(\d{2}))?$", start)
    return (m.group(1), m.group(2), m.group(3), "/" in v) if m else None


def classify(kind: str, a, b) -> str:
    """-> 'match' | 'cosmetic' | 'conflict' for two non-empty values."""
    x, y = _s(a), _s(b)
    if x == y:
        return "match"
    if x.casefold() == y.casefold():
        return "cosmetic"
    if kind != "num" and re.sub(r"[\W_]+", "", x.casefold()) == re.sub(r"[\W_]+", "", y.casefold()):
        return "cosmetic"   # only spaces / punctuation differ ("skin, lesion" vs "skin lesion", "Homosapiens")
    if kind == "num":
        try:
            return "match" if float(x) == float(y) else "conflict"
        except ValueError:
            return "conflict"
    if kind == "date":
        p, q = _date_parts(x), _date_parts(y)
        if not p or not q:
            return "conflict"
        n = min(sum(1 for t in p[:3] if t), sum(1 for t in q[:3] if t))
        same = p[:n] == q[:n]
        return "cosmetic" if same else "conflict"   # also covers ranges (2024-01-01/2024-12-09 vs 2024)
    if kind == "country":
        return "cosmetic" if _country(x) == _country(y) else "conflict"
    if kind == "place":
        if _words(x) == _words(y):
            return "cosmetic"
        return "cosmetic" if _country(re.split(r"[:,]", x)[0]) == _country(re.split(r"[:,]", y)[0]) else "conflict"
    if kind == "host":   # "Sus scrofa (pig)" is "Sus scrofa" with its common name appended; "Homo sapiens; 32; male" carries age and sex
        bare = lambda v: re.sub(r"\s*\([^)]*\)", "", v.split(";")[0]).strip().casefold()
        return "cosmetic" if bare(x) == bare(y) else "conflict"
    if kind == "strain":
        s1, s2 = x.casefold(), y.casefold()
        short, long_ = sorted((s1, s2), key=len)
        return "cosmetic" if len(short) >= 3 and re.search(r"(^|[^a-z0-9])" + re.escape(short) + r"($|[^a-z0-9])", long_) else "conflict"
    if kind == "contain":
        s1, s2 = x.casefold(), y.casefold()
        return "cosmetic" if s1 in s2 or s2 in s1 else "conflict"
    if kind == "tokens":
        t1, t2 = set(_words(x).split()), set(_words(y).split())
        return "cosmetic" if t1 and t2 and len(t1 & t2) / min(len(t1), len(t2)) >= 0.5 else "conflict"
    if kind == "set":
        return "cosmetic" if set(re.split(r"[;, ]+", x.casefold())) == set(re.split(r"[;, ]+", y.casefold())) else "conflict"
    return "conflict"


def analyze(per_source: dict[str, dict[str, dict]]) -> tuple[dict, list[tuple]]:
    """-> (report, rows) where rows = (accession, field, kind, src_a, val_a, src_b, val_b) for every non-exact pair."""
    present = [s for s in SOURCES if per_source.get(s)]
    accs = {s: set(per_source[s]) for s in present}
    ncbi = accs.get("ncbi", set())
    union = set().union(*accs.values()) if accs else set()
    cov = {"sources": {s: len(accs[s]) for s in present}, "union": len(union),
           "all_sources": len(set.intersection(*accs.values())) if accs else 0,
           "pairs": {f"{a}|{b}": len(accs[a] & accs[b]) for a, b in combinations(present, 2)},
           # not in NCBI = never enters the curated set (NCBI defines the accession set)
           "not_in_ncbi": {s: len(accs[s] - ncbi) for s in present if s != "ncbi"},
           "only_in_ncbi": {s: len(ncbi - accs[s]) for s in present if s != "ncbi"},
           "examples_not_in_ncbi": {s: sorted(accs[s] - ncbi)[:10] for s in present if s != "ncbi"}}
    fields: dict[str, dict] = {}
    rows: list[tuple] = []
    for a, b in combinations(present, 2):
        pair = f"{a}|{b}"
        for acc in accs[a] & accs[b]:
            ra, rb = per_source[a][acc], per_source[b][acc]
            for f, kind in FIELDS.items():
                va, vb = ra.get(f), rb.get(f)
                if va in (None, "", []) or vb in (None, "", []):
                    continue
                r = classify(kind, va, vb)
                d = fields.setdefault(f, {}).setdefault(pair, {"match": 0, "cosmetic": 0, "conflict": 0})
                d[r] += 1
                if r != "match":
                    rows.append((acc, f, r, a, _s(va)[:300], b, _s(vb)[:300]))
    report = {"coverage": cov, "fields": fields, "excluded": EXCLUDED, "checked_at": datetime.now().isoformat(timespec="seconds"),
              "totals": {"cosmetic": sum(1 for r in rows if r[2] == "cosmetic"), "conflict": sum(1 for r in rows if r[2] == "conflict")}}
    return report, rows


def store(con, watch_id: int, report: dict, rows: list[tuple]) -> None:
    con.execute("DELETE FROM source_conflict WHERE watch_id=?", (watch_id,))
    con.executemany("INSERT INTO source_conflict(watch_id,accession,field,kind,src_a,val_a,src_b,val_b) VALUES(?,?,?,?,?,?,?,?)",
                    [(watch_id, r[0], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows])
    con.execute("INSERT OR REPLACE INTO reconcile_report(watch_id,checked_at,report_json) VALUES(?,?,?)",
                (watch_id, report["checked_at"], json.dumps(report, ensure_ascii=False)))


def load(watch_id: int) -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT report_json FROM reconcile_report WHERE watch_id=?", (watch_id,)).fetchone()
        return json.loads(r[0]) if r else None
    finally:
        con.close()


def conflicts(watch_id: int, field: str = "", kind: str = "conflict", limit: int = 200) -> list[dict]:
    con = connect()
    try:
        q, a = "SELECT accession,field,kind,src_a,val_a,src_b,val_b FROM source_conflict WHERE watch_id=?", [watch_id]
        if field:
            q += " AND field=?"
            a.append(field)
        if kind:
            q += " AND kind=?"
            a.append(kind)
        q += " ORDER BY field, accession LIMIT ?"
        a.append(limit)
        return [dict(r) for r in con.execute(q, a)]
    finally:
        con.close()
