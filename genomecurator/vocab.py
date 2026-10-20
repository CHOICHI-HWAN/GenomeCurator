"""The app's own dictionary of terms and the standardisation of incoming values against it.

Sources spell the same thing in many ways (`Homo sapiens`, `Homosapiens`, `H. sapiens`, `Patent`; `USA` / `United States`; `Zaire` / `Democratic Republic of the Congo`; collection dates as
`2025-05-01/2025-05-28`, `15-May-2025`, ...). Compared as they come, such differences look like conflicts and end up as different values in the curated table. Here every incoming value of
a standardised field is looked up in **our dictionary** and replaced by the standard term, **before** the sources are compared:

    host_name, isolation_country, collection_date, clade

The dictionary lives in `resources/vocab/*.tsv` (host, country, clade: term, synonyms separated by `|`, and a flag for terms that are not what the field means, e.g. a cell line in the
host field) and can be extended or overridden with files of the same name in `<data>/vocab/` (kept out of git, survives updates). Dates follow rules instead of a list.

Nothing is guessed and nothing is deleted:
* a value found in the dictionary becomes its standard term (`mapped`), already-standard values are `same`;
* a format rule that is unambiguous gives `normalized` (`15-May-2025` -> `2025-05-15`, `2025-05-01/2025-05-28` -> `2025-05` with the original kept in the record);
* an empty marker (`missing`, `null`, `n/a` ...) becomes an empty value (`empty`);
* a value the dictionary knows to be *not what the field means* (cell line, environmental sample, genus-level `Wild monkey`) or that cannot be read safely (an ambiguous `05/06/2024`,
  a broad clade `II`, a recombinant) is **kept as it is** and marked `flagged`;
* a value the dictionary does not know is kept and marked `unmapped` - the list of those is what a person decides about.

Every change is written to the table `standardization` (raw value, standard value, status, rule), so the raw value is never lost and the review screen can show the reason."""
import csv
import re
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from . import config, log

_log = log.get("vocab")
FIELDS = ("host_name", "isolation_country", "collection_date", "clade")
FILES = {"host_name": "host", "isolation_country": "country", "clade": "clade"}
EMPTY = {"", "missing", "unknown", "unspecified", "n/a", "na", "none", "null", "nan", "not applicable", "not available", "not collected", "not provided", "not determined", "-", "--", "?"}
MONTHS = {m: i for i, m in enumerate(["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"], 1)}
_CACHE: dict = {}


@dataclass
class Std:
    value: str
    status: str            # same | mapped | normalized | empty | flagged | unmapped | nodict
    rule: str = ""
    note: str = ""


def key(s: str) -> str:
    """Comparison key: no accents, no case, punctuation as a space (`Cote d'Ivoire` = `Côte d’Ivoire`, `H. sapiens` = `h sapiens`)."""
    s = str(s).translate({0x2018: "'", 0x2019: "'", 0x02BC: "'", 0x2032: "'", 0x60: "'", 0xB4: "'", 0x2010: "-", 0x2011: "-", 0x2012: "-", 0x2013: "-", 0x2014: "-"})   # typographic quotes / dashes are the plain ones
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode().lower()
    return " ".join(re.sub(r"[^a-z0-9]+", " ", s).split())


def vocab_dirs() -> list[Path]:
    return [Path(__file__).parent / "resources" / "vocab", config.DATA_DIR / "vocab"]


def _rows(fname: str):
    for d in vocab_dirs():
        f = d / f"{fname}.tsv"
        if f.exists():
            for r in csv.DictReader((l for l in open(f, encoding="utf-8") if not l.startswith("#")), delimiter="\t"):
                if (r.get("term") or "").strip():
                    yield d, r


def load(field: str) -> dict:
    """{(scope, key): entry} of the field's dictionary: shipped terms first, then the user's file on top (its synonyms win)."""
    if field in _CACHE and _CACHE[field][0] == _stamp(field):
        return _CACHE[field][1]
    out: dict = {}
    for _d, r in _rows(FILES[field]):
        scope = (r.get("scope") or "*").strip() or "*"
        term = r["term"].strip()
        e = {"term": term, "flag": (r.get("flag") or "").strip(), "note": (r.get("note") or "").strip(), "iso3": (r.get("iso3") or "").strip(), "taxon_id": (r.get("taxon_id") or "").strip(), "scope": scope}
        for k in [term] + [x for x in (r.get("synonyms") or "").split("|") if x.strip()]:
            out[(scope, key(k))] = e
    _CACHE[field] = (_stamp(field), out)
    return out


def _stamp(field: str):
    return tuple((str(d / f"{FILES[field]}.tsv"), (d / f"{FILES[field]}.tsv").stat().st_mtime if (d / f"{FILES[field]}.tsv").exists() else 0) for d in vocab_dirs())


def terms(field: str) -> dict[str, dict]:
    """The standard terms of a field: {term: entry}."""
    return {e["term"]: e for e in load(field).values()}


def scope_for(*names: str) -> str:
    """The dictionary scope of a watch (`mpox` for monkeypox, `*` otherwise): clade terms only exist for mpox."""
    t = " ".join(names).lower()
    return "mpox" if ("monkeypox" in t or "mpox" in t) else "*"


# ------------------------------------------------------------------------------------------------ dates
def _valid(y: int, m: int | None = None, d: int | None = None) -> bool:
    try:
        if m is None:
            return 1000 <= y <= 2100
        date(y, m, d or 1)
        return True
    except ValueError:
        return False


def _one_date(s: str):
    """-> (value, kind) with kind 'exact' (already standard), 'format' (re-written) or None (cannot be read); ('', 'amb') when the order of day and month is ambiguous."""
    s = s.strip()
    m = re.fullmatch(r"(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?", s)
    if m:
        y, mo, d = int(m[1]), m[2] and int(m[2]), m[3] and int(m[3])
        return (s, "exact") if _valid(y, mo or None, d or None) else (None, None)
    m = re.fullmatch(r"(\d{4})-(\d{2}|XX)-(\d{2}|XX)", s, re.I)
    if m:
        y, mo = int(m[1]), m[2].upper()
        if mo == "XX":
            return (f"{y:04d}", "format") if _valid(y) else (None, None)
        if m[3].upper() == "XX":
            return (f"{y:04d}-{mo}", "format") if _valid(y, int(mo)) else (None, None)
    m = re.fullmatch(r"(\d{4})[/.](\d{1,2})(?:[/.](\d{1,2}))?", s)
    if m:
        y, mo, d = int(m[1]), int(m[2]), m[3] and int(m[3])
        if _valid(y, mo, d):
            return (f"{y:04d}-{mo:02d}" + (f"-{d:02d}" if d else ""), "format")
        return (None, None)
    m = re.fullmatch(r"(\d{4})(\d{2})(\d{2})", s)
    if m and _valid(int(m[1]), int(m[2]), int(m[3])):
        return f"{m[1]}-{m[2]}-{m[3]}", "format"
    m = re.fullmatch(r"(\d{1,2})[-\s/.]([A-Za-z]{3,9})[-\s/.,]+(\d{4})", s)
    if m and m[2][:3].lower() in MONTHS and _valid(int(m[3]), MONTHS[m[2][:3].lower()], int(m[1])):
        return f"{int(m[3]):04d}-{MONTHS[m[2][:3].lower()]:02d}-{int(m[1]):02d}", "format"
    m = re.fullmatch(r"([A-Za-z]{3,9})[-\s/.,]+(\d{1,2})[-\s/.,]+(\d{4})", s)                      # May 15, 2025
    if m and m[1][:3].lower() in MONTHS and _valid(int(m[3]), MONTHS[m[1][:3].lower()], int(m[2])):
        return f"{int(m[3]):04d}-{MONTHS[m[1][:3].lower()]:02d}-{int(m[2]):02d}", "format"
    m = re.fullmatch(r"([A-Za-z]{3,9})[-\s/.,]+(\d{4})", s)
    if m and m[1][:3].lower() in MONTHS:
        return f"{int(m[2]):04d}-{MONTHS[m[1][:3].lower()]:02d}", "format"
    m = re.fullmatch(r"(\d{1,2})[/.-](\d{1,2})[/.-](\d{4})", s)
    if m:
        a, b, y = int(m[1]), int(m[2]), int(m[3])
        if a > 12 and b <= 12 and _valid(y, b, a):
            return f"{y:04d}-{b:02d}-{a:02d}", "format"
        if b > 12 and a <= 12 and _valid(y, a, b):
            return f"{y:04d}-{a:02d}-{b:02d}", "format"
        if a == b and _valid(y, a, b):
            return f"{y:04d}-{a:02d}-{b:02d}", "format"
        return ("", "amb")
    return (None, None)


def standardize_date(raw: str) -> Std:
    s = raw.strip()
    if s.lower() in EMPTY:
        return Std("", "empty", "date.empty", "비어 있음을 뜻하는 표기")
    if "/" in s and re.fullmatch(r"\d{4}-\d{2}(-\d{2})?/\d{4}-\d{2}(-\d{2})?", s):                  # ISO interval a/b
        a, b = s.split("/")
        (va, ka), (vb, kb) = _one_date(a), _one_date(b)
        if not va or not vb:
            return Std(raw, "flagged", "date.range", "구간의 한쪽을 읽을 수 없음")
        if va[:7] == vb[:7]:
            return Std(va[:7], "normalized", "date.range->month", f"구간 {s}을(를) 월 단위로 (원본은 standardization 기록에 있음)")
        if va[:4] == vb[:4]:
            return Std(va[:4], "normalized", "date.range->year", f"구간 {s}을(를) 연 단위로 (월이 다름)")
        return Std(raw, "flagged", "date.range", "여러 해에 걸친 구간 — 하나의 날짜로 줄이지 않고 그대로 둠")
    v, kind = _one_date(s)
    if kind == "exact":
        return Std(s if s == raw else s, "same" if s == raw else "normalized", "" if s == raw else "date.trim")
    if kind == "format":
        return Std(v, "normalized", "date.format", f"{s} -> {v}")
    if kind == "amb":
        return Std(raw, "flagged", "date.ambiguous", "일/월 순서가 모호한 표기(예: 05/06/2024) — 그대로 둠")
    return Std(raw, "flagged", "date.unreadable", "날짜로 읽을 수 없는 표기 — 그대로 둠")


# ------------------------------------------------------------------------------------------------ term lookup
def standardize(field: str, raw, scope: str = "*") -> Std:
    """The standard form of one incoming value (see the module text for what the statuses mean)."""
    if raw is None:
        return Std("", "empty")
    s = str(raw)
    if field == "collection_date":
        return standardize_date(s)
    clean = " ".join(s.split())
    if clean.lower() in EMPTY:
        return Std("", "empty", f"{FILES.get(field, field)}.empty", "비어 있음을 뜻하는 표기") if clean else Std("", "empty")
    d = load(field)
    if field == "clade" and not any(k[0] == scope for k in d):                                       # no clade dictionary for this virus: leave it alone
        return Std(s, "nodict")
    probe = clean
    if field == "isolation_country" and ":" in clean:                                                  # INSDC style `Country: region`
        probe = clean.split(":", 1)[0].strip()
    hit = d.get((scope, key(probe))) or d.get(("*", key(probe)))
    if not hit:
        return Std(clean if clean != s else s, "unmapped", f"{FILES.get(field, field)}.unmapped", "사전에 없는 값 — 그대로 두고 검토 대상으로 표시")
    term = hit["term"]
    note = hit["note"]
    if hit["flag"]:
        return Std(term, "flagged", f"{FILES[field]}.flag:{hit['flag']}", note or hit["flag"])
    if probe != clean:
        return Std(term, "normalized", f"{FILES[field]}.split", f"'{clean}'에서 나라만 사용")
    if s == term:
        return Std(term, "same")
    if clean == term:
        return Std(term, "normalized", f"{FILES[field]}.trim")
    return Std(term, "mapped", f"{FILES[field]}.synonym", note)


# ------------------------------------------------------------------------------------------------ applying it
def apply_to_sources(per_source: dict[str, dict[str, dict]], scope: str = "*") -> tuple[list[dict], Counter]:
    """Standardise the fields of every record of every source in place. -> (changes, counts by (field, status)); `changes` lists what is not `same` (and keeps the raw value)."""
    changes: list[dict] = []
    counts: Counter = Counter()
    for src, recs in per_source.items():
        for acc, rec in recs.items():
            for field in FIELDS:
                v = rec.get(field)
                if v in (None, "") or not isinstance(v, str):
                    continue
                r = standardize(field, v, scope)
                counts[(field, r.status)] += 1
                if r.status in ("same", "nodict"):
                    continue
                if r.value != v:
                    rec[field] = r.value
                changes.append({"source": src, "accession": acc, "field": field, "raw": v, "standard": r.value, "status": r.status, "rule": r.rule, "note": r.note})
    return changes, counts


def store(con, watch_id: int, changes: list[dict]) -> None:
    """Replace the watch's rows of the table `standardization` with the changes of the latest check."""
    con.execute("DELETE FROM standardization WHERE watch_id=?", (watch_id,))
    con.executemany("INSERT OR REPLACE INTO standardization(watch_id,accession,source,field,raw,standard,status,rule,note,checked_at) VALUES(?,?,?,?,?,?,?,?,?,datetime('now'))",
                    [(watch_id, c["accession"], c["source"], c["field"], c["raw"], c["standard"], c["status"], c["rule"], c["note"]) for c in changes])


# ------------------------------------------------------------------------------------------------ what it would do to the data we hold
def preview(con, watch_id: int, scope: str = "*") -> dict:
    """Run the standardisation over the curated table as it is now (nothing is written): per field the distinct raw values with how many genomes carry them, the standard value
    and the status; plus the count per status. The same values come back through the next check as proposed changes."""
    out: dict = {}
    for field in FIELDS:
        rows = con.execute(f"SELECT {field} v, COUNT(*) n FROM genome WHERE watch_id=? AND IFNULL(genome_status,'')!='Removed' AND {field} IS NOT NULL AND {field}!='' GROUP BY {field} ORDER BY n DESC", (watch_id,)).fetchall()
        items, by = [], Counter()
        for r in rows:
            s = standardize(field, r["v"], scope)
            by[s.status] += r["n"]
            items.append({"raw": r["v"], "standard": s.value, "status": s.status, "rule": s.rule, "note": s.note, "n": r["n"]})
        out[field] = {"items": items, "by_status": dict(by), "distinct": len(items), "genomes": sum(r["n"] for r in rows),
                      "changes": sum(n for st, n in by.items() if st in ("mapped", "normalized", "empty")), "review": sum(n for st, n in by.items() if st in ("flagged", "unmapped")),
                      "distinct_after": len({i["standard"] for i in items})}
    return out
