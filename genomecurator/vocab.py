"""The app's own dictionary of terms and the standardisation of incoming values against it.

Sources spell the same thing in many ways (`Homo sapiens`, `Homosapiens`, `H. sapiens`, `Patent`; `USA` / `United States`; `Zaire` / `Democratic Republic of the Congo`; collection dates as
`2025-05-01/2025-05-28`, `15-May-2025`, ...). Compared as they come, such differences look like conflicts and end up as different values in the curated table. Here every incoming value of
a standardised field is looked up in **our dictionary** and replaced by the standard term, **before** the sources are compared:

    host_name, isolation_country, collection_date, clade, sequencing_platform (model-level standard names; the technology goes into the derived column `platform`), isolation_source (specimen type)

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
FIELDS = ("host_name", "isolation_country", "collection_date", "clade", "sequencing_platform", "isolation_source")
FILES = {"host_name": "host", "isolation_country": "country", "clade": "clade", "collection_date": "date", "sequencing_platform": "platform", "isolation_source": "specimen"}
COLUMNS = {"host_name": ["term", "synonyms", "taxon_id", "flag", "note"], "isolation_country": ["term", "synonyms", "iso3", "note"], "clade": ["scope", "term", "synonyms", "flag", "note"],
           "collection_date": ["raw", "standard", "note"], "sequencing_platform": ["term", "platform", "synonyms", "flag", "note"], "isolation_source": ["term", "synonyms", "flag", "note"]}
MULTI = {"sequencing_platform"}        # several values in one field (`Illumina; Oxford Nanopore`): every part is looked up, the parts are joined again
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
        e = {"term": term, "flag": (r.get("flag") or "").strip(), "note": (r.get("note") or "").strip(), "iso3": (r.get("iso3") or "").strip(), "taxon_id": (r.get("taxon_id") or "").strip(), "platform": (r.get("platform") or "").strip(), "scope": scope}
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


def date_overrides() -> dict[str, dict]:
    """{raw date string: {standard, note}} decided by a person (`<data>/vocab/date.tsv`): for the writings the rules leave alone, e.g. `05/06/2024`."""
    out: dict = {}
    for _d, r in _date_rows():
        out[r["raw"].strip()] = {"standard": r["standard"].strip(), "note": (r.get("note") or "").strip()}
    return out


def _date_rows():
    for d in vocab_dirs():
        f = d / "date.tsv"
        if f.exists():
            for r in csv.DictReader((l for l in open(f, encoding="utf-8") if not l.startswith("#")), delimiter="\t"):
                if (r.get("raw") or "").strip() and (r.get("standard") or "").strip():
                    yield d, r


def standardize_date(raw: str) -> Std:
    s = raw.strip()
    ov = date_overrides().get(s)
    if ov:
        return Std(ov["standard"], "mapped", "date.decided", ov["note"] or "사용자가 정한 값")
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
    if field in MULTI:
        return _standardize_multi(field, s, clean, d)
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


def platform_of(value) -> str:
    """The technology / manufacturer of a standardised `sequencing_platform` value (`Illumina MiSeq; Oxford Nanopore MinION` -> `Illumina; Oxford Nanopore`); empty when a part is not in the dictionary."""
    d = load("sequencing_platform")
    out = []
    for x in (p.strip() for p in str(value or "").split(";")):
        h = d.get(("*", key(x))) if x else None
        if not h or not h.get("platform"):
            return ""
        out.append(h["platform"])
    return "; ".join(dict.fromkeys(out))


def _standardize_multi(field: str, s: str, clean: str, d: dict) -> Std:
    """A field that may hold several values (separated by `;`, `,` or `/`): each part is looked up; the standard value is the distinct standard terms in order of first appearance joined by `; `.
    One unknown part makes the whole value `unmapped` (kept as it is), one flagged part makes it `flagged`."""
    parts = [x.strip() for x in re.split(r"[;,/+]", clean) if x.strip()]
    hits = []
    whole = d.get(("*", key(clean)))
    if whole:
        hits = [whole]
    else:
        for x in parts:
            h = d.get(("*", key(x)))
            if not h:
                return Std(clean if clean != s else s, "unmapped", f"{FILES[field]}.unmapped", f"사전에 없는 부분: '{x}'" if len(parts) > 1 else "사전에 없는 값 — 그대로 두고 검토 대상으로 표시")
            hits.append(h)
    terms = list(dict.fromkeys(h["term"] for h in hits))
    val = "; ".join(terms)
    bad = [h for h in hits if h["flag"]]
    if bad:
        return Std(val, "flagged", f"{FILES[field]}.flag:{bad[0]['flag']}", bad[0]["note"] or bad[0]["flag"])
    return Std(val, "same" if val == s else "mapped", f"{FILES[field]}.synonym" if val != s else "")


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
                if field == "sequencing_platform" and r.status not in ("unmapped", "nodict"):
                    rec["platform"] = platform_of(r.value)                                 # derived column: the technology, while the model stays in sequencing_platform
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


# ------------------------------------------------------------------------------------------------ editing the user's dictionary (<data>/vocab/*.tsv)
def user_file(field: str) -> Path:
    return config.DATA_DIR / "vocab" / f"{FILES[field]}.tsv"


def user_entries(field: str) -> list[dict]:
    f = user_file(field)
    if not f.exists():
        return []
    return [dict(r) for r in csv.DictReader((l for l in open(f, encoding="utf-8") if not l.startswith("#")), delimiter="\t")]


def _write_user(field: str, rows: list[dict]) -> None:
    f = user_file(field)
    f.parent.mkdir(parents=True, exist_ok=True)
    cols = COLUMNS[field]
    tmp = f.with_suffix(".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        fh.write("# 사용자 사전 — 이 앱의 용어 탭에서 더한 항목입니다. 같은 표준 용어는 앱에 들어 있는 사전을 덮어씁니다. 손으로 고쳐도 됩니다.\n")
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(f)
    _CACHE.pop(field, None)


def _owner(field: str, k: str, scope: str):
    d = load(field)
    return d.get((scope, k)) or d.get(("*", k))


def add_term(field: str, term: str, synonyms=(), flag: str = "", note: str = "", scope: str = "*", extra: dict | None = None) -> dict:
    """Add a standard term (and synonyms) to the user's dictionary, or update it. Refuses a writing that already belongs to another term, so one writing never means two things."""
    if field not in COLUMNS or field == "collection_date":
        raise ValueError("이 필드는 용어 사전이 아닙니다.")
    term = " ".join(str(term).split())
    if not term:
        raise ValueError("표준 용어가 비었습니다.")
    syn = [" ".join(x.split()) for x in synonyms if x and x.strip()]
    for k in [term] + syn:
        o = _owner(field, key(k), scope)
        if o and o["term"] != term:
            raise ValueError(f"'{k}'은(는) 이미 '{o['term']}'의 표기입니다. 먼저 그쪽에서 빼거나, 그 용어의 동의어로 두세요.")
    rows = user_entries(field)
    cur = next((r for r in rows if r["term"] == term and (r.get("scope") or "*") == scope), None)
    if cur is None:
        cur = {"term": term, "scope": scope}
        rows.append(cur)
        # a term the app's dictionary already has: the user's row replaces it for this term's own spellings, so carry those over (else its synonyms would keep the old flag)
        cur["synonyms"] = "|".join(kk for (sc, kk), e in load(field).items() if e["term"] == term and e["scope"] == scope and kk != key(term))
    have = [x for x in (cur.get("synonyms") or "").split("|") if x]
    cur["synonyms"] = "|".join(dict.fromkeys(have + syn))
    cur["flag"], cur["note"] = flag.strip(), note.strip()
    for k, v in (extra or {}).items():
        cur[k] = v
    _write_user(field, rows)
    return {"term": term, "synonyms": cur["synonyms"].split("|") if cur["synonyms"] else [], "flag": cur["flag"]}


def add_synonym(field: str, raw: str, term: str, scope: str = "*") -> dict:
    """The writing `raw` is another way of writing the standard term `term` (which must exist)."""
    ts = {e["term"]: e for k, e in load(field).items() if k[0] in (scope, "*")}
    if term not in ts:
        raise ValueError(f"'{term}'은(는) 사전에 없는 표준 용어입니다. 먼저 새 용어로 추가하세요.")
    e = ts[term]
    return add_term(field, term, [raw], e.get("flag", ""), e.get("note", ""), e.get("scope", "*"), {k: e.get(k, "") for k in ("taxon_id", "iso3") if e.get(k)})


def set_date(raw: str, standard: str, note: str = "") -> dict:
    """Decide the standard date of one raw writing (e.g. an ambiguous `05/06/2024`): it must itself be a valid YYYY / YYYY-MM / YYYY-MM-DD."""
    raw, standard = raw.strip(), standard.strip()
    v, kind = _one_date(standard)
    if kind != "exact":
        raise ValueError("표준 날짜는 YYYY, YYYY-MM, YYYY-MM-DD 형식의 올바른 날짜여야 합니다.")
    if not raw:
        raise ValueError("원본 표기가 비었습니다.")
    rows = [r for r in user_entries("collection_date") if r["raw"].strip() != raw]
    rows.append({"raw": raw, "standard": standard, "note": note.strip()})
    _write_user("collection_date", rows)
    return {"raw": raw, "standard": standard}


def remove_entry(field: str, term: str, scope: str = "*") -> bool:
    """Remove an entry the user added (the app's own dictionary is not touched)."""
    key_col = "raw" if field == "collection_date" else "term"
    rows = user_entries(field)
    keep = [r for r in rows if not (r.get(key_col) == term and (field == "collection_date" or (r.get("scope") or "*") == scope))]
    if len(keep) == len(rows):
        return False
    _write_user(field, keep)
    return True


def suggest(field: str, raw: str, scope: str = "*", n: int = 3) -> list[dict]:
    """The standard terms a value is most likely a writing of (closest by spelling): [{term, score}] best first, only reasonably close ones."""
    import difflib
    k = key(raw)
    best: dict[str, float] = {}
    for (sc, kk), e in load(field).items():
        if sc not in (scope, "*") or e.get("flag"):
            continue
        r = difflib.SequenceMatcher(None, k, kk).ratio()
        if r > best.get(e["term"], 0):
            best[e["term"]] = r
    return [{"term": t, "score": round(r, 2)} for t, r in sorted(best.items(), key=lambda x: -x[1])[:n] if r >= 0.72]


def review(con, watch_id: int, scope: str = "*") -> dict:
    """What a person has to decide: per field the values the dictionary does not know (`unmapped`, with the closest terms) and the ones it keeps but marks (`flagged`), with how many genomes carry them."""
    pv = preview(con, watch_id, scope)
    out: dict = {}
    for f, d in pv.items():
        un = [{"raw": i["raw"], "n": i["n"], "suggestions": suggest(f, i["raw"], scope) if f != "collection_date" else []} for i in d["items"] if i["status"] == "unmapped"]
        fl = [{"raw": i["raw"], "standard": i["standard"], "n": i["n"], "rule": i["rule"], "note": i["note"]} for i in d["items"] if i["status"] == "flagged"]
        out[f] = {"unmapped": un, "flagged": fl, "terms": sorted(t for t in {e["term"] for k, e in load(f).items() if k[0] in (scope, "*")}) if f != "collection_date" else []}
    return out


def search_terms(field: str, q: str = "", scope: str = "*", limit: int = 80) -> list[dict]:
    """The dictionary of a field for browsing: [{term, synonyms, flag, note, user}] filtered by a text; `user` is true for entries from the user's own file."""
    if field == "collection_date":
        return [{"term": r["raw"], "synonyms": [r["standard"]], "flag": "", "note": r.get("note", ""), "user": True} for r in user_entries(field)]
    users = {(r["term"], (r.get("scope") or "*")) for r in user_entries(field)}
    by: dict = {}
    for (sc, kk), e in load(field).items():
        if sc not in (scope, "*"):
            continue
        b = by.setdefault((e["term"], e["scope"]), {"term": e["term"], "synonyms": set(), "flag": e["flag"], "note": e["note"], "user": (e["term"], e["scope"]) in users, "keys": set()})
        b["keys"].add(kk)
    rows = {}
    for src in (Path(__file__).parent / "resources" / "vocab", config.DATA_DIR / "vocab"):
        f = src / f"{FILES[field]}.tsv"
        if f.exists():
            for r in csv.DictReader((l for l in open(f, encoding="utf-8") if not l.startswith("#")), delimiter="\t"):
                rows[((r.get("term") or "").strip(), (r.get("scope") or "*").strip() or "*")] = [x for x in (r.get("synonyms") or "").split("|") if x.strip()]
    qk = key(q)
    out = []
    for (t, sc), b in sorted(by.items()):
        syn = rows.get((t, sc), [])
        if qk and qk not in key(t) and not any(qk in key(x) for x in syn):
            continue
        out.append({"term": t, "synonyms": syn, "flag": b["flag"], "note": b["note"], "user": b["user"], "scope": sc})
    return out[:limit]
