"""Schema dry-run: have the sources (or our own database) changed shape in a way the code does not expect?

For each source a small live sample is fetched (read-only, no sequences kept), and our real adapter functions
(`normalize_ncbi`, `normalize_bvbrc`, `normalize_ppx`) are run on it through a recording wrapper, so the list of fields the
code *reads* is measured, not hand-maintained. That is compared with
  * what the source really offers now (new / vanished / retyped fields), and
  * the previous snapshot (stored in `schema_snapshot`), so a field that used to be filled and is now gone is caught.
The local database is checked against the schema the code expects (columns, adapter output columns, the analysis view).
The report says whether the source code needs an update, and where."""
import difflib
import json
from datetime import datetime

from . import config, curate, db, log, reconcile, sources
from .db import connect

_log = log.get("schema")
SAMPLE = 200

WHERE = {"ncbi": "sources.normalize_ncbi (genomecurator/sources.py)",
         "bvbrc": "sources.normalize_bvbrc (genomecurator/sources.py)",
         "pathoplexus": "sources.normalize_ppx / _PPX_FIELDS (genomecurator/sources.py)"}
LABEL = {"ncbi": "NCBI Virus", "bvbrc": "BV-BRC", "pathoplexus": "Pathoplexus"}


class _Track(dict):
    """dict that records every key path the code touches (even when the key is absent), recursively."""

    def __init__(self, d, seen, prefix=""):
        super().__init__(d)
        self._seen, self._p = seen, prefix

    def _wrap(self, k, v):
        path = self._p + k
        self._seen.add(path)
        if isinstance(v, dict):
            return _Track(v, self._seen, path + ".")
        if isinstance(v, list):
            return [_Track(x, self._seen, path + "[].") if isinstance(x, dict) else x for x in v]
        return v

    def get(self, k, default=None):
        return self._wrap(k, super().get(k)) if (k in self and super().get(k) is not None) else (self._seen.add(self._p + k) or default)

    def __getitem__(self, k):
        return self._wrap(k, super().__getitem__(k))


def _kind(v) -> str:
    return "null" if v is None else {bool: "bool", int: "int", float: "float", str: "str", list: "list", dict: "object"}.get(type(v), type(v).__name__)


def _walk(rec: dict, prefix: str, acc: dict):
    for k, v in rec.items():
        path = prefix + k
        a = acc.setdefault(path, {"n": 0, "nonnull": 0, "types": set(), "ex": None})
        a["n"] += 1
        a["types"].add(_kind(v))
        if v not in (None, "", [], {}):
            a["nonnull"] += 1
            if a["ex"] is None and not isinstance(v, (dict, list)):
                a["ex"] = str(v)[:60]
        if isinstance(v, dict):
            _walk(v, path + ".", acc)
        elif isinstance(v, list):
            for x in v[:5]:
                if isinstance(x, dict):
                    _walk(x, path + "[].", acc)


def profile(records: list[dict]) -> dict:
    acc: dict = {}
    for r in records:
        one: dict = {}
        _walk(r, "", one)   # a list of objects contributes several elements to one record: count the record once
        for p, a in one.items():
            t = acc.setdefault(p, {"nonnull": 0, "types": set(), "ex": None})
            t["nonnull"] += 1 if a["nonnull"] else 0
            t["types"] |= a["types"]
            t["ex"] = t["ex"] or a["ex"]
    n = max(len(records), 1)
    return {p: {"fill": round(a["nonnull"] / n, 3), "types": sorted(t for t in a["types"] if t != "null"), "ex": a["ex"]} for p, a in sorted(acc.items())}


# ----------------------------------------------------------------------------- live samples
def _ncbi_sample(c, taxon, n=SAMPLE):
    r = sources._ncbi_get(c, f"{config.NCBI_BASE}/taxon/{sources.quote(taxon)}/dataset_report", page_size=n)
    j = r.json()
    return j.get("reports", []), {"top": sorted(j), "status": r.status_code}


def _bvbrc_sample(c, taxon_id, n=SAMPLE):
    r = sources._get(c, f"{config.BVBRC_BASE}?eq(taxon_lineage_ids,{taxon_id})&limit({n})", headers={"Accept": "application/json"})
    r.raise_for_status()
    rows = r.json()
    return (rows if isinstance(rows, list) else []), {"top": ["<list>"] if isinstance(rows, list) else ["<" + type(rows).__name__ + ">"], "status": r.status_code}


def _ppx_sample(c, slug, n=SAMPLE):
    r = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details", params={**sources._PPX_FILTER, "limit": n})
    r.raise_for_status()
    j = r.json()
    return j.get("data", []), {"top": sorted(j), "status": r.status_code}


def _read_set(source: str, rows: list[dict]) -> tuple[set, list[dict]]:
    """Run the real adapter on the sample; -> (paths the code touched, the records it produced)."""
    seen: set = set()
    produced: list[dict] = []
    fn = {"ncbi": sources.normalize_ncbi, "bvbrc": sources.normalize_bvbrc, "pathoplexus": sources.normalize_ppx}[source]
    for row in rows:
        try:
            res = fn(_Track(row, seen))
        except Exception as e:   # the adapter cannot even digest the live shape: that is itself the finding
            raise RuntimeError(f"어댑터가 실제 응답을 처리하지 못했습니다: {type(e).__name__}: {e}") from e
        if res:
            produced.append(res[1])
    if source == "pathoplexus":
        seen |= set(sources._PPX_FIELDS) | set(sources._PPX_TAXON)   # also requested by name (fields=...)
    return seen, produced


# ----------------------------------------------------------------------------- comparison
def _compare(source: str, live: dict, read: set, prev: dict | None, shape: dict, prev_shape: list | None, mode: str = "full", n: int = SAMPLE) -> list[dict]:
    """mode: full = everything; virus = daily per-virus (shape and type only: three records of one virus cannot prove a field is gone);
    pooled = daily, all viruses together (only fields that were (almost) always filled and are now absent everywhere)."""
    out: list[dict] = []
    where = WHERE[source]

    def add(sev, kind, path, detail, hint=""):
        out.append({"source": source, "severity": sev, "kind": kind, "path": path, "detail": detail, "hint": hint})

    if mode != "pooled" and prev_shape is not None and shape["top"] != prev_shape:
        add("fail", "shape", ", ".join(shape["top"]), f"응답의 최상위 구조가 바뀌었습니다: {', '.join(prev_shape)} → {', '.join(shape['top'])}", f"{where}의 응답 파싱 확인")
    live_paths = set(live)
    for p in sorted(read):
        info = live.get(p)
        was = (prev or {}).get(p)
        was_fill = min(was.get("fill", 0), 1.0) if was else 0.0   # a baseline saved by older code could hold impossible values (> 100%)
        if info is None or info["fill"] == 0:
            if mode == "virus":
                continue
            if mode == "pooled" and not was_fill >= 0.9:
                continue
            if was and was_fill > 0 and was_fill * n < 5:
                continue   # a sparse field (expected < 5 hits in this sample): not turning up is likely chance, not evidence it is gone
            if was and was_fill > 0:
                close = difflib.get_close_matches(p.split(".")[-1], [q.split(".")[-1] for q in live_paths], n=2, cutoff=0.6)
                alt = [q for q in live_paths if q.split(".")[-1] in close and q not in read][:3]
                add("fail", "removed", p, f"코드가 읽는 필드인데 지난번엔 값이 {was_fill:.0%} 있었고 지금은 샘플 {n}건에서 모두 비어 있습니다.",
                    f"{where}에서 '{p}'를 읽는 부분을 확인하세요." + (f" 이름이 바뀐 것일 수 있습니다: {', '.join(alt)}" if alt else ""))
            elif prev is None:   # first run: nothing to compare with, so rare optional fields are only noted
                add("info", "unseen", p, f"코드가 읽는 필드인데 샘플 {SAMPLE}건에서 값이 없었습니다. 원래 드문 필드일 수 있어 기준으로 저장하고, 다음 점검부터 변화를 비교합니다.", "")
        elif mode != "pooled" and was and set(was.get("types", [])) and set(info["types"]) - set(was["types"]):
            add("warn", "type", p, f"값의 형태가 바뀌었습니다: {'/'.join(was['types'])} → {'/'.join(info['types'])}", f"{where}에서 '{p}' 처리 확인")
        elif mode == "pooled":
            continue
        elif was and was_fill >= 0.5 and info["fill"] < was_fill * 0.5:
            add("warn", "dropped", p, f"채워진 비율이 크게 줄었습니다: {was_fill:.0%} → {info['fill']:.0%}", "소스가 이 필드를 덜 채우기 시작했을 수 있습니다.")
    for p, info in ({} if mode != "full" else live).items():
        if p in read or any(r.startswith(p + ".") or r.startswith(p + "[].") for r in read):   # read, or a container whose children are read
            continue
        if info["fill"] == 0:
            continue
        was = (prev or {}).get(p)
        if prev is not None and was is None:
            add("info", "new", p, f"새로 나타난 필드(채움 {info['fill']:.0%}, 예: {info['ex']}). 코드는 아직 쓰지 않습니다.", f"쓰고 싶다면 {where}에 매핑을 추가하세요.")
        elif prev is None and info["fill"] >= 0.25:
            add("info", "unused", p, f"코드가 읽지 않는 필드(채움 {info['fill']:.0%}, 예: {info['ex']}).", "")
    return out


# ----------------------------------------------------------------------------- local schema
def local_checks(produced: dict[str, list[dict]]) -> list[dict]:
    out: list[dict] = []

    def add(sev, kind, path, detail, hint=""):
        out.append({"source": "local", "severity": sev, "kind": kind, "path": path, "detail": detail, "hint": hint})

    import sqlite3
    ref = sqlite3.connect(":memory:")
    ref.executescript(db.SCHEMA)
    db.migrate(ref)
    expect = {t: {r[1] for r in ref.execute(f"PRAGMA table_info({t})")} for (t,) in
              ref.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
    con = sqlite3.connect(config.DB_PATH)   # not connect(): that would silently re-apply migrations and hide a missing column
    try:
        actual = {t: {r[1] for r in con.execute(f"PRAGMA table_info({t})")} for (t,) in
                  con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
        for t, cols in expect.items():
            if t not in actual:
                add("fail", "table", t, f"코드가 기대하는 테이블 '{t}'가 로컬 DB에 없습니다.", "앱을 한 번 다시 시작하면 만들어집니다(마이그레이션).")
            elif cols - actual[t]:
                add("fail", "column", t, f"테이블 '{t}'에 열이 없습니다: {', '.join(sorted(cols - actual[t]))}", "마이그레이션(db.migrate)이 적용되지 않았습니다.")
        for t in sorted(set(actual) - set(expect)):
            add("info", "extra_table", t, f"코드가 모르는 테이블 '{t}'가 있습니다.", "")
        for t, cols in actual.items():
            if t in expect and cols - expect[t]:
                add("info", "extra_column", t, f"테이블 '{t}'에 코드가 모르는 열이 있습니다: {', '.join(sorted(cols - expect[t]))}", "")
        try:
            con.execute("SELECT * FROM v_genome_analysis LIMIT 0")
        except Exception as e:
            add("fail", "view", "v_genome_analysis", f"분석용 뷰가 열리지 않습니다: {e}", "앱을 다시 시작하면 뷰가 다시 만들어집니다.")
    finally:
        con.close()
    g = expect.get("genome", set())
    for name, cols in (("curate.GENOME_COLS", set(curate.GENOME_COLS)), ("config.TRACKED_FIELDS", set(config.TRACKED_FIELDS)),
                       ("reconcile.FIELDS", set(reconcile.FIELDS) - {"accession_version"})):
        gone = cols - g
        if gone:
            add("fail", "mapping", name, f"{name}에 있는 필드가 genome 테이블에 없습니다: {', '.join(sorted(gone))}", "db.SCHEMA/migrate 또는 해당 목록을 맞추세요.")
    for src, recs in produced.items():
        emitted = {k for r in recs for k in r}
        lost = emitted - set(curate.GENOME_COLS) - {"accession_version", "bvbrc_genome_id", "pathoplexus_accession"}
        if lost:
            add("fail", "dropped_on_apply", src, f"{LABEL[src]} 어댑터가 내는 필드가 genome 열에 없어 적용 때 버려집니다: {', '.join(sorted(lost))}", "curate.GENOME_COLS와 db 스키마에 열을 추가하세요.")
    return out


# ----------------------------------------------------------------------------- field -> column map ("which column does this field feed?")
_PROBES = ["§S§", "RefSeq", "COMPLETE"]       # the last two flip the few places where code compares against a known value
_TEMPLATE = {"ncbi": {"accession": "ZZ000001.1"}, "bvbrc": {}, "pathoplexus": {"insdcAccessionFull": "ZZ000001.1"}}   # the adapters drop rows without an accession


def _put(obj: dict, parts: list[str], val):
    """Create the nested path (lists of objects get one element) holding only `val`."""
    p = parts[0]
    many = p.endswith("[]")
    k = p[:-2] if many else p
    if len(parts) == 1:
        obj[k] = [val] if many else val
        return
    child = obj.setdefault(k, [{}] if many else {})
    _put(child[0] if many else child, parts[1:], val)


def field_map(source: str, read: set[str]) -> dict[str, list[str]]:
    """{path the code reads: [local columns it feeds]}, measured, not declared: the real adapter runs on a row that contains only
    this field (so "A, else B" fallbacks show up too) and the columns that change against an empty row are the ones it feeds."""
    fn = {"ncbi": sources.normalize_ncbi, "bvbrc": sources.normalize_bvbrc, "pathoplexus": sources.normalize_ppx}[source]
    fillers = [p for p in read if "[]" in p]   # siblings inside one list must exist, or the adapter's own x["name"] would raise

    def flat(res):
        if res is None:
            return {}
        if source == "bvbrc":
            accs, rec = res
            return {**rec, "accession": tuple(accs)}
        base, rec = res
        return {**rec, "accession": base}

    def row(path=None, val=None):
        r = json.loads(json.dumps(_TEMPLATE[source]))
        for f in fillers:
            _put(r, f.split("."), "~")
        if path:
            _put(r, path.split("."), val)
        return r

    try:
        base = flat(fn(row()))
    except Exception:
        base = {}
    out: dict[str, list[str]] = {}
    for path in sorted(read):
        cols: set[str] = set()
        for v in _PROBES:
            try:
                got = flat(fn(row(path, v)))
            except Exception:
                continue            # a container such as `host` (an object): its children carry the mapping
            cols |= {k for k in set(got) | set(base) if got.get(k) != base.get(k)}
        out[path] = sorted(cols)
    return out


# ----------------------------------------------------------------------------- the structure board
RANK = {"ok": 0, "unused": 0, "internal": 0, "new": 1, "unseen": 1, "warn": 2, "fail": 3}


def _field_status(path: str, read: bool, fill: float, fx: dict) -> tuple[str, str]:
    """-> (status, message) for one source field from this run's findings (own path first, then its parent containers)."""
    chain = [path] + [path.rsplit(sep, 1)[0] for sep in (".", "[]") if sep in path] + [path.split(".")[0], path.split("[]")[0]]
    best = None
    for q in chain:
        f = fx.get(q)
        if f and (best is None or RANK.get(f[0], 0) > RANK.get(best[0], 0)):
            best = f
    if best:
        return best
    if not read:
        return "unused", "코드가 읽지 않는 필드"
    return ("ok", "정상") if fill > 0 else ("unseen", "샘플에 값이 없었습니다")


def _groups(path: str) -> str:
    return (path.split(".")[0].replace("[]", "") if "." in path else "") or ""


def build_board(per: dict, findings: list[dict], local_missing: set[str], local_cols: list[str]) -> dict:
    """Everything the Schema tab needs to draw "source fields -> my genome columns" with problems highlighted."""
    KIND_SEV = {"removed": "fail", "type": "warn", "dropped": "warn", "new": "new", "unseen": "unseen"}
    board = {"sources": {}, "local": {}}
    feeders: dict[str, list[dict]] = {}
    for src, d in per.items():
        fx: dict[str, tuple[str, str]] = {}
        src_level = None
        for f in findings:
            if f["source"] != src:
                continue
            if f["kind"] in ("shape", "unreachable"):
                src_level = (f["severity"], f["detail"])
            sev = KIND_SEV.get(f["kind"])
            if sev and f["path"] != "-":
                old = fx.get(f["path"])
                if not old or RANK[sev] > RANK[old[0]]:
                    fx[f["path"]] = (sev if sev in ("fail", "warn") else sev, f["detail"])
        fields = []
        live = d["live"]
        paths = set(live) | set(d["read"])
        for path in sorted(paths):
            info = live.get(path, {"fill": 0.0, "types": [], "ex": None})
            if info["types"] in (["object"], ["list"]) and any(q.startswith(path + ".") or q.startswith(path + "[].") for q in paths):
                continue          # a container (host, virus.lineage ...): shown as a group heading, its children are the columns
            read = path in d["read"]
            st, msg = _field_status(path, read, info["fill"], fx)
            cols = d["fmap"].get(path, []) if read else []
            grp = _groups(path)
            fields.append({"path": path, "name": (path[len(grp) + 1:] if grp and path.startswith(grp + ".") else path), "group": grp, "read": read, "status": st, "msg": msg,
                           "fill": info["fill"], "types": info["types"], "ex": info["ex"], "cols": cols})
            for c in cols:
                feeders.setdefault(c, []).append({"source": src, "path": path, "status": st})
        fields.sort(key=lambda x: (x["group"], not x["read"], x["name"]))
        worst = max((RANK[x["status"]] for x in fields if x["read"]), default=0)
        board["sources"][src] = {"label": LABEL[src], "fields": fields, "n_read": sum(x["read"] for x in fields), "n_fields": len(fields),
                                 "status": ("fail" if src_level and src_level[0] == "fail" else "fail" if worst >= 3 else "warn" if (worst == 2 or src_level) else "ok"),
                                 "banner": src_level[1] if src_level else ""}
    order = {s: i for i, s in enumerate(config.SOURCE_PRIORITY)}
    cols = []
    for name in local_cols:
        fs = sorted(feeders.get(name, []), key=lambda x: (order.get(x["source"], 9), x["path"]))
        if name in local_missing:
            st, msg = "fail", "로컬 DB에 이 열이 없습니다 (마이그레이션 필요)"
        elif not fs:
            st, msg = "internal", "소스에서 오지 않는 열(앱이 채움)"
        else:
            bad = [x for x in fs if x["status"] == "fail"]
            part = [x for x in fs if x["status"] in ("warn",)]
            good = [x for x in fs if x["status"] in ("ok", "new", "unseen", "unused")]
            if bad and not good:
                st, msg = "fail", "이 열을 채우는 모든 소스 필드에 문제가 있습니다"
            elif bad or part:
                st, msg = "warn", "이 열을 채우는 소스 필드 중 일부에 문제가 있습니다" + (" (다른 소스가 대신 채울 수 있음)" if good else "")
            else:
                st, msg = "ok", "정상"
        cols.append({"name": name, "status": st, "msg": msg, "feeders": fs})
    board["local"] = {"table": "genome", "columns": cols}
    return board


def local_columns() -> tuple[list[str], set[str]]:
    """The genome table's columns as the code expects them (data columns first) and the ones missing from the real database."""
    import sqlite3
    ref = sqlite3.connect(":memory:")
    ref.executescript(db.SCHEMA)
    db.migrate(ref)
    expect = [r[1] for r in ref.execute("PRAGMA table_info(genome)")]
    con = sqlite3.connect(config.DB_PATH)
    try:
        have = {r[1] for r in con.execute("PRAGMA table_info(genome)")}
    finally:
        con.close()
    first = [c for c in curate.GENOME_COLS if c in expect] + [c for c in ("accession", "accession_version") if c in expect and c not in curate.GENOME_COLS]
    ordered = list(dict.fromkeys(first + expect))
    return ordered, {c for c in expect if c not in have}


# ----------------------------------------------------------------------------- run
def _sequence_probe(c, watch_id: int) -> list[dict]:
    """Tiny live probe of each sequence endpoint using the smallest local genome (a few KB at most)."""
    out = []
    con = connect()
    try:
        r = con.execute("SELECT accession, accession_version, genome_length, bvbrc_genome_id, pathoplexus_accession FROM genome "
                        "WHERE watch_id=? AND genome_length>0 ORDER BY genome_length LIMIT 1", (watch_id,)).fetchone()
        slug = con.execute("SELECT pathoplexus_organism FROM watchlist WHERE id=?", (watch_id,)).fetchone()[0]
    finally:
        con.close()
    if not r:
        return [{"source": s, "severity": "info", "kind": "sequence", "path": "-", "detail": "로컬 유전체가 없어 서열 엔드포인트 점검을 건너뜀", "hint": ""} for s in LABEL]

    def rec(src, ok, detail):
        out.append({"source": src, "severity": "ok" if ok else "fail", "kind": "sequence", "path": "서열 엔드포인트", "detail": detail,
                    "hint": "" if ok else f"{WHERE[src]} / 서열 내려받기 코드 확인"})
    try:
        got = sources.download_ncbi_fasta([f"{r['accession']}.{r['accession_version']}"])
        rec("ncbi", bool(got) and len(next(iter(got.values()))) == r["genome_length"], f"{r['accession']} {r['genome_length']} bp 내려받기 {'성공' if got else '실패'}")
    except Exception as e:
        rec("ncbi", False, f"내려받기 오류: {type(e).__name__}: {e}")
    if r["bvbrc_genome_id"]:
        try:
            rr = sources._get(c, f"https://www.bv-brc.org/api/genome_sequence/?eq(genome_id,{r['bvbrc_genome_id']})&select(genome_id,length,sequence)&limit(5)", headers={"Accept": "application/json"})
            rows = rr.json() if rr.status_code == 200 else []
            rec("bvbrc", bool(rows) and bool(rows[0].get("sequence")), f"HTTP {rr.status_code}, 서열 {len(rows[0].get('sequence', '')) if rows else 0} bp")
        except Exception as e:
            rec("bvbrc", False, f"오류: {type(e).__name__}: {e}")
    if slug and r["pathoplexus_accession"]:
        try:
            rr = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/unalignedNucleotideSequences", params={"accessionVersion": r["pathoplexus_accession"]})
            rec("pathoplexus", rr.status_code == 200 and rr.text.startswith(">"), f"HTTP {rr.status_code}, 응답 {len(rr.text)}바이트")
        except Exception as e:
            rec("pathoplexus", False, f"오류: {type(e).__name__}: {e}")
    return out


def run(watch_id: int | None = None, progress=None) -> dict:
    p = progress or (lambda *a: None)
    con = connect()
    try:
        w = (con.execute("SELECT * FROM watchlist WHERE id=?", (watch_id,)).fetchone() if watch_id
             else con.execute("SELECT * FROM watchlist ORDER BY (SELECT COUNT(*) FROM genome g WHERE g.watch_id=watchlist.id) DESC LIMIT 1").fetchone())
        snaps = {r[0]: json.loads(r[2]) for r in con.execute("SELECT source, taken, snapshot_json FROM schema_snapshot")}
        prev_taken = {r[0]: r[1] for r in con.execute("SELECT source, taken FROM schema_snapshot")}
    finally:
        con.close()
    if not w:
        raise ValueError("등록된 바이러스가 없어 샘플을 가져올 수 없습니다. 먼저 바이러스를 추가하세요.")
    findings: list[dict] = []
    produced: dict[str, list[dict]] = {}
    snapshots: dict[str, dict] = {}
    info: dict[str, dict] = {}
    per: dict[str, dict] = {}
    con = connect()
    try:
        allw = [dict(r) for r in con.execute("SELECT * FROM watchlist ORDER BY name")]
    finally:
        con.close()
    targets = [x for x in allw if watch_id is None or x["id"] == watch_id] or [dict(w)]
    ppx_slugs = sorted({x["pathoplexus_organism"] for x in targets if x["pathoplexus_organism"]}) or \
        [("mpox" if "mpox" in sources.PPX_KNOWN else sources.PPX_KNOWN[0])]
    with sources._client() as c:
        plans = [("ncbi", "NCBI Virus", [(x["name"], (lambda x=x: _ncbi_sample(c, x["taxon"]))) for x in targets]),
                 ("bvbrc", "BV-BRC", [(x["name"], (lambda x=x: _bvbrc_sample(c, x["taxon_id"]))) for x in targets]),
                 ("pathoplexus", "Pathoplexus", [(sl, (lambda sl=sl: _ppx_sample(c, sl))) for sl in ppx_slugs])]
        for i, (src, label, getters) in enumerate(plans):
            p(10 + i * 25, f"{label} 응답 구조 확인 중…")
            rows: list[dict] = []
            shape = None
            tags: list[str] = []
            err = None
            for tag, getter in getters:   # every watched virus contributes, so the baseline is not one virus' peculiarities
                try:
                    r, sh = getter()
                except Exception as e:
                    _log.error("schema check: %s (%s) failed: %s: %s", src, tag, type(e).__name__, e, exc_info=True)
                    err = f"{type(e).__name__}: {e}"
                    continue
                if r:
                    rows += r
                    shape = shape or sh
                    tags.append(tag)
            try:
                if not rows:
                    raise RuntimeError(err or "샘플 응답이 비어 있습니다")
                read, prod = _read_set(src, rows)
            except Exception as e:
                findings.append({"source": src, "severity": "fail", "kind": "unreachable", "path": "-", "detail": f"{label} 점검 실패: {type(e).__name__}: {e}",
                                 "hint": f"{WHERE[src]} 와 소스 상태 확인"})
                continue
            live = profile(rows)
            prev = snaps.get(src)
            # `read` depends on the data (an absent `host` object hides `host.organism_name`), so the paths the code read at the last good
            # baseline are added: a field whose parent vanished is then reported and mapped like any other.
            read_u = read | set(prev.get("read", [])) if prev else read
            findings += _compare(src, live, read_u, prev["fields"] if prev else None, shape, prev["shape"] if prev else None, n=len(rows))
            produced[src] = prod
            per[src] = {"live": live, "read": read_u, "fmap": field_map(src, read_u)}
            snapshots[src] = {"fields": live, "shape": shape["top"], "read": sorted(read)}
            info[src] = {"label": label, "sample": len(rows), "fields": len(live), "read": len(read), "shape": shape["top"],
                         "previous": prev_taken.get(src), "sample_of": ", ".join(tags)}
            _log.info("schema check %s: %d sample rows from %s, %d fields, code reads %d", src, len(rows), tags, len(live), len(read))
        _wid = w["id"]
        p(88, "서열 엔드포인트 확인 중…")
        findings += _sequence_probe(c, _wid)
    p(94, "로컬 DB 스키마 확인 중…")
    findings += local_checks(produced)
    lcols, lmissing = local_columns()
    board = build_board(per, findings, lmissing, lcols)
    sev = {"fail": 0, "warn": 1, "info": 2, "ok": 3}
    findings.sort(key=lambda f: (sev[f["severity"]], f["source"], f["path"]))
    n_fail = sum(f["severity"] == "fail" for f in findings)
    n_warn = sum(f["severity"] == "warn" for f in findings)
    verdict = ("update" if n_fail else ("review" if n_warn else "ok"))
    report = {"checked_at": datetime.now().isoformat(timespec="seconds"), "watch": (w["name"] if len(targets) == 1 else f"{len(targets)}종 합산"), "verdict": verdict, "fail": n_fail, "warn": n_warn,
              "info": sum(f["severity"] == "info" for f in findings), "sources": info, "findings": findings, "board": board}
    dirty = {f["source"] for f in findings if f["severity"] in ("fail", "warn")}
    kept = sorted(src for src in snapshots if src in snaps and src in dirty)
    report["baseline_kept"] = kept          # sources whose old baseline was NOT replaced because of findings
    report["_snapshots"] = snapshots        # what "accept as new baseline" will store
    con = connect()
    try:
        now = report["checked_at"]
        for src, snap in snapshots.items():
            if src in kept:   # a drifted state must not become the baseline, or the next run would call it normal
                continue
            con.execute("INSERT OR REPLACE INTO schema_snapshot(source, taken, snapshot_json) VALUES(?,?,?)", (src, now, json.dumps(snap, ensure_ascii=False)))
        con.execute("INSERT INTO schema_report(checked_at, report_json) VALUES(?,?)", (now, json.dumps(report, ensure_ascii=False)))
        con.execute("DELETE FROM schema_report WHERE id NOT IN (SELECT id FROM schema_report ORDER BY id DESC LIMIT 20)")
        db.audit(con, "schema_check", f"{w['name']}: verdict {verdict} (fail {n_fail}, warn {n_warn})" + (f"; baseline kept for {', '.join(kept)}" if kept else ""))
        con.commit()
    finally:
        con.close()
    _log.info("schema check done: verdict=%s fail=%d warn=%d", verdict, n_fail, n_warn)
    p(100, "스키마 점검 완료")
    return report


def accept() -> list[str]:
    """After the code was updated (or the change judged harmless): make the latest live profile the new baseline."""
    con = connect()
    try:
        r = next(((i, rj) for i, rj in con.execute("SELECT id, report_json FROM schema_report ORDER BY id DESC LIMIT 40")
                  if json.loads(rj).get("mode", "full") != "daily"), None)
        if not r:
            raise ValueError("먼저 스키마 점검을 실행하세요.")
        rep = json.loads(r[1])
        for src, snap in (rep.get("_snapshots") or {}).items():
            con.execute("INSERT OR REPLACE INTO schema_snapshot(source, taken, snapshot_json) VALUES(?,?,?)", (src, rep["checked_at"], json.dumps(snap, ensure_ascii=False)))
        rep["baseline_kept"] = []
        con.execute("UPDATE schema_report SET report_json=? WHERE id=?", (json.dumps(rep, ensure_ascii=False), r[0]))
        db.audit(con, "schema_accept", f"새 기준으로 저장: {', '.join(rep.get('_snapshots', {}))}")
        con.commit()
        _log.info("schema baseline accepted for %s", list(rep.get("_snapshots", {})))
        return list(rep.get("_snapshots", {}))
    finally:
        con.close()


def run_light(progress=None, sample: int | None = None) -> dict:
    """Daily check: every watched virus, `sample` records (default 3) per source. Cheap and strict: it reports only what
    three records can prove - a changed response shape, a field that is always filled and is now gone, a changed value type,
    or a source that cannot be read. Never overwrites the baseline. Without a baseline it runs the full check first."""
    p = progress or (lambda *a: None)
    n = sample or config.SCHEMA_CHECK_SAMPLE
    con = connect()
    try:
        watches = [dict(r) for r in con.execute("SELECT * FROM watchlist ORDER BY name")]
        snaps = {r[0]: json.loads(r[1]) for r in con.execute("SELECT source, snapshot_json FROM schema_snapshot")}
    finally:
        con.close()
    if not watches:
        raise ValueError("등록된 바이러스가 없습니다.")
    if not all(s in snaps for s in LABEL):
        _log.info("daily schema check: no baseline yet, running the full check first")
        p(5, "기준이 없어 전체 점검을 먼저 실행합니다…")
        rep = run(None, progress=lambda pc, m: p(5 + pc * 0.9, m))
        rep["mode"] = "full-baseline"
        return rep
    findings: list[dict] = []
    matrix: dict[str, dict[str, dict]] = {}
    pooled: dict[str, dict] = {s: {"rows": [], "read": set(), "shape": None} for s in LABEL}
    with sources._client() as c:
        for i, w in enumerate(watches):
            p(5 + 90 * i / len(watches), f"{w['name']}: 소스별 {n}건 확인 중…")
            row = matrix.setdefault(w["name"], {})
            jobs = [("ncbi", lambda w=w: _ncbi_sample(c, w["taxon"], n)), ("bvbrc", lambda w=w: _bvbrc_sample(c, w["taxon_id"], n))]
            if w["pathoplexus_organism"]:
                jobs.append(("pathoplexus", lambda w=w: _ppx_sample(c, w["pathoplexus_organism"], n)))
            for src, getter in jobs:
                with log.context(watch=w["name"]):
                    try:
                        rows, shape = getter()
                        if not rows:
                            row[src] = {"status": "empty", "n": 0}   # e.g. a virus BV-BRC has no genomes for: not a schema problem
                            continue
                        read, _ = _read_set(src, rows)
                        fs = _compare(src, profile(rows), read, snaps[src]["fields"], shape, snaps[src]["shape"], mode="virus", n=len(rows))
                        pooled[src]["rows"] += rows
                        pooled[src]["read"] |= read | set(snaps[src].get("read", []))
                        pooled[src]["shape"] = shape
                    except Exception as e:
                        _log.error("daily schema check %s/%s failed: %s: %s", w["name"], src, type(e).__name__, e, exc_info=True)
                        fs = [{"source": src, "severity": "warn", "kind": "unreachable", "path": "-",
                               "detail": f"점검 실패(일시적인 장애일 수 있음): {type(e).__name__}: {e}", "hint": f"{WHERE[src]} 와 소스 상태 확인"}]
                        rows = []
                    for f in fs:
                        f["virus"] = w["name"]
                    findings += fs
                    worst = "ok" if not fs else min((f["severity"] for f in fs), key=lambda x: {"fail": 0, "warn": 1, "info": 2}.get(x, 3))
                    row[src] = {"status": worst, "n": len(rows)}
            _log.info("daily schema check %r: %s", w["name"], {k: v["status"] for k, v in row.items()})
    for src, pl in pooled.items():   # "a field is gone" is judged over all viruses together, never from one virus' 3 records
        if pl["rows"]:
            for f in _compare(src, profile(pl["rows"]), pl["read"], snaps[src]["fields"], pl["shape"], snaps[src]["shape"], mode="pooled", n=len(pl["rows"])):
                f["virus"] = "(전체 합산)"
                findings.append(f)
    n_fail = sum(f["severity"] == "fail" for f in findings)
    n_warn = sum(f["severity"] == "warn" for f in findings)
    report = {"mode": "daily", "checked_at": datetime.now().isoformat(timespec="seconds"), "sample": n, "verdict": "update" if n_fail else ("review" if n_warn else "ok"),
              "fail": n_fail, "warn": n_warn, "info": 0, "matrix": matrix, "findings": findings}
    con = connect()
    try:
        con.execute("INSERT INTO schema_report(checked_at, report_json) VALUES(?,?)", (report["checked_at"], json.dumps(report, ensure_ascii=False)))
        con.execute("DELETE FROM schema_report WHERE id NOT IN (SELECT id FROM schema_report ORDER BY id DESC LIMIT 40)")
        db.audit(con, "schema_daily", f"{len(watches)}종 × 소스당 {n}건: verdict {report['verdict']} (fail {n_fail}, warn {n_warn})")
        con.commit()
    finally:
        con.close()
    if n_fail or n_warn:
        _log.error("daily schema check: verdict=%s fail=%d warn=%d", report["verdict"], n_fail, n_warn)
    else:
        _log.info("daily schema check: all %d viruses ok", len(watches))
    p(100, "매일 점검 완료")
    return report


def _latest_of(modes: tuple[str, ...]) -> dict | None:
    con = connect()
    try:
        for (rj,) in con.execute("SELECT report_json FROM schema_report ORDER BY id DESC LIMIT 40"):
            rep = json.loads(rj)
            if rep.get("mode", "full") in modes:
                rep.pop("_snapshots", None)
                return rep
        return None
    finally:
        con.close()


def latest_daily() -> dict | None:
    return _latest_of(("daily",))


def alert() -> dict | None:
    """The newest report of either kind, if it says something needs attention (for the badge / banner)."""
    rep = _latest_of(("full", "daily", "full-baseline"))
    return rep if rep and rep.get("verdict") in ("update", "review") else None


def latest() -> dict | None:
    con = connect()
    try:
        for (rj,) in con.execute("SELECT report_json FROM schema_report ORDER BY id DESC LIMIT 40"):
            rep = json.loads(rj)
            if rep.get("mode", "full") != "daily":
                rep.pop("_snapshots", None)
                return rep
        return None
    finally:
        con.close()
