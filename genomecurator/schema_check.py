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
        _walk(r, "", acc)
    n = max(len(records), 1)
    return {p: {"fill": round(a["nonnull"] / n, 3), "types": sorted(t for t in a["types"] if t != "null"), "ex": a["ex"]} for p, a in sorted(acc.items())}


# ----------------------------------------------------------------------------- live samples
def _ncbi_sample(c, taxon):
    r = sources._ncbi_get(c, f"{config.NCBI_BASE}/taxon/{sources.quote(taxon)}/dataset_report", page_size=SAMPLE)
    j = r.json()
    return j.get("reports", []), {"top": sorted(j), "status": r.status_code}


def _bvbrc_sample(c, taxon_id):
    r = sources._get(c, f"{config.BVBRC_BASE}?eq(taxon_lineage_ids,{taxon_id})&limit({SAMPLE})", headers={"Accept": "application/json"})
    r.raise_for_status()
    rows = r.json()
    return (rows if isinstance(rows, list) else []), {"top": ["<list>"] if isinstance(rows, list) else ["<" + type(rows).__name__ + ">"], "status": r.status_code}


def _ppx_sample(c, slug):
    r = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details", params={**sources._PPX_FILTER, "limit": SAMPLE})
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
def _compare(source: str, live: dict, read: set, prev: dict | None, shape: dict, prev_shape: list | None) -> list[dict]:
    out: list[dict] = []
    where = WHERE[source]

    def add(sev, kind, path, detail, hint=""):
        out.append({"source": source, "severity": sev, "kind": kind, "path": path, "detail": detail, "hint": hint})

    if prev_shape is not None and shape["top"] != prev_shape:
        add("fail", "shape", ", ".join(shape["top"]), f"응답의 최상위 구조가 바뀌었습니다: {', '.join(prev_shape)} → {', '.join(shape['top'])}", f"{where}의 응답 파싱 확인")
    live_paths = set(live)
    for p in sorted(read):
        info = live.get(p)
        was = (prev or {}).get(p)
        if info is None or info["fill"] == 0:
            if was and was.get("fill", 0) > 0:
                close = difflib.get_close_matches(p.split(".")[-1], [q.split(".")[-1] for q in live_paths], n=2, cutoff=0.6)
                alt = [q for q in live_paths if q.split(".")[-1] in close and q not in read][:3]
                add("fail", "removed", p, f"코드가 읽는 필드인데 지난번엔 값이 {was['fill']:.0%} 있었고 지금은 샘플 {SAMPLE}건에서 모두 비어 있습니다.",
                    f"{where}에서 '{p}'를 읽는 부분을 확인하세요." + (f" 이름이 바뀐 것일 수 있습니다: {', '.join(alt)}" if alt else ""))
            elif prev is None:   # first run: nothing to compare with, so rare optional fields are only noted
                add("info", "unseen", p, f"코드가 읽는 필드인데 샘플 {SAMPLE}건에서 값이 없었습니다. 원래 드문 필드일 수 있어 기준으로 저장하고, 다음 점검부터 변화를 비교합니다.", "")
        elif was and set(was.get("types", [])) and set(info["types"]) - set(was["types"]):
            add("warn", "type", p, f"값의 형태가 바뀌었습니다: {'/'.join(was['types'])} → {'/'.join(info['types'])}", f"{where}에서 '{p}' 처리 확인")
        elif was and was.get("fill", 0) >= 0.5 and info["fill"] < was["fill"] * 0.5:
            add("warn", "dropped", p, f"채워진 비율이 크게 줄었습니다: {was['fill']:.0%} → {info['fill']:.0%}", "소스가 이 필드를 덜 채우기 시작했을 수 있습니다.")
    for p, info in live.items():
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
    slug = w["pathoplexus_organism"] or ("mpox" if "mpox" in sources.PPX_KNOWN else sources.PPX_KNOWN[0])
    findings: list[dict] = []
    produced: dict[str, list[dict]] = {}
    snapshots: dict[str, dict] = {}
    info: dict[str, dict] = {}
    with sources._client() as c:
        for i, (src, label, getter) in enumerate((
                ("ncbi", "NCBI Virus", lambda: _ncbi_sample(c, w["taxon"])),
                ("bvbrc", "BV-BRC", lambda: _bvbrc_sample(c, w["taxon_id"])),
                ("pathoplexus", "Pathoplexus", lambda: _ppx_sample(c, slug)))):
            p(10 + i * 25, f"{label} 응답 구조 확인 중…")
            try:
                rows, shape = getter()
                if not rows:
                    raise RuntimeError("샘플 응답이 비어 있습니다")
                read, prod = _read_set(src, rows)
            except Exception as e:
                _log.error("schema check: %s failed: %s: %s", src, type(e).__name__, e, exc_info=True)
                findings.append({"source": src, "severity": "fail", "kind": "unreachable", "path": "-", "detail": f"{label} 점검 실패: {type(e).__name__}: {e}",
                                 "hint": f"{WHERE[src]} 와 소스 상태 확인"})
                continue
            live = profile(rows)
            prev = snaps.get(src)
            findings += _compare(src, live, read, prev["fields"] if prev else None, shape, prev["shape"] if prev else None)
            produced[src] = prod
            snapshots[src] = {"fields": live, "shape": shape["top"], "read": sorted(read)}
            info[src] = {"label": label, "sample": len(rows), "fields": len(live), "read": len(read), "shape": shape["top"],
                         "previous": prev_taken.get(src), "sample_of": (w["taxon"] if src != "pathoplexus" else slug)}
            _log.info("schema check %s: %d sample rows, %d fields, code reads %d", src, len(rows), len(live), len(read))
        p(88, "서열 엔드포인트 확인 중…")
        findings += _sequence_probe(c, w["id"])
    p(94, "로컬 DB 스키마 확인 중…")
    findings += local_checks(produced)
    sev = {"fail": 0, "warn": 1, "info": 2, "ok": 3}
    findings.sort(key=lambda f: (sev[f["severity"]], f["source"], f["path"]))
    n_fail = sum(f["severity"] == "fail" for f in findings)
    n_warn = sum(f["severity"] == "warn" for f in findings)
    verdict = ("update" if n_fail else ("review" if n_warn else "ok"))
    report = {"checked_at": datetime.now().isoformat(timespec="seconds"), "watch": w["name"], "verdict": verdict, "fail": n_fail, "warn": n_warn,
              "info": sum(f["severity"] == "info" for f in findings), "sources": info, "findings": findings}
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
        r = con.execute("SELECT id, report_json FROM schema_report ORDER BY id DESC LIMIT 1").fetchone()
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


def latest() -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT report_json FROM schema_report ORDER BY id DESC LIMIT 1").fetchone()
        if not r:
            return None
        rep = json.loads(r[0])
        rep.pop("_snapshots", None)
        return rep
    finally:
        con.close()
