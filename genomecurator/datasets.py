"""Ad hoc datasets for analysis: a saved definition - which genomes belong - that can be previewed and exported again and again.

A dataset is rules over the curated genomes (and, if asked, the private store) plus manual additions and removals:
* rules (all must hold; an empty list means "no restriction"): clade, lineage, outbreak, region, country, host, completeness, year range, length range,
  QC statuses to leave out, and whether a sequence file is required. They are judged on the same row the export writes (so *clade* means the
  curated clade, else Nextclade's), which keeps "what the preview counts" and "what the export contains" the same thing;
* manual include: genomes that belong whatever the rules say; manual exclude: genomes that never belong to this dataset;
* the global *exclude from dataset* marks (exclusions.py) still apply to every dataset.
A dataset is only a definition - nothing is copied or deleted. Each export is a snapshot: its manifest records the rules and a hash of the exact
member ids, and the export history keeps them, so a result can always be traced to the dataset as it was then."""
import hashlib
import json
from collections import Counter
from datetime import datetime

from .db import audit, connect

RULE_LISTS = ("clade", "lineage", "outbreak", "region", "country", "host", "completeness", "qc_not")
RULE_NUMS = ("year_min", "year_max", "min_length", "max_length")
BUILDS = ("mpxv", "hmpxv1", "clade-i")        # which nextstrain/mpox default config the exported overlay suggests
PRESETS = {
    "orthopoxvirus monkeypox": [
        {"name": "clade Ib", "build": "clade-i", "rules": {"clade": ["Ib"]}},
        {"name": "clade Ia", "build": "clade-i", "rules": {"clade": ["Ia"]}},
        {"name": "clade I (Ia + Ib)", "build": "clade-i", "rules": {"clade": ["Ia", "Ib"]}},
        {"name": "clade IIb", "build": "hmpxv1", "rules": {"clade": ["IIb"]}},
        {"name": "all clades, complete genomes", "build": "mpxv", "rules": {"completeness": ["Complete", "Near complete"]}},
    ]}


def normalize(rules: dict | None) -> dict:
    """Clean a rules dict: known keys only, lists of strings, numbers as ints (or None)."""
    r = rules or {}
    out = {k: [str(x).strip() for x in (r.get(k) or []) if str(x).strip()] for k in RULE_LISTS}
    for k in RULE_NUMS:
        v = r.get(k)
        out[k] = int(v) if str(v or "").strip().lstrip("-").isdigit() else None
    out["require_sequence"] = bool(r.get("require_sequence", True))
    out["only_manual"] = bool(r.get("only_manual"))          # a dataset made from a selection in the genome table: exactly the genomes listed in manual_include
    return out


def describe(rules: dict) -> str:
    r = normalize(rules)
    names = {"clade": "clade", "lineage": "lineage", "outbreak": "outbreak", "region": "대륙", "country": "국가", "host": "숙주", "completeness": "완성도", "qc_not": "QC 제외"}
    parts = ["고른 유전체만"] if r.get("only_manual") else []
    parts += [f"{names[k]} {'/'.join(r[k])}" for k in RULE_LISTS if r[k]]
    if r["year_min"] or r["year_max"]:
        parts.append(f"수집연도 {r['year_min'] or ''}–{r['year_max'] or ''}")
    if r["min_length"] or r["max_length"]:
        parts.append(f"길이 {r['min_length'] or ''}–{r['max_length'] or ''}bp")
    if not r["require_sequence"]:
        parts.append("서열 없는 것 포함")
    return " · ".join(parts) or "제한 없음 (전체)"


def _year(row: dict):
    d = str(row.get("date") or "")
    return int(d[:4]) if d[:4].isdigit() else None


def matches(ds: dict, row: dict, genome: dict | None = None, key: str = "") -> bool:
    """Is this export row (nextstrain profile) a member? `genome` is the curated genome record (completeness, QC) when there is one."""
    ident = key or row.get("PPX_accession") or ""
    ids = {ident, row.get("INSDC_accession", ""), row.get("accession", "")} - {""}
    if ids & set(ds.get("manual_exclude") or []):
        return False
    if ids & set(ds.get("manual_include") or []):
        return True
    r = normalize(ds.get("rules"))
    if r["only_manual"]:
        return False
    for k, field in (("clade", "clade"), ("lineage", "lineage"), ("outbreak", "outbreak"), ("region", "region"), ("country", "country"), ("host", "host")):
        if r[k] and str(row.get(field) or "") not in r[k]:
            return False
    g = genome or {}
    if r["completeness"] and g.get("completeness") not in r["completeness"]:
        return False
    y = _year(row)
    if (r["year_min"] or r["year_max"]) and (y is None or (r["year_min"] and y < r["year_min"]) or (r["year_max"] and y > r["year_max"])):
        return False
    ln = int(row.get("length") or 0)
    if (r["min_length"] and ln < r["min_length"]) or (r["max_length"] and ln > r["max_length"]):
        return False
    if r["qc_not"] and g.get("_qc"):
        from . import qc
        if qc.overall({t: v["status"] for t, v in g["_qc"].items()}) in set(r["qc_not"]):
            return False
    return True


def members_hash(ids: list[str]) -> str:
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()


# ------------------------------------------------------------------------------------------------ storage
def _row(r) -> dict:
    return {"id": r["id"], "watch_id": r["watch_id"], "name": r["name"], "description": r["description"] or "", "rules": normalize(json.loads(r["rules_json"] or "{}")),
            "include_private": bool(r["include_private"]), "manual_include": json.loads(r["manual_include_json"] or "[]"),
            "manual_exclude": json.loads(r["manual_exclude_json"] or "[]"), "build": r["build"] or "", "created": r["created"], "updated": r["updated"],
            "params_yaml": r["params_yaml"] or "", "squirrel_json": r["squirrel_json"] or "", "mutprofile_json": r["mutprofile_json"] or ""}


def save(watch_id: int, name: str, rules: dict | None = None, description: str = "", include_private: bool = False, manual_include: list | None = None,
         manual_exclude: list | None = None, build: str = "", dataset_id: int | None = None) -> dict:
    """Create a dataset (or update the one with dataset_id). Names are unique per virus."""
    name = (name or "").strip()
    if not name:
        raise ValueError("데이터셋 이름이 필요합니다.")
    if build and build not in BUILDS:
        raise ValueError(f"빌드는 {', '.join(BUILDS)} 중 하나입니다.")
    clean = lambda xs: [str(x).strip() for x in (xs or []) if str(x).strip()]
    con = connect()
    try:
        if not con.execute("SELECT 1 FROM watchlist WHERE id=?", (watch_id,)).fetchone():
            raise ValueError("바이러스를 찾을 수 없습니다.")
        now = datetime.now().isoformat(timespec="seconds")
        vals = (name, description.strip(), json.dumps(normalize(rules)), 1 if include_private else 0, json.dumps(clean(manual_include)), json.dumps(clean(manual_exclude)), build or "")
        if dataset_id:
            if not con.execute("SELECT 1 FROM dataset WHERE id=? AND watch_id=?", (dataset_id, watch_id)).fetchone():
                raise ValueError("데이터셋을 찾을 수 없습니다.")
            clash = con.execute("SELECT 1 FROM dataset WHERE watch_id=? AND name=? AND id!=?", (watch_id, name, dataset_id)).fetchone()
            if clash:
                raise ValueError(f"같은 이름의 데이터셋이 이미 있습니다: {name}")
            con.execute("UPDATE dataset SET name=?, description=?, rules_json=?, include_private=?, manual_include_json=?, manual_exclude_json=?, build=?, updated=? WHERE id=?", vals + (now, dataset_id))
            did = dataset_id
        else:
            if con.execute("SELECT 1 FROM dataset WHERE watch_id=? AND name=?", (watch_id, name)).fetchone():
                raise ValueError(f"같은 이름의 데이터셋이 이미 있습니다: {name}")
            did = con.execute("INSERT INTO dataset(watch_id,name,description,rules_json,include_private,manual_include_json,manual_exclude_json,build,created,updated) VALUES(?,?,?,?,?,?,?,?,?,?)",
                              (watch_id,) + vals + (now, now)).lastrowid
        audit(con, "dataset_save", f"{name} ({'수정' if dataset_id else '생성'}): {describe(rules or {})}")
        con.commit()
        return _row(con.execute("SELECT * FROM dataset WHERE id=?", (did,)).fetchone())
    finally:
        con.close()


def set_params(dataset_id: int, text: str) -> None:
    """Store the dataset's Nextstrain settings (YAML overrides of the workflow's config). Validation is nextstrain.parse_params's job."""
    con = connect()
    try:
        con.execute("UPDATE dataset SET params_yaml=?, updated=? WHERE id=?", (text, datetime.now().isoformat(timespec="seconds"), dataset_id))
        con.commit()
    finally:
        con.close()


def set_mutprofile(dataset_id: int, text: str) -> None:
    """Store the dataset's mutation-profile settings (JSON; validated by mutprofile.parse_settings)."""
    con = connect()
    try:
        con.execute("UPDATE dataset SET mutprofile_json=?, updated=? WHERE id=?", (text, datetime.now().isoformat(timespec="seconds"), dataset_id))
        con.commit()
    finally:
        con.close()


def set_squirrel(dataset_id: int, text: str) -> None:
    """Store the dataset's Squirrel settings (JSON; validated by squirrel.parse_settings)."""
    con = connect()
    try:
        con.execute("UPDATE dataset SET squirrel_json=?, updated=? WHERE id=?", (text, datetime.now().isoformat(timespec="seconds"), dataset_id))
        con.commit()
    finally:
        con.close()


def get(dataset_id: int) -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT * FROM dataset WHERE id=?", (dataset_id,)).fetchone()
        return _row(r) if r else None
    finally:
        con.close()


def find(watch_id: int, name: str) -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT * FROM dataset WHERE watch_id=? AND name=?", (watch_id, name)).fetchone()
        return _row(r) if r else None
    finally:
        con.close()


def listing(watch_id: int | None = None) -> list[dict]:
    con = connect()
    try:
        q = "SELECT d.*, (SELECT MAX(exported) FROM dataset_export e WHERE e.dataset_id=d.id) last_export FROM dataset d" + (" WHERE d.watch_id=?" if watch_id else "") + " ORDER BY d.name"
        out = []
        for r in con.execute(q, (watch_id,) if watch_id else ()):
            d = _row(r)
            d["last_export"] = r["last_export"]
            d["summary"] = describe(d["rules"])
            out.append(d)
        return out
    finally:
        con.close()


def delete(dataset_id: int) -> bool:
    """Delete the DEFINITION only (its export history too); genomes, sequences and exported folders are untouched."""
    con = connect()
    try:
        r = con.execute("SELECT name FROM dataset WHERE id=?", (dataset_id,)).fetchone()
        if not r:
            return False
        con.execute("DELETE FROM dataset_export WHERE dataset_id=?", (dataset_id,))
        con.execute("DELETE FROM dataset WHERE id=?", (dataset_id,))
        audit(con, "dataset_delete", f"{r['name']} (정의만 삭제, 유전체·내보낸 파일은 그대로)")
        con.commit()
        return True
    finally:
        con.close()


def record_export(dataset_id: int, directory: str, rows: int, private_rows: int, sha: str, rules: dict) -> None:
    con = connect()
    try:
        con.execute("INSERT INTO dataset_export(dataset_id,exported,dir,rows,private_rows,members_sha256,rules_json) VALUES(?,?,?,?,?,?,?)",
                    (dataset_id, datetime.now().isoformat(timespec="seconds"), directory, rows, private_rows, sha, json.dumps(rules, ensure_ascii=False)))
        con.commit()
    finally:
        con.close()


def export_history(dataset_id: int, limit: int = 10) -> list[dict]:
    con = connect()
    try:
        return [dict(r) for r in con.execute("SELECT exported, dir, rows, private_rows, members_sha256 FROM dataset_export WHERE dataset_id=? ORDER BY id DESC LIMIT ?", (dataset_id, limit))]
    finally:
        con.close()


def member_ids(watch_id: int, ds: dict) -> tuple[set[str], set[str]]:
    """(public accessions, private ids) that the dataset's rules and manual lists pick - for browsing and the dashboard. Unlike an export this does
    NOT require a sequence file (the view shows who lacks one); genomes marked excluded are not members."""
    from . import export, exclusions, private
    _w, rows = export._load(watch_id)
    pub = {g["accession"] for g in rows if not g.get("_excluded") and matches(ds, export.nextstrain_row(g), g)}
    priv: set[str] = set()
    if ds.get("include_private"):
        ex = exclusions.marked(watch_id, "private")
        for r in private.rows_for_export(watch_id):
            if r["id"] not in ex and matches(ds, export.private_row(r), None, key=r["id"]):
                priv.add(r["id"])
    return pub, priv


def run_qc(ds: dict, progress=None, force: bool = False) -> dict:
    """QC (sequence statistics, Nextclade) for exactly this dataset's members - the public genomes that have a sequence file and, when the
    dataset includes them, the private ones - instead of the whole virus."""
    from . import qc
    p = progress or (lambda *a: None)
    pub, priv = member_ids(ds["watch_id"], ds)
    res = {"dataset": ds["name"], "public": None, "private": None}
    res["public"] = qc.run(ds["watch_id"], lambda pct, m: p(pct * (0.8 if priv else 1), m), force, accessions=pub)
    if priv:
        res["private"] = qc.run_private(ds["watch_id"], lambda pct, m: p(80 + pct * 0.2, m), force, ids=priv)
    return res


# ------------------------------------------------------------------------------------------------ what values exist (for the rule builder)
def facets(watch_id: int) -> dict:
    """Value counts over the genomes the dataset rules can pick from (curated, not Removed, not marked excluded) - and the private store -
    split into 'with a sequence file' (what an export needs) and all."""
    from . import config, export, exclusions, private
    w, rows = export._load(watch_id)
    out: dict[str, Counter] = {k: Counter() for k in ("clade", "lineage", "outbreak", "region", "country", "host", "completeness", "year")}
    seq: dict[str, Counter] = {k: Counter() for k in out}
    n = n_seq = 0
    for g in rows:
        if g.get("_excluded"):
            continue
        r = export.nextstrain_row(g)
        has = bool(g.get("file_path") and config.data_path(g["file_path"]).exists())
        n += 1
        n_seq += has
        for k in out:
            v = g.get("completeness") if k == "completeness" else (_year(r) if k == "year" else r.get(k))
            v = "" if v is None else str(v)
            out[k][v] += 1
            if has:
                seq[k][v] += 1
    ex_priv = exclusions.marked(watch_id, "private")
    pn = sum(1 for p in private.rows_for_export(watch_id) if p["id"] not in ex_priv)
    return {"virus": w["name"], "taxon": w["taxon"], "genomes": n, "with_sequence": n_seq, "private": pn,
            "fields": {k: [{"value": v, "n": c, "with_seq": seq[k].get(v, 0)} for v, c in sorted(out[k].items(), key=lambda x: (-x[1], x[0]))][:200] for k in out},
            "presets": PRESETS.get((w["taxon"] or "").lower(), [])}
