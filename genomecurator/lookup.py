"""Pre-flight check for a new Watchlist entry: look the typed virus up in ICTV, NCBI Virus, BV-BRC
and Pathoplexus, report what each one knows, and propose a better (official) name when it differs.

Read-only: nothing is written to the database or disk. Every source is isolated, so one failing
service never hides the others."""
from urllib.parse import quote

from . import config, sources

OLS_SEARCH = "https://www.ebi.ac.uk/ols4/api/search"


def _src(status: str, title: str, lines: list[str], **data) -> dict:
    return {"status": status, "title": title, "lines": lines, **data}


# ------------------------------------------------------------------------------- ICTV
def _ictv(q: str) -> dict:
    """ICTV taxonomy via the official EVORA/OLS4 ontology API. Species ids look like .../MSL41/ICTV...,
    virus names (e.g. 'monkeypox virus') like .../VMR...."""
    with sources._client() as c:
        r = sources._get(c, OLS_SEARCH, params={"ontology": "ictv", "q": q, "rows": 12,
                                                "queryFields": "label,synonym,short_form,obo_id"})
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", [])
    species, names = [], []
    for d in docs:
        iri, label = d.get("iri", ""), d.get("label", "")
        if not label or d.get("is_obsolete"):
            continue
        if "/MSL" in iri:
            species.append({"name": label, "release": iri.split("/MSL")[1].split("/")[0].replace("MSL", ""),
                            "id": iri.rsplit("/", 1)[-1]})
        elif "/VMR" in iri:
            names.append(label)
    if not species and not names:
        return _src("none", "ICTV", ["ICTV 분류에서 찾을 수 없음"], species=[], virus_names=[])
    lines = [f"종(species): {s['name']} (MSL{s['release']}, {s['id']})" for s in species[:5]]
    if names:
        lines.append("바이러스 명칭: " + ", ".join(sorted(set(names))[:5]))
    return _src("ok" if species else "warn", "ICTV", lines, species=species, virus_names=sorted(set(names)))


# ------------------------------------------------------------------------------- NCBI
def _ncbi(q: str, ictv_species: list[str]) -> dict:
    """Resolve to an NCBI *virus* taxon (trying ICTV's species names too), widen to species, count records."""
    with sources._client() as c:
        def node(text):
            try:
                nodes = sources._ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{quote(text)}").json().get("taxonomy_nodes", [])
                return nodes[0].get("taxonomy") if nodes else None
            except Exception:
                return None

        found, tried = None, []
        for text in [q] + [s for s in ictv_species if s.lower() != q.lower()]:
            n = node(text)
            if n and 10239 in n.get("lineage", []):
                found = n
                break
            tried.append(text if not n else f"{text} (→ {n.get('organism_name')}, 바이러스 아님)")
        hints = []
        if not found:
            try:
                sug = sources._ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest/{quote(q)}").json()
                hints = [(x["sci_name"], int(x["tax_id"]), x.get("matched_term", ""))
                         for x in sug.get("sci_name_and_ids", []) if x.get("group_name") == "viruses"][:5]
                if hints:
                    found = node(str(hints[0][1]))
            except Exception:
                pass
        if not found:
            return _src("none", "NCBI Virus", ["NCBI Taxonomy에서 바이러스 분류군을 찾을 수 없음"] + [f"시도: {t}" for t in tried[:3]], taxon=None)
        tid, name, rank = int(found["tax_id"]), found["organism_name"], found.get("rank") or "no rank"
        scope = sources.ncbi_species_scope(tid)
        sid, sname = scope if scope else (tid, name)
        total = sources._ncbi_get(c, f"{config.NCBI_BASE}/taxon/{sid}/dataset_report", page_size=1).json().get("total_count") or 0
    lines = [f"분류군: {name} (taxon {tid}, {rank})"]
    if scope:
        lines.append(f"종 노드로 확장: {sname} (taxon {sid}) — 종 노드에 붙은 레코드까지 모두 수집")
    lines.append(f"NCBI Virus 유전체 레코드: {total:,}건")
    status = "ok" if total else "warn"
    if hints and (not tried or tried[0].lower() != name.lower()):
        lines.insert(0, f"입력과 일치하는 항목이 없어 가장 유사한 후보를 사용: {hints[0][0]}" + (f" ({hints[0][2]})" if hints[0][2] else ""))
        status = "warn"
    return _src(status, "NCBI Virus", lines, taxon={"id": sid, "name": sname, "input_id": tid, "input_name": name}, count=total)


# ------------------------------------------------------------------------------- BV-BRC
def _bvbrc(q: str, taxon_id: int | None) -> dict:
    with sources._client() as c:
        tax = []
        try:
            r = sources._get(c, f"https://www.bv-brc.org/api/taxonomy/?eq(taxon_name,{quote(q)}*)&select(taxon_id,taxon_name,taxon_rank,genomes)&sort(-genomes)&limit(5)",
                             headers={"Accept": "application/json"})
            r.raise_for_status()
            tax = r.json()
        except Exception:
            pass
        count = None
        if taxon_id:
            r = sources._get(c, f"{config.BVBRC_BASE}?eq(taxon_lineage_ids,{taxon_id})&limit(1)&select(genome_id)",
                             headers={"Accept": "application/solr+json"})
            r.raise_for_status()
            count = r.json().get("response", {}).get("numFound")
    lines = []
    if count is not None:
        lines.append(f"taxon {taxon_id} 계통에 속한 BV-BRC 유전체: {count:,}건")
    if tax:
        lines.append("이름이 비슷한 BV-BRC 분류군: " + "; ".join(f"{t['taxon_name']} (taxon {t['taxon_id']}, {t.get('taxon_rank', '?')}, {t.get('genomes', 0):,}건)" for t in tax[:3]))
    if not lines:
        return _src("none", "BV-BRC", ["BV-BRC에서 찾을 수 없음"], count=0)
    return _src("ok" if count else "warn", "BV-BRC", lines, count=count or 0, taxa=tax)


# ------------------------------------------------------------------------------- Pathoplexus
def _pathoplexus(q: str, taxon: str, taxon_id: int | None) -> dict:
    res = sources.pathoplexus_resolve(q, taxon or q, taxon_id)
    if not res:
        return _src("none", "Pathoplexus", ["지원 organism에 없음 — 이 소스는 건너뜁니다 (오류 아님, NCBI·BV-BRC로 수집)"], organism=None)
    slug, flt = res
    with sources._client() as c:
        def count(params):
            r = sources._get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/aggregated", params=params)
            r.raise_for_status()
            return r.json()["data"][0]["count"]
        total = count(dict(sources._PPX_FILTER))
        mine = sources.pathoplexus_matched(slug, taxon_id) or count(dict(flt))
    lines = [f"organism: {slug} — 전체 {total:,}건을 받아 NCBI accession과 연결 (종 노드로 직접 표기된 건 {mine:,}건)"]
    return _src("ok", "Pathoplexus", lines, organism=slug, count=mine, total=total)


# ------------------------------------------------------------------------------- combine
_cache: dict[str, dict] = {}


def lookup_virus(text: str) -> dict:
    q = " ".join(text.split())
    if not q:
        raise ValueError("바이러스 이름을 입력하세요.")
    if q.lower() in _cache:
        return _cache[q.lower()]
    out: dict = {"query": q, "sources": {}, "suggestions": [], "add": None}

    def run(key, fn):
        try:
            out["sources"][key] = fn()
        except Exception as e:  # a failing service must not hide the others
            out["sources"][key] = _src("error", key, [f"조회 실패: {e}"])
        return out["sources"][key]

    ictv = run("ictv", lambda: _ictv(q))
    species = [s["name"] for s in ictv.get("species", [])]
    ncbi = run("ncbi", lambda: _ncbi(q, species))
    tx = ncbi.get("taxon")
    run("bvbrc", lambda: _bvbrc(q, tx["id"] if tx else None))
    run("pathoplexus", lambda: _pathoplexus(q, tx["name"] if tx else q, tx["id"] if tx else None))

    if tx:
        out["add"] = {"taxon": tx["name"], "taxon_id": tx["id"]}
        names = []
        official = next((s for s in species if s.lower() == tx["name"].lower()), None) or (species[0] if len(species) == 1 else None)
        if official and official.lower() != q.lower():
            names.append({"name": official, "reason": f"ICTV 공식 종명 · NCBI taxon {tx['id']}와 일치"})
        if tx["name"].lower() not in {q.lower(), (official or "").lower()}:
            names.append({"name": tx["name"], "reason": "NCBI 종 노드 학명"})
        if tx["input_name"].lower() not in {n["name"].lower() for n in names} | {q.lower()}:
            names.append({"name": tx["input_name"], "reason": "입력과 가장 가까운 NCBI 이름"})
        out["suggestions"] = names
        vn = ictv.get("virus_names") or []
        if official and any(v.lower() == q.lower() for v in vn):
            out["note"] = f"'{q}'는 ICTV에서 바이러스 명칭(약칭)이고, 공식 종명은 '{official}'입니다. 이름이 섞여 쓰이므로 종명으로 등록하는 것을 권합니다."
    if not any(s["status"] == "error" for s in out["sources"].values()):
        _cache[q.lower()] = out
    return out
