"""Source adapters. Each fetch_* returns {base_accession: partial record} using
BV-BRC style field names. A missing source (e.g. organism not on Pathoplexus)
yields {} instead of an error."""
import io
import re
import json
import time
import zipfile
from urllib.parse import quote

import httpx

from . import config

_last_call = 0.0


def _client() -> httpx.Client:
    return httpx.Client(timeout=180, follow_redirects=True,
                        headers={"User-Agent": "GenomeCurator/0.1"})


def _throttle():
    """Space out every outgoing web request (default 300 ms, GC_REQUEST_DELAY_MS)."""
    global _last_call
    wait = config.REQUEST_DELAY_S - (time.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.monotonic()


def _get(c: httpx.Client, url: str, retries: int = 3, **kw) -> httpx.Response:
    """GET with the shared delay; retries 429/5xx and network errors with backoff."""
    for attempt in range(retries + 1):
        _throttle()
        try:
            r = c.get(url, **kw)
            if r.status_code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return r
        except httpx.TransportError:
            if attempt == retries:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError("unreachable")

def split_version(acc: str) -> tuple[str, int | None]:
    if "." in acc:
        base, v = acc.rsplit(".", 1)
        if v.isdigit():
            return base, int(v)
    return acc, None


def _ncbi_get(c: httpx.Client, url: str, **params):
    if config.NCBI_API_KEY:
        params["api_key"] = config.NCBI_API_KEY
    r = _get(c, url, params=params)
    r.raise_for_status()
    return r


# ---------------------------------------------------------------- NCBI Virus
def ncbi_resolve_taxon_id(taxon: str) -> int | None:
    try:
        with _client() as c:
            r = _ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{quote(taxon)}")
        return int(r.json()["taxonomy_nodes"][0]["taxonomy"]["tax_id"])
    except Exception:
        return None


def ncbi_resolve_virus(text: str) -> tuple[int, str]:
    """Resolve user text to an NCBI *virus* taxon -> (tax_id, official name).

    Raises ValueError with suggestions when the text is not a virus (e.g. 'andes' matches the genus
    Andes, not Andes virus) or when NCBI Virus holds no genome records for it."""
    text = text.strip()
    with _client() as c:
        node = None
        try:
            nodes = _ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{quote(text)}").json().get("taxonomy_nodes", [])
            node = nodes[0]["taxonomy"] if nodes and "taxonomy" in nodes[0] else None
        except Exception:
            pass
        if not node or 10239 not in node.get("lineage", []):
            hint = ""
            try:
                sug = _ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon_suggest/{quote(text)}").json()
                names = [f"{x['sci_name']} ({x['matched_term']})" if x["matched_term"] != x["sci_name"] else x["sci_name"]
                         for x in sug.get("sci_name_and_ids", []) if x.get("group_name") == "viruses"][:5]
                if names:
                    hint = " 후보: " + ", ".join(names)
            except Exception:
                pass
            what = f"NCBI에서 '{node['organism_name']}'(바이러스 아님)로 인식됨." if node else "NCBI에서 찾을 수 없음."
            raise ValueError(f"'{text}'는 바이러스 분류군이 아닙니다. {what}{hint}")
        tid, name = int(node["tax_id"]), node["organism_name"]
        total = _ncbi_get(c, f"{config.NCBI_BASE}/taxon/{tid}/dataset_report", page_size=1).json().get("total_count") or 0
        if not total:
            raise ValueError(f"'{name}'(taxon {tid})에는 NCBI Virus 유전체 레코드가 없습니다.")
        return tid, name


def ncbi_species_scope(taxon_id: int | None) -> tuple[int, str] | None:
    """If taxon_id is a rank-less virus name / strain / sub-species node, return its SPECIES ancestor.

    NCBI tags many records at the species node (e.g. 'Orthopoxvirus monkeypox' 3431483) while the
    well-known name ('Monkeypox virus' 10244) is a child, so querying the child misses records."""
    if not taxon_id:
        return None
    try:
        with _client() as c:
            node = _ncbi_get(c, f"https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/{taxon_id}").json()[
                "taxonomy_nodes"][0]["taxonomy"]
            rank = node.get("rank")
            if rank in ("SPECIES", "GENUS", "FAMILY", "SUBFAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "REALM"):
                return None
            ids = [i for i in node.get("lineage", []) if i > 1][-8:]
            if not ids:
                return None
            anc = _ncbi_get(c, "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(map(str, ids))).json()
            for n in reversed(anc["taxonomy_nodes"]):
                tx = n["taxonomy"]
                if tx.get("rank") == "SPECIES":
                    return int(tx["tax_id"]), tx["organism_name"]
    except Exception:
        pass
    return None


_STATUS = {"COMPLETE": "Complete", "PARTIAL": "Partial"}


def normalize_ncbi(rep: dict) -> tuple[str, dict]:
    base, ver = split_version(rep["accession"])
    iso = rep.get("isolate") or {}
    loc = rep.get("location") or {}
    virus = rep.get("virus") or {}
    sub = rep.get("submitter") or {}
    geo = loc.get("geographic_location") or ""
    refseq = base if rep.get("source_database") == "RefSeq" else None
    rec = {
        "accession_version": ver,
        "genbank_accessions": None if refseq else base,
        "refseq_accessions": refseq,
        "taxon_id": virus.get("tax_id"),
        "genome_name": virus.get("organism_name"),
        "taxon_lineage_ids": ",".join(str(x["tax_id"]) for x in virus.get("lineage", [])) or None,
        "taxon_lineage_names": ",".join(x["name"] for x in virus.get("lineage", [])) or None,
        "genome_status": _STATUS.get(rep.get("completeness"), rep.get("completeness")),
        "genome_length": rep.get("length"),
        "strain": iso.get("name"),
        "host_name": (rep.get("host") or {}).get("organism_name"),
        "isolation_country": geo.split(":")[0].strip() or None,
        "geographic_location": geo or None,
        "collection_date": iso.get("collection_date"),
        "isolation_source": iso.get("source"),
        "submitter": sub.get("affiliation"),
        "sequence_hash": (rep.get("nucleotide") or {}).get("sequence_hash"),
        "sra_accession": "; ".join(rep.get("sra_accessions") or []),
        "bioproject_accession": "; ".join(rep.get("bioprojects") or []),
        "biosample_accession": rep.get("biosample"),
        "release_date": (rep.get("release_date") or "")[:10] or None,
        "update_date": (rep.get("update_date") or "")[:10] or None,
    }
    return base, {k: v for k, v in rec.items() if v not in (None, "")}


def fetch_ncbi(taxon: str, filters: dict | None = None) -> dict[str, dict]:
    out: dict[str, dict] = {}
    params = {"page_size": 1000, **(filters or {})}
    with _client() as c:
        while True:
            data = _ncbi_get(c, f"{config.NCBI_BASE}/taxon/{quote(taxon)}/dataset_report", **params).json()
            for rep in data.get("reports", []):
                base, rec = normalize_ncbi(rep)
                out[base] = rec
            tok = data.get("next_page_token")
            if not tok:
                return out
            params["page_token"] = tok


def download_ncbi_fasta(accessions: list[str], progress=None) -> dict[str, str]:
    """accessions may carry .version. Returns {accession-with-version: sequence}."""
    seqs: dict[str, str] = {}
    with _client() as c:
        for i in range(0, len(accessions), 100):
            batch = accessions[i:i + 100]
            if progress:
                progress(i, len(accessions))
            r = _ncbi_get(c, f"{config.NCBI_BASE}/accession/{','.join(batch)}/genome/download",
                          include_sequence="GENOME")
            with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                name = next((n for n in z.namelist() if n.endswith("genomic.fna")), None)
                if not name:
                    continue
                text = z.read(name).decode()
            seqs.update(parse_fasta(text))
    return seqs


def parse_fasta(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    cur, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                out[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf)
    return out


# -------------------------------------------------------------------- BV-BRC
def fetch_bvbrc(taxon_id: int | None) -> dict[str, dict]:
    if not taxon_id:
        return {}
    out: dict[str, dict] = {}
    offset, page = 0, 1000
    with _client() as c:
        while True:
            r = _get(c, f"{config.BVBRC_BASE}?eq(taxon_lineage_ids,{taxon_id})&sort(+genome_id)&limit({page},{offset})",
                      headers={"Accept": "application/json"})
            r.raise_for_status()
            rows = r.json()
            for g in rows:
                host = g.get("host_name") or g.get("host_common_name")
                rec = {
                    "genome_name": g.get("genome_name"),
                    "genome_status": g.get("genome_status"),
                    "genome_length": g.get("genome_length"),
                    "strain": g.get("strain"),
                    "host_name": host,
                    "isolation_country": g.get("isolation_country"),
                    "geographic_location": g.get("geographic_location"),
                    "collection_date": g.get("collection_date"),
                    "isolation_source": g.get("isolation_source"),
                    "subtype": g.get("subtype"),
                    "genotype": g.get("genotype"),
                    "segment": g.get("segment"),
                    "publication": g.get("publication"),
                    "genome_quality": g.get("genome_quality"),
                    "sra_accession": g.get("sra_accession"),
                    "bioproject_accession": g.get("bioproject_accession"),
                    "biosample_accession": g.get("biosample_accession"),
                    "assembly_method": g.get("assembly_method"),
                    "sequencing_platform": g.get("sequencing_platform"),
                    "taxon_id": g.get("taxon_id"),
                    "taxon_lineage_ids": ",".join(g.get("taxon_lineage_ids") or []) or None,
                    "taxon_lineage_names": ",".join(g.get("taxon_lineage_names") or []) or None,
                }
                rec = {k: v for k, v in rec.items() if v not in (None, "")}
                for acc in (g.get("genbank_accessions") or "").replace(";", ",").split(","):
                    acc = acc.strip()
                    if acc:
                        out[split_version(acc)[0]] = {**rec, "bvbrc_genome_id": g.get("genome_id")}
            if len(rows) < page:
                return out
            offset += page


# --------------------------------------------------------------- Pathoplexus
# LAPIS has no organism list endpoint; these slugs were verified live (2026-10).
PPX_KNOWN = ["mpox", "ebola-zaire", "ebola-sudan", "ebola-bdbv", "marburg", "cchf", "west-nile", "hmpv", "rsv-a", "rsv-b",
             "measles", "dengue", "zika", "chikungunya", "yellow-fever", "andv"]
_PPX_FIELDS = ["insdcAccessionFull", "insdcAccessionBase", "accessionVersion", "specimenCollectorSampleId",
               "hostNameScientific", "hostNameCommon", "geoLocCountry", "geoLocAdmin1", "sampleCollectionDate",
               "authorAffiliations", "earliestReleaseDate", "anatomicalPart", "bodyProduct",
               "bioprojectAccession", "biosampleAccession", "insdcRawReadsAccession"]
# Organism-specific classification fields (mpox: clade/lineage, dengue: serotype/lineage, measles: genotype ...).
# LAPIS rejects unknown field names, so the ones an organism really has are discovered first.
_PPX_TAXON = ["clade", "lineage", "genotype", "subtype", "serotype", "outbreak"]
_PPX_FILTER = {"versionStatus": "LATEST_VERSION", "isRevocation": "false"}
_ppx_cache: dict[tuple, tuple | None] = {}


PPX_HOME = "https://pathoplexus.org/"
_ppx_index: dict[str, dict[int, int]] | None = None
_lineage_cache: dict[int, set[int]] = {}
_ppx_match: dict[tuple[str, int], int] = {}


def pathoplexus_organisms() -> list[str]:
    """Organisms Pathoplexus serves. LAPIS has no listing endpoint, so the slugs are read from the
    website's own navigation (new organisms such as ebola-bdbv or andv appear without a code change);
    PPX_KNOWN is only the fallback when the site cannot be read."""
    slugs = set(PPX_KNOWN)
    try:
        with _client() as c:
            r = _get(c, PPX_HOME, headers={"User-Agent": "Mozilla/5.0 GenomeCurator"})
            if r.status_code == 200:
                slugs |= set(re.findall(r'href="/([a-z0-9][a-z0-9\-]+)/search', r.text))
    except Exception:
        pass
    return sorted(slugs)


def pathoplexus_index(refresh: bool = False) -> dict[str, dict[int, int]]:
    """{organism slug: {ncbiVirusTaxId: number of sequences}} for every organism Pathoplexus serves."""
    global _ppx_index
    if _ppx_index is not None and not refresh:
        return _ppx_index
    index: dict[str, dict[int, int]] = {}
    with _client() as c:
        for slug in pathoplexus_organisms():
            try:
                r = _get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/aggregated", params={**_PPX_FILTER, "fields": "ncbiVirusTaxId"})
                if r.status_code != 200:
                    continue
                index[slug] = {int(x["ncbiVirusTaxId"]): x["count"] for x in r.json()["data"] if x.get("ncbiVirusTaxId") is not None}
            except Exception:
                continue
    _ppx_index = index
    return index


def _lineages(taxids: list[int]) -> dict[int, set[int]]:
    """NCBI lineage (ancestor ids, incl. itself) per taxon id, cached."""
    todo = [t for t in taxids if t not in _lineage_cache]
    with _client() as c:
        for i in range(0, len(todo), 50):
            chunk = todo[i:i + 50]
            try:
                nodes = _ncbi_get(c, "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/" + ",".join(map(str, chunk))).json().get("taxonomy_nodes", [])
            except Exception:
                continue
            for n in nodes:
                t = n.get("taxonomy") or {}
                if t.get("tax_id"):
                    _lineage_cache[int(t["tax_id"])] = {int(x) for x in t.get("lineage", [])} | {int(t["tax_id"])}
    return {t: _lineage_cache.get(t, {t}) for t in taxids}


def pathoplexus_matched(slug: str, taxon_id: int | None) -> int:
    """Sequences of this organism that belong to the taxon (the taxon itself or any descendant)."""
    return _ppx_match.get((slug, taxon_id or 0), 0)


def pathoplexus_resolve(watch_name: str, taxon: str, taxon_id: int | None):
    """-> (slug, filter) or None. Config override wins; otherwise the organism whose sequences are tagged with the
    watch's NCBI taxon or one of its descendants (Pathoplexus tags at varying ranks, e.g. a strain or the species)."""
    key = (watch_name, taxon_id)
    if key in _ppx_cache:
        return _ppx_cache[key]
    base = {**_PPX_FILTER}
    forced = config.PATHOPLEXUS_ORGANISMS.get(watch_name)
    found = None
    if not forced and taxon_id:
        index = pathoplexus_index()
        best: dict[str, int] = {}
        lin = _lineages([t for ids in index.values() for t in ids])
        for slug, ids in index.items():
            n = sum(cnt for t, cnt in ids.items() if t == taxon_id or taxon_id in lin.get(t, ()))
            if n:
                best[slug] = n
                _ppx_match[(slug, taxon_id)] = n
        if best:
            slug = max(best, key=best.get)
            found = (slug, {**base, "ncbiVirusTaxId": taxon_id} if taxon_id in index[slug] else base)
    if not found and not forced:  # last resort: match by virus name
        with _client() as c:
            for slug in pathoplexus_index():
                try:
                    r = _get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/aggregated", params={**base, "ncbiVirusName": taxon})
                    if r.status_code == 200 and r.json()["data"][0]["count"] > 0:
                        found = (slug, {**base, "ncbiVirusName": taxon})
                        break
                except Exception:
                    pass
    if not found and forced:  # forced organism: take everything
        found = (forced, base)
    _ppx_cache[key] = found
    return found


def pathoplexus_slug_for(accession_version: str) -> str | None:
    """Which Pathoplexus organism holds this accession (e.g. PP_006WPN1.1)? Probes the known organisms."""
    with _client() as c:
        for slug in pathoplexus_organisms():
            try:
                r = _get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/aggregated", params={"accessionVersion": accession_version})
                if r.status_code == 200 and r.json()["data"][0]["count"] > 0:
                    return slug
            except Exception:
                continue
    return None


def fetch_pathoplexus(watch_name: str, taxon: str = "", taxon_id: int | None = None) -> dict[str, dict]:
    res = pathoplexus_resolve(watch_name, taxon or watch_name, taxon_id)
    if not res:
        return {}
    slug, _ = res
    flt = dict(_PPX_FILTER)  # whole organism: NCBI defines the accession set, extra rows simply don't join
    out: dict[str, dict] = {}
    with _client() as c:
        probe = _get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details", params={**flt, "limit": 1})
        probe.raise_for_status()
        have = set((probe.json().get("data") or [{}])[0])
        fields = [f for f in _PPX_FIELDS if f in have] + [f for f in _PPX_TAXON if f in have]
        r = _get(c, f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details",
                  params={**flt, "fields": ",".join(fields), "limit": 200000}, timeout=300)
        r.raise_for_status()
        for row in r.json().get("data", []):
            full = row.get("insdcAccessionFull") or row.get("insdcAccessionBase")
            if not full:
                continue  # Pathoplexus-native sequence without INSDC accession: cannot join
            base, ver = split_version(full)
            geo = ", ".join(x for x in (row.get("geoLocCountry"), row.get("geoLocAdmin1")) if x)
            rec = {"accession_version": ver,
                   "pathoplexus_accession": row.get("accessionVersion"),
                   "strain": row.get("specimenCollectorSampleId"),
                   "host_name": row.get("hostNameScientific") or row.get("hostNameCommon"),
                   "isolation_country": row.get("geoLocCountry"),
                   "geographic_location": geo or None,
                   "collection_date": row.get("sampleCollectionDate"),
                   "clade": row.get("clade"), "lineage": row.get("lineage"), "serotype": row.get("serotype"),
                   "genotype": row.get("genotype"), "subtype": row.get("subtype"), "outbreak": row.get("outbreak"),
                   "isolation_source": row.get("anatomicalPart") or row.get("bodyProduct"),
                   "submitter": row.get("authorAffiliations"),
                   "release_date": row.get("earliestReleaseDate"),
                   "bioproject_accession": row.get("bioprojectAccession"),
                   "biosample_accession": row.get("biosampleAccession"),
                   "sra_accession": row.get("insdcRawReadsAccession")}
            out[base] = {k: v for k, v in rec.items() if v not in (None, "", "None")}
    return out

def fetch_all(watch_name: str, taxon: str, taxon_id: int | None, filters: dict | None, progress=None):
    """Returns ({source: {acc: rec}}, notes). Pathoplexus/BV-BRC failures are non-fatal."""
    notes: list[str] = []
    p = progress or (lambda *a: None)
    p(8, "NCBI Virus 메타데이터 가져오는 중…")
    res = {"ncbi": fetch_ncbi(taxon, filters)}
    for name, fn, pct, label in (("pathoplexus", lambda: fetch_pathoplexus(watch_name, taxon, taxon_id), 45, "Pathoplexus 확인 중…"),
                                 ("bvbrc", lambda: fetch_bvbrc(taxon_id), 55, "BV-BRC 메타데이터 가져오는 중…")):
        p(pct, label)
        try:
            res[name] = fn()
            if not res[name]:
                notes.append(f"{name}: no records (skipped)")
        except Exception as e:  # keep NCBI result usable
            res[name] = {}
            notes.append(f"{name}: failed ({e})")
    return res, notes
