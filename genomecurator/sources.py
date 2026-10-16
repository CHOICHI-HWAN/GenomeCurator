"""Source adapters. Each fetch_* returns {base_accession: partial record} using
BV-BRC style field names. A missing source (e.g. organism not on Pathoplexus)
yields {} instead of an error."""
import io
import json
import time
import zipfile
from urllib.parse import quote

import httpx

from . import config

_last_call = 0.0


def _client() -> httpx.Client:
    return httpx.Client(timeout=60, follow_redirects=True,
                        headers={"User-Agent": "GenomeCurator/0.1"})


def _throttle():
    """NCBI: 3 req/s without key, 10 with."""
    global _last_call
    gap = 0.11 if config.NCBI_API_KEY else 0.34
    wait = gap - (time.monotonic() - _last_call)
    if wait > 0:
        time.sleep(wait)
    _last_call = time.monotonic()


def split_version(acc: str) -> tuple[str, int | None]:
    if "." in acc:
        base, v = acc.rsplit(".", 1)
        if v.isdigit():
            return base, int(v)
    return acc, None


def _ncbi_get(c: httpx.Client, url: str, **params):
    _throttle()
    if config.NCBI_API_KEY:
        params["api_key"] = config.NCBI_API_KEY
    r = c.get(url, params=params)
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
            r = c.get(f"{config.BVBRC_BASE}?eq(taxon_lineage_ids,{taxon_id})&sort(+genome_id)&limit({page},{offset})",
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
PPX_KNOWN = ["mpox", "ebola-zaire", "ebola-sudan", "marburg", "cchf", "west-nile", "hmpv", "rsv-a", "rsv-b",
             "measles", "dengue", "zika", "chikungunya", "yellow-fever"]
_PPX_FIELDS = ["insdcAccessionFull", "insdcAccessionBase", "accessionVersion", "specimenCollectorSampleId",
               "hostNameScientific", "hostNameCommon", "geoLocCountry", "geoLocAdmin1", "sampleCollectionDate",
               "clade", "lineage", "authorAffiliations", "earliestReleaseDate", "anatomicalPart", "bodyProduct"]
_PPX_FILTER = {"versionStatus": "LATEST_VERSION", "isRevocation": "false"}
_ppx_cache: dict[tuple, tuple | None] = {}


def pathoplexus_resolve(watch_name: str, taxon: str, taxon_id: int | None):
    """-> (slug, filter) or None. Config override wins; otherwise probe known organisms by NCBI taxon."""
    key = (watch_name, taxon_id)
    if key in _ppx_cache:
        return _ppx_cache[key]
    base = {**_PPX_FILTER}
    found = None
    slug = config.PATHOPLEXUS_ORGANISMS.get(watch_name)
    cands = [(slug, None)] if slug else [(s, None) for s in PPX_KNOWN]
    with _client() as c:
        for s, _ in cands:
            for flt in ([{"ncbiVirusTaxId": taxon_id}] if taxon_id else []) + [{"ncbiVirusName": taxon}]:
                try:
                    r = c.get(f"{config.PATHOPLEXUS_BASE}/{s}/sample/aggregated", params={**base, **flt})
                    if r.status_code == 200 and r.json()["data"][0]["count"] > 0:
                        found = (s, {**base, **flt})
                        break
                except Exception:
                    pass
            if found:
                break
        if not found and slug:  # forced organism with no taxon match: take everything
            found = (slug, base)
    _ppx_cache[key] = found
    return found


def fetch_pathoplexus(watch_name: str, taxon: str = "", taxon_id: int | None = None) -> dict[str, dict]:
    res = pathoplexus_resolve(watch_name, taxon or watch_name, taxon_id)
    if not res:
        return {}
    slug, flt = res
    out: dict[str, dict] = {}
    with _client() as c:
        r = c.get(f"{config.PATHOPLEXUS_BASE}/{slug}/sample/details",
                  params={**flt, "fields": ",".join(_PPX_FIELDS), "limit": 200000}, timeout=300)
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
                   "genotype": row.get("clade"), "subtype": row.get("lineage"),
                   "isolation_source": row.get("anatomicalPart") or row.get("bodyProduct"),
                   "submitter": row.get("authorAffiliations"),
                   "release_date": row.get("earliestReleaseDate")}
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
