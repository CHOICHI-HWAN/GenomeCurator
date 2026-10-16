"""Read-only public snapshot (a single self-contained HTML page) built from the local database.

Only metadata that is already public (NCBI Virus / Pathoplexus / BV-BRC) is exported.
No sequences, file paths, curation state, audit log or collection code is included.
"""
import json
from datetime import date
from pathlib import Path

from .db import connect

TEMPLATE = Path(__file__).parent / "templates" / "static_template.html"
FIELDS = [
    "genome_id", "accession", "accession_version", "genome_name", "strain", "genome_status", "genome_length",
    "collection_date", "isolation_country", "geographic_location", "host_name", "isolation_source",
    "sra_accession", "bioproject_accession", "biosample_accession", "sequencing_platform", "assembly_method",
    "genome_quality", "bvbrc_genome_id", "release_date", "taxon_id", "genotype", "subtype", "segment", "publication",
    "pathoplexus_accession", "clade", "lineage", "serotype", "outbreak",
]


def snapshot(repo: str = "") -> dict:
    con = connect()
    try:
        viruses = []
        for w in con.execute("SELECT * FROM watchlist ORDER BY id").fetchall():
            rows = con.execute(
                "SELECT g.*, s.gc_content FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession "
                "WHERE g.watch_id=? AND IFNULL(g.genome_status,'')!='Removed' ORDER BY g.genome_id", (w["id"],)).fetchall()
            if not rows:
                continue
            prov: dict[str, dict] = {}
            for p in con.execute("SELECT p.genome_id, p.field, p.source FROM field_provenance p JOIN genome g "
                                 "ON g.genome_id=p.genome_id WHERE g.watch_id=?", (w["id"],)):
                prov.setdefault(p["genome_id"], {})[p["field"]] = p["source"]
            out = []
            for r in rows:
                d = {k: r[k] for k in FIELDS if r[k] not in (None, "")}
                d["genbank"] = f"{r['accession']}.{r['accession_version']}" if r["accession_version"] else r["accession"]
                if r["gc_content"] is not None:
                    d["gc_content"] = r["gc_content"]
                if r["genome_id"] in prov:
                    d["p"] = prov[r["genome_id"]]
                out.append(d)
            names = rows[0]["taxon_lineage_names"]
            viruses.append({"id": w["id"], "name": w["name"], "taxon": w["taxon"], "taxon_id": w["taxon_id"],
                            "lineage": names.split(",") if names else [w["name"]], "ppx": w["pathoplexus_organism"], "rows": out})
        return {"generated": date.today().isoformat(), "repo": repo, "viruses": viruses}
    finally:
        con.close()


def build_html(repo: str = "") -> str:
    payload = snapshot(repo)
    if not payload["viruses"]:
        raise ValueError("no genomes in the database; nothing to publish")
    blob = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    return TEMPLATE.read_text(encoding="utf-8").replace("/*DATA*/", blob)
