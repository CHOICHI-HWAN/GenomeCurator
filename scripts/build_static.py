"""Build a public, read-only snapshot (docs/index.html) from the local database.

Only metadata that is already public (NCBI Virus / Pathoplexus / BV-BRC) is exported.
No sequences, file paths, curation state, audit log or collection code is included.

    python scripts/build_static.py [--out docs/index.html] [--repo https://github.com/<user>/<repo>]
"""
import argparse
import json
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomecurator.db import connect  # noqa: E402  (read-only use)

FIELDS = [
    "genome_id", "accession", "accession_version", "genome_name", "strain", "genome_status", "genome_length",
    "collection_date", "isolation_country", "geographic_location", "host_name", "isolation_source",
    "sra_accession", "bioproject_accession", "biosample_accession", "sequencing_platform", "assembly_method",
    "genome_quality", "bvbrc_genome_id", "release_date", "taxon_id", "genotype", "subtype", "segment", "publication",
]


def build(out: Path, repo: str) -> None:
    con = connect()
    viruses = []
    for w in con.execute("SELECT * FROM watchlist ORDER BY id"):
        rows = con.execute(
            "SELECT g.*, s.gc_content, s.md5 FROM genome g LEFT JOIN genome_sequence s ON s.accession=g.accession "
            "WHERE g.watch_id=? AND IFNULL(g.genome_status,'')!='Removed' ORDER BY g.genome_id", (w["id"],)).fetchall()
        if not rows:
            continue
        prov: dict[str, dict] = {}
        for p in con.execute("SELECT p.genome_id, p.field, p.source FROM field_provenance p JOIN genome g "
                             "ON g.genome_id=p.genome_id WHERE g.watch_id=?", (w["id"],)):
            prov.setdefault(p["genome_id"], {})[p["field"]] = p["source"]
        out_rows = []
        for r in rows:
            d = {k: r[k] for k in FIELDS if r[k] not in (None, "")}
            d["genbank"] = f"{r['accession']}.{r['accession_version']}" if r["accession_version"] else r["accession"]
            if r["gc_content"] is not None:
                d["gc_content"] = r["gc_content"]
            if r["genome_id"] in prov:
                d["p"] = prov[r["genome_id"]]
            out_rows.append(d)
        lineage = (rows[0]["taxon_lineage_names"] or "").split(",") if rows[0]["taxon_lineage_names"] else [w["name"]]
        viruses.append({"id": w["id"], "name": w["name"], "taxon": w["taxon"], "taxon_id": w["taxon_id"],
                        "lineage": lineage, "rows": out_rows})
    if not viruses:
        sys.exit("no genomes in the database; nothing to publish")
    payload = {"generated": date.today().isoformat(), "repo": repo, "viruses": viruses}
    blob = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    html = (ROOT / "scripts" / "static_template.html").read_text(encoding="utf-8").replace("/*DATA*/", blob)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    n = sum(len(v["rows"]) for v in viruses)
    print(f"wrote {out} ({out.stat().st_size / 1e6:.2f} MB, {len(viruses)} viruses, {n} genomes)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "docs" / "index.html"))
    ap.add_argument("--repo", default="https://github.com/CHOICHI-HWAN/GenomeCurator")
    a = ap.parse_args()
    build(Path(a.out), a.repo)
