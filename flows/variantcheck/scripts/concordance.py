"""Agreement of the mutation-profile run's own calls, snippy and lofreq (logic in genomecurator/flowlogic.py, shared with the app)."""
import json
import sys

cfg = snakemake.config
sys.path.insert(0, cfg["repo"])
from genomecurator import flowlogic          # noqa: E402  (pure Python, no third-party imports)

judged = set(open(snakemake.input.judged, encoding="utf-8").read().split()) & set(open(snakemake.input.genomes, encoding="utf-8").read().split())
lo, hi = cfg.get("itr_lo"), cfg.get("itr_hi")
inwin = (lambda p: lo <= p <= hi) if lo and hi else None
callers = {"mutprofile": flowlogic.read_mutations(snakemake.input.mutations, judged)}
if inwin:
    callers["mutprofile"] = {g: {v for v in vs if inwin(v[0])} for g, vs in callers["mutprofile"].items()}
callers.update(flowlogic.read_calls(snakemake.input.tsv, judged, inwin))
rep = flowlogic.sites_report(callers, 300)
all_sites = rep.pop("all_sites")
json.dump({"available": True, **rep}, open(snakemake.output.json, "w", encoding="utf-8"), ensure_ascii=False)
open(snakemake.output.bed, "w", encoding="utf-8").write(flowlogic.review_bed(all_sites, cfg.get("seqid", "ref"), int(cfg.get("min_genomes", 3))))
