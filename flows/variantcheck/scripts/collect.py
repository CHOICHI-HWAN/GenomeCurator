"""Join the per-genome call tables into `variantcall.tsv` (same columns as genomecurator/variantcall.py writes) and list the genomes that were called."""
rows, genomes = [], []
for f in sorted(snakemake.input):
    text = open(f, encoding="utf-8").read()
    rows.append(text)
    first = text.split("\t", 1)[0]
    if first:
        genomes.append(first)
ids = {}
for f in sorted(snakemake.input):
    sid = f.split("/")[-1][:-4]
    ids[sid] = None
with open(snakemake.output.tsv, "w", encoding="utf-8", newline="") as out:
    out.write("genome\tpos\tref\talt\tsnippy\tlofreq\tdp\taf\n")
    out.write("".join(rows))
# a genome with no calls at all has an empty table: its id is what the flow was asked about, so list every input genome
all_ids = []
for line in open(snakemake.config["fasta"], encoding="utf-8"):
    if line.startswith(">"):
        all_ids.append(line[1:].split()[0])
open(snakemake.output.genomes, "w", encoding="utf-8").write("\n".join(sorted(all_ids)) + "\n")
