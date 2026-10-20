"""Write the record `params.gid` of the multi-FASTA `input[0]` to `output[0]` (plain Python; Snakemake's `script:` gives the `snakemake` object)."""
gid = snakemake.params.gid
out, keep = [], False
for line in open(snakemake.input[0], encoding="utf-8"):
    if line.startswith(">"):
        keep = line[1:].split()[0] == gid
        if keep:
            out.append(f">{gid}\n")
        elif out:
            break
    elif keep:
        out.append(line)
open(snakemake.output[0], "w", encoding="utf-8").write("".join(out))
