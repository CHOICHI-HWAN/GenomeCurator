"""Single-base substitutions of snippy's snps.vcf and lofreq's VCF of one genome -> `genome pos ref alt snippy lofreq dp af` rows."""


def parse(path):
    out = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 8 or len(p[3]) != 1 or len(p[4]) != 1 or p[4] == ".":
            continue
        info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
        out[(int(p[1]), p[3].upper(), p[4].upper())] = (int(info["DP"]) if info.get("DP", "").isdigit() else "", float(info["AF"]) if "AF" in info else "")
    return out


gid = snakemake.params.gid
sn = parse(f"{snakemake.input.snippy}/snps.vcf")
lf = parse(snakemake.input.lofreq)
rows = []
for v in sorted(set(sn) | set(lf)):
    dp, af = lf.get(v, ("", ""))
    rows.append((gid, *v, int(v in sn), int(v in lf), dp, af))
open(snakemake.output[0], "w", encoding="utf-8").write("".join("\t".join(map(str, r)) + "\n" for r in rows))
