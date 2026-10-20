"""Pure-Python logic shared by the app and the Snakemake flows (`flows/*/scripts`): no imports from the rest of the app, no third-party packages, so a flow script can load it with a plain
`sys.path` entry from inside the small `gc-snakemake` environment and the app can import it as usual - one implementation, not two copies."""
import csv
from collections import Counter


def read_calls(path, judged: set, in_window=None) -> dict:
    """{caller: {genome: {(pos, ref, alt)}}} for snippy and lofreq from a `variantcall.tsv` (columns genome, pos, ref, alt, snippy, lofreq, ...); only genomes in `judged`, and, with
    `in_window(pos)`, only positions it accepts (the ITR is left out that way)."""
    out = {t: {g: set() for g in judged} for t in ("snippy", "lofreq")}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        if r["genome"] in judged and (in_window is None or in_window(int(r["pos"]))):
            v = (int(r["pos"]), r["ref"], r["alt"])
            for t in ("snippy", "lofreq"):
                if r[t] == "1":
                    out[t][r["genome"]].add(v)
    return out


def read_mutations(path, judged: set) -> dict:
    """{genome: {(pos, ref, alt)}} from a mutation-profile `mutations.tsv` (the run's own calls)."""
    out = {g: set() for g in judged}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        if r["sample"] in judged:
            out[r["sample"]].add((int(r["pos"]), r["ref"].upper(), r["alt"].upper()))
    return out


def concordance(callers: dict) -> dict:
    """How well callers agree. `callers` = {caller: {genome id: {(pos, ref, alt)}}} -> {"callers", "rows": [{id, pos, ref, alt, by, n, of}], "summary"}.
    A genome is only judged by the callers that looked at it; a variant is called by `n` of those."""
    names = list(callers)
    genomes = sorted(set().union(*(set(c) for c in callers.values())))
    rows, agree = [], 0
    for g in genomes:
        have = [c for c in names if g in callers[c]]
        allv = set().union(*(set(callers[c][g]) for c in have))
        for v in sorted(allv):
            by = [c for c in have if v in callers[c][g]]
            rows.append({"id": g, "pos": v[0], "ref": v[1], "alt": v[2], "by": by, "n": len(by), "of": len(have)})
            agree += len(by) == len(have)
    return {"callers": names, "rows": rows, "summary": {"genomes": len(genomes), "variants": len(rows), "all_agree": agree,
                                                         "share_all_agree": round(agree / len(rows), 4) if rows else None,
                                                         "only_one": sum(1 for r in rows if r["n"] == 1 and r["of"] > 1)}}


def sites_report(callers: dict, limit: int = 300) -> dict:
    """The agreement per site: the sites where the callers do not all agree (worst first), the pattern counts and the summary."""
    c = concordance(callers)
    by_site: dict[tuple, dict] = {}
    for r in c["rows"]:
        s = by_site.setdefault((r["pos"], r["ref"], r["alt"]), {"pos": r["pos"], "ref": r["ref"], "alt": r["alt"], "genomes": 0, "all": 0, **{t: 0 for t in callers}})
        s["genomes"] += 1
        s["all"] += r["n"] == r["of"]
        for t in r["by"]:
            s[t] += 1
    sites = sorted((s for s in by_site.values() if s["all"] < s["genomes"]), key=lambda s: (-(s["genomes"] - s["all"]), s["pos"]))
    patterns = Counter(",".join(r["by"]) for r in c["rows"])
    return {"genomes": len(set().union(*(set(v) for v in callers.values()))), "summary": c["summary"], "patterns": dict(patterns.most_common()), "sites_total": len(by_site),
            "sites_disagree": len(sites), "sites": sites[:limit], "callers": list(callers), "all_sites": sites}


def review_bed(sites: list, seqid: str, min_genomes: int = 3) -> str:
    """BED (0-based, end-exclusive) of the positions where at least `min_genomes` genomes get different calls from the callers - a list of places to look at, **not** a mask."""
    pos = sorted({s["pos"] for s in sites if s["genomes"] - s["all"] >= min_genomes})
    return "".join(f"{seqid}\t{p - 1}\t{p}\n" for p in pos)
