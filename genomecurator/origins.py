"""How many times did each mutation arise independently? (evidence E7 of docs/ADAPTIVE.md)

From the dataset's Nextstrain tree: every branch that carries a mutation (position + new base) is one origin of it. An origin counts as **settled** when at least two genomes descend from that
branch; an origin with a single descendant genome cannot be told from a sequencing error or a change inside one patient, so it is counted apart. A mutation that arose several times on
unrelated branches and spread each time is a stronger sign of being kept than one that arose once and was sampled many times.

The tree comes from augur's ancestral reconstruction, so the counts are only as good as the tree and its placement of mutations; a recombinant, a long-branch artefact or a sequencing error can
make a mutation look recurrent. The positions are those of the tree's reference: the counts are returned only when that is the run's own reference (otherwise `None`, and the caller keeps the
lineage-based approximation)."""
import json
import math
import sys
from collections import Counter

from . import continents, insight, lineages, mutprofile, nextstrain, outbreak

_CACHE: dict = {}


def independent_origins(run: dict, min_tips: int = 2) -> dict | None:
    """{(position, new base): {events, settled, singles, tips, lineages, countries}} or None when the run has no usable tree (no Nextstrain run for the dataset, or a different reference)."""
    try:
        tree, _root, ns = outbreak._tree(run)
    except (KeyError, OSError, ValueError):                   # a run without a dataset / a tree file that cannot be read: no counts, the caller falls back
        return None
    if not tree:
        return None
    refname, refseq, cds = outbreak._tree_reference(ns)
    run_ref = mutprofile.summary(run).get("reference", "")
    if not refname or refname != run_ref:
        return None
    f = nextstrain.auspice_dir() / ns["auspice_file"]
    key = (run["id"], str(f), f.stat().st_mtime, min_tips)
    if key in _CACHE:
        return _CACHE[key]
    meta = {g["id"]: g for g in lineages._genomes(run)}
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    out: dict = {}

    def walk(node):
        kids = node.get("children") or []
        tips = [node["name"]] if not kids else [t for c in kids for t in walk(c)]
        for _ref, pos, alt in outbreak._nuc_muts(node):
            e = out.setdefault((pos, alt), {"events": 0, "settled": 0, "singles": 0, "tips": 0, "lineages": Counter(), "countries": Counter()})
            e["events"] += 1
            if len(tips) >= min_tips:
                e["settled"] += 1
                e["tips"] += len(tips)
                for t in tips:
                    g = meta.get(t)
                    if g:
                        e["lineages"][g["lineage"] or "(없음)"] += 1
                        if g["country"]:
                            e["countries"][g["country"]] += 1
            else:
                e["singles"] += 1
        return tips
    walk(tree)
    for e in out.values():
        e["lineages"], e["countries"] = dict(e["lineages"]), dict(e["countries"])
    neutral_baseline(out, refseq, cds)
    _CACHE[key] = out
    return out


def mutation_class(pos: int, alt: str, refseq: str) -> str:
    """APOBEC3-type (TC>TT, GA>AA), else transition / transversion: the classes a neutral mutation is compared within."""
    ref = refseq[pos - 1] if 0 < pos <= len(refseq) else "?"
    if insight._apo_branch(ref, alt, pos, refseq):
        return "apobec"
    return "transition" if {ref, alt} in ({"A", "G"}, {"C", "T"}) else "transversion"


def _count_sf(k: int, mean: float, var: float) -> float:
    """P(X >= k) for a count with the given mean and variance: Poisson, or negative binomial (method of moments) when the counts are overdispersed - recurrence clusters on a tree, so they usually are."""
    if k <= 0:
        return 1.0
    if mean <= 0:
        return 0.0 if k > 0 else 1.0
    if var <= mean * 1.0001:
        cdf = sum(math.exp(-mean + i * math.log(mean) - math.lgamma(i + 1)) for i in range(k))
    else:
        r = mean * mean / (var - mean)
        pr = r / (r + mean)
        cdf = sum(math.exp(math.lgamma(i + r) - math.lgamma(r) - math.lgamma(i + 1) + r * math.log(pr) + i * math.log(1 - pr)) for i in range(k))
    return max(0.0, min(1.0, 1.0 - cdf))


def neutral_baseline(o: dict, refseq: str, cds: list) -> None:
    """Evidence E8 of docs/ADAPTIVE.md, in place. Synonymous and non-coding mutations are taken as neutral; within the same mutation class their settled-origin counts give the chance distribution of
    "k independent origins" (a Poisson, or a negative binomial when overdispersed, fitted by moments to the neutral counts of the class). For a protein-changing mutation p = P(neutral count >= its settled
    origins), BH-corrected to q over all protein-changing mutations of the tree. Adds `effect`, `class`, `p_neutral`, `q_neutral`, `neutral_n` to every entry (p / q only for protein-changing ones).
    Caveats: the tail is a fitted model, extrapolated beyond what the neutral mutations show (an empirical tail would need thousands of neutral mutations to reach q < 0.05); the null is the recurrence among
    mutations seen at all, ignores how many genomes carry a mutation (more carriers allow more origins) and assumes synonymous sites are neutral."""
    neutral: dict[str, list[int]] = {}
    for (pos, alt), e in o.items():
        ref = refseq[pos - 1] if 0 < pos <= len(refseq) else "?"
        an = mutprofile.annotate(pos, ref, alt, cds, refseq) if ref in "ACGT" else {"effect": ""}
        e["effect"], e["class"] = an["effect"], mutation_class(pos, alt, refseq)
        if e["effect"] in ("synonymous", "intergenic"):
            neutral.setdefault(e["class"], []).append(e["settled"])
    moments = {}
    for c, v in neutral.items():
        m = sum(v) / len(v)
        moments[c] = (m, sum((x - m) ** 2 for x in v) / max(1, len(v) - 1))
    test = [(k, e) for k, e in o.items() if e["effect"] in ("nonsynonymous", "stop_gained")]
    for e in o.values():
        e["p_neutral"] = e["q_neutral"] = None
        e["neutral_n"] = len(neutral.get(e["class"], []))
    ps = []
    for k, e in test:
        m, v = moments.get(e["class"], (0.0, 0.0))
        e["p_neutral"] = _count_sf(e["settled"], m, v) if e["neutral_n"] >= 20 else None        # too few neutral mutations of this class to say
        ps.append(1.0 if e["p_neutral"] is None else e["p_neutral"])
    for (k, e), q in zip(test, continents.bh(ps)):
        e["q_neutral"] = None if e["p_neutral"] is None else q
