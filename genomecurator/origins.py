"""How many times did each mutation arise independently? (evidence E7 of docs/ADAPTIVE.md)

From the dataset's Nextstrain tree: every branch that carries a mutation (position + new base) is one origin of it. An origin counts as **settled** when at least two genomes descend from that
branch; an origin with a single descendant genome cannot be told from a sequencing error or a change inside one patient, so it is counted apart. A mutation that arose several times on
unrelated branches and spread each time is a stronger sign of being kept than one that arose once and was sampled many times.

The tree comes from augur's ancestral reconstruction, so the counts are only as good as the tree and its placement of mutations; a recombinant, a long-branch artefact or a sequencing error can
make a mutation look recurrent. The positions are those of the tree's reference: the counts are returned only when that is the run's own reference (otherwise `None`, and the caller keeps the
lineage-based approximation)."""
import json
import sys
from collections import Counter

from . import lineages, mutprofile, nextstrain, outbreak

_CACHE: dict = {}


def independent_origins(run: dict, min_tips: int = 2) -> dict | None:
    """{(position, new base): {events, settled, singles, tips, lineages, countries}} or None when the run has no usable tree (no Nextstrain run for the dataset, or a different reference)."""
    try:
        tree, _root, ns = outbreak._tree(run)
    except (KeyError, OSError, ValueError):                   # a run without a dataset / a tree file that cannot be read: no counts, the caller falls back
        return None
    if not tree:
        return None
    refname, _seq, _cds = outbreak._tree_reference(ns)
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
    _CACHE[key] = out
    return out
