"""A Nextstrain `clades.tsv` for the lineages of clade Ib, from the mpxv-lineages definitions.

Nextstrain's mpox build labels the clade I tree with only clade I, Ia, sh2024 and Ib (one site for Ib), so the lineages of the Ib outbreak (sh2023: A, A.1 ... A.4.1 ...) do not show on it, whereas the
hMPXV-1 build has its lineages in `clades.tsv`. This writes the missing part in the same format: a lineage is a clade whose parent is its parent lineage (a row `<lineage> clade <parent>`), followed by the
sites that *introduced* it (`<lineage> nuc <position> <base>`), in the coordinates of the build's reference (KJ642613.1 for the clade I build); the definitions are in NC_063383.1 coordinates and are
moved with a minimap2 liftover. `augur clades` then gives a node the deepest lineage whose sites it all carries (and whose parent it carries).

`validate` replays that rule on the mutation calls of a run and compares it with the lineage the genome table holds (the lineage Pathoplexus / Nextclade assigned), so a wrong coordinate or a definition
that the data do not support shows before the build is run."""
from collections import Counter, defaultdict
from pathlib import Path

from . import config, lineages, mutprofile, nextstrain

UPSTREAM = "defaults/clade-i/clades.tsv"
OUTBREAK = "sh2023"


def _order(defs: dict) -> list[str]:
    """Lineages with parents first."""
    out, seen = [], set()

    def visit(name):
        if name in seen or name not in defs:
            return
        seen.add(name)
        par = defs[name]["parent"] if defs[name]["parent"] not in ("None", "none", "") else ""
        if par:
            visit(par)
        out.append(name)
    for n in sorted(defs):
        if n:
            visit(n)
    return out


def own_sites(outbreak: str = OUTBREAK) -> dict[str, dict]:
    """{lineage: {"parent": parent lineage or the outbreak, "sites": [(NC position, base)] introduced by it}} for one outbreak, plus the outbreak's own root sites under the key ''."""
    defs = lineages.definitions().get(outbreak, {})
    out = {}
    root = defs.get("", {})
    out[""] = {"parent": "", "sites": [(p, n) for p, n in root.get("snps", [])]}
    for name in _order(defs):
        par = defs[name]["parent"]
        par = "" if par in ("None", "none", "") or par not in defs else par
        intro = [(p, n) for p, n, by in lineages.defining_snps(outbreak, name) if by == name]
        out[name] = {"parent": par, "sites": intro}
    return out


def to_run_coordinates(run: dict) -> dict[int, int] | None:
    ref = mutprofile.summary(run).get("reference", "")
    fa, _gff = mutprofile.reference_by_name(ref)
    return lineages.defining_in_run_coordinates(run, fa)


def build(run: dict, outbreak: str = OUTBREAK, base_clades: Path | None = None) -> dict:
    """-> {"text": the clades.tsv, "lineages": n, "sites": n, "unmapped": [lineage: NC position]}: the build's own file (clade I, Ia, sh2024, Ib) followed by the outbreak's lineages."""
    base = base_clades or (nextstrain.workflow_dir() / "defaults" / "clade-i" / "clades.tsv")
    text = base.read_text(encoding="utf-8").rstrip("\n") + "\n"
    cmap = to_run_coordinates(run)
    sites = own_sites(outbreak)
    ref = mutprofile.summary(run).get("reference", "")
    lines = ["", f"# Lineages of {outbreak} (clade Ib), from mpxv-lineages {lineages.commit()}; sites moved from NC_063383.1 to {ref} with a minimap2 liftover", f"{outbreak}\tclade\tclade Ib"]
    unmapped, n_sites = [], 0

    def rows(name, label):
        nonlocal n_sites
        for p, nuc in sorted(sites[name]["sites"]):
            q = cmap.get(p) if cmap is not None else p
            if not q:
                unmapped.append((label, p))
                continue
            lines.append(f"{label}\tnuc\t{q}\t{nuc}")
            n_sites += 1
    rows("", outbreak)
    for name in [n for n in sites if n]:
        lines += ["", f"{name}\tclade\t{sites[name]['parent'] or outbreak}"]
        rows(name, name)
    return {"text": text + "\n".join(lines) + "\n", "lineages": len(sites) - 1, "sites": n_sites, "unmapped": unmapped}


def validate(run: dict, outbreak: str = OUTBREAK) -> dict:
    """Replay `augur clades` on the SNV calls of a run: a genome gets the deepest lineage whose own sites it all carries with its parent also carried. Compare with the lineage in the genome table."""
    cmap = to_run_coordinates(run)
    sites = own_sites(outbreak)
    kids = defaultdict(list)
    for n, d in sites.items():
        if n:
            kids[d["parent"]].append(n)
    carried: dict[str, set] = defaultdict(set)
    snv_at: dict[str, set] = defaultdict(set)
    for r in mutprofile._load_mutations(run):
        carried[r["sample"]].add((int(r["pos"]), r["alt"].upper()))
        snv_at[r["sample"]].add(int(r["pos"]))
    meta = mutprofile._run_meta(run)
    fa, _gff = mutprofile.reference_by_name(mutprofile.summary(run).get("reference", ""))
    refseq = mutprofile.read_fasta(fa)[1].upper()

    unknown = lineages.nc_unknown(run)                                           # per genome, the NC_063383.1 intervals that are N, deleted or not covered
    win = mutprofile.run_itr(run)                                                # the mutation table leaves the ITR out: a site there cannot be read from it

    def need(name):
        return [(p, (cmap.get(p) if cmap is not None else p), n) for p, n in sites[name]["sites"] if (cmap.get(p) if cmap is not None else p)]

    def state(sample, site):
        """True / False whether the genome carries base n at the site, None when that cannot be told (N, deleted, outside the window): such a site does not count for or against."""
        p, q, n = site
        if (win and not (win[0] <= q <= win[1])) or lineages.is_unknown(unknown.get(sample, []), p):
            return None
        return (q, n) in carried.get(sample, ()) or (q not in snv_at.get(sample, ()) and 0 < q <= len(refseq) and refseq[q - 1] == n)

    def ok(sample, name):
        return all(state(sample, x) is not False for x in need(name))               # an unreadable site does not count against

    def known(sample, name):
        return sum(1 for x in need(name) if state(sample, x))

    def assign(sample):
        if not ok(sample, ""):
            return ""
        cur = ""
        while True:
            nxt = [k for k in kids[cur] if ok(sample, k) and (known(sample, k) > 0 or not need(k))]   # a lineage with sites needs one readable (an all-N genome is not put in it); one with no sites of its own (A) follows its parent
            if not nxt:
                return cur
            cur = max(nxt, key=lambda k: known(sample, k))
    res = Counter()
    bad = []
    for s, m in meta.items():
        if m.get("outbreak") not in (outbreak, "") or not m["lineage"]:
            continue
        got = assign(s)
        res[(m["lineage"], got)] += 1
        if got != m["lineage"] and len(bad) < 30:
            bad.append({"genome": s, "table": m["lineage"], "replayed": got or "(sh2023 root)"})
    total = sum(res.values())
    same = sum(v for (a, b), v in res.items() if a == b)
    return {"genomes": total, "same": same, "share": same / total if total else None, "pairs": {f"{a} -> {b or '(root)'}": v for (a, b), v in res.most_common() if a != b}, "examples": bad}


def fix_auspice(path: Path, outbreak: str = OUTBREAK, clade: str = "Ib") -> dict:
    """Make an Auspice file of a build run with the lineage `clades.tsv` read correctly. The workflow's `clades_renaming.py` knows only the hMPXV-1 lineage names: any name that is not clade I / Ia /
    sh2024 / outgroup is taken as an `sh2017` lineage of clade IIb, so the Ib lineages come out as `IIb` (clade_membership) with the branch label `sh2017/<lineage>`. This reads the branch labels (the
    lineage each branch starts), gives every node the lineage of its nearest labelled ancestor, and sets for the nodes of the outbreak's lineages `clade_membership` = Ib, `outbreak` = sh2023, `lineage` =
    the lineage, with a `lineage` colouring and the branch label `sh2023/<lineage>`. The original file is kept as `<name>.orig.json`. -> {"nodes": changed nodes, "tips": tips, "lineages": {lineage: tips}}"""
    import json
    import sys
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    names = {n for n in own_sites(outbreak) if n}
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    stats, tips = Counter(), [0, 0]

    def walk(node, cur):
        lab = (node.get("branch_attrs", {}).get("labels") or {}).get("clade")
        name = (lab or "").split("/")[-1]
        if name in names:
            cur = name
            node.setdefault("branch_attrs", {}).setdefault("labels", {})["clade"] = f"{outbreak}/{name}"
        elif lab:
            cur = None if lab in ("Ib", "IIb", "unassigned") else cur
        kids = node.get("children") or []
        if cur:
            na = node.setdefault("node_attrs", {})
            na["clade_membership"] = {"value": clade}
            na["outbreak"] = {"value": outbreak}
            na["lineage"] = {"value": cur}
            tips[0] += 1
            if not kids:
                stats[cur] += 1
                tips[1] += 1
        for c in kids:
            walk(c, cur)
    walk(data["tree"], None)
    cols = data.setdefault("meta", {}).setdefault("colorings", [])
    if not any(c.get("key") == "lineage" for c in cols):
        cols.insert(1, {"key": "lineage", "title": "Lineage", "type": "categorical"})
    orig = path.with_name(path.stem + ".orig.json")
    if not orig.exists():
        orig.write_bytes(path.read_bytes())
    path.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    return {"nodes": tips[0], "tips": tips[1], "lineages": dict(stats)}
