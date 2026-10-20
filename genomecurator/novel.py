"""Discovery of new variants in a dataset: groups of genomes that share mutations no designated lineage is defined by, and how fast such mutations keep appearing.

On the dataset's Nextstrain tree, a **candidate** is a monophyletic group of at least `min_size` genomes whose branch carries at least one mutation that is *not* a defining
mutation of any lineage in the lineage-designation definitions (positions are compared in NC_063383.1 coordinates, so the tree's reference does not matter). Candidates are
dated by the collection dates of their genomes and given a status: **new** (first genome within `window_days` of the newest date in the dataset), **growing** (at least half of
the genomes are that recent, or the last window has more genomes than the one before) or **steady**. Next to them: the **discovery curve** (how many new, non-defining clade
mutations appeared per month, and the running total) and, per such mutation, when it was first seen and how many genomes and countries carry it.

Limits: the tree is one IQ-TREE run without branch support, so a candidate is a hypothesis; a mutation seen once in a patient (private to one genome) is not listed (sequencing
errors and within-host changes cannot be told apart); sequencing effort differs by country and time, so *first seen* is not *emerged*; without the lineage definitions every
branch mutation counts as new (the page says so)."""
import re
from collections import Counter, defaultdict
from datetime import date, timedelta

from . import insight, lineages, log, mutprofile, nextstrain, outbreak

_log = log.get("novel")


def _to_date(s: str):
    try:
        y, m, *rest = s.split("-")
        return date(int(y), int(m), int(rest[0]) if rest else 15)
    except (ValueError, IndexError):
        return None


def _known_in_tree_coordinates(tree_fa) -> tuple[set, bool]:
    """(position in the tree's reference, new base) of every defining SNP of every designated lineage; False when there are no definitions."""
    if not lineages.status()["ok"]:
        return set(), False
    to_tree = lineages.defining_in_run_coordinates({}, tree_fa)               # NC_063383.1 -> the tree's reference (None when they are the same)
    out = set()
    for ob, d in lineages.definitions().items():
        for lin in d:
            for pos, nuc, _by in lineages.defining_snps(ob, lin):
                p = pos if to_tree is None else to_tree.get(pos)
                if p:
                    out.add((p, nuc))
    return out, True


def discover(run: dict, window_days: int = 120, min_size: int = 3, max_n: float | None = None, allow_homoplasy: bool = False) -> dict:
    tree, rootseq, nsrun = outbreak._tree(run)
    if not tree:
        return {"available": False, "reason": "이 데이터셋의 Nextstrain 분석(트리)이 필요합니다."}
    refname, refseq, cds = outbreak._tree_reference(nsrun)
    base = nextstrain.workflow_dir() / "defaults"
    d = base / (nsrun.get("build") or "")
    tree_fa = (d if (d / "reference.fasta").exists() else base) / "reference.fasta"
    known, have_defs = _known_in_tree_coordinates(tree_fa)
    meta = {g["id"]: g for g in lineages._genomes(run)}
    ok = mutprofile.good_samples(run, max_n)
    dates = {}
    for i, g in meta.items():
        kind, _w, _m, shown = lineages.parse_date(g["date"])
        if kind and (ok is None or i in ok):
            dates[i] = (_to_date(shown), shown)
    newest = max((v[0] for v in dates.values()), default=None)
    if not newest:
        return {"available": False, "reason": "수집일이 있는 유전체가 없습니다."}
    recent_from = newest - timedelta(days=window_days)
    prev_from = recent_from - timedelta(days=window_days)

    nodes: list[dict] = []
    branch_count: Counter = Counter()

    def walk(n, parent):
        rec = {"node": n, "parent": parent, "tips": [], "muts": outbreak._nuc_muts(n), "children": []}
        nodes.append(rec)
        for m in rec["muts"]:
            branch_count[m] += 1
        kids = n.get("children") or []
        if not kids:
            rec["tips"] = [n["name"]]
        for c in kids:
            walk(c, rec)
            rec["tips"] += tips_of[id(c)]
        tips_of[id(n)] = rec["tips"]
    tips_of: dict = {}
    import sys
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    walk(tree, None)
    L = len(refseq)

    def novel_muts(r):
        out = []
        for ref, pos, alt in r["muts"]:
            if mutprofile.in_itr(pos, L) or (pos, alt) in known:
                continue
            if branch_count[(ref, pos, alt)] > 1 and not allow_homoplasy:
                continue
            out.append((ref, pos, alt))
        return out

    cands = []
    for r in nodes:
        if r["parent"] is None or not r["node"].get("children"):
            continue
        tips = [t for t in r["tips"] if t in dates]
        if len(tips) < min_size:
            continue
        nm = novel_muts(r)
        if not nm:
            continue
        cands.append((r, tips, nm))
    # a chain of single-child nodes has the same genomes: keep one candidate and gather the mutations
    merged: dict[frozenset, dict] = {}
    for r, tips, nm in cands:
        key = frozenset(tips)
        if key in merged:
            merged[key]["muts"] += nm
        else:
            merged[key] = {"r": r, "tips": tips, "muts": list(nm)}
    items = list(merged.values())
    by_node = {id(it["r"]): it for it in items}

    def parent_candidate(it):
        p = it["r"]["parent"]
        while p is not None:
            if id(p) in by_node:
                return by_node[id(p)]
            p = p["parent"]
        return None

    out = []
    for it in items:
        tips = it["tips"]
        ds = sorted(dates[t][0] for t in tips)
        shown = sorted((dates[t][1] for t in tips), key=lambda s: (s[:7], s))
        n_recent = sum(1 for t in tips if dates[t][0] >= recent_from)
        n_prev = sum(1 for t in tips if prev_from <= dates[t][0] < recent_from)
        first = ds[0]
        if first >= recent_from:
            status = "new"
        elif n_recent / len(tips) >= 0.5 or (n_recent > n_prev and n_recent >= min_size):
            status = "growing"
        else:
            status = "steady"
        months = Counter(dates[t][0].strftime("%Y-%m") for t in tips)
        gs = [meta[t] for t in tips if t in meta]
        muts = []
        for ref, pos, alt in it["muts"]:
            an = mutprofile.annotate(pos, ref, alt, cds, refseq) if refseq else {"gene": "", "aa_change": "", "effect": ""}
            muts.append({"label": f"{ref}{pos}{alt}", "gene": an["gene"], "aa": an["aa_change"], "effect": an["effect"], "apobec": insight._apo_branch(ref, alt, pos, rootseq or refseq),
                         "branches": branch_count[(ref, pos, alt)]})
        pc = parent_candidate(it)
        out.append({"id": "", "n": len(tips), "n_recent": n_recent, "n_prev": n_prev, "status": status, "first": shown[0], "last": shown[-1], "first_ord": first.isoformat(),
                    "countries": dict(Counter(g["country"] for g in gs if g["country"]).most_common(8)), "n_countries": len({g["country"] for g in gs if g["country"]}),
                    "lineages": dict(Counter(g["lineage"] for g in gs if g["lineage"]).most_common(3)), "mutations": muts, "months": dict(sorted(months.items())),
                    "parent": id(pc["r"]) if pc else None, "_key": id(it["r"]), "apobec": sum(1 for m in muts if m["apobec"])})
    order = {"new": 0, "growing": 1, "steady": 2}
    out.sort(key=lambda c: (order[c["status"]], -c["n"], c["first_ord"]))
    names = {c["_key"]: f"V{i + 1}" for i, c in enumerate(out)}
    for c in out:
        c["id"] = names[c["_key"]]
        c["parent"] = names.get(c["parent"])
        del c["_key"]

    # discovery curve: a clade's novel mutations are first seen at the clade's earliest genome
    by_month: Counter = Counter()
    scatter = []
    seen = set()
    for r in nodes:
        if r["parent"] is None or not r["node"].get("children"):
            continue
        tips = [t for t in r["tips"] if t in dates]
        if len(tips) < min_size:
            continue
        for ref, pos, alt in novel_muts(r):
            if (ref, pos, alt) in seen:
                continue
            seen.add((ref, pos, alt))
            f = min(dates[t][0] for t in tips)
            by_month[f.strftime("%Y-%m")] += 1
            gs = [meta[t] for t in tips if t in meta]
            scatter.append({"label": f"{ref}{pos}{alt}", "first": f.isoformat(), "n": len(tips), "countries": len({g["country"] for g in gs if g["country"]}), "apobec": bool(insight._apo_branch(ref, alt, pos, rootseq or refseq)),
                            "recent": f >= recent_from})
    months = []
    if by_month:
        y, m = (int(x) for x in min(by_month).split("-"))
        ey, em = newest.year, newest.month
        total = 0
        while (y, m) <= (ey, em):
            k = f"{y}-{m:02d}"
            total += by_month.get(k, 0)
            months.append({"month": k, "new": by_month.get(k, 0), "total": total})
            y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    return {"available": True, "tree": nsrun["name"], "reference": refname, "newest": newest.isoformat(), "recent_from": recent_from.isoformat(), "window_days": window_days, "min_size": min_size,
            "definitions": have_defs, "max_n": max_n, "candidates": out, "counts": dict(Counter(c["status"] for c in out)), "curve": months, "scatter": scatter, "n_dated": len(dates),
            "n_novel_mutations": len(seen), "recent_novel_mutations": sum(1 for s in scatter if s["recent"])}
