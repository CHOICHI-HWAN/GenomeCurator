"""Outbreak analysis of a mutation-profile run, after the approach of Lee, Choi et al., *Unraveling the transmission dynamics and evolutionary trajectory of mpox in South Korea,
2022-2025* (Emerg Microbes Infect 2026, doi 10.1080/22221751.2026.2746255): the same questions on any dataset.

1. **Genomic subclusters** from the tree topology: monophyletic groups of at least `min_size` genomes whose branch carries at least one defining mutation (no SNP-distance
   threshold), with the branch-defining mutations (gene, amino-acid change, APOBEC3-type, homoplasy = the same mutation on other branches too), a cluster x mutation matrix
   (the paper's SNIPIT panel), first / last date and countries, and a bubble chart over time. Needs the dataset's Nextstrain tree (the newest finished run).
2. **Cluster size distribution**: Poisson against negative binomial (variance-to-mean ratio, Pearson dispersion, theta, AIC / BIC).
3. **APOBEC3-type mutations accumulate over time**: per genome the number of APOBEC3-type mutations (against NC_063383.1) by collection date, per year, with a trend.
4. **Mutation types per genome** (APOBEC3-type / other x nonsynonymous / synonymous), in chronological order and by cluster.
5. **Genes**: mutation burden per gene (hotspots above a threshold, ITRs excluded), APOBEC3 share, functional category (the paper's Supplementary Table 1, `resources/opg_categories.tsv`)
   and an *exploratory* dN/dS per gene (Nei-Gojobori, distinct mutations of the dataset against the run's reference) - not a test of selection.

The per-genome counts are against the run's reference; sites that almost every genome has (differences between the reference and the whole dataset) can be left out (`fixed`)."""
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from . import insight, lineages, log, mutprofile, nextstrain

_log = log.get("outbreak")
ITR = 6500                                                  # bp at each end counted as inverted terminal repeat (hotspots there are left out)
_CAT: dict = {}


def categories() -> dict[str, dict]:
    if not _CAT:
        f = Path(__file__).parent / "resources" / "opg_categories.tsv"
        if f.exists():
            for r in csv.DictReader((l for l in open(f, encoding="utf-8") if not l.startswith("#")), delimiter="\t"):
                _CAT.setdefault(r["OPG"], r)
    return _CAT


# ------------------------------------------------------------------------------------------------ tree
def _tree(run: dict):
    runs = [r for r in nextstrain.list_runs(run["dataset_id"]) if r["tool"] == "nextstrain" and r["status"] == "done" and r["auspice_file"]]
    if not runs:
        return None, None, None
    f = nextstrain.auspice_dir() / runs[0]["auspice_file"]
    rootf = nextstrain.auspice_dir() / runs[0]["auspice_file"].replace(".json", "_root-sequence.json")
    if not f.exists():
        return None, None, None
    root = json.loads(rootf.read_text(encoding="utf-8")).get("nuc", "") if rootf.exists() else ""
    return json.loads(f.read_text(encoding="utf-8"))["tree"], root, runs[0]


def _tree_reference(nsrun: dict):
    base = nextstrain.workflow_dir() / "defaults"
    d = base / (nsrun.get("build") or "")
    fa, gff = (d / "reference.fasta", d / "genome_annotation.gff3") if (d / "reference.fasta").exists() else (base / "reference.fasta", base / "genome_annotation.gff3")
    if not fa.exists():
        return "", "", []
    name, seq = mutprofile.read_fasta(fa)
    return name, seq, mutprofile.load_cds(gff if gff.exists() else None)


def _nuc_muts(node: dict) -> list[tuple[str, int, str]]:
    out = []
    for m in (node.get("branch_attrs", {}).get("mutations", {}) or {}).get("nuc", []):
        g = re.fullmatch(r"([ACGT])(\d+)([ACGT])", m)
        if g:
            out.append((g[1], int(g[2]), g[3]))
    return out


# ------------------------------------------------------------------------------------------------ statistics
def _lgamma(x: float) -> float:
    return math.lgamma(x)


def pois_ll(y: list[int], mu: float) -> float:
    return sum(v * math.log(mu) - mu - _lgamma(v + 1) for v in y) if mu > 0 else float("-inf")


def nb_ll(y: list[int], mu: float, theta: float) -> float:
    return sum(_lgamma(v + theta) - _lgamma(theta) - _lgamma(v + 1) + theta * math.log(theta / (theta + mu)) + v * math.log(mu / (theta + mu)) for v in y)


def size_distribution(sizes: list[int]) -> dict | None:
    """Poisson against negative binomial (intercept only: mu = the mean) for the cluster sizes: VMR, Pearson dispersion, theta (maximum likelihood, golden-section; its standard error is
    approximate, from the curvature of the log-likelihood in theta at fixed mu), AIC and BIC."""
    n = len(sizes)
    if n < 3 or sum(sizes) == 0:
        return None
    mu = sum(sizes) / n
    var = sum((v - mu) ** 2 for v in sizes) / (n - 1)
    pearson = sum((v - mu) ** 2 / mu for v in sizes) / (n - 1)
    f = lambda t: nb_ll(sizes, mu, math.exp(t))
    a, b = math.log(1e-3), math.log(1e4)
    g = (math.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    for _ in range(120):
        if f(c) > f(d):
            b = d
        else:
            a = c
        c, d = b - g * (b - a), a + g * (b - a)
    theta = math.exp((a + b) / 2)
    h = theta * 1e-3
    d2 = (nb_ll(sizes, mu, theta + h) - 2 * nb_ll(sizes, mu, theta) + nb_ll(sizes, mu, theta - h)) / h ** 2
    se = math.sqrt(-1 / d2) if d2 < 0 else None
    lp, ln_ = pois_ll(sizes, mu), nb_ll(sizes, mu, theta)
    return {"n": n, "mean": round(mu, 2), "variance": round(var, 2), "vmr": round(var / mu, 2), "pearson": round(pearson, 2), "theta": round(theta, 3), "theta_se": round(se, 3) if se else None,
            "aic_poisson": round(2 * 1 - 2 * lp, 2), "aic_nb": round(2 * 2 - 2 * ln_, 2), "bic_poisson": round(math.log(n) * 1 - 2 * lp, 2), "bic_nb": round(math.log(n) * 2 - 2 * ln_, 2),
            "overdispersed": var / mu > 1.5, "nb_better": (2 * 2 - 2 * ln_) < (2 * 1 - 2 * lp)}


def _rank(v: list[float]) -> list[float]:
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def trend(x: list[float], y: list[float]) -> dict | None:
    """Least-squares slope of y on x (per unit of x) and Spearman's rho (p from the Fisher z approximation)."""
    n = len(x)
    if n < 5:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx if sxx else 0.0
    rx, ry = _rank(x), _rank(y)
    mrx, mry = sum(rx) / n, sum(ry) / n
    num = sum((a - mrx) * (b - mry) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mrx) ** 2 for a in rx) * sum((b - mry) ** 2 for b in ry))
    rho = num / den if den else 0.0
    z = math.atanh(max(-0.999999, min(0.999999, rho))) * math.sqrt(max(n - 3, 1))
    return {"n": n, "slope_per_year": round(slope, 3), "rho": round(rho, 3), "p": math.erfc(abs(z) / math.sqrt(2))}


# ------------------------------------------------------------------------------------------------ 1-2: subclusters
def clusters(run: dict, min_size: int = 5, min_muts: int = 1, tips: bool = False) -> dict:
    tree, rootseq, nsrun = _tree(run)
    if not tree:
        return {"available": False, "reason": "이 데이터셋의 Nextstrain 분석(트리)이 필요합니다."}
    refname, refseq, cds = _tree_reference(nsrun)
    meta = {g["id"]: g for g in lineages._genomes(run)}
    # collect: tips below every node, mutations on every branch, how many branches carry the same mutation (homoplasy)
    nodes: list[dict] = []
    branch_count: Counter = Counter()

    def walk(n, parent):
        rec = {"node": n, "parent": parent, "tips": [], "muts": _nuc_muts(n), "children": []}
        nodes.append(rec)
        if parent is not None:
            parent["children"].append(rec)
        for m in rec["muts"]:
            branch_count[m] += 1
        kids = n.get("children") or []
        if not kids:
            rec["tips"] = [n["name"]]
        for c in kids:
            walk(c, rec)
            rec["tips"] += nodes_tips[id(c)]
        nodes_tips[id(n)] = rec["tips"]
    nodes_tips: dict = {}
    import sys
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20000))
    walk(tree, None)
    qual = [r for r in nodes if r["parent"] is not None and len(r["tips"]) >= min_size and len(r["muts"]) >= min_muts and r["node"].get("children")]
    qset = {id(r) for r in qual}
    minimal = [r for r in qual if not any(id(c) in qset or any(id(g) in qset for g in _desc(c)) for c in r["children"])]

    def cum(r):
        out, cur = [], r
        while cur is not None:
            out = cur["muts"] + out
            cur = cur["parent"]
        return out

    def describe(r, cid):
        tips = r["tips"]
        gs = [meta[t] for t in tips if t in meta]
        dated = sorted(((lineages.parse_date(g["date"])[3] or "9999"), g["id"]) for g in gs if lineages.parse_date(g["date"])[3])
        muts = []
        for ref, pos, alt in r["muts"]:
            an = mutprofile.annotate(pos, ref, alt, cds, refseq) if refseq else {"gene": "", "aa_change": "", "effect": ""}
            muts.append({"label": f"{ref}{pos}{alt}", "gene": an["gene"], "aa": an["aa_change"], "effect": an["effect"], "apobec": insight._apo_branch(ref, alt, pos, rootseq) if rootseq else False,
                         "branches": branch_count[(ref, pos, alt)], "homoplasy": branch_count[(ref, pos, alt)] > 1})
        months = Counter(lineages.parse_date(g["date"])[2] for g in gs if lineages.parse_date(g["date"])[2])
        return {"id": cid, "n": len(tips), "mutations": muts, "countries": dict(Counter(g["country"] for g in gs if g["country"]).most_common(8)),
                "lineages": dict(Counter(g["lineage"] for g in gs if g["lineage"]).most_common(4)), "first": dated[0][0] if dated else "", "last": dated[-1][0] if dated else "",
                "months": dict(sorted(months.items())), "tips": tips, "path": [f"{a}{b}{c}" for a, b, c in cum(r)]}

    def first_date(r):
        ds = [lineages.parse_date(meta[t]["date"])[3] for t in r["tips"] if t in meta and lineages.parse_date(meta[t]["date"])[3]]
        return min(ds) if ds else "9999"
    minimal.sort(key=first_date)
    cl = [describe(r, f"C{i + 1}") for i, r in enumerate(minimal)]
    inside = {t for c in cl for t in c["tips"]}
    sizes = [c["n"] for c in cl]
    cols = []
    for c in cl:
        for m in c["mutations"]:
            if m["label"] not in cols and not m["homoplasy"]:
                cols.append(m["label"])
    cols.sort(key=lambda s: int(re.search(r"\d+", s)[0]))
    matrix = [{"id": c["id"], "cells": {lab: (lab[-1] if lab in set(c["path"]) else "") for lab in cols}} for c in cl]
    ann = {m["label"]: m for c in cl for m in c["mutations"]}
    return {"available": True, "tree": nsrun["name"], "reference": refname, "min_size": min_size, "min_muts": min_muts, "tips": len(nodes_tips[id(tree)]), "clusters": [c if tips else {k: v for k, v in c.items() if k != "tips"} for c in cl],
            "unassigned": len(nodes_tips[id(tree)]) - len(inside), "size_distribution": size_distribution(sizes), "columns": [{"label": lab, "gene": ann[lab]["gene"], "aa": ann[lab]["aa"], "apobec": ann[lab]["apobec"]} for lab in cols],
            "matrix": matrix, "months": sorted({m for c in cl for m in c["months"]}), "n_nested": len(qual)}


def _desc(r):
    out, stack = [], list(r["children"])
    while stack:
        x = stack.pop()
        out.append(x)
        stack.extend(x["children"])
    return out


def within_lineage_trend(rows: list[dict], key: str = "apobec_var", min_n: int = 5) -> dict | None:
    """The trend of `key` over time **inside lineages**: both the date and the count are taken relative to the lineage's own mean, so a rise that only reflects later samples belonging to
    lineages that already carry more inherited mutations (phylogenetic depth) does not count. Lineages with fewer than `min_n` dated genomes are left out."""
    by = defaultdict(list)
    for r in rows:
        by[r["lineage"] or "-"].append(r)
    x, y, used = [], [], 0
    for lin, v in by.items():
        if len(v) < min_n:
            continue
        used += 1
        mt, my = sum(r["t"] for r in v) / len(v), sum(r[key] for r in v) / len(v)
        x += [r["t"] - mt for r in v]
        y += [r[key] - my for r in v]
    if len(x) < 10:
        return None
    t = trend(x, y)
    if t:
        t["lineages"] = used
    return t


# ------------------------------------------------------------------------------------------------ 3-4: per genome
def _rows(run: dict, max_n: float | None = None) -> list[dict]:
    return insight._mutations(run, max_n)


def _fixed_sites(run: dict, share: float, max_n: float | None = None) -> set:
    by = defaultdict(set)
    rows = _rows(run, max_n)
    for r in rows:
        by[(int(r["pos"]), r["alt"])].add(r["sample"])
    n = len({r["sample"] for r in rows}) or 1
    return {k for k, s in by.items() if len(s) / n >= share}


def accumulation(run: dict, fixed: float = 0.9, max_n: float | None = None) -> dict:
    """APOBEC3-type mutations per genome against NC_063383.1 over collection date (own rule: G>A followed by A, C>T preceded by T), all and without the sites nearly every genome has."""
    ncname, ncseq, _cds = lineages.nc_reference()
    nc = lineages._nc_snps(run)
    ok = mutprofile.good_samples(run, max_n)
    if ok is not None:
        nc = {i: v for i, v in nc.items() if i in ok}
    gs = {g["id"]: g for g in lineages._genomes(run)}
    allc = Counter(pos for sn in nc.values() for pos in sn)
    n = len(nc) or 1
    rows = []
    for i, sn in nc.items():
        g = gs.get(i)
        if not g:
            continue
        kind, wk, mo, shown = lineages.parse_date(g["date"])
        if kind != "day":
            continue
        y, m, d = (int(x) for x in shown.split("-"))
        from datetime import date
        t = date(y, m, d)
        dec = y + (t.timetuple().tm_yday - 1) / 365.25
        a_all = a_var = n_all = n_var = 0
        for pos, alt in sn.items():
            ref = ncseq[pos - 1] if 0 < pos <= len(ncseq) else ""
            if pos <= ITR or pos > len(ncseq) - ITR or ref not in "ACGT":
                continue
            apo = insight._apo_branch(ref, alt, pos, ncseq)
            n_all += 1
            a_all += apo
            if allc[pos] / n < fixed:
                n_var += 1
                a_var += apo
        rows.append({"id": i, "date": shown, "year": y, "t": dec, "apobec_all": a_all, "snps_all": n_all, "apobec_var": a_var, "snps_var": n_var, "lineage": g["lineage"]})
    rows.sort(key=lambda r: r["t"])
    years = {}
    for r in rows:
        years.setdefault(r["year"], []).append(r)

    def q(v, p):
        v = sorted(v)
        k = (len(v) - 1) * p
        lo, hi = math.floor(k), math.ceil(k)
        return v[lo] + (v[hi] - v[lo]) * (k - lo)
    by_year = [{"year": y, "n": len(v), "all_median": q([r["apobec_all"] for r in v], .5), "all_q1": q([r["apobec_all"] for r in v], .25), "all_q3": q([r["apobec_all"] for r in v], .75),
                "var_median": q([r["apobec_var"] for r in v], .5), "var_q1": q([r["apobec_var"] for r in v], .25), "var_q3": q([r["apobec_var"] for r in v], .75)} for y, v in sorted(years.items())]
    t = [r["t"] for r in rows]
    return {"reference": ncname, "fixed_share": fixed, "max_n": max_n, "n": len(rows), "by_year": by_year, "trend_all": trend(t, [r["apobec_all"] for r in rows]), "trend_variable": trend(t, [r["apobec_var"] for r in rows]),
            "trend_within_lineage": within_lineage_trend(rows),
            "points": [{k: r[k] for k in ("id", "date", "t", "apobec_all", "apobec_var", "snps_all", "snps_var", "lineage")} for r in rows]}


def spectrum(run: dict, fixed: float = 0.9, max_n: float | None = None) -> dict:
    """Per genome (chronological): APOBEC3-type / other x nonsynonymous / synonymous (stop gained counts as nonsynonymous; intergenic apart), against the run's reference."""
    gs = {g["id"]: g for g in lineages._genomes(run)}
    skip = _fixed_sites(run, fixed, max_n) if fixed < 1 else set()
    per: dict[str, Counter] = defaultdict(Counter)
    for r in _rows(run, max_n):
        if (int(r["pos"]), r["alt"]) in skip:
            continue
        a = "a3" if insight._is_apo(r) else "other"
        e = r["effect"]
        per[r["sample"]][f"{a}_{'syn' if e == 'synonymous' else 'inter' if e == 'intergenic' else 'nonsyn'}"] += 1
    rows = []
    for i, c in per.items():
        g = gs.get(i, {})
        kind, _w, _m, shown = lineages.parse_date(g.get("date", ""))
        rows.append({"id": i, "date": shown if kind else "", "lineage": g.get("lineage", ""), "a3_nonsyn": c["a3_nonsyn"], "a3_syn": c["a3_syn"], "other_nonsyn": c["other_nonsyn"], "other_syn": c["other_syn"],
                     "a3_inter": c["a3_inter"], "other_inter": c["other_inter"]})
    rows.sort(key=lambda r: (r["date"] or "9999", r["id"]))
    return {"fixed_share": fixed, "rows": rows}


# ------------------------------------------------------------------------------------------------ 5: genes
def _sites(codon: str) -> float:
    """Synonymous sites of a codon (Nei-Gojobori): the share of the 9 single changes that keep the amino acid (changes to / from a stop count as nonsynonymous)."""
    a = mutprofile.CODON.get(codon)
    if a is None:
        return 0.0
    s = 0.0
    for i in range(3):
        k = 0
        for b in "ACGT":
            if b != codon[i]:
                k += mutprofile.CODON.get(codon[:i] + b + codon[i + 1:]) == a and a != "*"
        s += k / 3
    return s


def _pois_sf(k: int, mu: float) -> float:
    """P(X >= k) for X ~ Poisson(mu)."""
    if k <= 0:
        return 1.0
    term = math.exp(-mu)
    cdf = 0.0
    for i in range(k):
        cdf += term
        term *= mu / (i + 1)
    return max(0.0, 1.0 - cdf)


def independent_by_gene(run: dict) -> dict | None:
    """Mutation **events** on the tree's branches (each event once, however many genomes inherit it) placed in the genes of the tree's reference, with the ITRs left out. A gene is a
    hotspot when it has more events than its length predicts at the genome-wide rate (Poisson, Bonferroni over the genes), not merely more mutations summed over genomes - that sum
    counts one ancestral mutation once per descendant and mostly shows which genes carry the mutations that define big clades."""
    tree, rootseq, nsrun = _tree(run)
    if not tree:
        return None
    refname, refseq, cds = _tree_reference(nsrun)
    b = insight._tree_branches(run)
    if not b or not refseq:
        return None
    _root, L, muts = b
    by = {c["gene"]: c for c in cds}
    glen = {g: sum(e - s + 1 for s, e in c["exons"]) for g, c in by.items()}
    cnt: dict[str, Counter] = defaultdict(Counter)
    sites: dict[str, set] = defaultdict(set)
    for internal, (ref, pos, alt) in muts:
        an = mutprofile.annotate(pos, ref, alt, cds, refseq)
        g = an["gene"] or "(intergenic)"
        c = cnt[g]
        c["events"] += 1
        c["internal" if internal else "terminal"] += 1
        c["a3"] += insight._apo_branch(ref, alt, pos, rootseq or refseq)
        c["nonsyn" if an["effect"] in ("nonsynonymous", "stop_gained") else "syn" if an["effect"] == "synonymous" else "inter"] += 1
        sites[g].add((pos, alt))
    interior = L - 2 * mutprofile.ITR_BP
    total = sum(c["events"] for c in cnt.values())
    rate = total / interior if interior > 0 else 0
    ng = max(1, len(by))
    out = {}
    for g, c in cnt.items():
        ln = glen.get(g)
        exp = rate * ln if ln else None
        p = _pois_sf(c["events"], exp) if exp else None
        out[g] = {"events": c["events"], "per_kb": round(c["events"] / ln * 1000, 2) if ln else None, "expected": round(exp, 2) if exp is not None else None, "p": p,
                  "p_bonferroni": min(1.0, p * ng) if p is not None else None, "distinct": len(sites[g]), "a3": c["a3"], "nonsyn": c["nonsyn"], "syn": c["syn"], "internal": c["internal"], "terminal": c["terminal"], "length": ln}
    return {"genes": out, "total_events": total, "rate_per_kb": round(rate * 1000, 3), "reference": refname,
            "recurrent_sites": sorted(({"label": f"{r}{p}{a}", "branches": n} for (r, p, a), n in Counter(k for _i, k in muts).items() if n >= 3), key=lambda x: -x["branches"])}


def _pois_cdf(k: int, mu: float) -> float:
    """P(X <= k) for X ~ Poisson(mu)."""
    if k < 0:
        return 0.0
    if mu <= 0:
        return 1.0
    term = math.exp(-mu)
    cdf = term
    for i in range(1, k + 1):
        term *= mu / i
        cdf += term
    return min(1.0, cdf)


def coldspot_model(events: list[tuple[str, int, str]], seq: str, regions: list[tuple[str, int, int, list[tuple[int, int]]]], itr: int, min_expected: float = 5.0, alpha: float = 0.05) -> list[dict]:
    """The core of the coldspot search, free of files so that it can be tested. `events` = [(ref, pos, alt)] independent mutation events, `seq` the sequence whose context is used (the tree's
    root), `regions` = [(name, first, last, exons)] to test (exons = the spans counted; for a window just [(first, last)]), `itr` the bp at each end that are left out.

    Mutations are not spread evenly: APOBEC3 events (TC>TT, GA>AA) land only on sites with that context. A region with few such sites would look cold under an even rate, so the expectation
    uses two rates - events of APOBEC3 type per APOBEC3-context site and all other events per other site - and counts the sites of each class inside the region.
    A region is a coldspot when it has fewer events than expected (lower-tail Poisson, Bonferroni over the regions tested), at least `min_expected` were expected (otherwise zero is not
    surprising) and p < alpha."""
    L = len(seq)
    lo, hi = itr + 1, L - itr                                   # 1-based interior
    a3site = lambda p: (seq[p - 1] == "C" and p >= 2 and seq[p - 2] == "T") or (seq[p - 1] == "G" and p < L and seq[p] == "A")
    ev = [(r, p, a) for r, p, a in events if lo <= p <= hi]
    ev_a3 = sum(1 for r, p, a in ev if insight._apo_branch(r, a, p, seq))
    ev_other = len(ev) - ev_a3
    n_a3 = sum(1 for p in range(lo, hi + 1) if a3site(p))
    n_other = (hi - lo + 1) - n_a3
    rate_a3 = ev_a3 / n_a3 if n_a3 else 0.0
    rate_other = ev_other / n_other if n_other else 0.0
    out = []
    for name, first, last, exons in regions:
        spans = [(max(a, lo), min(b, hi)) for a, b in exons if min(b, hi) >= max(a, lo)]
        ln = sum(b - a + 1 for a, b in spans)
        if ln < 1:
            continue
        sa3 = sum(1 for a, b in spans for p in range(a, b + 1) if a3site(p))
        exp = rate_a3 * sa3 + rate_other * (ln - sa3)
        k = sum(1 for r, p, a in ev if any(x <= p <= y for x, y in spans))
        p_low = _pois_cdf(k, exp)
        out.append({"name": name, "first": first, "last": last, "length": ln, "events": k, "expected": round(exp, 2), "per_kb": round(k / ln * 1000, 2), "a3_sites": sa3, "p": p_low})
    n = max(1, len(out))
    for x in out:
        x["p_bonferroni"] = min(1.0, x["p"] * n)
        low = x["expected"] >= min_expected and x["events"] < x["expected"]
        x["coldspot"] = bool(low and x["p_bonferroni"] < alpha)
        x["nominal"] = bool(low and x["p"] < alpha and not x["coldspot"])                  # below chance without the correction for the number of regions tested: a lead, not a finding
    return out


def coldspots(run: dict, window: int = 10000) -> dict | None:
    """Regions with significantly **fewer** independent mutation events than expected: genes (CDS) and, to catch non-coding stretches, windows of `window` bp. Same events as the hotspot
    search (each event once per tree branch, ITRs left out); the expectation takes the APOBEC3 sequence context into account (see `coldspot_model`). A coldspot is where mutations are
    rarer than chance - consistent with constraint (purifying selection) or a region that is hard to mutate - not proof of function."""
    tree, rootseq, nsrun = _tree(run)
    if not tree:
        return None
    refname, refseq, cds = _tree_reference(nsrun)
    b = insight._tree_branches(run)
    if not b or not refseq:
        return None
    _root, L, muts = b
    seq = rootseq or refseq
    events = [m for _internal, m in muts]
    genes = [(c["gene"], min(a for a, _ in c["exons"]), max(e for _, e in c["exons"]), c["exons"]) for c in cds]
    wins = [(f"{s}-{min(s + window - 1, L)}", s, min(s + window - 1, L), [(s, min(s + window - 1, L))]) for s in range(1, L + 1, window)]
    g_rows = coldspot_model(events, seq, genes, mutprofile.ITR_BP)
    w_rows = coldspot_model(events, seq, wins, mutprofile.ITR_BP)
    cat = categories()
    for x in g_rows:
        k = cat.get(x["name"], {})
        x["function"], x["category"] = k.get("function", ""), k.get("category", "")
    g_rows.sort(key=lambda x: (x["p"], -x["expected"]))
    w_rows.sort(key=lambda x: x["first"])
    return {"reference": refname, "total_events": len(events), "genes": g_rows, "windows": w_rows, "n_cold_genes": sum(1 for x in g_rows if x["coldspot"]),
            "n_cold_windows": sum(1 for x in w_rows if x["coldspot"]), "n_nominal_genes": sum(1 for x in g_rows if x["nominal"]), "window": window}


def genes(run: dict, hotspot: int = 100, fixed: float = 0.9, max_n: float | None = None) -> dict:
    ref = mutprofile.summary(run).get("reference", "")
    fa, gff = mutprofile.reference_by_name(ref) if ref else (None, None)
    refname, refseq = mutprofile.read_fasta(fa) if fa else ("", "")
    cds = mutprofile.load_cds(gff)
    cat = categories()
    skipfixed = _fixed_sites(run, fixed, max_n) if fixed < 1 else set()
    ev: dict[str, Counter] = defaultdict(Counter)
    distinct: dict[str, dict] = defaultdict(dict)
    L = len(refseq)
    for r in _rows(run, max_n):
        pos = int(r["pos"])
        if pos <= ITR or pos > L - ITR:
            continue
        g = r["gene"] or "(intergenic)"
        c = ev[g]
        isfixed = (pos, r["alt"]) in skipfixed
        e = r["effect"]
        a = insight._is_apo(r)
        c["all"] += 1
        c["a3" if a else "other"] += 1
        c["nonsyn" if e in ("nonsynonymous", "stop_gained") else "syn" if e == "synonymous" else "inter"] += 1
        if not isfixed:
            c["var_all"] += 1
            c["var_a3"] += a
            c["var_nonsyn"] += e in ("nonsynonymous", "stop_gained")
        distinct[g][(pos, r["alt"])] = (e, a)
    byname = {c["gene"]: c for c in cds}
    out = []
    for g, c in ev.items():
        d = distinct[g]
        nd = sum(1 for e, _a in d.values() if e in ("nonsynonymous", "stop_gained"))
        sd = sum(1 for e, _a in d.values() if e == "synonymous")
        dnds = None
        cd = byname.get(g)
        if cd and refseq and (nd + sd):
            seq = mutprofile._coding(cd, refseq)
            S = sum(_sites(seq[i:i + 3]) for i in range(0, len(seq) - 2, 3))
            N = len(seq) - S
            pn, ps = nd / N if N else 0, sd / S if S else 0
            jc = lambda p: -0.75 * math.log(1 - 4 * p / 3) if 0 < p < 0.75 else (0.0 if p == 0 else None)
            dn, ds_ = jc(pn), jc(ps)
            dnds = {"dN": round(dn, 5) if dn is not None else None, "dS": round(ds_, 5) if ds_ is not None else None, "ratio": round(dn / ds_, 2) if dn is not None and ds_ else None,
                    "N_sites": round(N, 1), "S_sites": round(S, 1)}
        k = cat.get(g, {})
        out.append({"gene": g, "function": k.get("function", ""), "category": k.get("category", "Unclassified" if g != "(intergenic)" else ""), "homolog": k.get("vaccinia_homolog", ""),
                    "events": c["all"], "a3": c["a3"], "a3_share": round(c["a3"] / c["all"], 3) if c["all"] else None, "nonsyn": c["nonsyn"], "syn": c["syn"], "inter": c["inter"],
                    "var_events": c["var_all"], "var_a3": c["var_a3"], "var_nonsyn": c["var_nonsyn"], "distinct": len(d), "distinct_nonsyn": nd, "distinct_syn": sd, "dnds": dnds,
                    "hotspot": c["all"] > hotspot, "var_hotspot": c["var_all"] > hotspot})
    ind = independent_by_gene(run)
    if ind:
        for x in out:
            v = ind["genes"].get(x["gene"])
            x["ind"] = v
            x["hotspot_ind"] = bool(v and v["p_bonferroni"] is not None and v["p_bonferroni"] < 0.05 and v["events"] >= 4 and x["gene"] != "(intergenic)")
        out.sort(key=lambda x: (-(x["ind"]["events"] if x.get("ind") else 0), -x["events"]))
    else:
        out.sort(key=lambda x: -x["events"])
    bycat = defaultdict(lambda: Counter())
    for x in out:
        if x["gene"] != "(intergenic)":
            bycat[x["category"]]["events"] += x["events"]
            bycat[x["category"]]["a3"] += x["a3"]
            bycat[x["category"]]["var_events"] += x["var_events"]
            bycat[x["category"]]["genes"] += 1
    return {"reference": refname, "hotspot": hotspot, "fixed_share": fixed, "genes": out[:120], "n_genes": len(out),
            "independent": ({k: ind[k] for k in ("total_events", "rate_per_kb", "reference", "recurrent_sites")} if ind else None),
            "categories": [{"category": k, **dict(v)} for k, v in sorted(bycat.items(), key=lambda kv: -kv[1]["events"])]}


def analysis(run: dict, min_size: int = 5, hotspot: int = 100, fixed: float = 0.9, max_n: float | None = None) -> dict:
    out = {"max_n": max_n}
    for name, fn in (("clusters", lambda: clusters(run, min_size)), ("accumulation", lambda: accumulation(run, fixed, max_n)), ("spectrum", lambda: spectrum(run, fixed, max_n)),
                     ("genes", lambda: genes(run, hotspot, fixed, max_n)), ("coldspots", lambda: coldspots(run))):
        try:
            out[name] = fn()
        except Exception as e:
            _log.warning("outbreak analysis %s failed", name, exc_info=True)
            out[name] = {"available": False, "reason": str(e)}
    return out
