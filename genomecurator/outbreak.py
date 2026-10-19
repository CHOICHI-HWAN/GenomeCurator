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
def clusters(run: dict, min_size: int = 5, min_muts: int = 1) -> dict:
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
    return {"available": True, "tree": nsrun["name"], "reference": refname, "min_size": min_size, "min_muts": min_muts, "tips": len(nodes_tips[id(tree)]), "clusters": [{k: v for k, v in c.items() if k != "tips"} for c in cl],
            "unassigned": len(nodes_tips[id(tree)]) - len(inside), "size_distribution": size_distribution(sizes), "columns": [{"label": lab, "gene": ann[lab]["gene"], "aa": ann[lab]["aa"], "apobec": ann[lab]["apobec"]} for lab in cols],
            "matrix": matrix, "months": sorted({m for c in cl for m in c["months"]}), "n_nested": len(qual)}


def _desc(r):
    out, stack = [], list(r["children"])
    while stack:
        x = stack.pop()
        out.append(x)
        stack.extend(x["children"])
    return out


# ------------------------------------------------------------------------------------------------ 3-4: per genome
def _fixed_sites(run: dict, share: float) -> set:
    by = defaultdict(set)
    for r in mutprofile._load_mutations(run):
        by[(int(r["pos"]), r["alt"])].add(r["sample"])
    n = len({r["sample"] for r in mutprofile._load_mutations(run)}) or 1
    return {k for k, s in by.items() if len(s) / n >= share}


def accumulation(run: dict, fixed: float = 0.9) -> dict:
    """APOBEC3-type mutations per genome against NC_063383.1 over collection date (own rule: G>A followed by A, C>T preceded by T), all and without the sites nearly every genome has."""
    ncname, ncseq, _cds = lineages.nc_reference()
    nc = lineages._nc_snps(run)
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
    return {"reference": ncname, "fixed_share": fixed, "n": len(rows), "by_year": by_year, "trend_all": trend(t, [r["apobec_all"] for r in rows]), "trend_variable": trend(t, [r["apobec_var"] for r in rows]),
            "points": [{k: r[k] for k in ("id", "date", "t", "apobec_all", "apobec_var", "snps_all", "snps_var", "lineage")} for r in rows]}


def spectrum(run: dict, fixed: float = 0.9) -> dict:
    """Per genome (chronological): APOBEC3-type / other x nonsynonymous / synonymous (stop gained counts as nonsynonymous; intergenic apart), against the run's reference."""
    gs = {g["id"]: g for g in lineages._genomes(run)}
    skip = _fixed_sites(run, fixed) if fixed < 1 else set()
    per: dict[str, Counter] = defaultdict(Counter)
    for r in mutprofile._load_mutations(run):
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


def genes(run: dict, hotspot: int = 100, fixed: float = 0.9) -> dict:
    ref = mutprofile.summary(run).get("reference", "")
    fa, gff = mutprofile.reference_by_name(ref) if ref else (None, None)
    refname, refseq = mutprofile.read_fasta(fa) if fa else ("", "")
    cds = mutprofile.load_cds(gff)
    cat = categories()
    skipfixed = _fixed_sites(run, fixed) if fixed < 1 else set()
    ev: dict[str, Counter] = defaultdict(Counter)
    distinct: dict[str, dict] = defaultdict(dict)
    L = len(refseq)
    for r in mutprofile._load_mutations(run):
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
    out.sort(key=lambda x: -x["events"])
    bycat = defaultdict(lambda: Counter())
    for x in out:
        if x["gene"] != "(intergenic)":
            bycat[x["category"]]["events"] += x["events"]
            bycat[x["category"]]["a3"] += x["a3"]
            bycat[x["category"]]["var_events"] += x["var_events"]
            bycat[x["category"]]["genes"] += 1
    return {"reference": refname, "hotspot": hotspot, "fixed_share": fixed, "genes": out[:120], "n_genes": len(out),
            "categories": [{"category": k, **dict(v)} for k, v in sorted(bycat.items(), key=lambda kv: -kv[1]["events"])]}


def analysis(run: dict, min_size: int = 5, hotspot: int = 100, fixed: float = 0.9) -> dict:
    out = {}
    for name, fn in (("clusters", lambda: clusters(run, min_size)), ("accumulation", lambda: accumulation(run, fixed)), ("spectrum", lambda: spectrum(run, fixed)), ("genes", lambda: genes(run, hotspot, fixed))):
        try:
            out[name] = fn()
        except Exception as e:
            _log.warning("outbreak analysis %s failed", name, exc_info=True)
            out[name] = {"available": False, "reason": str(e)}
    return out
