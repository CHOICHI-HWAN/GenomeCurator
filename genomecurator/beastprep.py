"""BEAST for the tMRCA: a package prepared from a clock run (`clock.py`), an optional run, and the posterior read back.

TreeTime (clock.py) is fast and shows whether there is a clock at all, but its intervals are conditional on one tree and a fitted rate. For a tMRCA that is meant to be cited, a Bayesian analysis
(BEAST: tip dates, a clock model, a tree prior, all uncertainty sampled together) is the standard. BEAST cannot take a thousand 197 kb genomes in reasonable time, so the package holds a
**subsample** chosen to keep the whole sampling period and every lineage (dated by day, not a clock outlier, not flagged by QC), the alignment of those genomes from the Squirrel alignment, tip dates, and a
BEAST 1.10-style XML (HKY + gamma sites; strict or uncorrelated log-normal relaxed clock; constant-size or exponential coalescent) that BEAST 10 reads. Nothing is hidden in the XML: the model is written out
and can be opened and changed in BEAUti / a text editor.

Reading the result: the tMRCA of the root and of the lineages as dates (mean / median / 95% HPD), the clock rate, and the effective sample size (ESS) of each - an ESS below 200 means the chain has not
run long enough and the interval must not be trusted. The package is a starting point: a real analysis needs several chains, checks in Tracer, and a model comparison."""
import json
from collections import Counter
import math
import random
import re
import shutil
import subprocess
import time
from datetime import date, datetime, timedelta
from pathlib import Path

from . import aligned, clock, log, nextstrain, runs

_log = log.get("beast")
DEFAULTS = {"n_tips": 150, "chain_length": 30000000, "samples": 10000, "clock": "relaxed", "tree_prior": "exponential", "seed": 1, "threads": 0}
CLOCKS = ("strict", "relaxed")
PRIORS = ("constant", "exponential")                     # the parametric tree priors the marginal-likelihood comparison (BETS / GSS) uses
TREE_PRIORS = PRIORS + ("skygrid",)
SUBST_MODELS = ("hky", "gtr")                            # HKY (one transition/transversion ratio) or GTR (six exchange rates)
GTR_RATES = ("ac", "ag", "at", "cg", "gt")               # BEAST fixes the sixth (C<->T) at 1; the XML elements are rateAC, rateAG, ...                       # all tree priors a BEAST analysis accepts


def env_prefix() -> Path | None:
    """Our own environment `gc-beast` (BEAST 10 + BEAGLE, `variantcall/beast.yml`), else one called `beast` that has BEAGLE too (BEAST 10 needs the BEAGLE library to run the likelihood)."""
    import os
    for name in (os.environ.get("GC_BEAST_ENV"), "gc-beast", "beast"):
        pre = nextstrain.find_env(name, "beast") if name else None
        if pre and any((pre / "lib").glob("libhmsbeagle*")):
            return pre
    return None


def status() -> dict:
    pre = env_prefix()
    ver = ""
    if pre:
        try:
            r = subprocess.run([str(pre / "bin" / "beast"), "-version"], capture_output=True, text=True, timeout=60)
            m = re.search(r"BEAST v([0-9.]+)", r.stdout + r.stderr)
            ver = m[1] if m else ""
        except Exception:
            ver = ""
    return {"ok": bool(pre), "version": ver, "how": "mamba env create -f variantcall/beast.yml   (gc-beast: BEAST 10 and the BEAGLE library, which BEAST 10 needs)"}


def parse_settings(raw: dict | None) -> dict:
    s = {**DEFAULTS, **(raw or {})}
    try:
        out = {"n_tips": int(s["n_tips"]), "chain_length": int(s["chain_length"]), "samples": int(s["samples"]), "seed": int(s["seed"]), "threads": int(s["threads"]), "clock": str(s["clock"]), "tree_prior": str(s["tree_prior"])}
    except (TypeError, ValueError):
        raise ValueError("BEAST 설정의 숫자가 올바르지 않습니다.")
    if not 10 <= out["n_tips"] <= 600:
        raise ValueError("표본 유전체 수는 10~600입니다.")
    if out["chain_length"] < 10000 or not 100 <= out["samples"] <= 100000:
        raise ValueError("체인 길이는 10,000 이상, 표본 수는 100~100,000이어야 합니다.")
    if out["clock"] not in CLOCKS or out["tree_prior"] not in TREE_PRIORS:
        raise ValueError(f"시계는 {', '.join(CLOCKS)}, 트리 사전분포는 {', '.join(TREE_PRIORS)} 중 하나입니다.")
    out["subst_model"] = str(s.get("subst_model") or "hky").lower() if str(s.get("subst_model") or "hky").lower() in SUBST_MODELS else "hky"
    out["site_rates"] = str(s.get("site_rates") or "gamma") if str(s.get("site_rates") or "gamma") in ("gamma", "inv", "gamma_inv") else "gamma"   # rate variation among sites: gamma(4), a proportion of invariant sites, or both
    out["skygrid_points"] = max(8, min(200, int(s.get("skygrid_points") or 36)))      # SkyGrid: number of grid cells (one per month by default)
    out["month_dates"] = bool(s.get("month_dates"))                                # use month-only dates at the middle of the month
    out["partition"] = str(s.get("partition") or "") if str(s.get("partition") or "") in ("apobec", "apobec_only", "other_only") else ""   # apobec: separate rates for the APOBEC3-context sites and the others; *_only: a plain model on one class of site
    out["mask_columns"] = sorted({int(x) for x in (s.get("mask_columns") or [])})   # alignment columns (0-based) left out of the alignment
    out["tag"] = re.sub(r"[^A-Za-z0-9_.\-]", "_", str(s.get("tag") or ""))[:40]
    out["tip_names"] = [str(x) for x in (s.get("tip_names") or [])]            # an explicit set of tips instead of the balanced subsample
    out["groups"] = {str(k): [str(x) for x in v] for k, v in (s.get("groups") or {}).items()}   # taxon sets whose MRCA age is logged (with tip_names)
    out["stem_groups"] = [str(x) for x in (s.get("stem_groups") or [])]        # groups whose stem age (the split from the sister lineage) is logged too
    rp = s.get("rate_prior")
    out["rate_prior"] = [float(rp[0]), float(rp[1])] if rp else None            # [mean, sd] of a log-normal prior on the clock rate, in log space
    return out


# ------------------------------------------------------------------------------------------------ the subsample
def choose_tips(tips: list[dict], n: int, seed: int = 1, month_dates: bool = False) -> list[dict]:
    """`n` of the clock run's tips, keeping the whole sampling period and the lineages: only genomes dated by day and not clock outliers; one genome from every month first (the earliest-dated of a
    random order), then more rounds from the busiest months and lineages until `n` are taken. Deterministic for a seed. `month_dates`: genomes whose date is only a month (precision `range`) are used as well, at the
    middle of the month (the usual convention; the clock run gives that date): up to half a month of error, about half a change at the rate of the mpox lineages studied here."""
    rng = random.Random(seed)
    pool = [t for t in tips if not t["outlier"] and (t["precision"] == "day" or (month_dates and t["precision"] == "range"))]
    if len(pool) <= n:
        return sorted(pool, key=lambda t: t["date"])
    by_month: dict[int, list] = {}
    for t in pool:
        by_month.setdefault(int(t["date"] * 12), []).append(t)
    for v in by_month.values():
        rng.shuffle(v)
    chosen: list[dict] = []
    lineages_seen: dict[str, int] = {}
    rounds = 0
    while len(chosen) < n and any(by_month.values()):
        for m in sorted(by_month, key=lambda k: -len(by_month[k])):
            if len(chosen) >= n:
                break
            if not by_month[m]:
                continue
            # within the month prefer a lineage not yet represented, so that rare lineages are kept
            lst = by_month[m]
            i = min(range(len(lst)), key=lambda j: lineages_seen.get(lst[j]["lineage"] or "-", 0))
            t = lst.pop(i)
            chosen.append(t)
            lineages_seen[t["lineage"] or "-"] = lineages_seen.get(t["lineage"] or "-", 0) + 1
        rounds += 1
    return sorted(chosen, key=lambda t: t["date"])


# ------------------------------------------------------------------------------------------------ the XML
def _safe(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", x)


def apobec_columns(seqs: dict[str, str]) -> tuple[list[int], list[int]]:
    """Alignment columns split by the context of the consensus base: (APOBEC3-type sites, all the others). A site is APOBEC3-type when the consensus there is a C right after a T (TC, where APOBEC3 makes
    C>T) or a G right before an A (GA, G>A) - the two contexts of the mutations the host enzymes make in human mpox transmission. Gaps and N do not vote."""
    names = list(seqs)
    L = len(seqs[names[0]])
    cons = []
    for i in range(L):
        c = Counter(seqs[n][i] for n in names if seqs[n][i] in "ACGT")
        cons.append(c.most_common(1)[0][0] if c else "N")
    a, b = [], []
    for i in range(L):
        isa = (cons[i] == "C" and i > 0 and cons[i - 1] == "T") or (cons[i] == "G" and i + 1 < L and cons[i + 1] == "A")
        (a if isa else b).append(i)
    return a, b


def build_xml(seqs: dict[str, str], dates: dict[str, float], groups: dict[str, list[str]], s: dict, rate0: float, fileprefix: str = "beast", mle: dict | None = None, proper: bool = False, fixed_freq: bool = False, gss: dict | None = None, init: dict | None = None) -> str:
    """The BEAST 1.10-style XML. `proper`: only proper priors (the clock rate and the population size get log-normal priors instead of the vague gamma / 1/x ones): a marginal likelihood, and so a Bayes factor, exists only for proper priors. `fixed_freq`: the base frequencies are the alignment's own and not sampled (needed for generalized stepping-stone sampling, which has no working prior for a simplex). `gss` = {"log", "burnin", "chain_length", "path_steps", "alpha"}: instead of path sampling / stepping-stone sampling, generalized stepping-stone sampling against a working prior built from the posterior samples of an earlier ordinary run of this model (`log` is read from the folder BEAST runs in, the first `burnin` states are discarded). `init` = {parameter id: starting value}: the chain starts there (for GSS, from the middle of the earlier posterior - a start outside the working prior has zero density and BEAST refuses it). `mle` = {"chain_length", "path_steps", "alpha"}: after the initial chain, path sampling and stepping-stone sampling of the marginal likelihood (`marginalLikelihoodEstimator`, written to `<prefix>.mle.log`). `seqs` {taxon: aligned sequence}, `dates` {taxon: decimal year}, `groups` {name: [taxa]} for the tMRCA of each group, `rate0` the starting clock rate (substitutions / site / year)."""
    names = list(seqs)
    sr = s.get("site_rates", "gamma")
    sm = s.get("subst_model", "hky")
    chain, every = s["chain_length"], max(1, s["chain_length"] // s["samples"])
    L = []
    w = L.append
    w('<?xml version="1.0" standalone="yes"?>')
    w('<!-- Written by Virus Genome Metadata Curator (genomecurator/beastprep.py). Model: HKY + gamma(4) sites, ' + ("strict clock" if s["clock"] == "strict" else "uncorrelated log-normal relaxed clock") + ", " + ("constant-size" if s["tree_prior"] == "constant" else "exponential-growth") + " coalescent tree prior, tip dates in years. Open in BEAUti or a text editor to change it. -->")
    w('<beast version="1.10.4">')
    w('\t<taxa id="taxa">')
    for n in names:
        w(f'\t\t<taxon id="{_safe(n)}"><date value="{dates[n]:.5f}" direction="forwards" units="years"/></taxon>')
    w('\t</taxa>')
    for g, members in groups.items():
        w(f'\t<taxa id="group.{_safe(g)}">')
        for n in members:
            w(f'\t\t<taxon idref="{_safe(n)}"/>')
        w('\t</taxa>')
    part = s.get("partition") == "apobec"
    if (part or s.get("partition") in ("apobec_only", "other_only")) and gss:
        raise ValueError("분할 모델(APOBEC3 자리 대 그 외)은 한계우도 추정(GSS)과 함께 쓸 수 없습니다.")
    if s.get("partition") in ("apobec_only", "other_only"):                     # one class of site, a plain (unpartitioned) model on it
        ca, cb = apobec_columns(seqs)
        cols = ca if s["partition"] == "apobec_only" else cb
        sets, nsites = [("", {n: "".join(seqs[n][i] for i in cols) for n in names})], [len(cols)]
    elif part:
        ca, cb = apobec_columns(seqs)
        sets = [(".a", {n: "".join(seqs[n][i] for i in ca) for n in names}), (".b", {n: "".join(seqs[n][i] for i in cb) for n in names})]
        nsites = [len(ca), len(cb)]
    else:
        sets, nsites = [("", seqs)], [len(seqs[names[0]])]
    for suf, sq in sets:
        w(f'\t<alignment id="alignment{suf}" dataType="nucleotide">')
        for n in names:
            w(f'\t\t<sequence><taxon idref="{_safe(n)}"/>{sq[n]}</sequence>')
        w('\t</alignment>')
        w(f'\t<patterns id="patterns{suf}" from="1" strip="false"><alignment idref="alignment{suf}"/></patterns>')
    sky = s["tree_prior"] == "skygrid"
    if sky:                                                                   # the starting tree only needs a demographic model; the SkyGrid itself is the tree prior below
        w('\t<constantSize id="demographic" units="years"><populationSize><parameter id="demographic.popSize" value="1.0" lower="0.0"/></populationSize></constantSize>')
    elif s["tree_prior"] == "constant":
        w('\t<constantSize id="demographic" units="years"><populationSize><parameter id="demographic.popSize" value="1.0" lower="0.0"/></populationSize></constantSize>')
    else:
        w('\t<exponentialGrowth id="demographic" units="years"><populationSize><parameter id="demographic.popSize" value="1.0" lower="0.0"/></populationSize><growthRate><parameter id="demographic.growthRate" value="1.0"/></growthRate></exponentialGrowth>')
    w('\t<coalescentSimulator id="startingTree"><taxa idref="taxa"/><' + ("exponentialGrowth" if s["tree_prior"] == "exponential" else "constantSize") + ' idref="demographic"/></coalescentSimulator>')
    w('\t<treeModel id="treeModel"><coalescentTree idref="startingTree"/><rootHeight><parameter id="treeModel.rootHeight"/></rootHeight>'
      '<nodeHeights internalNodes="true"><parameter id="treeModel.internalNodeHeights"/></nodeHeights><nodeHeights internalNodes="true" rootNode="true"><parameter id="treeModel.allInternalNodeHeights"/></nodeHeights></treeModel>')
    w('\t<treeLengthStatistic id="treeLength"><treeModel idref="treeModel"/></treeLengthStatistic>')
    w('\t<tmrcaStatistic id="age(root)" absolute="true"><treeModel idref="treeModel"/><mrca><taxa idref="taxa"/></mrca></tmrcaStatistic>')
    for g in groups:
        w(f'\t<tmrcaStatistic id="age({_safe(g)})" absolute="true"><treeModel idref="treeModel"/><mrca><taxa idref="group.{_safe(g)}"/></mrca></tmrcaStatistic>')
        if g in s.get("stem_groups", []):
            w(f'\t<tmrcaStatistic id="stem({_safe(g)})" absolute="true" includeStem="true"><treeModel idref="treeModel"/><mrca><taxa idref="group.{_safe(g)}"/></mrca></tmrcaStatistic>')
            w(f'\t<monophylyStatistic id="monophyly({_safe(g)})"><mrca><taxa idref="group.{_safe(g)}"/></mrca><treeModel idref="treeModel"/></monophylyStatistic>')       # 1 in the samples where the group is a clade
    if sky:
        span = max(dates.values()) - min(dates.values()); cut = round(span + 1.0, 2); npts = int(s.get("skygrid_points", 36))   # the grid reaches one year beyond the sampled span
        w(f'\t<gmrfSkyGridLikelihood id="skygrid"><populationSizes><parameter id="skygrid.logPopSize" dimension="{npts + 1}" value="1.0"/></populationSizes>'
          f'<precisionParameter><parameter id="skygrid.precision" value="0.1" lower="0.0"/></precisionParameter><numGridPoints><parameter id="skygrid.numGridPoints" value="{npts}"/></numGridPoints>'
          f'<cutOff><parameter id="skygrid.cutOff" value="{cut}"/></cutOff><populationTree><treeModel idref="treeModel"/></populationTree></gmrfSkyGridLikelihood>')
    else:
        w('\t<coalescentLikelihood id="coalescent"><model><' + ("constantSize" if s["tree_prior"] == "constant" else "exponentialGrowth") + ' idref="demographic"/></model><populationTree><treeModel idref="treeModel"/></populationTree></coalescentLikelihood>')
    if s["clock"] == "strict":
        w(f'\t<strictClockBranchRates id="branchRates"><rate><parameter id="clock.rate" value="{rate0:.4e}" lower="0.0"/></rate></strictClockBranchRates>')
        rate_param, rate_ref = "clock.rate", '<strictClockBranchRates idref="branchRates"/>'
    else:
        w(f'\t<discretizedBranchRates id="branchRates"><treeModel idref="treeModel"/><distribution><logNormalDistributionModel meanInRealSpace="true"><mean><parameter id="ucld.mean" value="{rate0:.4e}" lower="0.0"/></mean>'
          '<stdev><parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0"/></stdev></logNormalDistributionModel></distribution><rateCategories><parameter id="branchRates.categories"/></rateCategories></discretizedBranchRates>')
        rate_param, rate_ref = "ucld.mean", '<discretizedBranchRates idref="branchRates"/>'
    w(f'\t<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true"><treeModel idref="treeModel"/>{rate_ref}</rateStatistic>')
    for suf, sq in sets:
        if fixed_freq:
            cnt = {b: sum(q.count(b) for q in sq.values()) for b in "ACGT"}
            tot = max(1, sum(cnt.values()))
            fr = [round(cnt[b] / tot, 10) for b in "ACG"]
            freq_value = " ".join(f"{x:.10f}" for x in fr) + f" {1 - sum(fr):.10f}"                 # the last one is the remainder: BEAST refuses frequencies that do not sum to 1 (6 decimals gave 1.000001)
        else:
            freq_value = "0.25 0.25 0.25 0.25"
        freqs = f'<frequencies><frequencyModel dataType="nucleotide"><frequencies><parameter id="frequencies{suf}" value="' + freq_value + '"/></frequencies></frequencyModel></frequencies>'
        if sm == "gtr":
            w(f'\t<gtrModel id="gtr{suf}">' + freqs + "".join(f'<rate{r.upper()}><parameter id="{r}{suf}" value="1.0" lower="0.0"/></rate{r.upper()}>' for r in GTR_RATES) + '</gtrModel>')
            w(f'\t<compoundParameter id="allRates{suf}">' + "".join(f'<parameter idref="{r}{suf}"/>' for r in GTR_RATES) + '</compoundParameter>')
        else:
            w(f'\t<HKYModel id="hky{suf}">' + freqs + f'<kappa><parameter id="kappa{suf}" value="2.0" lower="0.0"/></kappa></HKYModel>')
        w(f'\t<siteModel id="siteModel{suf}"><substitutionModel>' + (f'<gtrModel idref="gtr{suf}"/>' if sm == "gtr" else f'<HKYModel idref="hky{suf}"/>') + '</substitutionModel>' + (f'<relativeRate><parameter id="mu{suf}" value="1.0" lower="0.0"/></relativeRate>' if part else "")
          + (f'<gammaShape gammaCategories="4"><parameter id="alpha{suf}" value="0.5" lower="0.0"/></gammaShape>' if sr != "inv" else "")
          + (f'<proportionInvariant><parameter id="pInv{suf}" value="0.5" lower="0.0" upper="1.0"/></proportionInvariant>' if sr != "gamma" else "") + '</siteModel>')
    if part:
        w('\t<compoundParameter id="allMus"><parameter idref="mu.a"/><parameter idref="mu.b"/></compoundParameter>')
    w('\t<treeDataLikelihood id="treeLikelihood" useAmbiguities="false">' + "".join(f'<partition><patterns idref="patterns{suf}"/><siteModel idref="siteModel{suf}"/></partition>' for suf, _q in sets) + f'<treeModel idref="treeModel"/>{rate_ref}</treeDataLikelihood>')
    w('\t<operators id="operators" optimizationSchedule="default">')
    for suf, _q in sets:
        if sm == "gtr":
            for r in GTR_RATES:
                w(f'\t\t<scaleOperator scaleFactor="0.75" weight="0.1"><parameter idref="{r}{suf}"/></scaleOperator>')
        else:
            w(f'\t\t<scaleOperator scaleFactor="0.75" weight="0.1"><parameter idref="kappa{suf}"/></scaleOperator>')
        if sr != "inv":
            w(f'\t\t<scaleOperator scaleFactor="0.75" weight="0.1"><parameter idref="alpha{suf}"/></scaleOperator>')
        if sr != "gamma":
            w(f'\t\t<randomWalkOperator windowSize="0.1" weight="0.1"><parameter idref="pInv{suf}"/></randomWalkOperator>')
        if not fixed_freq:
            w(f'\t\t<deltaExchange delta="0.01" weight="0.1"><parameter idref="frequencies{suf}"/></deltaExchange>')
    if part:                                                                                     # the relative rates of the two partitions: their mean over the sites stays 1, so the clock rate is the genome-wide rate
        w(f'\t\t<deltaExchange delta="0.75" weight="3" parameterWeights="{nsites[0]} {nsites[1]}"><parameter idref="allMus"/></deltaExchange>')
    w(f'\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="{rate_param}"/></scaleOperator>')
    if s["clock"] == "relaxed":
        w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="ucld.stdev"/></scaleOperator>')
        w('\t\t<swapOperator size="1" weight="10" autoOptimize="false"><parameter idref="branchRates.categories"/></swapOperator>')
        w('\t\t<uniformIntegerOperator weight="10"><parameter idref="branchRates.categories"/></uniformIntegerOperator>')
    w('\t\t<subtreeSlide size="0.1" gaussian="true" weight="30"><treeModel idref="treeModel"/></subtreeSlide>')
    w('\t\t<narrowExchange weight="30"><treeModel idref="treeModel"/></narrowExchange>')
    w('\t\t<wideExchange weight="3"><treeModel idref="treeModel"/></wideExchange>')
    w('\t\t<wilsonBalding weight="3"><treeModel idref="treeModel"/></wilsonBalding>')
    w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="treeModel.rootHeight"/></scaleOperator>')
    w('\t\t<uniformOperator weight="30"><parameter idref="treeModel.internalNodeHeights"/></uniformOperator>')
    if sky:
        w('\t\t<gmrfGridBlockUpdateOperator scaleFactor="1.0" weight="30"><gmrfSkyrideLikelihood idref="skygrid"/></gmrfGridBlockUpdateOperator>')
        w('\t\t<scaleOperator scaleFactor="0.75" weight="1"><parameter idref="skygrid.precision"/></scaleOperator>')
    else:
        w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="demographic.popSize"/></scaleOperator>')
    if s["tree_prior"] == "exponential":
        w('\t\t<randomWalkOperator windowSize="1.0" weight="3"><parameter idref="demographic.growthRate"/></randomWalkOperator>')
    w('\t</operators>')
    w(f'\t<mcmc id="mcmc" chainLength="{chain}" autoOptimize="true" operatorAnalysis="{fileprefix}.ops">')
    w('\t\t<posterior id="posterior"><prior id="prior">')
    for suf, _q in sets:
        if sm == "gtr":
            w(f'\t\t\t<gammaPrior shape="0.05" scale="20.0" offset="0.0">' + "".join(f'<parameter idref="{r}{suf}"/>' for r in GTR_RATES) + '</gammaPrior>')
        else:
            w(f'\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false"><parameter idref="kappa{suf}"/></logNormalPrior>')
        if sr != "inv":
            w(f'\t\t\t<exponentialPrior mean="0.5" offset="0.0"><parameter idref="alpha{suf}"/></exponentialPrior>')
        if sr != "gamma":
            w(f'\t\t\t<uniformPrior lower="0.0" upper="1.0"><parameter idref="pInv{suf}"/></uniformPrior>')
        if not fixed_freq:
            w(f'\t\t\t<dirichletPrior alpha="1.0" sumsTo="1.0"><parameter idref="frequencies{suf}"/></dirichletPrior>')
    if part:
        w('\t\t\t<exponentialPrior mean="1.0" offset="0.0"><parameter idref="mu.a"/></exponentialPrior>')
        w('\t\t\t<exponentialPrior mean="1.0" offset="0.0"><parameter idref="mu.b"/></exponentialPrior>')
    if s.get("rate_prior"):
        w(f'\t\t\t<logNormalPrior mean="{s["rate_prior"][0]}" stdev="{s["rate_prior"][1]}" offset="0.0" meanInRealSpace="false"><parameter idref="{rate_param}"/></logNormalPrior>')      # an informative prior on the rate (see the settings)
    elif proper:
        w(f'\t\t\t<logNormalPrior mean="-11.5" stdev="2.0" offset="0.0" meanInRealSpace="false"><parameter idref="{rate_param}"/></logNormalPrior>')        # median 1e-5, four orders of magnitude wide
    else:
        w(f'\t\t\t<gammaPrior shape="0.001" scale="1000.0" offset="0.0"><parameter idref="{rate_param}"/></gammaPrior>')
    if s["clock"] == "relaxed":
        w('\t\t\t<exponentialPrior mean="0.3333333333333333" offset="0.0"><parameter idref="ucld.stdev"/></exponentialPrior>')
    if sky:
        w('\t\t\t<gammaPrior shape="0.001" scale="1000.0" offset="0.0"><parameter idref="skygrid.precision"/></gammaPrior>')       # the usual vague prior on the smoothing precision
    elif proper:
        w('\t\t\t<logNormalPrior mean="0.0" stdev="2.0" offset="0.0" meanInRealSpace="false"><parameter idref="demographic.popSize"/></logNormalPrior>')
    else:
        w('\t\t\t<oneOnXPrior><parameter idref="demographic.popSize"/></oneOnXPrior>')
    if s["tree_prior"] == "exponential":
        w('\t\t\t<laplacePrior mean="0.0" scale="100.0"><parameter idref="demographic.growthRate"/></laplacePrior>')
    w('\t\t\t<gmrfSkyGridLikelihood idref="skygrid"/>' if sky else '\t\t\t<coalescentLikelihood idref="coalescent"/>')
    w('\t\t</prior><likelihood id="likelihood"><treeDataLikelihood idref="treeLikelihood"/></likelihood></posterior>')
    w('\t\t<operators idref="operators"/>')
    w(f'\t\t<log id="screenLog" logEvery="{every * 10}"><column label="Posterior" dp="4" width="12"><posterior idref="posterior"/></column><column label="Prior" dp="4" width="12"><prior idref="prior"/></column>'
      '<column label="Likelihood" dp="4" width="12"><likelihood idref="likelihood"/></column><column label="rootAge" sf="6" width="12"><tmrcaStatistic idref="age(root)"/></column></log>')
    w(f'\t\t<log id="fileLog" logEvery="{every}" fileName="{fileprefix}.log" overwrite="false"><posterior idref="posterior"/><prior idref="prior"/><likelihood idref="likelihood"/>'
      '<parameter idref="treeModel.rootHeight"/><tmrcaStatistic idref="age(root)"/>' + "".join(f'<tmrcaStatistic idref="age({_safe(g)})"/>' + (f'<tmrcaStatistic idref="stem({_safe(g)})"/><monophylyStatistic idref="monophyly({_safe(g)})"/>' if g in s.get("stem_groups", []) else "") for g in groups) + '<treeLengthStatistic idref="treeLength"/>' +
      ('<parameter idref="skygrid.precision"/><parameter idref="skygrid.logPopSize"/>' if sky else '<parameter idref="demographic.popSize"/>') + ('<parameter idref="demographic.growthRate"/>' if s["tree_prior"] == "exponential" else "") + "".join(("".join(f'<parameter idref="{r}{suf}"/>' for r in GTR_RATES) if sm == "gtr" else f'<parameter idref="kappa{suf}"/>') + (f'<parameter idref="alpha{suf}"/>' if sr != "inv" else "") + (f'<parameter idref="pInv{suf}"/>' if sr != "gamma" else "") for suf, _q in sets) + ('<parameter idref="mu.a"/><parameter idref="mu.b"/>' if part else "") +
      f'<parameter idref="{rate_param}"/>' + ('<parameter idref="ucld.stdev"/>' if s["clock"] == "relaxed" else "") + '<rateStatistic idref="meanRate"/><treeDataLikelihood idref="treeLikelihood"/>' + ('<gmrfSkyGridLikelihood idref="skygrid"/>' if sky else '<coalescentLikelihood idref="coalescent"/>') + '</log>')
    w(f'\t\t<logTree id="treeFileLog" logEvery="{every}" nexusFormat="true" fileName="{fileprefix}.trees" sortTranslationTable="true"><treeModel idref="treeModel"/>'
      f'<trait name="rate" tag="rate">{rate_ref}</trait><posterior idref="posterior"/></logTree>')
    w('\t</mcmc>')
    w('\t<report><property name="timer"><mcmc idref="mcmc"/></property></report>')
    if gss:
        f_, bi = gss["log"], gss["burnin"]
        ref = lambda tag, col: f'<{tag} fileName="{f_}" parameterColumn="{col}" burnin="{bi}">'
        w('\t<!-- generalized stepping-stone sampling: the working prior is built from the posterior samples of an earlier ordinary run of this model (' + f_ + ') -->')
        if s["tree_prior"] == "constant":
            w(f'\t<constantSize id="demographicReference" units="years"><populationSize><parameter id="demographicReference.popSize" value="0.1" parameterColumn="demographic.popSize" fileName="{f_}" burnin="{bi}" dimension="1"/></populationSize></constantSize>')
            w('\t<coalescentLikelihood id="coalescentReference"><model><constantSize idref="demographicReference"/></model><populationTree><treeModel idref="treeModel"/></populationTree></coalescentLikelihood>')
        else:
            w(f'\t<exponentialGrowth id="demographicReference" units="years"><populationSize><parameter id="demographicReference.popSize" value="0.1" parameterColumn="demographic.popSize" fileName="{f_}" burnin="{bi}" dimension="1"/></populationSize>'
              f'<growthRate><parameter id="demographicReference.growthRate" value="0.1" parameterColumn="demographic.growthRate" fileName="{f_}" burnin="{bi}" dimension="1"/></growthRate></exponentialGrowth>')
            w('\t<coalescentLikelihood id="coalescentReference"><model><exponentialGrowth idref="demographicReference"/></model><populationTree><treeModel idref="treeModel"/></populationTree></coalescentLikelihood>')
        w(f'\t<marginalLikelihoodEstimator id="MLE" chainLength="{gss["chain_length"]}" pathSteps="{gss["path_steps"]}" pathScheme="betaquantile" alpha="{gss.get("alpha", 0.3)}">')
        w('\t\t<samplers><mcmc idref="mcmc"/></samplers>')
        w('\t\t<pathLikelihood id="pathLikelihood"><source><posterior idref="posterior"/></source><destination><workingPrior id="workingPrior">')
        for par in ((GTR_RATES if sm == "gtr" else ("kappa",)) + ("alpha", rate_param)) + (("ucld.stdev",) if s["clock"] == "relaxed" else ()) + ("demographic.popSize",):
            w(f'\t\t\t{ref("logTransformedNormalReferencePrior", par)}<parameter idref="{par}"/></logTransformedNormalReferencePrior>')
        if s["tree_prior"] == "exponential":
            w(f'\t\t\t{ref("normalReferencePrior", "demographic.growthRate")}<parameter idref="demographic.growthRate"/></normalReferencePrior>')
        w('\t\t\t<coalescentLikelihood idref="coalescentReference"/>')
        w('\t\t</workingPrior></destination></pathLikelihood>')
        w(f'\t\t<log id="MLELog" logEvery="{max(1, gss["chain_length"] // 100)}" fileName="{fileprefix}.mle.log"><pathLikelihood idref="pathLikelihood"/></log>')
        w('\t</marginalLikelihoodEstimator>')
        w(f'\t<generalizedSteppingStoneSamplingAnalysis fileName="{fileprefix}.mle.log"><sourceColumn name="pathLikelihood.source"/><destinationColumn name="pathLikelihood.destination"/><thetaColumn name="pathLikelihood.theta"/></generalizedSteppingStoneSamplingAnalysis>')
    if mle:
        w(f'\t<marginalLikelihoodEstimator id="MLE" chainLength="{mle["chain_length"]}" pathSteps="{mle["path_steps"]}" pathScheme="betaquantile" alpha="{mle.get("alpha", 0.3)}">')
        w('\t\t<samplers><mcmc idref="mcmc"/></samplers>')
        w('\t\t<pathLikelihood id="pathLikelihood"><source><posterior idref="posterior"/></source><destination><prior idref="prior"/></destination></pathLikelihood>')
        w(f'\t\t<log id="MLELog" logEvery="{max(1, mle["chain_length"] // 100)}" fileName="{fileprefix}.mle.log"><pathLikelihood idref="pathLikelihood"/></log>')
        w('\t</marginalLikelihoodEstimator>')
        for el in ("pathSamplingAnalysis", "steppingStoneSamplingAnalysis"):
            w(f'\t<{el} fileName="{fileprefix}.mle.log"><likelihoodColumn name="pathLikelihood.delta"/><thetaColumn name="pathLikelihood.theta"/></{el}>')
    w('</beast>')
    xml = "\n".join(L) + "\n"
    for name, val in (init or {}).items():
        xml = re.sub(rf'(<parameter id="{re.escape(name)}" value=")[^"]*(")', lambda m: f"{m[1]}{val:.8g}{m[2]}", xml, count=1)
    return xml


# ------------------------------------------------------------------------------------------------ package, run, results
def folder(clock_run: dict, tag: str = "") -> Path:
    """`<clock run>/beast/`, or `beast/<tag>/` for a further analysis of the same clock run (e.g. one lineage with its own tips)."""
    return Path(clock_run["run_dir"]) / "beast" / (_safe(tag) if tag else "")


def prepare(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """Write the BEAST package of a clock run into `<clock run>/beast/`: subsample.tsv, alignment.fasta, dates.tsv, beast.xml, README.txt. -> {n, taxa, groups, files}"""
    p = progress or (lambda *a: None)
    s = parse_settings(settings)
    S = clock.summary(clock_run)
    if not S.get("tips"):
        raise ValueError("이 분자시계 실행에는 결과가 없습니다.")
    src = next((r for r in nextstrain.list_runs(clock_run["dataset_id"]) if r["id"] == S["source_run"]), None)
    if not src or src["tool"] != "squirrel":
        raise ValueError("BEAST 패키지는 Squirrel 정렬이 필요합니다. 분자시계를 Squirrel 트리로 실행하세요.")
    p(5, "표본을 고르는 중…")
    if s["tip_names"]:
        names = set(s["tip_names"])
        chosen = sorted((t for t in S["tips"] if t["name"] in names and (t["precision"] == "day" or (s["month_dates"] and t["precision"] == "range"))), key=lambda t: t["date"])
        missing = names - {t["name"] for t in chosen}
        if missing:
            raise ValueError(f"고른 유전체 중 이 분자시계 실행에 없거나 날짜가 일 단위가 아닌 것이 있습니다: {', '.join(sorted(missing)[:8])}")
    else:
        chosen = choose_tips(S["tips"], s["n_tips"], s["seed"], s["month_dates"])
    if len(chosen) < 10:
        raise ValueError(f"BEAST에 쓸 수 있는 유전체(날짜가 일 단위이고 이상치가 아님)가 {len(chosen)}개뿐입니다.")
    want = {t["name"] for t in chosen}
    p(15, "정렬에서 서열을 읽는 중…")
    seqs = {}
    for name, seq in aligned.read_alignment(Path(src["run_dir"]) / "squirrel" / "alignment.aln.fasta.gz"):
        if name in want:
            seqs[name] = seq.upper()
    chosen = [t for t in chosen if t["name"] in seqs]
    if s["mask_columns"]:                                                      # columns left out (e.g. sites that mutated again and again): the alignment is shortened
        drop = set(s["mask_columns"])
        seqs = {n: "".join(c for i, c in enumerate(q) if i not in drop) for n, q in seqs.items()}
    dates = {t["name"]: t["date"] for t in chosen}
    groups: dict[str, list[str]] = {}
    for t in chosen:
        for kind, val in (("clade", t["clade"]), ("lineage", t["lineage"])):
            if val:
                groups.setdefault(f"{kind}.{val}", []).append(t["name"])
    groups = {g: m for g, m in groups.items() if 5 <= len(m) < len(chosen)}
    if s["tip_names"]:
        groups = {g: [n for n in m if n in seqs] for g, m in s["groups"].items()}
        groups = {g: m for g, m in groups.items() if len(m) >= 2}
    d = folder(clock_run, s["tag"])
    d.mkdir(parents=True, exist_ok=True)
    rate0 = (S.get("time") or {}).get("rate") or (S.get("clock") or {}).get("rate") or 1e-5
    p(40, "XML을 쓰는 중…")
    (d / "beast.xml").write_text(build_xml(seqs, dates, groups, s, rate0), encoding="utf-8")
    with open(d / "alignment.fasta", "w", encoding="utf-8") as f:
        for n, q in seqs.items():
            f.write(f">{n}\n{q}\n")
    (d / "dates.tsv").write_text("name\tdate\tlineage\tclade\tcountry\tprecision\n" + "".join(f"{t['name']}\t{t['date']:.5f}\t{t['lineage']}\t{t['clade']}\t{t['country']}\t{'month (middle of the month)' if t['precision'] == 'range' else 'day'}\n" for t in chosen), encoding="utf-8")
    (d / "subsample.tsv").write_text("name\tdate\tlineage\n" + "".join(f"{t['name']}\t{t['date']:.5f}\t{t['lineage']}\n" for t in chosen), encoding="utf-8")
    (d / "settings.json").write_text(json.dumps(s), encoding="utf-8")
    (d / "README.txt").write_text(
        "BEAST package written by Virus Genome Metadata Curator.\n\n"
        f"{len(chosen)} genomes subsampled from the {len(S['tips'])} of the clock run (day-precision dates, no clock outliers), alignment length {len(next(iter(seqs.values()))):,}.\n"
        "Run:   conda run -n beast beast -overwrite -seed 1 beast.xml     (or open beast.xml in BEAUti / BEAST)\n"
        "Read:  beast.log (Tracer), beast.trees (TreeAnnotator). Check that every ESS is above 200 before using an interval.\n"
        "Model: see the comment at the top of beast.xml. This is a starting point: run several chains with different seeds and compare models before you cite a tMRCA.\n", encoding="utf-8")
    p(100, "완료")
    return {"n": len(chosen), "taxa": len(seqs), "groups": list(groups), "settings": s, "folder": str(d), "files": sorted(f.name for f in d.iterdir())}


def ess(x: list[float]) -> float:
    """Effective sample size of a chain (initial-positive-sequence estimator, as Tracer does it)."""
    n = len(x)
    if n < 10:
        return float(n)
    m = sum(x) / n
    v = sum((a - m) ** 2 for a in x) / n
    if v == 0:
        return float(n)
    tot, prev = 0.0, None
    for lag in range(1, n // 2):
        c = sum((x[i] - m) * (x[i + lag] - m) for i in range(n - lag)) / n / v
        if prev is not None and (c + prev) < 0:
            break
        tot += c
        prev = c
        if lag > 2000:
            break
    return max(1.0, n / (1 + 2 * tot)) if tot > -0.5 else float(n)


def hpd(x: list[float], mass: float = 0.95) -> tuple[float, float]:
    s = sorted(x)
    k = max(1, int(math.floor(mass * len(s))))
    best = min(range(len(s) - k + 1), key=lambda i: s[i + k - 1] - s[i]) if len(s) > k else 0
    return s[best], s[min(len(s) - 1, best + k - 1)]


def read_log(path: Path, burnin: float = 0.1) -> dict:
    """The posterior of a BEAST log: for the root and each group's age, the clock rate and the posterior: mean, median, 95% HPD and ESS after discarding `burnin` of the chain."""
    rows, header = [], None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if header is None:
            header = p
            continue
        if len(p) == len(header):
            rows.append(p)
    if header is None or len(rows) < 20:
        return {"states": len(rows), "params": {}}
    cut = int(len(rows) * burnin)
    out = {}
    for j, h in enumerate(header):
        if h == "state":
            continue
        try:
            xs = [float(r[j]) for r in rows[cut:]]
        except ValueError:
            continue
        if not xs:
            continue
        lo, hi = hpd(xs)
        xs_s = sorted(xs)
        out[h] = {"mean": sum(xs) / len(xs), "median": xs_s[len(xs) // 2], "lo": lo, "hi": hi, "ess": ess(xs)}
    if "mu.a" in header and "mu.b" in header:                       # the partitioned model: the rates of the two kinds of site per sample (clock rate x relative rate) and their ratio
        rk = "ucld.mean" if "ucld.mean" in header else "clock.rate"
        ia, ib, ir = header.index("mu.a"), header.index("mu.b"), header.index(rk)
        sub = rows[cut:]
        for name, f in (("rate.apobec", lambda r: float(r[ir]) * float(r[ia])), ("rate.other", lambda r: float(r[ir]) * float(r[ib])), ("rate.ratio", lambda r: float(r[ia]) / max(1e-12, float(r[ib])))):
            xs = [f(r) for r in sub]
            lo, hi = hpd(xs)
            out[name] = {"mean": sum(xs) / len(xs), "median": sorted(xs)[len(xs) // 2], "lo": lo, "hi": hi, "ess": ess(xs)}
    return {"states": len(rows), "last_state": int(float(rows[-1][0])), "burnin": burnin, "params": out}


def run(clock_run: dict, settings: dict | None = None, progress=None, threads: int | None = None) -> dict:
    """Prepare the package and run BEAST on it (a long job). The result is written as `beast/summary.json` when the chain ends; the log can be read while it runs (`read_log`)."""
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"BEAST 환경이 없습니다: {st['how']}")
    p = progress or (lambda *a: None)
    info = prepare(clock_run, settings, lambda pct, m="": p(pct * 0.02, m))
    s = info["settings"]
    d = folder(clock_run, s["tag"])
    for f in ("beast.log", "beast.trees", "beast.ops", "beast.screen.log", "summary.json"):
        (d / f).unlink(missing_ok=True)
    pre = env_prefix()
    n_threads = threads or s["threads"]
    cmd = [str(pre / "bin" / "beast"), "-overwrite", "-seed", str(s["seed"])] + (["-threads", str(n_threads)] if n_threads else []) + ["beast.xml"]
    t0 = time.time()
    with open(d / "beast.screen.log", "w", encoding="utf-8") as lf:
        proc = subprocess.Popen(cmd, cwd=str(d), env=nextstrain.tool_env(pre), stdout=lf, stderr=subprocess.STDOUT, text=True)
        try:
            while proc.poll() is None:
                time.sleep(5)
                lg = d / "beast.log"
                if lg.exists():
                    try:
                        last = lg.read_text(encoding="utf-8").strip().splitlines()[-1].split("\t")[0]
                        done = float(last) / s["chain_length"]
                        p(2 + 96 * min(done, 1.0), f"{int(float(last)):,}/{s['chain_length']:,} 상태 ({time.time() - t0:.0f}초)")
                    except (ValueError, IndexError):
                        pass
        finally:
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    proc.kill()
    if proc.returncode != 0:
        tail = "".join(open(d / "beast.screen.log", encoding="utf-8").readlines()[-8:])
        raise RuntimeError(f"BEAST가 실패했습니다 (종료 코드 {proc.returncode}): {tail.strip()[-400:]}")
    res = summarize(clock_run, s["tag"])
    res["seconds"] = round(time.time() - t0)
    (d / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    p(100, "BEAST 완료")
    return res


def summarize(clock_run: dict, tag: str = "") -> dict:
    d = folder(clock_run, tag)
    lg = d / "beast.log"
    if not lg.exists():
        return {}
    s = json.loads((d / "settings.json").read_text(encoding="utf-8")) if (d / "settings.json").exists() else {}
    r = read_log(lg, float(s.get("burnin", 0.1)))           # a folder of combined runs has its burn-in removed already (burnin 0)
    P = r["params"]
    rate_key = "ucld.mean" if "ucld.mean" in P else "clock.rate"
    ages = {k[4:-1]: v for k, v in P.items() if k.startswith("age(") and k != "age(root)"}
    stems = {k[5:-1]: v for k, v in P.items() if k.startswith("stem(")}
    mono = {k[10:-1]: v["mean"] for k, v in P.items() if k.startswith("monophyly(")}              # the share of the samples in which the group is a clade
    return {"settings": s, "states": r["states"], "last_state": r.get("last_state"), "complete": r.get("last_state", 0) >= s.get("chain_length", 1), "root": P.get("age(root)"), "rate": P.get(rate_key),
            "rate_param": rate_key, "groups": ages, "stems": stems, "monophyly": mono, "partition": {k: P[k] for k in ("rate.apobec", "rate.other", "rate.ratio") if k in P}, "posterior": P.get("posterior"), "min_ess": min((v["ess"] for k, v in P.items() if k.startswith("age(") or k == rate_key), default=None)}


def mcc_tree(clock_run: dict, tag: str = "", kind: str = "hipstr", burnin: float = 0.1, force: bool = False) -> dict:
    """The summary tree of a finished chain, written as `mcc-<kind>.tree` (NEXUS, node heights as medians, with posterior and 95% HPD annotations) next to `beast.trees`,
    and `mcc-tips.tsv` (tip, date, lineage, clade, country, ...) to join onto the tips in a tree viewer such as PearTree. -> {tree, tips, trees, burnin_trees, newest_tip}"""
    pre = env_prefix()
    if not pre:
        raise RuntimeError(f"BEAST 환경이 없습니다: {status()['how']}")
    if kind not in ("hipstr", "mcc", "mrhipstr"):
        raise ValueError("요약 트리 종류는 hipstr, mcc, mrhipstr 중 하나입니다.")
    d = folder(clock_run, tag)
    src = d / "beast.trees"
    if not src.exists():
        raise ValueError("이 분석에는 beast.trees가 없습니다.")
    out = d / f"mcc-{kind}.tree"
    n = 0
    with open(src, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().lower().startswith("tree "):
                n += 1
    st = d / "settings.json"
    if st.exists():                                        # combined runs: the burn-in was removed when they were combined
        burnin = float(json.loads(st.read_text(encoding="utf-8")).get("burnin", burnin))
    skip = int(n * burnin)
    if force or not out.exists() or out.stat().st_mtime < src.stat().st_mtime:
        tmp = d / f".mcc-{kind}.tmp"
        for attempt in (1, 2):       # on a busy machine TreeAnnotator's worker threads sometimes die; a second try works
            r = subprocess.run([str(pre / "bin" / "treeannotator"), "-type", kind, "-heights", "median", "-burninTrees", str(skip), str(src), str(tmp)],
                               capture_output=True, text=True, timeout=1800)
            if r.returncode == 0 and tmp.exists():
                break
        if r.returncode != 0 or not tmp.exists():
            raise RuntimeError("TreeAnnotator가 실패했습니다: " + (r.stderr or r.stdout)[-400:])
        tmp.replace(out)
    tips = d / "mcc-tips.tsv"
    newest = None
    dates = d / "dates.tsv"
    if dates.exists():
        tips.write_text(dates.read_text(encoding="utf-8"), encoding="utf-8")
        try:
            rows = [ln.split("\t") for ln in dates.read_text(encoding="utf-8").splitlines()[1:] if ln.strip()]
            newest = max(float(x[1]) for x in rows)
        except (ValueError, IndexError):
            newest = None
    return {"tree": out.name, "tips": tips.name if tips.exists() else None, "trees": n, "burnin_trees": skip, "newest_tip": newest, "kind": kind}


def combine_runs(clock_run: dict, tags: list[str], out_tag: str, burnin: float = 0.1) -> dict:
    """Join independent runs of the same analysis (different seeds) with LogCombiner, removing `burnin` of each run first, into a new folder `out_tag`
    that is summarized like any other (its settings say `burnin: 0` and name the runs). -> the summary of the combined runs."""
    pre = env_prefix()
    if not pre:
        raise RuntimeError(f"BEAST 환경이 없습니다: {status()['how']}")
    if len(tags) < 2:
        raise ValueError("두 개 이상의 실행이 필요합니다.")
    src = [folder(clock_run, t) for t in tags]
    sets = [json.loads((d / "settings.json").read_text(encoding="utf-8")) for d in src]
    same = {k for k in sets[0] if k not in ("seed", "tag", "threads")}
    for k in same:
        if any(x.get(k) != sets[0][k] for x in sets[1:]):
            raise ValueError(f"실행들의 설정이 다릅니다 ({k}); 같은 분석의 다른 시드만 합칠 수 있습니다.")
    for d in src:
        if not (d / "beast.log").exists() or not (d / "beast.trees").exists():
            raise ValueError(f"{d.name}에 beast.log / beast.trees가 없습니다.")
    states = [int(read_log(d / "beast.log", 0).get("last_state") or 0) for d in src]
    if any(x < sets[0].get("chain_length", 0) for x in states):
        raise ValueError("끝나지 않은 실행은 합칠 수 없습니다.")
    cut = int(sets[0]["chain_length"] * burnin)
    out = folder(clock_run, out_tag)
    out.mkdir(parents=True, exist_ok=True)
    lc = str(pre / "bin" / "logcombiner")
    for args, name in ((["-burnin", str(cut)], "beast.log"), (["-trees", "-burnin", str(cut)], "beast.trees")):
        r = subprocess.run([lc, *args, *[str(d / name) for d in src], str(out / name)], capture_output=True, text=True, timeout=3600)
        if r.returncode != 0 or not (out / name).exists():
            raise RuntimeError("LogCombiner가 실패했습니다: " + (r.stderr or r.stdout)[-400:])
    for f in ("dates.tsv", "alignment.fasta", "subsample.tsv", "beast.xml"):          # the same XML for every run (only the seed differs)
        if (src[0] / f).exists():
            shutil.copyfile(src[0] / f, out / f)
    st = dict(sets[0], tag=out_tag, burnin=0, combined_from=tags, combined_burnin=burnin, seeds=[x.get("seed") for x in sets],
              chain_length=sets[0]["chain_length"])
    (out / "settings.json").write_text(json.dumps(st, ensure_ascii=False, indent=1), encoding="utf-8")
    res = summarize(clock_run, out_tag)
    res["complete"] = True
    (out / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    return res


def _short_country(c: str) -> str:
    return {"Democratic Republic of the Congo": "DR Congo", "United States of America": "USA", "United Kingdom": "UK", "Central African Republic": "CAR"}.get(c, c)


def _tip_label(name: str, a: dict) -> str:
    """`Country  date  lineage  id`: what a tip stands for, readable in the viewer; the country is first so the tips of one place sit together in the eye, and Korea is marked."""
    try:
        y = float(a.get("date", ""))
        yr = int(y)
        d = datetime(yr, 1, 1) + timedelta(days=round((y - yr) * (365 + (yr % 4 == 0))))
        date = d.strftime("%Y-%m-%d") if a.get("precision") == "day" else d.strftime("%Y-%m")
    except (TypeError, ValueError):
        date = "?"
    c = _short_country(a.get("country") or "")
    star = "\u2605 " if c == "South Korea" else ""
    return f"{star}{c or '?'}  {date}  {a.get('lineage') or ''}  {name}".replace('"', "")


def annotated_tree(clock_run: dict, tag: str = "", kind: str = "hipstr") -> Path:
    """`mcc-<kind>.annotated.tree`: the summary tree with a readable `label` (country, date, lineage, id) and date, lineage, clade and country on the tips (from dates.tsv), for the in-app PearTree viewer."""
    from . import peartree
    info = mcc_tree(clock_run, tag, kind)
    d = folder(clock_run, tag)
    src, out = d / info["tree"], d / f"mcc-{kind}.annotated.v2.tree"
    if out.exists() and out.stat().st_mtime >= src.stat().st_mtime and out.stat().st_mtime >= (d / "dates.tsv").stat().st_mtime:
        return out
    tips = {}
    lines = (d / "dates.tsv").read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    for ln in lines[1:]:
        c = ln.split("\t")
        if len(c) == len(head):
            row = dict(zip(head, c))
            a = {"label": _tip_label(c[0], row)}
            a.update({h: v for h, v in row.items() if h not in ("name", "precision")})
            tips[c[0]] = a
    out.write_text(peartree.annotate(src.read_text(encoding="utf-8"), tips), encoding="utf-8")
    return out


def result(clock_run: dict, tag: str = "") -> dict:
    """What the page shows: the saved summary, or - while a chain runs - the summary of the log so far."""
    d = folder(clock_run, tag)
    f = d / "summary.json"
    if f.exists():
        try:
            return {**json.loads(f.read_text(encoding="utf-8")), "running": False}
        except ValueError:
            pass
    s = summarize(clock_run, tag) if (d / "beast.log").exists() else {}
    return {**s, "running": bool(s)} if s else {}


def run_recorded(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """`run` with a row in `analysis_run` (tool `beast`, note `clock:<id>`): running / done / error, so that it shows with the other runs and a failed or interrupted chain is visible."""
    from . import datasets
    ds = datasets.get(clock_run["dataset_id"])
    st = status()
    s = parse_settings(settings)
    d = folder(clock_run, s["tag"])
    d.mkdir(parents=True, exist_ok=True)
    name = f"{clock_run['name']}-beast" + (f"-{s['tag']}" if s["tag"] else "")
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    members = nextstrain.run_members(clock_run)
    from datetime import datetime
    rid = runs.start(ds, "beast", name, datetime.now(), ds.get("build") or "", s, "beast-v1", members, clock_run["members_sha"] or "", clock_run.get("content_sha") or "", st["version"], d, d / "beast.screen.log",
                     f"clock:{clock_run['id']}", {"beast": st["version"], "settings": s})
    try:
        res = run(clock_run, s, progress)
    except Exception as e:
        runs.finish(rid, "beast", name, "error", str(e)[:400])
        raise
    runs.finish(rid, "beast", name, "done", "", "summary.json")
    return {**res, "run": rid}


# ------------------------------------------------------------------------------------------------ BETS: is there a temporal signal at all?
BETS_DEFAULTS = {"n_tips": 60, "clock": "relaxed", "tree_prior": "exponential", "initial_chain": 1000000, "path_steps": 50, "step_chain": 250000, "estimator": "both", "seed": 1, "threads": 0, "compare": False}
ESTIMATORS = ("both", "gss", "ps_ss")
_MLE_RE = {"ps": re.compile(r"log marginal likelihood \(using path sampling\)[^=]*=\s*(-?[0-9.eE+-]+)"), "ss": re.compile(r"log marginal likelihood \(using stepping stone sampling\)[^=]*=\s*(-?[0-9.eE+-]+)"),
           "gss": re.compile(r"log marginal likelihood \(using generalized stepping stone sampling\)[^=]*=\s*(-?[0-9.eE+-]+)")}
GSS_BURN = 2000           # states of the ordinary chain run before the path begins: kept short so that the chain stays where it started (the middle of the earlier posterior, inside the working prior)
GSS_PARAMS = ("kappa", "alpha", "clock.rate", "ucld.mean", "ucld.stdev", "demographic.popSize", "demographic.growthRate")


def parse_bets_settings(raw: dict | None) -> dict:
    s = {**BETS_DEFAULTS, **(raw or {})}
    try:
        out = {k: int(s[k]) for k in ("n_tips", "initial_chain", "path_steps", "step_chain", "seed", "threads")}
    except (TypeError, ValueError):
        raise ValueError("BETS 설정의 숫자가 올바르지 않습니다.")
    out["clock"], out["tree_prior"], out["estimator"] = str(s["clock"]), str(s["tree_prior"]), str(s["estimator"])
    out["compare"] = bool(s.get("compare"))
    out["name"] = re.sub(r"[^A-Za-z0-9_.\-]", "_", str(s.get("name") or ""))[:40]     # its own folder `bets/<name>/`, so that two estimations can run side by side
    if not 10 <= out["n_tips"] <= 300:
        raise ValueError("표본 유전체 수는 10~300입니다.")
    if out["initial_chain"] < 10000 or out["step_chain"] < 1000 or not 4 <= out["path_steps"] <= 500:
        raise ValueError("초기 체인은 10,000 이상, 단계당 체인은 1,000 이상, 경로 단계는 4~500이어야 합니다 (신뢰할 만한 추정은 단계 50 이상, 단계당 250,000 이상).")
    if out["clock"] not in CLOCKS or out["tree_prior"] not in PRIORS or out["estimator"] not in ESTIMATORS:
        raise ValueError(f"시계는 {', '.join(CLOCKS)}, 트리 사전분포는 {', '.join(PRIORS)}, 추정기는 {', '.join(ESTIMATORS)} 중 하나입니다.")
    return out


def bets_folder(clock_run: dict) -> Path:
    return Path(clock_run["run_dir"]) / "beast" / "bets"


def _root(clock_run: dict, s: dict) -> Path:
    """Where this BETS flavour keeps its analyses: `bets/` (dated vs undated) or `bets/models/` (the four-model comparison)."""
    return bets_folder(clock_run) / (s["name"] if s.get("name") else "models" if s.get("compare") else "")


def _models(s: dict) -> list[dict]:
    """The analyses of one run: BETS = the same model with real and with identical dates; comparison = the four clock x tree-prior models, all with the real dates (their marginal likelihoods are compared)."""
    if s.get("compare"):
        return [{"tag": f"{c}-{t}", "dated": True, "clock": c, "tree_prior": t} for c in CLOCKS for t in PRIORS]
    return [{"tag": "dated", "dated": True, "clock": s["clock"], "tree_prior": s["tree_prior"]}, {"tag": "undated", "dated": False, "clock": s["clock"], "tree_prior": s["tree_prior"]}]


def _model(s: dict, tag: str) -> dict:
    return next(m for m in _models(s) if m["tag"] == tag)


def _bets_data(clock_run: dict, s: dict) -> tuple[dict, dict, float, int, float]:
    """(sequences, real dates, starting rate, genomes, span in years) of the BETS subsample."""
    S = clock.summary(clock_run)
    if not S.get("tips"):
        raise ValueError("이 분자시계 실행에는 결과가 없습니다.")
    src = next((r for r in nextstrain.list_runs(clock_run["dataset_id"]) if r["id"] == S["source_run"]), None)
    if not src or src["tool"] != "squirrel":
        raise ValueError("BETS는 Squirrel 정렬이 필요합니다. 분자시계를 Squirrel 트리로 실행하세요.")
    chosen = choose_tips(S["tips"], s["n_tips"], s["seed"])
    if len(chosen) < 10:
        raise ValueError(f"BETS에 쓸 수 있는 유전체(날짜가 일 단위이고 이상치가 아님)가 {len(chosen)}개뿐입니다.")
    want = {t["name"] for t in chosen}
    seqs = {n: q.upper() for n, q in aligned.read_alignment(Path(src["run_dir"]) / "squirrel" / "alignment.aln.fasta.gz") if n in want}
    real = {t["name"]: t["date"] for t in chosen if t["name"] in seqs}
    rate0 = (S.get("time") or {}).get("rate") or (S.get("clock") or {}).get("rate") or 1e-5
    return seqs, real, rate0, len(real), max(real.values()) - min(real.values())


def bets_prepare(clock_run: dict, settings: dict | None = None) -> dict:
    """The XMLs of the BETS analyses: for each of 'dated' (the real sampling dates) and 'undated' (every genome given the same date, the isochronous model), the same model on the same genomes:
    * `psss/beast.xml` - an initial chain, then path sampling and stepping-stone sampling of the marginal likelihood (destination = the prior);
    * `mcmc/beast.xml` - an ordinary chain; `mcmc/gss.xml` is written when it has finished, from its posterior (see `run_bets`), for generalized stepping-stone sampling.
    Proper priors and fixed (empirical) base frequencies, since a Bayes factor needs proper priors and GSS has no working prior for the frequency simplex."""
    s = parse_bets_settings(settings)
    seqs, real, rate0, n, span = _bets_data(clock_run, s)
    same = {k: max(real.values()) for k in real}
    out = {}
    for m in _models(s):
        tag, dates = m["tag"], (real if m["dated"] else same)
        ms = {**s, "chain_length": s["initial_chain"], "samples": 1000, "clock": m["clock"], "tree_prior": m["tree_prior"]}
        base = _root(clock_run, s) / tag
        if s["estimator"] in ("both", "ps_ss"):
            (base / "psss").mkdir(parents=True, exist_ok=True)
            (base / "psss" / "beast.xml").write_text(build_xml(seqs, dates, {}, ms, rate0, "beast", mle={"chain_length": s["step_chain"], "path_steps": s["path_steps"], "alpha": 0.3}, proper=True, fixed_freq=True), encoding="utf-8")
        if s["estimator"] in ("both", "gss"):
            (base / "mcmc").mkdir(parents=True, exist_ok=True)
            (base / "mcmc" / "beast.xml").write_text(build_xml(seqs, dates, {}, ms, rate0, "beast", proper=True, fixed_freq=True), encoding="utf-8")
        out[tag] = str(base)
    _root(clock_run, s).mkdir(parents=True, exist_ok=True)
    (_root(clock_run, s) / "settings.json").write_text(json.dumps(s), encoding="utf-8")
    return {"n": n, "settings": s, "folders": out, "span_years": span}


def gss_xml(clock_run: dict, tag: str, s: dict) -> str:
    """The generalized stepping-stone XML of one BETS analysis, from the finished ordinary chain in `<tag>/mcmc/`: the working prior is made from that chain's posterior, and the new chain starts from its middle."""
    seqs, real, rate0, _n, _span = _bets_data(clock_run, s)
    m = _model(s, tag)
    dates = real if m["dated"] else {k: max(real.values()) for k in real}
    d = _root(clock_run, s) / tag / "mcmc"
    post = read_log(d / "beast.log", 0.5)["params"]
    init = {k: post[k]["median"] for k in GSS_PARAMS if k in post}
    flat = [k for k in GSS_PARAMS if k in post and post[k]["hi"] - post[k]["lo"] <= 0]
    if flat:                                                  # a parameter that never moved has no density to build a working prior from (BEAST would give NaN)
        raise RuntimeError(f"일반 체인에서 {', '.join(flat)} 값이 변하지 않았습니다 — 체인이 너무 짧거나 섞이지 않아 작업 사전분포를 만들 수 없습니다. 초기 체인을 늘려 다시 실행하세요.")
    ms = {**s, "chain_length": GSS_BURN, "samples": 100, "clock": m["clock"], "tree_prior": m["tree_prior"]}
    return build_xml(seqs, dates, {}, ms, rate0, "gss", proper=True, fixed_freq=True, init=init,
                     gss={"log": "beast.log", "burnin": s["initial_chain"] // 2, "chain_length": s["step_chain"], "path_steps": s["path_steps"], "alpha": 0.3})


def bets_verdict(log_bf: float | None) -> tuple[str, str]:
    """(level, text) for a log Bayes factor of 'sampling dates' over 'no dates' (Kass & Raftery's scale; BETS reads a log BF of 5 or more as strong support for a temporal signal)."""
    if log_bf is None:
        return "unknown", "추정하지 못했습니다."
    if log_bf >= 5:
        return "strong", f"로그 베이즈 인수 {log_bf:.1f} ≥ 5: 날짜가 모델을 뚜렷하게 개선합니다 — 시계 신호가 있다는 강한 증거입니다."
    if log_bf >= 3:
        return "positive", f"로그 베이즈 인수 {log_bf:.1f} (3~5): 시계 신호의 증거가 있지만 결정적이지는 않습니다."
    if log_bf >= 1:
        return "weak", f"로그 베이즈 인수 {log_bf:.1f} (1~3): 약한 증거입니다. 시계 신호가 있다고 하기에는 부족합니다."
    if log_bf > -1:
        return "none", f"로그 베이즈 인수 {log_bf:.1f}: 날짜가 모델을 개선하지 않습니다 — 시계 신호가 있다는 증거가 없습니다."
    return "against", f"로그 베이즈 인수 {log_bf:.1f} < -1: 날짜 없는 모델이 더 낫습니다 — 시계 신호가 없다는 쪽의 증거입니다."


def _log_states(path: Path) -> int:
    try:
        return int(float(path.read_text(encoding="utf-8").strip().splitlines()[-1].split("\t")[0]))
    except (OSError, ValueError, IndexError):
        return 0


def _mle_rows(path: Path) -> int:
    try:
        return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith(("#", "state")))
    except OSError:
        return 0


def run_bets(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """Prepare and run the BETS analyses (two models x the estimators asked for), then compare the marginal likelihoods. The analyses run side by side as separate BEAST processes: the path /
    stepping-stone runs and the ordinary chains start together, each generalized stepping-stone run starts when its ordinary chain has finished."""
    import os
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"BEAST 환경이 없습니다: {st['how']}")
    p = progress or (lambda *a: None)
    info = bets_prepare(clock_run, settings)
    s = info["settings"]
    root = _root(clock_run, s)
    tags = [m["tag"] for m in _models(s)]
    pre = env_prefix()
    use_psss, use_gss = s["estimator"] in ("both", "ps_ss"), s["estimator"] in ("both", "gss")
    jobs = []                                    # (tag, kind, folder, xml name, states)
    for tag in tags:
        if use_psss:
            jobs.append((tag, "psss", root / tag / "psss", "beast.xml", s["initial_chain"] + (s["path_steps"] + 1) * s["step_chain"]))
        if use_gss:
            jobs.append((tag, "mcmc", root / tag / "mcmc", "beast.xml", s["initial_chain"]))
            jobs.append((tag, "gss", root / tag / "mcmc", "gss.xml", GSS_BURN + (s["path_steps"] + 1) * s["step_chain"]))
    for _t, kind, d, xml, _n in jobs:
        for f in ("beast.log", "beast.trees", "beast.ops", "beast.mle.log", "gss.log", "gss.trees", "gss.ops", "gss.mle.log", "screen.txt", "gss.screen.txt"):
            if (kind == "gss") == f.startswith("gss") or kind == "psss" and not f.startswith("gss") or kind == "mcmc" and not f.startswith("gss"):
                (d / f).unlink(missing_ok=True)
    (root / "summary.json").unlink(missing_ok=True)
    n_threads = s["threads"] or max(2, (os.cpu_count() or 4) // (12 if s["compare"] else 4))
    total = sum(j[4] for j in jobs)
    procs: dict[tuple, tuple] = {}
    t0 = time.time()

    def start(tag, kind, d, xml):
        lf = open(d / ("gss.screen.txt" if kind == "gss" else "screen.txt"), "w", encoding="utf-8")
        procs[(tag, kind)] = (subprocess.Popen([str(pre / "bin" / "beast"), "-overwrite", "-seed", str(s["seed"]), "-threads", str(n_threads), xml], cwd=str(d), env=nextstrain.tool_env(pre),
                                               stdout=lf, stderr=subprocess.STDOUT, text=True), lf)

    def done(tag, kind, d, n):
        if (tag, kind) not in procs:
            return 0
        if kind == "mcmc":
            return min(n, _log_states(d / "beast.log"))
        if kind == "psss":
            rows = _mle_rows(d / "beast.mle.log")
            return min(n, s["initial_chain"] + rows * (s["step_chain"] // 100)) if rows else min(n, _log_states(d / "beast.log"))
        rows = _mle_rows(d / "gss.mle.log")
        return min(n, GSS_BURN + rows * (s["step_chain"] // 100)) if rows else min(n, _log_states(d / "gss.log"))
    try:
        for tag, kind, d, xml, _n in jobs:
            if kind in ("psss", "mcmc"):
                start(tag, kind, d, xml)
        pending_gss = [j for j in jobs if j[1] == "gss"]
        while True:
            time.sleep(10)
            for tag, kind, d, xml, n in list(pending_gss):
                pr = procs.get((tag, "mcmc"))
                if pr and pr[0].poll() is not None:
                    if pr[0].returncode != 0:
                        raise RuntimeError(f"BEAST({tag}, 일반 체인)이 실패했습니다: {(d / 'screen.txt').read_text(encoding='utf-8')[-300:]}")
                    (d / "gss.xml").write_text(gss_xml(clock_run, tag, s), encoding="utf-8")
                    start(tag, "gss", d, "gss.xml")
                    pending_gss.remove((tag, kind, d, xml, n))
            f = sum(done(t, k, d, n) for t, k, d, _x, n in jobs) / max(1, total)
            p(2 + 96 * f, f"BETS {f * 100:.0f}% ({time.time() - t0:.0f}초, 분석 {sum(1 for v in procs.values() if v[0].poll() is None)}개 실행 중)")
            if not pending_gss and all(v[0].poll() is not None for v in procs.values()):
                break
    finally:
        for pr, lf in procs.values():
            if pr.poll() is None:
                pr.terminate()
                try:
                    pr.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pr.kill()
            lf.close()
    res = {"settings": s, "seconds": round(time.time() - t0), "n_tips": info["n"], "span_years": info["span_years"], **{t: {} for t in tags}}
    for (tag, kind), (pr, _lf) in procs.items():
        d = root / tag / ("psss" if kind == "psss" else "mcmc")
        out = (d / ("gss.screen.txt" if kind == "gss" else "screen.txt")).read_text(encoding="utf-8")
        if pr.returncode != 0:
            raise RuntimeError(f"BEAST({tag}, {kind})가 실패했습니다 (종료 코드 {pr.returncode}): {out.strip()[-300:]}")
        if kind == "psss":
            for k in ("ps", "ss"):
                m = _MLE_RE[k].search(out)
                res[tag][k] = float(m[1]) if m else None
        elif kind == "gss":
            m = _MLE_RE["gss"].search(out)
            res[tag]["gss"] = float(m[1]) if m else None
    if s["compare"]:
        res["models"] = compare_models({t: res[t] for t in tags})
        (root / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
        p(100, "모델 비교 완료")
        return res
    for k in ("ps", "ss", "gss"):
        a, b = res["dated"].get(k), res["undated"].get(k)
        res["log_bf_" + k] = (a - b) if a is not None and b is not None else None
    best = next((res["log_bf_" + k] for k in ("gss", "ss", "ps") if res["log_bf_" + k] is not None), None)
    res["verdict"] = dict(zip(("level", "text"), bets_verdict(best)))
    res["verdict"]["estimator"] = next((k for k in ("gss", "ss", "ps") if res["log_bf_" + k] is not None), None)
    vals = [res["log_bf_" + k] for k in ("gss", "ss", "ps") if res["log_bf_" + k] is not None]
    if len(vals) > 1 and max(vals) - min(vals) > 2:
        res["verdict"]["text"] += " (추정기들의 값이 2 이상 다릅니다 — 경로 단계와 단계당 체인을 늘려 다시 돌리세요. 신뢰할 수 있는 결론이 아닙니다.)"
    (root / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    p(100, "BETS 완료")
    return res


def compare_models(res: dict) -> dict:
    """Rank the models by marginal likelihood (GSS if there is one, else SS, else PS) and give each one's log Bayes factor against the best; the estimators' disagreement is reported per model.
    Kass & Raftery: a log BF of 1-3 is weak, 3-5 positive, over 5 strong evidence for the better model."""
    est = next((k for k in ("gss", "ss", "ps") if all(r.get(k) is not None for r in res.values())), None)
    if est is None:
        return {"estimator": None, "ranking": [], "text": "한계우도를 추정하지 못했습니다."}
    ranking = sorted(({"model": t, "ml": r[est], "ps": r.get("ps"), "ss": r.get("ss"), "gss": r.get("gss")} for t, r in res.items()), key=lambda x: -x["ml"])
    best = ranking[0]["ml"]
    for x in ranking:
        x["log_bf_vs_best"] = round(x["ml"] - best, 2)
        vals = [v for v in (x["ps"], x["ss"], x["gss"]) if v is not None]
        x["estimator_spread"] = round(max(vals) - min(vals), 2) if len(vals) > 1 else None
    gap = ranking[0]["ml"] - ranking[1]["ml"] if len(ranking) > 1 else None
    level, _t = bets_verdict(gap)
    text = f"{ranking[0]['model']}이(가) 가장 높습니다(추정기 {est.upper()})."
    if gap is not None:
        text += f" 2위({ranking[1]['model']})와의 로그 베이즈 인수 {gap:.1f}: " + {"strong": "뚜렷한 차이", "positive": "차이의 증거가 있으나 결정적이지 않음", "weak": "약한 차이", "none": "차이가 없음", "against": "차이가 없음"}.get(level, "")
    if any((x["estimator_spread"] or 0) > 2 for x in ranking):
        text += " 추정기들의 값이 2 이상 다른 모델이 있습니다 — 경로 단계와 단계당 체인을 늘려 다시 돌리세요."
    return {"estimator": est, "ranking": ranking, "text": text}


def bets_result(clock_run: dict, compare: bool = False, name: str = "") -> dict:
    f = (bets_folder(clock_run) / (_safe(name) if name else "models" if compare else "") / "summary.json")
    try:
        return json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    except ValueError:
        return {}


def run_bets_recorded(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """`run_bets` with a row in `analysis_run` (tool `bets`, note `clock:<id>`)."""
    from . import datasets
    from datetime import datetime
    ds = datasets.get(clock_run["dataset_id"])
    st = status()
    s = parse_bets_settings(settings)
    name = f"{clock_run['name']}-{'models' if s['compare'] else 'bets'}" + (f"-{s['name']}" if s["name"] else "")
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    d = bets_folder(clock_run)
    d.mkdir(parents=True, exist_ok=True)
    rid = runs.start(ds, "bets", name, datetime.now(), ds.get("build") or "", s, "bets-v2", nextstrain.run_members(clock_run), clock_run["members_sha"] or "", clock_run.get("content_sha") or "", st["version"], d,
                     _root(clock_run, s) / _models(s)[0]["tag"] / "psss" / "screen.txt", f"clock:{clock_run['id']}" + (":models" if s["compare"] else ""), {"beast": st["version"], "settings": s})
    try:
        res = run_bets(clock_run, s, progress)
    except Exception as e:
        runs.finish(rid, "bets", name, "error", str(e)[:400])
        raise
    runs.finish(rid, "bets", name, "done", "", "summary.json")
    return {**res, "run": rid}
