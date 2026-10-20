"""BEAST for the tMRCA: a package prepared from a clock run (`clock.py`), an optional run, and the posterior read back.

TreeTime (clock.py) is fast and shows whether there is a clock at all, but its intervals are conditional on one tree and a fitted rate. For a tMRCA that is meant to be cited, a Bayesian analysis
(BEAST: tip dates, a clock model, a tree prior, all uncertainty sampled together) is the standard. BEAST cannot take a thousand 197 kb genomes in reasonable time, so the package holds a
**subsample** chosen to keep the whole sampling period and every lineage (dated by day, not a clock outlier, not flagged by QC), the alignment of those genomes from the Squirrel alignment, tip dates, and a
BEAST 1.10-style XML (HKY + gamma sites; strict or uncorrelated log-normal relaxed clock; constant-size or exponential coalescent) that BEAST 10 reads. Nothing is hidden in the XML: the model is written out
and can be opened and changed in BEAUti / a text editor.

Reading the result: the tMRCA of the root and of the lineages as dates (mean / median / 95% HPD), the clock rate, and the effective sample size (ESS) of each - an ESS below 200 means the chain has not
run long enough and the interval must not be trusted. The package is a starting point: a real analysis needs several chains, checks in Tracer, and a model comparison."""
import json
import math
import random
import re
import subprocess
import time
from datetime import date
from pathlib import Path

from . import aligned, clock, log, nextstrain, runs

_log = log.get("beast")
DEFAULTS = {"n_tips": 150, "chain_length": 30000000, "samples": 10000, "clock": "relaxed", "tree_prior": "exponential", "seed": 1, "threads": 0}
CLOCKS = ("strict", "relaxed")
PRIORS = ("constant", "exponential")


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
    if out["clock"] not in CLOCKS or out["tree_prior"] not in PRIORS:
        raise ValueError(f"시계는 {', '.join(CLOCKS)}, 트리 사전분포는 {', '.join(PRIORS)} 중 하나입니다.")
    return out


# ------------------------------------------------------------------------------------------------ the subsample
def choose_tips(tips: list[dict], n: int, seed: int = 1) -> list[dict]:
    """`n` of the clock run's tips, keeping the whole sampling period and the lineages: only genomes dated by day and not clock outliers; one genome from every month first (the earliest-dated of a
    random order), then more rounds from the busiest months and lineages until `n` are taken. Deterministic for a seed."""
    rng = random.Random(seed)
    pool = [t for t in tips if not t["outlier"] and t["precision"] == "day"]
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


def build_xml(seqs: dict[str, str], dates: dict[str, float], groups: dict[str, list[str]], s: dict, rate0: float, fileprefix: str = "beast", mle: dict | None = None, proper: bool = False) -> str:
    """The BEAST 1.10-style XML. `proper`: only proper priors (the clock rate and the population size get log-normal priors instead of the vague gamma / 1/x ones): a marginal likelihood, and so a Bayes factor, exists only for proper priors. `mle` = {"chain_length", "path_steps", "alpha"}: after the initial chain, path sampling and stepping-stone sampling of the marginal likelihood (`marginalLikelihoodEstimator`, written to `<prefix>.mle.log`). `seqs` {taxon: aligned sequence}, `dates` {taxon: decimal year}, `groups` {name: [taxa]} for the tMRCA of each group, `rate0` the starting clock rate (substitutions / site / year)."""
    names = list(seqs)
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
    w('\t<alignment id="alignment" dataType="nucleotide">')
    for n in names:
        w(f'\t\t<sequence><taxon idref="{_safe(n)}"/>{seqs[n]}</sequence>')
    w('\t</alignment>')
    w('\t<patterns id="patterns" from="1" strip="false"><alignment idref="alignment"/></patterns>')
    if s["tree_prior"] == "constant":
        w('\t<constantSize id="demographic" units="years"><populationSize><parameter id="demographic.popSize" value="1.0" lower="0.0"/></populationSize></constantSize>')
    else:
        w('\t<exponentialGrowth id="demographic" units="years"><populationSize><parameter id="demographic.popSize" value="1.0" lower="0.0"/></populationSize><growthRate><parameter id="demographic.growthRate" value="1.0"/></growthRate></exponentialGrowth>')
    w('\t<coalescentSimulator id="startingTree"><taxa idref="taxa"/><' + ("constantSize" if s["tree_prior"] == "constant" else "exponentialGrowth") + ' idref="demographic"/></coalescentSimulator>')
    w('\t<treeModel id="treeModel"><coalescentTree idref="startingTree"/><rootHeight><parameter id="treeModel.rootHeight"/></rootHeight>'
      '<nodeHeights internalNodes="true"><parameter id="treeModel.internalNodeHeights"/></nodeHeights><nodeHeights internalNodes="true" rootNode="true"><parameter id="treeModel.allInternalNodeHeights"/></nodeHeights></treeModel>')
    w('\t<treeLengthStatistic id="treeLength"><treeModel idref="treeModel"/></treeLengthStatistic>')
    w('\t<tmrcaStatistic id="age(root)" absolute="true"><treeModel idref="treeModel"/><mrca><taxa idref="taxa"/></mrca></tmrcaStatistic>')
    for g in groups:
        w(f'\t<tmrcaStatistic id="age({_safe(g)})" absolute="true"><treeModel idref="treeModel"/><mrca><taxa idref="group.{_safe(g)}"/></mrca></tmrcaStatistic>')
    w('\t<coalescentLikelihood id="coalescent"><model><' + ("constantSize" if s["tree_prior"] == "constant" else "exponentialGrowth") + ' idref="demographic"/></model><populationTree><treeModel idref="treeModel"/></populationTree></coalescentLikelihood>')
    if s["clock"] == "strict":
        w(f'\t<strictClockBranchRates id="branchRates"><rate><parameter id="clock.rate" value="{rate0:.4e}" lower="0.0"/></rate></strictClockBranchRates>')
        rate_param, rate_ref = "clock.rate", '<strictClockBranchRates idref="branchRates"/>'
    else:
        w(f'\t<discretizedBranchRates id="branchRates"><treeModel idref="treeModel"/><distribution><logNormalDistributionModel meanInRealSpace="true"><mean><parameter id="ucld.mean" value="{rate0:.4e}" lower="0.0"/></mean>'
          '<stdev><parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0"/></stdev></logNormalDistributionModel></distribution><rateCategories><parameter id="branchRates.categories"/></rateCategories></discretizedBranchRates>')
        rate_param, rate_ref = "ucld.mean", '<discretizedBranchRates idref="branchRates"/>'
    w(f'\t<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true"><treeModel idref="treeModel"/>{rate_ref}</rateStatistic>')
    w('\t<HKYModel id="hky"><frequencies><frequencyModel dataType="nucleotide"><frequencies><parameter id="frequencies" value="0.25 0.25 0.25 0.25"/></frequencies></frequencyModel></frequencies>'
      '<kappa><parameter id="kappa" value="2.0" lower="0.0"/></kappa></HKYModel>')
    w('\t<siteModel id="siteModel"><substitutionModel><HKYModel idref="hky"/></substitutionModel><gammaShape gammaCategories="4"><parameter id="alpha" value="0.5" lower="0.0"/></gammaShape></siteModel>')
    w(f'\t<treeDataLikelihood id="treeLikelihood" useAmbiguities="false"><partition><patterns idref="patterns"/><siteModel idref="siteModel"/></partition><treeModel idref="treeModel"/>{rate_ref}</treeDataLikelihood>')
    w('\t<operators id="operators" optimizationSchedule="default">')
    w('\t\t<scaleOperator scaleFactor="0.75" weight="0.1"><parameter idref="kappa"/></scaleOperator>')
    w('\t\t<scaleOperator scaleFactor="0.75" weight="0.1"><parameter idref="alpha"/></scaleOperator>')
    w('\t\t<deltaExchange delta="0.01" weight="0.1"><parameter idref="frequencies"/></deltaExchange>')
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
    w('\t\t<scaleOperator scaleFactor="0.75" weight="3"><parameter idref="demographic.popSize"/></scaleOperator>')
    if s["tree_prior"] == "exponential":
        w('\t\t<randomWalkOperator windowSize="1.0" weight="3"><parameter idref="demographic.growthRate"/></randomWalkOperator>')
    w('\t</operators>')
    w(f'\t<mcmc id="mcmc" chainLength="{chain}" autoOptimize="true" operatorAnalysis="{fileprefix}.ops">')
    w('\t\t<posterior id="posterior"><prior id="prior">')
    w('\t\t\t<logNormalPrior mean="1.0" stdev="1.25" offset="0.0" meanInRealSpace="false"><parameter idref="kappa"/></logNormalPrior>')
    w('\t\t\t<exponentialPrior mean="0.5" offset="0.0"><parameter idref="alpha"/></exponentialPrior>')
    w('\t\t\t<dirichletPrior alpha="1.0" sumsTo="1.0"><parameter idref="frequencies"/></dirichletPrior>')
    if proper:
        w(f'\t\t\t<logNormalPrior mean="-11.5" stdev="2.0" offset="0.0" meanInRealSpace="false"><parameter idref="{rate_param}"/></logNormalPrior>')        # median 1e-5, four orders of magnitude wide
    else:
        w(f'\t\t\t<gammaPrior shape="0.001" scale="1000.0" offset="0.0"><parameter idref="{rate_param}"/></gammaPrior>')
    if s["clock"] == "relaxed":
        w('\t\t\t<exponentialPrior mean="0.3333333333333333" offset="0.0"><parameter idref="ucld.stdev"/></exponentialPrior>')
    if proper:
        w('\t\t\t<logNormalPrior mean="0.0" stdev="2.0" offset="0.0" meanInRealSpace="false"><parameter idref="demographic.popSize"/></logNormalPrior>')
    else:
        w('\t\t\t<oneOnXPrior><parameter idref="demographic.popSize"/></oneOnXPrior>')
    if s["tree_prior"] == "exponential":
        w('\t\t\t<laplacePrior mean="0.0" scale="100.0"><parameter idref="demographic.growthRate"/></laplacePrior>')
    w('\t\t\t<coalescentLikelihood idref="coalescent"/>')
    w('\t\t</prior><likelihood id="likelihood"><treeDataLikelihood idref="treeLikelihood"/></likelihood></posterior>')
    w('\t\t<operators idref="operators"/>')
    w(f'\t\t<log id="screenLog" logEvery="{every * 10}"><column label="Posterior" dp="4" width="12"><posterior idref="posterior"/></column><column label="Prior" dp="4" width="12"><prior idref="prior"/></column>'
      '<column label="Likelihood" dp="4" width="12"><likelihood idref="likelihood"/></column><column label="rootAge" sf="6" width="12"><tmrcaStatistic idref="age(root)"/></column></log>')
    w(f'\t\t<log id="fileLog" logEvery="{every}" fileName="{fileprefix}.log" overwrite="false"><posterior idref="posterior"/><prior idref="prior"/><likelihood idref="likelihood"/>'
      '<parameter idref="treeModel.rootHeight"/><tmrcaStatistic idref="age(root)"/>' + "".join(f'<tmrcaStatistic idref="age({_safe(g)})"/>' for g in groups) + '<treeLengthStatistic idref="treeLength"/>'
      '<parameter idref="demographic.popSize"/>' + ('<parameter idref="demographic.growthRate"/>' if s["tree_prior"] == "exponential" else "") + '<parameter idref="kappa"/><parameter idref="alpha"/>'
      f'<parameter idref="{rate_param}"/>' + ('<parameter idref="ucld.stdev"/>' if s["clock"] == "relaxed" else "") + '<rateStatistic idref="meanRate"/><treeDataLikelihood idref="treeLikelihood"/><coalescentLikelihood idref="coalescent"/></log>')
    w(f'\t\t<logTree id="treeFileLog" logEvery="{every}" nexusFormat="true" fileName="{fileprefix}.trees" sortTranslationTable="true"><treeModel idref="treeModel"/>'
      f'<trait name="rate" tag="rate">{rate_ref}</trait><posterior idref="posterior"/></logTree>')
    w('\t</mcmc>')
    w('\t<report><property name="timer"><mcmc idref="mcmc"/></property></report>')
    if mle:
        w(f'\t<marginalLikelihoodEstimator id="MLE" chainLength="{mle["chain_length"]}" pathSteps="{mle["path_steps"]}" pathScheme="betaquantile" alpha="{mle.get("alpha", 0.3)}">')
        w('\t\t<samplers><mcmc idref="mcmc"/></samplers>')
        w('\t\t<pathLikelihood id="pathLikelihood"><source><posterior idref="posterior"/></source><destination><prior idref="prior"/></destination></pathLikelihood>')
        w(f'\t\t<log id="MLELog" logEvery="{max(1, mle["chain_length"] // 100)}" fileName="{fileprefix}.mle.log"><pathLikelihood idref="pathLikelihood"/></log>')
        w('\t</marginalLikelihoodEstimator>')
        for el in ("pathSamplingAnalysis", "steppingStoneSamplingAnalysis"):
            w(f'\t<{el} fileName="{fileprefix}.mle.log"><likelihoodColumn name="pathLikelihood.delta"/><thetaColumn name="pathLikelihood.theta"/></{el}>')
    w('</beast>')
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------------------------------------ package, run, results
def folder(clock_run: dict) -> Path:
    return Path(clock_run["run_dir"]) / "beast"


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
    chosen = choose_tips(S["tips"], s["n_tips"], s["seed"])
    if len(chosen) < 10:
        raise ValueError(f"BEAST에 쓸 수 있는 유전체(날짜가 일 단위이고 이상치가 아님)가 {len(chosen)}개뿐입니다.")
    want = {t["name"] for t in chosen}
    p(15, "정렬에서 서열을 읽는 중…")
    seqs = {}
    for name, seq in aligned.read_alignment(Path(src["run_dir"]) / "squirrel" / "alignment.aln.fasta.gz"):
        if name in want:
            seqs[name] = seq.upper()
    chosen = [t for t in chosen if t["name"] in seqs]
    dates = {t["name"]: t["date"] for t in chosen}
    groups: dict[str, list[str]] = {}
    for t in chosen:
        for kind, val in (("clade", t["clade"]), ("lineage", t["lineage"])):
            if val:
                groups.setdefault(f"{kind}.{val}", []).append(t["name"])
    groups = {g: m for g, m in groups.items() if 5 <= len(m) < len(chosen)}
    d = folder(clock_run)
    d.mkdir(parents=True, exist_ok=True)
    rate0 = (S.get("time") or {}).get("rate") or (S.get("clock") or {}).get("rate") or 1e-5
    p(40, "XML을 쓰는 중…")
    (d / "beast.xml").write_text(build_xml(seqs, dates, groups, s, rate0), encoding="utf-8")
    with open(d / "alignment.fasta", "w", encoding="utf-8") as f:
        for n, q in seqs.items():
            f.write(f">{n}\n{q}\n")
    (d / "dates.tsv").write_text("name\tdate\tlineage\tclade\tcountry\n" + "".join(f"{t['name']}\t{t['date']:.5f}\t{t['lineage']}\t{t['clade']}\t{t['country']}\n" for t in chosen), encoding="utf-8")
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
    return {"states": len(rows), "last_state": int(float(rows[-1][0])), "burnin": burnin, "params": out}


def run(clock_run: dict, settings: dict | None = None, progress=None, threads: int | None = None) -> dict:
    """Prepare the package and run BEAST on it (a long job). The result is written as `beast/summary.json` when the chain ends; the log can be read while it runs (`read_log`)."""
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"BEAST 환경이 없습니다: {st['how']}")
    p = progress or (lambda *a: None)
    info = prepare(clock_run, settings, lambda pct, m="": p(pct * 0.02, m))
    s = info["settings"]
    d = folder(clock_run)
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
    res = summarize(clock_run)
    res["seconds"] = round(time.time() - t0)
    (d / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    p(100, "BEAST 완료")
    return res


def summarize(clock_run: dict) -> dict:
    d = folder(clock_run)
    lg = d / "beast.log"
    if not lg.exists():
        return {}
    r = read_log(lg)
    s = json.loads((d / "settings.json").read_text(encoding="utf-8")) if (d / "settings.json").exists() else {}
    P = r["params"]
    rate_key = "ucld.mean" if "ucld.mean" in P else "clock.rate"
    ages = {k[4:-1]: v for k, v in P.items() if k.startswith("age(") and k != "age(root)"}
    return {"settings": s, "states": r["states"], "last_state": r.get("last_state"), "complete": r.get("last_state", 0) >= s.get("chain_length", 1), "root": P.get("age(root)"), "rate": P.get(rate_key),
            "rate_param": rate_key, "groups": ages, "posterior": P.get("posterior"), "min_ess": min((v["ess"] for k, v in P.items() if k.startswith("age(") or k == rate_key), default=None)}


def result(clock_run: dict) -> dict:
    """What the page shows: the saved summary, or - while a chain runs - the summary of the log so far."""
    d = folder(clock_run)
    f = d / "summary.json"
    if f.exists():
        try:
            return {**json.loads(f.read_text(encoding="utf-8")), "running": False}
        except ValueError:
            pass
    s = summarize(clock_run) if (d / "beast.log").exists() else {}
    return {**s, "running": bool(s)} if s else {}


def run_recorded(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """`run` with a row in `analysis_run` (tool `beast`, note `clock:<id>`): running / done / error, so that it shows with the other runs and a failed or interrupted chain is visible."""
    from . import datasets
    ds = datasets.get(clock_run["dataset_id"])
    st = status()
    s = parse_settings(settings)
    d = folder(clock_run)
    d.mkdir(parents=True, exist_ok=True)
    name = f"{clock_run['name']}-beast"
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
BETS_DEFAULTS = {"n_tips": 60, "clock": "relaxed", "tree_prior": "exponential", "initial_chain": 1000000, "path_steps": 50, "step_chain": 250000, "seed": 1, "threads": 0}
_MLE_RE = {"ps": re.compile(r"log marginal likelihood \(using path sampling\)[^=]*=\s*(-?[0-9.eE+-]+)"), "ss": re.compile(r"log marginal likelihood \(using stepping stone sampling\)[^=]*=\s*(-?[0-9.eE+-]+)")}


def parse_bets_settings(raw: dict | None) -> dict:
    s = {**BETS_DEFAULTS, **(raw or {})}
    try:
        out = {k: int(s[k]) for k in ("n_tips", "initial_chain", "path_steps", "step_chain", "seed", "threads")}
    except (TypeError, ValueError):
        raise ValueError("BETS 설정의 숫자가 올바르지 않습니다.")
    out["clock"], out["tree_prior"] = str(s["clock"]), str(s["tree_prior"])
    if not 10 <= out["n_tips"] <= 300:
        raise ValueError("표본 유전체 수는 10~300입니다.")
    if out["initial_chain"] < 10000 or out["step_chain"] < 1000 or not 4 <= out["path_steps"] <= 500:
        raise ValueError("초기 체인은 10,000 이상, 단계당 체인은 1,000 이상, 경로 단계는 4~500이어야 합니다 (신뢰할 만한 추정은 단계 50 이상, 단계당 250,000 이상).")
    if out["clock"] not in CLOCKS or out["tree_prior"] not in PRIORS:
        raise ValueError(f"시계는 {', '.join(CLOCKS)}, 트리 사전분포는 {', '.join(PRIORS)} 중 하나입니다.")
    return out


def bets_folder(clock_run: dict) -> Path:
    return Path(clock_run["run_dir"]) / "beast" / "bets"


def bets_prepare(clock_run: dict, settings: dict | None = None) -> dict:
    """Two XMLs of the same model on the same genomes - one with the real sampling dates, one with every genome given the same date (the isochronous model) - each with path-sampling and
    stepping-stone estimation of the marginal likelihood after an initial chain. Proper priors, since a Bayes factor needs them."""
    s = parse_bets_settings(settings)
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
    chosen = [t for t in chosen if t["name"] in seqs]
    real = {t["name"]: t["date"] for t in chosen}
    same = {n: max(real.values()) for n in real}
    rate0 = (S.get("time") or {}).get("rate") or (S.get("clock") or {}).get("rate") or 1e-5
    ms = {**s, "chain_length": s["initial_chain"], "samples": 1000}
    out = {}
    for tag, dates in (("dated", real), ("undated", same)):
        d = bets_folder(clock_run) / tag
        d.mkdir(parents=True, exist_ok=True)
        (d / "beast.xml").write_text(build_xml(seqs, dates, {}, ms, rate0, "beast", mle={"chain_length": s["step_chain"], "path_steps": s["path_steps"], "alpha": 0.3}, proper=True), encoding="utf-8")
        out[tag] = str(d)
    (bets_folder(clock_run) / "settings.json").write_text(json.dumps(s), encoding="utf-8")
    return {"n": len(chosen), "settings": s, "folders": out, "span_years": max(real.values()) - min(real.values())}


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


def run_bets(clock_run: dict, settings: dict | None = None, progress=None) -> dict:
    """Prepare and run the two analyses side by side (two BEAST processes), then compare their marginal likelihoods. Long: initial chain + (path steps + 1) x chain per step, for each."""
    import os
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"BEAST 환경이 없습니다: {st['how']}")
    p = progress or (lambda *a: None)
    info = bets_prepare(clock_run, settings)
    s = info["settings"]
    root = bets_folder(clock_run)
    for tag in ("dated", "undated"):
        for f in ("beast.log", "beast.trees", "beast.ops", "beast.mle.log", "screen.txt"):
            (root / tag / f).unlink(missing_ok=True)
    (root / "summary.json").unlink(missing_ok=True)
    pre = env_prefix()
    n_threads = s["threads"] or max(2, (os.cpu_count() or 4) // 2)
    total_states = s["initial_chain"] + (s["path_steps"] + 1) * s["step_chain"]
    procs = {}
    t0 = time.time()
    for tag in ("dated", "undated"):
        lf = open(root / tag / "screen.txt", "w", encoding="utf-8")
        procs[tag] = (subprocess.Popen([str(pre / "bin" / "beast"), "-overwrite", "-seed", str(s["seed"]), "-threads", str(n_threads), "beast.xml"], cwd=str(root / tag),
                                       env=nextstrain.tool_env(pre), stdout=lf, stderr=subprocess.STDOUT, text=True), lf)

    def done_states(tag):
        d = root / tag
        n = 0
        lg = d / "beast.log"
        if lg.exists():
            try:
                n = int(float(lg.read_text(encoding="utf-8").strip().splitlines()[-1].split("\t")[0]))
            except (ValueError, IndexError):
                n = 0
        ml = d / "beast.mle.log"
        if ml.exists():
            rows = sum(1 for line in ml.read_text(encoding="utf-8").splitlines() if line and not line.startswith(("#", "state")))
            n = s["initial_chain"] + rows * (s["step_chain"] // 100)
        return min(n, total_states)
    try:
        while any(pr.poll() is None for pr, _lf in procs.values()):
            time.sleep(10)
            f = sum(done_states(t) for t in procs) / (2 * total_states)
            p(2 + 96 * f, f"BETS {f * 100:.0f}% ({time.time() - t0:.0f}초, 날짜 있음·없음 두 분석을 동시에)")
    finally:
        for pr, lf in procs.values():
            if pr.poll() is None:
                pr.terminate()
                try:
                    pr.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    pr.kill()
            lf.close()
    res = {"settings": s, "seconds": round(time.time() - t0), "n_tips": info["n"], "span_years": info["span_years"]}
    for tag, (pr, _lf) in procs.items():
        out = (root / tag / "screen.txt").read_text(encoding="utf-8")
        if pr.returncode != 0:
            raise RuntimeError(f"BEAST({tag})가 실패했습니다 (종료 코드 {pr.returncode}): {out.strip()[-300:]}")
        res[tag] = {k: (float(m[1]) if (m := rx.search(out)) else None) for k, rx in _MLE_RE.items()}
    for k in ("ps", "ss"):
        a, b = res["dated"][k], res["undated"][k]
        res["log_bf_" + k] = (a - b) if a is not None and b is not None else None
    res["verdict"] = dict(zip(("level", "text"), bets_verdict(res["log_bf_ss"] if res["log_bf_ss"] is not None else res["log_bf_ps"])))
    if res["log_bf_ps"] is not None and res["log_bf_ss"] is not None and abs(res["log_bf_ps"] - res["log_bf_ss"]) > 2:
        res["verdict"]["text"] += " (경로 표집과 디딤돌 표집의 값이 2 이상 다릅니다 — 경로 단계를 늘려 다시 돌리세요.)"
    (root / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
    p(100, "BETS 완료")
    return res


def bets_result(clock_run: dict) -> dict:
    f = bets_folder(clock_run) / "summary.json"
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
    name = f"{clock_run['name']}-bets"
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    d = bets_folder(clock_run)
    d.mkdir(parents=True, exist_ok=True)
    rid = runs.start(ds, "bets", name, datetime.now(), ds.get("build") or "", s, "bets-v1", nextstrain.run_members(clock_run), clock_run["members_sha"] or "", clock_run.get("content_sha") or "", st["version"], d,
                     d / "dated" / "screen.txt", f"clock:{clock_run['id']}", {"beast": st["version"], "settings": s})
    try:
        res = run_bets(clock_run, s, progress)
    except Exception as e:
        runs.finish(rid, "bets", name, "error", str(e)[:400])
        raise
    runs.finish(rid, "bets", name, "done", "", "summary.json")
    return {**res, "run": rid}
