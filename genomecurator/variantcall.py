"""A second, independent opinion on SNVs: snippy on the genome's FASTA, then lofreq on snippy's BAM (the way the Lee/Choi paper called them).

`snippy --ctgs` cuts each assembled genome into simulated 20x reads, maps them to the reference (BWA) and calls with FreeBayes; its BAM goes to `lofreq call`. The two callers use
different statistics, and neither depends on the alignment Nextclade / Squirrel make, so a variant that all of them find is well supported and one that only some find is a variant to look at -
nothing is removed or masked here. For an assembly the "reads" are simulated: depth and allele frequency say nothing about the original sample (AF is 0 or 1); real low-frequency calls need real reads.

snippy and lofreq each live in their own conda environment (`gc-snippy`, `gc-lofreq`; snippy 4.6.0 cannot be installed together with lofreq), both GPL and run as separate programs."""
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from collections import Counter

from . import flowlogic, log, mutprofile, nextstrain

_log = log.get("variantcall")
ENVS = {"snippy": ("GC_SNIPPY_ENV", "gc-snippy"), "lofreq": ("GC_LOFREQ_ENV", "gc-lofreq")}


def env_prefix(tool: str) -> Path | None:
    var, default = ENVS[tool]
    return nextstrain.find_env(os.environ.get(var, default), tool)


def status() -> dict:
    out = {}
    for t in ENVS:
        pre = env_prefix(t)
        out[t] = {"ok": bool(pre), "env": ENVS[t][1], "prefix": str(pre) if pre else ""}
    out["ok"] = all(v["ok"] for v in out.values() if isinstance(v, dict))
    return out


def tool_versions() -> dict:
    """{"snippy": "4.6.0", "lofreq": "2.1.5"} as the tools report them (empty when a tool cannot be run)."""
    out = {}
    for t, args in (("snippy", ["--version"]), ("lofreq", ["version"])):
        pre = env_prefix(t)
        v = ""
        if pre:
            try:
                r = subprocess.run([str(pre / "bin" / t)] + args, capture_output=True, text=True, timeout=60, env=nextstrain.tool_env(pre))
                text = (r.stdout + r.stderr)
                m = re.search(r"(\d+\.\d+(?:\.\d+)?)", text)
                v = m.group(1) if m else ""
            except Exception:
                v = ""
        out[t] = v
    return out


def _run(tool: str, cmd: list[str], cwd: Path, timeout: int = 1800):
    pre = env_prefix(tool)
    if not pre:
        raise RuntimeError(f"{tool} 환경({ENVS[tool][1]})이 없습니다: mamba create -n {ENVS[tool][1]} -c conda-forge -c bioconda {tool}")
    r = subprocess.run([str(pre / "bin" / cmd[0])] + cmd[1:], cwd=str(cwd), env=nextstrain.tool_env(pre), capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"{cmd[0]} 실패: {(r.stderr or r.stdout).strip()[-400:]}")
    return r


def parse_vcf(path: Path) -> dict:
    """{(pos, ref, alt): {"dp": depth, "af": allele frequency or None}} of the single-base substitutions of a VCF (indels and multi-base changes are not SNVs and are left out)."""
    out = {}
    for line in open(path, encoding="utf-8"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 8 or len(p[3]) != 1 or len(p[4]) != 1 or p[4] == ".":
            continue
        info = dict(kv.split("=", 1) for kv in p[7].split(";") if "=" in kv)
        out[(int(p[1]), p[3].upper(), p[4].upper())] = {"dp": int(info["DP"]) if info.get("DP", "").isdigit() else None, "af": float(info["AF"]) if "AF" in info else None}
    return out


def call_genome(name: str, seq: str, ref_fa: Path, work: Path, cpus: int = 2) -> dict:
    """{"snippy": {snv: info}, "lofreq": {snv: info}} for one genome. Work files go to `work` (a scratch directory)."""
    d = work / name
    d.mkdir(parents=True, exist_ok=True)
    fa = d / "genome.fa"
    fa.write_text(f">{name}\n{seq}\n", encoding="utf-8")
    _run("snippy", ["snippy", "--cpus", str(cpus), "--outdir", str(d / "s"), "--ref", str(Path(ref_fa).resolve()), "--ctgs", str(fa), "--force"], d)
    sn = parse_vcf(d / "s" / "snps.vcf")
    _run("lofreq", ["lofreq", "call", "-f", str(d / "s" / "reference" / "ref.fa"), "-o", str(d / "lofreq.vcf"), str(d / "s" / "snps.bam")], d)
    lf = parse_vcf(d / "lofreq.vcf")
    return {"snippy": sn, "lofreq": lf}


def call_many(fasta: Path, ref_fa: Path, ids: set | None = None, limit: int | None = None, workers: int = 4, progress=None) -> dict:
    """{genome id: {"snippy": ..., "lofreq": ...}} for the genomes of a FASTA file (all, those in `ids`, or the first `limit`)."""
    recs = []
    for name, seq in mutprofile.iter_fasta(fasta) if hasattr(mutprofile, "iter_fasta") else _iter(fasta):
        if ids is not None and name not in ids:
            continue
        recs.append((name, seq))
        if limit and len(recs) >= limit:
            break
    out, n = {}, 0
    work = Path(tempfile.mkdtemp(prefix="gc-variantcall-"))
    try:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for (name, _s), res in zip(recs, ex.map(lambda r: call_genome(r[0], r[1], ref_fa, work), recs)):
                out[name] = res
                n += 1
                if progress:
                    progress(n, len(recs))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return out


def _iter(fasta: Path):
    name, parts = None, []
    for line in open(fasta, encoding="utf-8"):
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts)


concordance = flowlogic.concordance                         # one implementation, shared with the Snakemake flow (flowlogic.py)


# ------------------------------------------------------------------------------------------------ on a finished mutation-profile run
def _vc_file(run: dict, source: str = "direct") -> Path:
    """The call table of a run: `direct` = this module (`mutprofile/`), `flow` = the Snakemake flow (`variantflow/`, see variantflow.py)."""
    return Path(run["run_dir"]) / ("variantflow" if source == "flow" else "mutprofile") / "variantcall.tsv"


def run_for_mutprofile(run: dict, workers: int | None = None, progress=None) -> dict:
    """snippy + lofreq on every genome of a finished mutation-profile run, against the run's reference; the calls are kept as `mutprofile/variantcall.tsv` next to the run
    (genome, pos, ref, alt, snippy, lofreq, dp, af). The run's own SNVs (from the Squirrel / minimap2 alignment) are not touched."""
    from . import datasets, squirrel
    ref = mutprofile.summary(run).get("reference", "")
    fa, _gff = mutprofile.reference_by_name(ref)
    ds = datasets.get(run["dataset_id"])
    members = set(nextstrain.run_members(run))
    tmp = Path(run["run_dir"]) / "vc_input.fasta"
    squirrel.write_input(ds, tmp)
    try:
        res = call_many(tmp, fa, members, workers=workers or max(1, (os.cpu_count() or 4) // 2), progress=progress)
    finally:
        tmp.unlink(missing_ok=True)
    rows = []
    for g, r in res.items():
        for v in sorted(set(r["snippy"]) | set(r["lofreq"])):
            lf = r["lofreq"].get(v) or {}
            rows.append((g, *v, int(v in r["snippy"]), int(v in r["lofreq"]), lf.get("dp") or "", lf.get("af") if lf.get("af") is not None else ""))
    out = _vc_file(run)
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("genome\tpos\tref\talt\tsnippy\tlofreq\tdp\taf\n")
        f.write("".join("\t".join(map(str, r)) + "\n" for r in rows))
    ids = Path(run["run_dir"]) / "mutprofile" / "variantcall.genomes.txt"
    ids.write_text("\n".join(sorted(res)) + "\n", encoding="utf-8")
    return {"genomes": len(res), "calls": len(rows), "reference": ref}


def available(run: dict, source: str = "direct") -> bool:
    return _vc_file(run, source).exists()


def concordance_for_run(run: dict, max_n: float | None = 0.05, limit: int = 300, source: str = "direct") -> dict:
    """Agreement of three callers on the SNVs of a mutation-profile run: the run's own (`aligner`, Squirrel alignment by default), snippy and lofreq. The ITR is left out like everywhere
    else. A site is listed when the callers do not all agree; for each: how many genomes each caller called it in, and in how many of them it is called by all."""
    if not available(run, source):
        return {"available": False}
    win = mutprofile.run_itr(run)
    inwin = (lambda p: win[0] <= p <= win[1]) if win else (lambda p: True)
    judged = set(_vc_file(run, source).with_name("variantcall.genomes.txt").read_text(encoding="utf-8").split())
    ok = mutprofile.good_samples(run, max_n)
    if ok is not None:
        judged &= ok
    callers = {"mutprofile": {g: set() for g in judged}}
    for r in mutprofile._load_mutations(run):                       # the run's own calls, with the ITR already left out (mutprofile._load_mutations)
        if r["sample"] in judged:
            callers["mutprofile"][r["sample"]].add((int(r["pos"]), r["ref"].upper(), r["alt"].upper()))
    callers.update(flowlogic.read_calls(_vc_file(run, source), judged, inwin))
    rep = flowlogic.sites_report(callers, limit)
    rep.pop("all_sites")
    return {"available": True, "max_n": max_n, **rep}
