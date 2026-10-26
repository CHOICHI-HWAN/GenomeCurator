"""The Snakemake version of the independent variant calls (`flows/variantcheck/Snakefile`): snippy -> lofreq for every genome of a mutation-profile run.

It runs *beside* `variantcall.py`, which stays until the flow is complete: both write the same table (`variantcall.tsv`: genome, pos, ref, alt, snippy, lofreq, dp, af), the direct one into
`<run>/mutprofile/`, this one into `<run>/variantflow/`, and `parity()` says whether the two agree. What Snakemake adds: an interrupted run continues where it stopped (finished genomes are not
redone), genomes run in parallel with retries, and the order of the steps is written in the Snakefile instead of in Python.

Snakemake, snippy and lofreq live in their own conda environments (`gc-snakemake`, `gc-snippy`, `gc-lofreq`); Snakemake activates the tool environments itself (`--use-conda` with their folders)."""
import filecmp
import os
import re
import subprocess
import time
from pathlib import Path

from . import callcache, log, mutprofile, nextstrain, squirrel, variantcall

_log = log.get("variantflow")
FLOW = Path(__file__).resolve().parent.parent / "flows" / "variantcheck" / "Snakefile"
STEPS = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")


def env_prefix() -> Path | None:
    return nextstrain.find_env(os.environ.get("GC_SNAKEMAKE_ENV", "gc-snakemake"), "snakemake")


def status() -> dict:
    sm = env_prefix()
    st = variantcall.status()
    ver = ""
    if sm:
        try:
            ver = subprocess.run([str(sm / "bin" / "snakemake"), "--version"], capture_output=True, text=True, timeout=60, env=nextstrain.tool_env(sm)).stdout.strip()
        except Exception:
            pass
    return {"ok": bool(sm) and st["ok"] and FLOW.exists(), "snakemake": {"ok": bool(sm), "version": ver, "env": os.environ.get("GC_SNAKEMAKE_ENV", "gc-snakemake")}, "tools": st,
            "how": "mamba create -n gc-snakemake -c conda-forge -c bioconda snakemake   (+ variantcall/snippy.yml, variantcall/lofreq.yml)"}


def _running(wd: Path) -> bool:
    """Is a Snakemake process working in this folder (its command line names the folder)?"""
    r = subprocess.run(["pgrep", "-f", f"snakemake.* -d {wd}"], capture_output=True, text=True)
    return bool(r.stdout.split()) and any(int(x) != os.getpid() for x in r.stdout.split())


def workdir(run: dict) -> Path:
    return Path(run["run_dir"]) / "variantflow"


def result_file(run: dict) -> Path:
    return workdir(run) / "variantcall.tsv"


def available(run: dict) -> bool:
    return result_file(run).exists()


def member_names(run: dict) -> set:
    """The run's members as the sequence names in the input FASTA: a private genome is `private:<id>` in the member list and `<id>` in the FASTA, the mutation table and the call tables."""
    return {m.split(":", 1)[1] if m.startswith("private:") else m for m in nextstrain.run_members(run)}


def run_for_mutprofile(run: dict, cores: int | None = None, progress=None, max_n: float | None = 0.05) -> dict:
    """Run the flow for a finished mutation-profile run. Calling it again after an interruption continues in the same folder."""
    from . import datasets
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"Snakemake / snippy / lofreq 환경이 없습니다: {st['how']}")
    p = progress or (lambda *a: None)
    ds = datasets.get(run["dataset_id"])
    members = member_names(run)
    ref = mutprofile.summary(run).get("reference", "")
    fa, _gff = mutprofile.reference_by_name(ref)
    wd = workdir(run)
    wd.mkdir(parents=True, exist_ok=True)
    tmp = wd / "input.new.fasta"
    squirrel.write_input(ds, tmp)
    keep = wd / "input.fasta"                                    # the members only, in a stable order; an identical file is left alone so that Snakemake does not see it as changed
    with open(wd / "input.sel.fasta", "w", encoding="utf-8") as out:
        for name, seq in variantcall._iter(tmp):
            if name in members:
                out.write(f">{name}\n{seq}\n")
    tmp.unlink()
    sel = wd / "input.sel.fasta"
    if keep.exists() and filecmp.cmp(sel, keep, shallow=False):
        sel.unlink()
    else:
        sel.replace(keep)
    # calls already made for the same sequence, reference and tool versions are taken from the shared cache (callcache.py); only the new genomes are called
    records = list(variantcall._iter(keep))
    tv = variantcall.tool_versions()
    ck = callcache.key(mutprofile.read_fasta(fa)[0], mutprofile.read_fasta(fa)[1], tv["snippy"], tv["lofreq"])
    from_cache = callcache.prefill(records, ck, wd / "calls") if tv["snippy"] and tv["lofreq"] else 0
    ok = mutprofile.good_samples(run, max_n)
    judged = wd / "judged.txt"
    new = "\n".join(sorted(members if ok is None else members & ok)) + "\n"
    if not judged.exists() or judged.read_text(encoding="utf-8") != new:                  # rewritten only when it changes, so that Snakemake does not redo the report for nothing
        judged.write_text(new, encoding="utf-8")
    win = mutprofile.run_itr(run)
    rname = mutprofile.read_fasta(fa)[0]
    cfg = wd / "config.yaml"
    cfg.write_text(f"fasta: {keep}\nreference: {Path(fa).resolve()}\nsnippy_env: {variantcall.env_prefix('snippy')}\nlofreq_env: {variantcall.env_prefix('lofreq')}\ncpus_per_job: 2\n"
                   f"repo: {Path(__file__).resolve().parent.parent}\nmutations: {Path(run['run_dir']) / 'mutprofile' / 'mutations.tsv'}\njudged: {judged}\n"
                   f"itr_lo: {win[0] if win else 0}\nitr_hi: {win[1] if win else 0}\nseqid: {rname}\nmin_genomes: 3\n", encoding="utf-8")
    sm = env_prefix()
    n = cores or max(2, (os.cpu_count() or 4))
    cmd = [str(sm / "bin" / "snakemake"), "-s", str(FLOW), "--configfile", str(cfg), "-d", str(wd), "--cores", str(n), "--use-conda", "--conda-frontend", "conda",
           "--rerun-incomplete", "--keep-going", "--retries", "1", "--printshellcmds"]
    logp = wd / "snakemake.log"
    t0 = time.time()
    if any((wd / ".snakemake" / "locks").glob("*")) if (wd / ".snakemake" / "locks").exists() else False:
        if _running(wd):
            raise RuntimeError("같은 폴더에서 Snakemake가 이미 실행 중입니다.")
        subprocess.run(cmd[:1] + cmd[1:7] + ["--unlock"], cwd=str(wd), env=nextstrain.tool_env(sm), capture_output=True, text=True)      # a run that was killed leaves its lock behind
    proc = None
    try:
        with open(logp, "a", encoding="utf-8") as lf:
            lf.write(f"\n$ {' '.join(cmd)}\n")
            proc = subprocess.Popen(cmd, cwd=str(wd), env=nextstrain.tool_env(sm), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
            for line in proc.stdout:
                lf.write(line)
                m = STEPS.search(line)
                if m:
                    p(float(m.group(3)), f"{m.group(1)}/{m.group(2)} 단계")
            rc = proc.wait()
    finally:
        if proc and proc.poll() is None:                           # the app's job was cancelled or crashed: do not leave Snakemake running on its own
            proc.terminate()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
    if rc != 0:
        raise RuntimeError(f"Snakemake가 실패했습니다 (종료 코드 {rc}); 다시 실행하면 끝난 유전체는 건너뛰고 이어서 합니다. 로그: {logp}")
    stored = callcache.store(records, ck, wd / "calls") if tv["snippy"] and tv["lofreq"] else 0
    return {"genomes": len((wd / "variantcall.genomes.txt").read_text(encoding="utf-8").split()), "seconds": round(time.time() - t0), "dir": str(wd),
            "from_cache": from_cache, "added_to_cache": stored}


def run_recorded(run: dict, cores: int | None = None, progress=None, max_n: float | None = 0.05) -> dict:
    """`run_for_mutprofile` with a row in `analysis_run` (tool `variantflow`): started / finished / failed, the Snakemake and tool versions, the folder and log. The row's note names the
    mutation-profile run it belongs to (`mutprofile:<id>`), and a failed or interrupted flow is simply run again - the folder keeps what was finished."""
    from . import datasets, runs
    ds = datasets.get(run["dataset_id"])
    members = nextstrain.run_members(run)
    wd = workdir(run)
    wd.mkdir(parents=True, exist_ok=True)
    st = status()
    name = f"{run['name']}-variantflow"
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    versions = {"snakemake": st["snakemake"]["version"], "flow": str(FLOW), "snippy": "", "lofreq": ""}
    rid = runs.start(ds, "variantflow", name, __import__("datetime").datetime.now(), ds.get("build") or "", {"max_n": max_n, "itr": list(mutprofile.run_itr(run) or [])}, "flow-v1",
                     members, run["members_sha"] or "", run.get("content_sha") or "", st["snakemake"]["version"], wd, wd / "snakemake.log", f"mutprofile:{run['id']}", versions)
    try:
        res = run_for_mutprofile(run, cores, progress, max_n)
    except Exception as e:
        runs.finish(rid, "variantflow", name, "error", str(e)[:400])
        raise
    runs.finish(rid, "variantflow", name, "done", "", str(wd / "concordance.json"))
    return {**res, "run": rid, "name": name}


def pending() -> list[dict]:
    """The finished mutation-profile runs (of every dataset) that have no variant-call result from the flow yet, oldest first."""
    return [r for r in reversed(nextstrain.list_runs(limit=100000)) if r["tool"] == "mutprofile" and r["status"] == "done" and not available(r)]


def run_all(cores: int | None = None, progress=None, max_n: float | None = 0.05) -> dict:
    """The flow for every run in `pending()`, one after the other (each run is recorded and can be resumed; one failing run does not stop the others)."""
    p = progress or (lambda *a: None)
    todo = pending()
    done, failed = [], []
    for i, r in enumerate(todo):
        try:
            run_recorded(r, cores, lambda pct, m, i=i, r=r: p((i + pct / 100) * 100 / len(todo), f"[{i + 1}/{len(todo)}] {r['dataset_name'] or r['name']} - {m}"), max_n)
            done.append(r["id"])
        except Exception as e:
            _log.warning("variant flow failed for run %s: %s", r["id"], e)
            failed.append({"run": r["id"], "error": str(e)[:200]})
    p(100, f"{len(done)} runs done" + (f", {len(failed)} failed" if failed else ""))
    return {"runs_done": done, "runs_failed": failed, "pending_before": len(todo)}


def _calls(path: Path) -> dict:
    out = {}
    if not path.exists():
        return out
    for line in list(open(path, encoding="utf-8"))[1:]:
        p = line.rstrip("\n").split("\t")
        out[(p[0], int(p[1]), p[2], p[3])] = (p[4], p[5])
    return out


def report(run: dict) -> dict:
    """The flow's own per-site report (`concordance.json`), or {available: False}."""
    f = workdir(run) / "concordance.json"
    if not f.exists():
        return {"available": False}
    import json
    return json.loads(f.read_text(encoding="utf-8"))


def parity(run: dict) -> dict:
    """Do the Snakemake flow and the direct implementation give the same calls? -> {available, direct, flow, same, only_direct, only_flow, differ}"""
    d, f = variantcall._vc_file(run), result_file(run)
    if not (d.exists() and f.exists()):
        return {"available": False, "direct": d.exists(), "flow": f.exists()}
    a, b = _calls(d), _calls(f)
    return {"available": True, "direct": len(a), "flow": len(b), "same": a == b, "only_direct": len(set(a) - set(b)), "only_flow": len(set(b) - set(a)),
            "differ": sum(1 for k in set(a) & set(b) if a[k] != b[k])}
