"""Checks that the installation works - not that the code is right (that is `pytest`): environments exist, each tool answers, and a tiny end-to-end run of the Snakemake variant flow
finds a planted substitution. Run through `./test.sh check`. Exit code 0 = everything required works.

Each line is OK / WARN (optional part missing) / FAIL (a required part is missing or broken)."""
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
results: list[tuple[str, str, str]] = []


def rec(level: str, what: str, detail: str = "") -> None:
    results.append((level, what, detail))
    mark = {"OK": "\033[32mOK  \033[0m", "WARN": "\033[33mWARN\033[0m", "FAIL": "\033[31mFAIL\033[0m"}[level]
    print(f"[{mark}] {what}" + (f" - {detail}" if detail else ""), flush=True)


def version(cmd: list[str], env=None) -> str:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=env)
        out = (r.stdout + r.stderr).strip().splitlines()
        return out[-1].strip() if out and r.returncode == 0 else ""
    except Exception:
        return ""


def main() -> int:
    try:
        from genomecurator import __version__, nextstrain, qc, squirrel, variantcall, variantflow
    except Exception as e:
        rec("FAIL", "import genomecurator", str(e))
        return 1
    rec("OK", "app package", f"genomecurator {__version__}, Python {sys.version.split()[0]}")
    try:
        import fastapi, jinja2, numpy, uvicorn  # noqa: F401
        rec("OK", "app libraries", f"fastapi {fastapi.__version__}, numpy {numpy.__version__}")
    except Exception as e:
        rec("FAIL", "app libraries", str(e))
    nc = qc.nextclade_bin()
    rec("OK" if nc else "WARN", "Nextclade (QC)", version([nc, "--version"]) if nc else "not installed - sequence statistics still work (mamba install -c bioconda nextclade)")

    sq = squirrel.status()
    rec("OK" if sq["ok"] else "WARN", "Squirrel environment", f"{sq.get('version') or ''} ({squirrel.env_name()})" if sq["ok"] else f"missing - {sq['how']}")
    if sq["ok"]:
        pre = squirrel.env_prefix()
        for t in ("minimap2", "mutation-profile", "iqtree3"):
            ok = (pre / "bin" / t).exists()
            rec("OK" if ok else "WARN", f"  {t}", "found" if ok else "missing in the Squirrel environment")
    ns = nextstrain.status()
    env_ok, wf_ok = ns["env"]["ok"], ns["workflow"]["ok"]
    rec("OK" if env_ok else "WARN", "Nextstrain environment", ns["env"]["prefix"] if env_ok else "missing - " + ns["env"]["how"])
    rec("OK" if wf_ok else "WARN", "Nextstrain mpox workflow", f"commit {ns['workflow']['commit']}" if wf_ok else "not downloaded - use the button on the Nextstrain tab")

    vc = variantcall.status()
    for t in ("snippy", "lofreq"):
        pre = variantcall.env_prefix(t)
        v = version([str(pre / "bin" / t), "--version" if t == "snippy" else "version"], nextstrain.tool_env(pre)) if pre else ""
        rec("OK" if pre and v else "WARN", f"{t} environment", v if v else "missing - see variantcall/*.yml (./setup.sh --only gc-snippy,gc-lofreq)")
    vf = variantflow.status()
    rec("OK" if vf["snakemake"]["ok"] else "WARN", "Snakemake environment", vf["snakemake"]["version"] if vf["snakemake"]["ok"] else "missing - ./setup.sh --only gc-snakemake")

    from genomecurator import beastprep, clock
    ck = clock.status()
    rec("OK" if ck["ok"] else "WARN", "TreeTime (molecular clock, tMRCA)", ck["version"] if ck["ok"] else "missing - " + ck["how"])
    bs = beastprep.status()
    rec("OK" if bs["ok"] else "WARN", "BEAST + BEAGLE (Bayesian tMRCA)", bs["version"] if bs["ok"] else "missing - " + bs["how"])
    if bs["ok"]:
        beast_smoke()
    if vf["ok"]:
        end_to_end()
    else:
        rec("WARN", "end-to-end variant flow", "skipped (needs snakemake + snippy + lofreq)")
    bad = [r for r in results if r[0] == "FAIL"]
    warn = [r for r in results if r[0] == "WARN"]
    print(f"\n{len(results) - len(bad) - len(warn)} ok, {len(warn)} warnings, {len(bad)} failures")
    if warn and not bad:
        print("Warnings are optional parts: curation, QC statistics and the app itself work without them; each analysis needs its own environment.")
    return 1 if bad else 0


def beast_smoke() -> None:
    """The BEAST XML the app writes, with random sequences and dates, run for a short chain: BEAST and BEAGLE must start, accept the XML and write a log."""
    import random
    from genomecurator import beastprep, nextstrain
    rng = random.Random(3)
    names = [f"s{i}" for i in range(12)]
    seqs = {n: "".join(rng.choice("ACGT") for _ in range(300)) for n in names}
    dates = {n: 2024.0 + i * 0.1 for i, n in enumerate(names)}
    s = {"n_tips": 12, "chain_length": 20000, "samples": 200, "clock": "strict", "tree_prior": "constant", "seed": 1, "threads": 0}
    tmp = Path(tempfile.mkdtemp(prefix="gc-beast-"))
    try:
        (tmp / "beast.xml").write_text(beastprep.build_xml(seqs, dates, {}, s, 1e-3), encoding="utf-8")
        pre = beastprep.env_prefix()
        r = subprocess.run([str(pre / "bin" / "beast"), "-overwrite", "-seed", "1", "beast.xml"], cwd=tmp, capture_output=True, text=True, timeout=600, env=nextstrain.tool_env(pre))
        log = tmp / "beast.log"
        if r.returncode == 0 and log.exists() and len(log.read_text().splitlines()) > 20:
            rec("OK", "BEAST runs the XML the app writes", "12 random sequences, 20,000 states")
        else:
            rec("FAIL", "BEAST runs the XML the app writes", (r.stdout + r.stderr)[-300:])
    except Exception as e:
        rec("FAIL", "BEAST runs the XML the app writes", f"{type(e).__name__}: {e}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def end_to_end() -> None:
    """A 4 kb random 'genome' with three planted substitutions, through the real Snakemake flow (snippy -> lofreq -> report)."""
    from genomecurator import nextstrain, variantcall, variantflow
    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(4000))
    planted = [(500, "A"), (2000, "C"), (3500, "G")]                      # position (1-based) and the new base; chosen to differ from the reference
    seq = list(ref)
    muts = []
    for pos, new in planted:
        old = ref[pos - 1]
        if old == new:
            new = "T" if old != "T" else "A"
        seq[pos - 1] = new
        muts.append((pos, old, new))
    tmp = Path(tempfile.mkdtemp(prefix="gc-check-"))
    try:
        (tmp / "ref.fa").write_text(f">ref\n{ref}\n")
        (tmp / "in.fa").write_text(f">sampleA\n{''.join(seq)}\n>sampleB\n{ref}\n")
        (tmp / "mut.tsv").write_text("sample\tpos\tref\talt\n" + "".join(f"sampleA\t{p}\t{o}\t{n}\n" for p, o, n in muts))
        (tmp / "judged.txt").write_text("sampleA\nsampleB\n")
        sm = variantflow.env_prefix()
        (tmp / "config.yaml").write_text(
            f"fasta: {tmp / 'in.fa'}\nreference: {tmp / 'ref.fa'}\nsnippy_env: {variantcall.env_prefix('snippy')}\nlofreq_env: {variantcall.env_prefix('lofreq')}\ncpus_per_job: 1\n"
            f"repo: {ROOT}\nmutations: {tmp / 'mut.tsv'}\njudged: {tmp / 'judged.txt'}\nitr_lo: 0\nitr_hi: 0\nseqid: ref\nmin_genomes: 1\n")
        r = subprocess.run([str(sm / "bin" / "snakemake"), "-s", str(variantflow.FLOW), "--configfile", str(tmp / "config.yaml"), "-d", str(tmp), "--cores", "2", "--use-conda",
                            "--conda-frontend", "conda", "--quiet", "rules"], capture_output=True, text=True, env=nextstrain.tool_env(sm), timeout=1800)
        if r.returncode != 0:
            rec("FAIL", "end-to-end variant flow", (r.stderr or r.stdout).strip()[-300:])
            return
        called = set()
        for line in list(open(tmp / "variantcall.tsv", encoding="utf-8"))[1:]:
            p = line.rstrip("\n").split("\t")
            if p[0] == "sampleA" and p[4] == "1" and p[5] == "1":
                called.add((int(p[1]), p[2], p[3]))
        missing = [m for m in muts if m not in called]
        if missing:
            rec("FAIL", "end-to-end variant flow", f"planted substitutions not found by snippy and lofreq: {missing}")
        elif (tmp / "concordance.json").exists() and (tmp / "review_sites.bed").exists():
            rec("OK", "end-to-end variant flow", f"snippy and lofreq both found the {len(muts)} planted substitutions; report written")
        else:
            rec("FAIL", "end-to-end variant flow", "the report files were not written")
    except Exception as e:
        rec("FAIL", "end-to-end variant flow", f"{type(e).__name__}: {e}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
