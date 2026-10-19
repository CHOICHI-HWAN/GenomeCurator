"""Run Squirrel (github.com/aineniamh/squirrel - MPXV alignment, sequence QC and phylogenetics) on a dataset, public and private genomes together.

Squirrel takes one FASTA and produces a reference-based alignment, a list of suspicious SNPs / sequences (QC), a maximum-likelihood tree (IQ-TREE) and an
HTML report. This module feeds it a dataset exactly like the Nextstrain runner does and keeps every run in the same table (`analysis_run`, tool = 'squirrel'):
members hash and counts, settings, tool versions, note, log; a run with the same members, settings and version is not repeated unless forced; each run can be
compared with the dataset as it is today. The tools live in their own conda environment (`squirrel/environment.yml`, name `gc-squirrel`); Squirrel is GPL-3.0
and runs as a separate program. The report is opened from the app (it loads its table / style libraries from CDNs in the browser; no data is sent anywhere).

Input names are the genome ids the export uses (the Pathoplexus accession for public genomes, the private id for private ones): plain characters only, as Squirrel
and IQ-TREE require. The input FASTA is rebuilt from the stored sequences and is not kept (members.txt records exactly what went in); the alignment is kept gzipped."""
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path

from . import datasets, export, log, nextstrain
from .db import audit, connect

_log = log.get("squirrel")
CLADES = ("cladei", "cladeia", "cladeib", "cladeii", "cladeiia", "cladeiib", "split")
BUILD_CLADE = {"clade-i": "cladei", "hmpxv1": "cladeiib", "hmpxv1_big": "cladeiib", "mpxv": "split"}       # the dataset's Nextstrain build hints at the clade
FLAGS = {"--no-mask", "--no-itr-mask", "--extract-cds", "--concatenate", "--bm-separate-dimers", "-bm", "--binary-partition-mask", "--interactive-tree"}
FLAGS_WITH_VALUE = {"--outgroups": r"[A-Za-z0-9_.|\-]+", "--point-style": r"circle|square", "--point-justify": r"left|right"}


def env_name() -> str:
    return os.environ.get("GC_SQUIRREL_ENV", "gc-squirrel")


def env_prefix() -> Path | None:
    return nextstrain.find_env(env_name(), "squirrel")


def tool_env() -> dict:
    return nextstrain.tool_env(env_prefix())


def make_shim(rdir: Path) -> Path | None:
    """IQ-TREE 3 refuses alignments with many N's ('cannot confidently detect the sequence data type') unless it is told `-st DNA`, and Squirrel does not tell it.
    A one-line `iqtree3` wrapper placed first on PATH for this run adds the option; nothing else is changed."""
    pre = env_prefix()
    real = (pre / "bin" / "iqtree3") if pre else None
    if not real or not real.exists():
        return None
    d = rdir / "shim"
    d.mkdir(exist_ok=True)
    f = d / "iqtree3"
    f.write_text(f'#!/bin/sh\nexec "{real}" -st DNA "$@"\n', encoding="utf-8")
    f.chmod(0o755)
    return d


def version() -> str:
    pre = env_prefix()
    if not pre:
        return ""
    try:
        r = subprocess.run(["squirrel", "-v"], capture_output=True, text=True, timeout=60, env=tool_env())
        return (r.stdout or r.stderr).strip().splitlines()[-1].strip()
    except Exception:
        return ""


def tool_versions() -> dict:
    pre = env_prefix()
    if not pre:
        return {}
    iq = ""
    try:
        r = subprocess.run(["iqtree", "--version"], capture_output=True, text=True, timeout=30, env=tool_env())
        iq = (r.stdout or r.stderr).strip().splitlines()[0].split(" built")[0]
    except Exception:
        pass
    return {"squirrel": version(), "iqtree": iq}


def status() -> dict:
    pre = env_prefix()
    return {"name": env_name(), "ok": bool(pre), "prefix": str(pre) if pre else "", "version": version() if pre else "",
            "how": f"mamba env create -f squirrel/environment.yml   (환경 이름 {env_name()})"}


# ------------------------------------------------------------------------------------------------ settings
def defaults_for(ds: dict) -> dict:
    return {"clade": BUILD_CLADE.get(ds.get("build") or "", "split"), "seq_qc": True, "run_phylo": True, "apobec3": False, "include_background": False, "exclude_flagged": True, "extra_args": ""}


def parse_settings(text) -> dict:
    """The dataset's Squirrel settings (JSON object or text of one): clade, seq_qc, run_phylo, apobec3, include_background, extra_args. Unknown keys and unsafe
    flags are refused; the output location and temp folders are managed by the tool."""
    if isinstance(text, dict):
        d = dict(text)
    elif not (text or "").strip():
        return {}
    else:
        try:
            d = json.loads(text)
        except ValueError as e:
            raise ValueError(f"설정(JSON)을 읽을 수 없습니다: {e}") from None
    if not isinstance(d, dict):
        raise ValueError("설정은 JSON 객체여야 합니다.")
    unknown = sorted(set(d) - {"clade", "seq_qc", "run_phylo", "apobec3", "include_background", "exclude_flagged", "extra_args", "threads"})
    if unknown:
        raise ValueError(f"알 수 없는 설정: {', '.join(unknown)}")
    out = {}
    if d.get("clade") not in (None, ""):
        if d["clade"] not in CLADES:
            raise ValueError(f"clade는 {', '.join(CLADES)} 중 하나입니다.")
        out["clade"] = d["clade"]
    for k in ("seq_qc", "run_phylo", "apobec3", "include_background", "exclude_flagged"):
        if k in d:
            out[k] = bool(d[k])
    if str(d.get("threads") or "").strip():
        if not str(d["threads"]).isdigit() or int(d["threads"]) < 1:
            raise ValueError("threads는 1 이상의 정수입니다.")
        out["threads"] = int(d["threads"])
    extra = str(d.get("extra_args") or "").strip()
    if extra:
        try:
            toks = shlex.split(extra)
        except ValueError as e:
            raise ValueError(f"추가 옵션을 읽을 수 없습니다: {e}") from None
        i = 0
        while i < len(toks):
            t = toks[i]
            if t in FLAGS:
                i += 1
            elif t in FLAGS_WITH_VALUE and i + 1 < len(toks) and re.fullmatch(FLAGS_WITH_VALUE[t], toks[i + 1]):
                i += 2
            else:
                raise ValueError(f"허용되지 않는 추가 옵션: {t} (허용: {', '.join(sorted(FLAGS | set(FLAGS_WITH_VALUE)))}; 입력·출력 위치는 도구가 관리합니다)")
        out["extra_args"] = extra
    return out


def effective(ds: dict) -> dict:
    return {**defaults_for(ds), **parse_settings(ds.get("squirrel_json", ""))}


def settings_sha(ds: dict) -> str:
    return hashlib.sha256(json.dumps(effective(ds), sort_keys=True).encode()).hexdigest()[:16]


def wants_phylo(s: dict) -> bool:
    return bool(s.get("apobec3") or s.get("run_phylo"))


def command(ds: dict, fasta: str, outdir: str, tmpdir: str, cores: int | None = None, qc_only: bool = False, exclude: str | None = None) -> list[str]:
    """The squirrel command line. qc_only is the first pass of a phylogeny run: alignment + QC, no tree (its suggested_to_exclude.csv feeds the second pass)."""
    s = effective(ds)
    cmd = ["squirrel", fasta, "--outdir", outdir, "--outfile", "alignment.fasta", "--tempdir", tmpdir, "--clade", s["clade"], "-t", str(cores or s.get("threads") or os.cpu_count() or 4), "--verbose"]
    if s.get("seq_qc") or qc_only:
        cmd.append("--seq-qc")
    if exclude:
        cmd += ["--exclude", exclude]
    if qc_only:
        return cmd + (shlex.split(s["extra_args"]) if s.get("extra_args") else [])
    if s.get("apobec3"):
        cmd.append("--run-apobec3-phylo")
    elif s.get("run_phylo"):
        cmd.append("--run-phylo")
    if s.get("include_background") and (s.get("apobec3") or s.get("run_phylo")):
        cmd.append("--include-background")
    if s.get("extra_args"):
        cmd += shlex.split(s["extra_args"])
    return cmd


# ------------------------------------------------------------------------------------------------ running
_STEPS = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
_RULE = re.compile(r"^(?:local)?rule (\w+):")


def write_input(ds: dict, path: Path) -> tuple[int, int, list[str]]:
    """The dataset's sequences as one FASTA (public first, then private), exactly the members an export would have. -> (public, private, member ids)"""
    _c, _r, sel, _w = export.table(ds["watch_id"], "nextstrain", ds)
    n_pub = n_priv = 0
    ids = []
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        for chunk in export.iter_fasta("nextstrain", sel):
            f.write(chunk)
            n_pub += 1
        ids += list(sel)
        if ds.get("include_private"):
            _pc, _pr, files = export.private_table(ds["watch_id"], ds)
            for chunk in export.iter_private_fasta(files):
                f.write(chunk)
                n_priv += 1
            ids += ["private:" + i for i in files]
    return n_pub, n_priv, ids


def run(ds: dict, cores: int | None = None, force: bool = False, note: str = "", progress=None) -> dict:
    """One Squirrel run on the dataset. Raises with the reason when something is missing; an identical finished run is returned instead unless force."""
    p = progress or (lambda *a: None)
    st = status()
    if not st["ok"]:
        raise RuntimeError(f"Squirrel 환경({env_name()})이 없습니다. 먼저: {st['how']}")
    settings = effective(ds)
    sha, ver = settings_sha(ds), st["version"]
    p(2, "데이터셋 멤버를 확인하는 중…")
    members = nextstrain.current_members(ds)
    if not members:
        raise ValueError("데이터셋에 멤버가 없습니다 (서열이 있는 유전체가 있어야 합니다).")
    msha = datasets.members_hash(members)
    if not force:
        for r in nextstrain.list_runs(ds["id"]):
            if r["tool"] == "squirrel" and r["status"] == "done" and r["members_sha"] == msha and r["params_sha"] == sha and r["workflow_commit"] == ver:
                p(100, f"같은 데이터·설정·버전으로 이미 분석했습니다 ({r['name']})")
                return {**r, "reused": True}
    now = datetime.now()
    virus = nextstrain.datasets_virus_name(ds)
    name = nextstrain.slug_name(virus, ds["name"], now.strftime("%Y%m%d-%H%M%S")).replace("-20", "-squirrel-20", 1)
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    rdir = nextstrain.runs_dir() / export.slug(virus) / name
    (rdir / "squirrel").mkdir(parents=True, exist_ok=True)
    tmpdir = rdir / "tmp"
    tmpdir.mkdir(exist_ok=True)
    logp = rdir / "run.log"
    versions = tool_versions()
    con = connect()
    try:
        rid = con.execute("INSERT INTO analysis_run(dataset_id,watch_id,name,status,started,runtime,build,rules_json,params_yaml,params_sha,members_sha,n_public,n_private,workflow_commit,run_dir,log_path,note,versions_json,dataset_name,tool)"
                          " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                          (ds["id"], ds["watch_id"], name, "running", now.isoformat(timespec="seconds"), "ambient", settings["clade"], json.dumps(ds["rules"], ensure_ascii=False),
                           json.dumps(settings, ensure_ascii=False), sha, msha, sum(1 for m in members if not m.startswith("private:")), sum(1 for m in members if m.startswith("private:")),
                           ver, str(rdir), str(logp), note.strip(), json.dumps(versions), ds["name"], "squirrel")).lastrowid
        con.commit()
    finally:
        con.close()

    def finish(status: str, error: str = "", result: str = ""):
        con = connect()
        try:
            con.execute("UPDATE analysis_run SET status=?, finished=?, error=?, auspice_file=? WHERE id=?", (status, datetime.now().isoformat(timespec="seconds"), error, result, rid))
            audit(con, "squirrel_run", f"{name}: {status}" + (f" ({error[:120]})" if error else ""))
            con.commit()
        finally:
            con.close()
    try:
        p(6, "서열을 하나의 FASTA로 모으는 중…")
        fasta = rdir / "input.fasta"
        write_input(ds, fasta)
        (rdir / "members.txt").write_text("\n".join(members) + "\n", encoding="utf-8")
        logf = open(logp, "w", encoding="utf-8")
        run_env = tool_env()
        shim = make_shim(rdir)
        if shim:
            run_env["PATH"] = f"{shim}{os.pathsep}{run_env['PATH']}"
            logf.write("! iqtree3 is called as `iqtree3 -st DNA ...` (IQ-TREE 3 cannot detect the data type of an alignment with many N's by itself).\n")
        ansi = re.compile(r"\x1b\[[0-9;]*m")
        deadline = time.time() + float(os.environ.get("GC_SQUIRREL_TIMEOUT_H", "12")) * 3600
        state = {"pct": 12.0, "last": "Squirrel 시작"}

        def stream(cmd: list[str], lo: float, hi: float, label: str) -> None:
            """One squirrel call: output to the log, progress between lo and hi."""
            logf.write("\n$ " + " ".join(shlex.quote(c) for c in cmd) + f"\n(in {rdir})\n\n")
            proc = subprocess.Popen(cmd, cwd=str(rdir), env=run_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True)
            state["pct"] = lo
            for line in proc.stdout:
                logf.write(line)
                logf.flush()
                s_ = ansi.sub("", line).strip()
                m, r = _STEPS.search(s_), _RULE.match(s_)
                if r:
                    state["last"] = r.group(1)
                if m:
                    state["pct"] = lo + (hi - lo) * int(m.group(3)) / 100
                    p(state["pct"], f"{label}{state['last']} ({m.group(1)}/{m.group(2)} 단계)")
                elif s_ and not s_.startswith(("[", "    ", "input", "output", "jobid", "reason", "resources", "localrule", "rule")):
                    p(state["pct"], label + s_[:80])
                if time.time() > deadline:
                    import signal
                    os.killpg(proc.pid, signal.SIGTERM)
                    raise RuntimeError("제한 시간(GC_SQUIRREL_TIMEOUT_H)을 넘겨 중단했습니다.")
            code = proc.wait()
            if code != 0:
                tail = "".join(ansi.sub("", logp.read_text(encoding="utf-8")).splitlines(True)[-25:])
                raise RuntimeError(f"Squirrel이 실패했습니다 (종료 코드 {code}). 로그 끝부분:\n{tail}")
        out = rdir / "squirrel"
        excluded: list[str] = []
        if wants_phylo(settings):
            # IQ-TREE cannot call the data type of an alignment with many N's, so: pass 1 = alignment + QC, then leave out what Squirrel's own QC suggests
            # (N content > 20% ...) and build the tree in pass 2 - the way Squirrel's QC is meant to be used.
            qc_out = rdir / "squirrel_qc"
            qc_out.mkdir(exist_ok=True)
            stream(command(ds, str(fasta), str(qc_out), str(tmpdir), cores, qc_only=True), 12, 50, "[1/2 QC] ")
            sug = qc_out / "suggested_to_exclude.csv"
            exclude_file = None
            if settings.get("exclude_flagged", True) and sug.exists():
                import csv as _csv
                rows = list(_csv.DictReader(open(sug, encoding="utf-8")))
                excluded = [r["name"] for r in rows if r.get("name")]
                if excluded:
                    exclude_file = str(sug)
                    logf.write(f"\n! QC가 제안한 제외 {len(excluded)}건을 계통수 분석에서 뺍니다 ({sug.name}).\n")
            stream(command(ds, str(fasta), str(out), str(tmpdir), cores, exclude=exclude_file), 50, 94, "[2/2 계통수] ")
            if sug.exists():
                shutil.copy2(sug, out / "qc_suggested_to_exclude.csv")
            shutil.rmtree(qc_out, ignore_errors=True)
        else:
            stream(command(ds, str(fasta), str(out), str(tmpdir), cores), 12, 94, "")
        logf.close()
        out = rdir / "squirrel"
        rep = next(iter(sorted(out.glob("*.report.html"))), None)
        if wants_phylo(settings) and not list(out.glob("*.tree")):
            raise RuntimeError("계통수를 만들지 못했습니다 (Squirrel은 정상 종료했지만 *.tree가 없습니다). 로그의 Full Traceback을 확인하세요 — 서열 N 비율이 너무 높거나 정렬이 비어 있을 수 있습니다.")
        if not rep:
            raise RuntimeError("Squirrel이 보고서를 만들지 못했습니다 (모든 서열이 정렬되지 않았을 수 있습니다). 로그를 확인하세요.")
        p(96, "결과를 정리하는 중…")
        for aln in out.glob("*.aln.fasta"):                      # the alignment is big: keep it gzipped
            with open(aln, "rb") as src, gzip.open(str(aln) + ".gz", "wb", compresslevel=6) as dst:
                shutil.copyfileobj(src, dst)
            aln.unlink()
        fasta.unlink(missing_ok=True)                            # rebuilt from the stored sequences + members.txt when needed
        shutil.rmtree(tmpdir, ignore_errors=True)
        shutil.rmtree(rdir / "shim", ignore_errors=True)
        (rdir / "versions.json").write_text(json.dumps({"tools": versions, "settings": settings, "squirrel": ver, "excluded_by_qc": excluded}, indent=2), encoding="utf-8")
        finish("done", result=f"squirrel/{rep.name}")
        p(100, "분석 완료")
    except BaseException as e:
        finish("error", str(e))
        shutil.rmtree(tmpdir, ignore_errors=True)
        shutil.rmtree(rdir / "shim", ignore_errors=True)
        raise
    return nextstrain.get_run(rid)


def result_files(run: dict) -> list[dict]:
    """What a finished run left (name, relative path, size) - for the page's file list."""
    out = Path(run["run_dir"] or "") / "squirrel"
    if not out.exists():
        return []
    return [{"name": f.name, "path": f"squirrel/{f.name}", "size": f.stat().st_size} for f in sorted(out.iterdir()) if f.is_file()]
