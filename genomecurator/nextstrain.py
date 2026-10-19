"""Run the Nextstrain mpox workflow (github.com/nextstrain/mpox, `phylogenetic/`) on a dataset and keep every result.

What this module does, in the order it happens
1. **Environment**: the workflow's tools (nextstrain-cli, augur, snakemake, nextclade, iqtree, mafft, auspice) live in their own conda environment
   (`nextstrain/environment.yml`, name `gc-nextstrain`, like Superset). The app only finds it and puts its `bin/` first on PATH.
2. **Workflow checkout**: a shallow clone of nextstrain/mpox (default `data/nextstrain/mpox`, a button does the clone); its commit is recorded.
3. **Inputs + config**: the dataset is exported straight into the checkout (`export.write(nextstrain_dir=...)`): local metadata / FASTA, the private files
   when the dataset has them, and a config *overlay*. The overlay replaces the workflow's web inputs, names this run's build, and appends the dataset's
   own YAML settings (filters, subsampling, root, clock ...). It is passed after the build's own `defaults/<build>/config.yaml`, which the workflow's
   Snakefile would not load otherwise.
4. **Run**: `nextstrain build --ambient . --configfile defaults/<build>/config.yaml build-configs/genomecurator/config.yaml` with the progress and log
   streamed into the job.
5. **Keep**: every run is recorded in `analysis_run` (members hash, rules, settings, workflow commit, tool versions, log) and its Auspice files are kept under
   `data/nextstrain/auspice/<run name>.json`, its inputs / manifest / merged config under `data/nextstrain/runs/<virus>/<run name>/`. A new run never replaces an
   old one; a run with exactly the same members, settings and workflow is not repeated unless forced; each run can be compared with the dataset as it is now.
6. **View**: a local Auspice server (`auspice view`) serves the kept results, so any run can be opened by name; nothing leaves this computer."""
import atexit
import copy
import hashlib
import json
import os
import re
import shutil
import signal
import socket
import subprocess
import time
from datetime import datetime
from pathlib import Path

from . import config, datasets, export, log
from .db import audit, connect

_log = log.get("nextstrain")
REPO = "https://github.com/nextstrain/mpox"
RESERVED = {"inputs", "additional_inputs", "private_metadata", "private_sequences", "use_remote_files", "build_name", "auspice_name", "strain_id_field"}
_viewer: subprocess.Popen | None = None


# ------------------------------------------------------------------------------------------------ places and tools
def env_name() -> str:
    return os.environ.get("GC_NEXTSTRAIN_ENV", "gc-nextstrain")


def auspice_port() -> int:
    return int(os.environ.get("GC_AUSPICE_PORT", "4040"))


def checkout_dir() -> Path:
    return Path(os.environ.get("GC_NEXTSTRAIN_DIR") or (config.DATA_DIR / "nextstrain" / "mpox")).expanduser()


def workflow_dir() -> Path:
    return checkout_dir() / "phylogenetic"


def auspice_dir() -> Path:
    return config.DATA_DIR / "nextstrain" / "auspice"


def runs_dir() -> Path:
    return config.DATA_DIR / "nextstrain" / "runs"


def find_env(name: str, binary: str) -> Path | None:
    """A conda environment by name that has `binary` (found next to the app's own environment, else through conda)."""
    import sys
    cands = [Path(sys.prefix).parent / name, Path(sys.prefix).parent.parent / "envs" / name, Path(os.environ.get("CONDA_PREFIX", "/nonexistent")).parent / name]
    for c in cands:
        if (c / "bin" / binary).exists():
            return c
    for exe in ("mamba", "conda"):
        if shutil.which(exe):
            try:
                out = subprocess.run([exe, "env", "list", "--json"], capture_output=True, text=True, timeout=30).stdout
                for pth in json.loads(out).get("envs", []):
                    if Path(pth).name == name and (Path(pth) / "bin" / binary).exists():
                        return Path(pth)
            except Exception:
                pass
    return None


def env_prefix() -> Path | None:
    """The conda environment that holds the Nextstrain workflow's tools."""
    return find_env(env_name(), "nextstrain")


def tool_env(pre: Path | None = None) -> dict:
    pre = pre or env_prefix()
    env = dict(os.environ)
    if pre:
        env["PATH"] = f"{pre / 'bin'}{os.pathsep}{env.get('PATH', '')}"
        env["PYTHONNOUSERSITE"] = "1"                       # the workflow must see its own packages, not the user's site-packages
    env["PYTHONUNBUFFERED"] = "1"
    return env


def _version(cmd: list[str]) -> str:
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=30, env=tool_env())
        return (r.stdout or r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr).strip() else ""
    except Exception:
        return ""


def tool_versions() -> dict:
    if not env_prefix():
        return {}
    return {"nextstrain": _version(["nextstrain", "--version"]), "augur": _version(["augur", "--version"]), "snakemake": _version(["snakemake", "--version"]),
            "nextclade": _version(["nextclade", "--version"]), "iqtree": _version(["iqtree", "--version"]).split(" built")[0], "auspice": _version(["auspice", "--version"])}


def workflow_commit() -> str:
    d = checkout_dir()
    if not (d / ".git").exists():
        return ""
    try:
        return subprocess.run(["git", "-C", str(d), "rev-parse", "--short=10", "HEAD"], capture_output=True, text=True, timeout=30, env=tool_env()).stdout.strip()
    except Exception:
        return ""


def status() -> dict:
    pre = env_prefix()
    wf = workflow_dir()
    return {"env": {"name": env_name(), "ok": bool(pre), "prefix": str(pre) if pre else "",
                    "how": f"mamba env create -f nextstrain/environment.yml   (환경 이름 {env_name()})"},
            "workflow": {"dir": str(checkout_dir()), "ok": (wf / "Snakefile").exists(), "commit": workflow_commit(), "repo": REPO,
                         "builds": sorted(p.name for p in (wf / "defaults").iterdir() if (p / "config.yaml").exists()) if (wf / "defaults").exists() else []},
            "viewer": viewer_status(), "docker": bool(shutil.which("docker"))}


def clone_workflow(progress=None) -> dict:
    """git clone --depth 1 of nextstrain/mpox into the checkout folder (or a fast-forward pull when it is there)."""
    p = progress or (lambda *a: None)
    d = checkout_dir()
    git = shutil.which("git") or str((env_prefix() or Path("/x")) / "bin" / "git")
    if (d / ".git").exists():
        p(20, "워크플로 저장소를 업데이트하는 중…")
        r = subprocess.run([git, "-C", str(d), "pull", "--ff-only"], capture_output=True, text=True, timeout=600)
    else:
        if d.exists() and any(d.iterdir()):
            raise RuntimeError(f"{d} 이(가) 비어 있지 않습니다. 비우거나 GC_NEXTSTRAIN_DIR로 다른 위치를 정하세요.")
        d.parent.mkdir(parents=True, exist_ok=True)
        p(20, f"{REPO} 를 받는 중…")
        r = subprocess.run([git, "clone", "--depth", "1", REPO, str(d)], capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f"git 실패: {(r.stderr or r.stdout).strip()[-300:]}")
    con = connect()
    try:
        audit(con, "nextstrain_workflow", f"{REPO} @ {workflow_commit()} -> {d}")
        con.commit()
    finally:
        con.close()
    p(100, f"워크플로 {workflow_commit()}")
    return {"dir": str(d), "commit": workflow_commit()}


# ------------------------------------------------------------------------------------------------ the dataset's own settings
def _yaml():
    try:
        import yaml
        return yaml
    except ImportError:
        raise ValueError("YAML 설정을 읽으려면 pyyaml이 필요합니다 (`mamba install pyyaml`).") from None


def parse_params(text: str) -> dict:
    """The dataset's Nextstrain settings: a YAML mapping of config keys to override (nested keys merge into the workflow's config)."""
    if not (text or "").strip():
        return {}
    try:
        d = _yaml().safe_load(text)
    except ValueError:
        raise
    except Exception as e:
        raise ValueError(f"YAML을 읽을 수 없습니다: {e}") from None
    if d is None:
        return {}
    if not isinstance(d, dict):
        raise ValueError("설정은 `키: 값` 형태의 YAML이어야 합니다.")
    bad = sorted(set(d) & RESERVED)
    if bad:
        raise ValueError(f"이 키는 도구가 관리합니다(데이터셋의 입력·이름): {', '.join(bad)}")
    return d


def params_sha(text: str) -> str:
    return hashlib.sha256(json.dumps(parse_params(text), sort_keys=True, default=str).encode()).hexdigest()[:16]


def _deep(a: dict, b: dict) -> dict:
    out = copy.deepcopy(a)
    for k, v in b.items():
        out[k] = _deep(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def default_config(build: str) -> dict:
    f = workflow_dir() / "defaults" / (build or "mpxv") / "config.yaml"
    if not f.exists():
        return {}
    return _yaml().safe_load(f.read_text(encoding="utf-8")) or {}


def merged_config(ds: dict, run_name: str = "<실행 이름>") -> dict:
    """What the workflow will effectively run with: the build's defaults, then the dataset's inputs, then its own settings."""
    base = default_config(ds.get("build"))
    inputs = {"use_remote_files": False, "inputs": [{"name": "genomecurator", "metadata": "data/genomecurator/metadata.tsv", "sequences": "data/genomecurator/sequences.fasta"}],
              "build_name": run_name, "auspice_name": run_name}
    if ds.get("include_private"):
        inputs.update(private_metadata="data/genomecurator/metadata-private.tsv", private_sequences="data/genomecurator/sequences-private.fasta")
    return _deep(_deep(base, inputs), parse_params(ds.get("params_yaml", "")))


# ------------------------------------------------------------------------------------------------ state of a dataset
def current_members(ds: dict) -> list[str]:
    """The exact member ids an export of the dataset would contain now (public ids, then private ones as 'private:<id>')."""
    _c, _r, sel, _w = export.table(ds["watch_id"], "nextstrain", ds)
    ids = list(sel)
    if ds.get("include_private"):
        ids += ["private:" + i for i in export.private_table(ds["watch_id"], ds)[2]]
    return ids


def content_sha(ds: dict) -> str:
    """A fingerprint of what an analysis would read for the dataset: per member the exported metadata row and the sequence's md5 (private members: their version and md5).
    The members hash only says *which* genomes; this also changes when a genome's sequence version, curated values or QC change, so an old result is never presented as current
    for data that has changed under the same ids."""
    _c, rows, sel, _w = export.table(ds["watch_id"], "nextstrain", ds)
    con = connect()
    try:
        md5 = {r["accession"]: r["md5"] for r in con.execute("SELECT accession, md5 FROM genome_sequence")}
        priv = {r["id"]: (r["version"], r["seq_md5"]) for r in con.execute("SELECT id, version, seq_md5 FROM private_genome WHERE watch_id=?", (ds["watch_id"],))}
    finally:
        con.close()
    lines = []
    for row in rows:
        g = sel.get(row.get("PPX_accession")) or {}
        lines.append(json.dumps([row, md5.get(g.get("accession"))], sort_keys=True, ensure_ascii=False, default=str))
    if ds.get("include_private"):
        _pc, prow, _f = export.private_table(ds["watch_id"], ds)
        for row in prow:
            lines.append("private:" + json.dumps([row, priv.get(row.get("PPX_accession"))], sort_keys=True, ensure_ascii=False, default=str))
    return hashlib.sha256("\n".join(sorted(lines)).encode()).hexdigest()


def qc_missing(ds: dict) -> int:
    """Members the workflow's default filter will drop because they have no Nextclade QC (`QC_rare_mutations`): run the dataset's QC first.
    Private members are not filtered by the workflow, so they are not counted."""
    _c, rows, _s, _w = export.table(ds["watch_id"], "nextstrain", ds)
    return sum(1 for r in rows if not r.get("QC_rare_mutations"))


def slug_name(virus: str, dataset: str, stamp: str) -> str:
    """The run's name: lower case, hyphens only (Auspice reads '_' in a file name as a folder separator)."""
    clean = lambda s: re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-") or "x"
    return f"{clean(virus)}-{clean(dataset)}-{stamp}"


def _row(r) -> dict:
    d = dict(r)
    d["tool"] = d.get("tool") or "nextstrain"
    d["versions"] = json.loads(d.pop("versions_json") or "{}")
    d["rules"] = json.loads(d.pop("rules_json") or "{}")
    return d


def list_runs(dataset_id: int | None = None, watch_id: int | None = None, limit: int = 100) -> list[dict]:
    con = connect()
    try:
        q, a = "SELECT * FROM analysis_run WHERE 1=1", []
        if dataset_id:
            q += " AND dataset_id=?"
            a.append(dataset_id)
        if watch_id:
            q += " AND watch_id=?"
            a.append(watch_id)
        return [_row(r) for r in con.execute(q + " ORDER BY id DESC LIMIT ?", a + [limit])]
    finally:
        con.close()


def get_run(run_id: int) -> dict | None:
    con = connect()
    try:
        r = con.execute("SELECT * FROM analysis_run WHERE id=?", (run_id,)).fetchone()
        return _row(r) if r else None
    finally:
        con.close()


def run_members(run: dict) -> list[str]:
    f = Path(run["run_dir"] or "") / "members.txt"
    return f.read_text(encoding="utf-8").split() if f.exists() else []


def compare(run: dict, now_ids: list[str], now_params_sha: str, commit: str, now_content: str | None = None) -> dict:
    """How a past run relates to the dataset as it is today. `same_content` is None for runs made before the content fingerprint existed (only the ids can be compared)."""
    then, now = set(run_members(run)), set(now_ids)
    same_members = run["members_sha"] == datasets.members_hash(now_ids)
    same_content = (run.get("content_sha") == now_content) if (run.get("content_sha") and now_content) else None
    return {"same_members": same_members, "same_content": same_content, "added": len(now - then), "removed": len(then - now), "same_params": run["params_sha"] == now_params_sha,
            "same_workflow": (run["workflow_commit"] or "") == commit, "current": same_members and same_content is not False and run["params_sha"] == now_params_sha}


def mark_interrupted() -> int:
    con = connect()
    try:
        n = con.execute("UPDATE analysis_run SET status='interrupted', finished=?, error=? WHERE status='running'",
                        (datetime.now().isoformat(timespec="seconds"), "서버가 다시 시작되어 분석이 중단되었습니다. 같은 데이터셋으로 다시 실행하세요.")).rowcount
        con.commit()
        return n
    finally:
        con.close()


def delete_run(run_id: int) -> bool:
    """Forget a run and delete its kept files (an explicit decision: runs are otherwise kept)."""
    r = get_run(run_id)
    if not r:
        return False
    for f in (auspice_dir() / f"{r['name']}.json", auspice_dir() / f"{r['name']}_root-sequence.json"):
        f.unlink(missing_ok=True)
    if r["run_dir"] and Path(r["run_dir"]).exists():
        shutil.rmtree(r["run_dir"], ignore_errors=True)
    con = connect()
    try:
        con.execute("DELETE FROM analysis_run WHERE id=?", (run_id,))
        audit(con, "nextstrain_delete", f"{r['name']}")
        con.commit()
    finally:
        con.close()
    return True


def set_note(run_id: int, note: str) -> None:
    con = connect()
    try:
        con.execute("UPDATE analysis_run SET note=? WHERE id=?", (note.strip(), run_id))
        con.commit()
    finally:
        con.close()


# ------------------------------------------------------------------------------------------------ running
_STEPS = re.compile(r"(\d+) of (\d+) steps \((\d+)%\) done")
_RULE = re.compile(r"^(?:local)?rule (\w+):")


def run(ds: dict, runtime: str = "ambient", cores: int | None = None, force: bool = False, note: str = "", progress=None) -> dict:
    """One analysis of the dataset. Raises with the reason when something is missing; a run with the same members / settings / workflow as a finished
    one is not repeated unless force. Returns the run record."""
    p = progress or (lambda *a: None)
    st = status()
    if not st["env"]["ok"]:
        raise RuntimeError(f"Nextstrain 환경({env_name()})이 없습니다. 먼저: {st['env']['how']}")
    if not st["workflow"]["ok"]:
        raise RuntimeError(f"워크플로 저장소가 없습니다({st['workflow']['dir']}). Nextstrain 탭의 '워크플로 받기'를 누르세요.")
    if runtime not in ("ambient", "docker"):
        raise ValueError("runtime은 ambient 또는 docker입니다.")
    build = ds.get("build") or "mpxv"
    if build not in st["workflow"]["builds"]:
        raise ValueError(f"워크플로에 빌드 `{build}`가 없습니다 (있는 것: {', '.join(st['workflow']['builds'])}).")
    params = parse_params(ds.get("params_yaml", ""))
    psha, commit = params_sha(ds.get("params_yaml", "")), st["workflow"]["commit"]
    p(2, "데이터셋 멤버를 확인하는 중…")
    members = current_members(ds)
    if not members:
        raise ValueError("데이터셋에 멤버가 없습니다 (서열이 있는 유전체가 있어야 합니다).")
    msha = datasets.members_hash(members)
    csha = content_sha(ds)
    if not force:
        for r in list_runs(ds["id"]):
            if r["tool"] == "nextstrain" and r["status"] == "done" and r["members_sha"] == msha and r.get("content_sha") == csha and r["params_sha"] == psha and r["workflow_commit"] == commit and r["build"] == build:
                p(100, f"같은 데이터·설정·워크플로로 이미 분석했습니다 ({r['name']})")
                return {**r, "reused": True}
    now = datetime.now()
    name = slug_name(datasets_virus_name(ds), ds["name"], now.strftime("%Y%m%d-%H%M%S"))
    taken = {r["name"] for r in list_runs(limit=100000)}
    base, k = name, 1
    while name in taken or (auspice_dir() / f"{name}.json").exists():       # two runs in the same second (a double click) still get their own name
        k += 1
        name = f"{base}-{k}"
    rdir = runs_dir() / export.slug(datasets_virus_name(ds)) / name
    rdir.mkdir(parents=True, exist_ok=True)
    logp = rdir / "run.log"
    versions = tool_versions()
    con = connect()
    try:
        rid = con.execute("INSERT INTO analysis_run(dataset_id,watch_id,name,status,started,runtime,build,rules_json,params_yaml,params_sha,members_sha,n_public,n_private,workflow_commit,run_dir,log_path,note,versions_json,dataset_name)"
                          " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                          (ds["id"], ds["watch_id"], name, "running", now.isoformat(timespec="seconds"), runtime, build, json.dumps(ds["rules"], ensure_ascii=False), ds.get("params_yaml", ""), psha,
                           msha, sum(1 for m in members if not m.startswith("private:")), sum(1 for m in members if m.startswith("private:")), commit, str(rdir), str(logp), note.strip(),
                           json.dumps(versions), ds["name"])).lastrowid
        con.execute("UPDATE analysis_run SET content_sha=? WHERE id=?", (csha, rid))
        con.commit()
    finally:
        con.close()

    def finish(status: str, error: str = "", auspice: str = ""):
        con = connect()
        try:
            con.execute("UPDATE analysis_run SET status=?, finished=?, error=?, auspice_file=? WHERE id=?", (status, datetime.now().isoformat(timespec="seconds"), error, auspice, rid))
            audit(con, "nextstrain_run", f"{name}: {status}" + (f" ({error[:120]})" if error else ""))
            con.commit()
        finally:
            con.close()
    try:
        p(5, "입력과 설정을 워크플로 폴더에 쓰는 중…")
        res = export.write(ds["watch_id"], None, "nextstrain", lambda pct, m: p(5 + pct * 0.1, m), nextstrain_dir=workflow_dir(), dataset=ds, run_name=name,
                           extra_config=ds.get("params_yaml", ""))
        (rdir / "members.txt").write_text("\n".join(members) + "\n", encoding="utf-8")
        cfg_overlay = Path(res["config"])
        wf = workflow_dir()
        cmd = ["nextstrain", "build", f"--{runtime}", ".", "--configfile", f"defaults/{build}/config.yaml", str(cfg_overlay.relative_to(wf)),
               "--cores", str(cores or os.cpu_count() or 4), "--forceall"]
        logf = open(logp, "w", encoding="utf-8")
        logf.write("$ " + " ".join(cmd) + f"\n(in {wf})\n")
        miss = qc_missing(ds)
        if miss:
            logf.write(f"! Nextclade QC가 없는 공개 멤버 {miss:,}건은 워크플로의 기본 필터(QC_rare_mutations)가 걸러냅니다. 데이터셋 QC를 먼저 실행하면 포함됩니다.\n")
            p(14, f"주의: QC가 없는 {miss:,}건은 워크플로가 걸러냅니다")
        logf.write("\n")
        proc = subprocess.Popen(cmd, cwd=str(wf), env=tool_env(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1, start_new_session=True)
        pct, last, deadline = 15.0, "워크플로 시작", time.time() + float(os.environ.get("GC_NEXTSTRAIN_TIMEOUT_H", "12")) * 3600
        for line in proc.stdout:
            logf.write(line)
            logf.flush()
            s = line.strip()
            m = _STEPS.search(s)
            r = _RULE.match(s)
            if r:
                last = r.group(1)
            if m:
                pct = 15 + 80 * int(m.group(3)) / 100
                p(pct, f"{last} ({m.group(1)}/{m.group(2)} 단계)")
            elif r:
                p(pct, f"{last} 진행 중")
            if time.time() > deadline:
                os.killpg(proc.pid, signal.SIGTERM)
                raise RuntimeError("제한 시간(GC_NEXTSTRAIN_TIMEOUT_H)을 넘겨 중단했습니다.")
        code = proc.wait()
        logf.close()
        if code != 0:
            tail = "".join(logp.read_text(encoding="utf-8").splitlines(True)[-25:])
            raise RuntimeError(f"워크플로가 실패했습니다 (종료 코드 {code}). 로그 끝부분:\n{tail}")
        src = wf / "auspice" / f"{name}.json"
        if not src.exists():
            raise RuntimeError(f"결과 파일이 없습니다: {src}")
        p(96, "결과를 보관하는 중…")
        auspice_dir().mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, auspice_dir() / f"{name}.json")
        root = wf / "auspice" / f"{name}_root-sequence.json"
        if root.exists():
            shutil.copy2(root, auspice_dir() / f"{name}_root-sequence.json")
        gd = wf / "data" / "genomecurator"
        for f in ("metadata.tsv", "metadata-private.tsv", "manifest.json", "excluded.tsv"):
            if (gd / f).exists():
                shutil.copy2(gd / f, rdir / f)           # the inputs' tables and manifest are kept; the FASTA is rebuilt from the stored sequences + members.txt
        shutil.copy2(cfg_overlay, rdir / "config-overlay.yaml")
        rc = wf / "results" / "run_config.yaml"
        if rc.exists():
            shutil.copy2(rc, rdir / "run_config.yaml")
        (rdir / "versions.json").write_text(json.dumps({"tools": versions, "workflow_commit": commit, "runtime": runtime, "build": build}, indent=2), encoding="utf-8")
        shutil.rmtree(wf / "results" / name, ignore_errors=True)             # intermediate files of this run; the kept results are what matter
        for f in (wf / "auspice").glob(f"{name}*"):
            f.unlink(missing_ok=True)
        finish("done", auspice=f"{name}.json")
        p(100, "분석 완료")
    except BaseException as e:
        finish("error", str(e))
        raise
    return get_run(rid)


def datasets_virus_name(ds: dict) -> str:
    con = connect()
    try:
        r = con.execute("SELECT name FROM watchlist WHERE id=?", (ds["watch_id"],)).fetchone()
        return r["name"] if r else "virus"
    finally:
        con.close()


# ------------------------------------------------------------------------------------------------ the viewer (a local Auspice server)
def viewer_listening() -> bool:
    try:
        with socket.create_connection(("127.0.0.1", auspice_port()), timeout=0.5):
            return True
    except OSError:
        return False


def viewer_url(name: str = "") -> str:
    return f"http://127.0.0.1:{auspice_port()}/{name}"


def viewer_status() -> dict:
    return {"running": viewer_listening(), "port": auspice_port(), "started_by_app": _viewer is not None and _viewer.poll() is None, "url": viewer_url(), "dir": str(auspice_dir())}


def start_viewer(progress=None, wait: int = 60) -> str:
    """`auspice view` over the folder that holds every kept result (started by the app, stopped with it)."""
    global _viewer
    p = progress or (lambda *a: None)
    if viewer_listening():
        return viewer_url()
    pre = env_prefix()
    if not pre or not (pre / "bin" / "auspice").exists():
        raise RuntimeError(f"Auspice가 없습니다. Nextstrain 환경({env_name()})을 만드세요: mamba env create -f nextstrain/environment.yml")
    auspice_dir().mkdir(parents=True, exist_ok=True)
    lf = config.DATA_DIR / "logs"
    lf.mkdir(parents=True, exist_ok=True)
    out = open(lf / "auspice.log", "ab")
    p(10, "Auspice 서버를 시작하는 중…")
    env = {**tool_env(), "PORT": str(auspice_port()), "HOST": "127.0.0.1"}        # auspice 3 reads the port from the environment; local only
    _viewer = subprocess.Popen(["auspice", "view", str(auspice_dir())], stdout=out, stderr=out, env=env, start_new_session=True)
    t0 = time.time()
    while time.time() - t0 < wait:
        if _viewer.poll() is not None:
            raise RuntimeError(f"Auspice가 바로 종료되었습니다 (코드 {_viewer.returncode}). 로그: {lf / 'auspice.log'}")
        if viewer_listening():
            p(100, "Auspice가 응답합니다")
            return viewer_url()
        time.sleep(1)
    raise RuntimeError(f"{wait}초 안에 Auspice가 응답하지 않았습니다. 로그: {lf / 'auspice.log'}")


def stop_viewer() -> bool:
    global _viewer
    if _viewer is None or _viewer.poll() is not None:
        return False
    try:
        os.killpg(_viewer.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        _viewer.terminate()
    _viewer = None
    return True


atexit.register(stop_viewer)
