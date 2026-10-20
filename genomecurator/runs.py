"""The bookkeeping every analysis run shares (Squirrel, mutation profile, and the new flow-based pipelines): reuse of an identical finished run, a unique run name and folder, the
`analysis_run` row (started / finished) and its audit entry. Kept in one place so that a change - say, what makes two runs "the same" - is made once.

A finished run is reused when the members, the exported content, the settings and the tool version are all the same (members_sha / content_sha / params_sha / workflow_commit)."""
import json
from datetime import datetime
from pathlib import Path

from . import datasets, export, nextstrain
from .db import audit, connect


def reusable(ds: dict, tool: str, msha: str, csha: str, sha: str, ver: str) -> dict | None:
    """The newest finished run of `tool` on this dataset with the same members, content, settings and version, or None."""
    for r in nextstrain.list_runs(ds["id"]):
        if r["tool"] == tool and r["status"] == "done" and r["members_sha"] == msha and r.get("content_sha") == csha and r["params_sha"] == sha and r["workflow_commit"] == ver:
            return r
    return None


def new_run_dir(ds: dict, tool: str, subdir: str = "") -> tuple[str, Path, datetime]:
    """(unique run name, its folder - created, now). The name is `<virus>-<dataset>-<tool>-<stamp>` (with -2, -3 ... when taken)."""
    now = datetime.now()
    virus = nextstrain.datasets_virus_name(ds)
    name = nextstrain.slug_name(virus, ds["name"], now.strftime("%Y%m%d-%H%M%S")).replace("-20", f"-{tool}-20", 1)
    taken = {r["name"] for r in nextstrain.list_runs(limit=100000)}
    base, k = name, 1
    while name in taken:
        k += 1
        name = f"{base}-{k}"
    rdir = nextstrain.runs_dir() / export.slug(virus) / name
    (rdir / subdir if subdir else rdir).mkdir(parents=True, exist_ok=True)
    return name, rdir, now


def start(ds: dict, tool: str, name: str, now: datetime, build: str, settings: dict, sha: str, members: list[str], msha: str, csha: str, ver: str, rdir: Path, logp: Path,
          note: str, versions: dict) -> int:
    """Insert the `running` row of a run and return its id."""
    con = connect()
    try:
        rid = con.execute("INSERT INTO analysis_run(dataset_id,watch_id,name,status,started,runtime,build,rules_json,params_yaml,params_sha,members_sha,n_public,n_private,workflow_commit,run_dir,log_path,note,versions_json,dataset_name,tool)"
                          " VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                          (ds["id"], ds["watch_id"], name, "running", now.isoformat(timespec="seconds"), "ambient", build, json.dumps(ds["rules"], ensure_ascii=False),
                           json.dumps(settings, ensure_ascii=False), sha, msha, sum(1 for m in members if not m.startswith("private:")), sum(1 for m in members if m.startswith("private:")),
                           ver, str(rdir), str(logp), note.strip(), json.dumps(versions), ds["name"], tool)).lastrowid
        con.execute("UPDATE analysis_run SET content_sha=? WHERE id=?", (csha, rid))
        con.commit()
        return rid
    finally:
        con.close()


def finish(rid: int, tool: str, name: str, status: str, error: str = "", result: str = "") -> None:
    con = connect()
    try:
        con.execute("UPDATE analysis_run SET status=?, finished=?, error=?, auspice_file=? WHERE id=?", (status, datetime.now().isoformat(timespec="seconds"), error, result, rid))
        audit(con, f"{tool}_run", f"{name}: {status}" + (f" ({error[:120]})" if error else ""))
        con.commit()
    finally:
        con.close()
