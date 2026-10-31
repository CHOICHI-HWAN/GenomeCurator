"""Cluster: worker servers on the same network take jobs from this (main) server.

The main server stays the only source of truth. A job is a fixed snapshot: its input files are listed with their SHA-256 and copied to the worker, which
runs one fixed command (no free-form commands, only the kinds in `KINDS`) and sends the result files back. The main server checks every hash, installs the
files in the folder the job came from and only then tells the worker it may delete its copy. A worker pulls its work (it is never contacted by the main
server) and shows it is alive by a heartbeat; a job whose worker goes silent is put back in the queue, so the same input gives the same result on another worker.

Private genomes stay on the main server unless the job is marked private and the worker was set up to accept them.
Nothing here reads genome content; it moves files and records who ran what.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
import secrets
import shutil
import threading
import uuid
from datetime import datetime, timedelta
from pathlib import Path

from . import config, log
from .db import audit_now, connect

_log = log.get("cluster")
LEASE = int(os.environ.get("GC_CLUSTER_LEASE", "600"))          # grace: seconds without a heartbeat before a running job is taken as failed there and handed out again
STALE = int(os.environ.get("GC_CLUSTER_STALE", "90"))
RESERVE_CPUS = int(os.environ.get("GC_CLUSTER_RESERVE_CPUS", "4"))   # cores left free on every worker, so the machine stays usable
DEFAULT_THREADS = 8                                              # assumed when a job or a worker does not say           # seconds without a heartbeat before a running job is shown as "not answering" (grace starts)
LOG_MAX = 300
# the only commands a worker runs; {x} are filled from the job's spec by the worker, never from free text
KINDS = {"beast": {"label": "BEAST", "inputs": ["beast.xml"], "outputs": ["beast.log", "beast.trees", "beast.ops", "beast.screen.log"], "required": ["beast.log", "beast.trees"]}}


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def root() -> Path:
    return config.DATA_DIR / "cluster"


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _safe_name(name: str) -> str:
    if not name or name != Path(name).name or name in (".", "..") or "\x00" in name:
        raise ValueError(f"not a plain file name: {name!r}")
    return name


# ---------------------------------------------------------------- nodes
def create_node(name: str, allow_private: bool = False, note: str = "") -> dict:
    """Register a worker. The key is returned once and only its hash is kept; it is set up on the worker with `genomecurator worker --key`."""
    name = (name or "").strip()
    if not name or len(name) > 60:
        raise ValueError("worker name: 1-60 characters")
    key = secrets.token_urlsafe(32)
    nid = uuid.uuid4().hex[:8]
    con = connect()
    try:
        if con.execute("SELECT 1 FROM cluster_node WHERE name=? AND revoked=0", (name,)).fetchone():
            raise ValueError(f"a worker named {name} already exists")
        con.execute("DELETE FROM cluster_node WHERE name=? AND revoked=1", (name,))
        con.execute("INSERT INTO cluster_node(id,name,key_sha,created,allow_private,note) VALUES(?,?,?,?,?,?)",
                    (nid, name, hashlib.sha256(key.encode()).hexdigest(), _now(), int(bool(allow_private)), note))
        con.commit()
    finally:
        con.close()
    audit_now("cluster_node", f"worker {name} registered (private genomes {'allowed' if allow_private else 'not allowed'})")
    return {"id": nid, "name": name, "key": key}


def revoke_node(nid: str) -> bool:
    con = connect()
    try:
        r = con.execute("SELECT name FROM cluster_node WHERE id=?", (nid,)).fetchone()
        if not r:
            return False
        con.execute("UPDATE cluster_node SET revoked=1 WHERE id=?", (nid,))
        con.commit()
    finally:
        con.close()
    _requeue_node(nid, "the worker was revoked")
    audit_now("cluster_node", f"worker {r['name']} revoked")
    return True


def authenticate(key: str | None) -> dict | None:
    if not key:
        return None
    digest = hashlib.sha256(key.encode()).hexdigest()
    con = connect()
    try:
        for r in con.execute("SELECT * FROM cluster_node WHERE revoked=0"):
            if hmac.compare_digest(r["key_sha"] or "", digest):
                return dict(r)
    finally:
        con.close()
    return None


def _touch(con, nid: str, caps: dict | None = None, metrics: dict | None = None) -> None:
    if metrics:
        con.execute("UPDATE cluster_node SET metrics_json=? WHERE id=?", (json.dumps(metrics), nid))
    if caps is not None:
        con.execute("UPDATE cluster_node SET last_seen=?, caps_json=? WHERE id=?", (_now(), json.dumps(caps), nid))
    else:
        con.execute("UPDATE cluster_node SET last_seen=? WHERE id=?", (_now(), nid))


def hello(node: dict, caps: dict) -> dict:
    con = connect()
    try:
        _touch(con, node["id"], caps, caps.get("metrics") if isinstance(caps, dict) else None)
        con.commit()
    finally:
        con.close()
    return {"node": node["name"], "lease": LEASE, "allow_private": bool(node["allow_private"]), "kinds": list(KINDS)}


def _online(last_seen: str | None) -> bool:
    try:
        return bool(last_seen) and datetime.now() - datetime.fromisoformat(last_seen) < timedelta(seconds=STALE)
    except ValueError:
        return False


def list_nodes() -> list[dict]:
    con = connect()
    try:
        rows = [dict(r) for r in con.execute("SELECT * FROM cluster_node WHERE revoked=0 ORDER BY name")]
        busy: dict[str, list] = {}
        for r in con.execute("SELECT * FROM cluster_job WHERE state='running' ORDER BY claimed"):
            busy.setdefault(r["node_id"], []).append(_public(dict(r)))
    finally:
        con.close()
    for r in rows:
        r.pop("key_sha", None)
        r["caps"] = json.loads(r.pop("caps_json") or "{}")
        r["metrics"] = json.loads(r.pop("metrics_json", None) or "{}")
        r["paused"] = bool(r["metrics"].get("paused"))
        r["online"] = _online(r["last_seen"])
        jobs_here = busy.get(r["id"]) or []
        r["running"] = jobs_here[0] if jobs_here else None     # kept for older pages
        r["running_jobs"] = jobs_here
        cpus = int(r["caps"].get("cpus") or 0)
        r["cpus_busy"] = sum(int((j.get("spec") or {}).get("threads") or DEFAULT_THREADS) for j in jobs_here) or len(jobs_here) * DEFAULT_THREADS
        r["cpus_free"] = max(0, cpus - RESERVE_CPUS - r["cpus_busy"]) if cpus else None
        env = r["caps"].get("env") or {}
        r["env_ok"] = env.get("ok")
    return rows


def local_node() -> dict | None:
    """The main server is a worker of itself by default: there is no reason for the machine that holds the data to sit idle while others compute.
    Created once (private genomes allowed, since nothing leaves the machine) and remembered in `app_meta`; if it is revoked, it is not created again."""
    import socket
    con = connect()
    try:
        r = con.execute("SELECT value FROM app_meta WHERE key='cluster_local_worker'").fetchone()
        if r:
            n = con.execute("SELECT id, revoked FROM cluster_node WHERE id=?", (r[0],)).fetchone()
            if n and not n["revoked"]:
                return None                      # already there and working: nothing to do
            if n:                                # the person revoked it: respect that
                return None
    finally:
        con.close()
    name = socket.gethostname() or "main"
    try:
        node = create_node(name, allow_private=True, note="the main server itself (local worker)")
    except ValueError:                           # a worker of that name exists already (e.g. added by hand)
        return None
    con = connect()
    try:
        con.execute("INSERT OR REPLACE INTO app_meta(key, value) VALUES('cluster_local_worker', ?)", (node["id"],))
        con.commit()
    finally:
        con.close()
    return node


# ---------------------------------------------------------------- join requests: a sub-server asks, the main server approves, the key is handed over once
REQUEST_TTL = 24 * 3600
MAX_PENDING = 20


def _expire(con) -> None:
    cut = (datetime.now() - timedelta(seconds=REQUEST_TTL)).isoformat(timespec="seconds")
    con.execute("DELETE FROM cluster_request WHERE state IN ('pending','rejected','approved') AND updated < ?", (cut,))
    con.execute("DELETE FROM cluster_request WHERE state='delivered' AND updated < ?", (cut,))


def request_join(name: str, ip: str, secret: str, info: dict | None = None) -> str:
    """A server asks to become a worker. It proves later that it is the same server with `secret` (only its hash is kept).
    A server whose environment check failed is refused here as well as on its own page."""
    name = (name or "").strip()[:60]
    if not name or not secret or len(secret) < 20:
        raise ValueError("a name and a secret are needed")
    info = info or {}
    if not (info.get("env") or {}).get("ok"):
        raise ValueError("the environment check of this server has not passed; install what is missing first")
    con = connect()
    try:
        _expire(con)
        if con.execute("SELECT COUNT(*) FROM cluster_request WHERE state='pending'").fetchone()[0] >= MAX_PENDING:
            raise ValueError("too many requests are waiting; ask the administrator of the main server")
        rid = uuid.uuid4().hex[:10]
        con.execute("INSERT INTO cluster_request(id,secret_sha,name,ip,created,updated,state,info_json) VALUES(?,?,?,?,?,?,?,?)",
                    (rid, hashlib.sha256(secret.encode()).hexdigest(), name, ip, _now(), _now(), "pending", json.dumps(info)))
        con.commit()
    finally:
        con.close()
    audit_now("cluster_request", f"{name} ({ip}) asked to join as a worker")
    return rid


def list_requests() -> list[dict]:
    con = connect()
    try:
        _expire(con)
        con.commit()
        rows = [dict(r) for r in con.execute("SELECT * FROM cluster_request WHERE state IN ('pending','approved') ORDER BY created")]
    finally:
        con.close()
    for r in rows:
        r.pop("secret_sha", None)
        r["info"] = json.loads(r.pop("info_json") or "{}")
    return rows


def decide_request(rid: str, approve: bool, allow_private: bool = False) -> bool:
    con = connect()
    try:
        r = con.execute("SELECT name, ip FROM cluster_request WHERE id=? AND state='pending'", (rid,)).fetchone()
        if not r:
            return False
        con.execute("UPDATE cluster_request SET state=?, allow_private=?, updated=? WHERE id=?", ("approved" if approve else "rejected", int(bool(allow_private)), _now(), rid))
        con.commit()
    finally:
        con.close()
    audit_now("cluster_request", f"{r['name']} ({r['ip']}) {'approved' if approve else 'rejected'}" + (" (private genomes allowed)" if approve and allow_private else ""))
    return True


def collect_request(rid: str, secret: str) -> dict:
    """The asking server checks its request. Once approved, the worker key is made now and returned this one time."""
    con = connect()
    try:
        r = con.execute("SELECT * FROM cluster_request WHERE id=?", (rid,)).fetchone()
        if not r or not hmac.compare_digest(r["secret_sha"] or "", hashlib.sha256((secret or "").encode()).hexdigest()):
            raise PermissionError("unknown request")
        state = r["state"]
        con.execute("UPDATE cluster_request SET updated=? WHERE id=? AND state='pending'", (_now(), rid))
        con.commit()
    finally:
        con.close()
    if state != "approved":
        return {"state": state}
    name = r["name"]
    con = connect()
    try:
        taken = con.execute("SELECT 1 FROM cluster_node WHERE name=? AND revoked=0", (name,)).fetchone()
    finally:
        con.close()
    if taken:
        name = f"{name}-{rid[:4]}"
    node = create_node(name, bool(r["allow_private"]), note=f"joined from {r['ip']} (request {rid})")
    con = connect()
    try:
        n = con.execute("UPDATE cluster_request SET state='delivered', node_id=?, updated=? WHERE id=? AND state='approved'", (node["id"], _now(), rid)).rowcount
        con.commit()
    finally:
        con.close()
    if not n:                                   # collected twice at the same moment: keep only the first key
        revoke_node(node["id"])
        return {"state": "delivered"}
    return {"state": "approved", "key": node["key"], "name": node["name"], "allow_private": bool(r["allow_private"])}


# ---------------------------------------------------------------- jobs
def _addlog(con, jid: str, msg: str) -> None:
    r = con.execute("SELECT log_json FROM cluster_job WHERE id=?", (jid,)).fetchone()
    lg = json.loads(r["log_json"] or "[]") if r else []
    if len(lg) < LOG_MAX:
        lg.append({"t": _now(), "msg": msg})
        con.execute("UPDATE cluster_job SET log_json=? WHERE id=?", (json.dumps(lg, ensure_ascii=False), jid))


def _private_ids() -> set[str]:
    con = connect()
    try:
        return {r[0] for r in con.execute("SELECT id FROM private_genome")}
    finally:
        con.close()


def submit_beast(folder: str | Path, label: str = "", max_attempts: int = 3, private: bool | None = None, group: dict | None = None,
                 extra_inputs: tuple[str, ...] = ()) -> str:
    """Queue the BEAST package in `folder` (its beast.xml and settings.json). The result files are put back into the same folder.
    `group` marks the job as one part of a larger piece of work (a split analysis, a BETS run) and is carried in the job's spec.
    `extra_inputs` are further files of the folder the run reads (a BETS stepping-stone job reads `prior.log`)."""
    d = Path(folder).resolve()
    try:
        d.relative_to(config.DATA_DIR.resolve())
    except ValueError:
        raise ValueError("the BEAST folder must be inside the data folder")
    xml = d / "beast.xml"
    if not xml.is_file():
        raise ValueError(f"{d} has no beast.xml (prepare the BEAST package first)")
    st = json.loads((d / "settings.json").read_text(encoding="utf-8")) if (d / "settings.json").is_file() else {}
    if private is None:
        private = bool(_private_ids() & set(st.get("tip_names", [])))
    jid = uuid.uuid4().hex[:8]
    jd = root() / "jobs" / jid / "in"
    jd.mkdir(parents=True)
    shutil.copy2(xml, jd / "beast.xml")
    inputs = [{"name": "beast.xml", "sha256": sha256(jd / "beast.xml"), "bytes": (jd / "beast.xml").stat().st_size}]
    for name in extra_inputs:
        f = d / _safe_name(name)
        if not f.is_file():
            raise ValueError(f"{d}에 {name}이 없습니다")
        shutil.copy2(f, jd / f.name)
        inputs.append({"name": f.name, "sha256": sha256(jd / f.name), "bytes": (jd / f.name).stat().st_size})
    spec = {"seed": int(st.get("seed", 1)), "threads": int(st.get("threads") or 0), "chain_length": int(st.get("chain_length") or 0), "tag": st.get("tag", "")}
    if group:
        spec["group"] = dict(group)
    con = connect()
    try:
        con.execute("INSERT INTO cluster_job(id,kind,label,state,created,max_attempts,private,spec_json,inputs_json,target_dir,log_json) VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                    (jid, "beast", label or f"BEAST {st.get('tag') or d.name}", "queued", _now(), max_attempts, int(bool(private)), json.dumps(spec),
                     json.dumps(inputs), str(d), "[]"))
        _addlog(con, jid, "queued" + (" (private genomes: only workers that accept them)" if private else ""))
        con.commit()
    finally:
        con.close()
    audit_now("cluster_job", f"queued {jid} ({d.name})")
    return jid


# ------------------------------------------------------------------------------------------------ one analysis over several machines
# A chain cannot be cut in two: each state follows from the one before. What is handed out in parallel is several shorter chains of the
# same analysis under different seeds, which hold together as many usable samples as the one long chain (beastprep.split_plan explains
# the arithmetic). When the last part comes back the main server joins them with LogCombiner, after checking that the chains agree.
def submit_beast_split(folder: str | Path, parts: int, label: str = "", max_attempts: int = 3, seed0: int = 1, mode: str = "replicate") -> dict:
    """Run the prepared BEAST package in `folder` on several machines at once and join the results into `<folder>-combined`.

    `mode="replicate"` (the default) gives every machine the analysis at its **full** length under its own seed; `mode="split"`
    divides one total between them. Measured on the 170-genome BDBV analysis, splitting a 50-million-state total four ways did not
    reproduce the undivided answer (the deep nodes and the clock rate came out wider and shifted, R-hat 1.07) while one undivided
    20-million-state chain did, so splitting is offered but is not the default. -> {"group", "jobs", "parts", "out", "mode"}"""
    from . import beastprep
    d = Path(folder).resolve()
    try:
        d.relative_to(config.DATA_DIR.resolve())
    except ValueError:
        raise ValueError("the BEAST folder must be inside the data folder")
    gid = uuid.uuid4().hex[:8]
    out = d.parent / f"{d.name}-combined"
    if mode not in ("replicate", "split"):
        raise ValueError("mode는 replicate 또는 split입니다.")
    dirs = beastprep.write_split(d, parts, seed0, mode)
    base = label or f"BEAST {json.loads((d / 'settings.json').read_text(encoding='utf-8')).get('tag') or d.name}"
    jobs = []
    for i, pd in enumerate(dirs, 1):
        jobs.append(submit_beast(pd, f"{base} [{i}/{parts}]", max_attempts,
                                 group={"id": gid, "part": i, "parts": parts, "out": str(out), "source": str(d), "mode": mode}))
    audit_now("cluster_job", f"{d.name}: {parts} chains over the cluster, mode {mode} (group {gid})")
    return {"group": gid, "jobs": jobs, "parts": [str(p) for p in dirs], "out": str(out), "mode": mode}


def _groups(rows: list[dict], kind: str = "split") -> dict[str, dict]:
    out: dict[str, dict] = {}
    for j in rows:
        g = (json.loads(j.get("spec_json") or "{}") or {}).get("group")
        if not g or g.get("kind", "split") != kind:
            continue
        e = out.setdefault(g["id"], {"id": g["id"], "kind": kind, "jobs": [], "dirs": {},
                                     **{k: g.get(k) for k in ("parts", "out", "source", "root", "tags", "stages", "name", "mode")}})
        e["jobs"].append({"id": j["id"], "part": g.get("part"), "tag": g.get("tag"), "stage": g.get("stage"),
                          "state": j["state"], "pct": j.get("pct"), "node_id": j.get("node_id"), "dir": j.get("target_dir")})
        key = g.get("part") if kind == "split" else f'{g.get("tag")}/{g.get("stage")}'      # a string, not a tuple: this goes out as JSON
        e["dirs"][key] = j.get("target_dir")
    for e in out.values():
        e["jobs"].sort(key=lambda x: (x["tag"] or "", x["part"] or 0, x["stage"] or ""))
        e["done"] = sum(1 for x in e["jobs"] if x["state"] == "done")
        e["failed"] = sum(1 for x in e["jobs"] if x["state"] in ("error", "cancelled"))
        if kind == "split":
            e["ready"] = e["done"] == (e["parts"] or 0) and bool(e["out"])
            ag = Path(e["out"] or "/nonexistent") / "agreement.json"
            e["joined"] = ag.is_file()
            if e["joined"]:
                try:
                    e["agreement"] = json.loads(ag.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    e["agreement"] = None
        else:
            e["expected"] = len(e["tags"] or []) * len(e["stages"] or [])
            e["ready"] = e["done"] == e["expected"] and e["expected"] > 0
            sm = Path(e["root"] or "/nonexistent") / "summary.json"
            e["joined"] = sm.is_file()
            if e["joined"]:
                try:
                    e["summary"] = json.loads(sm.read_text(encoding="utf-8"))
                except (OSError, ValueError):
                    e["summary"] = None
    return out


def _all_jobs() -> list[dict]:
    con = connect()
    try:
        return [dict(r) for r in con.execute("SELECT * FROM cluster_job ORDER BY created DESC, rowid DESC LIMIT 600")]
    finally:
        con.close()


def split_groups() -> list[dict]:
    return list(_groups(_all_jobs(), "split").values())           # newest first, as the rows came


def join_group(gid: str) -> dict:
    """Join the parts of a finished split analysis. Safe to call again: a group that is already joined is left alone."""
    from . import beastprep
    g = next((x for x in split_groups() if x["id"] == gid), None)
    if not g:
        raise ValueError(f"split group {gid} not found")
    if g["joined"]:
        return {"group": gid, "joined": True, "already": True, "out": g["out"]}
    if not g["ready"]:
        raise ValueError(f"{g['done']}/{g['parts']} parts are done; nothing to join yet")
    dirs = [g["dirs"][i] for i in sorted(g["dirs"]) if g["dirs"].get(i)]
    out = beastprep.combine_folders(dirs, Path(g["out"]), 0.1)
    agree = json.loads((out / "agreement.json").read_text(encoding="utf-8"))
    why = "the chains agree" if agree["ok"] else ("the chains do NOT agree: " + (", ".join(agree["disjoint"]) or f"R-hat {agree['worst'][1]}"))
    con = connect()
    try:
        for x in g["jobs"]:
            _addlog(con, x["id"], f"joined into {Path(g['out']).name} ({why})")
        con.commit()
    finally:
        con.close()
    audit_now("cluster_job", f"split group {gid} joined into {Path(g['out']).name} (agreement: {'ok' if agree['ok'] else 'check'})")
    return {"group": gid, "joined": True, "out": str(out), "agreement": agree}


def finish_split_groups() -> list[dict]:
    """Join every split analysis whose parts have all come back. Called from the background; errors are recorded on the jobs, not raised."""
    out = []
    for g in split_groups():
        if not g["ready"] or g["joined"]:
            continue
        try:
            out.append(join_group(g["id"]))
        except Exception as ex:                                  # a join that fails must not stop the server or the other groups
            _log.warning("split group %s could not be joined: %s", g["id"], ex)
            con = connect()
            try:
                for x in g["jobs"]:
                    _addlog(con, x["id"], f"joining the parts failed: {str(ex)[:200]}")
                con.commit()
            finally:
                con.close()
            out.append({"group": g["id"], "joined": False, "error": str(ex)})
    return out


# ------------------------------------------------------------------------------------------------ BETS over several machines
# A BETS run is already several independent BEAST analyses: per model a path/stepping-stone run and an ordinary chain, and after
# that chain a generalized stepping-stone run built from its posterior. Nothing is split here - each analysis goes whole to one
# machine - so there is no burn-in question and no agreement check; what the cluster buys is that the 6 (BETS) or 12 (model
# comparison) analyses run at the same time instead of fighting over the cores of one machine.
BETS_STAGES = {"psss": "path / stepping-stone", "mcmc": "ordinary chain (for the working prior)", "gss": "generalized stepping-stone"}


def _bets_settings(folder: Path, s: dict, tag: str, stage: str, chain: int, tips: list[str]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "settings.json").write_text(json.dumps({"seed": s["seed"], "threads": s["threads"], "chain_length": chain,
                                                      "tag": f"{s.get('name') or 'bets'}-{tag}-{stage}", "tip_names": tips}, ensure_ascii=False), encoding="utf-8")


def submit_bets(clock_run: dict, settings: dict | None = None, max_attempts: int = 3) -> dict:
    """Queue every analysis of a BETS run (or of a model comparison) as its own cluster job. The generalized stepping-stone analyses
    are queued later, each when its own ordinary chain comes back. -> {"group", "jobs", "root", "tags", "stages"}"""
    from . import beastprep
    info = beastprep.bets_prepare(clock_run, settings)
    s = info["settings"]
    tags = list(info["folders"])
    root_dir = Path(info["folders"][tags[0]]).parent
    try:
        root_dir.relative_to(config.DATA_DIR.resolve())
    except ValueError:
        raise ValueError("the BETS folder must be inside the data folder")
    stages = (["psss"] if s["estimator"] in ("both", "ps_ss") else []) + (["mcmc", "gss"] if s["estimator"] in ("both", "gss") else [])
    tips = list(info.get("tips") or [])                           # so that a run with private genomes goes only to workers that accept them
    gid = uuid.uuid4().hex[:8]
    what = "model comparison" if s["compare"] else "BETS"
    jobs = []
    for tag in tags:
        base = Path(info["folders"][tag])
        if "psss" in stages:
            _bets_settings(base / "psss", s, tag, "psss", s["initial_chain"], tips)
            jobs.append(_submit_bets_stage(gid, s, root_dir, tags, stages, tag, "psss", base / "psss", what, max_attempts))
        if "gss" in stages:
            (base / "gssrun").mkdir(parents=True, exist_ok=True)
            (base / "gssrun" / "template.xml").write_text(beastprep.bets_gss_template(clock_run, tag, s), encoding="utf-8")
            _bets_settings(base / "mcmc", s, tag, "mcmc", s["initial_chain"], tips)
            jobs.append(_submit_bets_stage(gid, s, root_dir, tags, stages, tag, "mcmc", base / "mcmc", what, max_attempts))
    (root_dir / "summary.json").unlink(missing_ok=True)
    audit_now("cluster_job", f"{what} {root_dir.name} queued as {len(tags) * len(stages)} analyses (group {gid})")
    return {"group": gid, "jobs": jobs, "root": str(root_dir), "tags": tags, "stages": stages}


def _submit_bets_stage(gid: str, s: dict, root_dir: Path, tags: list[str], stages: list[str], tag: str, stage: str,
                       folder: Path, what: str, max_attempts: int = 3) -> str:
    return submit_beast(folder, f"{what} [{tag} · {BETS_STAGES[stage]}]", max_attempts,
                        group={"id": gid, "kind": "bets", "tag": tag, "stage": stage, "root": str(root_dir),
                               "tags": tags, "stages": stages, "name": s.get("name") or ""},
                        extra_inputs=("prior.log",) if stage == "gss" else ())


def bets_groups() -> list[dict]:
    return list(_groups(_all_jobs(), "bets").values())


def finish_bets_groups() -> list[dict]:
    """Queue the generalized stepping-stone analysis of every ordinary chain that has come back, and write the summary of a run
    whose analyses are all in. Called from the background; a failure is recorded on the jobs, not raised."""
    from . import beastprep
    out = []
    for g in bets_groups():
        root_dir = Path(g["root"] or "/nonexistent")
        s = json.loads((root_dir / "settings.json").read_text(encoding="utf-8")) if (root_dir / "settings.json").is_file() else None
        if not s:
            continue
        have = {(x["tag"], x["stage"]) for x in g["jobs"]}
        for x in [j for j in g["jobs"] if j["stage"] == "mcmc" and j["state"] == "done"]:
            if (x["tag"], "gss") in have:
                continue
            base = Path(x["dir"]).parent                         # <root>/<tag>/mcmc -> <root>/<tag>
            gd = base / "gssrun"
            try:
                gd.mkdir(parents=True, exist_ok=True)
                gd.joinpath("beast.xml").write_text(beastprep.bets_gss_start((gd / "template.xml").read_text(encoding="utf-8"), Path(x["dir"]) / "beast.log"), encoding="utf-8")
                shutil.copyfile(Path(x["dir"]) / "beast.log", gd / "prior.log")
                mcmc_set = json.loads((Path(x["dir"]) / "settings.json").read_text(encoding="utf-8"))
                _bets_settings(gd, s, x["tag"], "gss", 0, mcmc_set.get("tip_names", []))   # 0: the path sampling has no state count to show progress against
                jid = _submit_bets_stage(g["id"], s, root_dir, g["tags"], g["stages"], x["tag"], "gss", gd,
                                         "model comparison" if s.get("compare") else "BETS")
                out.append({"group": g["id"], "queued": jid, "tag": x["tag"], "stage": "gss"})
            except Exception as ex:
                _log.warning("BETS %s/%s: the stepping-stone analysis could not be queued: %s", g["id"], x["tag"], ex)
                con = connect()
                try:
                    _addlog(con, x["id"], f"the generalized stepping-stone analysis could not be prepared: {str(ex)[:200]}")
                    con.commit()
                finally:
                    con.close()
        if not g["ready"] or g["joined"]:
            continue
        try:
            screens = {t: {} for t in g["tags"]}
            for x in g["jobs"]:
                if x["stage"] in ("psss", "gss") and x["state"] == "done" and x["dir"]:
                    screens[x["tag"]][x["stage"]] = Path(x["dir"]) / "beast.screen.log"
            res = beastprep.bets_summary(s, screens)
            (root_dir / "summary.json").write_text(json.dumps(res, ensure_ascii=False), encoding="utf-8")
            audit_now("cluster_job", f"BETS group {g['id']} summarized in {root_dir.name}")
            out.append({"group": g["id"], "summary": str(root_dir / "summary.json")})
        except Exception as ex:
            _log.warning("BETS group %s could not be summarized: %s", g["id"], ex)
            out.append({"group": g["id"], "error": str(ex)})
    return out


def _job(con, jid: str) -> dict | None:
    r = con.execute("SELECT * FROM cluster_job WHERE id=?", (jid,)).fetchone()
    return dict(r) if r else None


def _public(j: dict, with_log: bool = False) -> dict:
    out = {k: j.get(k) for k in ("id", "kind", "label", "state", "node_id", "created", "claimed", "finished", "heartbeat", "attempts", "max_attempts", "pct", "msg", "error", "private")}
    if j.get("state") == "running":
        try:
            silent = (datetime.now() - datetime.fromisoformat(j.get("heartbeat") or j.get("claimed"))).total_seconds()
        except (TypeError, ValueError):
            silent = 0
        out["silent"] = round(silent)
        out["stale"] = silent >= STALE
        out["requeue_in"] = max(0, round(LEASE - silent)) if silent >= STALE else None
    out["spec"] = json.loads(j.get("spec_json") or "{}")
    out["inputs"] = json.loads(j.get("inputs_json") or "[]")
    out["outputs"] = json.loads(j.get("outputs_json") or "[]")
    if with_log:
        out["log"] = json.loads(j.get("log_json") or "[]")
    return out


def get_job(jid: str) -> dict | None:
    con = connect()
    try:
        j = _job(con, jid)
        if not j:
            return None
        out = _public(j, True)
        n = con.execute("SELECT name FROM cluster_node WHERE id=?", (j["node_id"],)).fetchone() if j["node_id"] else None
        out["node"] = n["name"] if n else None
        return out
    finally:
        con.close()


def list_jobs(limit: int = 100) -> list[dict]:
    con = connect()
    try:
        names = {r["id"]: r["name"] for r in con.execute("SELECT id,name FROM cluster_node")}
        rows = [_public(dict(r)) for r in con.execute("SELECT * FROM cluster_job ORDER BY created DESC, rowid DESC LIMIT ?", (limit,))]
    finally:
        con.close()
    for r in rows:
        r["node"] = names.get(r["node_id"])
    return rows


_UP_SINCE = datetime.now()                                       # when this main server started: workers need time to check in again after it was away


def reap() -> int:
    """Running jobs whose worker has not been heard from for LEASE seconds go back to the queue (or to error after the last attempt).
    Silence while this server was down does not count: after a restart every worker gets the same LEASE to report in, so that a job
    that is still computing on its worker is not handed out a second time."""
    if datetime.now() - _UP_SINCE < timedelta(seconds=LEASE):
        return 0                                                 # this server has just started: every worker still has its full lease to report in
    cut = (datetime.now() - timedelta(seconds=LEASE)).isoformat(timespec="seconds")
    n = 0
    con = connect()
    try:
        for r in con.execute("SELECT id FROM cluster_job WHERE state='running' AND COALESCE(heartbeat, claimed) < ?", (cut,)).fetchall():
            _release(con, r["id"], "no heartbeat from the worker")
            n += 1
        con.commit()
    finally:
        con.close()
    return n


def _release(con, jid: str, why: str, count: bool = True) -> None:
    """Take a running job back. `count=False` when the worker gave it back on purpose (stopped by a person): that attempt is not counted."""
    j = _job(con, jid)
    if not j or j["state"] != "running":
        return
    if not count:
        con.execute("UPDATE cluster_job SET state='queued', node_id=NULL, claimed=NULL, pct=0, msg=NULL, attempts=MAX(attempts-1,0) WHERE id=?", (jid,))
        _addlog(con, jid, f"back in the queue: {why} (not counted as an attempt)")
        return
    if j["attempts"] >= j["max_attempts"]:
        con.execute("UPDATE cluster_job SET state='error', finished=?, error=?, node_id=NULL WHERE id=?", (_now(), f"{why}; no attempts left", jid))
        _addlog(con, jid, f"failed: {why}; no attempts left")
    else:
        con.execute("UPDATE cluster_job SET state='queued', node_id=NULL, claimed=NULL, pct=0, msg=NULL WHERE id=?", (jid,))
        _addlog(con, jid, f"back in the queue: {why} (attempt {j['attempts']} of {j['max_attempts']})")


def _requeue_node(nid: str, why: str) -> None:
    con = connect()
    try:
        for r in con.execute("SELECT id FROM cluster_job WHERE state='running' AND node_id=?", (nid,)).fetchall():
            _release(con, r["id"], why)
        con.commit()
    finally:
        con.close()


def claim(node: dict, kinds: list[str] | None = None, metrics: dict | None = None) -> dict | None:
    """Give the worker the oldest queued job it may run (its kind is known to the worker; a private job only to a worker that accepts private genomes)."""
    reap()
    kinds = [k for k in (list(KINDS) if kinds is None else kinds) if k in KINDS]      # an empty list (stopped, or environment broken) takes nothing
    con = connect()
    try:
        _touch(con, node["id"], metrics=metrics)
        if not kinds:
            con.commit()
            return None
        caps = json.loads(con.execute("SELECT caps_json FROM cluster_node WHERE id=?", (node["id"],)).fetchone()[0] or "{}")
        cpus = int(caps.get("cpus") or 0)
        running = [json.loads(r[0] or "{}") for r in con.execute("SELECT spec_json FROM cluster_job WHERE state='running' AND node_id=?", (node["id"],))]
        used = sum(int(x.get("threads") or DEFAULT_THREADS) for x in running)
        free = (cpus - RESERVE_CPUS - used) if cpus else (0 if running else DEFAULT_THREADS)
        if running and free < 1:
            con.commit()
            return None                     # this worker's cores are taken
        q = "SELECT * FROM cluster_job WHERE state='queued' AND kind IN (%s)" % ",".join("?" * len(kinds))
        if not node["allow_private"]:
            q += " AND private=0"
        for j in con.execute(q + " ORDER BY created, rowid", kinds).fetchall():
            want = int((json.loads(j["spec_json"] or "{}")).get("threads") or DEFAULT_THREADS)
            if running and want > free:
                continue                    # it would oversubscribe this machine; another worker can take it
            n = con.execute("UPDATE cluster_job SET state='running', node_id=?, claimed=?, heartbeat=?, attempts=attempts+1, pct=0, msg='assigned' WHERE id=? AND state='queued'",
                            (node["id"], _now(), _now(), j["id"])).rowcount
            if n:
                _addlog(con, j["id"], f"assigned to {node['name']}")
                con.commit()
                return {"id": j["id"], "kind": j["kind"], "label": j["label"], "spec": json.loads(j["spec_json"] or "{}"), "inputs": json.loads(j["inputs_json"] or "[]"),
                        "outputs": KINDS[j["kind"]]["outputs"], "required": KINDS[j["kind"]]["required"], "lease": LEASE}
        con.commit()
        return None
    finally:
        con.close()


def ask_recheck(nid: str) -> bool:
    """Ask a worker to check its environment again (it answers at its next contact; the result appears in the worker list)."""
    con = connect()
    try:
        n = con.execute("UPDATE cluster_node SET recheck=1 WHERE id=? AND revoked=0", (nid,)).rowcount
        con.commit()
    finally:
        con.close()
    return bool(n)


def pop_recheck(node: dict) -> bool:
    con = connect()
    try:
        n = con.execute("UPDATE cluster_node SET recheck=0 WHERE id=? AND recheck=1", (node["id"],)).rowcount
        con.commit()
    finally:
        con.close()
    return bool(n)


def heartbeat(node: dict, jid: str, pct: float | None = None, msg: str | None = None, metrics: dict | None = None) -> dict:
    """-> {"continue": False} when the job was cancelled or given to another worker (the worker stops and removes its copy)."""
    con = connect()
    try:
        _touch(con, node["id"], metrics=metrics)
        j = _job(con, jid)
        ok = bool(j) and j["state"] == "running" and j["node_id"] == node["id"]
        if ok:
            con.execute("UPDATE cluster_job SET heartbeat=?, pct=COALESCE(?,pct), msg=COALESCE(?,msg) WHERE id=?", (_now(), pct, msg, jid))
        con.commit()
        return {"continue": ok}
    finally:
        con.close()


def input_path(node: dict, jid: str, name: str) -> Path:
    con = connect()
    try:
        j = _job(con, jid)
    finally:
        con.close()
    if not j or j["state"] != "running" or j["node_id"] != node["id"]:
        raise PermissionError("this job is not assigned to this worker")
    name = _safe_name(name)
    if name not in {i["name"] for i in json.loads(j["inputs_json"] or "[]")}:
        raise FileNotFoundError(name)
    return root() / "jobs" / jid / "in" / name


def out_dir(jid: str) -> Path:
    return root() / "jobs" / jid / "out"


def store_output(node: dict, jid: str, name: str, tmp: Path, sha: str) -> dict:
    """Keep an uploaded result file when its SHA-256 is the one the worker declared; otherwise delete it."""
    try:
        con = connect()
        try:
            j = _job(con, jid)
        finally:
            con.close()
        if not j or j["state"] != "running" or j["node_id"] != node["id"]:
            raise PermissionError("this job is not assigned to this worker")
        name = _safe_name(name)
        if name not in KINDS[j["kind"]]["outputs"]:
            raise ValueError(f"{name} is not an output of a {j['kind']} job")
        got = sha256(tmp)
        if got != sha:
            raise ValueError(f"checksum mismatch for {name}")
        od = out_dir(jid)
        od.mkdir(parents=True, exist_ok=True)
        os.replace(tmp, od / name)
        return {"name": name, "sha256": got, "bytes": (od / name).stat().st_size}
    finally:
        Path(tmp).unlink(missing_ok=True)


def complete(node: dict, jid: str, outputs: list[dict], seconds: float | None = None) -> dict:
    """The worker says it is finished and lists its files. Every listed file must have arrived with the declared hash; the files are then installed in
    the job's folder and the answer `{"ack": true}` lets the worker delete its copy. Missing files are named so the worker sends them again."""
    con = connect()
    try:
        j = _job(con, jid)
        if not j or j["node_id"] != node["id"] or j["state"] not in ("running", "received"):
            raise PermissionError("this job is not assigned to this worker")
        spec = KINDS[j["kind"]]
        od = out_dir(jid)
        missing, listed = [], []
        for o in outputs:
            nm = _safe_name(o.get("name", ""))
            f = od / nm
            if nm not in spec["outputs"]:
                raise ValueError(f"{nm} is not an output of a {j['kind']} job")
            if not f.is_file() or sha256(f) != o.get("sha256"):
                missing.append(nm)
            else:
                listed.append({"name": nm, "sha256": o["sha256"], "bytes": f.stat().st_size})
        names = {o["name"] for o in listed}
        for r in spec["required"]:
            if r not in names and r not in missing:
                missing.append(r)
        if missing:
            con.execute("UPDATE cluster_job SET state='running', heartbeat=? WHERE id=?", (_now(), jid))
            con.commit()
            return {"ack": False, "missing": sorted(set(missing))}
        target = Path(j["target_dir"]) if j["target_dir"] else None
        installed = False
        if target and target.is_dir():
            for o in listed:
                tmp = target / (o["name"] + ".part")
                shutil.copyfile(od / o["name"], tmp)
                os.replace(tmp, target / o["name"])
            (target / "cluster_job.json").write_text(json.dumps({"job": jid, "worker": node["name"], "seconds": seconds, "attempts": j["attempts"], "outputs": listed}, indent=2), encoding="utf-8")
            installed = True
        con.execute("UPDATE cluster_job SET state='done', finished=?, pct=100, msg='done', outputs_json=?, error=NULL WHERE id=?", (_now(), json.dumps(listed), jid))
        _addlog(con, jid, f"results received from {node['name']} and checked ({len(listed)} files)" + ("; installed in the folder" if installed else "; the folder is gone, files kept in the cluster folder"))
        con.commit()
    finally:
        con.close()
    if installed:
        shutil.rmtree(root() / "jobs" / jid, ignore_errors=True)       # the folder now holds the only copy
    audit_now("cluster_job", f"{jid} done on {node['name']}")
    if json.loads(j.get("spec_json") or "{}").get("group"):
        _join_in_background()                                          # the last part of a split analysis: join without waiting for anyone to open the page
    return {"ack": True}


def fail(node: dict, jid: str, error: str) -> dict:
    con = connect()
    try:
        j = _job(con, jid)
        if not j or j["node_id"] != node["id"] or j["state"] != "running":
            return {"ack": True}                       # nothing to keep: the worker may delete its copy
        _addlog(con, jid, f"failed on {node['name']}: {error[:300]}")
        _release(con, jid, f"worker error: {error[:200]}")
        con.commit()
        return {"ack": True}
    finally:
        con.close()


def give_back(node: dict, jid: str, why: str = "stopped on the sub-server") -> dict:
    """The worker stopped the job (its Stop button): hand it out again at once."""
    con = connect()
    try:
        j = _job(con, jid)
        if j and j["node_id"] == node["id"] and j["state"] == "running":
            _release(con, jid, f"{why} ({node['name']})", count=False)
            con.commit()
    finally:
        con.close()
    return {"ack": True}


def sync_running(node: dict, running: list[str]) -> dict:
    """The worker says which jobs it runs. Any other job the main server still counts as running there is handed out again at once
    (not counted as an attempt: the job did not fail, the worker lost it, e.g. it was restarted)."""
    keep = set(running or [])
    released = []
    con = connect()
    try:
        for r in con.execute("SELECT id FROM cluster_job WHERE state='running' AND node_id=?", (node["id"],)).fetchall():
            if r["id"] not in keep:
                _release(con, r["id"], f"{node['name']} no longer runs it", count=False)
                released.append(r["id"])
        con.commit()
    finally:
        con.close()
    return {"released": released}


def cancel(jid: str) -> bool:
    con = connect()
    try:
        n = con.execute("UPDATE cluster_job SET state='cancelled', finished=?, node_id=NULL WHERE id=? AND state IN ('queued','running')", (_now(), jid)).rowcount
        if n:
            _addlog(con, jid, "cancelled")
        con.commit()
    finally:
        con.close()
    if n:
        shutil.rmtree(root() / "jobs" / jid, ignore_errors=True)
    return bool(n)


_joining = threading.Lock()


def _join_in_background() -> None:
    """What happens when a job of a group comes back - LogCombiner over gigabytes of trees, or writing the next BETS analysis -
    runs off the request thread, one at a time."""
    if _joining.locked():
        return
    def work():
        with _joining:
            finish_split_groups()
            finish_bets_groups()
            try:
                from . import beastplan
                beastplan.advance_all()                          # a plan whose step is now in starts its next one
            except Exception:
                _log.warning("could not advance the analysis plans", exc_info=True)
    threading.Thread(target=work, name="gc-group-finish", daemon=True).start()


def state() -> dict:
    reap()
    rows = _all_jobs()
    groups, bets = list(_groups(rows, "split").values()), list(_groups(rows, "bets").values())
    if any(g["ready"] and not g["joined"] for g in groups + bets) or any(
            x["stage"] == "mcmc" and x["state"] == "done" and (x["tag"], "gss") not in {(y["tag"], y["stage"]) for y in g["jobs"]}
            for g in bets for x in g["jobs"]):
        _join_in_background()
    try:
        from . import beastplan
        plans = beastplan.unfinished()
    except Exception:
        plans = []
    return {"nodes": list_nodes(), "requests": list_requests(), "jobs": list_jobs(100), "groups": groups, "bets": bets,
            "plans": plans, "joining": _joining.locked(), "lease": LEASE, "stale": STALE, "now": _now()}
