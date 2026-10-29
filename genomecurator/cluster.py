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
import uuid
from datetime import datetime, timedelta
from pathlib import Path

from . import config, log
from .db import audit_now, connect

_log = log.get("cluster")
LEASE = int(os.environ.get("GC_CLUSTER_LEASE", "300"))          # seconds without a heartbeat before a running job goes back to the queue
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


def _touch(con, nid: str, caps: dict | None = None) -> None:
    if caps is not None:
        con.execute("UPDATE cluster_node SET last_seen=?, caps_json=? WHERE id=?", (_now(), json.dumps(caps), nid))
    else:
        con.execute("UPDATE cluster_node SET last_seen=? WHERE id=?", (_now(), nid))


def hello(node: dict, caps: dict) -> dict:
    con = connect()
    try:
        _touch(con, node["id"], caps)
        con.commit()
    finally:
        con.close()
    return {"node": node["name"], "lease": LEASE, "allow_private": bool(node["allow_private"]), "kinds": list(KINDS)}


def _online(last_seen: str | None) -> bool:
    try:
        return bool(last_seen) and datetime.now() - datetime.fromisoformat(last_seen) < timedelta(seconds=LEASE)
    except ValueError:
        return False


def list_nodes() -> list[dict]:
    con = connect()
    try:
        rows = [dict(r) for r in con.execute("SELECT * FROM cluster_node WHERE revoked=0 ORDER BY name")]
        busy = {r["node_id"]: r["label"] for r in con.execute("SELECT node_id,label FROM cluster_job WHERE state='running'")}
    finally:
        con.close()
    for r in rows:
        r.pop("key_sha", None)
        r["caps"] = json.loads(r.pop("caps_json") or "{}")
        r["online"] = _online(r["last_seen"])
        r["running"] = busy.get(r["id"])
    return rows


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


def submit_beast(folder: str | Path, label: str = "", max_attempts: int = 3, private: bool | None = None) -> str:
    """Queue the BEAST package in `folder` (its beast.xml and settings.json). The result files are put back into the same folder."""
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
    spec = {"seed": int(st.get("seed", 1)), "threads": int(st.get("threads") or 0), "chain_length": int(st.get("chain_length") or 0), "tag": st.get("tag", "")}
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


def _job(con, jid: str) -> dict | None:
    r = con.execute("SELECT * FROM cluster_job WHERE id=?", (jid,)).fetchone()
    return dict(r) if r else None


def _public(j: dict, with_log: bool = False) -> dict:
    out = {k: j.get(k) for k in ("id", "kind", "label", "state", "node_id", "created", "claimed", "finished", "heartbeat", "attempts", "max_attempts", "pct", "msg", "error", "private")}
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


def reap() -> int:
    """Running jobs whose worker has not been heard from for LEASE seconds go back to the queue (or to error after the last attempt)."""
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


def _release(con, jid: str, why: str) -> None:
    j = _job(con, jid)
    if not j or j["state"] != "running":
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


def claim(node: dict, kinds: list[str] | None = None) -> dict | None:
    """Give the worker the oldest queued job it may run (its kind is known to the worker; a private job only to a worker that accepts private genomes)."""
    reap()
    kinds = [k for k in (kinds or list(KINDS)) if k in KINDS]
    if not kinds:
        return None
    con = connect()
    try:
        _touch(con, node["id"])
        busy = con.execute("SELECT 1 FROM cluster_job WHERE state='running' AND node_id=?", (node["id"],)).fetchone()
        if busy:
            con.commit()
            return None                     # one job at a time per worker
        q = "SELECT * FROM cluster_job WHERE state='queued' AND kind IN (%s)" % ",".join("?" * len(kinds))
        if not node["allow_private"]:
            q += " AND private=0"
        for j in con.execute(q + " ORDER BY created, rowid", kinds).fetchall():
            n = con.execute("UPDATE cluster_job SET state='running', node_id=?, claimed=?, heartbeat=?, attempts=attempts+1, pct=0, msg='assigned' WHERE id=? AND state='queued'",
                            (node["id"], _now(), _now(), j["id"])).rowcount
            if n:
                _addlog(con, j["id"], f"assigned to {node['name']}")
                con.commit()
                return {"id": j["id"], "kind": j["kind"], "spec": json.loads(j["spec_json"] or "{}"), "inputs": json.loads(j["inputs_json"] or "[]"),
                        "outputs": KINDS[j["kind"]]["outputs"], "required": KINDS[j["kind"]]["required"], "lease": LEASE}
        con.commit()
        return None
    finally:
        con.close()


def heartbeat(node: dict, jid: str, pct: float | None = None, msg: str | None = None) -> dict:
    """-> {"continue": False} when the job was cancelled or given to another worker (the worker stops and removes its copy)."""
    con = connect()
    try:
        _touch(con, node["id"])
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


def state() -> dict:
    reap()
    return {"nodes": list_nodes(), "jobs": list_jobs(100), "lease": LEASE, "now": _now()}
