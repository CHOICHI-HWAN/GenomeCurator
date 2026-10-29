"""The worker of the cluster (see cluster.py): `genomecurator worker --main http://main:8000 --key KEY`.

It pulls one job at a time from the main server, copies the input files (checking their SHA-256), runs the fixed command of the job kind, sends the result
files back and deletes its own copy only after the main server has said it holds them. It keeps no data of its own. While the main server is down a running
job goes on; the finished files wait in the work folder (with their state on disk) and are sent as soon as the main server answers again, also after a
restart of the worker. Only the standard library is used, so a worker needs nothing but this package and the tools of the job kinds it offers.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from . import log

_log = log.get("worker")
KINDS = {"beast": {"outputs": ["beast.log", "beast.trees", "beast.ops", "beast.screen.log"]}}   # what each kind produces (the main server checks the same list)


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Gone(Exception):
    """The main server no longer wants this job here (cancelled, given to another worker, or unknown)."""


class Main:
    def __init__(self, base: str, key: str, timeout: int = 60):
        self.base, self.key, self.timeout = base.rstrip("/"), key, timeout

    def _open(self, method: str, path: str, data=None, headers: dict | None = None, length: int | None = None):
        h = {"Authorization": f"Bearer {self.key}", **(headers or {})}
        if length is not None:
            h["Content-Length"] = str(length)
        req = urllib.request.Request(self.base + path, data=data, method=method, headers=h)
        try:
            return urllib.request.urlopen(req, timeout=self.timeout)
        except urllib.error.HTTPError as e:
            if e.code in (403, 404, 409):
                raise Gone(f"{e.code} {e.read()[:200]!r}")
            raise

    def json(self, method: str, path: str, body: dict | None = None) -> dict:
        data = json.dumps(body or {}).encode()
        with self._open(method, path, data, {"Content-Type": "application/json"}, len(data)) as r:
            return json.loads(r.read() or b"{}")

    def download(self, path: str, dest: Path) -> None:
        with self._open("GET", path) as r, open(dest, "wb") as f:
            shutil.copyfileobj(r, f, 1 << 20)

    def upload(self, path: str, src: Path) -> dict:
        with open(src, "rb") as f:
            with self._open("PUT", path, f, {"Content-Type": "application/octet-stream"}, src.stat().st_size) as r:
                return json.loads(r.read() or b"{}")


def capabilities(beast: str | None) -> dict:
    caps = {"kinds": [], "cpus": os.cpu_count(), "platform": platform.platform(), "tools": {}}
    if beast:
        caps["kinds"].append("beast")
        try:
            out = subprocess.run([beast, "-version"], capture_output=True, text=True, timeout=60)
            import re
            m = re.search(r"v([0-9][0-9.]*)", out.stdout + out.stderr)
            caps["tools"]["beast"] = m[1] if m else "unknown"
        except Exception:
            caps["tools"]["beast"] = "unknown"
    return caps


def find_beast(explicit: str | None = None) -> str | None:
    if explicit:
        return explicit if Path(explicit).exists() else None
    try:
        from . import beastprep
        pre = beastprep.env_prefix()
        if pre and (pre / "bin" / "beast").exists():
            return str(pre / "bin" / "beast")
    except Exception:
        pass
    return shutil.which("beast")


class Worker:
    def __init__(self, base: str, key: str, workdir: str | Path, beast: str | None = None, poll: float = 10.0, beat: float = 20.0, retry: float = 15.0):
        self.main, self.work = Main(base, key), Path(workdir)
        self.beast, self.poll, self.beat, self.retry = find_beast(beast), poll, beat, retry
        self.stop = threading.Event()
        self.status = {"state": "starting", "job": None, "error": None, "since": time.strftime("%Y-%m-%dT%H:%M:%S")}
        self.work.mkdir(parents=True, exist_ok=True)

    def _set(self, state: str, job: str | None = None, error: str | None = None) -> None:
        self.status = {"state": state, "job": job, "error": error, "since": time.strftime("%Y-%m-%dT%H:%M:%S")}

    # ---- one job
    def _state(self, jd: Path, **kw) -> dict:
        p = jd / "state.json"
        st = json.loads(p.read_text()) if p.exists() else {}
        st.update(kw)
        p.write_text(json.dumps(st))
        return st

    def _progress(self, jd: Path, chain: int) -> tuple[float | None, str | None]:
        try:
            last = (jd / "run" / "beast.log").read_text(encoding="utf-8").strip().splitlines()[-1].split("\t")[0]
            if chain:
                return round(98 * min(float(last) / chain, 1.0), 1), f"{int(float(last)):,}/{chain:,} states"
        except (OSError, ValueError, IndexError):
            pass
        return None, None

    def run_job(self, job: dict) -> None:
        jid = job["id"]
        jd = self.work / "jobs" / jid
        shutil.rmtree(jd, ignore_errors=True)
        (jd / "in").mkdir(parents=True)
        (jd / "run").mkdir()
        self._state(jd, phase="running", job=job)
        self._set("running", jid)
        try:
            for i in job["inputs"]:
                dest = jd / "in" / i["name"]
                self.main.download(f"/api/cluster/w/jobs/{jid}/in/{urllib.parse.quote(i['name'])}", dest)
                if sha256(dest) != i["sha256"]:
                    raise RuntimeError(f"checksum mismatch for the input {i['name']}")
                shutil.copy2(dest, jd / "run" / i["name"])
        except Gone:
            return self._discard(jd)
        except Exception as e:
            return self._fail(jd, jid, f"could not fetch the inputs: {e}")
        if job["kind"] != "beast" or not self.beast:
            return self._fail(jd, jid, f"this worker cannot run {job['kind']} jobs")
        spec = job["spec"]
        cmd = [self.beast, "-overwrite", "-seed", str(int(spec.get("seed", 1)))] + (["-threads", str(int(spec["threads"]))] if spec.get("threads") else []) + ["beast.xml"]
        env = os.environ.copy()
        try:
            from . import beastprep, nextstrain
            pre = beastprep.env_prefix()
            if pre:
                env = nextstrain.tool_env(pre)
        except Exception:
            pass
        t0 = time.time()
        gone = threading.Event()
        with open(jd / "run" / "beast.screen.log", "w", encoding="utf-8") as lf:
            proc = subprocess.Popen(cmd, cwd=str(jd / "run"), env=env, stdout=lf, stderr=subprocess.STDOUT, text=True)
            last_beat = 0.0
            while proc.poll() is None and not self.stop.is_set():
                time.sleep(1)
                if time.time() - last_beat < self.beat:
                    continue
                last_beat = time.time()
                pct, msg = self._progress(jd, int(spec.get("chain_length") or 0))
                try:
                    if not self.main.json("POST", f"/api/cluster/w/jobs/{jid}/heartbeat", {"pct": pct, "msg": msg}).get("continue", False):
                        gone.set()
                except Gone:
                    gone.set()
                except Exception:
                    pass                            # the main server is away: the job goes on
                if gone.is_set():
                    break
            if proc.poll() is None:
                proc.terminate()
                try:
                    proc.wait(30)
                except subprocess.TimeoutExpired:
                    proc.kill()
        if gone.is_set():
            return self._discard(jd)
        if self.stop.is_set():
            return                                  # stopped on request: the job is left to the lease of the main server
        if proc.returncode != 0:
            tail = "".join(open(jd / "run" / "beast.screen.log", encoding="utf-8").readlines()[-6:])[-300:]
            return self._fail(jd, jid, f"BEAST exited with code {proc.returncode}: {tail.strip()}")
        outs = [{"name": n, "sha256": sha256(jd / "run" / n), "bytes": (jd / "run" / n).stat().st_size} for n in KINDS["beast"]["outputs"] if (jd / "run" / n).is_file()]
        self._state(jd, phase="finished", outputs=outs, seconds=round(time.time() - t0))
        self.deliver(jd)

    def _fail(self, jd: Path, jid: str, why: str) -> None:
        _log.warning("job %s failed: %s", jid, why)
        try:
            self.main.json("POST", f"/api/cluster/w/jobs/{jid}/fail", {"error": why})
        except Gone:
            pass
        except Exception:
            self._state(jd, phase="failed", error=why)      # tell the main server later
            return
        shutil.rmtree(jd, ignore_errors=True)

    def _discard(self, jd: Path) -> None:
        shutil.rmtree(jd, ignore_errors=True)

    def deliver(self, jd: Path) -> bool:
        """Send the finished files and ask the main server to take over; the copy is deleted only after `ack`. -> True when this job is settled."""
        st = json.loads((jd / "state.json").read_text())
        jid = st["job"]["id"]
        try:
            if st["phase"] == "failed":
                self.main.json("POST", f"/api/cluster/w/jobs/{jid}/fail", {"error": st.get("error", "")})
                self._discard(jd)
                return True
            for _ in range(3):
                for o in st["outputs"]:
                    self.main.upload(f"/api/cluster/w/jobs/{jid}/out/{urllib.parse.quote(o['name'])}?sha={o['sha256']}", jd / "run" / o["name"])
                r = self.main.json("POST", f"/api/cluster/w/jobs/{jid}/done", {"outputs": st["outputs"], "seconds": st.get("seconds")})
                if r.get("ack"):
                    self._discard(jd)
                    return True
            return False
        except Gone:
            self._discard(jd)
            return True
        except Exception as e:
            _log.warning("job %s: results not delivered yet (%s); they stay here and are sent again", jid, e)
            return False

    def resume(self) -> None:
        for jd in sorted((self.work / "jobs").glob("*")):
            try:
                st = json.loads((jd / "state.json").read_text())
            except (OSError, ValueError):
                shutil.rmtree(jd, ignore_errors=True)
                continue
            if st.get("phase") in ("finished", "failed"):
                self.deliver(jd)
            else:
                shutil.rmtree(jd, ignore_errors=True)        # interrupted while running: the main server hands the job out again

    # ---- loop
    def serve(self) -> None:
        caps = capabilities(self.beast)
        _log.info("worker ready: %s", caps)
        known = False
        while not self.stop.is_set():
            try:
                if not known:
                    self.main.json("POST", "/api/cluster/w/hello", caps)
                    known = True
                self.resume()
                job = self.main.json("POST", "/api/cluster/w/claim", {"kinds": caps["kinds"]}).get("job")
                if job:
                    self.run_job(job)
                    self._set("idle")
                    continue
                if self.status["state"] != "idle":
                    self._set("idle")
            except Gone as e:
                _log.error("the main server refused this worker (%s); check the key", e)
                self._set("refused", error="the main server did not accept the key or the address")
                self.stop.wait(60)
                continue
            except Exception as e:
                _log.warning("main server not reachable (%s)", e)
                self._set("unreachable", error=str(e)[:200])
                known = False
                self.stop.wait(self.retry)
                continue
            self.stop.wait(self.poll)
        _log.info("worker stopped")


# ---- the "activate this server as a sub-server" switch of the Cluster page: the worker runs inside this program on a background thread
def _cfg_path() -> Path:
    from . import config
    return config.DATA_DIR / "cluster" / "worker.json"


def load_config() -> dict:
    try:
        return json.loads(_cfg_path().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def save_config(cfg: dict) -> None:
    p = _cfg_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(cfg), encoding="utf-8")
    try:
        os.chmod(p, 0o600)                      # holds the worker key; kept out of the database and so out of backups
    except OSError:
        pass


_sub: dict = {"worker": None, "thread": None}


def sub_status() -> dict:
    cfg, w = load_config(), _sub["worker"]
    running = bool(w and _sub["thread"] and _sub["thread"].is_alive())
    return {"enabled": bool(cfg.get("enabled")), "running": running, "main": cfg.get("main", ""), "has_key": bool(cfg.get("key")),
            **({"status": w.status, "beast": w.beast} if running else {})}


def sub_start(main: str | None = None, key: str | None = None, beast: str | None = None) -> dict:
    """Activate: remember the address and key (for the next start of the program) and start the worker thread."""
    cfg = load_config()
    if main:
        cfg["main"] = main.strip()
    if key:
        cfg["key"] = key.strip()
    if not cfg.get("main", "").startswith(("http://", "https://")) or not cfg.get("key"):
        raise ValueError("the address of the main server (http://host:port) and the worker key are needed")
    sub_stop(wait=True)
    cfg["enabled"] = True
    save_config(cfg)
    from . import config
    w = Worker(cfg["main"], cfg["key"], config.DATA_DIR / "cluster" / "worker", beast or cfg.get("beast"))
    t = threading.Thread(target=w.serve, daemon=True, name="cluster-worker")
    _sub["worker"], _sub["thread"] = w, t
    t.start()
    return sub_status()


def sub_stop(wait: bool = False, forget: bool = False) -> dict:
    """Deactivate. A job in progress is left to the main server (it goes back to the queue after the lease)."""
    w, t = _sub["worker"], _sub["thread"]
    if w:
        w.stop.set()
        if wait and t:
            t.join(15)
    _sub["worker"] = _sub["thread"] = None
    cfg = load_config()
    if cfg:
        cfg["enabled"] = False
        if forget:
            cfg.pop("key", None)
        save_config(cfg)
    return sub_status()


def sub_autostart() -> None:
    """At the start of the program: a server that was activated stays a sub-server."""
    cfg = load_config()
    if cfg.get("enabled") and cfg.get("main") and cfg.get("key"):
        try:
            sub_start()
        except Exception:
            _log.warning("could not restart the sub-server worker", exc_info=True)
