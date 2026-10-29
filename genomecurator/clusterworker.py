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

from . import log, sysinfo

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
    import sys
    from . import __version__
    caps = {"kinds": [], "cpus": os.cpu_count(), "platform": platform.platform(), "host": platform.node(), "python": sys.version.split()[0], "app_version": __version__, "tools": {}}
    mem = sysinfo.memory()
    if mem:
        caps["mem_total"] = mem[0]
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


MIN_FREE_GB = float(os.environ.get("GC_WORKER_MIN_FREE_GB", "10"))
ENV_EVERY = 600                                  # seconds between environment checks while the worker runs


def env_check(beast: str | None, workdir: str | Path) -> dict:
    """Can this machine run the jobs? BEAST starts, the BEAGLE library is next to it, and the work folder has room. -> {"ok", "checks": [...]}"""
    checks = []
    ver = None
    if beast:
        try:
            r = subprocess.run([beast, "-version"], capture_output=True, text=True, timeout=120)
            import re
            m = re.search(r"v([0-9][0-9.]*)", r.stdout + r.stderr)
            ver = m[1] if m else None
        except Exception:
            ver = None
    checks.append({"name": "BEAST", "ok": bool(ver), "detail": f"BEAST {ver}" if ver else ("found but does not start" if beast else "not installed")})
    libs = []
    if beast:
        pre = Path(beast).resolve().parent.parent
        libs = list((pre / "lib").glob("libhmsbeagle*")) if (pre / "lib").is_dir() else []
        if not libs:
            try:
                from . import beastprep
                ep = beastprep.env_prefix()
                libs = list((ep / "lib").glob("libhmsbeagle*")) if ep else []
            except Exception:
                libs = []
    checks.append({"name": "BEAGLE", "ok": bool(libs), "detail": "library found" if libs else "the BEAGLE library BEAST 10 needs is missing"})
    try:
        p = Path(workdir)
        while not p.exists() and p != p.parent:
            p = p.parent
        free = shutil.disk_usage(p).free / 2**30
        checks.append({"name": "disk", "ok": free >= MIN_FREE_GB, "detail": f"{free:.0f} GB free (at least {MIN_FREE_GB:.0f} GB needed)"})
    except OSError as e:
        checks.append({"name": "disk", "ok": False, "detail": str(e)})
    return {"ok": all(c["ok"] for c in checks), "checks": checks, "beast_version": ver}


def _end(proc: subprocess.Popen, grace: float = 30) -> None:
    """Stop a job's process and everything it started (SIGTERM to the group, SIGKILL after `grace` seconds)."""
    import signal
    def group(sig):
        try:
            if os.name == "posix":
                os.killpg(proc.pid, sig)
            else:
                proc.terminate() if sig != getattr(signal, "SIGKILL", None) else proc.kill()
        except (ProcessLookupError, PermissionError):
            pass
    group(signal.SIGTERM)
    try:
        proc.wait(grace)
    except subprocess.TimeoutExpired:
        group(getattr(signal, "SIGKILL", signal.SIGTERM))
        proc.wait(10)


class Worker:
    def __init__(self, base: str, key: str, workdir: str | Path, beast: str | None = None, poll: float = 10.0, beat: float = 20.0, retry: float = 15.0):
        self.main, self.work = Main(base, key), Path(workdir)
        self.beast, self.poll, self.beat, self.retry = find_beast(beast), poll, beat, retry
        self.stop = threading.Event()
        self.status = {"state": "starting", "job": None, "error": None, "since": time.strftime("%Y-%m-%dT%H:%M:%S")}
        self.history: list[dict] = []
        self.paused = threading.Event()          # Stop button of the sub-server page: give the running job back, take no new work           # the last jobs this worker ran (shown on its page)
        self.work.mkdir(parents=True, exist_ok=True)

    def metrics(self) -> dict:
        try:
            m = sysinfo.snapshot(self.work)
        except Exception:
            m = {}
        m["paused"] = self.paused.is_set()
        return m

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
        self._outcome = None
        self._set("running", jid)
        self.status["label"] = job.get("label") or jid
        self._t0 = time.time()
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
            # its own process group: `beast` is a shell script that starts Java, and stopping must end both
            proc = subprocess.Popen(cmd, cwd=str(jd / "run"), env=env, stdout=lf, stderr=subprocess.STDOUT, text=True, **({"start_new_session": True} if os.name == "posix" else {}))
            last_beat = 0.0
            while proc.poll() is None and not self.stop.is_set():
                time.sleep(1)
                if self.paused.is_set():
                    break
                if time.time() - last_beat < self.beat:
                    continue
                last_beat = time.time()
                pct, msg = self._progress(jd, int(spec.get("chain_length") or 0))
                self.status.update({"pct": pct, "msg": msg})
                try:
                    if not self.main.json("POST", f"/api/cluster/w/jobs/{jid}/heartbeat", {"pct": pct, "msg": msg, "metrics": self.metrics()}).get("continue", False):
                        gone.set()
                except Gone:
                    gone.set()
                except Exception:
                    pass                            # the main server is away: the job goes on
                if gone.is_set():
                    break
            if proc.poll() is None:
                _end(proc)
        if gone.is_set():
            return self._discard(jd)
        if self.paused.is_set():
            self._outcome = "stopped"
            try:
                self.main.json("POST", f"/api/cluster/w/jobs/{jid}/release", {})
            except Exception:
                pass                                # the main server hands it out again after the grace period anyway
            return self._discard(jd)
        if self.stop.is_set():
            return                                  # stopped on request: the job is left to the lease of the main server
        if proc.returncode != 0:
            tail = "".join(open(jd / "run" / "beast.screen.log", encoding="utf-8").readlines()[-6:])[-300:]
            return self._fail(jd, jid, f"BEAST exited with code {proc.returncode}: {tail.strip()}")
        outs = [{"name": n, "sha256": sha256(jd / "run" / n), "bytes": (jd / "run" / n).stat().st_size} for n in KINDS["beast"]["outputs"] if (jd / "run" / n).is_file()]
        self._state(jd, phase="finished", outputs=outs, seconds=round(time.time() - t0))
        self._outcome = "done"
        self.deliver(jd)

    def _fail(self, jd: Path, jid: str, why: str) -> None:
        self._outcome = "failed"
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
        if getattr(self, "_outcome", None) not in ("failed", "done"):
            self._outcome = "stopped"
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
        known, env_at = False, 0.0
        while not self.stop.is_set():
            try:
                if time.time() - env_at > ENV_EVERY:
                    env_at = time.time()
                    env = env_check(self.beast, self.work)
                    if env["ok"] != caps.get("env", {}).get("ok"):
                        known = False                   # tell the main server about the change
                    caps["env"] = env
                if not known:
                    self.main.json("POST", "/api/cluster/w/hello", {**caps, "metrics": self.metrics()})
                    known = True
                self.resume()
                if self.paused.is_set():
                    if self.status["state"] != "paused":
                        self._set("paused")
                    self.main.json("POST", "/api/cluster/w/claim", {"kinds": [], "metrics": self.metrics()})
                    self.stop.wait(self.poll)
                    continue
                if not caps["env"]["ok"]:
                    # the environment broke after approval: take no work, keep reporting
                    bad = ", ".join(f"{c['name']}: {c['detail']}" for c in caps["env"]["checks"] if not c["ok"])
                    if self.status["state"] != "env_broken":
                        self._set("env_broken", error=bad)
                    if self.main.json("POST", "/api/cluster/w/claim", {"kinds": [], "metrics": self.metrics()}).get("recheck"):
                        env_at = 0.0
                    self.stop.wait(self.poll)
                    continue
                got = self.main.json("POST", "/api/cluster/w/claim", {"kinds": caps["kinds"], "metrics": self.metrics()})
                if got.get("recheck"):
                    env_at = 0.0                        # the main server asked for a fresh environment check
                    continue
                job = got.get("job")
                if job:
                    self.run_job(job)
                    self.history = ([{"id": job["id"], "label": job.get("label") or job["id"], "seconds": round(time.time() - getattr(self, "_t0", time.time())),
                                      "result": getattr(self, "_outcome", "?"), "ended": time.strftime("%Y-%m-%dT%H:%M:%S")}] + self.history)[:10]
                    self._set("idle")
                    continue
                if self.status["state"] != "idle":
                    self._set("idle")
            except Gone as e:
                _log.error("the main server refused this worker (%s); check the key", e)
                self._set("refused", error="the main server did not accept this worker (revoked or unknown key)")
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


_sub: dict = {"worker": None, "thread": None, "join": None}


def _workdir() -> Path:
    from . import config
    return config.DATA_DIR / "cluster" / "worker"


_env_cache: dict = {"at": 0.0, "env": None}


def local_env(beast: str | None = None, fresh: bool = False) -> dict:
    """The environment check of this server (kept for a minute: the page asks every few seconds and starting BEAST takes a moment)."""
    if fresh or not _env_cache["env"] or time.time() - _env_cache["at"] > 60:
        _env_cache["env"], _env_cache["at"] = env_check(find_beast(beast or load_config().get("beast")), _workdir()), time.time()
    return _env_cache["env"]


def sub_status() -> dict:
    cfg, w = load_config(), _sub["worker"]
    running = bool(w and _sub["thread"] and _sub["thread"].is_alive())
    jt = _sub["join"]
    return {"enabled": bool(cfg.get("enabled")), "running": running, "paused": bool(cfg.get("paused")), "main": cfg.get("main", ""), "has_key": bool(cfg.get("key")),
            "join": {"state": cfg.get("join_state"), "id": cfg.get("join_id"), "error": cfg.get("join_error"), "waiting": bool(jt and jt.is_alive())},
            **({"status": w.status, "beast": w.beast, "history": w.history} if running else {})}


def _start_worker(cfg: dict) -> None:
    w = Worker(cfg["main"], cfg["key"], _workdir(), cfg.get("beast"))
    if cfg.get("paused"):
        w.paused.set()
    t = threading.Thread(target=w.serve, daemon=True, name="cluster-worker")
    _sub["worker"], _sub["thread"] = w, t
    t.start()


def sub_request(main: str, name: str | None = None) -> dict:
    """Ask the main server to take this server on as a worker. Only a server whose environment check passes may ask; the worker starts by itself
    once the main server has approved (the key is handed over then, nobody copies it)."""
    import secrets as _secrets
    main = (main or "").strip().rstrip("/")
    if not main.startswith(("http://", "https://")):
        raise ValueError("the address of the main server, e.g. http://10.0.0.1:8000")
    env = local_env(fresh=True)
    if not env["ok"]:
        raise ValueError("the environment check failed: " + "; ".join(f"{c['name']} ({c['detail']})" for c in env["checks"] if not c["ok"]))
    sub_stop(wait=True)
    secret = _secrets.token_urlsafe(32)
    caps = capabilities(find_beast(load_config().get("beast")))
    info = {**caps, "env": env, "metrics": sysinfo.snapshot(_workdir())}
    host = name or platform.node() or "worker"
    try:
        r = Main(main, "").json("POST", "/api/cluster/w/join", {"name": host, "secret": secret, "info": info})
    except Gone as e:
        raise ValueError(f"the main server refused the request ({e})")
    except Exception as e:
        raise ValueError(f"the main server cannot be reached at {main} ({e})")
    cfg = load_config()
    cfg.update({"main": main, "join_id": r["id"], "join_secret": secret, "join_state": "pending", "join_error": None, "enabled": False})
    cfg.pop("key", None)
    save_config(cfg)
    _wait_for_approval()
    return sub_status()


def _wait_for_approval(every: float = 5.0) -> None:
    if _sub["join"] and _sub["join"].is_alive():
        return
    stop = threading.Event()
    _sub["join_stop"] = stop

    def loop():
        while not stop.is_set():
            cfg = load_config()
            if not cfg.get("join_id") or cfg.get("join_state") != "pending":
                return
            try:
                r = Main(cfg["main"], "").json("POST", f"/api/cluster/w/join/{cfg['join_id']}", {"secret": cfg["join_secret"]})
            except Gone:
                r = {"state": "unknown"}
            except Exception as e:
                cfg["join_error"] = f"main server not reachable ({str(e)[:120]})"
                save_config(cfg)
                stop.wait(every * 3)
                continue
            if r.get("state") == "approved" and r.get("key"):
                cfg.update({"key": r["key"], "join_state": "approved", "join_error": None, "enabled": True, "name": r.get("name")})
                cfg.pop("join_secret", None)
                save_config(cfg)
                _log.info("approved by the main server as %s; starting the worker", r.get("name"))
                _start_worker(cfg)
                return
            if r.get("state") in ("rejected", "unknown", "delivered"):
                cfg.update({"join_state": r["state"], "join_error": None})
                cfg.pop("join_secret", None)
                save_config(cfg)
                return
            cfg["join_error"] = None
            save_config(cfg)
            stop.wait(every)
    t = threading.Thread(target=loop, daemon=True, name="cluster-join")
    _sub["join"] = t
    t.start()


def sub_start(main: str | None = None, key: str | None = None, beast: str | None = None) -> dict:
    """Start with a key that is already known (handed over at approval, or typed for a worker added by hand on the main server)."""
    cfg = load_config()
    if main:
        cfg["main"] = main.strip().rstrip("/")
    if key:
        cfg["key"] = key.strip()
    if beast:
        cfg["beast"] = beast
    if not cfg.get("main", "").startswith(("http://", "https://")) or not cfg.get("key"):
        raise ValueError("this server has not been approved by a main server yet")
    env = env_check(find_beast(cfg.get("beast")), _workdir())
    if not env["ok"]:
        raise ValueError("the environment check failed: " + "; ".join(f"{c['name']} ({c['detail']})" for c in env["checks"] if not c["ok"]))
    sub_stop(wait=True)
    cfg["enabled"] = True
    save_config(cfg)
    _start_worker(cfg)
    return sub_status()


def sub_pause(on: bool) -> dict:
    """Stop (on=True): the running job is stopped and given back to the main server at once, and no new work is taken until Resume. Kept across restarts."""
    cfg = load_config()
    cfg["paused"] = bool(on)
    save_config(cfg)
    w = _sub["worker"]
    if w:
        (w.paused.set if on else w.paused.clear)()
    return sub_status()


def sub_stop(wait: bool = False, forget: bool = False) -> dict:
    """Deactivate (a pending request is withdrawn too). A job in progress is left to the main server (it goes back to the queue after the lease)."""
    w, t = _sub["worker"], _sub["thread"]
    if w:
        w.stop.set()
        if wait and t:
            t.join(15)
    if _sub.get("join_stop"):
        _sub["join_stop"].set()
    _sub["worker"] = _sub["thread"] = None
    cfg = load_config()
    if cfg:
        cfg["enabled"] = False
        if cfg.get("join_state") == "pending":
            cfg["join_state"] = "withdrawn"
            cfg.pop("join_secret", None)
        if forget:
            for k in ("key", "join_id", "join_secret", "join_state"):
                cfg.pop(k, None)
        save_config(cfg)
    return sub_status()


def sub_autostart() -> None:
    """At the start of the program: an approved server stays a worker, and a request still waiting goes on waiting."""
    cfg = load_config()
    try:
        if cfg.get("enabled") and cfg.get("main") and cfg.get("key"):
            sub_start()
        elif cfg.get("join_state") == "pending" and cfg.get("join_secret"):
            _wait_for_approval()
    except Exception:
        _log.warning("could not restart the sub-server worker", exc_info=True)


# ---- "install": make the environment a sub-server needs (BEAST 10 and its BEAGLE library) with one button
BEAST_PACKAGES = ["beast", "beagle-lib"]
_install: dict = {"job": None}


def install_state() -> dict:
    """What is installed on this machine for the job kinds, and the state of the install running now (if any)."""
    from . import beastprep, jobs
    st = beastprep.status()
    pkgmgr = shutil.which("mamba") or shutil.which("micromamba") or shutil.which("conda")
    job = jobs.get(_install["job"]) if _install["job"] else None
    return {"beast": {"ok": st["ok"], "version": st["version"]}, "installer": pkgmgr, "job": job and {k: job[k] for k in ("id", "status", "pct", "msg", "error", "log")}}


def install_beast() -> str:
    """Create the conda environment `gc-beast` (the one the app looks for) in a background job; the progress is the job's log. Returns the job id."""
    from . import beastprep, jobs
    exe = shutil.which("mamba") or shutil.which("micromamba") or shutil.which("conda")
    if not exe:
        raise RuntimeError("mamba or conda is needed to install BEAST (install miniforge first)")
    if beastprep.status()["ok"]:
        raise ValueError("BEAST is already installed")
    if _install["job"]:
        cur = jobs.get(_install["job"])
        if cur and not cur["done"]:
            raise ValueError("an install is already running")

    def work(progress):
        progress(2, "creating the environment gc-beast (BEAST 10 and BEAGLE), this takes a few minutes")
        cmd = [exe, "create", "-y", "-n", "gc-beast", "-c", "conda-forge", "-c", "bioconda"] + BEAST_PACKAGES
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        n = 0
        for line in p.stdout:
            line = line.rstrip()
            if line:
                n += 1
                progress(min(95, 2 + n * 0.5), line[-160:])
        if p.wait() != 0:
            raise RuntimeError(f"{Path(exe).name} exited with code {p.returncode}; see the log")
        st = beastprep.status()
        if not st["ok"]:
            raise RuntimeError("the environment was created but BEAST was not found in it")
        local_env(fresh=True)
        progress(100, f"BEAST {st['version']} installed")
        return "/cluster/sub"
    jid = jobs.start(work, "BEAST environment install", "install")
    _install["job"] = jid
    return jid
