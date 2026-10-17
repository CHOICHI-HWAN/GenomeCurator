"""Apache Superset next to the app: is it set up, is it running, start/stop it (it stays a separate program in its own conda env).

The app only launches `scripts/superset.sh run` and watches the port; Superset reads the database read-only (see docs/SUPERSET.md)."""
import atexit
import os
import socket
import subprocess
import time
from pathlib import Path

from . import config, log

_log = log.get("analysis")
ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "superset.sh"
HOME = ROOT / "superset" / "home"
_proc: subprocess.Popen | None = None


def port() -> int:
    return int(os.environ.get("SUPERSET_PORT", "8089"))


def listening() -> bool:
    try:
        with socket.create_connection(("127.0.0.1", port()), timeout=0.5):
            return True
    except OSError:
        return False


def status() -> dict:
    """installed = the repo scripts exist and `init` has been run (metadata DB + admin password); running = something answers on the port."""
    installed = SCRIPT.exists() and (HOME / "admin_password.txt").exists()
    return {"installed": installed, "running": listening(), "port": port(), "url": f"http://127.0.0.1:{port()}/",
            "started_by_app": _proc is not None and _proc.poll() is None, "password_file": str(HOME / "admin_password.txt"),
            "log": str(config.DATA_DIR / "logs" / "superset.log")}


def start(progress=None, wait: int = 120) -> str:
    """Launch Superset (detached from the web request) and wait until it answers. Returns its URL; raises with the reason when it cannot."""
    global _proc
    p = progress or (lambda *a: None)
    st = status()
    if st["running"]:
        p(100, "Superset이 이미 실행 중입니다")
        return st["url"]
    if not st["installed"]:
        raise RuntimeError("Superset이 아직 설치되지 않았습니다. 먼저 `mamba env create -f superset/environment.yml` 와 `scripts/superset.sh init` 을 실행하세요 (docs/SUPERSET.md).")
    lf = config.DATA_DIR / "logs"
    lf.mkdir(parents=True, exist_ok=True)
    out = open(lf / "superset.log", "ab")
    env = {**os.environ, "GC_DATA_DIR": str(config.DATA_DIR), "GC_DB": str(config.DB_PATH)}
    p(10, "Superset 프로세스를 시작하는 중…")
    _proc = subprocess.Popen(["bash", str(SCRIPT), "run"], stdout=out, stderr=out, env=env, start_new_session=True)
    _log.info("superset started (pid %s) on port %s", _proc.pid, port())
    t0 = time.time()
    while time.time() - t0 < wait:
        if _proc.poll() is not None:
            raise RuntimeError(f"Superset이 바로 종료되었습니다 (코드 {_proc.returncode}). 로그: {lf / 'superset.log'}")
        if listening():
            p(100, "Superset이 응답합니다")
            return status()["url"]
        p(10 + 80 * (time.time() - t0) / wait, "Superset이 올라오는 중… (처음엔 1분 가까이 걸립니다)")
        time.sleep(1.5)
    raise RuntimeError(f"{wait}초 안에 Superset이 응답하지 않았습니다. 로그: {lf / 'superset.log'}")


def stop() -> bool:
    """Stop the Superset this app started (one that was started by hand is left alone)."""
    global _proc
    if _proc is None or _proc.poll() is not None:
        return False
    import signal
    try:
        os.killpg(_proc.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        _proc.terminate()
    _log.info("superset stopped (pid %s)", _proc.pid)
    _proc = None
    return True


atexit.register(stop)      # Superset the app launched goes down with the app
