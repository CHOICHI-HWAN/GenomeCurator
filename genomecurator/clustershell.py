"""Web terminal of the cluster: the main server opens a shell on a sub-server in the browser (Linux and macOS sub-servers).

* A sub-server allows it only when its own page has switched it on (off by default; `worker.json` remembers the choice).
* No password: the main server makes a one-minute, one-time token signed with the SHA-256 of the sub-server's worker key. The main server keeps only that
  hash, the sub-server can compute it from its key, and nobody else has either. The token names the sub-server, so it is useless on another one.
* The shell runs in a pseudo-terminal (top, vim, tail -f work) as the user who started the sub-server program, in its home folder.
* Every session is written to the audit log of both servers (who, when, from where, how long); what is typed is not recorded.
"""
from __future__ import annotations

import base64
import hashlib
import hmac
import json
import os
import secrets
import time

from . import log

_log = log.get("shell")
TOKEN_SECONDS = 60
_used: dict[str, float] = {}           # nonces already used (kept until they would have expired anyway)
sessions: dict[str, dict] = {}         # open sessions on this sub-server: id -> {"from", "since"}


def supported() -> bool:
    return os.name == "posix"


def _b64(b: bytes) -> str:
    return base64.urlsafe_b64encode(b).decode().rstrip("=")


def _unb64(s: str) -> bytes:
    return base64.urlsafe_b64decode(s + "=" * (-len(s) % 4))


def make_token(key_sha: str, node_id: str, who: str = "main") -> str:
    """On the main server: a token for one shell session on the sub-server `node_id` (signed with that sub-server's key hash)."""
    body = _b64(json.dumps({"n": node_id, "exp": int(time.time()) + TOKEN_SECONDS, "nonce": secrets.token_urlsafe(12), "who": who}).encode())
    sig = _b64(hmac.new(key_sha.encode(), body.encode(), hashlib.sha256).digest())
    return f"{body}.{sig}"


def check_token(token: str, key: str | None, node_name: str | None = None) -> dict:
    """On the sub-server: is this a fresh, unused token signed for us? -> its content, or ValueError."""
    if not key:
        raise ValueError("this server has not been approved by a main server")
    try:
        body, sig = token.split(".", 1)
        want = _b64(hmac.new(hashlib.sha256(key.encode()).hexdigest().encode(), body.encode(), hashlib.sha256).digest())
        if not hmac.compare_digest(want, sig):
            raise ValueError("bad signature")
        t = json.loads(_unb64(body))
    except ValueError:
        raise ValueError("the token is not valid for this server")
    now = time.time()
    for n, exp in list(_used.items()):
        if exp < now:
            _used.pop(n, None)
    if t.get("exp", 0) < now:
        raise ValueError("the token has expired; open the terminal again from the main server")
    if t.get("nonce") in _used:
        raise ValueError("the token was used already")
    _used[t["nonce"]] = t["exp"]
    return t


_CHILD = ("import os, sys, fcntl, termios\n"
          "os.setsid()\n"
          "fcntl.ioctl(0, termios.TIOCSCTTY, 0)\n"          # the pseudo-terminal becomes the shell's controlling terminal (job control works)
          "sh = sys.argv[1]\n"
          "try:\n    os.execvp(sh, [sh, '-l'])\nexcept OSError:\n    os.execvp('/bin/sh', ['/bin/sh'])\n")


def spawn(cols: int = 100, rows: int = 30):
    """Start a login shell in a pseudo-terminal. -> (process, master fd). Started through subprocess (not fork) because the server runs threads."""
    import pty
    import subprocess
    import sys
    master, slave = pty.openpty()
    resize(master, cols, rows)
    home = os.path.expanduser("~")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GC_", "UVICORN"))}
    env.update({"TERM": "xterm-256color", "HOME": home, "GC_WEBSHELL": "1"})
    proc = subprocess.Popen([sys.executable, "-c", _CHILD, os.environ.get("SHELL") or "/bin/bash"], stdin=slave, stdout=slave, stderr=slave,
                            cwd=home if os.path.isdir(home) else None, env=env, close_fds=True)
    os.close(slave)
    return proc, master


def resize(fd: int, cols: int, rows: int) -> None:
    import fcntl
    import struct
    import termios
    try:
        fcntl.ioctl(fd, termios.TIOCSWINSZ, struct.pack("HHHH", max(2, min(rows, 500)), max(10, min(cols, 1000)), 0, 0))
    except OSError:
        pass


def end(proc, fd: int) -> None:
    """Hang up the shell and everything it started (its own session), then close the terminal."""
    import signal
    import subprocess
    for sig in (signal.SIGHUP, signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(proc.pid, sig)
        except (ProcessLookupError, PermissionError):
            break
        try:
            proc.wait(1.5)
            break
        except subprocess.TimeoutExpired:
            continue
    try:
        os.close(fd)
    except OSError:
        pass
