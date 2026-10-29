"""Netdata (netdata.cloud, GPL-3.0 agent) as an optional source of machine metrics for the Cluster pages.

When a server runs a Netdata agent, its API gives per-second history (CPU, memory, disk and network traffic) that this program does not collect itself.
Only the agent's HTTP API is read; nothing of Netdata is bundled or installed by this program. Without Netdata the pages fall back to the current values
of `sysinfo` and show no history.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request

PORT = 19999
LOCAL = f"http://127.0.0.1:{PORT}"


def _get(url: str, timeout: float = 3.0) -> dict:
    with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=timeout) as r:
        return json.loads(r.read())


def probe(base: str = LOCAL, timeout: float = 2.0) -> dict:
    """Is a Netdata agent answering at `base`? -> {"ok", "version", "claimed" (connected to Netdata Cloud)}."""
    try:
        d = _get(f"{base}/api/v1/info", timeout)
        return {"ok": True, "version": d.get("version"), "claimed": bool(d.get("agent-claimed")), "port": urllib.parse.urlparse(base).port or PORT}
    except (OSError, ValueError, urllib.error.URLError):
        return {"ok": False}


def _series(base: str, chart: str, minutes: int, points: int, host: str | None, timeout: float) -> tuple[list[str], list[list[float]]]:
    prefix = f"{base}/host/{urllib.parse.quote(host)}" if host else base
    q = urllib.parse.urlencode({"chart": chart, "after": -60 * minutes, "points": points, "group": "average", "format": "json", "options": "seconds"})
    d = _get(f"{prefix}/api/v1/data?{q}", timeout)
    rows = sorted(d.get("data") or [], key=lambda r: r[0])
    return d.get("labels") or [], rows


def history(base: str = LOCAL, minutes: int = 60, points: int = 120, host: str | None = None, timeout: float = 4.0) -> dict:
    """CPU use (%), memory use (%), load, disk read/write (MB/s) and network in/out (Mbit/s) over the last `minutes`, oldest first.
    `host` reads a child that streams to the Netdata agent at `base` (a parent)."""
    minutes = max(5, min(int(minutes), 7 * 24 * 60))
    out: dict = {"ok": True, "minutes": minutes, "source": base + (f"/host/{host}" if host else "")}
    lab, rows = _series(base, "system.cpu", minutes, points, host, timeout)
    out["t"] = [r[0] for r in rows]
    out["cpu"] = [round(sum(x or 0 for x in r[1:]), 1) for r in rows]

    def take(chart, fn):
        try:
            lab, rows = _series(base, chart, minutes, points, host, timeout)
            idx = {n: i for i, n in enumerate(lab)}
            by_t = {r[0]: fn(idx, r) for r in rows}
            return [by_t.get(t) for t in out["t"]]
        except (OSError, ValueError, urllib.error.URLError, ZeroDivisionError, KeyError, IndexError, TypeError, StopIteration):
            return None

    def ram(idx, r):
        tot = sum(abs(r[i]) for n, i in idx.items() if n != "time" and r[i] is not None)
        return round(100 * abs(r[idx["used"]]) / tot, 1) if tot else None
    out["mem"] = take("system.ram", ram)
    out["load1"] = take("system.load", lambda idx, r: r[idx["load1"]])
    def dim(idx, *names):                       # dimension names differ between Netdata versions (in/out, reads/writes)
        return next(idx[n] for n in names if n in idx)
    out["disk_read"] = take("system.io", lambda idx, r: round(abs(r[dim(idx, "reads", "in")] or 0) / 1024, 2))       # KiB/s -> MiB/s
    out["disk_write"] = take("system.io", lambda idx, r: round(abs(r[dim(idx, "writes", "out")] or 0) / 1024, 2))
    out["net_in"] = take("system.net", lambda idx, r: round(abs(r[dim(idx, "received", "InOctets")] or 0) / 1000, 2))   # kbit/s -> Mbit/s
    out["net_out"] = take("system.net", lambda idx, r: round(abs(r[dim(idx, "sent", "OutOctets")] or 0) / 1000, 2))
    return out


def history_or_none(*a, **kw) -> dict:
    try:
        return history(*a, **kw)
    except (OSError, ValueError, urllib.error.URLError, KeyError, IndexError, TypeError) as e:
        return {"ok": False, "error": str(e)[:200]}


# ---- streaming: sub-servers send their metrics to the Netdata of the main server (a "parent"), so one place holds all history, without Netdata Cloud
def stream_key() -> str:
    """The API key the parent accepts (made once, kept in app_meta)."""
    import uuid
    from .db import connect
    con = connect()
    try:
        r = con.execute("SELECT value FROM app_meta WHERE key='netdata_stream_key'").fetchone()
        if r:
            return r[0]
        k = str(uuid.uuid4())
        con.execute("INSERT INTO app_meta(key, value) VALUES('netdata_stream_key', ?)", (k,))
        con.commit()
        return k
    finally:
        con.close()


def address_towards(ip: str | None) -> str | None:
    """The address of this machine on the link that reaches `ip` (no packet is sent)."""
    import socket
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            s.connect((ip or "10.255.255.255", 9))
            return s.getsockname()[0]
        finally:
            s.close()
    except OSError:
        return None


def parent_conf(key: str, allow: str = "*") -> str:
    """`allow` limits which addresses may stream (e.g. "10.0.0.*" for the direct link)."""
    return (f"# on the MAIN server: /etc/netdata/stream.conf (in Docker: the netdataconfig volume), then restart Netdata\n"
            f"[{key}]\n    enabled = yes\n    default memory mode = dbengine\n    health enabled by default = auto\n    allow from = {allow}\n")


def child_conf(key: str, parent: str) -> str:
    return (f"# on each SUB-server: /etc/netdata/stream.conf, then restart Netdata\n"
            f"[stream]\n    enabled = yes\n    destination = {parent}:{PORT}\n    api key = {key}\n")


def streaming_hosts(base: str = LOCAL) -> list[str]:
    """Hosts whose metrics this agent holds (itself and the children streaming to it)."""
    try:
        return list(_get(f"{base}/api/v1/info").get("mirrored_hosts") or [])
    except (OSError, ValueError, urllib.error.URLError):
        return []
