"""CPU, memory and disk use of this machine, for the Cluster page (no extra library: /proc on Linux, with a smaller answer elsewhere)."""
from __future__ import annotations

import os
import shutil
import time
from pathlib import Path

_prev: tuple[int, int] | None = None


def _cpu_times() -> tuple[int, int] | None:
    try:
        v = [int(x) for x in Path("/proc/stat").read_text().splitlines()[0].split()[1:]]
        return sum(v), v[3] + (v[4] if len(v) > 4 else 0)       # total, idle (+ iowait)
    except (OSError, ValueError, IndexError):
        return None


def cpu_percent() -> float | None:
    """Share of the CPU in use since the previous call (the first call measures over a short moment)."""
    global _prev
    cur = _cpu_times()
    if cur is None:
        return None
    if _prev is None:
        time.sleep(0.2)
        first, cur = cur, _cpu_times() or cur
        _prev = first
    dt, di = cur[0] - _prev[0], cur[1] - _prev[1]
    _prev = cur
    return round(100.0 * (1 - di / dt), 1) if dt > 0 else None


def memory() -> tuple[int, int] | None:
    """(total, used) bytes; 'used' leaves out what the system can give back (cache)."""
    try:
        m = {l.split(":")[0]: int(l.split()[1]) * 1024 for l in Path("/proc/meminfo").read_text().splitlines() if ":" in l}
        return m["MemTotal"], m["MemTotal"] - m["MemAvailable"]
    except (OSError, KeyError, ValueError):
        return None


def snapshot(path: str | Path = ".") -> dict:
    out: dict = {"cpus": os.cpu_count(), "cpu_pct": cpu_percent(), "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}
    try:
        out["load1"] = round(os.getloadavg()[0], 2)
    except (OSError, AttributeError):
        out["load1"] = None
    mem = memory()
    if mem:
        out["mem_total"], out["mem_used"] = mem
    try:
        p = Path(path)
        while not p.exists() and p != p.parent:
            p = p.parent
        d = shutil.disk_usage(p)
        out["disk_total"], out["disk_used"], out["disk_free"] = d.total, d.used, d.free
    except OSError:
        pass
    return out
