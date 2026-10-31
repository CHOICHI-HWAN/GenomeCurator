"""CPU, memory and disk use of this machine, for the Cluster page.

No extra library. /proc gives the full answer on Linux; macOS has no /proc, so memory comes from `sysctl` and `vm_stat` (cheap, and
cached for a moment because the page asks every few seconds). There is no equally cheap way to read the CPU share on macOS -
`top -l 2` needs two samples and a second of waiting - so `cpu_pct` is left empty there and the load average, which every platform
gives, is reported beside the core count. The page shows whichever of the two it has rather than passing a load off as a percentage.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

_prev: tuple[int, int] | None = None
_mem_cache: dict = {"at": 0.0, "v": None}


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


def _mac_memory() -> tuple[int, int] | None:
    """(total, used) on macOS: `sysctl` for the size, `vm_stat` for what is free. Free, inactive and speculative pages are the
    ones the system can hand back, so 'used' is the rest - the same meaning as MemAvailable on Linux."""
    try:
        total = int(subprocess.run(["sysctl", "-n", "hw.memsize"], capture_output=True, text=True, timeout=5).stdout.strip())
        out = subprocess.run(["vm_stat"], capture_output=True, text=True, timeout=5).stdout
        page = int(re.search(r"page size of (\d+)", out)[1])
        get = lambda k: int(m[1]) * page if (m := re.search(rf"^{k}:\s+(\d+)", out, re.M)) else 0
        free = get("Pages free") + get("Pages inactive") + get("Pages speculative")
        return total, max(0, total - free)
    except (OSError, ValueError, TypeError, AttributeError, subprocess.SubprocessError):
        return None


def memory() -> tuple[int, int] | None:
    """(total, used) bytes; 'used' leaves out what the system can give back (cache)."""
    try:
        m = {l.split(":")[0]: int(l.split()[1]) * 1024 for l in Path("/proc/meminfo").read_text().splitlines() if ":" in l}
        return m["MemTotal"], m["MemTotal"] - m["MemAvailable"]
    except (OSError, KeyError, ValueError):
        pass
    if sys.platform != "darwin":
        return None
    if time.time() - _mem_cache["at"] > 2:                      # the page asks every few seconds; two subprocesses each time is wasteful
        _mem_cache["v"], _mem_cache["at"] = _mac_memory(), time.time()
    return _mem_cache["v"]


def snapshot(path: str | Path = ".") -> dict:
    out: dict = {"cpus": os.cpu_count(), "cpu_pct": cpu_percent(), "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}
    try:
        out["load1"] = round(os.getloadavg()[0], 2)
    except (OSError, AttributeError):
        out["load1"] = None
    mem = memory()
    if mem:
        out["mem_total"], out["mem_used"] = mem
        out["mem_pct"] = round(100.0 * mem[1] / mem[0], 1) if mem[0] else None
    if out["cpu_pct"] is None and out["load1"] is not None and out["cpus"]:
        out["load_pct"] = round(min(100.0, 100.0 * out["load1"] / out["cpus"]), 1)   # no /proc here: the load per core instead, shown as a load, not as a CPU share
    try:
        p = Path(path)
        while not p.exists() and p != p.parent:
            p = p.parent
        d = shutil.disk_usage(p)
        out["disk_total"], out["disk_used"], out["disk_free"] = d.total, d.used, d.free
    except OSError:
        pass
    return out
