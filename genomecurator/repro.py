"""A reproducibility manifest for a figure or table: what it was made from and with what.

`manifest()` collects the program version and git commit, the Python and library versions, the versions of the external tools that made the inputs (BEAST,
TreeAnnotator when their environment exists), the SHA-256 and size of every input file, the parameters of the figure and the SHA-256 of the files written, and
writes it as `<figure>.manifest.json` next to the figure. Anyone with the inputs can then check that they hold the same files and re-run the same call.

Nothing here reads genome content; the manifest names files and numbers only. Paths are written as given (relative paths stay relative), so a manifest can be
shared without the folder layout of the machine it was made on when the caller passes relative paths.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

from . import __version__


def sha256(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def git_state(root: str | Path | None = None) -> dict:
    """The commit of the program and whether files were modified since (uncommitted changes make a run not exactly reproducible from the commit alone)."""
    root = Path(root or Path(__file__).resolve().parent.parent)
    def run(*a):
        try:
            return subprocess.run(["git", *a], cwd=root, capture_output=True, text=True, timeout=20).stdout.strip()
        except Exception:
            return ""
    commit = run("rev-parse", "HEAD")
    return {"commit": commit or None, "dirty": bool(run("status", "--porcelain", "--untracked-files=no")) if commit else None}


def library_versions(names=("numpy", "matplotlib", "pyyaml", "openpyxl")) -> dict:
    from importlib import metadata
    out = {}
    for n in names:
        try:
            out[n] = metadata.version(n)
        except metadata.PackageNotFoundError:
            pass
    return out


def tool_versions() -> dict:
    """Versions of the external tools of the BEAST environment, when it is installed (no error when it is not)."""
    out = {}
    try:
        from . import beastprep
        pre = beastprep.env_prefix()
        if pre:
            for name, exe in (("beast", "beast"), ("treeannotator", "treeannotator")):
                try:
                    r = subprocess.run([str(pre / "bin" / exe), "-version"], capture_output=True, text=True, timeout=60)
                    import re
                    m = re.search(r"v([0-9][0-9.]*)", r.stdout + r.stderr)
                    if m:
                        out[name] = m[1]
                except Exception:
                    pass
    except Exception:
        pass
    return out


def manifest(outputs: list[str | Path], inputs: list[str | Path] = (), parameters: dict | None = None, extra: dict | None = None, write: bool = True) -> dict:
    """Build (and by default write, next to the first output) the manifest of `outputs` made from `inputs` with `parameters`."""
    ins = []
    for p in inputs:
        p = Path(p)
        ins.append({"file": str(p), "sha256": sha256(p), "bytes": p.stat().st_size} if p.exists() else {"file": str(p), "missing": True})
    outs = [{"file": str(Path(o)), "sha256": sha256(o), "bytes": Path(o).stat().st_size} for o in outputs if Path(o).exists()]
    m = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "program": {"name": "genomecurator", "version": __version__, **git_state()},
        "python": sys.version.split()[0], "platform": platform.platform(),
        "libraries": library_versions(), "tools": tool_versions(),
        "inputs": ins, "parameters": parameters or {}, "outputs": outs, **(extra or {}),
    }
    if write and outputs:
        first = Path(outputs[0])
        first.with_suffix(".manifest.json").write_text(json.dumps(m, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")
    return m


def verify(manifest_file: str | Path, base: str | Path = ".") -> list[str]:
    """Check the inputs and outputs named in a manifest against the files on disk. -> list of problems (empty when everything matches)."""
    m = json.loads(Path(manifest_file).read_text(encoding="utf-8"))
    bad = []
    for kind in ("inputs", "outputs"):
        for e in m.get(kind, []):
            f = Path(base) / e["file"]
            if e.get("missing"):
                bad.append(f"{kind[:-1]} {e['file']}: was missing when the manifest was made")
            elif not f.exists():
                bad.append(f"{kind[:-1]} {e['file']}: not found")
            elif sha256(f) != e["sha256"]:
                bad.append(f"{kind[:-1]} {e['file']}: content differs")
    return bad
