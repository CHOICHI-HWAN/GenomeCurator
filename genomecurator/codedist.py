"""Sending the program's own code from the main server to the workers that allow it (see cluster.py, clusterworker.py).

The main server is already the single source of the data: a worker holds no genomes of its own, copies what a job needs and
deletes its copy once the main server has the result. The code was the exception. Three machines ran three checkouts, nobody
knew which, and the cluster page could only report that two of them were running something it could not name. A result is
only reproducible if you can say what produced it, so the code has to be as controlled as the data.

So: the main server publishes the package **at its HEAD commit** as one tar.gz with its SHA-256, and a worker that has been
allowed to do so fetches it, checks it, unpacks it into a folder of its own named after the commit, and restarts itself on it.

What makes this safe enough to do at all:
  * It is off until the person switches it on **on the worker's own machine**. Running code someone else sends you is the
    same power as the web terminal, and it is granted the same way: deliberately, once, on the machine that grants it.
  * Only what `git archive HEAD` holds travels - committed code, no data, no settings, no keys. Uncommitted changes on the
    main server are not distributed, and the manifest says so.
  * The bundle is checked against the SHA-256 the manifest declares before anything is unpacked, and unpacking takes regular
    files under `genomecurator/` only - no links, no absolute paths, no `..`, no device files, nothing above the size cap.
  * Nothing is overwritten. Each commit lands in `cluster/code/<commit>/`, and switching means pointing at another folder, so
    going back is pointing back. A tree that cannot even be imported is thrown away and the worker keeps the code it has.
  * A worker swaps code only when it is running nothing, so no analysis is ever killed by an update.
"""
from __future__ import annotations

import hashlib
import os
import subprocess
import tarfile
import tempfile
from pathlib import Path

from . import build, log

_log = log.get("codedist")

MAX_BUNDLE = 96 << 20          # the package is about 2 MB; this only has to be smaller than "a disk"
MAX_FILE = 16 << 20
PREFIX = "genomecurator/"      # what a bundle may contain, and nothing else


def _root() -> Path:
    return Path(__file__).resolve().parent.parent


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bundle_path(commit: str) -> Path:
    from . import config
    return config.DATA_DIR / "cluster" / "code" / f"{commit}.tar.gz"


def bundle() -> dict:
    """The package at this server's HEAD as one tar.gz -> its manifest. Built once per commit and kept, so asking costs nothing."""
    b = build()
    commit = b.get("commit")
    if not commit:
        raise RuntimeError("this server's code is not in a git checkout, so there is no commit to send")
    out = bundle_path(commit)
    if not out.is_file():
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_suffix(".tmp")
        r = subprocess.run(["git", "-C", str(_root()), "archive", "--format=tar.gz", "-o", str(tmp), commit, "genomecurator"],
                           capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            tmp.unlink(missing_ok=True)
            raise RuntimeError(f"the code could not be packed: {r.stderr.strip()[:200]}")
        os.replace(tmp, out)
        _log.info("packed the code of %s (%d bytes)", commit, out.stat().st_size)
    return {"commit": commit, "version": b.get("version"), "sha256": sha256(out), "size": out.stat().st_size,
            "dirty": bool(b.get("dirty"))}


def manifest() -> dict:
    """What a worker is told. A failure here is an answer, not an exception: a worker must go on working with the code it has."""
    try:
        return {"ok": True, **bundle()}
    except Exception as ex:
        _log.warning("the code bundle could not be made: %s", ex)
        return {"ok": False, "error": str(ex)[:200], "commit": build().get("commit"), "version": build().get("version")}


# ---- the worker's side

def _safe_members(tar: tarfile.TarFile):
    """Only the package's own regular files. Anything else in the archive is a reason to refuse the whole bundle, not to skip a file."""
    for m in tar.getmembers():
        name = m.name
        if name in (".", "./"):
            continue
        if name.startswith("./"):
            name = name[2:]
        p = Path(name)
        if p.is_absolute() or ".." in p.parts:
            raise ValueError(f"the bundle holds a path outside itself: {m.name}")
        if not (name == PREFIX.rstrip("/") or name.startswith(PREFIX)):
            raise ValueError(f"the bundle holds something that is not the package: {m.name}")
        if m.issym() or m.islnk() or m.isdev() or m.isfifo():
            raise ValueError(f"the bundle holds a link or a device: {m.name}")
        if not (m.isfile() or m.isdir()):
            raise ValueError(f"the bundle holds an unexpected kind of entry: {m.name}")
        if m.size > MAX_FILE:
            raise ValueError(f"the bundle holds a file that is too large: {m.name}")
        m.name = name
        m.mode = 0o755 if m.isdir() else 0o644
        yield m


def unpack(tar_path: Path, dest: Path, commit: str) -> Path:
    """Unpack a verified bundle into `dest` (made fresh). -> the folder to put on PYTHONPATH.

    A `BUILD` file goes in beside the package, because an unpacked bundle is not a git checkout and would otherwise report no
    commit at all - which is the very blindness this is here to end (`genomecurator.build`)."""
    import shutil
    if tar_path.stat().st_size > MAX_BUNDLE:
        raise ValueError("the code bundle is larger than this program ever is")
    staging = dest.with_name(dest.name + ".new")
    shutil.rmtree(staging, ignore_errors=True)
    staging.mkdir(parents=True)
    with tarfile.open(tar_path, "r:gz") as t:
        # `filter="data"` is Python 3.14's default and refuses the dangerous kinds of entry on its own; `_safe_members` is
        # the check that matters here, and runs first, but there is no reason to turn the second one off.
        t.extractall(staging, members=_safe_members(t), filter="data")      # noqa: S202 - every member is checked above
    if not (staging / "genomecurator" / "__init__.py").is_file():
        shutil.rmtree(staging, ignore_errors=True)
        raise ValueError("the code bundle holds no package")
    (staging / "genomecurator" / "BUILD").write_text(commit + "\n", encoding="utf-8")
    shutil.rmtree(dest, ignore_errors=True)
    os.replace(staging, dest)
    return dest


def runs_here(code_dir: Path, commit: str, version: str | None = None) -> bool:
    """Does that tree actually run, and does it say it is the commit it was sent as? Code that cannot import itself must never
    be started: the worker would exit on it and come back to the same tree for ever."""
    import json
    import sys
    env = dict(os.environ, PYTHONPATH=str(code_dir) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    try:
        r = subprocess.run([sys.executable, "-c", "import genomecurator as g, json; print(json.dumps(g.build()))"],
                           capture_output=True, text=True, timeout=180, env=env, cwd=str(code_dir))
    except Exception as ex:
        _log.warning("the new code could not be tried: %s", ex)
        return False
    if r.returncode != 0:
        _log.warning("the new code does not import: %s", r.stderr.strip()[:300])
        return False
    try:
        b = json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        _log.warning("the new code did not report what it is")
        return False
    if b.get("commit") != commit or (version and b.get("version") != version):
        _log.warning("the new code reports %s %s, not %s %s", b.get("version"), b.get("commit"), version, commit)
        return False
    return True


def fetch(main, commit: str, sha: str, dest_root: Path, version: str | None = None) -> Path:
    """Download the bundle from the main server, check it against the SHA-256 it declared, unpack it and try it.
    -> the folder to run from. Raises if anything does not add up; the caller then keeps the code it has."""
    dest_root.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=dest_root, suffix=".tar.gz", delete=False) as f:
        tmp = Path(f.name)
    try:
        main.download("/api/cluster/w/code/bundle", tmp)
        got = sha256(tmp)
        if got != sha:
            raise ValueError(f"the code bundle does not match what the main server declared ({got[:12]} vs {sha[:12]})")
        code_dir = unpack(tmp, dest_root / commit, commit)
    finally:
        tmp.unlink(missing_ok=True)
    if not runs_here(code_dir, commit, version):
        import shutil
        shutil.rmtree(code_dir, ignore_errors=True)
        raise ValueError("the code the main server sent does not run here; keeping the code this worker has")
    return code_dir
