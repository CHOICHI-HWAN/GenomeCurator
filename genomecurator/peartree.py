"""PearTree (Andrew Rambaut, MIT; https://github.com/artic-network/peartree) as the tree viewer for BEAST summary trees.

Its code is fetched into the data folder on request and served from there; it is not part of this repository (its `pearcore` library carries no licence file, so it is not
redistributed). Nothing is built: the browser app is plain static files. Everything runs in the page, the tree is not sent anywhere."""
from __future__ import annotations

import json
import re
import shutil
import subprocess

from . import config

APP_REPO, APP_COMMIT = "https://github.com/artic-network/peartree", "7a6ce610e35e87eff4ebbf9c70285c3b8b8769ed"
CORE_REPO, CORE_COMMIT = "https://github.com/rambaut/pearcore-lib", "4abb653d654b38fca470b602bf160f4265655e2e"
KEEP_OUT = ("manual", "data", "src-tauri", "node_modules", "docs")      # not needed to view a tree; the example data and the manual are tens of MB


def root():
    return config.DATA_DIR / "tools" / "peartree"


def installed() -> bool:
    r = root()
    return (r / "peartree" / "peartree.html").exists() and (r / "pearcore" / "js" / "pearcore-ui.js").exists()


def status() -> dict:
    meta = root() / "INSTALLED.json"
    return {"installed": installed(), "app_commit": APP_COMMIT, "core_commit": CORE_COMMIT,
            "url": "/tools/peartree/peartree/peartree.html" if installed() else None,
            "meta": json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else None}


def _git(*args, cwd=None):
    r = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, timeout=900)
    if r.returncode != 0:
        raise RuntimeError("git " + " ".join(args[:2]) + ": " + (r.stderr or r.stdout)[-300:])


def install(progress=None, force: bool = False) -> dict:
    """Fetch PearTree and pearcore at the pinned commits into `<data>/tools/peartree/` (about 5 MB after removing the manual and example data)."""
    p = progress or (lambda *a: None)
    if installed() and not force:
        return status()
    if not shutil.which("git"):
        raise RuntimeError("git이 필요합니다.")
    r = root()
    if r.exists():
        shutil.rmtree(r)
    r.parent.mkdir(parents=True, exist_ok=True)
    for repo, commit, dest, pct, label in ((APP_REPO, APP_COMMIT, r, 10, "PearTree 받는 중…"), (CORE_REPO, CORE_COMMIT, r / "pearcore", 60, "pearcore 받는 중…")):
        p(pct, label)
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        _git("init", "-q", cwd=dest)
        _git("remote", "add", "origin", repo, cwd=dest)
        _git("fetch", "-q", "--depth", "1", "origin", commit, cwd=dest)
        _git("checkout", "-q", "FETCH_HEAD", cwd=dest)
        shutil.rmtree(dest / ".git", ignore_errors=True)
    for name in KEEP_OUT:
        for d in (r / name, r / "peartree" / name):
            shutil.rmtree(d, ignore_errors=True)
    if not installed():
        raise RuntimeError("PearTree 파일이 예상한 위치에 없습니다 (저장소 구조가 바뀌었을 수 있습니다).")
    (r / "INSTALLED.json").write_text(json.dumps({"app": APP_COMMIT, "core": CORE_COMMIT, "licence": "PearTree: MIT (Andrew Rambaut); pearcore: no licence file in its repository"}), encoding="utf-8")
    p(100, "완료")
    return status()


def annotate(nexus: str, tips: dict[str, dict]) -> str:
    """Put each tip's annotations (country, lineage, clade, date ...) into the tip's own `[&...]` of a TreeAnnotator tree, so the viewer can colour by them.
    tips: {tip name: {annotation: value}}. Tips are the numbers of the Translate block; internal nodes follow ')' and are left alone."""
    m = re.search(r"Translate\s*(.*?);\s*tree\s", nexus, re.S | re.I)
    if not m:
        return nexus
    num = {}
    for line in m.group(1).split(","):
        parts = line.split(None, 1)
        if len(parts) == 2:
            num[parts[0].strip()] = parts[1].strip().strip("'\"")
    def piece(k, v):
        try:
            return f"{k}={float(v)}" if k == "date" else f'{k}="{str(v).replace(chr(34), "")}"'
        except (TypeError, ValueError):
            return f'{k}="{str(v).replace(chr(34), "")}"'
    def sub(mm):
        name = num.get(mm.group(2))
        a = tips.get(name or "")
        if not a:
            return mm.group(0)
        return f"{mm.group(1)}{mm.group(2)}[&" + ",".join(piece(k, v) for k, v in a.items() if v not in (None, "")) + ","
    head, sep, tree = nexus.partition("\ntree ")
    return head + sep + re.sub(r"([(,])(\d+)\[&", sub, tree)
