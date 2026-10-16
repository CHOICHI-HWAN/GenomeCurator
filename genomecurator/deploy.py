"""One-click publishing of the public read-only page (docs/index.html -> GitHub Pages).

Only docs/index.html is committed and pushed (``git commit -- docs/index.html``), so unrelated work in
the repository is never swept into the publish commit. Credentials are whatever git already has
(Git Credential Manager); prompts are disabled so a missing login fails fast instead of hanging."""
import os
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

from .static_site import build_html, snapshot

ROOT = Path(__file__).resolve().parent.parent
REL = "docs/index.html"
_GEN = re.compile(r'"generated":"\d{4}-\d{2}-\d{2}"')


def _git() -> str:
    exe = shutil.which("git") or next((p for p in (r"C:\Program Files\Git\cmd\git.exe",
                                                    r"C:\Program Files (x86)\Git\cmd\git.exe") if os.path.exists(p)), None)
    if not exe:
        raise RuntimeError("git을 찾을 수 없습니다. git을 설치하고 PATH에 추가하세요.")
    return exe


def _run(root: Path, *args: str, timeout: int = 180) -> str:
    r = subprocess.run([_git(), *args], cwd=root, capture_output=True, text=True, encoding="utf-8", errors="replace",
                       timeout=timeout, env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})
    if r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args[:2])} 실패: {(r.stderr or r.stdout).strip()[:400]}")
    return r.stdout.strip()


def _check_repo(root: Path) -> str:
    if not (root / ".git").exists():
        raise RuntimeError("이 폴더는 git 저장소가 아닙니다. 저장소 폴더에서 앱을 실행하세요.")
    try:
        return _run(root, "remote", "get-url", "origin")
    except RuntimeError:
        raise RuntimeError("git remote 'origin'이 설정되어 있지 않습니다.") from None


def pages_url(remote: str) -> str | None:
    m = re.search(r"github\.com[:/]+([^/]+)/([^/]+?)(?:\.git)?/?$", remote)
    return f"https://{m.group(1).lower()}.github.io/{m.group(2)}/" if m else None


def _repo_web(remote: str) -> str:
    m = re.search(r"github\.com[:/]+([^/]+)/([^/]+?)(?:\.git)?/?$", remote)
    return f"https://github.com/{m.group(1)}/{m.group(2)}" if m else ""


def _norm(html: str) -> str:
    return _GEN.sub('"generated":""', html)


def preview(root: Path = ROOT) -> dict:
    """What would be published, and does it differ from what is already in docs/index.html?"""
    remote = _check_repo(root)
    snap = snapshot(_repo_web(remote))
    if not snap["viruses"]:
        raise RuntimeError("게시할 유전체가 없습니다.")
    html = build_html(_repo_web(remote))
    target = root / REL
    old = target.read_text(encoding="utf-8") if target.exists() else ""
    return {"viruses": len(snap["viruses"]), "genomes": sum(len(v["rows"]) for v in snap["viruses"]),
            "changed": _norm(html) != _norm(old), "pages_url": pages_url(remote),
            "names": [f"{v['name']} {len(v['rows'])}건" for v in snap["viruses"]]}


def publish(progress=None, root: Path = ROOT) -> str:
    p = progress or (lambda *a: None)
    p(5, "저장소 확인 중…")
    remote = _check_repo(root)
    p(15, "공개 화면 생성 중…")
    html = build_html(_repo_web(remote))
    target = root / REL
    old = target.read_text(encoding="utf-8") if target.exists() else ""
    if _norm(html) == _norm(old):
        p(100, "변경 없음")
        return "이미 최신 상태라 게시할 변경이 없습니다."
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(html, encoding="utf-8")
    (target.parent / ".nojekyll").touch()
    p(40, "변경 사항 커밋 중…")
    _run(root, "add", "--", REL, "docs/.nojekyll")
    _run(root, "commit", "-m", f"Update public snapshot ({datetime.now():%Y-%m-%d %H:%M})", "--", REL, "docs/.nojekyll")
    p(65, "GitHub로 업로드 중…")
    try:
        _run(root, "push", "origin", "HEAD")
    except RuntimeError as e:
        raise RuntimeError(f"{e} — 원격에 새 커밋이 있으면 먼저 `git pull --rebase` 후 다시 시도하세요. "
                           "(생성된 파일은 로컬에 커밋되어 있습니다)") from None
    p(100, "업로드 완료")
    url = pages_url(remote)
    return "배포 완료. GitHub Pages 반영까지 1~2분 걸립니다." + (f" {url}" if url else "")
