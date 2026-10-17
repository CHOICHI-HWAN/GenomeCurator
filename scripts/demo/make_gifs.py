"""Record the README demo GIFs against throw-away copies of demo data (the real data/ is never touched).

    python scripts/demo/make_gifs.py explore            # one scenario
    python scripts/demo/make_gifs.py all

Scenarios that need data expect a prepared folder per scenario under $GC_DEMO_ROOT (see DEMO_DIRS below); `quickstart`
starts from an empty folder and talks to the real NCBI / BV-BRC / Pathoplexus (it needs internet, ~3 minutes)."""
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from recorder import Chrome  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "demo"
DEMO_ROOT = Path(os.environ.get("GC_DEMO_ROOT", Path.home() / ".cache" / "genomecurator-demo"))


class Server:
    """A GenomeCurator server on its own port and data folder."""

    def __init__(self, data_dir: Path, port: int, script: list[str] | None = None, env: dict | None = None):
        self.port = port
        cmd = script or [sys.executable, "-m", "genomecurator", "serve", "--port", str(port), "--no-scheduler"]
        self.proc = subprocess.Popen(cmd, cwd=ROOT, env={**os.environ, "GC_DATA_DIR": str(data_dir), "GC_LOG_CONSOLE_LEVEL": "CRITICAL", **(env or {})},
                                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(60):
            try:
                urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=2)
                return
            except Exception:
                time.sleep(0.5)
        raise SystemExit(f"server on {port} did not start")

    @property
    def url(self):
        return f"http://127.0.0.1:{self.port}"

    def stop(self):
        self.proc.terminate()


def fresh_copy(name: str) -> Path:
    """Copy a prepared demo folder to a scratch place so a recording can change it freely."""
    src, dst = DEMO_ROOT / name, DEMO_ROOT / f"_run_{name}"
    if not src.exists():
        raise SystemExit(f"demo data '{name}' not found under {DEMO_ROOT}")
    shutil.rmtree(dst, ignore_errors=True)
    shutil.copytree(src, dst)
    return dst


def finish(c: Chrome, name: str):
    OUT.mkdir(parents=True, exist_ok=True)
    n, size = c.save_gif(OUT / name)
    print(f"{name}: {n} frames, {size / 1024:.0f} KB")


# ---------------------------------------------------------------------------------------------------- scenarios
def quickstart():
    """From an EMPTY folder: register Henipavirus nipahense, check the sources, approve, apply. Real internet, ~3 minutes.
    The folder it leaves behind (_run_quickstart) is what the other scenarios show."""
    d = DEMO_ROOT / "_run_quickstart"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    srv, c = Server(d, 8020), Chrome()
    try:
        c.goto(srv.url + "/")
        c.caption("<b>① 지켜볼 바이러스</b>를 등록합니다 — 이름만 넣으면 됩니다")
        c.hold(1800)
        c.click("#lkq", after=0.2)
        c.type("#lkq", "Henipavirus nipahense", step=3, ms=150)
        c.hold(700)
        c.click("#lkb", after=0.4)
        c.caption("ICTV · NCBI · BV-BRC · Pathoplexus에서 <b>먼저 확인</b>합니다")
        c.record_until("document.querySelector('#lkres table')", every=5, timeout=240, ms=600, max_frames=5)
        c.caption("소스마다 <b>확인된 결과</b>가 나옵니다 (바이러스가 아니면 정중히 거절)")
        c.hold(2800)
        c.click(text="입력한 이름 그대로 추가", after=1.2)
        c.wait("document.querySelector('table')")
        c.caption("등록 완료 — 이제 <b>[지금 확인]</b>으로 소스 3곳을 살펴봅니다")
        c.hold(1800)
        c.click(text="지금 확인", after=0.8)
        c.caption("진행 상황이 <b>로그</b>와 함께 실시간으로 보입니다")
        c.record_until("location.pathname.indexOf('/changeset')===0", every=1.4, timeout=240, ms=550, max_frames=14)
        c.wait("document.readyState==='complete'")
        c.settle(1.0)
        c.caption("<b>② 변경 사항 442건</b> — 아직 아무것도 내 데이터가 아닙니다")
        c.hold(2400)
        c.click(text="Select all", after=0.6)
        c.click(text="Approve", after=1.0)
        c.wait("document.readyState==='complete'")
        c.settle(0.8)
        c.caption("<b>승인한 것만</b> 반영됩니다 (거절·보류도 가능)")
        c.hold(2000)
        c.click(text="Apply", after=0.8)
        c.caption("서열을 받아 <b>검증하고</b> 표에 넣습니다")
        c.record_until("location.pathname.indexOf('/history')===0", every=1.2, timeout=240, ms=550, max_frames=10)
        c.wait("document.readyState==='complete'")
        c.goto(srv.url + "/genomes")
        c.caption("<b>③ 반영 완료</b> — 유전체 442건이 표에 들어왔습니다")
        c.hold(3000)
        finish(c, "01-quickstart.gif")
    finally:
        c.close()
        srv.stop()


def explore():
    d = DEMO_ROOT / "_run_quickstart"
    if not (d / "genomecurator.db").exists():
        raise SystemExit("run the 'quickstart' scenario first (the explore demo shows its result)")
    d2 = DEMO_ROOT / "_run_explore"
    shutil.rmtree(d2, ignore_errors=True)
    shutil.copytree(d, d2)
    srv, c = Server(d2, 8021), Chrome()
    try:
        c.goto(srv.url + "/genomes")
        c.caption("<b>유전체 표</b> — 엑셀처럼 키워드·필터·정렬")
        c.hold(1600)
        c.click("#gq", after=0.2)
        c.type("#gq", "Bangladesh", step=2, ms=150)
        c.settle(0.5)
        c.caption("키워드를 치면 <b>바로 걸러집니다</b>")
        c.hold(1500)
        c.js("(()=>{const e=document.querySelector('#gq');e.value='';e.dispatchEvent(new Event('input',{bubbles:true}))})()")
        c.click(css="#fp .fv", text="Pteropus", after=0.6)
        c.caption("왼쪽 <b>필터</b>는 값을 누르면 겹쳐서 적용됩니다")
        c.hold(1600)
        c.click(css="#gb tr", nth=0, after=0.6)
        c.caption("행을 누르면 오른쪽 <b>상세 패널</b>에 값과 <b>출처</b>가 나옵니다")
        c.hold(2000)
        pos = c.locate("#rz")
        c.drag(pos["x"], pos["y"], pos["x"] - 260, pos["y"])
        c.caption("패널 너비는 <b>끌어서 조절</b>할 수 있습니다")
        c.hold(2200)
        finish(c, "04-explore.gif")
    finally:
        c.close()
        srv.stop()


def integrity():
    d = fresh_copy("mpox")
    srv, c = Server(d, 8022), Chrome()
    try:
        c.goto(srv.url + "/integrity")
        c.caption("<b>정합성</b> — 소스 3곳이 같은 유전체를 서로 얼마나 다르게 말하나")
        c.hold(2600)
        c.js("document.querySelectorAll('h3')[1].scrollIntoView()")
        c.caption("필드마다 <b style='color:#6fd48c'>일치</b> · <b style='color:#ffd25e'>표기만 다름</b> · <b style='color:#ff8a80'>진짜 다름</b>")
        c.hold(2800)
        c.click("td.cell a", nth=1, after=1.0)
        c.wait("document.readyState==='complete'")
        c.js("document.getElementById('ex').scrollIntoView()")
        c.caption("칸을 누르면 <b>실제로 다른 값</b>이 나열됩니다 (예: 수집일 <b>2202-07</b> ↔ 2022-07)")
        c.hold(3200)
        c.js("document.getElementById('seq').scrollIntoView()")
        c.caption("애매하면 <b>서열 본문</b>으로 확인 — 먼저 <b>dry-run</b>으로 계획을 봅니다")
        c.hold(2600)
        c.click(text="계획 확인", after=0.6)
        c.record_until("location.hash==='#seq' && document.body.innerText.includes('실행 계획')", every=1.0, timeout=120, ms=500, max_frames=8)
        c.wait("document.readyState==='complete'")
        c.js("document.querySelector('#seq').scrollIntoView()")
        c.settle(0.6)
        c.caption("받을 서열은 <b>1~2MB</b>, 길이 비교로 <b>수백 MB를 생략</b> — 확인하고 실행합니다")
        c.hold(3400)
        finish(c, "02-integrity.gif")
    finally:
        c.close()
        srv.stop()


def schema():
    d = fresh_copy("schema")
    flag = d / "drift.flag"
    srv = Server(d, 8023, script=[sys.executable, str(Path(__file__).parent / "demo_server.py"), "8023"], env={"GC_DEMO_FLAG": str(flag)})
    c = Chrome()
    try:
        c.goto(srv.url + "/schema")
        c.caption("<b>스키마 점검</b> — NCBI·BV-BRC·Pathoplexus의 구조가 바뀌지 않았나 감시합니다")
        c.hold(2200)
        c.click(text="스키마 점검 실행", after=0.6)
        c.caption("소스마다 실제 응답 <b>샘플</b>을 읽습니다 (데이터는 건드리지 않음)")
        c.record_until("location.pathname==='/schema' && document.body.innerText.includes('필요하지 않습니다')", every=0.9, timeout=120, ms=500, max_frames=8)
        c.wait("document.readyState==='complete'")
        c.settle(0.8)
        c.caption("지금은 <b>이상 없음</b> — 코드는 그대로 써도 됩니다")
        c.hold(2600)
        flag.write_text("pretend the sources changed")
        c.caption("<b>(시연)</b> 소스가 필드 이름을 바꾸고 응답 구조도 바꿨다고 가정해 봅니다")
        c.hold(2200)
        c.js("document.getElementById('daily').scrollIntoView()")
        c.click(text="지금 바이러스별 점검 실행", after=0.6)
        c.record_until("location.pathname==='/schema' && document.body.innerText.includes('확인 필요')", every=0.9, timeout=120, ms=500, max_frames=8)
        c.wait("document.readyState==='complete'")
        c.js("document.getElementById('daily').scrollIntoView()")
        c.settle(0.8)
        c.caption("<b>⚠ 확인 필요</b> — 어느 함수를 고쳐야 하는지까지 알려 줍니다. 상단 탭에는 <b>!</b> 배지가 켜집니다")
        c.hold(3600)
        c.js("window.scrollBy(0, 250)")
        c.hold(2600)
        finish(c, "03-schema.gif")
    finally:
        c.close()
        srv.stop()


SCENARIOS = {"quickstart": quickstart, "explore": explore, "integrity": integrity, "schema": schema}

if __name__ == "__main__":
    which = sys.argv[1:] or ["all"]
    for name in (SCENARIOS if which == ["all"] else which):
        print("== recording", name, flush=True)
        SCENARIOS[name]()
