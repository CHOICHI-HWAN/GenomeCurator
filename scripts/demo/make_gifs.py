"""Record the README demo GIFs against throw-away copies of demo data (the real data/ is never touched).

    python scripts/demo/make_gifs.py explore            # one scenario
    python scripts/demo/make_gifs.py all

Scenarios that need data expect a prepared folder per scenario under $GC_DEMO_ROOT (see DEMO_DIRS below); `quickstart`
starts from an empty folder and talks to the real NCBI / BV-BRC / Pathoplexus (it needs internet, ~3 minutes).
The analysis scenarios (spread_map, snipit, outbreak, novel_pca) run on a scratch copy of the database in ./data (or $GC_DEMO_REAL_DATA) that has a finished mutation-profile run;
they open only pages that show no genome ids."""
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
        c.caption("<b>스키마 점검</b> — 소스 3곳이 내 테이블의 어떤 열을 채우는지 <b>한 표</b>로 봅니다")
        c.js("document.querySelector('.vtabs').scrollIntoView()")
        c.hold(2800)
        c.js("window.scrollBy(0, 330)")
        c.caption("초록이면 정상 · 색 알약은 그 열을 채우는 <b>소스 필드</b>입니다")
        c.hold(2400)
        c.js("document.querySelector('.vtabs').scrollIntoView()")
        c.click(css=".vt[data-v=board]", after=0.8)
        c.js("document.querySelector('#v-board').scrollIntoView()")
        c.caption("필드 단위로 선을 따라가는 <b>구조도</b>도 있습니다 — 열을 누르면 채우는 필드가 우선순위 순으로 나옵니다")
        c.hold(1800)
        c.click(css=".col[data-col=host_name]", after=0.4)
        c.hold(2800)
        c.click(css=".col[data-col=host_name]", after=0.3, show=False)   # unpin
        c.js("document.querySelector('.vtabs').scrollIntoView()")
        c.click(css=".vt[data-v=simple]", after=0.6)
        flag.write_text("pretend the sources changed")
        c.caption("<b>(시연)</b> 소스가 필드를 바꿨다고 가정 — NCBI: <b>host</b> 삭제, <b>length</b> 이름 변경 / Pathoplexus: 응답 구조 변경")
        c.hold(2800)
        c.js("window.scrollTo(0, 0)")
        c.click(text="스키마 점검 실행", after=0.6)
        c.caption("소스마다 실제 응답 <b>샘플</b>을 읽습니다 (데이터는 건드리지 않음)")
        c.record_until("location.pathname==='/schema' && document.body.innerText.includes('업데이트가 필요합니다')", every=0.9, timeout=120, ms=500, max_frames=6)
        c.wait("document.readyState==='complete'")
        c.settle(0.8)
        c.js("document.querySelector('.vtabs').scrollIntoView()")
        c.caption("<b>⚠ 확인이 필요한 열</b>이 맨 위에 모이고, 소스 카드도 <b style='color:#ff8a80'>빨갛게</b> 바뀝니다")
        c.hold(3400)
        c.js("window.scrollBy(0, 470)")
        c.caption("표에서도 영향받는 열이 <b style='color:#ffd25e'>주황</b>, 사라진 필드가 <b style='color:#ff8a80'>빨강</b> — 다른 소스가 받쳐 주는지도 보입니다")
        c.hold(3600)
        c.js("document.querySelector('.vtabs').scrollIntoView()")
        c.click(css=".vt[data-v=board]", after=0.8)
        c.js("document.querySelector('#v-board').scrollIntoView()")
        c.click(css=".col[data-col=genome_length]", after=0.4)
        c.caption("구조도에서는 같은 문제가 <b>연결선</b>까지 칠해져 어느 필드가 어느 열에 영향을 주는지 따라갈 수 있습니다")
        c.hold(3600)
        finish(c, "03-schema.gif")
    finally:
        c.close()
        srv.stop()


def resolve():
    """Needs the demo data 'resolve': a virus that has been checked (an open changeset) but not applied."""
    d = fresh_copy("resolve")
    srv, c = Server(d, 8024), Chrome(height=760)
    try:
        c.goto(srv.url + "/changeset/1")
        c.js("(()=>{const r=[...document.querySelectorAll('#gb tr')].find(x=>/불일치/.test(x.innerText));r.scrollIntoView({block:'center'})})()")
        c.caption("<b>검토</b> — 소스끼리 값이 다른 항목은 <b>정합성</b> 열에 <b>불일치 · 표기차이</b>로 표시됩니다")
        c.hold(2600)
        c.js("(()=>{const r=[...document.querySelectorAll('#gb tr')].find(x=>/불일치/.test(x.innerText));r.click()})()")
        c.wait("document.querySelector('#rslv .rf')", 15)
        c.settle(0.6)
        c.caption("항목을 누르면 <b>무엇이 다른지</b> 소스별 값으로 보여 주고, <b>★ 추천</b>과 그 이유·신뢰도를 알려 줍니다")
        c.hold(3600)
        c.click(css="#rslv input[name=rf0][value='1']", after=0.3)
        c.caption("마음에 드는 값을 고르고 <b>이대로 선택 저장</b> — 항목에 바로 반영되고 다음 확인에서도 유지됩니다")
        c.hold(1800)
        c.click(css="#rslv .btns button.p", after=1.4)
        c.hold(2400)
        c.caption("행의 정합성이 <b>선택 완료</b>로 바뀌었습니다")
        c.hold(1800)
        c.click(text="Auto-pick", after=0.8)
        c.caption("신뢰도가 높은 추천은 <b>Auto-pick</b>으로 한 번에 — <b>낮음</b>은 직접 확인하도록 남겨 둡니다")
        c.wait("document.readyState==='complete'")
        c.settle(1.6)
        c.hold(3600)
        finish(c, "05-resolve.gif")
    finally:
        c.close()
        srv.stop()


def _analysis_setup() -> tuple[Path, dict, int]:
    """A scratch data folder for the analysis demos: a consistent copy of the database and of the kept Auspice files of the real data folder ($GC_DEMO_REAL_DATA, default ./data), read-only
    access to the Nextstrain checkout and the lineage definitions. The scenes only open pages that show no genome ids (private genomes stay unnamed). The folder is deleted afterwards."""
    import sqlite3
    real = Path(os.environ.get("GC_DEMO_REAL_DATA", ROOT / "data"))
    d = DEMO_ROOT / "_run_analysis"
    shutil.rmtree(d, ignore_errors=True)
    (d / "nextstrain" / "auspice").mkdir(parents=True)
    src, dst = sqlite3.connect(real / "genomecurator.db"), sqlite3.connect(d / "genomecurator.db")
    src.backup(dst)
    rid = dst.execute("SELECT id FROM analysis_run WHERE tool='mutprofile' AND status='done' ORDER BY id DESC LIMIT 1").fetchone()
    src.close(), dst.close()
    if not rid:
        raise SystemExit("the analysis demos need a finished mutation-profile run in the data folder")
    if (real / "genomes").exists():       # a dataset counts only genomes that have a sequence file: link the real folder (read only, a link is removed by rmtree without touching its target)
        os.symlink(real / "genomes", d / "genomes", target_is_directory=True)
    for f in (real / "nextstrain" / "auspice").glob("*.json"):
        shutil.copy2(f, d / "nextstrain" / "auspice" / f.name)
    env = {"GC_NEXTSTRAIN_DIR": str(real / "nextstrain" / "mpox"), "GC_LINEAGE_DIR": str(real / "lineage-designation")}
    return d, env, rid[0]


def _scene(c: Chrome, js_target: str, caption: str, ms: int = 1700, block: str = "start"):
    c.js(f"(()=>{{const e={js_target};if(e)e.scrollIntoView({{block:'{block}'}})}})()")
    c.settle(0.35)
    c.caption(caption)
    c.hold(ms)


def _heading(prefix: str) -> str:
    return f"[...document.querySelectorAll('h3')].find(h=>h.textContent.trim().startsWith('{prefix}'))"


def spread_map():
    d, env, rid = _analysis_setup()
    srv, c = Server(d, 8022, env=env), Chrome(width=1180, height=720)
    try:
        c.goto(srv.url + f"/analysis/{rid}/spread", wait=1.0)
        c.wait("typeof D!=='undefined' && D && document.querySelector('#mapsvg')", 90)
        _scene(c, "document.getElementById('mapcard')", "<b>전파 지도</b> — 나라별 유전체가 <b>보고된 시기대로</b> 나타납니다", 1800)
        for frac in (0.12, 0.25, 0.38, 0.5, 0.62, 0.75, 0.88, 1.0):
            c.js(f"stepTo(Math.min(D.periods.length-1,Math.floor(D.periods.length*{frac})))")
            c.settle(0.45)
            c.caption("<b>" + str(c.js("D.periods[CUR]")) + "</b> — 원 크기 = 유전체 수, 색 = lineage, 분홍 고리 = <b>처음 관찰된 나라</b>")
            c.hold(750)
        _scene(c, "document.getElementById('fig')", "<b>타임라인</b> — 변이가 처음 나타난 시기에 라벨이 <b>겹치지 않게</b> 놓입니다", 2200)
        _scene(c, "document.getElementById('lcard')", "아래 <b>라인차트</b>는 같은 가로축으로 보고 유전체 수를 보여 줍니다", 2000, "center")
        c.js("document.getElementById('v').value='mutation';load()")
        c.wait("typeof D!=='undefined' && D.view==='mutation' && D.rows.length>0", 120)
        c.settle(1.0)
        _scene(c, "document.getElementById('fig')", "<b>대표 변이</b> 보기 — lineage가 나뉘어도 변이 기준 그림은 유지됩니다", 2400)
        finish(c, "06-spread-map.gif")
    finally:
        c.close()
        srv.stop()
        shutil.rmtree(d, ignore_errors=True)


def snipit():
    d, env, rid = _analysis_setup()
    srv, c = Server(d, 8023, env=env), Chrome(width=1180, height=720)
    try:
        c.goto(srv.url + f"/analysis/{rid}/snipit", wait=1.0)
        c.wait("document.querySelector('#fig svg')", 90)
        c.js("document.getElementById('f_ids').checked=false;render()")
        c.settle(0.5)
        c.caption("<b>변이 그림</b> — 참조 위에 서열마다 한 줄, 색 = 바뀐 염기, 분홍 점 = <b>APOBEC3형</b>")
        c.hold(2400)
        _scene(c, "document.getElementById('fig')", "왼쪽은 lineage·clade, 굵은 눈금은 그 lineage의 <b>정의 변이</b>입니다", 2200)
        c.js("document.getElementById('f_gene').value='OPG175';draw()")
        c.wait("typeof D!=='undefined' && D.start>1", 30)
        c.settle(0.6)
        c.js("document.getElementById('f_ids').checked=false;render()")
        c.caption("<b>유전자로 확대</b> — 아래 유전자 트랙에서 위치와 방향이 보입니다")
        c.hold(2400)
        c.js("document.getElementById('f_gene').value='';draw();document.getElementById('f_ap').checked=true;")
        c.wait("typeof D!=='undefined' && D.start===1", 30)
        c.js("draw()")
        c.settle(1.0)
        c.js("document.getElementById('f_ids').checked=false;render()")
        c.caption("<b>APOBEC3형만</b> 보면 편집 변이의 분포가 한눈에 들어옵니다")
        c.hold(2400)
        finish(c, "07-snipit.gif")
    finally:
        c.close()
        srv.stop()
        shutil.rmtree(d, ignore_errors=True)


def outbreak_analysis():
    d, env, rid = _analysis_setup()
    srv, c = Server(d, 8024, env=env), Chrome(width=1180, height=720)
    try:
        c.goto(srv.url + f"/analysis/{rid}/outbreak", wait=1.0)
        c.wait("document.querySelector('#root table')", 120)
        c.caption("<b>유행 분석</b> — 트리에서 <b>유전체 서브클러스터</b>와 가지 정의 변이를 자동으로 찾습니다")
        c.hold(2400)
        _scene(c, "document.querySelector('#root svg')", "<b>클러스터 × 변이</b> 행렬 — 색 = 바뀐 염기", 2200)
        _scene(c, _heading("2."), "<b>APOBEC3형 변이</b>는 시간이 지나며 늘어납니다 (N 많은 서열을 거르고, lineage 안에서도 확인)", 2600)
        _scene(c, _heading("4."), "<b>유전자 hotspot</b> — 길이로 기대되는 것보다 독립 사건이 유의하게 많은 유전자만", 2600)
        c.goto(srv.url + f"/analysis/{rid}/insight", wait=1.0)
        c.wait("document.querySelector('#root table')", 120)
        c.caption("<b>변이 분석</b> — APOBEC3형 변이도 전파되는가?")
        c.hold(2400)
        _scene(c, _heading("4."), "<b>APOBEC3 대 배경 변이</b> — 독립 사건의 약 70%가 APOBEC3형, 우연의 2배 이상 몰려 있습니다", 3000)
        finish(c, "08-outbreak.gif")
    finally:
        c.close()
        srv.stop()
        shutil.rmtree(d, ignore_errors=True)


def novel_pca():
    d, env, rid = _analysis_setup()
    srv, c = Server(d, 8025, env=env), Chrome(width=1180, height=720)
    try:
        c.goto(srv.url + f"/analysis/{rid}/novel", wait=1.0)
        c.wait("document.querySelector('#root table')", 120)
        c.caption("<b>신규 변이 탐색</b> — 정의에 없는 변이를 공유하는 새 계통군 후보")
        c.hold(2400)
        _scene(c, _heading("신규 변이가 계속"), "월별로 <b>새로 보인 변이</b>(막대)와 누적(선) — 붉은 영역은 검출 지연 구간", 2400)
        _scene(c, _heading("후보"), "후보마다 <b>신규 · 성장 중 · 정체</b> 상태와 월별 분포가 붙습니다", 2800)
        c.goto(srv.url + f"/analysis/{rid}/pca", wait=1.0)
        c.wait("document.querySelector('#root svg')", 120)
        c.caption("<b>변이 묶음</b> — 유전체 PCA: 계통군이 덩어리로 갈라집니다")
        c.hold(2600)
        c.js("BY='country';render()")
        c.settle(0.5)
        c.caption("색 기준을 <b>나라</b>로 바꿔 봅니다")
        c.hold(1800)
        _scene(c, _heading("2."), "<b>함께 나타나는 변이</b>는 한 묶음으로 — 신규 후보와도 연결됩니다", 3000)
        finish(c, "09-novel-pca.gif")
    finally:
        c.close()
        srv.stop()
        shutil.rmtree(d, ignore_errors=True)


SCENARIOS = {"quickstart": quickstart, "explore": explore, "integrity": integrity, "schema": schema, "resolve": resolve,
             "spread_map": spread_map, "snipit": snipit, "outbreak": outbreak_analysis, "novel_pca": novel_pca}

if __name__ == "__main__":
    which = sys.argv[1:] or ["all"]
    for name in (SCENARIOS if which == ["all"] else which):
        print("== recording", name, flush=True)
        SCENARIOS[name]()
