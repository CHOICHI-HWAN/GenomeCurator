"""Write docs/DEPENDENCIES.md: every library GenomeCurator uses, with the version it was tested with and its licence, read from
the installed packages (so the list cannot drift from reality). Re-run it when pyproject.toml or superset/requirements.txt change:

    python scripts/dependencies.py --superset-python /path/to/gc-superset/bin/python

The purposes below are written by hand; versions and licences are measured."""
import argparse
import json
import re
import subprocess
import sys
from datetime import date
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# group, package, purpose (한국어 / English)
DIRECT = [
    ("runtime", "fastapi", "웹 서버·API 프레임워크 / Web framework (UI and JSON endpoints)"),
    ("runtime", "uvicorn", "ASGI 서버 / ASGI server that runs the web app"),
    ("runtime", "jinja2", "HTML 템플릿 / HTML templates"),
    ("runtime", "httpx", "NCBI·BV-BRC·Pathoplexus 호출(재시도·시간제한) / HTTP client for all source requests"),
    ("runtime", "apscheduler", "매일 확인·스키마 점검 스케줄러 / Scheduler for the daily checks"),
    ("runtime", "python-multipart", "웹 폼(승인·적용 등) 해석 / Parses web form posts"),
    ("dev", "pytest", "테스트 / Test runner"),
    ("demo", "pillow", "README GIF 조립 / Assembles the README GIFs"),
    ("demo", "websocket-client", "헤드리스 Chrome 조종(DevTools 프로토콜) / Talks to headless Chrome over the DevTools protocol"),
    ("superset", "apache-superset", "다차원 분석 화면 / Multidimensional analysis UI"),
    ("superset", "rich", "Superset CLI가 필요로 하는데 패키지 의존성에서 빠진 것 / Needed by Superset's CLI but missing from its own dependency list"),
    ("superset", "cachetools", "위와 같음 (Superset DB 엔진 모듈) / Same: needed by a Superset DB engine module"),
    ("superset", "flask-caching", "Superset과 맞는 <2.4로 고정 / Pinned <2.4 to match Superset"),
    ("superset", "requests", "superset/bootstrap.py가 Superset REST API를 호출 / Used by superset/bootstrap.py to call Superset's REST API"),
]
GROUPS = {"runtime": "실행에 필요 / Runtime (`pip install -e .`)", "dev": "개발·테스트 / Development (`.[dev]`)",
          "demo": "README GIF 녹화(선택) / Demo recording, optional (`.[demo]`)", "superset": "Superset 분석(선택) / Superset analysis, optional (`superset/environment.yml`)"}
SERVICES = [
    ("NCBI Datasets v2 (`api.ncbi.nlm.nih.gov`)", "accession 집합·메타데이터·서열 / The accession set, metadata and sequences", "https://www.ncbi.nlm.nih.gov/datasets/"),
    ("BV-BRC (`www.bv-brc.org`)", "보강 메타데이터·서열 / Enriching metadata and sequences", "https://www.bv-brc.org/"),
    ("Pathoplexus LAPIS (`lapis.pathoplexus.org`, `pathoplexus.org`)", "보강 메타데이터·변이·서열, organism 목록 / Enriching metadata, mutations, sequences; organism list", "https://pathoplexus.org/"),
    ("ICTV via EBI OLS4 (`www.ebi.ac.uk`)", "워치 등록 전 종명 확인 / Species-name check before adding a watch", "https://ictv.global/"),
    ("GitHub Pages / `git`", "공개 페이지 게시(선택) / Publishing the public page (optional)", "https://pages.github.com/"),
]


def licence(dist) -> str:
    md = dist.metadata
    expr = md.get("License-Expression")
    if expr:
        return expr
    lic = (md.get("License") or "").strip()
    if lic and len(lic) <= 60 and "\n" not in lic:
        return lic
    cls = [c.split("::")[-1].strip() for c in md.get_all("Classifier") or [] if c.startswith("License ::") and "OSI Approved" in c]
    cls = [c for c in cls if c != "OSI Approved"]
    return " / ".join(cls) or (lic.splitlines()[0][:60] if lic else "?")


def info(name: str):
    try:
        d = metadata.distribution(name)
        return d.version, licence(d)
    except metadata.PackageNotFoundError:
        return None


def info_in(python: str, names: list[str]) -> dict:
    code = ("import json,sys;from importlib import metadata;sys.path.insert(0,%r);import dependencies as D;"
            "print(json.dumps({n:D.info(n) for n in %r}))") % (str(Path(__file__).parent), names)
    out = subprocess.run([python, "-c", code], capture_output=True, text=True, check=True).stdout
    return json.loads(out.strip().splitlines()[-1])


def closure(roots: list[str]) -> dict[str, tuple[str, str]]:
    """Runtime dependencies pulled in by the direct runtime libraries (name -> (version, licence))."""
    seen, todo = {}, list(roots)
    while todo:
        n = todo.pop()
        key = re.sub(r"[-_.]+", "-", n).lower()
        if key in seen:
            continue
        i = info(n)
        if not i:
            continue
        seen[key] = i
        for r in metadata.distribution(n).requires or []:
            if "extra ==" in r:
                continue
            todo.append(re.split(r"[ ;<>=!~\[(]", r, maxsplit=1)[0])
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--superset-python", help="python of the gc-superset environment (to read its versions/licences)")
    a = ap.parse_args()
    here = {n: info(n) for g, n, _ in DIRECT if g != "superset"}
    sup = info_in(a.superset_python, [n for g, n, _ in DIRECT if g == "superset"]) if a.superset_python else {}
    pyver = sys.version.split()[0]
    L = ["# 사용한 라이브러리 / Libraries used", "",
         f"_이 문서는 `scripts/dependencies.py`가 설치된 패키지에서 뽑아 만듭니다 (버전·라이선스는 측정값, 역할은 손으로 씀). 생성: {date.today()}, Python {pyver}._  ",
         f"_Generated by `scripts/dependencies.py` from the installed packages (versions and licences are measured, purposes are written by hand). Generated {date.today()}, Python {pyver}._", "",
         "버전은 **이 프로젝트를 테스트한 버전**이고, 허용 범위는 `pyproject.toml`·`superset/requirements.txt`가 정합니다. ",
         "Versions are the ones **this project was tested with**; the allowed ranges live in `pyproject.toml` and `superset/requirements.txt`.", ""]
    for g, title in GROUPS.items():
        L += [f"## {title}", "", "| 라이브러리 / Library | 테스트한 버전 / Tested | 라이선스 / Licence | 역할 / Purpose |", "|---|---|---|---|"]
        for gg, n, purpose in DIRECT:
            if gg != g:
                continue
            i = (sup if g == "superset" else here).get(n)
            L.append(f"| `{n}` | {i[0] if i else '(not measured)'} | {i[1] if i else '?'} | {purpose} |")
        L.append("")
    L += ["## 브라우저에서 불러오는 것 / Loaded in the browser", "",
          "| 라이브러리 / Library | 버전 / Version | 라이선스 / Licence | 역할 / Purpose |", "|---|---|---|---|",
          "| [Chart.js](https://www.chartjs.org/) | 4.4.1 (jsDelivr CDN) | MIT | 대시보드와 공개 페이지의 차트(보는 쪽에 인터넷 필요) / Charts on the dashboard and the public page (viewer needs internet) |", "",
          "## 도구 / Tools (not Python libraries)", "",
          "- **conda / mamba** — 환경 생성 / environment creation (`environment.yml`).",
          "- **Google Chrome (headless)** — README GIF 녹화에만 사용(선택) / only for recording the README GIFs (optional).",
          "- **git** — 공개 페이지 게시(선택) / publishing the public page (optional).", "",
          "## 데이터를 가져오는 외부 서비스 / External services the app contacts", "",
          "| 서비스 / Service | 용도 / Use | 홈페이지 / Home |", "|---|---|---|"]
    L += [f"| {n} | {u} | {h} |" for n, u, h in SERVICES]
    L += ["", "데이터는 각 서비스의 공개 데이터입니다. 재배포·인용 전에 각 서비스의 이용 정책을 확인하세요. ",
          "The data come from public services; check each service's data-use policy before redistributing or citing.", ""]
    tr = closure([n for g, n, _ in DIRECT if g == "runtime"])
    direct_keys = {re.sub(r"[-_.]+", "-", n).lower() for _, n, _ in DIRECT}
    L += ["## 간접 의존성 (실행) / Transitive runtime dependencies", "",
          "<details><summary>실행에 함께 설치되는 패키지 %d개 / %d packages installed with the runtime libraries</summary>" % (len(tr) - len([k for k in tr if k in direct_keys]), len(tr) - len([k for k in tr if k in direct_keys])), "",
          "| 패키지 / Package | 버전 / Version | 라이선스 / Licence |", "|---|---|---|"]
    L += [f"| `{k}` | {v[0]} | {v[1]} |" for k, v in sorted(tr.items()) if k not in direct_keys]
    L += ["", "</details>", "",
          "Superset은 자체 의존성(Flask 등)을 다수 설치합니다. 목록은 [Superset 문서](https://superset.apache.org/)와 `pip show apache-superset`를 보세요. ",
          "Superset installs many dependencies of its own (Flask and others); see its documentation and `pip show apache-superset`.", ""]
    (ROOT / "docs" / "DEPENDENCIES.md").write_text("\n".join(L), encoding="utf-8", newline="\n")
    print("wrote docs/DEPENDENCIES.md")


if __name__ == "__main__":
    main()
