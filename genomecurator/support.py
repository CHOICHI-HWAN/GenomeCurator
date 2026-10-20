"""How far each part of the app is built for each virus - shown as badges so that nobody mistakes "works for monkeypox" for "works for every virus".

Levels: `done` (built and checked on real data), `wip` (works in part or only for some data; read the note), `none` (not built for this virus). The matrix is written down here, in one place, and
is updated by hand when something is finished - it is a statement by the developers, not something the app measures."""
AREAS = [("collect", "수집·큐레이션"), ("vocab", "용어 사전"), ("qc", "서열 QC"), ("analysis", "분석(Nextstrain·Squirrel·변이)")]
LABEL = {"done": "완료", "wip": "개발중", "none": "미지원"}

MPOX = {"collect": ("done", "NCBI Virus · BV-BRC · Pathoplexus에서 가져와 변경을 승인하며 큐레이션합니다."),
        "vocab": ("done", "숙주·나라·수집일·clade·플랫폼·검체 사전이 원숭이두창 데이터로 만들어졌습니다. 모호한 값은 검토 대상으로 남깁니다."),
        "qc": ("done", "서열 통계 + Nextclade(mpox 데이터셋) + Squirrel QC."),
        "analysis": ("done", "Nextstrain 트리, Squirrel 정렬, 변이 프로파일과 결과 화면 전부가 원숭이두창 clade I(Ib)로 검증되었습니다. clade II 빌드는 같은 코드지만 실제 데이터로 충분히 점검하지 못했습니다.")}
OTHER = {"collect": ("done", "수집·비교·승인 흐름은 바이러스에 상관없이 같습니다. 아직 이 바이러스를 큐레이션하지는 않았습니다."),
         "vocab": ("wip", "숙주·나라·수집일·플랫폼·검체 사전은 공통이지만 clade 사전이 없습니다. 이 바이러스의 값으로 점검하지 않았습니다."),
         "qc": ("wip", "서열 통계는 공통입니다. Nextclade는 해당 바이러스의 데이터셋이 있을 때만 돌고, 그 결과를 이 바이러스로 검증하지는 않았습니다."),
         "analysis": ("none", "Nextstrain·Squirrel·변이 프로파일은 원숭이두창 워크플로 기준입니다. 다른 바이러스에는 아직 준비되지 않았습니다.")}


def for_virus(*names: str) -> list[dict]:
    """[{area, label, level, text, note}] for a virus given its name / taxon."""
    t = " ".join(n or "" for n in names).lower()
    m = MPOX if ("monkeypox" in t or "mpox" in t) else OTHER
    return [{"area": a, "label": lbl, "level": m[a][0], "text": LABEL[m[a][0]], "note": m[a][1]} for a, lbl in AREAS]
