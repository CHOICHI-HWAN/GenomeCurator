# GenomeCurator

로컬에서 실행하는 바이러스 유전체 수집·큐레이션 도구. NCBI Virus(서열·accession 기준)를 받고,
Pathoplexus > NCBI Virus > BV-BRC 우선순위로 메타데이터를 병합해 BV-BRC 스타일 테이블로 저장합니다.
새 데이터는 **changeset**으로만 쌓이고, 웹 페이지에서 승인해야 로컬에 반영됩니다.

## 실행 (Windows · macOS · Linux)

Python 3.10 이상과 git이 필요합니다. 저장소를 clone한 폴더에서 실행하세요.

| OS | 한 번에 시작 (가상환경 생성 → 설치 → 실행 → 브라우저 열기) |
|---|---|
| Windows | `scripts\start.bat` 더블클릭 (또는 터미널에서 실행) |
| macOS | Finder에서 `scripts/start.command` 더블클릭 (처음 한 번 `chmod +x scripts/start.*` 필요) |
| Linux | `./scripts/start.sh` |

직접 설치하려면 (모든 OS 공통):
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate     macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
python -m genomecurator watch-add "Henipavirus nipahense"
python -m genomecurator serve --open          # http://127.0.0.1:8000
python -m pytest                              # 테스트
```
- 웹(큐레이션 페이지)에서 항목 승인/거절/보류 → Apply.
- 데이터 위치: 기본 `<저장소>/data` (어느 폴더에서 실행해도 같은 곳). 바꾸려면 `GC_DATA_DIR` (예: Windows `set GC_DATA_DIR=D:\gc`, macOS/Linux `export GC_DATA_DIR=~/gc`).
- 환경변수: `NCBI_API_KEY`, `GC_CHECK_INTERVAL_HOURS`(기본 24), `GC_DRYRUN_SAMPLE`(기본 3), `GC_REQUEST_DELAY_MS`(기본 300).
- 공개 페이지 "배포하기"는 git이 PATH 또는 기본 설치 위치에 있고 `origin`에 push 권한(자격 증명)이 있어야 합니다.
- 파일은 모든 OS에서 UTF-8 / LF로 저장됩니다. FASTA도 동일합니다. CI(GitHub Actions)가 Windows·macOS·Linux × Python 3.10/3.12에서 테스트를 돌립니다.

## 변경 판정
메타데이터 값 변경, 또는 GenBank accession `.version` 증가 → changed. 신규 → new, NCBI에서 사라짐 → removed.
거절한 항목은 같은 내용이면 다시 올라오지 않습니다.

## Pathoplexus
`genomecurator/config.py`의 `PATHOPLEXUS_ORGANISMS`에 `{워치 이름: LAPIS organism slug}`를 넣으면 활성화됩니다.
(Nipah/Henipavirus는 현재 Pathoplexus에 organism이 없어 건너뜁니다. 필드 매핑은 미검증.)

## 공개용 정적 페이지 (읽기 전용)

`docs/index.html` 한 파일로 Overview(차트)와 Genomes(필터·정렬·상세·TSV)를 볼 수 있는 스냅샷이다. 수집·검토·반영(다운로드) 코드와 서열, 파일 경로, 이력은 들어 있지 않고 공개 메타데이터만 담긴다.

```bash
python scripts/build_static.py      # 로컬 DB에서 docs/index.html 재생성
```

로컬 앱을 실행 중이면 http://127.0.0.1:8000/public (상단 "공개 보기")에서 같은 화면을 현재 DB 기준으로 바로 볼 수 있다. `docs/index.html`은 파일을 더블클릭해서 열어도 된다.

GitHub Pages로 공개하려면 저장소 Settings → Pages → Source를 `main` 브랜치의 `/docs`로 지정한다. 차트는 CDN(Chart.js)을 쓰므로 보는 쪽에 인터넷이 필요하다.

## 백업 / 복원 (PC 교체·OS 변경)
`data/`(수집 DB·FASTA)는 용량이 커서 git에 넣지 않습니다. 옮기거나 보관할 때는 압축 백업을 만드세요.

```bash
python -m genomecurator backup my-backup.zip      # 큐레이션 결과·결정 이력·FASTA만 담는 압축 파일 (원본 응답·대기 중 항목은 제외)
python -m genomecurator restore my-backup.zip     # 새 환경의 data 폴더로 복원 (기존 DB가 있으면 --force 필요)
```
백업은 `/` 경로와 UTF-8만 써서 Windows·macOS·Linux 어디서나 복원됩니다. 제외된 원본 응답·대기 중 항목은 다음 "지금 확인"에서 다시 만들어집니다.
