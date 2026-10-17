<div align="center">

```
 ╔══════════════════════════════════════════════════════════════╗
 ║                                                              ║
 ║        V I R U S   G E N O M E   C U R A T O R               ║
 ║            ~ Metadata Curator for virus genomes ~            ║
 ║                                                              ║
 ╚══════════════════════════════════════════════════════════════╝
        바이러스 유전체를 "승인받고" 들이는 깐깐한 사서
```

[![tests](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml/badge.svg)](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.10%20%7C%203.12-blue)
![version](https://img.shields.io/github/v/tag/CHOICHI-HWAN/GenomeCurator?label=version&sort=semver)
![license](https://img.shields.io/badge/license-MIT-green)
![os](https://img.shields.io/badge/Windows%20%7C%20macOS%20%7C%20Linux-ok-lightgrey)

</div>

> 데이터베이스 세 곳이 같은 바이러스를 서로 다르게 말할 때, 누구 말을 믿을까요?
> Virus Genome Curator는 **셋 다 듣고, 서로 얼마나 다른지 숫자로 보여 주고, 당신이 도장을 찍은 것만** 내 컴퓨터에 들입니다.

> **범위**: 이 도구는 **바이러스 전용**입니다 (NCBI Virus · BV-BRC · Pathoplexus). 이름의 "Genome"은 *바이러스 유전체*를 뜻합니다.
> 세균은 레코드 단위(어셈블리, contig 여러 개)와 서열 크기(수 MB)가 달라서 이 도구에 억지로 넣지 않고 **별도 프로젝트**로 만들 예정입니다.
> 저장소·패키지·명령어 이름은 그대로 `GenomeCurator` / `genomecurator`입니다.

## 한 줄 요약 (Synopsis)

```bash
mamba env create -f environment.yml && mamba activate genomecurator

python -m genomecurator watch-add "Henipavirus nipahense"   # 이 바이러스를 지켜보세요
python -m genomecurator serve --open                         # 그리고 웹에서 승인하세요
```

그게 전부입니다. 나머지는 이 문서가 아니라 웹 화면이 알려 줍니다.

<div align="center">

![빈 화면에서 바이러스를 등록하고, 확인하고, 승인해서 반영하기까지](docs/demo/01-quickstart.gif)

<sub>▶ 빈 화면에서 Henipavirus nipahense를 등록 → 소스 3곳 확인 → 442건 승인 → 반영까지 (실제 NCBI·BV-BRC·Pathoplexus에 접속한 화면, 약 2분을 30초로 줄였습니다)</sub>

</div>

## 왜 만들었나 (그리고 왜 이렇게 깐깐한가)

유전체 메타데이터를 NCBI에서 받고, BV-BRC에서 받고, Pathoplexus에서도 받다 보면 이런 일이 생깁니다.

* 같은 accession인데 **길이가 197,179 bp**라는 곳과 **197,157 bp**라는 곳이 있다.
* 수집일이 `2022-07`이라는 곳, `2022-07-11`이라는 곳, 그리고 **`1922-05-30`** 이라는 (아마도 오타인) 곳이 있다.
* 어제 있던 필드가 오늘 사라졌는데, 아무도 알려 주지 않는다.
* 스크립트가 조용히 실패하고, 빈 칸으로 덮어쓴 걸 석 달 뒤에 발견한다.

이 도구는 이 사고들을 **한 번도 겪지 않기 위해** 만들어졌습니다. 그래서 의심이 많습니다.

| 철칙 | 뜻 |
|---|---|
| **승인 전에는 아무것도 안 바뀐다** | 새 소식은 전부 *changeset*(변경 사항 묶음)으로만 쌓입니다. 승인하고 `Apply`를 눌러야 내 데이터가 됩니다. |
| **거절한 건 다시 안 온다** | 같은 내용이면 다음 확인 때 또 올라오지 않습니다. |
| **소스가 하나라도 아프면 멈춘다** | 일부 소스만으로 병합하면 값이 비워질 수 있어서, 조용히 넘어가지 않고 확인 자체를 중단합니다. |
| **말이 다르면 서열로 물어본다** | 메타데이터가 애매하면 서열 본문(헤더 뺀)이 같은지로 *같은 레코드인지* 확인합니다. |
| **소스도 바뀐다** | NCBI·BV-BRC·Pathoplexus의 응답 구조가 바뀌었는지 매일 감시합니다. |

## 소스 삼형제

```
                ┌────────────┐
   NCBI Virus ──┤  어느 accession이 있는지 "정함" (기준)
                └─────┬──────┘
                      │  같은 accession끼리 묶어서
   Pathoplexus ───────┤  메타데이터를 채움 (1순위)
   BV-BRC ────────────┘  그다음 NCBI, 마지막 BV-BRC
                      ▼
            필드마다 "누구 값을 썼는지" 기록 (field_provenance)
```

* **NCBI Virus**가 accession 집합을 정합니다. 나머지는 그 목록에 *보강*할 뿐, 목록을 늘리지 않습니다.
* 값이 겹치면 **Pathoplexus > NCBI Virus > BV-BRC** 순으로 채택합니다.
* 내부 표는 BV-BRC 스타일 컬럼(`genome`, `genome_sequence`)입니다.

## 설치 (Installation)

### 추천: conda / mamba

```bash
git clone https://github.com/CHOICHI-HWAN/GenomeCurator.git
cd GenomeCurator
mamba env create -f environment.yml      # conda 쓰시면 mamba 대신 conda
mamba activate genomecurator
python -m pytest -q                      # 선택: 전부 초록불이면 정상
```

### 한 번에 시작 (환경 생성 → 실행 → 브라우저 열기)

| OS | 실행 |
|---|---|
| Linux | `./scripts/start.sh` |
| macOS | Finder에서 `scripts/start.command` 더블클릭 (처음 한 번 `chmod +x scripts/start.*`) |
| Windows | `scripts\start.bat` 더블클릭 |

스크립트는 `mamba`/`conda`가 있으면 `genomecurator` 환경을 (없으면 만들어서) 쓰고, 둘 다 없으면 `.venv`로 대신합니다.
포트를 바꾸려면 `./scripts/start.sh --port 8001`.

### conda 없이

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

Python 3.10 이상, 인터넷(NCBI·BV-BRC·Pathoplexus 접속)이 필요합니다. "공개 페이지 배포"를 쓰려면 git과 `origin` push 권한도 필요합니다.

## 빠른 시작 (Quick start)

```bash
# 1) 지켜볼 바이러스 등록 (이름을 NCBI에서 확인하고, 바이러스가 아니면 정중히 거절합니다)
python -m genomecurator watch-add "Henipavirus nipahense"

# 2) 웹 열기  →  http://127.0.0.1:8000
python -m genomecurator serve --open

# 3) 큐레이션 화면에서 [지금 확인] → [검토하기] → 항목 선택 → Approve → Apply
```

도중에 길을 잃으면 각 화면 상단의 파란 안내 상자를 읽어 보세요. 사람 말로 적어 두었습니다.

> 포트 8000이 이미 쓰이고 있다면(`Address already in use`) `--port 8001`처럼 바꾸세요.

## 화면 구경 (The tour)

| 탭 | 하는 일 |
|---|---|
| **큐레이션** | 바이러스 목록(추가·편집·삭제) + *검토할 변경 사항* 카드. "새로 발견 / 정보 바뀜 / 사라짐" 숫자, 승인·거절·보류 진행 막대, 다음 할 일을 한 문장으로 알려 줍니다. |
| **대시보드** | 월별 등록, 수집 연도, 국가, 숙주 같은 차트. |
| **유전체** | 엑셀처럼 필터·정렬·열 선택·TSV 내보내기. 오른쪽 상세 패널은 **끌어서 너비를 조절**할 수 있고, Pathoplexus 변이 정보도 불러옵니다. |
| **정합성** | 세 소스가 서로 얼마나 일치하는지 한눈에. 일치(초록) / 표기만 다름(노랑) / 진짜 다름(빨강). 서열로 검증하는 버튼도 여기. |
| **스키마** | 소스와 내 DB의 *구조*가 코드의 기대와 맞는지 점검. **간단히 보기**(내 열 × 소스 3곳 표, 확인이 필요한 열이 맨 위에 모임)와 **자세히 보기**(소스 필드가 내 `genome` 열에 이어지는 선까지 보이는 구조도)를 탭으로 오가며, 이상이 있는 필드·열을 색으로 강조합니다. 문제가 있으면 탭에 `!` 배지가 뜹니다. |
| **작업** | 서버에서 도는 모든 작업(추가·조회·확인·적용·다운로드·검증·점검·배포·자동 실행)의 상태·진행률·로그·오류. 다른 메뉴로 가도 계속되고 여기서 다시 봅니다. |
| **이력** | 누가 언제 무엇을 했는지 + 서버 로그(경고만 / 오류만 필터). |
| **공개 보기 / 배포하기** | 읽기 전용 정적 페이지를 GitHub Pages로 게시. |

오래 걸리는 작업은 화면 가운데에 진행 모달(큰 진행 막대 + 로그)이 뜨고, 새로고침해도 이어집니다.

<div align="center">

![유전체 표: 키워드 검색, 필터, 상세 패널, 패널 너비 조절](docs/demo/04-explore.gif)

<sub>▶ 유전체 표 — 키워드·필터·상세 패널(값마다 <i>어느 소스에서 왔는지</i> 표시)·끌어서 너비 조절</sub>

</div>

## 승인 흐름 (changeset)

```
 [지금 확인] ──► 소스 3곳 조회 ──► 내 데이터와 비교 ──► changeset 생성 (아직 아무것도 안 바뀜)
                                                         │
                          ┌──────────────────────────────┘
                          ▼
              검토: Approve / Reject / Hold  (정합성 열로 "소스끼리 다른 항목"만 골라 볼 수 있음)
                          │
                          ▼
            적용 전 점검(Dry-run): 디스크 · 링크 · 샘플 서열 확인
                          │
                          ▼
                [Apply]  ──►  표와 FASTA에 반영 (중간에 예외가 나면 파일·DB 모두 되돌림)
```

판정 규칙: 메타데이터 값이 바뀌거나 GenBank `.version`이 올라가면 **changed**, 처음 보면 **new**, NCBI에서 사라지면 **removed**.

## 정합성 — 이 도구의 집착

### 1. 소스끼리 맞는지 (reconcile)

`check`할 때마다 같은 accession을 세 소스에서 필드별로 비교합니다. 결과는 세 부류입니다.

* **일치**: 그냥 같음.
* **표기 차이**: 같은 사실을 다르게 적음 (`USA: CA` vs `USA, CA`, `2022-07` vs `2022-07-11`, `Zaire` vs `Democratic Republic of the Congo`). 불일치로 세지 않습니다.
* **불일치**: 진짜 다름. 사람이 봐야 합니다.

mpox 1만 3천 건에서는 소스 간 일치율 99%대였고, 남은 불일치 중에는 오타(`1922-05-30`)와 국가 충돌(콩고민주공화국 ↔ 대한민국)이 있었습니다. 이런 걸 잡으려고 만들었습니다.

<div align="center">

![정합성 탭: 소스별 일치 막대, 실제로 다른 값 목록, 서열 검증 dry-run 계획](docs/demo/02-integrity.gif)

<sub>▶ mpox 1만 3천 건으로 본 정합성 탭 — 필드별 일치 막대 → 실제로 다른 값 → 서열 검증 계획(dry-run)</sub>

</div>

### 2. 서열로 확인 (seqcheck)

메타데이터로 판단이 안 되면 서열 본문을 비교합니다 (헤더·공백 제거, 대문자, SHA-256). 기준은 NCBI.
몇백 MB를 막 받지 않도록 **먼저 dry-run**을 합니다: 링크와 길이만 조회해서 "길이가 다르면 받을 것도 없이 다름 확정", "같으면 받아서 비교"로 나눕니다. mpox 404건에서 약 314MB를 아꼈습니다.
받은 서열은 저장하지 않고 결과만 남깁니다.

### 3. 소스의 구조 변화 (schema-check)

소스의 응답 구조는 예고 없이 바뀝니다. **우리 코드가 실제로 읽는 필드를 측정**해서 사라지거나(이름이 바뀐 후보도 알려 줌) 형태가 바뀌거나 응답 구조가 달라지면 "소스 코드 업데이트 필요"라고 말해 줍니다.
매일 새벽 3시에 바이러스마다 소스당 3건씩 가볍게 점검하고, 바뀐 게 있으면 **스키마** 탭에 `!`가 켜집니다.

결과는 두 가지 모양으로 볼 수 있습니다. **간단히 보기**(기본)는 내 `genome` 열마다 NCBI·BV-BRC·Pathoplexus 중 어느 소스의 어떤 필드가 채우는지를 한 표에 색 알약으로 보여 주고, 확인이 필요한 열을 맨 위에 모아 줍니다. **자세히 보기**는 **테이블 구조도**입니다. 왼쪽에 소스 3곳의 필드(NCBI는 중첩 객체 그대로), 오른쪽에 내 `genome` 테이블의 열을 놓고 *코드가 실제로 읽어서 채우는 연결*을 선으로 잇습니다. 연결은 손으로 적은 표가 아니라 어댑터에 값을 넣어 보며 **측정**한 것이라 코드가 바뀌어도 어긋나지 않습니다.
사라진 필드는 <b>빨강</b>, 그 필드가 채우던 열은 <b>주황</b>(다른 소스가 대신 채울 수 있으면 오류가 아니라 경고), 새로 생긴 필드는 <b>파랑</b>으로 표시되고, 열을 누르면 어느 소스의 어느 필드가 우선순위 순으로 채우는지 나옵니다.

<div align="center">

![스키마 탭: 간단히 보기 표와 구조도, 그리고 소스가 바뀌었다고 가정했을 때의 강조](docs/demo/03-schema.gif)

<sub>▶ 스키마 점검의 간단히 보기와 테이블 구조도 — 처음은 실제 결과(이상 없음), 뒤쪽은 <b>시연을 위해 "소스가 필드를 지우고 이름과 응답 구조를 바꿨다"고 가정</b>한 화면입니다 (실제로 바뀐 것이 아닙니다)</sub>

</div>

### 4. 다른 값 중에서 고르기 (resolve)

소스끼리 값이 다르면 **검토 화면에서 항목을 눌러** 무엇이 다른지 보고 직접 고릅니다.

* 오른쪽 패널에 다른 필드마다 **소스별 값**이 나란히 나오고, 표기 차이인지 진짜 불일치인지 표시됩니다.
* **★ 추천**이 미리 골라져 있고 **이유와 신뢰도**가 붙습니다. 규칙은 이렇습니다: 불가능한 날짜는 뺀다 → 길이·버전은 NCBI(서열을 NCBI에서 받아 검증하므로) → 같은 사실의 다른 표기는 가장 잘 쓴 쪽(가장 구체적인 날짜, 학명, 표준 국가명) = **높음** → 소스 다수가 같은 값 = **보통** → 모두 다르면 우선순위로 제안하되 **낮음**(직접 확인).
* 값을 고르거나 직접 입력해서 저장하면 그 항목에 **바로 반영**되고, **다음 `check`에서도 유지**됩니다(소스가 바뀌어도 내가 정한 값을 존중). 언제든 기본값으로 되돌릴 수 있습니다.
* 툴바의 **Auto-pick**은 신뢰도 "보통" 이상인 추천을 한 번에 선택하고, **"낮음"은 사람이 보도록 남겨 둡니다.** 어차피 소스가 둘뿐인 바이러스는 대부분 "낮음"이니, 근거 없는 일괄 적용을 막는 장치입니다.

<div align="center">

![검토 화면에서 항목을 눌러 소스별 값을 보고, 추천을 참고해 고르고, Auto-pick으로 일괄 선택하기](docs/demo/05-resolve.gif)

<sub>▶ 검토 화면: 항목을 누르면 소스별 값과 추천(이유·신뢰도) → 고르기 → "선택 완료" → Auto-pick</sub>

</div>

### 5. 내 데이터가 멀쩡한지 (verify)

FASTA 실제 길이·md5가 기록과 맞는지, 날짜가 미래가 아닌지, 고아 행이 없는지 등을 봅니다. 문제가 있으면 종료 코드 1.

## 명령어 모음 (Cheat sheet)

| 명령 | 설명 |
|---|---|
| `watch-add <이름>` | 바이러스 등록 (`--taxon`으로 NCBI 이름·ID 지정) |
| `watch-list` | 등록된 바이러스 |
| `check <이름>` | 소스 3곳 조회 → changeset 생성 (데이터는 안 바뀜) |
| `dry-run <번호>` | changeset의 다운로드 계획 점검 |
| `apply <번호>` | 승인된 항목 적용 (`--approve-all`은 CLI 지름길, 웹 검토를 권장) |
| `reconcile <이름>` | 마지막 check의 소스 간 대조 보고서 |
| `verify [--no-files]` | 로컬 데이터 무결성 검사 |
| `schema-check [--daily]` | 소스·로컬 스키마 점검 (`--daily`는 스케줄러가 쓰는 가벼운 판) |
| `backup [zip]` / `restore <zip>` | 이사할 때: 큐레이션 결과 + FASTA만 담은 이식 가능한 백업 |
| `serve [--port N] [--open] [--no-scheduler]` | 웹 서버 |

```bash
python -m genomecurator --help        # 정확한 옵션은 여기서
```

## 설정 (환경변수)

| 변수 | 기본값 | 뜻 |
|---|---|---|
| `GC_DATA_DIR` | `<저장소>/data` | 데이터 위치 (Windows `set`, macOS/Linux `export`) |
| `NCBI_API_KEY` | (없음) | 있으면 NCBI 요청 한도가 늘어납니다 |
| `GC_CHECK_INTERVAL_HOURS` | `24` | 자동 확인 주기 |
| `GC_REQUEST_DELAY_MS` | `300` | 모든 외부 요청 사이 최소 간격 (예의 바른 손님) |
| `GC_NEAR_COMPLETE_PCT` | `90` | Near complete 기준: NCBI가 Complete가 아니라고 한 유전체 중 길이가 기준 게놈의 이 % 이상이면 "준완전" |
| `GC_DRYRUN_SAMPLE` | `3` | 적용 전 점검에서 샘플로 받아 볼 서열 수 |
| `GC_SCHEMA_CHECK_HOUR` | `3` | 매일 스키마 점검 시각(시) |
| `GC_SCHEMA_CHECK_SAMPLE` | `3` | 매일 점검에서 바이러스·소스당 읽을 레코드 수 |
| `GC_SEQCHECK_MAX` | `2000` | 서열 검증 한 번에 다룰 최대 accession 수 |
| `GC_MIN_FREE_BYTES` | 500 MB | 이보다 디스크 여유가 적으면 적용 전 점검이 막습니다 |
| `GC_LOG_LEVEL` / `GC_LOG_CONSOLE_LEVEL` | `INFO` / `WARNING` | 로그 파일 / 콘솔 수준 |

## 데이터는 어디에 (Output files)

```
data/
├── genomecurator.db          # SQLite: 큐레이션된 표, 출처 기록, 이력, 대조 결과
├── genomes/<바이러스>/*.fna  # 받아 둔 서열 (UTF-8, LF, 모든 OS 동일)
├── snapshots/                # 적용 직전 파일 백업 (실패 시 되돌리는 용도)
├── staging/                  # 적용 중 임시 파일
└── logs/genomecurator.log    # 서버 로그 (5MB씩 5개 회전, job/watch 맥락 포함)
```

`data/`는 용량이 크고 공개하면 안 되는 것이 섞여 있어서 **git에 올라가지 않습니다**. 백업 zip도 마찬가지입니다.

## 로그 읽는 법

한 줄이 이렇게 생겼습니다.

```
2026-10-17 07:14:15,374 INFO  genomecurator.curate [job=add839d0 watch=Henipavirus nipahense] check done in 27.8s
```

시각, 수준, 모듈, `[작업ID 바이러스 changeset]`, 내용. 외부 호출은 URL(키는 가림)·상태·크기·소요 시간이 남고, 실패는 traceback까지 남습니다. 웹의 **이력 → 서버 로그**에서 같은 파일을 볼 수 있습니다.

## Superset으로 다차원 분석

메타데이터를 바이러스·숙주·국가·연도·clade로 이리저리 돌려 보고 싶다면 Apache Superset을 붙이세요. DB는 **읽기 전용**으로 연결되고, 평탄화된 뷰(`v_genome_analysis`)가 준비돼 있습니다.

```bash
mamba env create -f superset/environment.yml
scripts/superset.sh init && scripts/superset.sh run          # http://127.0.0.1:8089
scripts/superset.sh bootstrap                                # 연결·데이터셋·차트·대시보드 자동 생성
```

자세한 내용은 [docs/SUPERSET.md](docs/SUPERSET.md).

## 공개용 정적 페이지 (읽기 전용)

`docs/index.html` 한 파일에 개요 차트와 유전체 표(필터·정렬·상세·TSV)가 들어 있습니다. 공개 메타데이터만 담기고 서열·경로·이력·수집 코드는 들어 있지 않습니다.

```bash
python scripts/build_static.py     # 로컬 DB로 docs/index.html 다시 만들기
```

앱의 **배포하기** 버튼은 미리보기로 바이러스·유전체 수를 보여 주고 확인을 받은 뒤, 그 파일만 커밋·push합니다. GitHub Pages는 저장소 Settings → Pages → `main` 브랜치 `/docs`로 지정하세요. 차트는 CDN(Chart.js)을 쓰므로 보는 쪽에 인터넷이 필요합니다.

## 이사 가기 (백업 / 복원)

```bash
python -m genomecurator backup my-backup.zip      # 큐레이션 결과 · 결정 이력 · FASTA (원본 응답과 대기 중 항목은 제외)
python -m genomecurator restore my-backup.zip     # 새 컴퓨터의 data 폴더로 (기존 DB가 있으면 --force)
```

Windows ↔ macOS ↔ Linux 어디서든 풀립니다. 제외된 것은 다음 "지금 확인"에서 다시 만들어집니다.

## 자주 묻는 질문 (FAQ)

**Q. 승인을 안 하고 그냥 두면요?**
아무 일도 안 일어납니다. 그게 설계입니다. 상단 *큐레이션* 탭의 배지가 "검토 안 한 변경 사항이 있어요"라고 조용히 눈치를 줍니다.

**Q. "Pathoplexus: no records (skipped)"가 보여요.**
그 바이러스가 Pathoplexus에 없다는 뜻입니다 (예: Nipah). 정상입니다. 반면 Pathoplexus가 **응답하지 않았다면** 건너뛰지 않고 확인 전체를 중단합니다. 둘은 일부러 다르게 취급합니다.

**Q. 새 Pathoplexus organism은 어떻게 추가하나요?**
따로 설정할 필요가 없습니다. 사이트에서 organism 목록을 읽어 와서, 등록한 바이러스의 NCBI taxon과 맞는 것을 자동으로 연결합니다.

**Q. 불일치가 많은데 다 봐야 하나요?**
아니요. 정합성 탭에서 *표기 차이*는 이미 걸러져 있고, 검토 화면의 **정합성** 열로 "불일치"만 골라 볼 수 있습니다. 항목을 누르면 소스별 값과 **★ 추천**이 나오고, 툴바의 **Auto-pick**이 신뢰도가 높은 추천을 한 번에 골라 줍니다 ("낮음"만 직접 보면 됩니다). 서열 검증 버튼은 길이·버전이 다른 건만 대상으로 삼습니다.

**Q. 스키마 점검이 "코드 업데이트 필요"라고 해요.**
오류 항목의 "조치"에 고칠 함수 이름이 적혀 있습니다 (예: `sources.normalize_ncbi`). 고친 뒤 다시 점검하고 *새 기준으로 저장*하세요. 고치지 않은 채 기준을 갱신해 문제를 덮는 일은 막아 두었습니다.

**Q. 서버를 꺼 두면 매일 점검은요?**
그날은 건너뜁니다. 서버가 켜져 있어야 스케줄러가 돕니다.

## 한계 (솔직하게)

* 서버가 켜져 있어야 자동 확인과 매일 스키마 점검이 돕니다.
* 매일 점검의 NCBI·Pathoplexus 샘플은 "처음 N건"(오래된 레코드)이라 최신 레코드에만 생기는 변화는 늦게 발견될 수 있습니다. 가끔 *스키마* 탭에서 전체 점검을 눌러 주세요.
* 서열이 같다는 것은 "같은 레코드"라는 증명이지, 숙주·날짜 같은 메타데이터가 *옳다*는 증명은 아닙니다.
* 분류 이름(`genome_name`)은 소스 원본 그대로라 `Monkeypox virus`와 `Orthopoxvirus monkeypox`가 섞여 보일 수 있습니다 (수집 범위에는 영향 없음).

## 코드 지도

| 파일 | 역할 |
|---|---|
| `genomecurator/sources.py` | 세 소스 어댑터 (NCBI, BV-BRC, Pathoplexus) + 정규화 |
| `genomecurator/curate.py` | check → changeset → dry-run → apply |
| `genomecurator/reconcile.py` | 소스 간 일치/표기차이/불일치 분류 |
| `genomecurator/seqcheck.py` | 서열 본문 비교(dry-run 계획 + 실행) |
| `genomecurator/schema_check.py` | 소스·로컬 스키마 점검, 매일 점검 |
| `genomecurator/verify.py` | 로컬 데이터 무결성 |
| `genomecurator/web.py` + `templates/` | FastAPI 웹 UI |
| `genomecurator/log.py` | 로깅 (파일 회전 + job/watch 맥락) |
| `genomecurator/db.py`, `config.py` | 스키마·마이그레이션, 설정 |
| `scripts/` | 시작 스크립트, 정적 페이지 생성, Superset 도우미 |
| `CHANGELOG.md` · `CHANGELOG.ko.md` | 버전별 변경 이력 (영문 · 한글) |
| `docs/RELEASING.md` · `docs/DEPENDENCIES.md` | 릴리스 절차, 사용한 라이브러리 |
| `docs/STATUS.md` | 진행 상황과 결정 사항 (다른 곳에서 이어서 작업할 때 먼저 읽을 문서) |

개발할 때는 `python -m pytest -q`. CI는 Windows·macOS·Linux × Python 3.10/3.12, 그리고 `environment.yml`이 실제로 빌드되는지도 매번 확인합니다.

## 이 GIF들은 어떻게 만들었나

위 GIF는 손으로 찍은 게 아니라 **스크립트가 실제 앱을 조종해서** 만듭니다 (`scripts/demo/`). 헤드리스 Chrome을 DevTools 프로토콜로 움직이면서 프레임을 모으고, 커서·자막을 화면 위에 얹어 GIF로 조립합니다. 앱이 바뀌면 GIF도 다시 뽑으면 됩니다.

```bash
pip install -e ".[demo]"                       # Pillow + websocket-client (Google Chrome도 필요)
python scripts/demo/make_gifs.py quickstart    # 빈 폴더에서 시작, 실제 인터넷 접속, 3~4분
python scripts/demo/make_gifs.py explore       # quickstart가 만든 데이터를 보여 줌
```

시연은 임시 복사본에서 돌기 때문에 진짜 `data/`는 건드리지 않습니다. `integrity`와 `schema`는 미리 준비한 시연 데이터가 필요합니다 (`GC_DEMO_ROOT`, 자세한 내용은 `make_gifs.py` 맨 위 설명).

## 버전과 변경 이력

현재 버전은 **0.2.0** 입니다. 버전은 `genomecurator/__init__.py` 한 곳에서만 정하고, 화면 하단·`/api/version`·`--version`에 보입니다.
또한 **DB·changeset·백업**이 "어느 버전이 만들었는지"를 기록하므로, 나중에 이 데이터가 어느 코드에서 나왔는지 추적할 수 있습니다.

```bash
python -m genomecurator --version      # genomecurator 0.2.0 (DB schema 1)
```

* 변경 이력: [CHANGELOG.ko.md](CHANGELOG.ko.md) (한국어) · [CHANGELOG.md](CHANGELOG.md) (English)
* 릴리스 방법과 버전 규칙: [docs/RELEASING.md](docs/RELEASING.md)

## 사용한 라이브러리

| 라이브러리 | 하는 일 | 라이선스 |
|---|---|---|
| [FastAPI](https://fastapi.tiangolo.com/) · [Uvicorn](https://www.uvicorn.org/) | 웹 앱과 서버 | MIT · BSD-3 |
| [Jinja2](https://jinja.palletsprojects.com/) | HTML 템플릿 | BSD-3 |
| [HTTPX](https://www.python-httpx.org/) | 소스 3곳 호출 (재시도·시간제한) | BSD-3 |
| [APScheduler](https://apscheduler.readthedocs.io/) | 매일 확인·스키마 점검 | MIT |
| [python-multipart](https://github.com/Kludex/python-multipart) | 웹 폼 해석 | Apache-2.0 |
| [Chart.js](https://www.chartjs.org/) 4.4.1 (CDN) | 대시보드·공개 페이지 차트 | MIT |
| 선택: [pytest](https://pytest.org/) · [Pillow](https://python-pillow.org/) · [websocket-client](https://github.com/websocket-client/websocket-client) · [Apache Superset](https://superset.apache.org/) | 테스트 · GIF 만들기 · 분석 | MIT · MIT-CMU · Apache-2.0 · Apache-2.0 |

테스트한 버전, 전체 라이선스, 간접 의존성, 데이터를 가져오는 외부 서비스(NCBI·BV-BRC·Pathoplexus·ICTV)는 [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)에 있습니다
(설치된 패키지에서 `scripts/dependencies.py`가 뽑아 만든 문서라 코드와 어긋나지 않습니다).

## 라이선스 · 만든 사람 · 문의

MIT License. © 2026 CHOICHI-HWAN.
버그와 아이디어는 [Issues](https://github.com/CHOICHI-HWAN/GenomeCurator/issues)로 보내 주세요.

## 같이 쓰면 좋은 도구

* [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) · [BV-BRC](https://www.bv-brc.org/) · [Pathoplexus](https://pathoplexus.org/) — 데이터를 내주는 고마운 곳들입니다. 요청 사이에 간격을 두는 이유입니다.
* [Snippy](https://github.com/tseemann/snippy) — 이 README의 구성이 닮아 있다면 그 덕입니다. 변이 호출은 그쪽이 훨씬 잘합니다.
* [Apache Superset](https://superset.apache.org/) — 분석 화면을 맡습니다.
