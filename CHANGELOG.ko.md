# 변경 이력 (Changelog)

GenomeCurator의 눈에 띄는 변경은 모두 여기에 기록합니다. 형식은 [Keep a Changelog](https://keepachangelog.com/ko/1.1.0/)를,
버전은 [유의적 버전(SemVer)](https://semver.org/lang/ko/)을 따릅니다 (메이저가 0인 동안은 마이너 버전이 올라갈 때 동작이 바뀔 수 있고,
패치 버전은 버그 수정만 합니다). 릴리스 방법은 [docs/RELEASING.md](docs/RELEASING.md).
**English:** [CHANGELOG.md](CHANGELOG.md)

읽는 법: `추가`는 새 기능, `변경`은 이전부터 있던 동작의 변화, `수정`은 버그, `데이터·호환성`은 DB·백업·큐레이션 데이터에 닿는 변화입니다
(업그레이드 전에 미리 알 수 있도록 따로 적습니다), `의존성`은 추가·삭제·버전 고정이 바뀐 라이브러리입니다(라이선스까지 담은 전체 목록은
[docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)).

## [Unreleased]

### 추가
- **스키마** 탭의 **테이블 구조도**: 왼쪽에 NCBI(중첩 객체는 그룹으로 유지)·BV-BRC·Pathoplexus의 필드, 오른쪽에 내 `genome` 테이블의 열을 놓고,
  어느 소스 필드가 어느 열을 채우는지 선으로 잇습니다. 연결은 손으로 적은 것이 아니라 *측정*한 값입니다. 실제 어댑터에 필드를 하나씩만 넣고 실행해
  바뀌는 열을 읽기 때문에 "A가 없으면 B" 같은 대체 규칙과 한 필드가 두 열을 채우는 경우도 찾아냅니다. 사라진 필드는 빨강, 그 필드가 채우던 열은 주황
  (다른 소스가 여전히 채우면 경고, 채우는 소스가 하나도 없으면 오류), 새 필드는 파랑으로 표시합니다. 마우스를 올리면 연결이 강조되고, 열을 누르면 그 열을
  채우는 소스 필드가 우선순위 순으로 형태·채움 비율·예시와 함께 나옵니다. 응답 구조가 바뀐 소스는 빨간 배너가 붙습니다.

### 변경
- 스키마 탭은 이제 전체 점검과 그 구조도가 맨 위에 오고, 매일 바이러스별 점검은 그 아래로 옮겼습니다.
- 사라진 필드를 찾을 때 지난 정상 기준에서 코드가 읽던 필드도 함께 고려합니다. 코드가 읽는 필드는 데이터에 따라 달라서(`host` 객체가 없으면
  `host.organism_name`이 안 보임), 부모가 사라지면 자식이 가려지는 문제가 있었습니다.

### 수정
- 스키마 점검이 객체 목록(NCBI 레코드의 lineage 등)의 채움 비율을 100% 넘게 보고하던 문제: 목록의 원소마다 세던 것을 레코드당 한 번으로 고쳤습니다.

## [0.1.0] - 2026-10-17

첫 릴리스. 공개 소스 세 곳에서 바이러스 유전체를 모으고, 모든 변경을 사람이 검토·승인하게 하며, 소스끼리·서열과·코드의 기대와 맞는지 확인합니다.

### 추가

**수집과 병합**
- 바이러스 워치리스트: 이름으로 등록합니다. 먼저 ICTV · NCBI Virus · BV-BRC · Pathoplexus에서 찾아 보고, 바이러스가 아니면 후보와 함께 거절합니다.
  편집·삭제(이름을 직접 입력해 확인, 영향 범위 미리보기)가 됩니다. 균주·아종 이름은 종(species)으로 넓혀서 종 노드에 붙은 레코드를 놓치지 않습니다.
- 소스 세 곳. **NCBI Virus**(Datasets v2)가 어떤 accession이 있는지 정하고, **BV-BRC**와 **Pathoplexus**가 메타데이터를 보강합니다.
  Pathoplexus organism은 사이트에서 직접 찾아 NCBI taxon(하위 분류 포함)으로 연결하고, organism별 필드(clade, lineage, serotype, outbreak,
  genotype, subtype)도 자동으로 읽습니다.
- 필드 우선순위는 Pathoplexus > NCBI Virus > BV-BRC이고, 필드마다 어느 소스 값을 썼는지 기록합니다(`field_provenance`).
- SRA, BioProject, BioSample, 시퀀싱 플랫폼, 어셈블리 방법 열과 NCBI 링크.
- 모든 웹 요청은 최소 300ms 간격이고 실패하면 지수 백오프로 최대 3회 재시도합니다.

**검토와 적용 (changeset)**
- `check`는 로컬 데이터를 바꾸지 않습니다. 신규 / 변경 / 삭제된 레코드를 *changeset*으로 만들 뿐이고(변경은 메타데이터가 바뀌거나 GenBank `.version`이
  올라간 경우), 사람이 승인하기 전에는 아무것도 내 데이터가 되지 않습니다. 거절한 항목은 내용이 같으면 다시 올라오지 않습니다.
- 엑셀 같은 검토 표에서 승인 / 거절 / 보류. 1,000행이 넘는 선택과 "승인 + 적용" 한 번 클릭도 됩니다.
- 적용 전에는 항상 dry-run(디스크 여유, 링크, 샘플 서열)을 합니다. 적용은 원자적입니다: 파일을 먼저 백업하고 예외가 나면 모두 되돌리며,
  개별 accession 실패는 건너뛰고 보고하고, 끊긴 다운로드는 임시 파일에서 이어받습니다.
- 메타데이터 먼저 반영하고 FASTA는 나중에(*Get FASTA*), 메타데이터만 바뀐 항목은 FASTA를 받지 않습니다.
- changeset 삭제(영향 미리보기와 확인). 검토 화면은 단계별 안내와 함께 쉬운 말로 씌어 있습니다.

**데이터 정합성**
- `check`마다 소스 간 대조: 필드를 *일치*, *표기 차이*(같은 사실을 다르게 적음. 예 `USA: CA` ↔ `USA, CA`, `2022-07` ↔ `2022-07-11`), *불일치*로 분류하고,
  어느 소스가 어떤 accession을 갖는지도 봅니다. **정합성** 탭과, 검토 중 열·필터·경고로 보입니다.
- 서열 본문 검증: 길이·버전·해시가 다른 accession의 서열 본문(헤더와 공백을 빼고 대문자로 맞춘 뒤 SHA-256)을 NCBI·BV-BRC·Pathoplexus에서 받아 NCBI 기준으로 비교합니다.
  먼저 dry-run 계획을 보여 주고(가벼운 길이 조회만으로 대부분 판정해 다운로드를 줄임), 이어서 실행할 수 있으며, NCBI와 똑같은 서열을 가진 소스는 그 쌍의 길이·버전 불일치를 해소합니다.
- `verify`: FASTA 길이·md5가 기록과 맞는지, 날짜, 고아 행 등.
- 응답하지 않는 소스가 있으면 건너뛰지 않고 확인을 중단합니다. 그 소스를 빼고 병합하면 큐레이션된 필드가 비워질 수 있기 때문입니다.

**소스 스키마 감시**
- 스키마 dry-run(`schema-check`, **스키마** 탭): 어댑터를 실제 응답 샘플에 돌려 *코드가 실제로 읽는 필드*를 측정하고, 지난 기준과 비교해 사라진·형태가 바뀐·새로 생긴 필드와
  응답 구조 변화를 찾으며, 로컬 DB가 코드의 기대와 맞는지도 확인합니다. 소스 코드 업데이트가 필요한지, 어디를 고쳐야 하는지 알려 줍니다.
- 문제가 있는 기준은 사람이 새 상태를 수락하기 전까지 바꾸지 않습니다.
- 매일(03시) 바이러스별로 소스당 3건씩 가볍게 점검하고, 확인이 필요하면 **스키마** 탭에 배지가 켜집니다.

**웹 화면**
- 필터·정렬·열 선택·TSV 내보내기·전체 화면이 되는 엑셀형 표, 값마다 출처 배지와 Pathoplexus 변이 정보를 보여 주고 너비를 끌어서 조절하는 상세 패널,
  대시보드 차트, 모든 화면이 같은 틀.
- 큰 진행 막대와 서버가 보관하는 작업 로그를 보여 주는 진행 모달(새로고침해도 이어짐), 서버 로그 보기가 있는 이력 화면.

**분석과 공개**
- Apache Superset 연동: 읽기 전용 평탄화 뷰(`v_genome_analysis`), `scripts/superset.sh`, API로 만드는 차트 10개와 대시보드.
- 읽기 전용 공개 페이지(`docs/index.html`, 로컬 `/public`에서도 열림)와 GitHub Pages 원클릭 게시.

**운영**
- Windows·macOS·Linux: OS와 무관한 데이터 폴더, UTF-8 / LF 파일, 이식 가능한 상대 경로 FASTA, 안전한 폴더 이름, conda/mamba 환경을 우선 쓰는 시작 스크립트
  (`environment.yml`, `superset/environment.yml`). CI는 세 OS 행렬로 돌립니다.
- 이식 가능한 백업·복원(`backup` / `restore`): 큐레이션 결과, 결정 이력, FASTA. 더 새로운 DB 구조로 만든 백업은 거절합니다.
- 구조화된 로깅: job / watch / changeset 맥락이 붙는 회전 로그 파일, 모든 HTTP 시도(비밀값은 가림), 실패한 작업이 롤백돼도 남는 이력 기록.
- 버전 관리: 버전은 한 곳(`genomecurator.__version__`)에서만 정하고, `genomecurator --version`, 화면 하단, `/api/version`이 보여 주며,
  DB·모든 changeset·모든 백업이 어느 버전이 만들었는지 기록합니다.

**문서**
- 실제 앱을 스크립트로 녹화한 시연 GIF가 들어간 새 README(`scripts/demo/`), 이 변경 이력, 릴리스 절차([docs/RELEASING.md](docs/RELEASING.md)),
  사용한 라이브러리 목록([docs/DEPENDENCIES.md](docs/DEPENDENCIES.md), 설치된 패키지에서 뽑아 생성).

### 의존성
- Python **3.10 이상** (3.10과 3.12에서 테스트, conda 환경은 3.12).
- 실행: `fastapi>=0.110`, `uvicorn>=0.29`, `jinja2>=3.1`, `httpx>=0.27`, `apscheduler>=3.10,<4`, `python-multipart>=0.0.9`.
- 선택: `pytest>=8`(테스트), `pillow>=10`·`websocket-client>=1.6`(README GIF 녹화, Google Chrome 필요),
  Apache Superset `6.1.0`과 `rich`·`cachetools`·`flask-caching<2.4`(분석, 별도 환경).
- 브라우저: jsDelivr CDN의 Chart.js 4.4.1(차트에는 인터넷이 필요).
- 데이터는 NCBI Datasets, BV-BRC, Pathoplexus(LAPIS)에서, 이름 확인에는 EBI OLS4를 거친 ICTV에서 가져옵니다.

### 데이터·호환성
- DB 구조(스키마 버전 **1**): 지금까지의 변경은 모두 앱이 DB를 열 때 자동으로 적용되는 추가형 마이그레이션이라, 이전 개발 빌드의 DB도 그대로 업그레이드됩니다.
  0.1.0이 DB를 처음 열 때 이력(`version`)과 `app_meta`에 자신을 기록합니다.
- 이 버전이 만든 백업의 `manifest.json`에는 `app_version`과 `db_schema_version`이 들어 있습니다.

### 알려진 한계
- 매일 스키마 점검은 NCBI·Pathoplexus가 돌려주는 *처음* 레코드(오래된 것)를 샘플로 쓰므로 새 레코드에서만 나타나는 변화는 늦게 발견될 수 있습니다.
  가끔 **스키마** 탭에서 전체 점검을 실행하세요. BV-BRC는 최신순으로 뽑습니다.
- 서열 본문이 같다는 것은 두 레코드가 같은 대상이라는 증명이지, 숙주·날짜 같은 메타데이터가 맞다는 증명은 아닙니다.
- 스케줄러(매일 점검)는 서버가 켜져 있는 동안에만 돕니다.

[Unreleased]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/releases/tag/v0.1.0
