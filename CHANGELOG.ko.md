# 변경 이력 (Changelog)

Virus Genome Metadata Curator(패키지 `genomecurator`, 이전 표시 이름 "GenomeCurator")의 눈에 띄는 변경은 모두 여기에 기록합니다. 형식은 [Keep a Changelog](https://keepachangelog.com/ko/1.1.0/)를,
버전은 [유의적 버전(SemVer)](https://semver.org/lang/ko/)을 따릅니다 (메이저가 0인 동안은 마이너 버전이 올라갈 때 동작이 바뀔 수 있고,
패치 버전은 버그 수정만 합니다). 릴리스 방법은 [docs/RELEASING.md](docs/RELEASING.md).
**English:** [CHANGELOG.md](CHANGELOG.md)

읽는 법: `추가`는 새 기능, `변경`은 이전부터 있던 동작의 변화, `수정`은 버그, `데이터·호환성`은 DB·백업·큐레이션 데이터에 닿는 변화입니다
(업그레이드 전에 미리 알 수 있도록 따로 적습니다), `의존성`은 추가·삭제·버전 고정이 바뀐 라이브러리입니다(라이선스까지 담은 전체 목록은
[docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)).

## [Unreleased]

### 추가
- **다른 프로그램에 넘기는 내보내기 탭·CLI·API**: 같은 id·같은 순서의 `metadata.tsv` + `sequences.fasta`, 그리고 `manifest.json`(건수, 경고, 두 파일의 SHA-256, Nextclade 데이터셋 버전, 데이터 이용 안내)과 바로 쓸 수 있는 `nextstrain-config.yaml`. **nextstrain** 형식은 Nextstrain mpox 워크플로가 읽는 열(id는 `INSDC_accession_version`, `2023-XX-XX` 날짜, region/country/division/location, host, clade/lineage/outbreak, 기본 필터가 쓰는 Nextclade `QC_*` 열)에 맞추고, **plain**은 큐레이션된 모든 열을 줍니다. 필터: 최소 길이, 완성도, QC 불량·오류 제외, clade, 서열 유무. `python -m genomecurator export <바이러스>`, `GET /api/export/{virus}/metadata.tsv|sequences.fasta|preview`, `POST /api/export/{virus}/build`, `GET /api/exports`. 공개(Open) 데이터만 나갑니다.
- `python -m genomecurator qc <바이러스> [--force]`. Nextclade 결과가 Nextstrain 워크플로에서 쓰는 이름의 QC 열도 함께 저장합니다(이전에 만든 결과는 QC를 한 번 다시 실행하면 채워집니다).
- **QC 탭(서열 QC)**: 저장된 모든 유전체의 서열 통계(길이, N·모호 염기 비율, GC, 기록된 길이와의 일치)와, Nextclade 데이터셋이 있는 바이러스(mpox, 에볼라 종들, 홍역, 뎅기, 웨스트나일, 풍진, VZV, hMPV; `GC_NEXTCLADE_DATASETS`로 추가)는 **Nextclade**(clade, lineage, 자체 품질 점수). 경고만 합니다. 결과는 새 테이블 `qc_result`에 저장되는 파생 데이터입니다.
- Nextclade의 clade/lineage는 큐레이션에서 네 번째 *소스*로 참여합니다: 빈 값은 채우는 후보가 되고, NCBI / BV-BRC / Pathoplexus와 **다르면** 아무것도 덮어쓰지 않고 검토 화면에 추천(데이터베이스가 같은 값이면 높음, 아니면 보통)과 함께 올라옵니다. *QC 실행 + 확인*은 두 단계를 한 작업으로 합니다.
- Nextclade 데이터셋은 필요할 때 내려받고, 바이러스마다 데이터셋 버전과 새 버전 유무를 보여 주며 **데이터셋 업데이트** 버튼이 있습니다. 유전체 표에 QC / Nextclade clade / Nextclade lineage 열과 QC 필터가 생겼습니다.
- `environment.yml`에 `bioconda` 채널과 `nextclade`(MIT, x86-64·ARM)를 넣었습니다. Nextclade가 없으면 서열 통계만 동작합니다. 데이터셋은 라이선스 파일 없이 배포되어 재배포하지 않습니다.
- **소스 신호등**: 모든 페이지 맨 위에 NCBI — Pathoplexus — BV-BRC를 선으로 이어 최근 스키마 점검 결과를 불빛으로 보여 줍니다(파랑 = 변화 없음, 노랑 = 확인 필요, 빨강 = 소스 코드 업데이트 필요, 회색 = 점검 전). 마우스를 올리면 이유가, 누르면 스키마 탭이 열립니다.
- **선택한 항목에 추천값 일괄 적용**: 검토 화면에서 행을 고르고(체크박스, **마우스로 행을 끌어서**, **Shift+클릭**으로 범위, 또는 정합성 필터 후 *Select all*) *선택 항목에 추천 적용*을 누릅니다. 신뢰도 구성을 보여 주는 창에서 *보통 이상*(권장) 또는 *낮음까지 모두*를 고릅니다. 전체 버튼은 *Auto-pick (전체)*로 이름을 바꿨습니다.
- **분석 탭**: 앱에서 Apache Superset을 켜고 끄고(여전히 별도 conda 환경의 별도 프로그램), 설치·실행 여부를 보고 바로 엽니다. `serve --superset`은 앱과 함께 켭니다. 앱이 켠 Superset은 앱이 꺼질 때 같이 꺼집니다.
- `Dockerfile`(`environment.yml` 기반, linux/amd64·linux/arm64 공용), 두 아키텍처를 빌드하고 버전 태그 때 GHCR에 게시하는 `docker` 워크플로, CI의 ARM 작업.
- **작업 탭**: 오래 걸리는 모든 작업(바이러스 추가·조회, 확인, 적용, FASTA 받기, 서열 검증, 스키마 점검, 배포, 스케줄러 실행)이 서버 작업으로 기록됩니다: 진행 중 / 완료 / 오류 / 중단, 진행률, 시각이 붙은 로그, 오류와 발생 위치, 해당 작업의 서버 로그 링크. 진행 중인 작업 수나 못 본 오류가 탭에 표시되고, 다른 페이지에 있을 때 작업이 끝나면 알림이 뜹니다.
- 진행 창에 **백그라운드로** 버튼.
- 서버 로그 검색(예: `job=ab12cd34`).

### 변경
- 공개 페이지 푸터에 출처별 정책 링크(NCBI 고지, Pathoplexus 이용 약관, BV-BRC)와 Pathoplexus는 공개 데이터만 포함한다는 문구를 넣었습니다. 새 문서: `docs/LICENSES.md`(라이선스·데이터 이용 검토), `docs/WORKFLOW.md`.
- README를 짧고 실용적으로(기능, 화면, 설치, 실행, 출력) 다시 썼습니다. 영문 `README.md`와 한글 `README.ko.md`로 나누고, 긴 설명은 `docs/GUIDE.ko.md`로 옮겼습니다.
- 바이러스 추가와 조회가 서버 작업으로 실행되어, 다른 메뉴를 눌러도 사라지지 않습니다(조회는 돌아오면 이어서 보입니다).
- 자동 확인과 매일 스키마 점검도 로그뿐 아니라 작업으로 기록됩니다.

### 수정
- 동시에 DB 연결이 열릴 때 `view v_genome_analysis already exists` 오류가 날 수 있었습니다. 분석 뷰는 정의가 바뀐 때만 한 트랜잭션으로 다시 만듭니다.

### 데이터·호환성
- 새 파생 테이블 `qc_result`(자동 생성, 백업에 포함되지 않으며 QC 탭에서 다시 계산). 스키마 버전은 1 그대로입니다.
- 새 운영용 테이블 `job`(자동 생성, 백업에는 포함되지 않음). 스키마 버전은 1 그대로입니다.

## [0.2.0] - 2026-10-17

유전체의 완성도와 소스끼리의 불일치를 *바로 손볼 수 있는* 것으로 만든 버전입니다: Near complete 등급, 추천과 함께 소스의 값 중에서 고르는 기능, 스키마 점검을 읽는 두 가지 방법.
도구가 **바이러스 전용**임도 분명히 했습니다.

### 변경
- 호스트가 괄호 속 일반명만 다른 경우(`Sus scrofa` ↔ `Sus scrofa (pig)`)는 이제 불일치가 아니라 *표기 차이*입니다.
- **범위와 이름.** 이 도구는 **바이러스 전용**입니다(NCBI Virus, BV-BRC, Pathoplexus). 세균은 어셈블리(contig 여러 개)이고 서열이 수 Mb라서 **별도 프로젝트**로 만듭니다.
  "genome"이 *사람 유전체*로 읽히지 않도록 화면에 보이는 이름을 **Virus Genome Metadata Curator**(페이지 헤더·제목·안내문은 **Virus Genome Curator**)로 바꿨습니다.
  저장소, 파이썬 패키지 `genomecurator`, 명령어, 데이터 폴더, 환경변수, 백업 형식 식별자는 **그대로**라서 옮길 것이 없습니다.

### 추가
- 검토 화면에서 **소스의 값 중 고르기**. 항목을 누르면 오른쪽 패널에 소스끼리 다른 필드가 모두 나오고, 소스별 값과 표기만 다른 것인지 진짜 다른 것인지를 보여 줍니다.
  **★ 추천**이 이유·신뢰도와 함께 미리 골라져 있습니다: 불가능한 날짜는 뺀다 → 길이·버전은 NCBI(서열을 NCBI에서 받아 길이를 검증하므로) → 같은 사실의 다른 표기는 가장 잘 쓴 쪽
  (가장 구체적인 날짜, 학명, 표준 국가명) = *높음* → 소스 다수가 같은 값 = *보통* → 그 외에는 소스 우선순위로 제안하되 직접 확인하라고 안내 = *낮음*.
  값을 고르거나 직접 입력해 저장하면 그 항목에 바로 반영되고 **이후 모든 확인에서 유지**됩니다(고른 소스의 현재 값을 쓰고, 직접 입력한 값은 그대로 두며, 출처는 *manual*로 표시). 필드마다 되돌릴 수 있습니다.
  툴바의 **Auto-pick**은 신뢰도 *보통* 이상의 추천을 한 번에 선택하고 *낮음*은 사람에게 남깁니다. 선택은 새 테이블 `field_choice`에 저장되고 백업에 포함되며, 사람이 정한 필드는 더 이상 열린 불일치로 세지 않습니다.
- **Near complete(준완전)** 유전체. NCBI는 Complete와 Partial만 구분하는데, "Partial" 대부분은 양 끝 몇백 염기만 모자란 것입니다. 그 사이에 *완성도* 등급을 새로 두었습니다:
  **Complete**(NCBI가 Complete로 표시한 것. 신뢰하므로 기준 게놈보다 원래 짧은 계통이 강등되지 않음), **Near complete**(Complete가 아니면서 길이가 기준 게놈의 **90% 이상**),
  **Partial**(그 미만이거나 기준 없음). 바이러스의 기준 게놈은 수집된 RefSeq 게놈 중 가장 긴 것입니다. 기준 비율은 `GC_NEAR_COMPLETE_PCT`(기본 90, CheckV high-quality 기준과 같고
  실제 데이터에서도 길이 분포의 골짜기에 해당: mpox 13,462건 중 50~90% 구간은 33건뿐)로 바꿀 수 있습니다. 열·필터·대시보드 KPI·공개 페이지·분석용 뷰(`completeness`, `ref_coverage_pct`, `ref_length`,
  SQL에도 같은 규칙)에 나타납니다. mpox에서 NCBI가 Partial로 표시한 8,223건 중 7,583건이 Near complete입니다.
- **스키마** 탭의 **간단히 보기**(기본, 탭으로 아래 자세한 구조도와 전환): 내 `genome` 열을 행, NCBI / BV-BRC / Pathoplexus를 열로 놓은 표 하나에
  그 열을 채우는 소스 필드 이름을 상태 색과 함께 보여 줍니다. 확인이 필요한 열은 위에 모아 어느 소스 필드를 봐야 하는지 알려 주고, 열은 묶음(식별, 분류, 유전체, 시료, 날짜,
  외부 ID)으로 나뉩니다. 고른 보기는 기억합니다.
- **스키마** 탭의 **테이블 구조도**(자세히 보기): 왼쪽에 NCBI(중첩 객체는 그룹으로 유지)·BV-BRC·Pathoplexus의 필드, 오른쪽에 내 `genome` 테이블의 열을 놓고,
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

### 데이터·호환성
- 스키마 버전은 **1** 그대로입니다. DB 변경은 새 테이블 `field_choice`(소스 값 중 사람이 고른 선택) 하나뿐이고, 0.2.0이 DB를 처음 열 때 자동으로 만들어지며 기존 것은 바뀌지 않습니다.
  0.1.0 백업은 0.2.0에 복원되고, 0.2.0 백업에는 `field_choice`가 들어갑니다.
- 분석용 뷰 `v_genome_analysis`에 `completeness`, `ref_coverage_pct`, `ref_length` 열이 늘었습니다. Superset에서는 `scripts/superset.sh bootstrap`을 한 번 실행해 데이터셋을 새로고침하면 보입니다.
  뷰는 Near complete 기준 비율을 `app_meta.near_complete_pct`에서 읽으며, 앱이 `GC_NEAR_COMPLETE_PCT`와 맞춰 둡니다.
- 새 환경변수: `GC_NEAR_COMPLETE_PCT`(기본 90).
- NCBI가 *Partial*로 표시한 유전체가 *Near complete*로 보일 수 있습니다. 저장된 값을 고쳐 쓰는 것은 아니고 보여 줄 때 계산합니다. NCBI 상태는 별도 열로 그대로 남아 있습니다.

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

[Unreleased]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/releases/tag/v0.1.0
