# Virus Genome Curator 진행 상황 / 이어서 작업하기

임시 환경에서 작업하므로 이 문서와 저장소만 있으면 어디서든 이어갈 수 있습니다.

## 실행
```bash
mamba env create -f environment.yml && mamba activate genomecurator   # conda 환경 'genomecurator' (Python 3.12)
python -m genomecurator watch-add "Henipavirus nipahense"
python -m genomecurator serve        # http://127.0.0.1:8000 (8000이 쓰이면 --port 8001)
pytest
```
- 한 번에 시작: `scripts/start.sh | start.command | start.bat` (mamba/conda 환경을 쓰고, 없으면 만들고, 둘 다 없으면 .venv).
- Superset은 별도 환경: `mamba env create -f superset/environment.yml` (`docs/SUPERSET.md`).
- 이 머신의 쉘 기본값은 `base`다. 직접 실행하려면 먼저 `mamba activate genomecurator`.
데이터(DB, 서열)는 `data/`에 생기며 git에는 올리지 않습니다 (`GC_DATA_DIR`로 위치 변경).
Windows에서 `git`이 PATH에 없으면 `C:\Program Files\Git\cmd\git.exe`로 실행하세요.

## 확정된 결정
- 소스 우선순위: Pathoplexus > NCBI Virus > BV-BRC. NCBI가 accession 집합을 정하고 나머지는 메타데이터 보강. 필드별 채택 소스는 `field_provenance`에 기록.
- 변경 판정: 메타데이터 변경, GenBank `.version` 증가, 신규, NCBI에서 사라짐(removed).
- 승인 전에는 로컬 데이터를 바꾸지 않음 (changeset → 승인 → Apply). 거절한 항목은 같은 내용이면 재등장하지 않음.
- 내부 포맷: BV-BRC 스타일 컬럼명 (`genome`, `genome_sequence`).
- 스택: Python + FastAPI + SQLite + APScheduler + Jinja2 (승인은 안 받았지만 그대로 진행 중).

## 검증된 것 (2026-10-16, Windows 로컬)
- NCBI Datasets v2 `dataset_report` (442건), `genome/download?include_sequence=GENOME` 동작. HEAD는 405.
- BV-BRC `genome/?eq(taxon_lineage_ids,3052225)` 동작 (445건).
- Henipavirus nipahense taxon_id = 3052225.
- check → dry-run → 승인 → apply (442건) → 재확인(변경 없음), 버전/메타데이터/삭제 감지.
- Pathoplexus(LAPIS): organism 목록 API가 없어 `sources.PPX_KNOWN`의 slug(mpox, ebola-zaire, dengue 등, 실호출로 확인)를 NCBI taxon으로 자동 탐색한다. 필터는 `ncbiVirusTaxId`(없으면 `ncbiVirusName`). Monkeypox(10244) → `mpox` 연결 확인, 10,437건 join. Nipah/Henipavirus는 organism 없음 → 자동 건너뜀. 새 organism이 생기면 `PPX_KNOWN`에 slug 추가.

## 검증된 것 (2026-10-16, Linux, conda 환경 `genomecurator`, Python 3.12)
- `mamba create -n genomecurator -c conda-forge python=3.12 pip` 후 `pip install -e ".[dev]"`. pytest 12건 통과.
- 웹 클릭 검증: 바이러스 추가(조회 약 48초) → 지금 확인 → 전체 선택 → 승인 → Apply로 442건 적용, 유전체 화면 표시 정상.
- 스케줄러 동작 확인(간격 18초로 시험, audit_log에 check 기록). 백업→복원→재확인은 "no changes".
- mpox check: Pathoplexus `mpox` 연결, NCBI 13,462 / Pathoplexus 13,436 / BV-BRC 13,647건.
- 백업 버그 수정: 마이그레이션으로 추가된 열(assembly_method, sequencing_platform, clade 등)이 백업에서 빠지던 문제(`db.migrate`).

## 미완료 / 다음 작업
- Pathoplexus 대상 바이러스 추가 시 `config.PATHOPLEXUS_ORGANISMS` 설정과 필드 매핑 검증
- genome_feature(주석), 재개 가능한 다운로드, 로깅

## 외부 호출 규칙
- 모든 웹 요청(NCBI/BV-BRC/Pathoplexus)은 최소 300ms 간격 (GC_REQUEST_DELAY_MS), 429/5xx는 지수 백오프로 최대 3회 재시도 (sources._get).
- NCBI는 Datasets v2 REST API 사용. datasets CLI는 미사용(미설치). 대용량 서열 일괄 다운로드가 필요하면 CLI 백엔드를 선택 옵션으로 추가 검토.

## 분류군(taxon) 범위 규칙
- NCBI는 레코드를 종(species) 노드(예: `Orthopoxvirus monkeypox` 3431483)나 그 하위 이름(`Monkeypox virus` 10244, `Mpox virus` 등)에 섞어서 붙인다. 하위 노드로만 조회하면 종 노드에 붙은 레코드가 빠진다(예: PX379520, 652건).
- 워치 등록·확인 시 `sources.ncbi_species_scope()`가 rank 없는 노드의 SPECIES 조상을 찾아 워치의 taxon을 종으로 넓힌다(감사 로그 `watch_scope`, changeset 노트에 표시). NCBI 13,462건, BV-BRC는 같은 taxon_id의 lineage 조회.
- Pathoplexus는 taxon으로 organism(slug)만 찾고, 가져올 때는 organism 전체를 받아 INSDC accession으로 join한다(NCBI가 accession 집합을 정함).
- 표시되는 `genome_name`은 소스 원본 이름 그대로라 `Monkeypox virus`/`Orthopoxvirus monkeypox`가 섞여 보인다(수집 범위에는 영향 없음). 통일이 필요하면 종 이름 + strain으로 정규화하는 단계를 추가한다.

## 워치 등록 검증
- `sources.ncbi_resolve_virus()`가 입력 텍스트를 NCBI **바이러스** taxon으로 확인한다. 바이러스가 아니거나(예: `andes`는 속 Andes) 유전체 레코드가 없으면 등록을 거부하고 후보를 보여 준다(예: `Andes virus` → Orthohantavirus andesense).
- 공식 학명으로 저장하고, 종 아래 노드면 종으로 넓힌다. NCBI가 0건이면 check가 `refusing to compute removals`로 중단하며 원인 taxon을 알려 준다.

## 워치 추가 사전 확인 / Pathoplexus 변이
- `genomecurator/lookup.py`: 입력한 이름을 ICTV(EVORA/OLS4 API) · NCBI Virus · BV-BRC · Pathoplexus에서 확인하고 공식 종명을 제안한다(`/api/lookup?q=`, 읽기 전용, 소스별 오류 격리, 요청 간 300ms라 10~20초).
- 유전체 상세 패널: Pathoplexus accession이 있으면 LAPIS(CORS 허용)에서 clade/lineage, SNP·결실·삽입 합계, 아미노산·염기 변이를 바로 불러오고, "페이지 보기"로 pathoplexus.org 서열 페이지를 iframe 서랍에 임베드한다(X-Frame-Options 없음 확인).
- 워치별 Pathoplexus organism은 `watchlist.pathoplexus_organism`에 저장(check 때 갱신, 없으면 한 번 탐색).

## 배포하기 (로컬 → GitHub Pages)
- 상단 **배포하기** 버튼: 미리보기(`/api/deploy/preview`)로 바이러스·유전체 수를 보여 주고 확인을 받은 뒤, `docs/index.html`을 현재 DB로 다시 만들어 그 파일(과 `.nojekyll`)만 커밋·push한다(`genomecurator/deploy.py`). 이전 게시본과 데이터가 같으면 아무것도 하지 않는다.
- git 인증은 기존 자격 증명(Credential Manager)을 쓰며 프롬프트는 끈다. push가 거절되면(원격이 앞섬) 메시지로 알리고 로컬 커밋은 남는다.
- 다른 사이트에서 오는 POST는 Origin 검사로 거절한다.

## 서열(FASTA) 다운로드 최소화
- `curate._needs_sequence()`: 신규(new)는 항상 FASTA를 받고, 변경(changed)은 `accession_version`/`genome_length`/`sequence_hash`가 바뀌었거나 로컬 파일이 없을 때만 받는다. 그 외 메타데이터만 바뀐 항목은 표만 갱신하고 기존 FASTA·md5를 그대로 둔다.
- dry-run은 다운로드가 필요한 건수와 용량만 계산하고 샘플 검증하며(`sequence_plan` 항목에 건수 표시), apply 결과 메시지에 "서열을 다시 받지 않은" 건수를 알려 준다.
- 소스가 더 이상 주지 않는 값은 changed 항목 적용 시 비워진다(그렇지 않으면 같은 변경이 매번 다시 올라온다).

## 메타데이터만 먼저 반영 / 서열 나중에 받기
- 큐레이션 화면의 **서열(FASTA)도 받기** 체크를 끄고 Apply하면 신규 항목은 `genome` 표에만 들어가고 FASTA·`genome_sequence`는 만들지 않는다(`curate.apply(fetch_sequences=False)`). 서열이 바뀌어야 하는 변경 항목은 승인 상태로 남겨 둔다.
- 유전체 화면의 **Get FASTA**(`/genomes/sequences`)가 FASTA가 없는 유전체(선택한 것, 없으면 필터에 보이는 전체)를 받는다(`curate.fetch_missing_sequences`). 100건씩 받아 배치마다 커밋하므로 중단돼도 남은 것만 이어서 받는다. `FASTA` 필터/열(Yes/No)로 상태를 볼 수 있다.

## 대량 선택 (1000건 초과)
- 선택한 항목 id는 `items` 한 필드(쉼표 구분)로 보낸다. 항목마다 `item` 필드를 보내면 Starlette 폼 한도(1000 필드)에 걸려 `Too many fields` 오류가 난다(`web._item_ids`, 기존 `item` 방식도 호환).

## 크로스플랫폼 (Windows · macOS · Linux)
- 데이터 위치: `GC_DATA_DIR` → `<저장소>/data` → OS 사용자 데이터 폴더 순(`config._default_data_dir`). 현재 폴더에 의존하지 않는다.
- 모든 텍스트 파일 쓰기는 UTF-8 + LF(`curate._write_fasta`, deploy, build_static). `.gitattributes`로 저장소는 LF 고정. BOM 제거(테스트로 검사).
- 워치 이름은 `config.safe_name()`으로 폴더명에 안전하게 변환(Windows 금지 문자 대응). git은 PATH → 각 OS 기본 설치 위치 순으로 탐색(`deploy._git`).
- CLI 출력은 UTF-8로 재설정(cp949 등 콘솔 대응), `serve --open`으로 브라우저 자동 열기. 시작 스크립트: `scripts/start.bat | start.sh | start.command`.
- `.github/workflows/ci.yml`: ubuntu/windows/macos × Python 3.10/3.12에서 pytest.

## 워치 편집/삭제
- 워치 목록의 **편집**: 이름, 수집 분류군(NCBI로 바이러스·레코드 확인 후 종 노드로 확장), 자동 확인 켜기/끄기. 분류군을 바꾸면 이미 반영된 유전체가 다음 확인에서 '삭제됨' 후보로 올라올 수 있다고 알려 준다.
- **삭제**: 영향 범위(유전체·FASTA 파일·changeset 수)를 보여 주고 워치 이름을 그대로 입력해야 활성화된다. 유전체, 출처, 이력, FASTA 파일(`data/genomes` 안에 있는 것만), changeset을 함께 지운다(되돌릴 수 없음). 데이터를 남기고 수집만 멈추려면 편집에서 자동 확인을 끈다. 서버 쪽에서도 이름 일치와 Origin을 검사한다.
- FASTA 경로는 DB에 `data` 기준 상대경로(`/` 구분)로 저장한다(`config.rel_data_path/data_path`). 폴더를 옮기거나 OS를 바꿔도 링크가 유지되고, 예전 절대경로(Windows 경로 포함)도 읽힌다.

## Pathoplexus organism 자동 발견
- `sources.pathoplexus_organisms()`가 pathoplexus.org 메뉴에서 organism 목록을 읽는다(예: `ebola-bdbv`, `andv`). 실패하면 `PPX_KNOWN`으로 대체한다.
- 워치의 NCBI taxon이 그 organism 서열의 `ncbiVirusTaxId` 자신이거나 하위 노드이면 연결한다(`pathoplexus_index`, NCBI lineage로 판정). 이름으로도 한 번 더 찾는다.

## 백업 / 복원
- `python -m genomecurator backup <zip>` / `restore <zip> [--force]` (`genomecurator/backup.py`). 포함: watchlist, genome, genome_sequence, field_provenance, genome_version, changeset, audit_log, 결정이 있는 changeset 항목(rejected/hold/applied/approved), FASTA. 제외(재생성 가능): source_record, pending 항목. zip-slip 방지, 실행 중 DB는 잠금 메시지로 안내.
- 데이터(`data/`)와 백업 파일은 GitHub에 올리지 않는다(저장소가 공개). 옮길 때는 백업 zip을 직접 복사해 가져가 `restore`한다.

## Superset 분석 (2026-10-16)
- `scripts/superset.sh init|run|bootstrap`, 설명은 `docs/SUPERSET.md`. 평탄화 뷰 `v_genome_analysis`(앱이 열 때마다 재생성), 읽기 전용 연결, 포트 8089.
- 검증: Superset 6.1.0 설치·기동, bootstrap으로 연결/데이터셋/차트 10/대시보드 생성 성공. 차트 실제 렌더링은 브라우저 로그인 후 확인 필요(미검증).

## 소스 간 정합성 (일치/불일치) 검증
- `check`마다 `reconcile.analyze()`가 NCBI·BV-BRC·Pathoplexus를 accession 단위로 대조해 필드별로 **일치 / 표기 차이(같은 사실, 다르게 표기) / 불일치**로 분류한다. 요약은 `reconcile_report`, 건별 목록은 `source_conflict`(체크마다 재생성, 백업 제외)에 저장한다. 채택 값은 여전히 소스 우선순위가 정한다(대조는 측정·기록만).
- 비교 제외 필드: genome_name, genome_status, release_date, update_date(소스마다 의미가 달라 차이가 정상). 표기 차이로 흡수: 날짜 정밀도/범위/월 이름, 국가 별칭(Zaire 등), `USA: CA`/`USA, CA`, 숙주의 나이·성별 덧붙임, strain 포함 관계, 공백·구두점.
- 웹 **정합성** 탭(`/integrity`): 일치율·불일치 수·NCBI 밖 레코드 요약 → 소스별 보유 현황 → 필드×소스쌍 막대(칸 클릭 시 건별 목록) → 로컬 무결성 검사. CLI: `genomecurator reconcile <이름>`, `genomecurator verify [--no-files]`(오류 시 종료코드 1).
- `verify.py`: watch 없는 유전체, 고아 서열 행, 날짜(미래/형식/공개일 역전), 길이 이상, FASTA 실제 길이·md5 ≠ 기록값, Removed, 출처 누락.
- 실측(mpox, 2026-10-16): NCBI 13,462 / BV-BRC 13,647 / Pathoplexus 13,436. NCBI에 없어 수집 제외되는 레코드 233건(BV-BRC 208 = PDB 체인 등, Pathoplexus 25). 실제 불일치 약 2,000건(strain, genome_length, isolation_source, accession_version 156, 수집일 오타 `2202-07`/`1922-05-30`, 국가 DRC↔South Korea 등). nipah: 불일치 180건.
- 다음 단계 후보: 불일치 해소 규칙(필드별 '더 구체적인 값' 채택), 버전 불일치 항목에서 Pathoplexus 보강 제외, Superset에 `source_conflict` 연결.

## 검토 단계에서의 정합성 적용
- changeset 화면: 상단 경고 배너(불일치 항목 수), **정합성**(일치/표기차이/불일치)·**불일치 필드** 열과 필터. 불일치만 골라 보류/거절할 수 있다.
- 승인·Apply 시 선택 항목 중 불일치가 있으면 확인창에 건수를 보여 준다(차단은 아님, 사람이 판단).
- Dry-run에 `source_agreement` 항목(정보용, READY를 막지 않음)으로 불일치 건수와 예시를 남긴다.
- 정합성 기준은 마지막 check의 `source_conflict`이므로 check 후 만들어진 changeset 기준으로 정확하다.

## UI 일관성 / 상세 패널
- 오른쪽 상세 패널(Genome detail)은 가운데 손잡이를 드래그해 너비를 조절한다(200px~창 너비의 75%, 더블클릭=기본값, localStorage `gcRpW`에 기억).
- 필터(facet) 영역은 그룹당 높이 96px로 줄였다.
- 모든 페이지가 같은 틀을 쓴다: 전체 너비·같은 여백(`main`), 제목+컨트롤 한 줄(`.pagehead`), 공통 KPI 카드(`.cards/.kpi`), 같은 입력·버튼 높이(30px), 표 글자 13px. 이력은 긴 상세를 160자로 줄이고 "전체 보기"로 펼친다.

## 변경 사항(changeset) 화면 쉬운 말로
- 큐레이션 화면 "검토할 변경 사항": 안내문(변경 사항이란? ①확인→②검토·승인→③적용) + 변경 사항마다 카드(새로 발견/정보 바뀜/NCBI에서 사라짐 숫자, 승인·거절·보류·미결정 막대, 소스 간 불일치 경고, 다음 할 일 문장, 검토하기·삭제). 적용 완료분은 아래 한 줄 요약 표.
- 상세 화면: 한국어 상태/제목, "이렇게 쓰세요" 안내와 용어 설명(new/changed/removed, pending/approved/rejected/hold/applied), Dry-run → "적용 전 점검".

## 서열 본문 검증 (seqcheck)
- 메타데이터가 애매할 때 서열 본문(헤더·공백 제거, 대문자, SHA-256)이 같은지로 같은 레코드인지 확인한다. 기준은 NCBI. 대상은 길이/버전/해시가 소스마다 다른 항목(기본) 또는 검토 화면에서 고른 항목. 서열은 20건씩 받아 비교 후 버리고 결과만 `sequence_check`에 저장한다. 정합성 탭 "서열로 검증" 버튼, 검토 화면 Verify seq 버튼.
- 소스 쌍별로 해소한다: 비(非)기준 소스가 NCBI와 서열이 같으면 그 쌍의 길이/버전/해시 불일치는 "해소"로 바뀌고 정합성 열/배너/dry-run 집계에서 빠진다.
- 실측(mpox 404건): Pathoplexus는 392건이 NCBI와 동일(156건은 INSDC 버전 표기만 뒤처짐, 서열은 같음), 12건은 다름. BV-BRC는 404건 모두 NCBI와 다름(길이 다름 401, 길이 같고 내용 다름 3) → 길이 불일치 401건은 실제로 다른 서열.
- 개선 예정: 다운로드 전 dry-run(링크 유효성·예상 용량·길이 선비교로 다운로드 대상 축소).
- **dry-run 먼저**: ① "계획 확인"이 서열 없이 3개 소스에 링크·길이만 조회해 항목별 상태(받아서 비교 / 길이가 달라 다름 확정 / 이미 검증됨 / 소스에 없음 / 오류)와 용량·요청 수·시간을 보여 주고(`sequence_plan`), ② "이 계획대로 실행"이 필요한 것만 받는다. mpox 404건: 받을 서열 3건(1.2MB), 길이 비교로 약 314MB 생략. 이미 검증된 건은 재개 시 건너뛴다.
- NCBI `dataset_report`는 20건씩 페이징되므로 `page_size`를 줘야 한다(안 주면 대부분 "없음"으로 보임).

## 진행 모달
- 오래 걸리는 작업(확인·적용·서열 검증·배포)은 화면 중앙 모달로 보여 준다: 큰 진행 막대와 퍼센트, 경과 시간, 그 아래 시각이 찍히는 로그. 작업 중에는 뒤쪽 조작을 막고(오류일 때만 닫기), 새로고침해도 이어서 보인다.

## 로깅 (2026-10-17)
- `genomecurator/log.py`: `data/logs/genomecurator.log`(5MB×5 회전, UTF-8) + 콘솔(WARNING 이상). 모든 줄에 맥락이 붙는다: `시각 레벨 모듈 [job=작업ID watch=바이러스 changeset=번호] 내용`. `GC_LOG_LEVEL`(파일, 기본 INFO), `GC_LOG_CONSOLE_LEVEL`(콘솔, 기본 WARNING). CLI·웹 시작 시 `log.setup()`.
- 외부 HTTP는 `sources._get`이 매 시도마다 기록한다(URL은 api_key 마스킹, 상태코드, 크기, 소요 시간, 재시도·백오프). DEBUG=성공 호출, WARNING=재시도/4xx, ERROR=최종 실패.
- 단계별 이정표(INFO): check(소스별 건수·소요 시간·대조 결과·changeset), apply(적용/실패/보류, 실패 accession별 WARNING), 서열 검증 계획·판정, 스케줄러 실행, 작업 시작·진행·완료 시간.
- **조용히 삼키던 오류 제거**: Pathoplexus organism 색인·탐색, NCBI 계통/종 범위 조회가 실패하면 `SourceUnavailable`을 던진다. `fetch_all`은 소스가 실패하면 **확인 자체를 중단**한다(일부 소스만으로 병합하면 보강 필드가 비워질 수 있음): changeset·유전체에 아무것도 쓰지 않고 메시지와 로그만 남긴다. 소스가 단순히 비어 있는 경우만 "no records (skipped)" 노트. CLI는 traceback 대신 한 줄 오류.
- 웹 작업(`start_job`): 작업별 로그를 서버가 보관해 모달이 그대로 보여 준다(새로고침해도 유지). 실패 시 traceback은 로그 파일에, 요약은 모달과 `audit_log`(`job_failed`)에 남는다.
- `audit_log`: 상세는 600자로 줄이고, 실패 기록은 `audit_now`로 별도 트랜잭션에 써서 롤백돼도 남는다(`job_failed`, `scheduled_check_failed`, `command_failed`). `decide`는 id 목록 대신 건수와 범위만 저장.
- 이력 화면 아래 **서버 로그** 보기(전체/경고 이상/오류만, 최근 N줄).
- 알려진 점: Pathoplexus가 없는 바이러스(예: nipah)는 organism 탐색에 약 20초가 걸린다(색인 8초 + 이름 탐색 16회). 필요하면 색인 캐시를 디스크에 저장해 줄일 수 있다.

## 스키마 점검 (schema dry-run, 2026-10-17)
- 웹 **스키마** 탭의 "스키마 점검 실행" 버튼 / `python -m genomecurator schema-check`(소스 코드 업데이트가 필요하면 종료코드 1). 읽기 전용, 약 5~10초. 소스마다 실제 응답 샘플 200건만 읽고 서열은 받지 않는다(서열 엔드포인트는 가장 작은 로컬 유전체 1건으로만 확인).
- **코드가 읽는 필드를 측정한다**: 어댑터(`normalize_ncbi`/`normalize_bvbrc`/`normalize_ppx`, 이번에 분리)를 `_Track` 래퍼로 실제 샘플에 돌려 접근한 키 경로를 모은다. 손으로 쓴 목록이 아니라서 코드가 바뀌어도 어긋나지 않는다.
- 판정: **필드가 사라짐**(지난번엔 값이 있었는데 지금 샘플에서 전부 비어 있음, 이름이 바뀐 후보 제시) · **값 형태 변경** · **채움 비율 급감** · **응답 최상위 구조 변경**(예: `data`→`results`) · 접속/응답 실패 · 서열 엔드포인트 실패는 오류/경고. 새로 생긴 필드·읽지 않는 필드(채움 25%↑)는 참고. 첫 실행은 기준(baseline)만 저장한다.
- **로컬 점검**: 코드가 기대하는 테이블·열(마이그레이션 포함) 존재, `GENOME_COLS`/`TRACKED_FIELDS`/`reconcile.FIELDS`가 genome 열에 있는지, 어댑터가 내는 필드가 genome 열에 없어 적용 때 버려지는지, 분석 뷰가 열리는지. 원본 sqlite 연결로 본다(`connect()`는 마이그레이션을 다시 적용해 누락을 숨긴다).
- **기준 보존**: 문제(오류/경고)가 있는 소스의 기준은 갱신하지 않는다(드리프트 상태가 기준이 되어 다음번에 "정상"으로 보이는 것을 막음). 코드를 고쳤거나 영향이 없다고 판단하면 화면의 "현재 상태를 새 기준으로 저장"(`schema_check.accept`).
- 이력: `schema_report`(최근 20건), 기준: `schema_snapshot`. `fetch_pathoplexus`는 요청 필드가 사라지면 조용히 건너뛰던 것을 경고 로그로 남긴다.
- 실측(2026-10-17): 세 소스 모두 정상. NCBI 필드 41개 중 코드가 읽는 27개, BV-BRC 76개 중 24개, Pathoplexus(mpox) 136개 중 22개. 인위적으로 필드 삭제·이름 변경(`length`→`seqLength`)·형태 변경·응답 구조 변경을 주입했을 때 모두 잡혔다.

## 매일 스키마 점검 (바이러스별, 2026-10-17)
- 서버가 켜져 있으면 매일 `GC_SCHEMA_CHECK_HOUR`시(기본 03시)에 `schema_check.run_light()`가 돈다: 등록된 **모든 바이러스**에서 소스마다 `GC_SCHEMA_CHECK_SAMPLE`건(기본 3건)을 읽는다. 수동: 스키마 탭 "지금 바이러스별 점검 실행" 또는 `schema-check --daily`. 기준(baseline)이 없으면 전체 점검을 먼저 돌려 기준을 만든다. 기준은 바꾸지 않는다.
- **3건으로 증명할 수 있는 것만 판정**한다(오탐 방지): 바이러스별로는 응답 구조 변경·값 형태 변경·접속 실패(일시 장애일 수 있어 경고)·어댑터 예외, **"필드가 사라짐"은 모든 바이러스를 합쳐서** 지난 기준에서 거의 항상(≥90%) 채워지던 필드가 전부 비었을 때만. 바이러스 하나의 3건만으로는 판정하지 않는다(실제로 mpox의 오래된 레코드 3건에 `host`가 없어서 오탐이 났고, 이 규칙으로 고침).
- 전체 점검의 기준도 이제 등록된 모든 바이러스의 샘플(바이러스마다 200건, Pathoplexus는 organism별)을 합쳐서 만든다. 한 바이러스의 특성이 기준이 되지 않는다.
- 문제가 있으면 상단 **스키마** 탭에 `!` 배지(오류=빨강, 경고=주황)가 뜬다(1분 캐시). 결과는 스키마 탭 "매일 점검" 표(바이러스×소스)와 문제 목록, 이력(`schema_daily`, 실패는 `schema_daily_failed`)과 로그에 남는다.
- 한계: NCBI·Pathoplexus 샘플은 "처음 N건"(오래된 레코드)이다. NCBI는 정렬 파라미터가 먹지 않았고 LAPIS의 내림차순 지정은 확인하지 못했다(BV-BRC만 최신순 정렬 가능). 새 필드가 최신 레코드에만 생기는 경우는 매일 점검이 늦게 발견할 수 있어 전체 점검(수동)을 병행한다.
- 실측: 바이러스 2종(nipah, mpox)에서 오탐 없이 통과, 주입한 구조 변경·형태 변경·전체 소멸 필드는 모두 감지.

## 환경 파일 / README (2026-10-17)
- `environment.yml`(앱), `superset/environment.yml`(Superset) 추가. 새 환경을 만들어 다른 폴더에서 임포트·pytest 24건 통과로 검증했다(`-e .`는 yml이 있는 폴더 기준이라 저장소 루트의 파일로 만들어야 한다). pip 단계에서 `--no-deps` 줄은 mamba가 패키지 이름으로 읽어 실패하므로 쓰지 않는다.
- CI에 `conda-env` 작업 추가(environment.yml이 실제로 빌드되고 테스트가 통과하는지).
- `.claude/launch.json`은 `mamba run -n genomecurator ...`, 포트 8001. 이 컴퓨터 전용 절대 경로를 저장소에 넣지 않는다.
- README를 새로 썼다(철칙·소스 삼형제·정합성 3종·명령어·환경변수·FAQ·한계·코드 지도).

## README 시연 GIF (2026-10-17)
- `scripts/demo/`: `recorder.py`(헤드리스 Chrome을 DevTools로 조종, 커서·자막 오버레이, 프레임 수집, GIF 조립), `make_gifs.py`(시나리오 4개), `demo_server.py`(시연 전용: 스키마 화면에서 "소스가 바뀐 척"하는 서버, 실제 앱에는 쓰이지 않음). 의존성은 `pip install -e ".[demo]"`.
- 결과: `docs/demo/01-quickstart.gif`(빈 폴더에서 nipah 등록→확인→승인→적용, 실제 인터넷), `02-integrity.gif`(mpox), `03-schema.gif`(앞쪽 실제, 뒤쪽은 시연용 가정 — README에 명시), `04-explore.gif`. 합계 약 4MB.
- 시연은 임시 복사본에서 돌아 실제 `data/`를 건드리지 않는다. 녹화 도구는 Chrome 연결이 실패해도 Chrome을 남기지 않고 빈 포트를 골라 쓴다(이전 실행의 Chrome에 잘못 붙는 문제를 겪고 고침).
- 같은 날 사용자 요청으로 실제 `data/`에서 nipah(유전체 442·FASTA 90·changeset 1)를 삭제했다. 삭제 전 백업: `~/genomecurator-backup-before-nipah-delete-20261017-1008.zip`(`restore --force`로 복구 가능). mpox 13,462건은 그대로.

## 버전 관리 / 변경 이력 / 라이브러리 (2026-10-17, v0.1.0)
- **0.1.0이 첫 릴리스**다(지금까지 만든 것 전체). 이후 `0.1.1`(버그), `0.2.0`(기능·동작 변화)처럼 올린다. 규칙과 절차: `docs/RELEASING.md`(한·영).
- 버전은 `genomecurator/__init__.py`의 `__version__` 한 곳(pyproject는 `dynamic`으로 읽음). `DB_SCHEMA_VERSION`(=1)은 앱 버전과 별개로, 오래된 코드가 못 읽는 DB 구조 변경일 때만 올린다.
- 기록: DB `app_meta`(first/current version, schema version) + 버전이 바뀔 때 이력 `version`, 각 changeset `summary.app_version`, 백업 `manifest.json`(`app_version`, `db_schema_version`; 더 새로운 DB 구조의 백업은 복원 거부). 노출: `--version`, 화면 하단, `/api/version`, 로그 첫 줄.
- `CHANGELOG.md`(영문)·`CHANGELOG.ko.md`(한글), Keep a Changelog 형식(Added/Changed/Fixed/Data & compatibility/Dependencies). 테스트 `test_changelog_matches_the_code_version`이 코드 버전·날짜·순서·링크·영한 일치를 강제한다.
- `docs/DEPENDENCIES.md`는 `scripts/dependencies.py`가 설치된 패키지에서 버전·라이선스를 측정해 만든다(역할은 수기). `test_every_declared_dependency_is_documented`가 pyproject·superset requirements의 모든 라이브러리가 문서에 있는지 검사한다.
- 태그 `v0.1.0`. 태그를 push하면 그 커밋까지의 모든 커밋이 올라가므로, 워크플로 파일을 바꾼 커밋(`3223f0e`)이 앞에 섞이면 `workflow` 권한이 필요하다 — 그래서 CI 커밋은 릴리스 커밋 뒤에 둔다.

## 스키마 탭: 테이블 구조도 (2026-10-17, Unreleased)
- 스키마 탭 상단에 구조도: 소스 3곳의 필드(NCBI는 중첩 객체를 그룹으로) → 내 `genome` 열을 선으로 잇고 상태를 색으로 칠한다(`schema_check.build_board`, 리포트 JSON의 `board`).
- 필드→열 연결은 **측정값**이다(`schema_check.field_map`): 어댑터(`normalize_*`)를 필드 하나만 든 행에 돌려 바뀌는 열을 읽는다(탐침 값 3개: 임의 문자열 / `RefSeq` / `COMPLETE`로 값 비교 분기까지). 시험 행에는 어댑터가 버리지 않도록 기본 accession을 넣고, 목록 안 형제 필드는 채워 둔다. `bvbrc_genome_id`가 측정에 잡히도록 `normalize_bvbrc`로 옮겼다(동작 동일).
- 상태: 소스 필드 = 정상/경고/오류/새 필드/샘플에 값 없음/읽지 않음. 열 = 채우는 필드가 모두 오류면 오류, 일부만 문제면 경고(다른 소스가 대신 채울 수 있음), 로컬 DB에 열이 없으면 오류, 소스에서 안 오면 "앱이 채움". 부모 객체가 사라지면 자식도 같이 표시한다.
- `_compare`/`field_map`은 이번 실행의 읽기 기록 ∪ 지난 기준의 읽기 기록을 쓴다(읽기 기록은 데이터에 따라 달라짐). 채움 비율은 레코드당 한 번만 센다(이전에는 목록 원소마다 세어 500%가 나왔다).
- 시연용 가짜 드리프트는 `scripts/demo/demo_server.py`(깃 추적, 실제 앱 영향 없음). `03-schema.gif`를 구조도 중심으로 다시 녹화했다.
- 아직 릴리스하지 않았다(`CHANGELOG`의 [Unreleased]). 다음 릴리스는 새 기능이므로 0.2.0이 맞다.
- **간단히 보기**(스키마 탭 기본, 탭으로 구조도와 전환, 선택은 `localStorage gcSchemaView`): 열 × 소스 표(`schema.html`의 `#v-simple`), 확인이 필요한 열을 위에 모으고, 열 묶음은 `schema_check.COLUMN_GROUPS`(화면용 분류일 뿐 연결은 측정값). 구조도는 `#v-board`에 그대로 두었다. `03-schema.gif`는 간단히 → 구조도 → 드리프트 순서로 다시 녹화.

## 범위·이름 결정 (2026-10-17)
- **바이러스 전용으로 확정.** 세균(어셈블리 단위·contig 여러 개·수 Mb 서열·종당 수십만 건)은 레코드 단위와 정합성 검증 방식이 달라서 이 프로젝트에 넣지 않고 **새 프로젝트**로 만든다. 세균 설계 때 참고: NCBI Datasets 일반 유전체 엔드포인트, BV-BRC(MLST·AMR 강점), 어셈블리 수준·RefSeq 필터, 서열은 요청 시에만 받기.
- 표시 이름: 정식 `Virus Genome Metadata Curator`, 짧게 `Virus Genome Curator`(헤더·푸터·제목). **코드 이름은 그대로**(패키지 `genomecurator`, 명령어, 저장소 `GenomeCurator`, 데이터 폴더 `GenomeCurator`, `GC_*`, 백업 식별자 `genomecurator-backup`, Superset 대시보드 제목 `GenomeCurator 메타데이터 분석`). 코드까지 바꾸려면 데이터 폴더 이전·옛 명령어 별칭·백업 식별자 호환이 필요해 별도 릴리스로 미룬다(0.2.0 이후 검토).
- 공개 페이지(`docs/index.html`)의 제목은 다음 "배포하기" 때 새 이름으로 바뀐다.

## Near complete (준완전) 완성도 (2026-10-17, Unreleased)
- `completeness.py`: Complete(NCBI가 Complete — 신뢰) / Near complete(Complete 아님 + 길이 ≥ 기준 게놈의 `GC_NEAR_COMPLETE_PCT`% , 기본 90) / Partial. 기준 게놈 = 그 바이러스의 수집된 RefSeq 중 가장 긴 것(없으면 Near complete는 판정 안 함). 같은 규칙을 분석 뷰 SQL(`completeness`, `ref_coverage_pct`, `ref_length`)에도 썼고 임계값은 `app_meta.near_complete_pct`에서 읽는다. 파이썬과 SQL이 13,462건 전부 일치하는지, 임계값 변경 시에도 일치하는지 테스트한다.
- 90%의 근거: 실데이터 분포에 골짜기가 있다(mpox 50~90% 구간 33건/13,462, nipah 2건/442). CheckV high-quality(≥90%)와 같다. 더 엄격히(95%) 하려면 환경변수.
- 한계: ① 길이 비율일 뿐 정렬 기반 커버리지·N 비율은 아니다(FASTA가 로컬에 있을 때 N ≤ 5%를 추가할 수 있음). ② 기준이 하나라 계통별로 길이가 다른 종에서는 Complete가 아닌 짧은 계통이 Partial로 보일 수 있다(mpox clade I은 기준 대비 약 95~97%, 그래도 90% 이상이라 Near complete). ③ 사용자가 기준 길이를 직접 지정하는 화면은 아직 없다.
- 실측: mpox에서 NCBI Partial 8,223건 중 7,583건이 Near complete, 640건만 Partial.

## 소스 값 중 고르기 (resolve, 2026-10-17, Unreleased)
- `resolve.py`: `analyse`(필드별로 소스 값이 다른 곳 + 표기차이/불일치), `suggest`(추천 규칙·신뢰도·이유), `save`(선택 저장 + 열린 changeset 항목의 `new_record_json`/`provenance_json`/diff/fingerprint 갱신), `preview_auto`/`apply_auto`(Auto-pick). 선택은 `field_choice`(watch, accession, field → source|custom, value)에 저장되고 `curate.merge(per_source, choices)`가 이후 모든 check에서 존중한다(고른 소스의 **현재** 값을 쓰고, 그 소스에 값이 없어지면 우선순위 병합으로 되돌아감). 백업 대상 테이블에 포함. `source_agreement`는 선택한 필드를 열린 불일치에서 뺀다(정합성 열 "선택 완료").
- 추천 규칙: ① 불가능한 날짜(1900 이전·미래) 제외 ② 길이·버전·해시는 NCBI ③ 같은 사실의 다른 표기는 가장 잘 쓴 값(날짜: 구체성 > ISO > 길이, 지역/포함관계: 가장 자세한 것, 호스트: 학명, 국가: 표준명) = 높음(단, 한쪽이 다른 쪽을 포함할 뿐이면 보통) ④ 소스 다수가 같은 값 = 보통 ⑤ 다수 없음 = 낮음(우선순위로 제안, 직접 확인). Auto-pick 기본은 보통 이상.
- 실측(nipah, 소스 2개): 충돌 182항목 중 높음 11 · 보통 18 · 낮음 178 필드. 낮음은 `oronasopharynx` ↔ `throat swab`(어휘가 다름), `Pteropus medius` ↔ `Pteropus giganteus`(이명인지 다른 종인지 알 수 없음) 등이라 일부러 사람에게 남긴다. 근거 없이 일괄 적용하는 것이 더 위험하다고 판단.
- 검증: 선택→항목 반영→되돌리기→직접 입력→잘못된 입력 거부→일괄 적용→**새 check 후에도 29개 선택이 모두 유지**(0개 유실). 화면에서 라디오 클릭+저장으로 행이 "선택 완료"가 되는 것까지 확인.
- 고친 문제: 일괄 적용이 쓰기 트랜잭션을 연 채 서열 검증 결과를 읽으려고 다른 연결을 열어 "database is locked"가 났다(WAL에서 읽기 스냅샷이 낡은 연결은 쓰기로 올라갈 수 없음). 쓰기 전에 읽도록 순서를 바꿨다 — 중첩 `connect()`는 쓰기 트랜잭션 안에서 부르지 말 것.
- 한계: 필드별 도메인 지식은 최소만 넣었다(학명 동의어·어휘 매핑은 없음). 소스가 둘뿐인 바이러스는 대부분 "낮음"이다. 소스 3곳(Pathoplexus가 있는 바이러스)이면 다수결이 가능해 "보통"이 늘어난다.

## 릴리스 v0.2.0 (2026-10-17)
- 담은 것: Near complete 등급, 소스 값 중 고르기(추천·신뢰도·Auto-pick, `field_choice`), 스키마 탭의 간단히 보기와 테이블 구조도, 바이러스 전용 범위와 표시 이름(Virus Genome Metadata Curator). 자세한 내용은 `CHANGELOG`.
- DB 스키마 버전은 1 그대로(`field_choice`는 추가형). 의존성 변경 없음.
- 절차는 `docs/RELEASING.md`대로: CHANGELOG(한·영) 정리 → `__version__` → 테스트 → 릴리스 커밋 → 주석 태그 `v0.2.0` → push(커밋+태그). 워크플로 파일을 바꾼 커밋(CI)은 릴리스 커밋 **뒤**에 두어야 `workflow` 권한 없이 태그를 올릴 수 있다.
