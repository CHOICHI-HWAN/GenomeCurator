# GenomeCurator 진행 상황 / 이어서 작업하기

임시 환경에서 작업하므로 이 문서와 저장소만 있으면 어디서든 이어갈 수 있습니다.

## 실행
```bash
pip install -e ".[dev]"
python -m genomecurator watch-add "Henipavirus nipahense"
python -m genomecurator serve        # http://127.0.0.1:8000
pytest
```
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

## 미완료 / 다음 작업
- 웹에서 `전체 승인 + 적용` 실제 클릭 검증 (코드 경로는 CLI로 검증됨)
- 스케줄러(`GC_CHECK_INTERVAL_HOURS`) 실동작 확인
- 웹 UI 화면 스크린샷 확인 (스프레드시트형으로 변경 후 미확인)
- Pathoplexus 대상 바이러스 추가 시 `config.PATHOPLEXUS_ORGANISMS` 설정과 필드 매핑 검증
- genome_feature(주석), 재개 가능한 다운로드, 로깅

## �ܺ� ȣ�� ��Ģ
- ��� �� ��û(NCBI/BV-BRC/Pathoplexus)�� �ּ� 300ms ���� (GC_REQUEST_DELAY_MS), 429/5xx�� ���� ������� �ִ� 3ȸ ��õ� (sources._get).
- NCBI�� Datasets v2 REST API ���. datasets CLI�� �̻��(�̼�ġ). ��뷮 ���� �ϰ� �ٿ�ε尡 �ʿ��ϸ� CLI �鿣�带 ���� �ɼ����� �߰� ����.

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
- 2026-10-16 이 PC의 data(726MB)에서 만든 백업(21.7MB)은 저장소의 `data-backup` 브랜치에 올려 두었다(`genomecurator-backup.zip`). 복원: 새 환경에서 `git checkout data-backup -- genomecurator-backup.zip` 후 `restore`.
