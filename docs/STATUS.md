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
