# GenomeCurator 핸드오프 (새 세션용 컨텍스트)

## 목적
로컬에서 실행하는 웹앱 기반 바이러스 유전체 수집·큐레이션 도구.
외부 DB의 새 데이터를 주기적으로 확인하고, **사용자가 큐레이션 페이지에서 승인해야만** 로컬에 반영한다.
(무엇이 바뀌었는지 사람이 확인하는 것이 핵심 요구사항)

## 확정된 결정
- 첫 대상: **Henipavirus nipahense**. 바이러스는 이후 계속 추가 → watchlist(taxon + 필터) 구조.
- 데이터 소스 3개를 join해서 출력. **우선순위: Pathoplexus > NCBI Virus > BV-BRC**
  (같은 필드가 충돌하면 우선순위가 높은 소스 값을 채택하고, 출처를 필드별로 기록).
  - 서열/accession의 기준 수집은 NCBI Virus. Pathoplexus·BV-BRC는 메타데이터 보강.
  - Pathoplexus에 해당 바이러스가 없으면 그 소스는 건너뛴다(오류 아님).
- **변경 판정 기준**: (a) 메타데이터 값 변경, (b) GenBank accession의 `.version` 증가. 둘 다 changeset에 올린다.
- 내부 포맷: BV-BRC 스타일 테이블/컬럼명(genome, genome_sequence, 이후 genome_feature).
- 실행 환경: **사용자 로컬 PC** (로컬에서 실행 가능해야 함, 설치 쉬워야 함). 현재 세션은 외부이며 GitHub로 연결.
- 승인 전에는 로컬 데이터를 절대 변경하지 않는다(staging → 승인 → 트랜잭션 Apply → 이력 보존).

## 아직 미확정 (사용자에게 물어볼 것)
- 기술 스택 최종 확인: 제안은 Python + FastAPI + SQLite(WAL) + APScheduler + Jinja2/HTMX. (명시적 승인은 아직 없음)
- 단일 사용자 가정(인증 없음)으로 갈지.

## 업데이트 흐름
Watchlist → 스케줄러 → Check + **Dry-run** → Changeset(staging) → 큐레이션 페이지(건별/일괄 승인·거절·보류, 필드 diff) → Apply → 이력/changelog.
- Dry-run: 다운로드 엔드포인트 응답, accession 존재/삭제 여부, 예상 건수·용량·디스크 여유, 샘플 N건 파싱·md5 검증.
  통과 시 READY, 실패 시 BLOCKED+사유. READY일 때만 사용자에게 알림.
- 거절한 항목은 같은 내용으로 재등장하지 않게 기억.

## 데이터 모델 초안
- `genome`: BV-BRC 필드명 (genome_id, genome_name, taxon_id, taxon_lineage_*, genome_status, genome_length, genbank_accessions, refseq_accessions, segment, strain, subtype, genotype, host_name, isolation_country, collection_date, geographic_location, sequence_md5, date_inserted, date_modified …)
- `genome_sequence`: accession, **version**, length, gc_content, md5, 파일경로
- `source_record`: 소스별(pathoplexus/ncbi/bvbrc) 원본 응답 JSON + 수집시각 (join·diff의 근거)
- `field_provenance`: 필드별 채택 소스
- `genome_version`, `changeset`, `changeset_item`, `audit_log`
- 디스크: `data/genomes/<taxon>/…fna`, `data/staging/<changeset_id>/`, `data/snapshots/…`

## 소스별 접근 메모 (모두 **미검증** — 0단계에서 실호출로 확인할 것)
- NCBI: Datasets v2 Virus API (`/virus/taxon/{taxon}/dataset_report`, `/genome/download`). API key 선택(3→10 req/s). Henipavirus nipahense의 taxid는 직접 조회해 확정(추정 3052225, 미확인).
- BV-BRC: Data API (`https://www.bv-brc.org/api/genome/`, RQL 쿼리). 필드명이 내부 포맷의 기준.
- Pathoplexus: Loculus 기반, LAPIS API (`https://lapis.pathoplexus.org/<organism>/sample/details` 등). 지원 organism 목록에 Nipah/Henipavirus가 있는지 먼저 확인할 것. 없으면 어댑터는 비활성.
- join 키: INSDC accession(버전 제외한 base accession)으로 매칭. 소스별 accession 필드 이름 확인 필요.

## 환경 주의
- 클라우드 샌드박스는 `api.ncbi.nlm.nih.gov`, `bv-brc.org`에 403(네트워크 정책)이라 실호출 검증 불가.
  → 외부 호출은 인터페이스 뒤에 두고 **녹화 fixture**로 테스트. 실호출 검증은 사용자 로컬에서 하거나 정책에서 호스트 허용 후 진행.
- 저장소는 현재 README, LICENSE만 있음. 개발 브랜치: `claude/affectionate-einstein-419mes`.

## 구현 단계
0. API 실호출로 스키마 확정 + fixture 저장 (`docs/api-notes.md`) — 로컬/허용 환경 필요
1. 골격 + DB 스키마 + 3개 소스 클라이언트 + CLI(`check`, `dry-run`) + join/우선순위 병합
2. Dry-run 엔진 + changeset 모델 + 변경 판정(메타데이터 diff, accession.version)
3. 웹 UI: watchlist, 큐레이션 페이지, 승인/Apply, 이력
4. 스케줄러, 알림, 설정(API key, 주기, 저장 경로), 로컬 실행 패키징(`pip install -e .` 또는 `uv run`, Docker 선택)
5. 안정화: 재개 다운로드, 레이트리밋, 로깅, 테스트, BV-BRC형 TSV/FASTA export, feature 주석

## 다음 세션 시작 지침
1. 이 문서를 읽고 스택 확정 여부만 사용자에게 짧게 확인.
2. 1단계부터 구현(소스 어댑터는 fixture 기반 테스트 우선). 0단계 실호출은 사용자 로컬에서 돌릴 수 있는 스크립트(`scripts/probe_sources.py`)로 제공해 결과를 fixture로 가져오게 하면 효율적.
