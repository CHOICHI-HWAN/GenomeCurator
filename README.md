# GenomeCurator

로컬에서 실행하는 바이러스 유전체 수집·큐레이션 도구. NCBI Virus(서열·accession 기준)를 받고,
Pathoplexus > NCBI Virus > BV-BRC 우선순위로 메타데이터를 병합해 BV-BRC 스타일 테이블로 저장합니다.
새 데이터는 **changeset**으로만 쌓이고, 웹 페이지에서 승인해야 로컬에 반영됩니다.

## 실행
```bash
pip install -e .
python -m genomecurator watch-add "Henipavirus nipahense"
python -m genomecurator check "Henipavirus nipahense"   # changeset 생성 + dry-run
python -m genomecurator serve                           # http://127.0.0.1:8000
```
웹(큐레이션 페이지)에서 항목 승인/거절/보류 → Apply. 데이터는 `data/`(변경: `GC_DATA_DIR`).
환경변수: `NCBI_API_KEY`, `GC_CHECK_INTERVAL_HOURS`(기본 24), `GC_DRYRUN_SAMPLE`(기본 3).

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

GitHub Pages로 공개하려면 저장소 Settings → Pages → Source를 `main` 브랜치의 `/docs`로 지정한다. 차트는 CDN(Chart.js)을 쓰므로 보는 쪽에 인터넷이 필요하다.
