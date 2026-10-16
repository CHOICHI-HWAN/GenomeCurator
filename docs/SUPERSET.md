# Apache Superset으로 메타데이터 다차원 분석

GenomeCurator DB를 **읽기 전용**으로 연결해 Superset에서 바이러스·숙주·국가·수집 연도·clade 등 여러 축으로 분석한다.
Superset은 별도 conda 환경(`gc-superset`, Python 3.11)에서 돌고, GenomeCurator 환경과 섞이지 않는다.

## 구성
- `v_genome_analysis` (SQLite 뷰, `db.ANALYSIS_VIEW`): genome + watchlist + FASTA 유무를 한 줄로 평탄화.
  파생 열: `virus, host, country, collection_year/month, release_year, length_bin, has_sra, in_pathoplexus, has_fasta, n`.
  앱이 DB를 열 때마다 다시 만들어지므로 열 정의를 바꾸면 자동 반영된다.
- Superset 연결은 `mode=ro` URI라서 분석 중 큐레이션 데이터가 바뀔 수 없다.
- `superset/superset_config.py`: 메타DB와 비밀 키는 `superset/home/`(git 제외)에 저장.

## 처음 한 번
```bash
mamba create -y -n gc-superset -c conda-forge python=3.11 pip
mamba run -n gc-superset pip install -r superset/requirements.txt
python -m genomecurator serve --help >/dev/null   # (선택) 앱을 한 번 열어 뷰 생성; 이미 쓰던 DB면 불필요
scripts/superset.sh init        # 메타DB, admin 계정 (비밀번호: superset/home/admin_password.txt)
```
## 실행
```bash
scripts/superset.sh run         # http://127.0.0.1:8089  (SUPERSET_PORT로 변경)
scripts/superset.sh bootstrap   # 연결·데이터셋·차트 10개·대시보드 생성/갱신 (다시 실행해도 안전)
```
로그인: 사용자 `admin`, 비밀번호는 `superset/home/admin_password.txt`. 대시보드 "GenomeCurator 메타데이터 분석".
다른 DB 위치는 `GC_DB=/path/genomecurator.db scripts/superset.sh bootstrap`.

## 알아둘 점
- 포트 8000/8088은 이 머신의 다른 프로그램이 쓰고 있어 8089를 기본으로 했다.
- 데이터가 갱신(check/apply)돼도 Superset은 같은 파일을 읽으므로 대시보드 새로고침만 하면 된다. 캐시 때문에 안 바뀌면 차트의 "강제 새로고침".
- 차트는 API로 만들어 저장된 query context가 없다. 대시보드에서 렌더링되지만, 차트를 편집 후 저장하면 채워진다.
- 앱을 실행하지 않은 새 DB에는 뷰가 없으니 앱(`serve`/`check`)을 한 번 실행한다.
