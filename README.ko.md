# Virus Genome Metadata Curator

[![tests](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml/badge.svg)](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.10%20%7C%203.12-blue)
![version](https://img.shields.io/github/v/tag/CHOICHI-HWAN/GenomeCurator?label=version&sort=semver)
![license](https://img.shields.io/badge/license-GPL--3.0-blue)

[English](README.md)

[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/), [BV-BRC](https://www.bv-brc.org/), [Pathoplexus](https://pathoplexus.org/)에서 **바이러스** 유전체 메타데이터를 모아 같은 유전체에 대한 세 소스의 값을 비교하고, 사람이 승인한 것만 내 데이터베이스에 넣는 로컬 웹 도구입니다. 패키지·명령어 이름은 `genomecurator`입니다.

바이러스 전용입니다. 세균은 레코드 구조(여러 contig의 어셈블리)와 서열 크기가 달라 별도 프로젝트로 만들 예정입니다.

## 주요 기능

### 큐레이션: 모으고, 비교하고, 승인하기

1. **승인 전에는 아무것도 바뀌지 않습니다.** 새로 생기거나 바뀐 레코드는 *changeset*으로 쌓이고, 승인해서 적용한 것만 데이터가 됩니다. 거절한 것은 다시 올라오지 않습니다.
2. **세 소스, 하나의 레코드.** NCBI Virus가 어떤 accession이 있는지 정하고, 값은 Pathoplexus > NCBI > BV-BRC 순으로 채택하며, 값마다 출처가 기록됩니다.
3. **소스끼리의 일치를 수치로 봅니다.** 필드별로 일치 / 표기 차이(`2022-07` vs `2022-07-11`) / 진짜 불일치.
4. **서열 검증.** 애매하면 서열 본문(헤더 제거, SHA-256)을 비교하며, 큰 다운로드를 피하려고 dry-run을 먼저 합니다.
5. **소스의 값 중에서 고르기.** 추천값·이유·신뢰도가 붙고, 고른 값은 이후 확인에서도 유지됩니다.
5b. **우리 용어 사전.** 숙주·나라·수집일·clade를 소스끼리 비교하기 전에 앱의 사전으로 표준화합니다(`Homosapiens` → `Homo sapiens`, `Zaire` → `Democratic Republic of the Congo`, 날짜 구간). 모르거나 의심스러운 값은 그대로 두고 표시하며 원본 값은 항상 기록합니다.
6. **스키마 감시.** 어댑터가 읽는 필드의 이름이 바뀌거나 사라졌는지 매일, 그리고 원할 때 점검합니다.
7. **유전체 완성도 등급**: Complete / Near complete(가장 긴 RefSeq 게놈의 90% 이상, 조절 가능) / Partial.
8. **작업 탭.** 오래 걸리는 작업은 서버에서 돌고 상태·진행률·로그·오류가 남으며, 다른 메뉴로 가거나 서버를 재시작해도 기록이 유지됩니다.

### 품질 · 비공개 유전체 · 데이터셋

9. **서열 QC.** 모든 유전체의 서열 통계(길이, N·모호 염기, GC)와, Nextclade 데이터셋이 있는 바이러스(mpox, 에볼라 종들, 홍역, 뎅기 등)는 [Nextclade](https://github.com/nextstrain/nextclade)의 clade/lineage와 품질 점수. 소스와 다른 clade는 덮어쓰지 않고 검토 화면으로 보내며, 데이터셋은 버튼으로 업데이트합니다.
10. **다른 프로그램으로 내보내기.** 탭·CLI·내부 API에서 정렬된 `metadata.tsv` + `sequences.fasta`(체크섬이 든 manifest 포함)를 만들며, [Nextstrain mpox](https://github.com/nextstrain/mpox) 워크플로의 입력 열에 맞춘 형식이 있습니다.
11. **비공개 저장소.** 공개되지 않은 로컬 유전체를 하나씩 또는 스프레드시트(템플릿 제공)로 올려(올릴 때마다 버전 기록, dry-run 먼저) 공개 세트와 분리된 저장소에 두고, Nextstrain 워크플로 자체의 `private_metadata` / `private_sequences`에 연결되는 *별도 파일*로 내보내므로, 둘을 섞지 않고 최종 분석용 데이터세트를 만들 수 있습니다.
12. **삭제하지 않고 데이터셋에서 제외.** 최종 데이터셋에서 뺄 유전체를 사유와 함께 표시(되돌릴 수 있음)하면 내보내기가 건너뛰고 `excluded.tsv`에 기록합니다.
13. **데이터셋.** 이름 붙인 임의 선별(예: clade Ib만)을 저장하고 실시간 미리보기와 추적 가능한 내보내기(manifest에 규칙·멤버 해시)를 합니다.
    데이터셋마다 유전체 표처럼 보고 대시보드에서도 볼 수 있습니다.

### 데이터셋 분석 (공개+비공개 유전체를 함께)

모두 저장한 *데이터셋*을 대상으로 하고, 분석마다 무엇이 들어갔는지 함께 보관하며, 같은 id인데 내용이 바뀌면(id 목록만이 아니라 내용 지문으로) 알려 줍니다. 외부 도구는 각자의 conda 환경에서 별도 프로그램으로 실행됩니다.

- **분석 허브.** 데이터셋을 고르면 모든 분석 실행(과 오늘의 데이터와 같은지)을 보고 결과 화면을 한곳에서 엽니다. 메뉴는 큐레이션 / 유전체 / 데이터셋 / 분석 / 제출·공유 / 시스템으로 묶여 있습니다.
- **Nextstrain.** 데이터셋으로 Nextstrain mpox 워크플로를 돌리고(설정 overlay와 데이터셋별 YAML 설정은 도구가 처리), 모든 분석을 보관하며, 오늘의 데이터셋과 비교하고 로컬 Auspice로 엽니다.
- **Squirrel.** 정렬 · 서열 QC · IQ-TREE 계통수와 HTML 보고서; QC가 제안한 서열을 계통수에서 빼는 2단계 실행.
- **변이 프로파일.** 빌드 참조 서열(NC_063383.1, clade I 빌드는 Nextstrain의 KJ642613.1) 대비 모든 SNV에 **유전자 · 코돈 · 아미노산 변화**와 APOBEC3형 표시를 붙입니다([mutation_profile](https://github.com/insapathogenomics/mutation_profile)). 개수를 왜곡하는 양 끝 반복 영역(ITR)은 마스킹하고, 유전체의 N 비율 상한을 둘 수 있습니다.
- **변이 그림**([snipit](https://github.com/aineniamh/snipit) 스타일): 유전체마다 한 줄, lineage·clade 띠, APOBEC3형 점, 유전자 트랙, lineage 정의 변이는 굵게.
- **전파 타임라인과 지도.** lineage 또는 lineage 정의 변이([lineage-designation](https://github.com/mpxv-lineages/lineage-designation))를 나라 × ISO 주(주가 너무 얇으면 월)로, 겹치지 않는 첫 관찰 표시, 유전체 수 라인차트, **세계 지도 애니메이션**(오프라인, Natural Earth 경계).
- **변이 분석.** APOBEC3형 변이는 전달되는가? 변화 중 APOBEC3와 배경 변이는 얼마나 되는가(독립 사건, 우연 대비 농축)? 신규 변이와 lineage 시그니처.
- **유행 분석**([Lee, Choi 외 2026](https://doi.org/10.1080/22221751.2026.2746255)의 접근): 트리의 서브클러스터와 정의 변이, 클러스터 크기 과분산(포아송 대 음이항), APOBEC3 시간 누적, 독립 사건 기준 유전자 hotspot.
- **신규 변이 탐색.** 어떤 지정 lineage의 정의에도 없는 변이를 공유하는 계통군을 상태(신규 / 성장 중 / 정체)와 발견 곡선, 첫 관찰 산점도로 보이고, **변이 묶음(PCA)**은 함께 다니는 변이를 보여 줍니다.
- **Pathoplexus 제출 준비.** 비공개 유전체를 생물별 템플릿에 맞춰 점검하고 파일을 만듭니다. **앱은 아무것도 업로드하지 않고 인증 정보도 다루지 않습니다.**
- **위 분석 전체의 타당성 검토**가 [docs/VALIDITY.md](docs/VALIDITY.md)에 있습니다: 믿을 만한 것(SNV 호출과 lineage 정의의 일치, APOBEC3형 변이의 전달), 탐색용인 것, 이 데이터로는 알 수 없는 것.

### 공유

- **Apache Superset 분석**(선택, 읽기 전용)과 GitHub Pages용 **읽기 전용 정적 페이지**.

## 화면

![빈 화면에서 승인·반영까지](docs/demo/01-quickstart.gif)

*Henipavirus nipahense 등록 → 소스 3곳 확인 → 442건 승인 → 적용 (실제 데이터, 속도를 높인 화면).*

| | |
|---|---|
| ![정합성 탭](docs/demo/02-integrity.gif) | ![값 고르기](docs/demo/05-resolve.gif) |
| 정합성: 필드별 일치, 진짜 불일치, 서열 검증 계획 | 검토: 소스별 값, 추천, Auto-pick |
| ![스키마 탭](docs/demo/03-schema.gif) | ![유전체 표](docs/demo/04-explore.gif) |
| 스키마: 간단히 보기와 테이블 구조도 | 유전체 표와 상세 패널(값마다 출처) |

**분석** (Clade Ib 데이터셋의 실제 데이터, 유전체 id는 나오지 않음):

| | |
|---|---|
| ![전파 지도](docs/demo/06-spread-map.gif) | ![변이 그림](docs/demo/07-snipit.gif) |
| 전파: 세계 지도, 첫 관찰 라벨이 붙은 타임라인, 라인차트, 대표 변이 보기 | snipit 스타일 변이 그림: lineage 띠, APOBEC3 점, 유전자 확대 |
| ![유행 분석](docs/demo/08-outbreak.gif) | ![신규 변이와 PCA](docs/demo/09-novel-pca.gif) |
| 유행 분석, APOBEC3 대 배경 변이, 유전자 hotspot | 신규 변이 탐색과 변이 묶음(PCA) |


## 설치

Python 3.10 이상과 세 소스에 접속할 인터넷이 필요합니다. x86-64와 ARM(Apple Silicon, Linux aarch64)을 지원합니다.

**1. conda / mamba (권장)**

```bash
git clone https://github.com/CHOICHI-HWAN/GenomeCurator.git
cd GenomeCurator
mamba env create -f environment.yml      # conda라면 conda env create -f environment.yml (bioconda의 Nextclade 포함)
mamba activate genomecurator
python -m pytest -q                      # 선택: 전부 통과하면 정상
```

*한 번에 시작(환경이 없으면 만들고, 서버를 켜고, 브라우저를 엽니다):* `./scripts/start.sh`(Linux), `scripts/start.command`(macOS), `scripts\start.bat`(Windows).

**2. Docker** (linux/amd64, linux/arm64)

```bash
docker build -t genomecurator .
docker run -d -p 127.0.0.1:8001:8001 -v genomecurator-data:/data genomecurator
```

*로그인 기능이 없으므로 포트는 127.0.0.1로만 공개하세요.*

**3. pip**

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

**분석 도구(선택)** — 앱 환경이나 Docker 이미지에 들어가지 않는 별도 conda 환경입니다:

```bash
mamba env create -f nextstrain/environment.yml   # gc-nextstrain: nextstrain-cli, augur, snakemake, nextclade, iqtree, mafft, auspice
mamba env create -f squirrel/environment.yml     # gc-squirrel: squirrel, minimap2, mutation-profile
```

*Nextstrain 탭*이 두 환경을 점검하고, Nextstrain mpox 워크플로와 lineage 정의를 버튼으로 받습니다(`data/`에 `git clone`). GPL / AGPL 프로그램을 별도 프로그램으로 실행하는 것이며 [docs/LICENSES.md](docs/LICENSES.md)에 정리했습니다.

## 실행

```bash
python -m genomecurator watch-add "Henipavirus nipahense"   # 바이러스 등록 (이름은 NCBI에서 확인합니다)
python -m genomecurator serve --open                        # 웹 화면 http://127.0.0.1:8000
```

브라우저에서 **지금 확인 → 검토하기 → 항목 선택 → Approve → Apply** 순서로 진행합니다. 각 화면 위에 짧은 파란 안내 상자가 있습니다.

다른 명령(전체 옵션은 `python -m genomecurator --help`):

| 명령 | 하는 일 |
|---|---|
| `check <이름>` | 세 소스를 조회해 changeset 생성 (데이터는 그대로) |
| `dry-run <번호>` / `apply <번호>` | 적용 전 점검 / 승인된 항목 적용 |
| `reconcile <이름>` | 소스 간 일치 보고서 |
| `verify` | 로컬 데이터 무결성 검사 |
| `export <바이러스> [--profile nextstrain\|plain] [--out 폴더]` | 다른 프로그램용 `metadata.tsv`, `sequences.fasta`, `manifest.json` 만들기 |
| `dataset-create / dataset-list / dataset-export / dataset-delete` | 이름 붙여 저장한 선별과 추적 가능한 내보내기 |
| `exclude <바이러스> <id들> [--reason ..] [--undo]` | 삭제하지 않고 데이터셋에서 제외 |
| `private-check` / `private-import` | 비공개 유전체의 dry-run과 업로드(스프레드시트 또는 TSV + FASTA) |
| `analyze <바이러스> <데이터셋>` | 데이터셋으로 Nextstrain 워크플로 실행 (다른 분석은 *Nextstrain* 탭) |
| `pathoplexus <바이러스> [--check]` | 비공개 유전체를 점검하고 제출 파일 준비 (업로드는 하지 않음) |
| QC 탭 | 바이러스별 서열 통계·Nextclade, 데이터셋 업데이트, *QC 실행 + 확인* |
| `schema-check [--daily]` | 소스·로컬 스키마 점검 |
| `backup [zip]` / `restore <zip>` | 큐레이션 결과, 데이터셋 정의·제외 표시, FASTA의 이식 가능한 백업 (비공개는 `--with-private`일 때만) |
| `serve [--port N] [--open] [--superset] [--no-scheduler]` | 웹 서버 (Superset도 함께 켤 수 있음) |

주요 설정(환경변수): `GC_MIN_COVERAGE_PCT`(기본 80: 기준 게놈 대비 이보다 짧으면 큐레이션 안 함), `GC_DATA_DIR`(데이터 위치), `NCBI_API_KEY`(NCBI 요청 한도 증가), `GC_CHECK_INTERVAL_HOURS`(기본 24), `GC_NEAR_COMPLETE_PCT`(기본 90), `GC_SCHEMA_CHECK_HOUR`(기본 3). 전체 목록은 [가이드](docs/GUIDE.ko.md#설정-환경변수)에 있습니다.

## 출력 파일

```
data/
├── genomecurator.db          # SQLite: 큐레이션된 표, 출처 기록, 이력, 작업 기록
├── genomes/<바이러스>/*.fna  # 받아 둔 서열
├── private/                  # 비공개 유전체 (공개 산출물·일반 백업에는 나가지 않음)
├── nextstrain/               # 워크플로 체크아웃, 보관된 분석과 Auspice 파일
├── lineage-designation/      # lineage 정의 (버튼으로 받음)
├── pathoplexus/              # 제출 패키지 (파일만)
├── snapshots/  staging/      # 적용 중 안전 복사본과 임시 파일
└── logs/genomecurator.log    # 서버 로그 (회전, 작업·바이러스 맥락 포함)
```

`data/`는 git에 올라가지 않습니다. BI 도구용으로 `v_genome_analysis` 같은 분석 뷰도 DB에 들어 있습니다.

## 한계

- 자동 확인은 서버가 켜져 있어야 돕니다.
- 매일 스키마 점검은 소스당 몇 건만 표본으로 읽으므로, 가끔 스키마 탭에서 전체 점검을 눌러 주세요.
- 서열이 같다는 것은 "같은 레코드"의 증명이지 숙주·날짜가 옳다는 증명이 아닙니다.
- 분석: 합의 서열로는 환자 내 소수 변이를 알 수 없습니다. *첫 관찰*은 *유입*이 아닙니다(나라마다 시퀀싱 양이 다름). 신규 변이 후보와 서브클러스터는 지지도 없는 트리 하나에서 나온 가설이고 dN/dS는 탐색용입니다. 자세한 내용은 [docs/VALIDITY.md](docs/VALIDITY.md).

## 문서

자세한 한국어 가이드: [docs/GUIDE.ko.md](docs/GUIDE.ko.md) (탭별 설명, 정합성 점검 원리, FAQ, 로그 읽는 법, 코드 지도).
라이선스·데이터 이용 정책 검토: [docs/LICENSES.md](docs/LICENSES.md) · 전체 워크플로: [docs/WORKFLOW.md](docs/WORKFLOW.md) · 분석의 타당성: [docs/VALIDITY.md](docs/VALIDITY.md).
함께 보세요: [CHANGELOG.ko.md](CHANGELOG.ko.md) · [docs/SUPERSET.md](docs/SUPERSET.md) · [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md) · [docs/RELEASING.md](docs/RELEASING.md).

## 기여와 라이선스

버그와 아이디어는 [Issues](https://github.com/CHOICHI-HWAN/GenomeCurator/issues)로 보내 주세요. GNU General Public License v3.0 이상([LICENSE](LICENSE)), © 2026 CHOICHI-HWAN. v0.2.0까지의 릴리스는 MIT 라이선스였습니다.
데이터를 내주는 곳들에 예의를 지키려고 요청 사이에 간격을 둡니다(`GC_REQUEST_DELAY_MS`).
