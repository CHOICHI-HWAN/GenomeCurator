# Virus Genome Metadata Curator — 상세 가이드 (한국어)

> 처음이라면 [README](../README.ko.md)부터 보세요. 이 문서는 화면별 설명, 정합성·스키마 점검의 원리, FAQ, 로그 읽는 법, 코드 지도 등 **자세한 내용**을 담고 있습니다.

## 한 줄 요약 (Synopsis)

```bash
mamba env create -f environment.yml && mamba activate genomecurator

python -m genomecurator watch-add "Henipavirus nipahense"   # 이 바이러스를 지켜보세요
python -m genomecurator serve --open                         # 그리고 웹에서 승인하세요
```

그게 전부입니다. 나머지는 이 문서가 아니라 웹 화면이 알려 줍니다.

<div align="center">

![빈 화면에서 바이러스를 등록하고, 확인하고, 승인해서 반영하기까지](demo/01-quickstart.gif)

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

### Docker (x86 · ARM 공용 이미지)

```bash
docker build -t genomecurator .
docker run -d -p 127.0.0.1:8001:8001 -v genomecurator-data:/data genomecurator    # http://127.0.0.1:8001
```

이미지는 `environment.yml` 그대로 conda 환경을 만들며 `linux/amd64`, `linux/arm64`에서 돌아갑니다. 데이터는 `/data` 볼륨에 남습니다.
로그인 기능이 없으므로 위처럼 **127.0.0.1로만** 공개하세요.

### conda 없이

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

ARM(Apple Silicon, Raspberry Pi·Graviton 등)도 conda-forge에 모든 의존성이 있어 같은 방법으로 설치됩니다.

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

### 메뉴 구조

상단 메뉴는 6개 묶음이고, 묶음을 누르면 바로 아래에 그 묶음의 하위 탭이 나옵니다.

| 묶음 | 하위 탭 | 하는 일 |
|---|---|---|
| 큐레이션 | 변경 검토 · 정합성 · 스키마 · 이력 | 소스에서 모으고, 비교하고, 승인·적용 |
| 유전체 | 유전체 표 · 대시보드 · QC · 비공개 저장소 | 큐레이션된 유전체를 보고, 서열 QC, 비공개 유전체 관리 |
| 데이터셋 | 데이터셋 · 내보내기 | 분석에 쓸 선별(규칙 + 수동), 제외 표시, 다른 프로그램용 내보내기 |
| 분석 | **분석 허브** · 실행·이력 | 데이터셋별 분석 실행과 결과 화면 |
| 제출·공유 | Pathoplexus 제출 · Superset | 제출 파일 준비(업로드는 직접), BI 도구용 보기 |
| 시스템 | 작업 | 서버에서 도는 작업의 상태·로그 |

**분석 허브**(`/hub`)가 분석의 시작점입니다: 데이터셋을 고르면 Nextstrain·Squirrel·변이 프로파일 실행이 오늘의 데이터와 같은지 한눈에 보이고, 볼 실행을 골라 결과 화면 7개(*변이*: 변이 그림·히트맵·묶음, *계통·전파*: 타임라인·지도, *통계·가설*: 변이 분석·유행 분석·신규 변이)를 엽니다. 결과 화면 위의 위치 표시와 탭 줄로 같은 실행의 다른 화면으로 바로 옮겨 갈 수 있습니다. 공개 보기와 배포하기는 오른쪽 위에 그대로 있습니다.

## 화면 구경 (The tour)

| 탭 | 하는 일 |
|---|---|
| **큐레이션** | 바이러스 목록(추가·편집·삭제) + *검토할 변경 사항* 카드. "새로 발견 / 정보 바뀜 / 사라짐" 숫자, 승인·거절·보류 진행 막대, 다음 할 일을 한 문장으로 알려 줍니다. |
| **대시보드** | 월별 등록, 수집 연도, 국가, 숙주 같은 차트. |
| **유전체** | 엑셀처럼 필터·정렬·열 선택·TSV 내보내기. 오른쪽 상세 패널은 **끌어서 너비를 조절**할 수 있고, Pathoplexus 변이 정보도 불러옵니다. |
| **정합성** | 세 소스가 서로 얼마나 일치하는지 한눈에. 일치(초록) / 표기만 다름(노랑) / 진짜 다름(빨강). 서열로 검증하는 버튼도 여기. |
| **스키마** | 소스와 내 DB의 *구조*가 코드의 기대와 맞는지 점검. **간단히 보기**(내 열 × 소스 3곳 표, 확인이 필요한 열이 맨 위에 모임)와 **자세히 보기**(소스 필드가 내 `genome` 열에 이어지는 선까지 보이는 구조도)를 탭으로 오가며, 이상이 있는 필드·열을 색으로 강조합니다. 문제가 있으면 탭에 `!` 배지가 뜹니다. |
| **QC** | 서열 통계·Nextclade(clade/lineage, 품질) 결과, 데이터셋 버전과 업데이트 버튼, *QC 실행 + 확인*. |
| **분석** | Apache Superset 켜기/끄기/열기(설치 여부도 안내). 별도 프로그램이라 처음 한 번만 `superset.sh init`이 필요합니다. |
| **작업** | 서버에서 도는 모든 작업(추가·조회·확인·적용·다운로드·검증·점검·배포·자동 실행)의 상태·진행률·로그·오류. 다른 메뉴로 가도 계속되고 여기서 다시 봅니다. |
| **이력** | 누가 언제 무엇을 했는지 + 서버 로그(경고만 / 오류만 필터). |
| **공개 보기 / 배포하기** | 읽기 전용 정적 페이지를 GitHub Pages로 게시. |

오래 걸리는 작업은 화면 가운데에 진행 모달(큰 진행 막대 + 로그)이 뜨고, 새로고침해도 이어집니다.

<div align="center">

![유전체 표: 키워드 검색, 필터, 상세 패널, 패널 너비 조절](demo/04-explore.gif)

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

![정합성 탭: 소스별 일치 막대, 실제로 다른 값 목록, 서열 검증 dry-run 계획](demo/02-integrity.gif)

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

![스키마 탭: 간단히 보기 표와 구조도, 그리고 소스가 바뀌었다고 가정했을 때의 강조](demo/03-schema.gif)

<sub>▶ 스키마 점검의 간단히 보기와 테이블 구조도 — 처음은 실제 결과(이상 없음), 뒤쪽은 <b>시연을 위해 "소스가 필드를 지우고 이름과 응답 구조를 바꿨다"고 가정</b>한 화면입니다 (실제로 바뀐 것이 아닙니다)</sub>

</div>

### 4. 다른 값 중에서 고르기 (resolve)

소스끼리 값이 다르면 **검토 화면에서 항목을 눌러** 무엇이 다른지 보고 직접 고릅니다.

* 오른쪽 패널에 다른 필드마다 **소스별 값**이 나란히 나오고, 표기 차이인지 진짜 불일치인지 표시됩니다.
* **★ 추천**이 미리 골라져 있고 **이유와 신뢰도**가 붙습니다. 규칙은 이렇습니다: 불가능한 날짜는 뺀다 → 길이·버전은 NCBI(서열을 NCBI에서 받아 검증하므로) → 같은 사실의 다른 표기는 가장 잘 쓴 쪽(가장 구체적인 날짜, 학명, 표준 국가명) = **높음** → 소스 다수가 같은 값 = **보통** → 모두 다르면 우선순위로 제안하되 **낮음**(직접 확인).
* 값을 고르거나 직접 입력해서 저장하면 그 항목에 **바로 반영**되고, **다음 `check`에서도 유지**됩니다(소스가 바뀌어도 내가 정한 값을 존중). 언제든 기본값으로 되돌릴 수 있습니다.
* 표에서 행을 **끌어서**(드래그) 또는 **Shift+클릭**·체크박스로 여러 건을 고르고 **선택 항목에 추천 적용**을 누르면, 고른 항목만 추천값으로 한 번에 선택합니다("보통 이상" 또는 "낮음까지 모두"를 고르는 창이 뜹니다).
* 툴바의 **Auto-pick (전체)**는 신뢰도 "보통" 이상인 추천을 한 번에 선택하고, **"낮음"은 사람이 보도록 남겨 둡니다.** 어차피 소스가 둘뿐인 바이러스는 대부분 "낮음"이니, 근거 없는 일괄 적용을 막는 장치입니다.

<div align="center">

![검토 화면에서 항목을 눌러 소스별 값을 보고, 추천을 참고해 고르고, Auto-pick으로 일괄 선택하기](demo/05-resolve.gif)

<sub>▶ 검토 화면: 항목을 누르면 소스별 값과 추천(이유·신뢰도) → 고르기 → "선택 완료" → Auto-pick</sub>

</div>

### 5. 내 데이터가 멀쩡한지 (verify)

FASTA 실제 길이·md5가 기록과 맞는지, 날짜가 미래가 아닌지, 고아 행이 없는지 등을 봅니다. 문제가 있으면 종료 코드 1.

### 6. 서열 QC (서열 통계 · Nextclade)

저장된 서열을 한 번 더 검사합니다 (**QC** 탭). 경고만 하며 아무것도 막지 않고, 큐레이션 값을 몰래 바꾸지 않습니다. 결과는 `qc_result` 테이블에 쌓이는 *파생 데이터*라서 지우고 다시 계산해도 됩니다.

* **서열 통계**(모든 바이러스): 길이, N 비율(1% 이상 주의 · 5% 이상 불량), 모호 염기(0.1% 이상 주의), GC, **기록된 길이와 서열 길이가 같은지**(다르면 불량).
* **Nextclade**(데이터셋이 있는 바이러스만): clade, lineage와 자체 품질 점수(결측·혼합 부위·프레임시프트·정지코돈 등 → 양호/주의/불량). 데이터셋이 없는 바이러스(예: Nipah)는 서열 통계만 합니다.
* **데이터셋**: 필요할 때 Nextstrain에서 내려받아 `data/nextclade/`에 둡니다(저장소에는 넣지 않습니다). 바이러스마다 로컬 버전과 새 버전 유무가 보이고, 오래됐으면 **데이터셋 업데이트** 버튼이 나옵니다. 어떤 바이러스에 어떤 데이터셋을 쓰는지는 `qc.py`의 `DATASETS`(종명 → 데이터셋)이고, `GC_NEXTCLADE_DATASETS="종명=데이터셋경로;…"`로 추가할 수 있습니다.
* **큐레이션과의 연결**: Nextclade가 부른 clade/lineage는 NCBI·BV-BRC·Pathoplexus와 나란히 **네 번째 소스**로 참여합니다. 소스와 같으면 아무 일도 없고, 소스가 비어 있으면 채우는 후보, **다르면** 검토 화면의 "소스 간 불일치"에 Nextclade 값이 함께 올라옵니다(다른 데이터베이스가 같은 값이면 추천 신뢰도 *높음*, 아니면 *보통*). 고르는 건 사람이고, 선택은 이후 확인에서도 유지됩니다. **QC 실행 + 확인** 버튼이 QC와 지금 확인을 한 작업으로 합니다.
* 같은 서열(md5)과 같은 도구·데이터셋 버전이면 다시 계산하지 않습니다. 유전체 화면에는 QC · Nextclade Clade/Lineage 열과 QC 필터가 있습니다.
* CheckM2는 세균·고세균용이라 쓰지 않습니다. 바이러스 완성도는 Near complete 등급(참조 게놈 대비 길이)으로 봅니다.

### 7. 다른 프로그램에 넘기기 (export)

큐레이션한 유전체를 **내보내기** 탭(또는 `python -m genomecurator export`, 내부 API)으로 넘깁니다. 같은 id·같은 순서의 `metadata.tsv`와 `sequences.fasta`, 그리고 `manifest.json`(건수·경고·SHA-256·Nextclade 데이터셋 버전)과 `nextstrain-config.yaml`이 한 폴더(기본 `data/exports/<바이러스>/<시각>/`)에 만들어집니다.

* **형식 `nextstrain`**: [nextstrain/mpox](https://github.com/nextstrain/mpox)의 `phylogenetic/` 워크플로가 읽는 열에 맞춥니다. id는 **`PPX_accession`**(FASTA 헤더와 같음; 워크플로가 포함·제외 목록을 PPX accession으로 바꿔 `augur filter --exclude`와 맞추므로 INSDC id를 쓰면 제외가 조용히 무시됩니다. Pathoplexus accession이 없는 유전체는 빠지고 미리보기가 건수를 알려 줍니다), 날짜는 `YYYY-MM-DD`에 모르는 부분은 `XX`(`2023-XX-XX`), `region`(대륙)·`country`·`division`·`location`은 국가·`geographic_location`에서 만들고, `clade/lineage/outbreak`는 큐레이션 값이 우선이며 없으면 Nextclade 값, `coverage`·`QC_rare_mutations` 등 Nextclade QC 열은 QC 탭의 결과입니다. **QC를 먼저 실행하세요**: 워크플로의 기본 필터가 `QC_rare_mutations`를 보기 때문에 값이 없는 행은 걸러집니다(미리보기에서 경고로 알려 줍니다).
* **형식 `plain`**: `genome` 표의 모든 열과 QC 열을 이름 그대로 줍니다.
* **필터**: 최소 길이, 완성도(Complete/Near complete/Partial), QC 불량·오류 제외, clade, 서열이 있는 것만(기본).
* **공개(Open) 데이터만** 나갑니다. 수집 대상이 NCBI(INSDC) accession으로 정해지므로 Pathoplexus 제한 사용 데이터는 들어 있지 않고, `dataUseTerms`는 `OPEN`으로 적힙니다.
* **내부 API**: `GET /api/export/{바이러스id}/metadata.tsv`, `.../sequences.fasta`, `.../preview`, `POST .../build`(서버 작업), `GET /api/exports`. 다른 프로그램에서는 `curl`로 바로 받을 수 있습니다.
* **쓰는 법(mpox)**: QC 탭에서 QC 실행 → `git clone https://github.com/nextstrain/mpox` → 체크아웃 안으로 바로 내보내기:
  ```bash
  python -m genomecurator export "Orthopoxvirus monkeypox" --nextstrain-dir ~/mpox/phylogenetic
  cd ~/mpox/phylogenetic
  nextstrain build . --configfile defaults/mpxv/config.yaml build-configs/genomecurator/config.yaml   # 전체 clade
  ```
  * **config를 두 개** 줘야 합니다. 워크플로의 `Snakefile`은 `defaults/*/config.yaml`을 *config가 비어 있을 때만* 읽으므로, 우리 overlay만 주면 `build_name`이 없다는 `KeyError`가 납니다(확인함). overlay는 `inputs`(로컬 파일)와 `use_remote_files: false`만 덮어씁니다. 다른 빌드는 `defaults/hmpxv1/config.yaml`(clade IIb), `defaults/clade-i/config.yaml`.
  * 기본 입력이 매번 웹(`data.nextstrain.org`)에서 받는 URL이라 로컬 데이터를 쓰려면 이 `inputs` 교체가 필요하고, 워크플로 코드(`path_or_url`)는 로컬 경로를 그대로 쓰므로 **규칙을 고칠 필요는 없습니다**.
  * 파일은 체크아웃의 `data/genomecurator/`에, 설정은 `build-configs/genomecurator/config.yaml`에 상대 경로로 써집니다. `--nextstrain-dir` 없이 `--out`을 쓰면 절대 경로라 Docker 런타임(`nextstrain build --docker`)에서는 파일이 보이지 않으니 `--ambient`/`--conda`로 돌리세요.
* 확인한 것: 위 내보내기로 만든 파일에 Nextstrain mpox 저장소의 `scripts/map_accessions.py`(실제 `exclude.txt`), `augur filter`(기본 쿼리), Snakemake의 `decompress`·`map_accessions`·`filter` 규칙을 실제로 돌렸고, 전체 DAG(19단계)가 원격 입력 없이 풀리는 것까지 봤습니다. 정렬·트리 구축(`align`~`export`)은 서열이 충분해야 의미가 있어 돌려 보지 않았습니다. 워크플로 설정이 바뀌면 열이 달라질 수 있습니다.

### 8. 비공개 유전체와 최종 분석용 데이터세트 (private)

목표: 공개 데이터(위의 큐레이션 세트)와 **로컬 비공개 유전체**를 합친 분석용 데이터세트를 만들어 Nextstrain에서 돌리는 것입니다. 데이터세트를 준비하고, 필요하면 *Nextstrain 탭*(아래 11번)에서 워크플로까지 돌립니다.

* **분리 보관**: 비공개 유전체는 `genome` 표가 아니라 별도 저장소(`private_genome` 표 + `data/private/<바이러스>/`)에 둡니다. 큐레이션 세트는 NCBI accession으로 정해지므로 `genome`에 넣으면 다음 *지금 확인*에서 "NCBI에서 사라짐"이 되고, 공개 페이지·Superset·TSV 다운로드·백업에도 섞이기 때문입니다. 분리해 두면 그 어느 것에도 나타나지 않습니다.
* **올리는 두 가지 방법** (**비공개** 탭):
  1. **스프레드시트로 입력**: 페이지 안에 엑셀처럼 이어진 **진짜 스프레드시트**([Jspreadsheet CE](https://github.com/jspreadsheet/ce), MIT)가 있습니다. 한 줄이 유전체 하나이고 한 건만 올리려면 한 줄만 씁니다. 칸에 바로 쓰고, 끌어서 범위를 선택하고, **엑셀·구글 시트에서 복사해 Ctrl+V**, 모서리를 끌어 채우고, 오른쪽 클릭으로 줄·열을 삽입/삭제합니다(`+ 열 추가`로 Nextstrain에 그대로 전달할 열을 더합니다). **파일을 시트로 불러오기**로 기존 `.xlsx`/`.tsv`/`.csv`를 열어 고칠 수도 있습니다. **줄마다 서열 칸**이 있습니다: 시트 마지막 열 *서열*을 **더블클릭**하면 그 줄의 서열을 붙여넣거나 FASTA 파일을 고를 수 있고(FASTA 파일을 그 칸에 **끌어다 놓아도** 됩니다), 칸에는 `✔ 길이 bp`가 표시됩니다. 여러 건은 시트 아래 **여러 서열 한꺼번에 (FASTA)**에서 파일(여러 개 가능)이나 붙여넣기로 연결하며 헤더가 시트의 id와 같은 줄에 붙고, 시트에 없는 id는 새 줄이 됩니다. 서열은 id에 묶여 있으니 id 칸을 고치면 연결이 끊긴 것으로 표시됩니다. *검사*가 문제 있는 줄을 색으로 표시하고(빨강 = 오류, 노랑 = 경고), 입력 중인 시트 내용은 이 브라우저에 임시 저장됩니다(서열은 저장하지 않으니 새로고침하면 다시 연결). 시트 라이브러리는 CDN에서 불러오므로 이 탭에는 인터넷이 필요합니다.
  2. **스프레드시트 파일(.xlsx)로 올리기**: 탭에서 **`.xlsx` 템플릿**을 받아 `genomes` 시트에 한 줄에 하나씩 채웁니다(필수 열은 노란색, `열 설명` 시트에 열별 설명, mpox는 clade 드롭다운). `.tsv`/`.csv`도 됩니다. FASTA는 따로 올리며 헤더는 id와 같아야 합니다. 시트에 다른 열을 더 넣어도 보존됩니다.
  * **이 업로드의 기본 clade**: clade 칸이 빈 행에 채울 값(mpox는 *Ib*가 미리 선택되어 있음). 직접 쓴 clade는 덮어쓰지 않고, Nextclade가 다르게 판정하면 ⚠로 표시합니다.
* **올릴 때마다 버전**: 같은 id를 바뀐 내용으로 올리면 v2, v3…으로 기록됩니다(바뀐 열, 서열이 바뀌었는지, 원본 파일, 메모). 같은 내용이면 변경 없음이고 새 버전이 생기지 않습니다. 이전 버전(서열 포함)은 삭제한 뒤에도 남아 있고, *이 버전으로 복원*은 그 내용을 **새 버전**으로 만듭니다(이력은 고쳐 쓰지 않습니다). 서열은 `data/private/<바이러스>/`(현재)와 `_versions/`(버전별)에 있습니다. 내보내기 manifest에는 내보낸 버전이 기록됩니다.
* **커버리지 하한**: 비공개 업로드도 `GC_MIN_COVERAGE_PCT`(기본 80)를 따라, 기준 게놈 대비 짧은 서열은 검사에서 오류로 막힙니다.
* **가져오기**(**비공개** 탭 또는 CLI): Nextstrain에 준하는 메타데이터 TSV(id 열 `PPX_accession` 또는 `accession`, 그리고 `strain`, `date`, `country`, `division`, `location`, `host`, `clade` …)와 헤더가 그 id인 FASTA. 먼저 **검사(dry-run)**: 구조 오류(id 누락·중복·공개 accession과 충돌·서열 없는 행·행 없는 서열·서열 문자 오류)는 가져오기를 막고, 날짜(`2024-03` → `2024-03-XX`)·대륙(국가에서)·길이(서열에서)는 맞추고 경고로 알립니다. 같은 id를 다시 가져오면 갱신합니다.
  ```bash
  python -m genomecurator private-check  "Orthopoxvirus monkeypox" --metadata my.tsv --fasta my.fa
  python -m genomecurator private-import "Orthopoxvirus monkeypox" --metadata my.tsv --fasta my.fa
  python -m genomecurator qc "Orthopoxvirus monkeypox" --private        # 서열 통계 + Nextclade
  ```
* **QC**: 공개 유전체와 같은 서열 통계·Nextclade. Nextclade가 부른 clade/lineage/outbreak와 QC 열이 비어 있는 칸을 채우고, 직접 입력한 clade가 Nextclade와 다르면 비공개 탭에 ⚠로 표시합니다(덮어쓰지 않음).
* **내보내기**: `export --with-private`(탭의 *비공개 포함*). 공개 파일(`metadata.tsv`, `sequences.fasta`)에는 **비공개가 섞이지 않고**, 비공개는 **별도 파일**(`metadata-private.tsv`, `sequences-private.fasta`)로 나가 설정 overlay의 `private_metadata` / `private_sequences`에 연결됩니다. 워크플로가 이를 공개 데이터 필터 *뒤에* 붙이므로 비공개 행은 날짜·QC 필터를 받지 않습니다(`source_private` 열로 구분되어 Auspice에서 색을 입힐 수 있습니다). manifest에 `contains_private`, 폴더에 `.gitignore`(`*private*`)와 `PRIVATE_DATA.txt`가 생깁니다. 이 폴더는 공유하지 마세요.
* **`accession` 열**: 워크플로의 병합 스크립트(`combine_data_sources.py`)가 아직 `accession` 열로 행을 맞춰서, 현재의 PPX 기반 메타데이터에서는 이 열이 없으면 `metadata file missing 'accession'`으로 실패합니다(확인함). 내보내기는 공개·비공개 파일 모두에 id와 같은 `accession` 열을 넣어 이를 피합니다. 수정 없는 스크립트로 공개 7 + 비공개 3행이 병합되는 것까지 확인했습니다.
* **백업**: 일반 백업(`backup`)에는 비공개가 **들어가지 않습니다**(백업 파일을 공유할 수 있도록). 비공개까지 지키려면 `backup --with-private`를 쓰고, 그 파일은 비공개 데이터와 똑같이 다루세요.

### 9. 데이터셋에서 제외 (삭제하지 않고 보관)

최종 분석용 데이터셋에서 빼고 싶은 유전체(조립이 이상한 것, 오염 의심, 분석 대상이 아닌 것…)는 **삭제하지 말고 제외 표시**를 합니다.

* **유전체 탭**에서 유전체를 고르고(체크박스, 행 끌기, Shift+클릭, 필터 후 Select all) 툴바의 **데이터셋에서 제외**를 누릅니다. 사유는 선택이며 `excluded.tsv`에 기록됩니다. **제외 해제**로 되돌립니다. *데이터셋* 열(포함/제외)과 필터로 상태를 봅니다. 비공개 유전체는 **비공개 탭**의 줄마다 제외/풀기 버튼이 있습니다.
* 표시일 뿐이라 유전체·서열·이력·QC는 **그대로 보관**됩니다. *지금 확인*·적용·QC도 평소처럼 동작하고, 제외 표시는 새 changeset이 와도 유지됩니다.
* **내보내기**는 제외된 유전체를 건너뜁니다. 미리보기에 건수가 경고로 나오고, 내보낸 폴더에 `excluded.tsv`(accession, strain, 사유, 시각; 비공개 포함 내보내기에는 `excluded-private.tsv`도)와 manifest의 `excluded` 항목이 남아 어떤 데이터셋에서 무엇이 왜 빠졌는지 추적할 수 있습니다.
* `Removed`와 다릅니다: `Removed`는 NCBI에서 사라졌거나 커버리지 하한 미만이라 *큐레이션 세트에서 빠진* 것이고, 제외는 *특정 데이터셋에서만 빼는* 결정입니다.
* CLI: `python -m genomecurator exclude "Orthopoxvirus monkeypox" ACC1 ACC2 --reason "오염 의심"` (`--undo`로 해제, `--private`로 비공개 id, `--list`로 목록). API: `GET/POST /api/exclusions/{바이러스id}`.

### 10. 데이터셋 (분석할 유전체를 골라 저장하기)

분석마다 대상이 다릅니다(예: clade Ib만, 2024년 이후만, 완전한 게놈만). **데이터셋** 탭에서 대상을 *정의*해 이름으로 저장하고 필요할 때마다 내보냅니다.

* **정의일 뿐입니다**: 유전체·서열을 복사하거나 지우지 않습니다. 구성: **규칙**(모두 만족해야 하며 비우면 제한 없음) — clade, lineage, outbreak, 대륙, 국가, 숙주, 완성도, 수집 연도 범위, 길이 범위, 뺄 QC 상태, 서열이 있는 것만; **직접 넣기**(규칙과 상관없이 포함, accession 한 줄에 하나)와 **직접 빼기**(이 데이터셋에서만 제외); **비공개 포함**(별도 파일로).
* **규칙 고르기**: 탭이 데이터에 실제로 있는 값을 건수·서열 있는 수와 함께 보여 줍니다(많으면 찾기). 클릭할 때마다 **멤버 수**, clade·국가·연도 분포, 경고(서열이 없어 빠지는 수 등)가 바로 나옵니다. 처음이면 **프리셋**(mpox: clade Ib / Ia / I(Ia+Ib) / IIb / 완전한 게놈)으로 시작하세요.
* **판단 기준**: 규칙은 내보내기가 쓰는 *그 행*(큐레이션된 clade, 없으면 Nextclade 값 등)으로 판정하므로 "미리보기가 센 것"과 "내보낸 것"이 같습니다. 전역 **데이터셋에서 제외** 표시는 모든 데이터셋에서 빠지고, *직접 빼기*가 *직접 넣기*와 규칙보다 우선합니다.
* **내보내기**: `exports/<바이러스>/<데이터셋 이름>/<시각>/`에 `metadata.tsv`, `sequences.fasta`, `manifest.json`, 설정(비공개 포함이면 별도 파일도)을 씁니다. manifest의 `dataset`에 이름·규칙·**멤버 id 해시**가 기록되고 데이터셋의 *내보낸 기록*에도 남으므로, 규칙이나 데이터가 나중에 바뀌어도 어떤 멤버로 만든 결과인지 추적할 수 있습니다(같은 멤버면 해시가 같습니다). **Nextstrain 빌드**(`mpxv`/`hmpxv1`/`clade-i`)를 정해 두면 설정 파일이 그 빌드의 명령을 안내합니다(Ib는 `clade-i`).
* **정의 삭제**는 정의와 내보낸 기록만 지웁니다. 유전체·서열·이미 내보낸 파일은 그대로입니다.
* CLI: `python -m genomecurator dataset-create "Orthopoxvirus monkeypox" "clade Ib" --clade Ib --build clade-i` / `dataset-list` / `dataset-export ... [--nextstrain-dir ~/mpox/phylogenetic]` / `dataset-delete`. API: `GET/POST /api/datasets`, `POST /api/datasets/{id}/export`, `GET /api/export/{virus}/metadata.tsv?dataset={id}`(+ `sequences.fasta`; 다운로드에는 비공개가 실리지 않음).
* **데이터셋 QC**: 데이터셋 보기 위쪽의 **이 데이터셋 QC 실행**(또는 데이터셋 편집 화면의 *QC 실행*)이 바이러스 전체가 아니라 *그 데이터셋의 멤버만* 서열 통계·Nextclade로 검사합니다(서버 작업, 비공개 포함이면 비공개 멤버도). 바뀌지 않은 서열은 건너뛰고 **다시 계산**은 모두 다시 합니다. 같은 카드에 도구별 양호/주의/불량과 *아직 QC를 안 한 멤버 수*가 나오고, 표의 **QC** 열·필터로 주의/불량만 골라 볼 수 있으며 데이터셋 규칙의 *뺄 QC 상태*로 내보내기에서 뺄 수도 있습니다. 내보내기 전에 QC를 해 두세요(Nextstrain 기본 필터가 `QC_rare_mutations`를 봅니다).
* **유전체처럼 보기**: 데이터셋의 **목록 보기**(또는 `/dataset/<id>`)는 유전체 탭과 같은 표(필터, 정렬, 열 선택, 상세 패널, TSV)를 그 데이터셋의 멤버로만 보여 줍니다. 서열이 없는 멤버도 표시되고(내보내기에서는 빠짐) 위의 배너 **모두 받기**는 이 데이터셋의 서열만 받습니다. 데이터셋이 *비공개 포함*이면 비공개 멤버가 🔒 *출처: 비공개*로 같은 표에 섞여 보이고 *출처* 필터로 나눌 수 있습니다(공개 페이지·일반 유전체 탭에는 여전히 나타나지 않습니다).
* **대시보드에서 보기**: 대시보드 위의 *데이터셋* 선택에서 고르면 모든 차트와 지표가 그 데이터셋의 멤버만으로 그려지고, 차트를 눌러 들어가는 목록도 그 데이터셋 안에서 필터됩니다.
* 서열을 아직 다 받지 않았다면 멤버 수가 서열 있는 것만 셉니다(미리보기가 건수를 알려 줍니다) — 유전체 탭에서 서열을 받은 뒤 다시 보세요.

### 11. Nextstrain 분석 (데이터셋 → 워크플로 → 보관된 결과)

**Nextstrain** 탭에서 데이터셋으로 nextstrain/mpox 워크플로를 돌립니다. 도구가 *내보내기 → 설정(config) 만들기 → 실행 → 결과 보관*을 이어서 합니다.

* **준비물(한 번만)**: ① 도구 환경 — 터미널에서 `mamba env create -f nextstrain/environment.yml`(이름 `gc-nextstrain`; nextstrain-cli, augur, snakemake, nextclade, iqtree, mafft, auspice가 들어 있고 약 1.7 GB, **앱 환경과 별개**입니다). ② 워크플로 저장소 — 탭의 **워크플로 받기**(`git clone https://github.com/nextstrain/mpox` → `data/nextstrain/mpox`; **업데이트** 버튼은 `git pull`이며 분석마다 어느 커밋으로 돌렸는지 기록됩니다). 탭 위쪽이 둘의 상태(✔/✖)를 보여 줍니다.
* **config가 어떻게 만들어지나**: 워크플로의 `Snakefile`은 `defaults/*/config.yaml`을 *설정이 없을 때만* 읽기 때문에, 도구는 **빌드 기본값 → 도구가 만든 overlay** 순서로 두 설정 파일을 넘깁니다(`nextstrain build . --configfile defaults/<빌드>/config.yaml build-configs/genomecurator/config.yaml`). overlay는 ① 웹 주소 입력을 로컬 파일로 바꾸고(`inputs`, `use_remote_files: false`), ② 비공개 포함이면 `private_metadata`/`private_sequences`를 연결하고, ③ 이 분석만의 `build_name`/`auspice_name`(= 실행 이름)을 정하고, ④ **이 데이터셋의 설정**을 덧붙입니다.
* **데이터셋별 설정(YAML)**: 탭의 입력칸에 워크플로 config 중 바꿀 키만 씁니다. 중첩 키는 합쳐집니다 — 예:
  ```yaml
  filter:
    min_length: 150000     # 빌드 기본 170000(clade-i)을 낮춤
  timetree: true           # 시간 기반 트리
  ```
  오른쪽에 **실제로 적용될 설정**(빌드 기본값 + 입력 + 내 설정)이 바로 보이고, `inputs`·`private_*`·`build_name`·`auspice_name`·`strain_id_field`는 도구가 관리하므로 쓰면 거부됩니다. 설정은 데이터셋에 저장됩니다. **빌드**(`mpxv`/`hmpxv1`/`clade-i`)는 데이터셋 정의에서 고릅니다(Ib는 `clade-i`).
* **실행**: *분석 실행*(또는 `python -m genomecurator analyze "Orthopoxvirus monkeypox" "clade Ib"`)은 서버 작업이며 단계별 진행률이 나옵니다. QC가 없어 워크플로의 기본 필터(`QC_rare_mutations`)에서 빠질 멤버가 있으면 **실행 전에 경고**합니다 — 먼저 데이터셋 QC를 하세요. 실행 방식은 `ambient`(Nextstrain 환경의 도구, 기본) 또는 `docker`입니다.
* **분석은 쌓이고 골라 봅니다**: 분석마다 *그때의* 멤버 해시·건수, 규칙, 설정, 워크플로 커밋, 도구 버전, 메모, 로그가 기록되고 결과는 `data/nextstrain/auspice/<실행이름>.json`, 입력 표·manifest·최종 설정은 `data/nextstrain/runs/<바이러스>/<실행이름>/`에 **계속 보관**됩니다. 다시 돌려도 이전 결과를 덮어쓰지 않으며, **멤버·설정·워크플로가 완전히 같으면** 반복하지 않고 기존 분석을 알려 줍니다(*같은 입력이어도 다시 분석*으로 강제). 데이터가 바뀐 뒤(유전체 추가, 제외, 비공개 업로드, 규칙 변경)에 돌리면 **새 분석**이 하나 더 생깁니다.
* **오늘과 비교**: 분석 표는 지난 분석마다 현재 데이터셋과 어떻게 다른지 보여 줍니다 — *멤버 같음* 또는 *데이터가 바뀜(현재 기준 +추가 / −삭제)*, *설정이 다름*, *워크플로 버전이 다름*, 그리고 완전히 같으면 *← 현재와 같은 입력*. 어느 결과가 지금 데이터에 맞는지 한눈에 알 수 있습니다.
* **결과 보기**: 분석 줄의 **Auspice로 보기** — 이 컴퓨터에서 `auspice view`(127.0.0.1만)를 켜서(처음 한 번 확인 창) 보관된 결과 중 그 분석을 엽니다. 데이터가 밖으로 나가지 않습니다. **로그**로 실행 로그 끝부분을, 메모 칸으로 분석에 설명을 남기고, **삭제**는 그 분석의 결과와 기록만 지웁니다(확인 후, 되돌릴 수 없음).
* 라이선스 참고: augur·auspice는 AGPL-3.0, IQ-TREE는 GPL-2.0+라 **별도 환경의 별도 프로그램**으로만 실행하고 저장소·Docker 이미지에 넣지 않습니다. Auspice 서버를 네트워크에 공개하면 AGPL 의무가 생길 수 있으니 `docs/LICENSES.md`를 보세요.
* 한계: 아직 실행 중인 분석을 중간에 취소하는 버튼은 없습니다(서버를 끄면 *중단*으로 기록됨). 큰 데이터셋의 트리 구축은 서열 수에 따라 오래 걸릴 수 있습니다(예: Ib 213건 + 비공개 3건 → 약 1분).

**검증(예제 데이터)**: 업스트림이 같이 주는 `example_data`(63건)를 ① 업스트림 CI 설정으로 직접, ② 우리 방식(입력만 바꾸는 오버레이 + 같은 설정)으로 돌려 비교했습니다. 팁 53개 동일, 트리 토폴로지 완전 일치(Robinson-Foulds 거리 0), 팁 간 SNP 거리 1,378쌍 전부 동일, clade·lineage 전부 일치. 차이는 모든 팁이 똑같이 갖는 뿌리 가지 변이 4개(팁 사이 거리에 영향 없음)와 날짜가 불확실한 2개 팁의 추정 날짜(TreeTime의 확률적 추정)뿐이었습니다. 공식 공개 데이터셋과의 비교는 하지 않았습니다.

### 12. Squirrel · 변이 프로파일 (같은 데이터셋, 다른 분석)

**Nextstrain 탭**에서 데이터셋을 고르면 Nextstrain 아래에 두 패널이 더 있습니다(공개+비공개 멤버가 함께 들어갑니다). 도구는 별도 conda 환경 `gc-squirrel`에 있습니다: `mamba env create -f squirrel/environment.yml`.

- **Squirrel**: MPXV 정렬·서열 QC·IQ-TREE 계통수·HTML 보고서. clade와 옵션은 데이터셋마다 저장됩니다. 계통수를 만들기 전에 QC만 먼저 돌려 Squirrel이 제안한 서열(N 비율 20% 초과 등)을 빼고(끌 수 있음), IQ-TREE에는 `-st DNA`를 알려 줍니다. 보고서는 *보고서 열기*로 보며 브라우저가 표·스타일 라이브러리를 CDN에서 불러옵니다(데이터는 어디로도 나가지 않음).
- **변이 프로파일**([mutation_profile](https://github.com/insapathogenomics/mutation_profile)): 각 유전체를 참조 서열에 (기본 RefSeq **NC_063383.1**, Nextclade·Squirrel과 같은 기준. **clade I 빌드는 Nextstrain 빌드의 기본 참조 KJ642613.1**(Congo_8)을 자동으로 씁니다 — 그래야 Nextstrain 트리와 좌표가 맞습니다. RefSeq clade I 참조 NC_003310.1은 *다른 유전체*(약 99.65% 일치, SNV 150·indel 약 550염기 차이)이며 설정에서 고를 수 있습니다(처음 한 번 NCBI에서 GenBank 기록을 받아 `data/reference/`에 보관); 설정에서 고를 수 있고 쓴 참조는 결과에 기록됩니다) minimap2로 정렬해 SNP를 읽고, 앞뒤 염기 맥락과 APOBEC3형(GA>AA, TC>TT) 프로파일을 붙입니다. **각 변이가 어느 유전자에 있는지**는 그 참조의 유전자 주석(워크플로의 `genome_annotation.gff3`)으로 찾아 유전자·산물·코돈 위치·아미노산 변화·동의/비동의/종결/유전자 사이를 함께 기록합니다. 결과: `mutations.tsv`(변이 한 줄씩), `gene_summary.tsv`(유전자별 요약), `summary.json`, 원본 입력 `snps.tsv`. 이웃 서열이 N이거나 삽입·결손인 곳은 SNP로 세지 않습니다.
- **변이 그림**: 변이 프로파일 결과 줄의 *변이 그림* 버튼이 [snipit](https://github.com/aineniamh/snipit)과 비슷한 그림을 엽니다(snipit 자체는 쓰지 않고 `mutations.tsv`로 직접 그립니다). 위쪽 참조는 **그 분석에 설정된 참조 서열**이고(clade I 빌드는 Nextstrain 기본 KJ642613.1), 서열마다 한 줄에 SNP를 색 눈금으로 그리며 왼쪽에 lineage·clade, 눈금 위 분홍 점으로 **APOBEC3형 변이**, 아래 **유전자 트랙**이 있습니다. 서열이 많으니 lineage·clade를 고르거나 *그룹당 N건*으로 줄이고, 유전자·구간으로 확대하며, *공통 SNP 숨김*으로 서열 사이의 차이만 볼 수 있습니다. SVG/PNG로 저장합니다.
- **전파 타임라인**(*전파 타임라인* 버튼): 가로축 ISO 주(적으면 자동 월), 행은 **lineage별** 또는 **대표 변이별**(lineage-designation의 정의 SNP; 처음 한 번 *정의 받기*를 누르면 git으로 받습니다). 행을 눌러 나라별 줄을 펼치고, 선 = 처음~마지막 관찰 `(날짜 ── 날짜)`, ◆ = 첫 관찰, 옅은 띠 = 월까지만 있는 날짜. 대표 변이 보기는 lineage 이름과 무관하게 그 SNP를 가진 유전체를 세니 lineage가 새로 생겨도 그림이 유지되고, 두 보기를 나란히 보며 헷갈림을 줄일 수 있습니다. 정의 SNP는 NC_063383.1 좌표이므로 보유 여부는 이 좌표로 따로 정렬해 판정합니다(처음 열 때 1~2분, 결과는 분석 폴더에 보관). **첫 관찰은 유입 시점이 아닙니다.** 변이 그림에서는 서열의 lineage 정의 SNP가 굵게 보입니다.
- **애니메이션과 지도**: 전파 타임라인의 ▶로 주(월)별 재생, 슬라이더로 이동. 커서 시점까지 관찰된 것만 그리고, **세계 지도**에 나라별 원(면적 = 유전체 수, 색 = lineage 구성, 처음 관찰된 나라는 분홍 고리)이 움직이며 아래에 나라별 누적 막대가 있습니다. 지도 위 `표시`에서 lineage(또는 대표 변이) 하나만 골라 그 확산만 볼 수 있습니다. 지도는 앱에 들어 있는 Natural Earth 경계로 그려 인터넷 없이 됩니다. 분석 화면의 표에는 막대·스파크라인이 들어 있습니다.
- **변이 분석**(*변이 분석* 버튼): ① **APOBEC3형 변이도 전파되는가** — 합의 서열로는 *둘 이상이 가진 변이(전달됨)*와 *한 유전체에만 있는 변이*만 구분됩니다. 두 집단의 APOBEC3형 비율을 위치별·가지별(Nextstrain 트리의 끝/안쪽 가지)로 비교하고 검정합니다. 이 데이터셋(Ib, 1,145건, ITR 마스킹·N 5% 이하)에서는 안쪽(물려받은) 가지 변이의 79%가 APOBEC3형(끝 가지 68%, 4가지 계산 모두 p≈0.001~0.002)이라 APOBEC3형 변이는 전달됩니다. 위치 기준(한 유전체에만 62.7% vs 공유 64.3%, p=0.72)에서는 차이가 없었으니 "끝 가지에서 더 적다"는 차이는 작고 방법에 따라 갈립니다(`docs/VALIDITY.md`). **환자 안의 소수 변이(검체 부위마다 다른 변이)가 전파 때 전달되는지는 합의 서열로 알 수 없습니다** — 같은 환자의 여러 부위 검체나 전파 쌍의 원시 리드(iSNV)가 필요합니다. ② **신규 변이 후보** — 정의에 없는 변이를 첫 출현 순으로, 최근 것은 *최근* 표시, 같은 유전체들이 가진 변이 묶음은 새 lineage 후보(확정은 계통수 필요). ③ **lineage 시그니처** — 데이터에서 찾은 시그니처 변이를 공식 정의와 대조.
- **유행 분석**(*유행 분석* 버튼, [Lee, Choi 외 2026](https://doi.org/10.1080/22221751.2026.2746255)의 접근): ① **서브클러스터** — Nextstrain 트리에서 가지 정의 변이가 있는 5건 이상(조절 가능)의 단계통군을 자동 제안하고, 변이(유전자·아미노산·APOBEC3형, *는 다른 가지에도 있는 homoplasy), 클러스터 × 변이 행렬(논문 Fig 2B), 월별 거품 그림(Fig 3B)을 보여 줍니다. 서브클러스터 이름(C1…)은 첫 출현 순으로 붙는 임시 이름입니다. ② **크기 분포**: VMR·Pearson·음이항 θ·AIC/BIC로 포아송과 비교(논문 수치로 검증: 7개 클러스터에서 VMR≈21). ③ **APOBEC3형 변이의 시간 누적**(Fig 4C,D): NC_063383.1 대비 유전체당 개수와 추세. N 5% 이하(236건)에서는 새로 생긴 변이 기준 연 5.3건, lineage 안에서만 연 6.7건으로 알려진 APOBEC3 시계(연 약 6건)와 같은 범위입니다. **N 비율 상한을 풀면 lineage 안의 추세가 사라지니 N 필터를 켜 둔 채 해석하세요.** ④ **변이 유형**별 유전체 막대(Fig 4B,D). ⑤ **유전자**: 유전자별 **독립 돌연변이 사건**(트리 가지마다 한 번)과 길이 대비 기대(포아송, Bonferroni)로 정하는 hotspot(ITR 제외; 합산 변이 수는 상속이 중복되어 계통군 정의 변이가 있는 유전자가 부풀려지므로 판정에 쓰지 않음)·APOBEC3 비율·기능 범주(보충표 1)·**탐색용** dN/dS(선택압 검정 아님; 논문의 KaKs 쌍별 추정과 방법이 다릅니다). 역학 정보(접촉 유형 등)와 저빈도 변이는 넣지 않았습니다. 기능 범주 표는 논문 보충표 1을 옮겨 적은 것이므로 원문과 한 번 대조해 주세요.
- **신규 변이**(*신규 변이* 버튼): 정의(lineage-designation)에 없는 변이를 공유하는 N건 이상의 계통군을 후보로 찾아 *신규*(데이터의 가장 최근 날짜에서 N일 안에 처음 관찰) / *성장 중* / *정체*로 표시합니다. 위쪽 곡선은 정의 외 계통군 변이가 월마다 얼마나 새로 보였는지(붉은 영역은 후보가 최소 크기에 못 미쳐 적게 보이는 검출 지연 구간), 산점도는 변이마다 처음 관찰일과 퍼진 정도입니다. 후보는 가설이며 확정에는 계통수 지지도와 시퀀싱 오류 점검이 필요합니다.
- **변이 묶음(PCA)**(*변이 묶음* 버튼): 함께 나타나는 변이를 유전체 집합의 겹침(Jaccard)으로 묶고(같은 가지의 변이가 한 묶음), 유전체를 PCA로 배치합니다. 데이터셋 1(N 5% 이하 392건)에서는 PC1이 분산의 33%, PC1·PC2가 lineage로 76%·67% 설명되어 묶음의 대부분이 계통 구조입니다. 계통과 다르게 나타나는 점·변이(반복 변이, 새 후보)를 찾는 데 쓰세요. 함께 나타남은 같은 가지에 있다는 뜻이지 기능적 연관이 아닙니다.
- **정렬 도구**: 변이 프로파일의 *정렬* 설정은 기본(auto)이 **원숭이두창은 Squirrel 정렬**, 그 외 바이러스는 minimap2입니다. Squirrel의 정렬 전용 실행(트리 없음, 모든 유전체)을 한 번 돌려 거기서 SNV를 읽으므로 처음에는 몇 분 더 걸리고 이후 같은 조건이면 재사용합니다. N·마스킹·미정렬 위치는 *알 수 없음*으로 다룹니다. 이전 방식과의 비교는 `docs/VALIDITY.md`.
- **변이 히트맵**(*변이 히트맵* 버튼): 변이 그림이 전체를 보는 화면이라면 이쪽은 자세히 보는 화면입니다. 행 = lineage(또는 clade, 트리 서브클러스터), 열 = 변이, 칸 색 = 그 행의 유전체 중 보유 비율, 굵은 테두리 = 정의 변이, 분홍 열 = APOBEC3형. 유전자 막대를 누르면 그 유전자의 변이 목록, 검색·위치 이동·빈도 슬라이더·"정의 변이만" 등 스위치가 있고 칸을 누르면 행별 빈도와 유전자 기능이 나옵니다. *GVF 받기*는 lineage별 GVF를 zip으로 내려받아 VIRUS-MVP 같은 도구에 쓸 수 있게 합니다(VIRUS-MVP를 그대로 붙이지 않은 이유: SARS-CoV-2 전용 데이터, 낡은 의존성 고정, Nextflow 필요).
- 두 분석도 분석 이력에 **보관**되고(멤버 해시·설정·도구 버전·메모·로그), 같은 조건이면 다시 돌리지 않으며(*다시 분석*으로 강제), 오늘의 데이터셋과 비교됩니다. CLI는 없고 화면/API(`POST /api/datasets/{id}/squirrel`, `/mutprofile`)로 실행합니다.

### 13. Pathoplexus 제출 준비 (앱은 업로드하지 않음)

**Pathoplexus 탭**(또는 `python -m genomecurator pathoplexus <바이러스> [--check]`)은 비공개 유전체를 해당 생물의 제출 템플릿과 대조해 줄마다 *통과 / 경고 / 오류*를 보여 줍니다(필수 열, 날짜 형식, 저자 `성, 이름; 성, 이름`, 숙주, 공개 서열과 같은 서열인지, QC). 오류가 없는 행을 골라 **이용 조건**(OPEN 또는 RESTRICTED+만료일)과 그룹 ID를 정하고 *제출 파일 만들기*를 누르면 `data/pathoplexus/<바이러스>/<시각>/`에 `metadata.tsv`, `sequences.fasta`, `submission.json`, `submit.sh`, `README.txt`가 생기고 *zip 받기*로 내려받을 수 있습니다.
- **앱은 올리지 않고 인증 정보를 받지 않습니다.** Pathoplexus 사이트의 *Submit* 화면에서 두 파일을 올려 검토·승인하거나, 본인 터미널에서 `PATHOPLEXUS_TOKEN`을 설정하고 `submit.sh`를 실행하세요.
- **OPEN은 되돌릴 수 없습니다**(승인 즉시 공개, INSDC로 전달). RESTRICTED는 최대 1년이며 기간은 줄이기만 할 수 있습니다.
- 올린 뒤에는 패키지 상태를 *제출함 → 공개됨*으로 직접 바꿔 기록해 두면 비공개 탭에 상태가 표시되고 같은 유전체를 다시 고르지 않게 안내합니다.

### 용어 표준화 (우리 사전)

소스마다 같은 것을 다르게 적습니다(`Homosapiens` / `Homo sapiens`, `Zaire` / `Democratic Republic of the Congo`, 날짜 구간). 그대로 비교하면 가짜 충돌이 되고 큐레이션 표에도 여러 표기가 섞이므로, **우리 사전**으로 들어오는 값을 표준 용어로 바꿉니다. 소스를 가져온 직후, 비교하기 전에 적용합니다.

- **대상(1단계)**: 숙주, 나라, 수집일, clade. 사전 파일은 `genomecurator/resources/vocab/{host,country,clade}.tsv`(표준 용어, 동의어는 `|`로 구분, 숙주 아닌 값에 붙이는 flag)이고, **`data/vocab/`에 같은 이름의 파일을 두면 덮어쓰거나 보탤 수 있습니다**(업데이트해도 남음). 비교는 대소문자·구두점·악센트를 무시합니다. 날짜는 사전이 아니라 규칙입니다.
- **추측하지 않고 지우지 않습니다**: 동의어 → 표준 용어(`mapped`), 모호하지 않은 형식 정리(`normalized`: `15-May-2025`, 한 달 안의 구간 → 그 달, 한 해 안의 구간 → 그 해), 빈 값 표기(`missing`, `null`)는 빈 값(`empty`). **숙주 칸의 세포주·환경 시료, 상위 clade `II`, 재조합 `Ib/IIb`, 일/월 순서가 모호한 `05/06/2025`, 여러 해에 걸친 구간은 그대로 두고 `flagged`**, 사전에 없는 값은 그대로 두고 `unmapped`. 바뀐 값은 원본·표준·규칙이 `standardization` 표에 남습니다.
- **미리보기**: *용어 탭*(또는 `python -m genomecurator vocab <바이러스>`)은 지금 큐레이션된 표에 적용하면 무엇이 바뀌는지 보여 줍니다(저장 없음). 원숭이두창에서는 숙주 55건, 나라 6건, 날짜 83건이 바뀌고 11개 값이 검토 필요입니다. 실제 반영은 **다음 *지금 확인*에서 변경 제안**으로 올라오며 평소처럼 승인해야 합니다.
- **검토 화면(2단계)**: 용어 탭 위쪽 *결정할 값*에 사전에 없는 값과 표시된 값이 필드별로 나옵니다. 가까운 표준 용어 버튼(→)이나 목록으로 **동의어로 추가**, **새 표준 용어로 추가**, 표시된 용어를 **표준으로 인정**, 모호한 날짜에 **표준 날짜 지정**을 누르면 `data/vocab/`에 저장되고 같은 값은 다음 *지금 확인*부터 자동 표준화됩니다(그 변경은 평소처럼 승인). 아래 *사전 보기·고치기*에서 용어·동의어를 검색하고, ✎ 표시된 내가 더한 항목은 삭제할 수 있습니다(삭제하면 내장 사전 상태로 돌아갑니다). 한 표기가 두 용어에 속하게 하는 추가는 거절됩니다.
- **3단계(플랫폼·검체)**: `sequencing_platform`은 **모델을 유지**하고 철자만 표준화합니다(`novaseq` → `Illumina NovaSeq 6000`, `MinION` → `Oxford Nanopore MinION`). 기술은 **`platform` 열**에 따로 채웁니다(`Illumina`, `Oxford Nanopore`; 여러 개면 `Illumina; Oxford Nanopore`). 검체(isolation_source)는 병변 면봉·피부 병변(수포/가피)·혈액 등으로 모읍니다. `clinical`·`patient` 같은 모호한 값은 추측하지 않고 검토 대상으로 남깁니다. 조립 방법은 도구 버전이 정보라 바꾸지 않습니다. 원본 값은 항상 `standardization` 표에 남습니다.
- **다음 단계**: 검증 경고, Pathoplexus 제출 양식과의 연결.

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
| `backup [zip]` / `restore <zip>` | 이사할 때: 큐레이션 결과 + 데이터셋 정의·제외 표시 + FASTA를 담은 이식 가능한 백업 (비공개는 `--with-private`일 때만) |
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
| `GC_MIN_COVERAGE_PCT` | `80` | 커버리지 하한: 기준 게놈(가장 긴 RefSeq) 대비 이 % 미만이면 큐레이션 대상이 아님(제안 안 함, 이미 있으면 `removed`로 미리 승인; NCBI Complete는 예외; 0이면 끔). 비공개 가져오기에도 적용 |
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
scripts/superset.sh init && scripts/superset.sh run          # http://127.0.0.1:8089 (또는 앱의 '분석' 탭에서 켜기, 또는 python -m genomecurator serve --superset)
scripts/superset.sh bootstrap                                # 연결·데이터셋·차트·대시보드 자동 생성
```

자세한 내용은 [docs/SUPERSET.md](SUPERSET.md).

## 공개용 정적 페이지 (읽기 전용)

`docs/index.html` 한 파일에 개요 차트와 유전체 표(필터·정렬·상세·TSV)가 들어 있습니다. 공개 메타데이터만 담기고 서열·경로·이력·수집 코드는 들어 있지 않습니다.

```bash
python scripts/build_static.py     # 로컬 DB로 docs/index.html 다시 만들기
```

앱의 **배포하기** 버튼은 미리보기로 바이러스·유전체 수를 보여 주고 확인을 받은 뒤, 그 파일만 커밋·push합니다. GitHub Pages는 저장소 Settings → Pages → `main` 브랜치 `/docs`로 지정하세요. 차트는 CDN(Chart.js)을 쓰므로 보는 쪽에 인터넷이 필요합니다.

## 이사 가기 (백업 / 복원)

```bash
python -m genomecurator backup my-backup.zip      # 큐레이션 결과 · 결정 이력 · 데이터셋 정의와 제외 표시 · FASTA (원본 응답·대기 중 항목·분석 결과 파일은 제외)
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

* 변경 이력: [CHANGELOG.ko.md](../CHANGELOG.ko.md) (한국어) · [CHANGELOG.md](../CHANGELOG.md) (English)
* 릴리스 방법과 버전 규칙: [docs/RELEASING.md](RELEASING.md)

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

테스트한 버전, 전체 라이선스, 간접 의존성, 데이터를 가져오는 외부 서비스(NCBI·BV-BRC·Pathoplexus·ICTV)는 [docs/DEPENDENCIES.md](DEPENDENCIES.md)에 있습니다
(설치된 패키지에서 `scripts/dependencies.py`가 뽑아 만든 문서라 코드와 어긋나지 않습니다).

## 라이선스 · 만든 사람 · 문의

MIT License. © 2026 CHOICHI-HWAN.
버그와 아이디어는 [Issues](https://github.com/CHOICHI-HWAN/GenomeCurator/issues)로 보내 주세요.

## 같이 쓰면 좋은 도구

* [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) · [BV-BRC](https://www.bv-brc.org/) · [Pathoplexus](https://pathoplexus.org/) — 데이터를 내주는 고마운 곳들입니다. 요청 사이에 간격을 두는 이유입니다.
* [Snippy](https://github.com/tseemann/snippy) — 초기 문서의 구성이 닮아 있다면 그 덕입니다. 변이 호출은 그쪽이 훨씬 잘합니다.
* [Apache Superset](https://superset.apache.org/) — 분석 화면을 맡습니다.


### Nextstrain 정렬의 합의 마스크

원숭이두창에서 Nextclade 정렬에 오류가 있어도 워크플로(`nextclade3 run`)는 그대로 둡니다. 대신 Squirrel 정렬과 **어긋나는 열만 가립니다**.

```bash
python -m genomecurator mask-consensus <squirrel 정렬.aln.fasta.gz> <Squirrel 참조.fa> <nextclade 정렬.fasta> <워크플로 참조.fa> --min-genomes 3 -o mask.bed
```

그 BED를 데이터셋의 *추가 설정*에 `mask: {maskfile: /경로/mask.bed}`로 지정합니다. N·갭은 불일치로 세지 않습니다. clade I 1,145건에서 N≥3이면 14열(0.007%)이 가려집니다. 한계: `mpxv` 분할 빌드와 `QC_rare_mutations`는 계속 Nextclade 기준이고, 입력인 두 정렬 파일은 직접 준비해야 합니다(자동 실행 연결은 아직 없음).


### 독립 변이 호출 (snippy + lofreq)

정렬 기반 호출(Nextclade·Squirrel)과 **독립적인** 두 번째 의견입니다. 환경 만들기: `mamba env create -f variantcall/snippy.yml`, `mamba env create -f variantcall/lofreq.yml` (activation은 conda로).

```bash
python -m genomecurator variantcall 유전체.fasta 워크플로참조.fa --limit 30 --squirrel-aln … --squirrel-ref … --nextclade-aln … -o 변이.tsv
```

snippy `--ctgs`는 조립본을 20x 모의 리드로 자르므로 lofreq가 FASTA에서도 돌아갑니다(깊이 ≈17~20, AF는 0 또는 1). **조립본의 AF·깊이는 원래 시료의 값이 아니라** 소수 변이 판정에는 쓸 수 없습니다. 결과 TSV의 `called_by`가 변이를 찾은 도구이고, 도구들이 어긋나는 변이는 가리지 않고 눈으로 확인할 대상으로 표시합니다. 30건 점검에서 네 도구 일치가 81%였고, snippy/lofreq가 놓친 변이의 상당수는 Squirrel·Nextclade만 부른 것이라 ITR·반복 구간 등 모의 리드 정렬의 한계로 보이지만 원인 분류는 아직 안 했습니다. 아직 화면(변이 목록)에는 연결하지 않았습니다.


### QC: Squirrel과 Nextclade를 나란히

*서열 QC* 화면 아래쪽에 Squirrel QC 카드가 바이러스별로 나옵니다: Squirrel이 N 비율 등으로 **제외를 권한 유전체**, 변이 위치 마스킹 권고(`clustered_snps`, `N_adjacent`, `gap_adjacent`)가 있는 유전체, 그리고 Nextclade **불량**과의 대조(둘 다 / Squirrel만 / Nextclade만)입니다. Nextclade의 *주의* 등급은 원숭이두창 clade I의 대부분이라 비교에서 뺐습니다. 정보만 보여 주며 아무것도 제거하지 않습니다. Squirrel 실행이 끝난 데이터셋이 있어야 나옵니다.

### 변이 호출 일치도 화면

변이 프로파일 결과의 *호출 일치도* 탭에서 **snippy + lofreq 실행**을 누르면 모든 유전체를 두 도구로 다시 호출하고(1,000건 기준 십여 분), 이 실행의 호출과 비교해 도구가 어긋나는 위치를 보여 줍니다. 가리거나 지우지 않으며 어느 쪽이 맞는지 판정하지 않습니다. 사용 전 `variantcall/snippy.yml`·`variantcall/lofreq.yml`로 환경을 만드세요.


### 구현 상태 뱃지

워치리스트와 분석 허브에서 바이러스마다 *수집·큐레이션 / 용어 사전 / 서열 QC / 분석*이 **완료**(초록), **개발중**(노랑), **미지원**(회색)인지 표시합니다. 뱃지에 마우스를 올리면 이유가 나옵니다. 현재 원숭이두창만 완료이고, 다른 바이러스는 수집·큐레이션 흐름만 완료로 두었습니다(아직 큐레이션하지 않았습니다). 기준은 `genomecurator/support.py`에 손으로 적혀 있으며, 기능을 끝낼 때 같이 고칩니다.


### Snakemake 워크플로 (시험 적용, 기존 구현과 병행)

변이 호출(snippy → lofreq)을 `flows/variantcheck/Snakefile`로도 돌릴 수 있습니다. 호출 일치도 화면의 **Snakemake로 실행**을 누르면 되고, 결과 선택에서 *직접 구현*과 *Snakemake*를 바꿔 볼 수 있습니다. 두 결과가 같은지는 화면 맨 위의 “두 구현의 일치”가 알려 줍니다. 환경: `mamba env create -f variantcall/snakemake.yml`(그리고 snippy·lofreq 환경). **차이**: 중간에 멈춰도 다시 실행하면 끝난 유전체는 건너뛰고 이어서 합니다. **방침**: 워크플로가 완성되고 두 구현이 계속 같은 결과를 내면 직접 구현을 deprecated 합니다. 그 전에는 둘 다 유지합니다. Windows에서는 WSL2에서 실행하는 것을 전제로 합니다.


### 설치와 확인: `setup.sh`, `test.sh`

```bash
./setup.sh              # 환경 6개를 mamba로 만듭니다 (--core: 앱만, --only gc-squirrel,gc-snippy: 골라서, --update: 갱신, --dry-run: 계획만)
conda activate genomecurator
./test.sh unit          # 단위 테스트 (인터넷·분석 도구 불필요)
./test.sh check         # 설치가 실제로 동작하는지: 환경·도구 응답·작은 종단 실행(Snakemake 변이 호출이 심어 둔 치환 3개를 찾아야 함)
./test.sh               # 둘 다
```

환경은 mamba로 만들고 활성화는 conda로 합니다. `check`의 **WARN**은 없어도 되는 부분(예: Nextclade가 없으면 서열 통계만 동작)이고 **FAIL**은 필수 부분의 고장입니다. 분석마다 필요한 환경이 다르므로 쓰지 않을 분석의 환경은 만들지 않아도 됩니다.


### 호출 캐시와 레퍼런스 통일

snippy + lofreq 호출은 `data/variantcache/<레퍼런스>-<md5>-<snippy 버전>-<lofreq 버전>/` 아래에 **서열 md5별**로 저장되어 실행 사이에 공유됩니다. 같은 서열은 어느 데이터셋에서든 다시 호출하지 않으므로 **새 데이터를 추가하면 새 유전체만 호출**합니다(그래서 “모든 실행에 돌리기”는 처음 한 번만 오래 걸립니다). 레퍼런스가 같아야 캐시가 맞으므로 클레이드별로 레퍼런스를 하나로 정해 두는 것이 좋습니다: clade I → KJ642613.1, clade II → NC_063383.1(변이 프로파일의 `reference` 설정 `auto`가 이미 그렇게 고릅니다). 레퍼런스나 도구 버전이 바뀌면 새 캐시 폴더가 만들어지고 이전 캐시는 쓰이지 않습니다(지워도 됩니다).


### 유전체 표 → 데이터셋 → 분석 (한 흐름)

1. **유전체** 화면에서 표를 필터하고 행을 고릅니다(선택이 없으면 필터에 보이는 전체).
2. 도구줄의 **📦 데이터셋 만들기** → 이름을 정하면 그 유전체들만 담은 데이터셋이 만들어지고 데이터셋 화면이 열립니다(비공개가 섞여 있으면 비공개 포함이 자동으로 켜집니다).
3. 데이터셋 화면의 **분석** 패널에서 ① Nextstrain 트리 ② Squirrel ③ 변이 프로파일 ④ 변이 호출(snippy+lofreq)을 바로 실행하고, 끝난 것은 같은 줄의 링크로 결과를 엽니다. 단계마다 *현재 데이터와 같은 입력으로 완료*인지 *데이터·설정이 바뀜*인지 보이고, 같은 입력이면 다시 계산하지 않습니다.
4. 규칙으로 정의하는 데이터셋(clade 등)은 **데이터셋** 메뉴에서 만들고, 규칙 데이터셋이든 고른 데이터셋이든 3번은 같습니다. 고급 설정과 이력은 *실행·이력*과 *분석 허브*에 있습니다.

### QC 도구: 서열 통계 · Nextclade · Squirrel · CheckV

**QC 실행**은 서열 통계와 Nextclade(데이터셋이 있을 때)에 더해 최근 끝난 **Squirrel 실행의 판정**을 유전체별로 저장하고(제외 권고 = 불량, 위치 마스킹 권고 = 주의), **CheckV**가 설치되어 있으면 완성도·오염도 평가합니다. CheckV 준비: `mamba env create -f variantcall/checkv.yml` 후 `conda run -n gc-checkv checkv download_database data/checkv`(약 1.5 GB; 환경 안에서 실행해야 마지막 DIAMOND 색인 단계가 됩니다). 데이터베이스가 없으면 CheckV만 건너뛰고 나머지는 돕니다. 모두 정보일 뿐이며 아무것도 제거하지 않습니다.


### 한 번에 분석

데이터셋 화면의 분석 패널 맨 위 **한 번에 분석**에서 단계(QC, Squirrel, 변이 프로파일, 변이 호출, Nextstrain)를 체크하고 **▶ 선택한 단계 모두 실행**을 누르면 하나의 백그라운드 작업이 차례로 실행합니다(진행은 *작업* 탭). 입력이 그대로인 단계는 다시 계산하지 않고, 한 단계가 실패해도 나머지는 계속합니다. Nextstrain은 가장 오래 걸려 마지막입니다. 끝난 결과는 같은 패널의 줄마다 링크로 엽니다. 분석이 준비되지 않은 바이러스는 QC만 실행됩니다.


### QC 우선순위

도구마다 판정이 다를 때 **Squirrel > Nextclade > CheckV**(그다음 내장 서열 통계) 순서로 따릅니다. 유전체 표의 QC 열, 데이터셋의 *QC 제외* 규칙, 내보내기 필터가 쓰는 대표 상태는 *결과가 있는* 가장 우선순위 높은 도구의 판정입니다(오류가 난 도구는 판정이 없는 것으로 칩니다). 예: Squirrel이 양호이고 Nextclade가 불량이면 대표 상태는 양호입니다. 다른 도구의 판정은 QC 화면과 유전체의 상세에 그대로 남으니 의견이 갈리는 유전체를 따로 볼 수 있습니다.


### 결과 검증: `./test.sh verify`

`pytest`는 코드가 맞는지, `./test.sh check`는 도구가 도는지 봅니다. **`./test.sh verify`는 지금 데이터 폴더의 실제 결과가 맞는지** 봅니다(읽기 전용, 10초 안팎). 항목: 데이터 무결성, 용어 사전 멱등성, QC 우선순위, 비공개 유전체 유출 없음, Nextstrain 트리(멤버·clade), Squirrel(정렬·QC 제외·계통수 tip), 변이 프로파일(참조 염기·유전자·아미노산 변화·요약 개수), Nextclade 정렬과의 대조, snippy·lofreq와의 일치, Snakemake 워크플로와 직접 구현의 동일성, lineage 정의 변이 보유율, 지정 변이 타임라인, hotspot, 심은 치환 시험. 줄마다 PASS / WARN(살펴볼 것) / FAIL(나오면 안 되는 결과)이고 FAIL이 하나라도 있으면 종료 코드 1입니다. 입력이 없는 항목은 통과로 치지 않고 SKIP으로 표시합니다.


### 분자시계와 tMRCA (TempEst 방식 → TreeTime → BEAST)

데이터셋 분석 패널의 **분자시계·tMRCA** 버튼(또는 한 번에 분석)이 Squirrel 계통수(없으면 Nextstrain 트리)와 수집일로 두 단계를 실행하고 결과 화면(`/clock/<id>`)을 엽니다(1~2분).

1. **루트-투-팁 회귀**(TempEst과 같은 계산): 산점도에서 직선의 기울기 = 치환 속도, **R²** = 시계다운 정도, 직선이 0과 만나는 날짜 = 루트 날짜. 빨간 점은 시계에서 크게 벗어난 끝점(재조합·오염·날짜 오류·APOBEC3 폭발)이며 회귀에서 뺍니다.
2. **시간 보정 트리**: 모든 노드의 날짜와 **90% 구간**. 그룹(clade·lineage)의 **tMRCA**는 그 유전체들의 공통 조상 날짜입니다. *단계통 아님*은 그 그룹이 트리에서 한 가지로 묶이지 않는다는 뜻(조상 lineage에 흔함)이고, 그때의 tMRCA는 그룹의 시작이 아니라 더 위쪽 조상입니다.

**읽을 때**: 루트 날짜는 수집 기간 밖으로의 외삽이고, R²가 낮거나(0.3 미만) 수집 기간이 짧으면(5년 미만) 구간이 매우 넓습니다. 화면의 경고를 숫자와 함께 읽으세요. clade Ib 데이터는 R² 0.37, 2.9년이라 TreeTime의 루트 2020.4(2019.4~2022.3)는 인용할 값이 아닙니다.

**BEAST**(같은 화면 아래): 인용할 추정은 이쪽입니다. 표본 수·체인 길이·시계(strict/relaxed)·트리 사전분포를 정하고 **패키지만 만들기**(XML과 정렬, 날짜) 또는 **BEAST 실행**(백그라운드, 수 시간)을 누릅니다. 결과는 tMRCA의 중앙값·95% HPD와 **ESS**입니다. **ESS가 200 미만이면 체인이 덜 돈 것이니 구간을 믿지 마세요**(실행 중에도 중간 결과가 나오지만 ESS가 낮습니다). 환경: `mamba env create -f variantcall/beast.yml`(BEAST 10 + BEAGLE). 속도 참고: 40개 유전체에서 약 3,400상태/초(BEAGLE CPU), 유전체가 많을수록 느립니다. 정식 분석은 시드를 바꿔 체인을 여러 개 돌려 Tracer로 확인하고 모델(시계·트리 사전분포)을 비교하세요.


### BETS: BEAST 결과를 믿어도 되나

분자시계 화면의 **BETS**가 **같은 모델을 실제 수집일로 돌린 분석과 모든 유전체를 같은 날짜로 둔 분석** 두 개를 동시에 돌려 로그 베이즈 인수(날짜 있음 − 없음)를 냅니다. **5 이상**이면 날짜가 모델을 뚜렷하게 개선한다는(시계 신호가 있다는) 강한 증거입니다. **−1~1이거나 음수**면 날짜가 도움이 안 되므로 BEAST의 tMRCA는 날짜가 아니라 사전분포가 정한 값입니다 — 그 값을 인용하지 마세요. 추정기는 **GSS(일반화 디딤돌 표집, 가장 정확)**와 PS/SS를 함께 돌리는 것이 기본입니다. GSS는 모델마다 일반 체인을 먼저 돌려 작업 사전분포를 만든 뒤 경로를 걷습니다(그래서 단계가 하나 더 있고 프로세스가 4개까지 동시에 돕니다). 설정은 기본값(경로 50단계 × 단계당 250,000상태)을 줄이지 마세요: 짧게 돌리면 추정기들이 크게 어긋나고 베이즈 인수는 의미가 없습니다(화면이 어긋남을 알립니다). 한 번으로 단정하지 말고 시드를 바꿔 반복해 값이 안정적인지 보세요. 걸리는 시간: 표본 수와 단계 수에 비례하는 수 시간(두 분석이 병렬). 이 다음 단계는 같은 방식으로 시계 2종 × 트리 사전분포 2종의 네 모델을 비교하는 것입니다(아직 구현 전).


### 변이 간트 차트 (시간으로 보는 변이)

### 계통도·성장 (clade·lineage 정의도, Sankey, 성장 이점)

*계통·전파 › 계통도·성장*. **계통도**는 Mpox virus › clade › outbreak › lineage 순서로 이어 그리며 노드 옆 **+k**는 그 노드에서 새로 생긴 정의 변이 수, **Σm**은 뿌리부터 누적된 수입니다. 가로축을 *누적 정의 변이 수*로 두면 변이가 쌓이는 모양이 보입니다. **Sankey**는 띠 굵기가 유전체 수이며 *유전체 수 기준*으로 전체 또는 최근 몇 주를 고릅니다. 노드를 누르면 뿌리부터의 경로(새 변이·누적)와 그 노드의 정의 변이(유전자·아미노산)가 나옵니다.

**성장 이점**은 lineage 비율의 다항 로지스틱 회귀입니다. 읽는 법: 기준(또는 평균) lineage보다 *로그 오즈/주*가 양수이고 95% 구간이 0을 걸치지 않으면 표본에서 그 lineage의 비율이 더 빨리 늘고 있다는 뜻입니다. **전파력의 추정이 아닙니다.** 나라별 표집 차이가 크므로 가능하면 나라를 하나 고르고, 구간이 넓은 lineage는 해석하지 마세요(페이지에 주의사항이 함께 표시됩니다).

**lineage 선택·정의 변이**: lineage·clade는 검색 칩 선택기(건수 표시, *하위 lineage 포함*, 바꾸면 자동으로 다시 그림)로 고릅니다. 해당 lineage를 *정의*하는 변이는 줄에 태그가 붙고, *대표 변이* 필터로 정의 변이만/제외해 볼 수 있습니다.

**히트맵 시간별 보기**: (변이를 누르면 맨 위에 **호출 도구별 비교**(이 실행·snippy·lofreq·Nextclade)가 나옵니다. Nextclade 결과는 호출 일치도 페이지의 *Nextclade 실행*으로 만듭니다.) (그 아래에는 같은 변이의 **전파 타임라인**(나라 × 시간, 누적 곡선; 지도는 뺌)이 이어서 나옵니다.) 히트맵에서 변이를 누르면 행별 빈도 아래에 월/ISO 주별 보유 유전체 수(리니지별 누적)와 표집 대비 비율(점선)이 나옵니다. 막대의 증감은 표집량 때문일 수 있으니 비율 선을 같이 보세요. **변이 PCA**는 *변이 묶음* 페이지의 유전체 PCA 아래에 있습니다.

변이 그림 위쪽 **보기**에서 *시간 간트*를 고르면 **한 줄 = 변이, 가로 = 시간(월 또는 ISO 주)**로 바뀝니다. 칸은 그 기간에 수집된 보유 유전체가 많을수록 진하고, 회색 막대는 처음~마지막 관측, ◆는 처음 관측입니다. 줄은 **처음 관측된 시기 순**(먼저 나온 변이가 위)이 기본이며 *위치 순*, *보유 유전체 수 순*으로도 정렬합니다. 필터(lineage, clade, 유전자·구간, APOBEC3형만, 공통 변이 숨김)는 유전체별 그림과 같고 *최소 보유 수*와 *최대 행 수*가 더 있습니다. 변이를 우클릭하면 히트맵·전파 타임라인으로 갑니다.

**읽을 때**: *처음 관측*은 표본에서 가장 이른 수집일이지 변이가 생긴 시점이 아닙니다. 표본 추출이 시간·장소에 따라 고르지 않아서, 일찍 생긴 변이가 늦게 시퀀싱되면 뒤에 나옵니다. 수집 기간의 첫 유전체부터 있던 변이는 실제 나이와 상관없이 맨 위에 모입니다(그래서 거의 모두가 가진 변이는 기본으로 숨깁니다). 날짜가 월까지만 있는 유전체는 월 단위 보기에서만 칸에 들어가고 그 달의 시작으로 칩니다.
