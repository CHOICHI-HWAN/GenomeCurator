# 라이선스·이용 정책 검토 / Licences and data-use review

검토일 2026-10-19. 코드·라이브러리 라이선스는 설치된 패키지와 각 저장소의 LICENSE 파일에서 확인했고, 데이터 이용 정책은 각 기관의 공식 페이지에서 확인했습니다. 법률 자문이 아니며, 정책은 바뀔 수 있으니 공개·인용 전에 원문을 확인하세요.
Reviewed 2026-10-19 from the installed packages, the projects' LICENSE files and the providers' official pages. This is not legal advice; policies change, so check the originals before publishing or citing.

## 결론 / Verdict

| 질문 | 답 |
|---|---|
| 이 저장소의 **코드**를 GitHub에 공개해도 되나? | **예.** 쓰는 라이브러리는 모두 MIT·BSD·Apache·PSF·MPL 등 허용적 라이선스이고, 이 프로젝트는 MIT입니다. 강한 카피레프트(GPL·AGPL)를 코드에 끌어들이는 의존성은 없습니다. |
| **Docker 이미지**를 GHCR에 올려도 되나? | **예.** 베이스 이미지(miniforge3: Ubuntu + conda-forge)에 시스템 라이브러리(GPL/LGPL: readline, libgcc 등)가 들어 있으나, 별도 프로그램으로 동봉될 뿐 이 프로젝트 코드와 결합되지 않습니다. 일반적인 컨테이너 배포 관행과 같습니다. |
| **공개 페이지**(`docs/index.html`)에 데이터를 올려도 되나? | **조건부로 예.** 아래 "데이터" 표를 보세요. NCBI는 자체 제한이 없고, Pathoplexus는 *공개(Open)* 데이터만 들어가며, BV-BRC는 명시적 재배포 조항을 찾지 못해 출처 표기를 권합니다. 이번 검토로 페이지 푸터에 출처·면책 고지 링크를 넣었습니다. |
| **Nextclade 데이터셋**을 저장소에 넣어도 되나? | **넣지 않습니다.** 데이터셋에는 라이선스 파일이 없어 재배포 조건이 불명확합니다. 앱이 필요할 때 Nextstrain 서버에서 내려받아 `data/`(git 제외)에만 둡니다. |

## 1. 소프트웨어 / Software

### 실행 (pip / conda)
| 구성요소 | 라이선스 | 비고 |
|---|---|---|
| FastAPI, Uvicorn, Jinja2, HTTPX, APScheduler, python-multipart | MIT / BSD-3 / Apache-2.0 | 허용적 |
| 간접 의존성 16개 | MIT 8 · BSD 5 · Apache 1 · PSF 1 · **MPL-2.0 1 (`certifi`)** | `certifi`의 MPL-2.0은 파일 단위 카피레프트로, 수정 없이 의존성으로 설치만 하므로 문제 없습니다. |
| **Nextclade CLI** | **MIT** (저장소 LICENSE 확인) | `environment.yml`에 포함(bioconda). 별도 프로그램을 실행(subprocess)하는 방식이라 결합되지 않으며, 없으면 서열 통계만 동작합니다. |
| conda 기본 라이브러리 | GPL-3.0(+GCC 예외)·LGPL: `readline`, `libgcc`, `libgomp`, `ld_impl` 등 | 시스템 런타임입니다. 코드에 링크하지 않습니다. |
| 브라우저의 Chart.js 4.4.1 (jsDelivr CDN) | MIT | 공개 페이지는 CDN에서 불러옵니다(저장소에 포함 안 됨). |

### Nextstrain 분석 (별도 환경, 저장소에 포함하지 않음)
| 구성요소 | 라이선스 | 비고 |
|---|---|---|
| nextstrain/mpox 워크플로 | MIT | `git clone`으로 `data/nextstrain/mpox`에 받습니다(저장소에는 없음). 분석마다 커밋 해시를 기록합니다. |
| nextstrain-cli, snakemake, TreeTime, Nextclade, Node.js | MIT | 별도 conda 환경 `gc-nextstrain` |
| MAFFT | BSD | 〃 |
| **augur, auspice** | **AGPL-3.0** | 〃. 수정하지 않고 **별도 프로그램으로 실행**만 하며 이 저장소·Docker 이미지에는 들어가지 않습니다. Auspice 보기 서버는 `127.0.0.1`에만 열립니다. **다른 사람이 쓰도록 네트워크에 공개**하면 AGPL의 소스 제공 의무(수정하지 않았다면 업스트림 소스 링크)가 그 운영자에게 생기니, 공개하려면 따로 검토하세요. |
| **IQ-TREE** | **GPL-2.0+** | 〃 (별도 프로그램으로 실행) |
분석 결과(트리 JSON)와 입력 표에는 공개 데이터의 출처 조건과 *비공개* 데이터의 공유 제한이 그대로 따라갑니다: 비공개 포함 분석 결과(`data/nextstrain/auspice/`, `runs/`)는 공유하지 마세요.

### 선택 (개발·시연·분석)
| 구성요소 | 라이선스 | 비고 |
|---|---|---|
| pytest | MIT | 테스트 |
| Pillow | MIT-CMU(HPND) | GIF 조립 |
| websocket-client | Apache-2.0 | 헤드리스 Chrome 조종 (`DEPENDENCIES.md`에는 현재 환경에 없어 "측정 안 됨"으로 표시) |
| Apache Superset | Apache-2.0 | **별도 conda 환경**에서 별도 프로세스로 실행하며 저장소에는 코드가 없습니다(`superset/` 설정 파일만). |
| Google Chrome | 독점 | GIF 녹화에만 사용, 배포하지 않습니다. |

### CI·배포 도구
GitHub Actions(`checkout`, `setup-python`, `setup-miniconda`, `setup-qemu`, `setup-buildx`, `build-push`, `login`)는 MIT/Apache-2.0이며 서비스로 쓰는 것이라 배포물에 들어가지 않습니다.

## 2. 데이터 / Data

| 출처 | 조건 (공식 페이지 기준) | 이 도구의 상태 | 해야 할 일 |
|---|---|---|---|
| **NCBI** | 분자 데이터(서열 등)에 대해 NCBI 자체는 사용·배포 제한을 두지 않음. 단, 제출자가 특허·저작권 등을 주장했을 수 있고 NCBI가 판단·허락하지 않음. NCBI 데이터를 쓰는 서비스는 NCBI 면책·저작권 고지를 사용자에게 보여야 함. E-utilities는 초당 3건 이하. | 요청 간 300ms 간격(키 있으면 더 빠르게). 공개 페이지에는 메타데이터만, 서열은 없음. | 페이지 푸터에 NCBI 정책·면책 링크를 넣음(반영됨). README에도 안내. |
| **Pathoplexus** | *공개(Open)* 데이터: 사용 제한 없음, 출처 표시 권장·재공유 허용. *제한 사용(Restricted-Use)*(최대 1년): 재배포 시 `dataUseTerms` 열 유지·링크 등 조건. | 제한 사용 데이터는 **INSDC에 아직 공개되지 않아 accession이 없고**, 이 도구는 NCBI가 정한 accession 집합에만 값을 보강하므로 **들어오지 않습니다**(mpox 확인: 제한 1,644건 중 INSDC accession 보유 0건). | 푸터에 "공개 데이터만" 명시(반영됨). 인용 시 accession 나열. |
| **BV-BRC** | 공식 문서에서 재배포 조항을 찾지 못함(접근 제한으로 약관 페이지는 읽지 못함). 인용을 요청: Olson RD et al., *Nucleic Acids Res* (PMID 36350631). | 낮은 우선순위의 보강 소스. 일부 필드 값이 BV-BRC에서 올 수 있음. | **BV-BRC 약관을 직접 확인**하거나 문의하세요. 출처 표기(푸터)와 인용을 유지합니다. |
| **ICTV** | 사이트 푸터에 CC BY 4.0 표기(VMR은 Zenodo 기록 확인). | 이름 확인용으로 EBI OLS4를 호출할 뿐 데이터를 저장·재배포하지 않음. | 필요 없음(논문에서 쓰면 인용). |
| **Nextclade 데이터셋** | 데이터셋 파일에 라이선스가 없음. Nextstrain 매뉴얼은 데이터 제공자가 공개했더라도 자유 이용 허가는 아니니 인용하라고 안내(예: mpox). | 저장소에 포함하지 않고 앱이 내려받음. 공개 페이지에는 Nextclade 결과를 싣지 않음. | 논문·발표에서 쓸 때 Nextclade와 데이터셋 저자를 인용(Aksamentov et al., *JOSS* 2021, doi:10.21105/joss.03773). |

### 비공개 데이터 (사용자 소유)
`비공개` 저장소의 유전체는 이 도구의 공개 경로(공개 페이지 배포, Superset 뷰, TSV 다운로드, 일반 백업)에 들어가지 않도록 `genome` 표와 분리되어 있습니다. 내보내기에서 *비공개 포함*을 켠 폴더와 `backup --with-private`로 만든 백업에만 들어가며, 이 둘은 공유·커밋하면 안 됩니다(내보내기 폴더에는 `.gitignore`와 안내 파일이 자동으로 생깁니다).

## 3. 공개(GitHub Pages) 전 점검표 / Before publishing

- [x] 코드는 MIT, 저장소에 비밀값 없음 (`data/`, `superset/home/` 제외; 추적 파일·이력에서 토큰·키 패턴 검색 결과 없음).
- [x] 서열·로컬 경로·이력·수집 코드는 공개 페이지에 들어가지 않음.
- [x] 출처·면책 고지가 페이지 푸터에 있음.
- [x] Pathoplexus 제한 사용 데이터는 수집 구조상 포함되지 않음.
- [ ] BV-BRC 약관 직접 확인 (아직 못 함).
- [ ] 공개 페이지(`docs/index.html`)는 약 18MB이고 배포할 때마다 커밋되어 저장소 이력이 커집니다. 자주 배포한다면 파일 분할·압축, 또는 릴리스 자산으로 옮기는 것을 고려하세요.
- [ ] 논문·보고서에 쓸 때: NCBI, Pathoplexus(accession·SeqSet), BV-BRC, Nextclade와 데이터셋 저자를 인용하세요.
