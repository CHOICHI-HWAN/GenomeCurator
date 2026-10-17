# 릴리스 절차 / Releasing

## 한국어

### 버전 규칙
[유의적 버전](https://semver.org/lang/ko/) `MAJOR.MINOR.PATCH`를 따릅니다. 지금은 `0.x`(초기 단계)입니다.

| 올리는 자리 | 언제 |
|---|---|
| **PATCH** (0.1.**1**) | 버그 수정만. 동작·화면·DB 구조는 그대로. |
| **MINOR** (0.**2**.0) | 새 기능, 눈에 띄는 동작 변화. `0.x`인 동안은 호환이 깨지는 변화도 여기서 합니다 (변경 이력의 `변경`·`데이터·호환성`에 꼭 적습니다). |
| **MAJOR** (**1**.0.0) | 쓰임새와 DB 구조가 안정되었다고 판단할 때 첫 1.0.0. 이후 호환이 깨지면 MAJOR. |

**DB 스키마 버전**(`genomecurator.DB_SCHEMA_VERSION`)은 앱 버전과 별개입니다. 오래된 코드가 읽을 수 없는 DB 구조 변경을 할 때만 올립니다.
지금까지의 변경은 모두 `db.migrate`의 추가형 마이그레이션이라 1입니다. 올렸다면 변경 이력의 `데이터·호환성`에 업그레이드 방법을 적으세요.

### 버전은 어디에 있고 어디에 기록되나
- **정하는 곳은 한 곳**: `genomecurator/__init__.py`의 `__version__`. `pyproject.toml`은 그 값을 읽습니다(`dynamic`). 다른 곳에 손으로 적지 마세요.
- 보이는 곳: `genomecurator --version`, 화면 하단, `/api/version`, 로그 첫 줄.
- 기록되는 곳: DB(`app_meta`, 버전이 바뀔 때마다 이력의 `version`), 각 changeset의 `summary.app_version`, 백업의 `manifest.json`(`app_version`, `db_schema_version`).

### 평소: 변경은 [Unreleased]에
기능을 넣거나 버그를 고칠 때마다 `CHANGELOG.md`(영문)와 `CHANGELOG.ko.md`(한글)의 `## [Unreleased]` 아래에 `Added / Changed / Fixed / Data & compatibility / Dependencies` 중 맞는 항목으로 **두 파일 모두** 적습니다. 사용자가 알아야 하는 것만 적고 커밋 메시지를 베끼지 마세요.

### 릴리스할 때 (예: 0.2.0)
1. 두 변경 이력 파일에서 `## [Unreleased]` 아래 내용을 `## [0.2.0] - YYYY-MM-DD`로 옮기고, 맨 위에 빈 `## [Unreleased]`를 다시 둡니다. 맨 아래 링크도 갱신합니다
   (`[Unreleased]: .../compare/v0.2.0...HEAD`, `[0.2.0]: .../releases/tag/v0.2.0`).
2. `genomecurator/__init__.py`의 `__version__`을 `"0.2.0"`으로. DB 구조가 호환되지 않게 바뀌었으면 `DB_SCHEMA_VERSION`도 올립니다.
3. 의존성이 바뀌었으면 `pyproject.toml`·`environment.yml`·`superset/requirements.txt`를 맞추고 `python scripts/dependencies.py --superset-python <gc-superset의 python>`으로
   `docs/DEPENDENCIES.md`를 다시 만든 뒤 변경 이력의 `의존성`에 적습니다.
4. `python -m pytest -q` — `test_changelog_matches_the_code_version`이 변경 이력과 코드 버전이 맞는지(날짜, 순서, 링크, 영문·한글 일치) 검사합니다.
5. 커밋하고 주석 태그를 만들어 올립니다.
   ```bash
   git commit -am "Release v0.2.0"
   git tag -a v0.2.0 -m "GenomeCurator 0.2.0"
   git push origin main v0.2.0
   ```
   태그를 올리면 그 커밋까지의 모든 커밋이 함께 올라갑니다. 워크플로 파일(`.github/workflows/`)이 바뀐 커밋이 섞여 있으면 토큰에 `workflow` 권한이 필요합니다.
6. (선택) GitHub의 Releases에서 `v0.2.0` 태그로 릴리스를 만들고 변경 이력의 해당 절을 붙여 넣습니다.

### 되돌리기·주의
- 이미 올린 태그는 고쳐 쓰지 말고 다음 PATCH를 냅니다. 틀린 변경 이력은 다음 릴리스에서 `Fixed`가 아니라 `Changed`로 바로잡으세요.
- 데이터가 들어 있는 `data/`와 백업은 저장소에 올라가지 않습니다. 업그레이드 전에는 `genomecurator backup`을 권합니다.

## English

### Versioning
[Semantic Versioning](https://semver.org/) `MAJOR.MINOR.PATCH`. The project is in its `0.x` (initial) phase.

| Bump | When |
|---|---|
| **PATCH** (0.1.**1**) | Bug fixes only. No change to behaviour, screens or the database layout. |
| **MINOR** (0.**2**.0) | New features and noticeable behaviour changes. While the version is `0.x`, breaking changes also go here (always list them under `Changed` / `Data & compatibility` in the changelog). |
| **MAJOR** (**1**.0.0) | The first 1.0.0 when usage and the database layout are judged stable; after that, a breaking change means a MAJOR bump. |

The **database schema version** (`genomecurator.DB_SCHEMA_VERSION`) is separate from the app version. Raise it only for a database layout change that older code
cannot read. Every change so far is an additive migration in `db.migrate`, so it is 1. If you raise it, explain the upgrade under `Data & compatibility`.

### Where the version lives and where it is recorded
- **One source**: `__version__` in `genomecurator/__init__.py`. `pyproject.toml` reads it (`dynamic`). Never type it anywhere else.
- Shown by: `genomecurator --version`, the page footer, `/api/version`, the first log line.
- Recorded in: the database (`app_meta`, plus a `version` entry in the history whenever it changes), each changeset (`summary.app_version`), and every backup (`manifest.json`: `app_version`, `db_schema_version`).

### Day to day: write changes under [Unreleased]
Whenever you add a feature or fix a bug, add an entry under `## [Unreleased]` in **both** `CHANGELOG.md` (English) and `CHANGELOG.ko.md` (Korean), in the matching group
(`Added / Changed / Fixed / Data & compatibility / Dependencies`). Write what a user needs to know, not a copy of the commit message.

### Making a release (e.g. 0.2.0)
1. In both changelog files, move the `## [Unreleased]` content under a new `## [0.2.0] - YYYY-MM-DD` heading and leave an empty `## [Unreleased]` on top. Update the links at the bottom
   (`[Unreleased]: .../compare/v0.2.0...HEAD`, `[0.2.0]: .../releases/tag/v0.2.0`).
2. Set `__version__ = "0.2.0"` in `genomecurator/__init__.py`. If the database layout changed incompatibly, raise `DB_SCHEMA_VERSION` too.
3. If dependencies changed, update `pyproject.toml`, `environment.yml` and `superset/requirements.txt`, regenerate `docs/DEPENDENCIES.md` with
   `python scripts/dependencies.py --superset-python <python of gc-superset>`, and note it under `Dependencies` in the changelog.
4. Run `python -m pytest -q`. `test_changelog_matches_the_code_version` checks that the changelog and the code agree (dates, order, links, English/Korean in step).
5. Commit, create an annotated tag, push both.
   ```bash
   git commit -am "Release v0.2.0"
   git tag -a v0.2.0 -m "GenomeCurator 0.2.0"
   git push origin main v0.2.0
   ```
   Pushing a tag pushes every commit up to it. If a commit that changes a workflow file (`.github/workflows/`) is among them, the token needs the `workflow` scope.
6. (Optional) Create a GitHub Release from the `v0.2.0` tag and paste the matching changelog section.

### Notes
- Never rewrite a tag that has been pushed; ship the next PATCH instead. Correct a wrong changelog entry in the next release under `Changed`, not by editing history.
- `data/` and backups are not in the repository. Run `genomecurator backup` before upgrading.
