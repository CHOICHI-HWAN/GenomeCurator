# Changelog

All notable changes to GenomeCurator are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/) (while the major version is 0, a minor bump may change behaviour;
a patch bump only fixes bugs). How a release is made: [docs/RELEASING.md](docs/RELEASING.md).
**한국어:** [CHANGELOG.ko.md](CHANGELOG.ko.md)

How to read an entry: `Added` new capabilities, `Changed` behaviour that existed before, `Fixed` bugs, `Data & compatibility`
anything that touches the database, backups or curated data (so you know before you upgrade), `Dependencies` libraries added, removed or
re-pinned (the full list with licences is in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)).

## [Unreleased]

### Added
- **Table structure board** on the **Schema** tab: the fields of NCBI (nested objects kept as groups), BV-BRC and Pathoplexus on the left, the columns of
  the local `genome` table on the right, joined by lines that show which source field fills which column. The links are *measured*, not declared:
  the real adapters are run on a row holding one field at a time, and the columns that change are the ones it feeds (so fallbacks such as
  "A, else B" and one field feeding two columns are found too). Vanished fields turn red, the columns they fed turn orange (a warning when
  another source still fills the column, an error when none does), new fields blue; hovering highlights a connection, clicking shows the source
  fields behind a column in priority order, with type, fill rate and an example. A source whose response shape changed gets a red banner.

### Changed
- The Schema tab now starts with the full check and its board; the daily per-virus check moved below it.
- Vanished-field detection also considers the fields the code read at the last good baseline. What the code reads depends on the data
  (an absent `host` object hides `host.organism_name`), so a vanished parent used to hide its children.

### Fixed
- The schema check reported a fill rate above 100% for lists of objects (such as the lineage of an NCBI record) because each list element was counted;
  a record now counts once.

## [0.1.0] - 2026-10-17

First release. It collects virus genomes from three public sources, lets a person review and approve every change, and checks
that the sources agree with each other, with the sequences, and with the code's own expectations.

### Added

**Collecting and merging**
- Watchlist: register a virus by name. The name is first looked up in ICTV, NCBI Virus, BV-BRC and Pathoplexus; names that are not
  viruses are refused with candidate suggestions. Watches can be edited and deleted (typed confirmation and an impact preview).
  A strain or sub-species is widened to its species so records tagged at the species node are not missed.
- Three sources. **NCBI Virus** (Datasets v2) defines which accessions exist; **BV-BRC** and **Pathoplexus** enrich them.
  Pathoplexus organisms are discovered from the site itself and matched by NCBI taxon (including descendants); organism-specific
  fields (clade, lineage, serotype, outbreak, genotype, subtype) are read automatically.
- Field priority is Pathoplexus > NCBI Virus > BV-BRC, and the source that supplied each field is recorded (`field_provenance`).
- Columns include SRA, BioProject, BioSample, sequencing platform and assembly method, with NCBI links.
- Every web request is spaced at least 300 ms apart and retried with exponential backoff (up to 3 times).

**Review and apply (changesets)**
- `check` never changes local data. It produces a *changeset* of new / changed / removed records (a changed record means a metadata
  change or a higher GenBank `.version`). Nothing becomes local data until a person approves it; rejected items do not come back
  when their content is unchanged.
- Approve / Reject / Hold in a spreadsheet-style review grid, including selections above 1,000 rows and a one-click approve + apply.
- A dry-run precedes every apply: free disk space, links, and sample sequences. Apply is atomic: files are backed up first and
  everything is rolled back on an exception; individual failed accessions are skipped and reported; interrupted downloads resume from
  staged files.
- Metadata-first apply (FASTA can be fetched later with *Get FASTA*), and no FASTA download for metadata-only changes.
- Changesets can be deleted (impact preview, confirmation). The review screens are written in plain language with step-by-step guidance.

**Data integrity**
- Cross-source reconciliation at every check: each field is classified as *match*, *cosmetic* (same fact written differently, e.g.
  `USA: CA` vs `USA, CA`, `2022-07` vs `2022-07-11`) or *conflict*, plus coverage (which source has which accession). Shown on the
  **Integrity** tab and as a column, filter and warnings during review.
- Sequence-body verification: sequence bodies (header and whitespace removed, upper-cased, SHA-256) from NCBI, BV-BRC and Pathoplexus
  are compared for accessions whose length, version or hash disagree, with NCBI as the reference. A dry-run plan comes first (light
  length queries decide most cases without downloading), results are resumable, and a source that holds NCBI's exact sequence resolves
  that pair's length/version conflict.
- `verify`: FASTA length and md5 against the records, dates, orphan rows, and more.
- A source that fails to respond now stops the check instead of being skipped, because merging without it could blank curated fields.

**Source schema monitoring**
- Schema dry-run (`schema-check`, **Schema** tab): the fields the adapters really read are measured by running them over live samples,
  then compared with the previous baseline to find vanished, retyped and new fields and changed response shapes, and the local database
  is checked against what the code expects. The report says whether the source code needs an update and where.
- A baseline that shows problems is not replaced until a person accepts the new state.
- A daily per-virus check (3 records per source per virus, 03:00) with a badge on the **Schema** tab when something needs attention.

**Web interface**
- Spreadsheet-style grids with facet filters, sorting, column picker, TSV export and maximize; a detail panel with per-value source
  badges and Pathoplexus variant details, whose width can be dragged; dashboard charts; one consistent layout.
- A progress modal with a large bar and a server-kept job log that survives page reloads; a history page with the server log viewer.

**Analysis and publishing**
- Apache Superset integration: a read-only flat view (`v_genome_analysis`), `scripts/superset.sh`, and 10 charts plus a dashboard
  created through its API.
- A read-only public page (`docs/index.html`, also served locally at `/public`) and a one-click publish to GitHub Pages.

**Operation**
- Windows, macOS and Linux: OS-independent data folder, UTF-8 / LF files, portable relative FASTA paths, safe folder names, and
  launch scripts that prefer a conda/mamba environment (`environment.yml`, plus `superset/environment.yml`). CI runs a three-OS matrix.
- Portable backup and restore (`backup` / `restore`): curated data, decisions and FASTA; refuses a backup from a newer database layout.
- Structured logging: a rotating log file with job / watch / changeset context, every HTTP attempt (secrets masked), and failures kept
  in the audit log even when the failed work is rolled back.
- Versioning: the version has one source (`genomecurator.__version__`); `genomecurator --version`, the page footer and `/api/version`
  show it; the database, every changeset and every backup record which version produced them.

**Documentation**
- A rewritten README with demo GIFs that a script records from the real application (`scripts/demo/`), this changelog, the release process
  ([docs/RELEASING.md](docs/RELEASING.md)) and the list of libraries used ([docs/DEPENDENCIES.md](docs/DEPENDENCIES.md), generated from the installed packages).

### Dependencies
- Python **3.10 or newer** (tested on 3.10 and 3.12; the conda environment uses 3.12).
- Runtime: `fastapi>=0.110`, `uvicorn>=0.29`, `jinja2>=3.1`, `httpx>=0.27`, `apscheduler>=3.10,<4`, `python-multipart>=0.0.9`.
- Optional: `pytest>=8` (tests); `pillow>=10` and `websocket-client>=1.6` (recording the README GIFs, plus Google Chrome);
  Apache Superset `6.1.0` with `rich`, `cachetools` and `flask-caching<2.4` (analysis, separate environment).
- In the browser: Chart.js 4.4.1 from the jsDelivr CDN (charts need internet).
- Data come from NCBI Datasets, BV-BRC, Pathoplexus (LAPIS) and, for name checks, ICTV through EBI OLS4.

### Data & compatibility
- Database layout (schema version **1**): all changes so far are additive migrations applied automatically when the app opens a
  database, so a database from an earlier development build upgrades in place. The first time 0.1.0 opens a database it records
  itself in the audit log (`version`) and `app_meta`.
- Backups made by this version carry `app_version` and `db_schema_version` in `manifest.json`.

### Known limitations
- The daily schema check samples the *first* records NCBI and Pathoplexus return (older records), so a change that only appears on
  new records can be noticed late; run the full check from the **Schema** tab now and then. BV-BRC is sampled newest-first.
- Identical sequence bodies prove two records are the same entity; they do not prove that host, date or other metadata is correct.
- The scheduler (daily checks) runs only while the server is running.

[Unreleased]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/releases/tag/v0.1.0
