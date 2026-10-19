# Changelog

All notable changes to Virus Genome Metadata Curator (package `genomecurator`, formerly shown as "GenomeCurator") are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/) (while the major version is 0, a minor bump may change behaviour;
a patch bump only fixes bugs). How a release is made: [docs/RELEASING.md](docs/RELEASING.md).
**한국어:** [CHANGELOG.ko.md](CHANGELOG.ko.md)

How to read an entry: `Added` new capabilities, `Changed` behaviour that existed before, `Fixed` bugs, `Data & compatibility`
anything that touches the database, backups or curated data (so you know before you upgrade), `Dependencies` libraries added, removed or
re-pinned (the full list with licences is in [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md)).

## [Unreleased]

### Added
- **Nextstrain analysis, connected**: run the Nextstrain mpox workflow (`phylogenetic/`) on a dataset from the **Nextstrain** tab (or `python -m genomecurator analyze <virus> <dataset>`). The tool exports the dataset straight into the workflow checkout, writes the config overlay, passes it after the build's own defaults (the Snakefile only loads those when no config is given), runs `nextstrain build --ambient` as a server job with step progress and the log, and keeps the result. **The dataset's own Nextstrain settings** are a YAML box (nested keys such as `filter.min_length` merge into the workflow's config, the effective merged config is shown before running; inputs, names and the private files are managed by the tool). The tools live in a separate conda environment (`nextstrain/environment.yml`, name `gc-nextstrain`: nextstrain-cli, augur, snakemake, nextclade, iqtree, mafft, auspice); the tab checks it, clones / updates the workflow with a button and records the workflow commit and tool versions per run.
- **Every analysis is kept and can be chosen**: each run is recorded (members hash and counts, rules, settings, workflow commit, tool versions, note, log) with its Auspice files under `data/nextstrain/auspice/` and its input tables / manifest / merged config under `data/nextstrain/runs/`. Running again never replaces an older result; a run with exactly the same members, settings and workflow is not repeated unless forced. The runs table shows for each past run **how it relates to the dataset today** (same members, or +added / −removed, different settings, different workflow), so you can tell which result matches the current data. A local **Auspice** server (`auspice view`, 127.0.0.1 only) opens any kept run by name. QC-less members (dropped by the workflow's default filter) are announced before the run.
- `exclude` accepts a Pathoplexus accession (the id Nextstrain trees show) as well as the INSDC accession.
- **QC for a dataset**: *Run dataset QC* (the Dataset View, the Datasets editor, `POST /api/datasets/{id}/qc`) runs sequence statistics and Nextclade for exactly the dataset's members - the public ones with a sequence file and, if the dataset includes them, the private ones - not the whole virus; unchanged sequences are skipped (*Recalculate* forces it). The Dataset View shows a QC summary card (good / mediocre / bad per tool, how many members have not been through QC yet) and the table's QC column and filter pick out the genomes that need a look.
- **Browse a dataset like the genome table, and see it on the dashboard**: every dataset has a *Dataset View* (`/dataset/<id>`) - the same table with filters, sorting, columns, the detail panel with value provenance and TSV, restricted to the dataset's members (members without a sequence are shown and counted; private members, when the dataset includes them, appear marked *private* with a *Source* column and filter). The **Dashboard** has a dataset selector: all charts, KPIs and the click-through to the table then cover only that dataset. Dataset cards and the editor link to both. The genome table tells how many private genomes exist (they are kept apart on purpose), and the dataset preview always shows how many private genomes match the rules and whether they are included.
- **Datasets (ad hoc, named, saved)**: define which genomes belong to an analysis - e.g. only clade Ib - save it under a name, preview it live and export it again and again. A dataset is a *definition*: rules (clade, lineage, outbreak, region, country, host, completeness, collection years, length range, QC statuses to leave out, sequence required) plus manual adds / removals and an option to include the private store; the global *exclude from dataset* marks still apply. The **Datasets** tab builds the rules from the values that exist in the data (with counts and how many have a sequence), shows the member count with breakdown by clade / country / year and warnings as you click, and has presets for mpox (clade Ib, Ia, I, IIb, complete genomes). Export writes `exports/<virus>/<dataset>/<timestamp>/` and the manifest records the dataset, its rules and a **hash of the exact member ids** (the export history keeps them), so a result stays traceable even after the rules change; the dataset can name the nextstrain/mpox build (`mpxv`, `hmpxv1`, `clade-i`) that the generated config suggests. `dataset-create` / `dataset-list` / `dataset-export` / `dataset-delete` on the CLI; `GET/POST /api/datasets`, `/api/datasets/{id}/export`, and `?dataset=<id>` on `/api/export/{virus}/metadata.tsv|sequences.fasta` for other programs (downloads never carry private genomes).
- **Leave genomes out of the dataset - without deleting anything**: select genomes in the genome table (or a row in the Private tab) and press *Exclude from dataset* with an optional reason; *Undo exclusion* puts them back. It is only a mark (`dataset_exclusion`): the genome, its sequence, history and QC all stay. Exports skip marked genomes, the preview says how many, `excluded.tsv` (and `excluded-private.tsv` with a private export) records who / why / when, and the manifest counts them. A *Dataset* column and filter in the genome table show the state; CLI `exclude <virus> <ids> [--reason ..] [--private] [--undo] [--list]`; API `GET/POST /api/exclusions/{virus}`. Separate from `Removed` (dropped by NCBI or below the coverage floor): that is about the curated set, this is a per-dataset decision.
- **Private store for final analysis datasets**: local, unpublished genomes are imported (`private-check` dry run, then `private-import`; or the **Private** tab) from a Nextstrain-style metadata TSV + FASTA into a store that is **separate from the public set** - they never appear in the genome table, the public page, the Superset view, TSV downloads, `check`, or a general backup. Structural problems (missing / duplicate ids, ids that collide with public accessions, rows without a sequence, sequences without a row, non-nucleotide characters) are errors that stop the import; dates, region and length are normalised and reported. QC (sequence statistics, Nextclade) runs on them (`qc <virus> --private`) and a supplied clade that differs from Nextclade's is flagged.
- Private genomes go in as a **spreadsheet** (`.xlsx`, `.tsv` or `.csv`; a ready-to-fill workbook is one click away: a `genomes` sheet, a `열 설명` sheet explaining every column and an example sheet; the clade column is a drop-down for mpox) **or in a real in-page spreadsheet** ([Jspreadsheet CE](https://github.com/jspreadsheet/ce), MIT, loaded from the jsDelivr CDN: continuous cells, range selection, Ctrl+V from Excel / Google Sheets, fill handle, insert / delete rows and columns, your own extra columns; one line = one genome, so one line uploads one genome). Open an existing .xlsx / .tsv / .csv in the sheet to fix it first; every line has its own **sequence cell** (double-click it to paste the sequence or pick a FASTA file for that line, or drop a FASTA file on it; the cell shows `✔ 189,267 bp`), and many sequences can be attached at once from FASTA files or pasted FASTA by id; the draft is kept in the browser; the dry run highlights the lines with problems. An upload can set a default clade for its blank cells (e.g. *Ib* for a clade Ib batch; a supplied clade is never overwritten, a difference from Nextclade is flagged).
- **Every upload is a version** (what changed: which columns, whether the sequence; source file and a note). The same content again changes nothing; old versions and their sequences are kept, also after a delete, and any version can be restored as a new version. The export manifest records the exported versions. Private uploads respect the coverage floor too.
- `export --with-private` (or the Export tab's *include private*) writes the private genomes as **separate files** (`metadata-private.tsv`, `sequences-private.fasta`) wired to the Nextstrain mpox workflow's own `private_metadata` / `private_sequences`; the public files stay private-free, the manifest says `contains_private`, and the folder gets a `.gitignore` and a notice. Every metadata file now carries an `accession` column equal to the id: the workflow's private-data merge script (`combine_data_sources.py`) still matches rows on that name and fails with "metadata file missing 'accession'" on the current PPX-based metadata without it (verified; with the column the unmodified script merges public + private rows).
- `backup --with-private` stores the private store too (a general backup leaves it out on purpose so a backup can be shared).
- **Export tab, CLI and API for other programs**: `metadata.tsv` + `sequences.fasta` with the same ids in the same order, plus `manifest.json` (counts, warnings, SHA-256 of both files, Nextclade dataset versions, data-use note) and a ready `nextstrain-config.yaml`. The **nextstrain** profile matches the columns the Nextstrain mpox workflow reads (`PPX_accession` as the id, because the workflow maps its include/exclude lists to PPX accessions - genomes without one are left out and reported; `2023-XX-XX` dates, region/country/division/location, host, clade/lineage/outbreak, and the Nextclade `QC_*` columns its default filter uses). `--nextstrain-dir <checkout>/phylogenetic` writes the files and a config overlay straight into a nextstrain/mpox checkout with relative paths (works in its Docker runtime too); the overlay sets `use_remote_files: false` and only replaces `inputs`, so it is passed after the build's own `defaults/<build>/config.yaml`. Checked against that repository: its `map_accessions.py`, `augur filter` (default query) and the Snakemake rules up to `filter` ran on an export, and a dry run of the whole DAG resolves with no remote inputs; **plain** gives every curated column. Filters: minimum length, completeness, exclude QC bad/error, clade, sequence present. `python -m genomecurator export <virus>`, `GET /api/export/{virus}/metadata.tsv|sequences.fasta|preview`, `POST /api/export/{virus}/build`, `GET /api/exports`. Open data only.
- `python -m genomecurator qc <virus> [--force]`; Nextclade results now also keep the QC columns under the names Nextstrain workflows use (re-run QC once to fill them for results made before).
- **QC tab (sequence QC)**: sequence statistics (length, N and ambiguity code share, GC, agreement with the recorded length) for every stored genome, and **Nextclade** (clade, lineage, its own quality score) for viruses that have a Nextclade dataset (mpox, Ebola species, measles, dengue, West Nile, rubella, VZV, hMPV; add more with `GC_NEXTCLADE_DATASETS`). Warnings only. Results are derived data in the new table `qc_result`.
- Nextclade's clade/lineage take part in curation as a fourth *source*: they fill a blank, or - when they differ from NCBI / BV-BRC / Pathoplexus - appear in the review screen with a suggestion (high when a database agrees, medium otherwise) instead of overwriting anything. *Run QC + check* does both steps in one job.
- Nextclade datasets are downloaded on demand and each virus shows its dataset version and whether a newer one exists, with an **Update dataset** button. Genome table gets QC / Nextclade clade / Nextclade lineage columns and a QC filter.
- `environment.yml` now includes the `bioconda` channel and `nextclade` (MIT; x86-64 and ARM). Without Nextclade only the sequence statistics run. Datasets are not redistributed (they ship without a licence file).
- **Source traffic light** in the page header: NCBI — Pathoplexus — BV-BRC connected by a line, each with a light from the newest schema check (blue = no change, yellow = needs a look, red = source code update needed, grey = never checked). Hover for the reason; click to open the Schema tab.
- **Batch suggestions for selected rows**: in the review screen, select rows (checkboxes, **drag the mouse over rows**, **Shift+click** for a range, or filter by integrity and *Select all*) and press *Apply suggestions to selection*. A dialog shows the confidence breakdown and lets you take *medium and above* (recommended) or *everything, including low*. The whole-changeset button is now labelled *Auto-pick (all)*.
- **Analysis tab**: start and stop Apache Superset from the app (it stays a separate program in its own conda env), see whether it is installed/running, open it. `serve --superset` starts it together with the app; a Superset started by the app goes down with it.
- `Dockerfile` (multi-arch: linux/amd64 and linux/arm64, built from `environment.yml`), a `docker` workflow that builds both architectures and publishes to GHCR on version tags, and an ARM job in CI.
- **Jobs tab**: every long task (adding or looking up a virus, check, apply, FASTA download, sequence verification, schema check, publish, and the scheduler's own runs) is a server job with a durable record: running / done / error / interrupted, progress, a timestamped log, the error with where it happened, and a link to that job's lines in the server log. The tab shows a live count of running jobs or unseen failures, and a toast appears when a job ends while you are on another page.
- Progress window gets a **Run in background** button.
- Server-log viewer can be searched (e.g. `job=ab12cd34`).

### Changed
- Genome table: a banner tells how many genomes still lack a FASTA (and roughly how much that is) with an **Get all** button; the confirmation shows size and estimated time, the job shows the remaining time, and the finished message points to the QC tab.
- **Coverage floor**: genomes shorter than `GC_MIN_COVERAGE_PCT` (default **80**) % of the reference (the longest RefSeq) are not worth curating - `check` never proposes them, and ones already held are listed in the changeset as `removed` with the reason and **pre-approved** (soft: status `Removed`, history kept, the FASTA file stays on disk). NCBI *Complete* genomes are exempt; `0` turns the rule off. Removed genomes are hidden from the genome table, dashboard, TSV downloads and QC (the public page already hid them).
- Public page footer now links the sources' policies (NCBI notice, Pathoplexus data-use terms, BV-BRC) and states that only open Pathoplexus data is included. New docs: `docs/LICENSES.md` (licence and data-use review) and `docs/WORKFLOW.md`.
- README rewritten to be short and practical (features, screenshots, installation, running, output), in English with a Korean twin `README.ko.md`; the long explanations moved to `docs/GUIDE.ko.md`.
- Adding a virus and the virus lookup now run as server jobs, so clicking another menu no longer loses them (the lookup resumes where you left it).
- Scheduled checks and the daily schema check are recorded as jobs instead of only appearing in the log.

### Fixed
- Concurrent database connections could fail with `view v_genome_analysis already exists`; the analysis view is now rebuilt only when its definition changed, in one transaction.

### Data & compatibility
- New table `analysis_run` and column `dataset.params_yaml` (created / added automatically; not part of general backups); new dependency `pyyaml` (MIT). The Nextstrain tools are a separate environment and are not part of the Docker image.
- New tables `private_genome` and `private_version` (upload history); new dependency `openpyxl` (MIT) for the spreadsheet upload. Previously: `private_genome` (created automatically; not part of general backups, use `backup --with-private`). Sequences of private genomes live in `data/private/<virus>/`. Schema version stays 1.
- New derived table `qc_result` (created automatically; not part of backups, recomputed by the QC tab). Schema version stays 1.
- New operational table `job` (created automatically; not included in backups). Schema version stays 1.

## [0.2.0] - 2026-10-17

Makes the quality of a genome and the disagreement between sources something you can *act on*: a Near complete class, a way to choose between the
sources' values with a suggestion, and two ways to read the schema check. The tool is now explicitly for viruses only.

### Changed
- A host that differs only by a common name in parentheses (`Sus scrofa` vs `Sus scrofa (pig)`) is now a *cosmetic* difference, not a conflict.
- **Scope and name.** The tool is for **viruses only** (NCBI Virus, BV-BRC, Pathoplexus); bacteria will be a separate project because genomes there are
  assemblies made of many contigs and are megabases long. To stop "genome" from being read as *human genome*, the displayed name is now
  **Virus Genome Metadata Curator** (**Virus Genome Curator** in the page header, titles and messages). The repository, the Python package `genomecurator`, the
  command line, the data folder, environment variables and the backup format identifier are **unchanged**, so nothing needs migrating.

### Added
- **Choosing between the sources' values** in the review screen. Click an item and the side panel lists every field the sources disagree on, with each source's value and whether
  it is only written differently or really different. A **suggestion (★)** is pre-selected with its reason and confidence: an impossible date is left out; length/version follow NCBI
  (the sequence is downloaded from NCBI and its length verified); values that are the same fact in different words go to the best-written one (most specific date, scientific
  host name, standard country name) = *high*; a value most sources agree on = *medium*; otherwise the source priority decides and the person is told to look = *low*.
  Pick a value or type your own and save: it is applied to the item at once and **kept by every later check** (the chosen source's current value is used; a typed value stays as typed;
  *manual* is shown as the source). Choices can be undone per field. **Auto-pick** in the toolbar selects the suggestions of *medium* confidence or better in one go and leaves the *low* ones
  to the person. Choices are stored in the new table `field_choice` and are part of backups; fields a person has decided no longer count as open disagreements.
- **Near complete** genomes. NCBI only distinguishes Complete from Partial, but most of the "Partial" genomes are missing just a few hundred bases at the ends. A new
  *Completeness* class sits between them: **Complete** (NCBI says Complete; trusted, so a lineage that is legitimately shorter than the reference is not demoted),
  **Near complete** (not Complete, length at least **90%** of the reference genome) and **Partial** (below that, or no reference). The reference of a virus is its longest
  RefSeq genome among the collected records. The threshold is `GC_NEAR_COMPLETE_PCT` (default 90, the CheckV high-quality cut-off; on real data it sits in the valley of the length
  distribution: for mpox only 33 of 13,462 genomes lie between 50% and 90%). Shown as a column, a filter and dashboard KPIs, on the public page, and in the analysis view
  (`completeness`, `ref_coverage_pct`, `ref_length`), which applies the same rule in SQL. On the mpox data 7,583 of the 8,223 NCBI-Partial genomes are Near complete.
- **Simple view** on the **Schema** tab (the default; a tab switches to the detailed board below): one table with the local `genome` columns as rows and
  NCBI / BV-BRC / Pathoplexus as columns, each cell naming the source field that fills the column, coloured by status. Columns that need attention are
  listed on top with the source field to look at; columns are grouped (identity, taxonomy, genome, sample, dates, external IDs). The choice of view is remembered.
- **Table structure board** (the detailed view) on the **Schema** tab: the fields of NCBI (nested objects kept as groups), BV-BRC and Pathoplexus on the left, the columns of
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

### Data & compatibility
- Schema version stays **1**. The only database change is the new table `field_choice` (your choices between source values), created automatically the first time
  0.2.0 opens a database; nothing existing is altered. A 0.1.0 backup restores into 0.2.0, and 0.2.0 backups include `field_choice`.
- The analysis view `v_genome_analysis` gained the columns `completeness`, `ref_coverage_pct` and `ref_length`. In Superset run `scripts/superset.sh bootstrap` once
  to refresh the dataset and see them. The view reads the Near complete threshold from `app_meta.near_complete_pct`, which the app keeps in step with `GC_NEAR_COMPLETE_PCT`.
- New environment variable: `GC_NEAR_COMPLETE_PCT` (default 90).
- Genomes that NCBI marks *Partial* may now read *Near complete*; nothing is rewritten, the class is computed when shown. The NCBI status is still available as its own column.

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

[Unreleased]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/CHOICHI-HWAN/GenomeCurator/releases/tag/v0.1.0
