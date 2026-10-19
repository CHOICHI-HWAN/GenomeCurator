# Virus Genome Metadata Curator

[![tests](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml/badge.svg)](https://github.com/CHOICHI-HWAN/GenomeCurator/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.10%20%7C%203.12-blue)
![version](https://img.shields.io/github/v/tag/CHOICHI-HWAN/GenomeCurator?label=version&sort=semver)
![license](https://img.shields.io/badge/license-MIT-green)

[한국어](README.ko.md)

A local web tool that collects **virus** genome metadata from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/), [BV-BRC](https://www.bv-brc.org/) and [Pathoplexus](https://pathoplexus.org/), compares what the three sources say about the same genome, and adds records to your database only after you have approved them. Package and command name: `genomecurator`.

It is for viruses only. Bacteria differ in record structure (multi-contig assemblies) and sequence size, and will be a separate project.

## Features

1. **Nothing changes without approval.** New and changed records collect as a *changeset*; only what you approve and apply reaches your data. Rejected items do not come back.
2. **Three sources, one record.** NCBI Virus defines which accessions exist; Pathoplexus > NCBI > BV-BRC decide each field, and every value keeps a record of the source it came from.
3. **Agreement between sources is measured.** Per field: match / notation difference (`2022-07` vs `2022-07-11`) / real conflict.
4. **Sequence verification.** Disagreements can be settled by comparing the sequence body (SHA-256, headers removed), with a dry-run first so large downloads are avoided.
5. **Choose between sources' values**, with a suggested value, the reason and a confidence level; choices survive later checks.
6. **Schema watch.** Detects when a source renames or removes a field your adapters read, daily and on demand.
7. **Genome completeness classes**: Complete / Near complete (≥ 90% of the longest RefSeq genome, adjustable) / Partial.
8. **Jobs tab.** Every long task runs on the server with status, progress, log and error; it survives page changes and restarts.
9. **Sequence QC.** Sequence statistics (length, N and ambiguity, GC) for every genome, and [Nextclade](https://github.com/nextstrain/nextclade) clade/lineage plus its quality score for viruses that have a Nextclade dataset (mpox, Ebola species, measles, dengue, ...). A clade that differs from the sources goes to the review screen instead of being overwritten; datasets are updated from a button.
10. **Analysis with Apache Superset** (optional, read-only view of your data) and a **read-only static page** for GitHub Pages.

## Screenshots

![From an empty screen to an approved changeset](docs/demo/01-quickstart.gif)

*Register Henipavirus nipahense, check the three sources, approve 442 records, apply (real data, sped up).*

| | |
|---|---|
| ![Integrity tab](docs/demo/02-integrity.gif) | ![Resolve screen](docs/demo/05-resolve.gif) |
| Integrity: agreement per field, real conflicts, sequence-check plan | Review: per-source values, suggestion, Auto-pick |
| ![Schema tab](docs/demo/03-schema.gif) | ![Genome table](docs/demo/04-explore.gif) |
| Schema: simple view and table-structure board | Genome table with provenance in the detail panel |

## Installation

Requires Python 3.10+ and internet access to the three sources. x86-64 and ARM (Apple Silicon, Linux aarch64) are supported.

**1. conda / mamba (recommended)**

```bash
git clone https://github.com/CHOICHI-HWAN/GenomeCurator.git
cd GenomeCurator
mamba env create -f environment.yml      # or: conda env create -f environment.yml (includes Nextclade from bioconda)
mamba activate genomecurator
python -m pytest -q                      # optional: all green means it works
```

*One-step start (creates the environment if needed, starts the server, opens the browser):* `./scripts/start.sh` (Linux), `scripts/start.command` (macOS), `scripts\start.bat` (Windows).

**2. Docker** (linux/amd64 and linux/arm64)

```bash
docker build -t genomecurator .
docker run -d -p 127.0.0.1:8001:8001 -v genomecurator-data:/data genomecurator
```

*The app has no login. Publish the port to 127.0.0.1 only.*

**3. pip**

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Running

```bash
python -m genomecurator watch-add "Henipavirus nipahense"   # register a virus (the name is checked against NCBI)
python -m genomecurator serve --open                        # web UI at http://127.0.0.1:8000
```

Then in the browser: **Check now → Review → select items → Approve → Apply**. Each page has a short blue help box.

Other commands (`python -m genomecurator --help` for all options):

| Command | What it does |
|---|---|
| `check <name>` | Query the three sources and create a changeset (data unchanged) |
| `dry-run <id>` / `apply <id>` | Pre-apply checks / apply approved items |
| `reconcile <name>` | Source-vs-source agreement report |
| `verify` | Integrity checks on local data |
| QC tab | Sequence statistics and Nextclade per virus; update Nextclade datasets; *Run QC + check* |
| `schema-check [--daily]` | Check the sources' and local schema |
| `backup [zip]` / `restore <zip>` | Portable backup of curated data and FASTA |
| `serve [--port N] [--open] [--superset] [--no-scheduler]` | Web server (optionally with Superset) |

Useful settings (environment variables): `GC_DATA_DIR` (data location), `NCBI_API_KEY` (higher NCBI rate limit), `GC_CHECK_INTERVAL_HOURS` (default 24), `GC_NEAR_COMPLETE_PCT` (default 90), `GC_SCHEMA_CHECK_HOUR` (default 3). The full list is in the [guide](docs/GUIDE.ko.md#설정-환경변수).

## Output

```
data/
├── genomecurator.db          # SQLite: curated tables, provenance, history, job records
├── genomes/<virus>/*.fna     # downloaded sequences
├── snapshots/  staging/      # safety copies and temp files while applying
└── logs/genomecurator.log    # server log (rotated, with job/watch context)
```

`data/` is not committed to git. Analysis views such as `v_genome_analysis` are in the database for BI tools.

## Limitations

- The server must be running for the scheduled checks.
- The daily schema check samples a few records per source; use the full check in the Schema tab now and then.
- Identical sequences prove "same record", not that host or date are correct.

## Documentation

The detailed guide is in Korean: [docs/GUIDE.ko.md](docs/GUIDE.ko.md) (tabs, how integrity checks work, FAQ, reading logs, code map).
Also: [CHANGELOG.md](CHANGELOG.md) · [docs/SUPERSET.md](docs/SUPERSET.md) · [docs/DEPENDENCIES.md](docs/DEPENDENCIES.md) · [docs/RELEASING.md](docs/RELEASING.md).

## Contributing and license

Bug reports and ideas: [Issues](https://github.com/CHOICHI-HWAN/GenomeCurator/issues). MIT License, © 2026 CHOICHI-HWAN.
Please be kind to the data providers: requests are spaced out on purpose (`GC_REQUEST_DELAY_MS`).
