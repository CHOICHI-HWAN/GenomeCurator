import sqlite3
from contextlib import contextmanager

from . import config

SCHEMA = """
CREATE TABLE IF NOT EXISTS watchlist(
  id INTEGER PRIMARY KEY, name TEXT UNIQUE NOT NULL, taxon TEXT NOT NULL,
  taxon_id INTEGER, filters_json TEXT DEFAULT '{}', enabled INTEGER DEFAULT 1,
  last_checked TEXT, created TEXT DEFAULT CURRENT_TIMESTAMP);

-- BV-BRC style genome table (one row per INSDC accession; non-segmented viruses)
CREATE TABLE IF NOT EXISTS genome(
  genome_id TEXT PRIMARY KEY, watch_id INTEGER, accession TEXT UNIQUE NOT NULL,
  accession_version INTEGER, genbank_accessions TEXT, refseq_accessions TEXT,
  taxon_id INTEGER, genome_name TEXT, taxon_lineage_ids TEXT, taxon_lineage_names TEXT,
  genome_status TEXT, genome_length INTEGER, segment TEXT, strain TEXT, subtype TEXT,
  genotype TEXT, host_name TEXT, isolation_country TEXT, geographic_location TEXT,
  collection_date TEXT, isolation_source TEXT, publication TEXT, submitter TEXT,
  genome_quality TEXT, sequence_md5 TEXT, sequence_hash TEXT, bvbrc_genome_id TEXT, pathoplexus_accession TEXT, release_date TEXT, update_date TEXT,
  date_inserted TEXT, date_modified TEXT);

CREATE TABLE IF NOT EXISTS genome_sequence(
  accession TEXT PRIMARY KEY, genome_id TEXT NOT NULL, version INTEGER, length INTEGER,
  gc_content REAL, md5 TEXT, file_path TEXT);

CREATE TABLE IF NOT EXISTS source_record(
  source TEXT NOT NULL, accession TEXT NOT NULL, watch_id INTEGER,
  payload_json TEXT, fetched_at TEXT DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY(source, accession));

CREATE TABLE IF NOT EXISTS field_provenance(
  genome_id TEXT, field TEXT, source TEXT, PRIMARY KEY(genome_id, field));

CREATE TABLE IF NOT EXISTS genome_version(
  id INTEGER PRIMARY KEY, genome_id TEXT, snapshot_json TEXT, changeset_id INTEGER,
  changed_at TEXT DEFAULT CURRENT_TIMESTAMP);

CREATE TABLE IF NOT EXISTS changeset(
  id INTEGER PRIMARY KEY, watch_id INTEGER, created TEXT DEFAULT CURRENT_TIMESTAMP,
  status TEXT DEFAULT 'PENDING',  -- PENDING|READY|BLOCKED|APPLIED|DISCARDED
  summary_json TEXT, dryrun_json TEXT, applied_at TEXT);

CREATE TABLE IF NOT EXISTS changeset_item(
  id INTEGER PRIMARY KEY, changeset_id INTEGER, accession TEXT, kind TEXT,
  diff_json TEXT, new_record_json TEXT, provenance_json TEXT, fingerprint TEXT,
  decision TEXT DEFAULT 'pending');  -- pending|approved|rejected|hold

-- cross-source disagreements found at the last check (regenerated each check); see reconcile.py
CREATE TABLE IF NOT EXISTS source_conflict(
  id INTEGER PRIMARY KEY, watch_id INTEGER, accession TEXT, field TEXT, kind TEXT,  -- cosmetic|conflict
  src_a TEXT, val_a TEXT, src_b TEXT, val_b TEXT);
CREATE INDEX IF NOT EXISTS ix_conflict_watch ON source_conflict(watch_id, field, kind);
-- verdict of the sequence-body comparison against NCBI (seqcheck.py); no sequences are stored
CREATE TABLE IF NOT EXISTS sequence_check(
  watch_id INTEGER, accession TEXT, source TEXT, status TEXT,  -- ref|same|diff|missing|error
  length INTEGER, sha256 TEXT, ref_length INTEGER, prefix_same INTEGER, suffix_same INTEGER, detail TEXT, checked_at TEXT,
  PRIMARY KEY(watch_id, accession, source));
CREATE TABLE IF NOT EXISTS sequence_plan(watch_id INTEGER PRIMARY KEY, created TEXT, plan_json TEXT);  -- dry-run of seqcheck
CREATE TABLE IF NOT EXISTS reconcile_report(watch_id INTEGER PRIMARY KEY, checked_at TEXT, report_json TEXT);

CREATE TABLE IF NOT EXISTS audit_log(
  id INTEGER PRIMARY KEY, ts TEXT DEFAULT CURRENT_TIMESTAMP, action TEXT, detail TEXT);
"""

ANALYSIS_VIEW = """
CREATE VIEW v_genome_analysis AS
SELECT g.genome_id, g.accession, w.name AS virus, g.genome_name, g.strain, g.genome_status, g.genome_quality,
  g.genome_length,
  CASE WHEN g.genome_length IS NULL THEN NULL WHEN g.genome_length < 1000 THEN '<1 kb'
       WHEN g.genome_length < 5000 THEN '1-5 kb' WHEN g.genome_length < 10000 THEN '5-10 kb'
       WHEN g.genome_length < 50000 THEN '10-50 kb' ELSE '>=50 kb' END AS length_bin,
  COALESCE(NULLIF(g.host_name,''), '(unknown)') AS host,
  COALESCE(NULLIF(g.isolation_country,''), '(unknown)') AS country,
  g.geographic_location, g.isolation_source,
  g.collection_date,
  CASE WHEN g.collection_date GLOB '[12][0-9][0-9][0-9]*' THEN CAST(substr(g.collection_date,1,4) AS INTEGER) END AS collection_year,
  CASE WHEN g.collection_date GLOB '[12][0-9][0-9][0-9]-[01][0-9]*' THEN substr(g.collection_date,1,7) END AS collection_month,
  CASE WHEN g.release_date GLOB '[12][0-9][0-9][0-9]*' THEN CAST(substr(g.release_date,1,4) AS INTEGER) END AS release_year,
  g.release_date, g.update_date,
  g.clade, g.lineage, g.genotype, g.subtype, g.serotype, g.outbreak,
  g.sequencing_platform, g.assembly_method, g.submitter,
  CASE WHEN g.sra_accession IS NOT NULL AND g.sra_accession != '' THEN 'Yes' ELSE 'No' END AS has_sra,
  CASE WHEN g.pathoplexus_accession IS NOT NULL AND g.pathoplexus_accession != '' THEN 'Yes' ELSE 'No' END AS in_pathoplexus,
  CASE WHEN s.file_path IS NOT NULL THEN 'Yes' ELSE 'No' END AS has_fasta,
  g.bioproject_accession, g.biosample_accession,
  1 AS n
FROM genome g
LEFT JOIN watchlist w ON w.id = g.watch_id
LEFT JOIN genome_sequence s ON s.genome_id = g.genome_id
"""


def migrate(con: sqlite3.Connection) -> None:
    """Add columns introduced after the base SCHEMA (also used to build backup copies)."""
    have = {r[1] for r in con.execute("PRAGMA table_info(genome)")}
    for col in ("sra_accession", "bioproject_accession", "biosample_accession", "assembly_method", "sequencing_platform",
                "clade", "lineage", "serotype", "outbreak"):
        if col not in have:
            con.execute(f"ALTER TABLE genome ADD COLUMN {col} TEXT")
    wcols = {r[1] for r in con.execute("PRAGMA table_info(watchlist)")}
    if "pathoplexus_organism" not in wcols:
        con.execute("ALTER TABLE watchlist ADD COLUMN pathoplexus_organism TEXT")
    con.execute("DROP VIEW IF EXISTS v_genome_analysis")   # flat view for BI tools (Superset); rebuilt on every start
    con.execute(ANALYSIS_VIEW)


def connect() -> sqlite3.Connection:
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(config.DB_PATH, timeout=30)
    con.row_factory = sqlite3.Row
    try:
        con.execute("PRAGMA journal_mode=WAL")
    except sqlite3.OperationalError:  # e.g. some network / virtual file systems do not support WAL
        pass
    con.execute("PRAGMA foreign_keys=ON")
    con.executescript(SCHEMA)
    migrate(con)
    con.commit()
    return con


@contextmanager
def tx():
    con = connect()
    try:
        yield con
        con.commit()
    except BaseException:
        con.rollback()
        raise
    finally:
        con.close()


def audit(con, action: str, detail: str = ""):
    con.execute("INSERT INTO audit_log(action, detail) VALUES(?,?)", (action, detail))

