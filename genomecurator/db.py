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

CREATE TABLE IF NOT EXISTS audit_log(
  id INTEGER PRIMARY KEY, ts TEXT DEFAULT CURRENT_TIMESTAMP, action TEXT, detail TEXT);
"""


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
    have = {r[1] for r in con.execute("PRAGMA table_info(genome)")}
    for col in ("sra_accession", "bioproject_accession", "biosample_accession", "assembly_method", "sequencing_platform",
                "clade", "lineage", "serotype", "outbreak"):
        if col not in have:
            con.execute(f"ALTER TABLE genome ADD COLUMN {col} TEXT")
    wcols = {r[1] for r in con.execute("PRAGMA table_info(watchlist)")}
    if "pathoplexus_organism" not in wcols:
        con.execute("ALTER TABLE watchlist ADD COLUMN pathoplexus_organism TEXT")
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

