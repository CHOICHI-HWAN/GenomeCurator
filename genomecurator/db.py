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
-- every background job (jobs.py): status, progress, log, error, result - survives page reloads and server restarts
CREATE TABLE IF NOT EXISTS job(
  id TEXT PRIMARY KEY, kind TEXT, label TEXT, status TEXT,  -- running|done|error|interrupted
  pct REAL DEFAULT 0, msg TEXT, started TEXT, finished TEXT, error TEXT, redirect TEXT, log_json TEXT, result_json TEXT);
CREATE INDEX IF NOT EXISTS ix_job_started ON job(started);
-- a person's choice between the sources' values (resolve.py); honoured by every later check
CREATE TABLE IF NOT EXISTS field_choice(
  watch_id INTEGER, accession TEXT, field TEXT, source TEXT,  -- ncbi|bvbrc|pathoplexus|custom
  value TEXT, auto INTEGER DEFAULT 0, decided_at TEXT, PRIMARY KEY(watch_id, accession, field));
-- which version of the code wrote / last opened this database (see record_version)
CREATE TABLE IF NOT EXISTS app_meta(key TEXT PRIMARY KEY, value TEXT, updated TEXT DEFAULT CURRENT_TIMESTAMP);
-- schema dry-run (schema_check.py): last live field profile per source, and the report history
CREATE TABLE IF NOT EXISTS schema_snapshot(source TEXT PRIMARY KEY, taken TEXT, snapshot_json TEXT);
CREATE TABLE IF NOT EXISTS schema_report(id INTEGER PRIMARY KEY, checked_at TEXT, report_json TEXT);
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
  r.ref_length,
  CASE WHEN r.ref_length IS NOT NULL AND g.genome_length > 0 THEN ROUND(g.genome_length * 100.0 / r.ref_length, 1) END AS ref_coverage_pct,
  CASE WHEN g.genome_status = 'Removed' THEN 'Removed' WHEN g.genome_status = 'Complete' THEN 'Complete'
       WHEN r.ref_length IS NOT NULL AND g.genome_length > 0 AND g.genome_length * 100.0 / r.ref_length
            >= COALESCE((SELECT CAST(value AS REAL) FROM app_meta WHERE key = 'near_complete_pct'), 90) THEN 'Near complete'
       ELSE 'Partial' END AS completeness,
  1 AS n
FROM genome g
LEFT JOIN watchlist w ON w.id = g.watch_id
LEFT JOIN (SELECT watch_id, MAX(genome_length) AS ref_length FROM genome
           WHERE refseq_accessions IS NOT NULL AND refseq_accessions != '' AND genome_length > 0 AND IFNULL(genome_status,'') != 'Removed'
           GROUP BY watch_id) r ON r.watch_id = g.watch_id
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
    # flat view for BI tools (Superset); rebuilt only when its definition changed. Several threads connect at once (background jobs,
    # pollers), so the drop+create is one transaction and a lost race is harmless.
    cur = con.execute("SELECT sql FROM sqlite_master WHERE type='view' AND name='v_genome_analysis'").fetchone()
    if cur is None or (cur[0] or "").strip() != ANALYSIS_VIEW.strip():
        try:
            con.execute("BEGIN IMMEDIATE")
            con.execute("DROP VIEW IF EXISTS v_genome_analysis")
            con.execute(ANALYSIS_VIEW)
            con.commit()
        except sqlite3.OperationalError:
            con.rollback()
            if con.execute("SELECT 1 FROM sqlite_master WHERE name='v_genome_analysis'").fetchone() is None:
                raise


_recorded: set[str] = set()


def record_version(con: sqlite3.Connection) -> None:
    """Remember which GenomeCurator version works on this database: `app_meta` holds first/current version and the DB
    schema version, and every change of version is written to the audit log (so "which code produced this data" is answerable).
    Runs once per process and database file. A database written by a NEWER layout than this code knows is only warned about."""
    key = str(config.DB_PATH)
    if key in _recorded:
        return
    from . import DB_SCHEMA_VERSION, __version__, log
    meta = {r[0]: r[1] for r in con.execute("SELECT key, value FROM app_meta")}
    prev, prev_schema = meta.get("app_version"), int(meta.get("db_schema_version") or 0)
    if prev_schema > DB_SCHEMA_VERSION:
        log.get("db").error("this database uses schema v%s, newer than this code (v%s, GenomeCurator %s); update GenomeCurator before changing data",
                            prev_schema, DB_SCHEMA_VERSION, __version__)
    if prev != __version__:
        con.execute("INSERT OR IGNORE INTO app_meta(key, value) VALUES('first_version', ?)", (prev or __version__,))
        for k, v in (("app_version", __version__), ("db_schema_version", str(DB_SCHEMA_VERSION))):
            con.execute("INSERT INTO app_meta(key, value, updated) VALUES(?,?,CURRENT_TIMESTAMP) ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated=CURRENT_TIMESTAMP", (k, v))
        audit(con, "version", f"Virus Genome Curator {__version__} (DB schema {DB_SCHEMA_VERSION})" + (f" - 이전 {prev}" if prev else " - 이 데이터베이스를 처음 열었습니다"))
        log.get("db").info("database now handled by Virus Genome Curator %s (was %s)", __version__, prev or "new")
    pct = str(config.NEAR_COMPLETE_PCT)      # the analysis view reads the threshold from here (same rule as completeness.classify)
    con.execute("INSERT INTO app_meta(key, value, updated) VALUES('near_complete_pct', ?, CURRENT_TIMESTAMP) "
                "ON CONFLICT(key) DO UPDATE SET value=excluded.value, updated=CURRENT_TIMESTAMP WHERE value != excluded.value", (pct,))
    _recorded.add(key)


def version_info() -> dict:
    from . import DB_SCHEMA_VERSION, __version__
    con = connect()
    try:
        meta = {r[0]: r[1] for r in con.execute("SELECT key, value FROM app_meta")}
    finally:
        con.close()
    return {"version": __version__, "db_schema_version": DB_SCHEMA_VERSION, "first_version": meta.get("first_version"),
            "db_recorded_version": meta.get("app_version")}


def parse_version(v: str) -> tuple:
    """'0.10.2' -> (0, 10, 2); anything non-numeric sorts as 0 (good enough for "is this backup newer than me")."""
    import re
    return tuple(int(x) for x in (re.findall(r"\d+", v or "")[:3] + ["0", "0", "0"])[:3])


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
    record_version(con)
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


AUDIT_MAX = 600


def _short(detail: str) -> str:
    return detail if len(detail) <= AUDIT_MAX else detail[:AUDIT_MAX] + f"… (+{len(detail) - AUDIT_MAX}자)"


def audit(con, action: str, detail: str = ""):
    """Record a user-visible action inside the caller's transaction (it is rolled back with it)."""
    con.execute("INSERT INTO audit_log(action, detail) VALUES(?,?)", (action, _short(detail)))


def audit_now(action: str, detail: str = ""):
    """Record immediately in its own transaction - for failures, which must survive the rollback of the failed work."""
    try:
        con = connect()
        try:
            con.execute("INSERT INTO audit_log(action, detail) VALUES(?,?)", (action, _short(detail)))
            con.commit()
        finally:
            con.close()
    except Exception:   # never let bookkeeping hide the original error
        from . import log
        log.get("db").exception("could not write audit_log (%s)", action)

