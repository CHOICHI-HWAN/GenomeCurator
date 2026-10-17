"""Virus Genome Metadata Curator (package `genomecurator`). The version lives here and nowhere else (pyproject reads it); see CHANGELOG.md and docs/RELEASING.md."""
__version__ = "0.2.0"
# Bumped when the database layout changes in a way older code cannot read. Today every change is an additive migration
# (db.migrate), so it stays 1; it is recorded in the database and in backups so a newer data folder is never opened blindly.
DB_SCHEMA_VERSION = 1
