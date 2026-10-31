"""Virus Genome Metadata Curator (package `genomecurator`). The version lives here and nowhere else (pyproject reads it); see CHANGELOG.md and docs/RELEASING.md."""
__version__ = "0.3.0"

_build: dict | None = None


def build() -> dict:
    """`{"version", "commit", "dirty"}` — what code this actually is.

    The version alone names a release, and between releases it names nothing: two weeks and 160 commits all called themselves
    0.2.0, so a worker could not be told from a worker on a fortnight-old checkout, and a reproducibility manifest or a
    changeset stamped with the version alone could not be traced back to the code that made it. The commit closes that gap;
    `dirty` says whether anything was uncommitted when it ran, because then the commit does not describe it either.
    Read once per process - it shells out to git, and a checkout without git simply has no commit."""
    global _build
    if _build is None:
        import subprocess
        from pathlib import Path
        root = Path(__file__).resolve().parent.parent

        def git(*a: str) -> str:
            try:
                return subprocess.run(["git", "-C", str(root), *a], capture_output=True, text=True, timeout=10).stdout.strip()
            except Exception:
                return ""
        c = git("rev-parse", "--short=10", "HEAD")
        _build = {"version": __version__, "commit": c or None,
                  "dirty": bool(git("status", "--porcelain", "--untracked-files=no")) if c else None}
    return dict(_build)
# Bumped when the database layout changes in a way older code cannot read. Today every change is an additive migration
# (db.migrate), so it stays 1; it is recorded in the database and in backups so a newer data folder is never opened blindly.
DB_SCHEMA_VERSION = 1
