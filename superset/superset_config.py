"""Superset settings for GenomeCurator (local analysis). Loaded via SUPERSET_CONFIG_PATH."""
import os
from pathlib import Path

HOME = Path(os.environ.get("SUPERSET_HOME", Path(__file__).parent / "home"))
HOME.mkdir(parents=True, exist_ok=True)

_key = HOME / "secret_key"
if not _key.exists():
    import secrets
    _key.write_text(secrets.token_urlsafe(48), encoding="utf-8")
SECRET_KEY = _key.read_text(encoding="utf-8").strip()

SQLALCHEMY_DATABASE_URI = f"sqlite:///{(HOME / 'superset.db').as_posix()}"
# The analysed GenomeCurator database is SQLite; Superset blocks it unless this is off. Local use only.
PREVENT_UNSAFE_DB_CONNECTIONS = False
WTF_CSRF_ENABLED = True
TALISMAN_ENABLED = False
ROW_LIMIT = 100000
