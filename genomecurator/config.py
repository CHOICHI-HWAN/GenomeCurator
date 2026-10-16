"""Runtime settings. Override with environment variables."""
import os
import re
import sys
from pathlib import Path


def _default_data_dir() -> Path:
    """GC_DATA_DIR, else <repo>/data when running from a clone (the usual case), else the OS user-data folder.
    Never relative to the current directory, so the same data is used wherever the app is started from."""
    env = os.environ.get("GC_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    repo = Path(__file__).resolve().parent.parent
    if (repo / "pyproject.toml").exists():
        return repo / "data"
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local")
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_DATA_HOME") or Path.home() / ".local" / "share")
    return base / "GenomeCurator"


def safe_name(name: str) -> str:
    """Folder/file-safe version of a watch name (spaces, slashes, colons, reserved characters on Windows)."""
    return re.sub(r"[^\w.\-]+", "_", name, flags=re.UNICODE).strip("._") or "virus"


DATA_DIR = _default_data_dir()
DB_PATH = DATA_DIR / "genomecurator.db"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY") or None
# Minimum gap between any two outgoing web requests (NCBI, BV-BRC, Pathoplexus).
REQUEST_DELAY_S = float(os.environ.get("GC_REQUEST_DELAY_MS", "300")) / 1000
CHECK_INTERVAL_HOURS = float(os.environ.get("GC_CHECK_INTERVAL_HOURS", "24"))
DRYRUN_SAMPLE_SIZE = int(os.environ.get("GC_DRYRUN_SAMPLE", "3"))
MIN_FREE_BYTES = int(os.environ.get("GC_MIN_FREE_BYTES", str(500 * 1024 * 1024)))

# Source priority: earlier wins when the same field conflicts.
SOURCE_PRIORITY = ["pathoplexus", "ncbi", "bvbrc"]

# watchlist name -> Pathoplexus (Loculus/LAPIS) organism slug. Absent = skipped.
PATHOPLEXUS_ORGANISMS: dict[str, str] = {}
PATHOPLEXUS_BASE = "https://lapis.pathoplexus.org"
NCBI_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus"
BVBRC_BASE = "https://www.bv-brc.org/api/genome/"

# Metadata fields whose change counts as a "changed" record.
TRACKED_FIELDS = [
    "genome_name", "genome_status", "genome_length", "strain", "host_name",
    "isolation_country", "geographic_location", "collection_date",
    "isolation_source", "sequence_hash", "publication", "submitter",
    "genome_quality", "subtype", "genotype", "segment",
    "sra_accession", "bioproject_accession", "biosample_accession", "assembly_method", "sequencing_platform",
    "clade", "lineage", "serotype", "outbreak",
]

