"""Runtime settings. Override with environment variables."""
import os
from pathlib import Path

DATA_DIR = Path(os.environ.get("GC_DATA_DIR", "data")).resolve()
DB_PATH = DATA_DIR / "genomecurator.db"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY") or None
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
]

