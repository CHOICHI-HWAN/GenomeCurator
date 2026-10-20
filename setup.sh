#!/usr/bin/env bash
# Creates the conda environments of Virus Genome Metadata Curator with mamba (falls back to conda when mamba is missing).
# Creating uses mamba; ACTIVATING an environment is done with conda: `conda activate genomecurator`.
#
#   ./setup.sh                  all environments
#   ./setup.sh --core           only the app environment (enough for curation, QC statistics, Nextclade)
#   ./setup.sh --only gc-squirrel,gc-snippy     just those (names below)
#   ./setup.sh --update         update environments that already exist instead of skipping them
#   ./setup.sh --dry-run        print what would be done
#
# Environments (name <- file): the app and everything the curation tabs need, and one separate environment per analysis tool (they are run as separate programs):
#   genomecurator <- environment.yml            gc-nextstrain <- nextstrain/environment.yml   gc-squirrel <- squirrel/environment.yml
#   gc-snippy     <- variantcall/snippy.yml     gc-lofreq     <- variantcall/lofreq.yml        gc-snakemake <- variantcall/snakemake.yml   gc-checkv <- variantcall/checkv.yml (its database, about 1.5 GB, is downloaded separately: see the file)
# Check the result with ./test.sh check
set -euo pipefail
cd "$(dirname "$0")"

ENVS=(
  "genomecurator:environment.yml"
  "gc-nextstrain:nextstrain/environment.yml"
  "gc-squirrel:squirrel/environment.yml"
  "gc-snippy:variantcall/snippy.yml"
  "gc-lofreq:variantcall/lofreq.yml"
  "gc-snakemake:variantcall/snakemake.yml"
  "gc-checkv:variantcall/checkv.yml"
)
ONLY=""; UPDATE=0; DRY=0
while [ $# -gt 0 ]; do
  case "$1" in
    --core) ONLY="genomecurator" ;;
    --only) ONLY="${2:?--only needs a comma-separated list of environment names}"; shift ;;
    --update) UPDATE=1 ;;
    --dry-run) DRY=1 ;;
    -h|--help) sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1 (see ./setup.sh --help)" >&2; exit 2 ;;
  esac
  shift
done

TOOL=""
for t in mamba conda; do command -v "$t" >/dev/null 2>&1 && TOOL="$t" && break; done
[ -n "$TOOL" ] || { echo "mamba or conda is required (https://github.com/conda-forge/miniforge)" >&2; exit 1; }
[ "$TOOL" = "mamba" ] || echo "note: mamba not found, using conda (slower); install mamba with: conda install -n base -c conda-forge mamba"
command -v conda >/dev/null 2>&1 || echo "note: conda is not on PATH; activate environments with 'mamba activate' or put conda on PATH"

exists() { "$TOOL" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }
wanted() { [ -z "$ONLY" ] && return 0; case ",$ONLY," in *",$1,"*) return 0 ;; esac; return 1; }
run() { echo "+ $*"; [ "$DRY" = 1 ] || "$@"; }

made=0; skipped=0
for item in "${ENVS[@]}"; do
  name="${item%%:*}"; file="${item#*:}"
  wanted "$name" || continue
  [ -f "$file" ] || { echo "missing $file" >&2; exit 1; }
  if exists "$name"; then
    if [ "$UPDATE" = 1 ]; then echo "== updating $name"; run "$TOOL" env update -n "$name" -f "$file" --prune; made=$((made+1))
    else echo "== $name exists (use --update to update it)"; skipped=$((skipped+1)); fi
  else
    echo "== creating $name from $file"; run "$TOOL" env create -y -f "$file"; made=$((made+1))
  fi
done

echo
echo "done: $made created/updated, $skipped already present."
echo "next:  conda activate genomecurator   then   ./test.sh check   and   python -m genomecurator serve --open"
