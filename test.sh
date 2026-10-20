#!/usr/bin/env bash
# Tests for Virus Genome Metadata Curator.
#   ./test.sh unit       the unit tests (pytest; no internet, no analysis tools needed) - run in the `genomecurator` environment
#   ./test.sh check      does the installation WORK: environments, each tool answers, and a small real end-to-end run of the Snakemake variant flow
#   ./test.sh verify     are the RESULTS of the pipelines right, on the real data in the data folder: invariants, cross-checks between independent methods, a planted-mutation test, privacy (read-only)
#   ./test.sh            unit + check (verify is run on request: it needs finished analyses)
#   ./test.sh -h         this text
# Run it from an activated environment (`conda activate genomecurator`); if none is active it activates `genomecurator` itself with conda.
# Exit code 0 = passed. `check` reports optional parts that are missing as WARN (they do not fail the run), required ones as FAIL.
set -uo pipefail
cd "$(dirname "$0")"
MODE="${1:-all}"
case "$MODE" in -h|--help) sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;; unit|check|verify|all) ;; *) echo "usage: ./test.sh [unit|check]" >&2; exit 2 ;; esac

PY=python
if ! python -c "import genomecurator, pytest" >/dev/null 2>&1; then
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate genomecurator 2>/dev/null || true
  fi
  python -c "import genomecurator, pytest" >/dev/null 2>&1 || { echo "genomecurator (and pytest) not importable: run ./setup.sh --core, then 'conda activate genomecurator'" >&2; exit 1; }
fi

rc=0
if [ "$MODE" = unit ] || [ "$MODE" = all ]; then
  echo "== unit tests"; $PY -m pytest -q || rc=1
fi
if [ "$MODE" = check ] || [ "$MODE" = all ]; then
  echo "== installation check"; $PY scripts/check_install.py || rc=1
fi
if [ "$MODE" = verify ]; then
  echo "== pipeline results"; $PY scripts/verify_pipeline.py || rc=1
fi
[ $rc -eq 0 ] && echo "ALL GOOD" || echo "FAILED (see above)"
exit $rc
