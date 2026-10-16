#!/usr/bin/env sh
# GenomeCurator launcher for macOS / Linux.
# Prefers the conda/mamba environment "genomecurator" (created from environment.yml on first run);
# without conda/mamba it falls back to a plain .venv. Extra arguments go to `serve` (e.g. --port 8001).
set -e
cd "$(dirname "$0")/.."
TOOL=""
for t in mamba conda; do command -v "$t" >/dev/null 2>&1 && TOOL="$t" && break; done
if [ -n "$TOOL" ]; then
  ENV_PY="$("$TOOL" run -n genomecurator python -c 'import sys; print(sys.executable)' 2>/dev/null || true)"
  if [ -z "$ENV_PY" ]; then
    echo "Creating conda environment 'genomecurator' from environment.yml ..."
    "$TOOL" env create -f environment.yml
    ENV_PY="$("$TOOL" run -n genomecurator python -c 'import sys; print(sys.executable)')"
  fi
  exec "$ENV_PY" -m genomecurator serve --open "$@"
fi
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "conda/mamba or python3 (3.10+) is required: https://www.python.org/downloads/"; exit 1; }
if [ ! -x .venv/bin/python ]; then
  echo "No conda/mamba found: creating a virtual environment (.venv) ..."
  "$PY" -m venv .venv
  .venv/bin/python -m pip install --upgrade pip >/dev/null
  .venv/bin/python -m pip install -e .
fi
exec .venv/bin/python -m genomecurator serve --open "$@"
