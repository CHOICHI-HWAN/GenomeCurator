#!/usr/bin/env sh
# GenomeCurator launcher for macOS / Linux: creates .venv on first run, installs, starts the app.
set -e
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "python3 (3.10+) is required: https://www.python.org/downloads/"; exit 1; }
if [ ! -x .venv/bin/python ]; then
  echo "Creating virtual environment (.venv) ..."
  "$PY" -m venv .venv
  .venv/bin/python -m pip install --upgrade pip >/dev/null
  .venv/bin/python -m pip install -e .
fi
exec .venv/bin/python -m genomecurator serve --open "$@"
