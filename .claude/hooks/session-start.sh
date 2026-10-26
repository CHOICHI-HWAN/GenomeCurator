#!/usr/bin/env bash
# Cloud sessions only: install the app (core environment, no conda tools) so `python -m pytest` works.
set -euo pipefail
[ "${CLAUDE_CODE_REMOTE:-}" = "true" ] || exit 0
cd "${CLAUDE_PROJECT_DIR:-$(dirname "$0")/../..}"
pip install -q -e ".[dev]"
