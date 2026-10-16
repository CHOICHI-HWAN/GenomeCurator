#!/usr/bin/env bash
# Apache Superset for GenomeCurator. Needs a conda env with apache-superset (see docs/SUPERSET.md).
#   scripts/superset.sh init     one-time: metadata DB, admin user (password -> superset/home/admin_password.txt)
#   scripts/superset.sh run      start Superset on $SUPERSET_PORT (default 8089)
#   scripts/superset.sh bootstrap  create/refresh the connection, dataset, charts and dashboard
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ENV_NAME="${SUPERSET_ENV:-gc-superset}"
PREFIX="$(conda run -n "$ENV_NAME" python -c 'import sys;print(sys.prefix)' 2>/dev/null || mamba run -n "$ENV_NAME" python -c 'import sys;print(sys.prefix)')"
export PATH="$PREFIX/bin:$PATH"
export SUPERSET_HOME="$ROOT/superset/home"
export SUPERSET_CONFIG_PATH="$ROOT/superset/superset_config.py"
export FLASK_APP=superset
export GC_DB="${GC_DB:-${GC_DATA_DIR:-$ROOT/data}/genomecurator.db}"
PORT="${SUPERSET_PORT:-8089}"
mkdir -p "$SUPERSET_HOME"
case "${1:-run}" in
  init)
    superset db upgrade
    if [ ! -f "$SUPERSET_HOME/admin_password.txt" ]; then
      python -c 'import secrets;print(secrets.token_urlsafe(12))' > "$SUPERSET_HOME/admin_password.txt"
      chmod 600 "$SUPERSET_HOME/admin_password.txt"
    fi
    superset fab create-admin --username admin --firstname Admin --lastname User --email admin@localhost \
      --password "$(cat "$SUPERSET_HOME/admin_password.txt")" || true
    superset init ;;
  run)
    exec superset run -p "$PORT" -h 127.0.0.1 --with-threads ;;
  bootstrap)
    python "$ROOT/superset/bootstrap.py" --url "http://127.0.0.1:$PORT" --password-file "$SUPERSET_HOME/admin_password.txt" --db "$GC_DB" ;;
  *) echo "usage: $0 init|run|bootstrap"; exit 2 ;;
esac
