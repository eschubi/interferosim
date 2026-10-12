#!/usr/bin/env bash
# setup_venv.sh -- create the project virtual environment .venv (Linux/macOS)
# Usage:  ./setup_venv.sh            (PYTHON=python3.12 ./setup_venv.sh to choose)
#         ./setup_venv.sh --recreate --no-tests
# Version 1.0.0 (2026-10-12)
set -euo pipefail
cd "$(dirname "$0")"

RECREATE=0; TESTS=1
for a in "$@"; do
  case "$a" in
    --recreate) RECREATE=1 ;;
    --no-tests) TESTS=0 ;;
    *) echo "unknown option $a"; exit 2 ;;
  esac
done

PY="${PYTHON:-}"
if [ -z "$PY" ]; then
  for c in python3.13 python3.12 python3.11 python3; do
    if command -v "$c" >/dev/null 2>&1 && \
       "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)'; then
      PY="$c"; break
    fi
  done
fi
[ -n "$PY" ] || { echo "no Python >= 3.11 found (set PYTHON=...)"; exit 1; }
echo "Python: $PY ($("$PY" -c 'import sys; print(sys.version.split()[0])'))"

[ "$RECREATE" = 1 ] && rm -rf .venv
[ -x .venv/bin/python ] || "$PY" -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e ".[notebooks,test]"
.venv/bin/python -m ipykernel install --user --name interferosim --display-name "Python (interferosim)"
[ "$TESTS" = 1 ] && .venv/bin/python -m pytest -q

echo
echo "done. Activate with:  source .venv/bin/activate"
