#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PYTHON="${PAPER1_PYTHON:-python}"
exec "$PYTHON" "$ROOT/paper1_design_oneclick.py" --workers "${PAPER1_WORKERS:-36}" "$@"
