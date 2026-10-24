#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PARENT="$(dirname -- "$ROOT")"
OUTPUT="${PAPER1_OUTPUT:-$PARENT/paper1_complete_from_public_20261022}"
WORKERS="${PAPER1_WORKERS:-60}"
PYTHON="${PAPER1_PYTHON:-python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd -- "$ROOT"
CACHES=()
for candidate in "$PARENT/paper1_from_raw_results_v1" "$PARENT/paper1_public_evidence"; do
  if [[ -d "$candidate" ]]; then CACHES+=(--raw-cache "$candidate"); fi
done
exec "$PYTHON" "$ROOT/reproduce_from_public.py" --workers "$WORKERS" \
  --output "$OUTPUT" "${CACHES[@]}" "$@"
