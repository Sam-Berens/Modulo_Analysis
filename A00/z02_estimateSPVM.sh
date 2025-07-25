#!/usr/bin/env bash
# Run spvm over all subjects with 8-character IDs under ../../Data
#
# Usage:
#   ./z02_estimateSPVM.sh [DATA_DIR]
#
# Environment variables:
#   SPVM_BIN    Path to the spvm executable (default: ./spvm)
#   NUM_CHAINS  Number of chains for spvm (default: 8)
#   NUM_THREADS Number of threads for spvm (default: 8)
#   DRY_RUN     If set to 1, only print commands without running them (default: 0)

set -euo pipefail
IFS=$'\n\t'

DATA_DIR="${1:-../../Data}"
SPVM_BIN="${SPVM_BIN:-./spvm}"
NUM_CHAINS="${NUM_CHAINS:-8}"
NUM_THREADS="${NUM_THREADS:-8}"
DRY_RUN="${DRY_RUN:-0}"

# Verify prerequisites
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: DATA_DIR does not exist: $DATA_DIR" >&2
  exit 1
fi

if [[ ! -x "$SPVM_BIN" ]]; then
  echo "Warning: spvm binary not found or not executable at: $SPVM_BIN" >&2
  echo "         Proceeding anyway; execution will fail unless this is corrected." >&2
fi

# Collect subject IDs: directories whose names are exactly 8 characters
subjects=()
shopt -s nullglob
for path in "$DATA_DIR"/*; do
  name=$(basename "$path")
  if [[ -d "$path" && ${#name} -eq 8 ]]; then
    subjects+=("$name")
  fi
done
shopt -u nullglob

if (( ${#subjects[@]} == 0 )); then
  echo "No 8-character subject directories found under $DATA_DIR. Nothing to do." >&2
  exit 0
fi

echo "Found ${#subjects[@]} subject(s)."

# Loop over subjects
for sid in "${subjects[@]}"; do
  echo "=== Subject: $sid ==="
  INPUTS="$DATA_DIR/$sid/Behavioural/A00/InputData"
  OUTPUTS="$DATA_DIR/$sid/Behavioural/A00/Outputs"

  if [[ ! -d "$INPUTS" ]]; then
    echo "  Warning: inputs directory missing: $INPUTS — skipping subject." >&2
    continue
  fi

  mkdir -p "$OUTPUTS"

  # pairId: 0..35 (since 6^2 - 1 = 35)
  for pid in $(seq 0 35); do
    printf -v in_file  "%s/P%02d.json" "$INPUTS"  "$pid"
    printf -v out_file "%s/P%02d.csv"  "$OUTPUTS" "$pid"

    if [[ ! -f "$in_file" ]]; then
      echo "  [$(printf %02d "$pid")] Missing input: $in_file — skipping." >&2
      continue
    fi

    # spvm expects tokens like: data file=... output file=...
    cmd=("$SPVM_BIN" "sample" "num_chains=$NUM_CHAINS" "num_threads=$NUM_THREADS" \
         "data" "file=$in_file" "output" "file=$out_file")

    echo "  [$(printf %02d "$pid")] ${cmd[*]}"
    if [[ "$DRY_RUN" != "1" ]]; then
      "${cmd[@]}"
    fi
  done

echo "=== Done: $sid ==="
done

echo "All done."
