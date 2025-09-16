#!/usr/bin/env bash
# Run the hspvm model over specified subjects or all subjects with 8-character IDs
#
# Usage:
#   ./z02_estimateHSPVM.sh [subject_id...]
#
# Environment variables:
#   ADAPTDELTA   Target acceptance rate (default: 0.90)
#   DATA_DIR     Path to data directory (default: ../../Data)
#   HSPVM_BIN    Path to the hspvm executable (default: ./hspvm)
#   NUM_CHAINS   Number of chains for hspvm (default: 8)
#   NUM_THREADS  Number of threads for hspvm (default: 8)
#   INITS        Path to initialization file (default: ./Inits/Init.json)
#   DRY_RUN      If set to 1, only print commands without running them (default: 0)
#   OVERWRITE    If set to 1, overwrite existing output files (default: 0)

set -euo pipefail
IFS=$'\n\t'

ADAPTDELTA="${ADAPTDELTA:-0.90}"
DATA_DIR="${DATA_DIR:-../../Data}"
HSPVM_BIN="${HSPVM_BIN:-./hspvm}"
NUM_CHAINS="${NUM_CHAINS:-8}"
NUM_THREADS="${NUM_THREADS:-8}"
INITS="${INITS:-./Inits/Init.json}"
DRY_RUN="${DRY_RUN:-0}"
OVERWRITE="${OVERWRITE:-0}"

# Verify prerequisites
if [[ ! -d "$DATA_DIR" ]]; then
  echo "Error: DATA_DIR does not exist: $DATA_DIR" >&2
  exit 1
fi

if [[ ! -x "$HSPVM_BIN" ]]; then
  echo "Warning: hspvm binary not found or not executable at: $HSPVM_BIN" >&2
  echo "         Proceeding anyway; execution will fail unless this is corrected." >&2
fi

# Handle subject IDs
subjects=()
if (( $# > 0 )); then
  # Use provided subject IDs
  subjects=("$@")
else
  # Collect subject IDs: directories whose names are exactly 8 characters
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

  # Prompt for confirmation
  echo "Found ${#subjects[@]} subject(s) to process."
  read -p "Proceed with processing all subjects? [y/N] " -n 1 -r
  echo
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborted by user." >&2
    exit 1
  fi
fi

# Loop over subjects
for sid in "${subjects[@]}"; do
  echo "=== Subject: $sid ==="
  INPUT="$DATA_DIR/$sid/Analysis/A00/InputData.json"
  OUTPUT="$DATA_DIR/$sid/Analysis/A00/Output.csv"

  if [[ ! -f "$INPUT" ]]; then
    echo "  Missing input: $INPUT — skipping." >&2
    continue
  fi

  if [[ "$OVERWRITE" != "1" && -f "$DATA_DIR/$sid/Analysis/A00/Output_1.csv" ]]; then
    echo "  Output already exists: $OUTPUT — skipping." >&2
    continue
  fi

  cmd=("$HSPVM_BIN" "sample" "adapt" "delta=$ADAPTDELTA" \
    "num_chains=$NUM_CHAINS" "num_threads=$NUM_THREADS" \
    "random" "seed=1729" "init=$INITS" "data" "file=$INPUT" "output" "file=$OUTPUT")

  echo "  ${cmd[*]}"
  if [[ "$DRY_RUN" != "1" ]]; then
    "${cmd[@]}"
  fi
  
  echo "=== Done: $sid ==="
done

echo "All done."
