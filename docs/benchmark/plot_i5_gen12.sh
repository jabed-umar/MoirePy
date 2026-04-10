#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

python "$SCRIPT_DIR/plot.py" \
  --input-json "$SCRIPT_DIR/i5_gen12/benchmark_results.json" \
  --output-dir "$SCRIPT_DIR/../../paper/benchmark"

python "$SCRIPT_DIR/plot.py" \
  --input-json "$SCRIPT_DIR/i5_gen12/benchmark_results.json" \
  --output-dir "$SCRIPT_DIR" \
  --save-webp
