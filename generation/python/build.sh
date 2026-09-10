#!/usr/bin/env bash
# Builds libmsingen.so from msingen.c + fft.c for use by multisine.py (ctypes).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GEN_DIR="$(dirname "$SCRIPT_DIR")"

gcc -D__PYCLIB -fPIC -shared -O2 -o "$SCRIPT_DIR/libmsingen.so" \
    "$GEN_DIR/msingen.c" "$GEN_DIR/fft.c" -lm

echo "Built $SCRIPT_DIR/libmsingen.so"
