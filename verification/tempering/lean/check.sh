#!/bin/bash
# Compile the Lean files against a local Mathlib v4.33.0 checkout whose build cache
# has been fetched (`lake exe cache get` inside that checkout), without needing
# lake to resolve this project's own dependencies.
#   ./check.sh /path/to/mathlib4 SwapKernel ProductSwap Conjugate PTKernel Check
set -e
M=${1:?path to the mathlib4 checkout}; shift
LP=$M/.lake/build/lib/lean
for d in $M/.lake/packages/*; do LP=$LP:$d/.lake/build/lib/lean; done
OUT=${OUT:-/tmp/tempering-olean}; mkdir -p $OUT/Tempering
export LEAN_PATH=$LP:$OUT
cd "$(dirname "$0")"
for f in "$@"; do
  echo "== $f"
  lean -o $OUT/Tempering/$f.olean Tempering/$f.lean
done
