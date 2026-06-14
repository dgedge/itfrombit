#!/bin/bash
# Canonical-K04 equilibrium sweep driver (run on deep). Emits JSONL to k04_eq_results.jsonl.
PY=${PY:-/home/dave/tenpy-env/bin/python}
OUT=${OUT:-k04_eq_results.jsonl}
cd "$(dirname "$0")"
: > "$OUT"
JOBS=$(mktemp)
for W4 in 1.5 2 3; do
  S=$($PY -c "print(($W4+1)/3)")
  for N in 24 48; do
    SW=$([ "$N" = 24 ] && echo 3000 || echo 2000)
    for TF in 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.3 2.6 3.0; do
      T=$($PY -c "print(round($S*$TF,4))")
      for ST in cold hot; do
        for R in 1 2 3 4; do
          echo "$N $W4 1 $T $ST $SW $R" >> "$JOBS"
        done
      done
    done
  done
done
# N = 96 spine at fewer T for finite-size trend (ratio 2 only)
for TF in 1.0 1.4 1.8 2.2; do
  T=$($PY -c "print(round($TF,4))")
  for ST in cold hot; do
    for R in 1 2; do
      echo "96 2 1 $T $ST 1200 $R" >> "$JOBS"
    done
  done
done
echo "jobs: $(wc -l < "$JOBS")"
xargs -a "$JOBS" -L1 -P 22 sh -c "$PY k04_eq_sweep.py \$0 \$1 \$2 \$3 \$4 \$5 \$6" >> "$OUT"
echo "DONE $(wc -l < "$OUT") results"
