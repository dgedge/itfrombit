#!/bin/bash
# Canonical-K04 REFINED overnight sweep (deep). Fine T-windows inside each transition,
# ratio density around the preliminary 1.8, deep equilibration, 8 reps, FSS spine.
# T_c estimate per ratio from first light: T_c ~ 0.85 + 0.55 r; window +-0.35, 12 pts.
PY=${PY:-/home/dave/tenpy-env/bin/python}
OUT=${OUT:-k04_eq2_results.jsonl}
cd "$(dirname "$0")"
: > "$OUT"
JOBS=$(mktemp)
for W4 in 1.5 1.6 1.7 1.8 1.9 2.0 2.25 2.5; do
  TC=$($PY -c "print(0.85+0.55*$W4)")
  for N in 24 48; do
    SW=$([ "$N" = 24 ] && echo 100000 || echo 60000)
    for i in 0 1 2 3 4 5 6 7 8 9 10 11; do
      T=$($PY -c "print(round($TC-0.35+0.0636*$i,4))")
      for ST in cold hot; do
        for R in 1 2 3 4 5 6 7 8; do
          echo "$N $W4 1 $T $ST $SW $R" >> "$JOBS"
        done
      done
    done
  done
done
# N = 96 FSS spine at the three ratios nearest the preliminary pin
for W4 in 1.7 1.8 1.9; do
  TC=$($PY -c "print(0.85+0.55*$W4)")
  for i in 0 2 4 6 8 10; do
    T=$($PY -c "print(round($TC-0.35+0.0636*$i,4))")
    for ST in cold hot; do
      for R in 1 2 3 4; do
        echo "96 $W4 1 $T $ST 30000 $R" >> "$JOBS"
      done
    done
  done
done
echo "jobs: $(wc -l < "$JOBS")"
xargs -a "$JOBS" -L1 -P 22 sh -c "$PY k04_eq_sweep.py \$0 \$1 \$2 \$3 \$4 \$5 \$6" >> "$OUT"
echo "DONE $(wc -l < "$OUT") results"
