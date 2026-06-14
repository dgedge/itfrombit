#!/bin/bash
# Coarse ordering-window bracket for the Z3-EMBEDDED K04 ensemble (replaces the
# superseded configuration-model toy). L in {4,6}; r in {1.5, 2.0, 2.5}; 10 T points
# spanning the expected window; hot+cold; 2 reps. Fine grid follows once bracketed.
cd "$(dirname "$0")"
PY=/home/dave/tenpy-env/bin/python
OUT=k04_embedded_results.jsonl
JOBS=k04_embedded_jobs.txt
: > "$JOBS"
for r in 1.5 2.0 2.5; do
  for T in 1.5 2.0 2.5 3.0 3.5 4.0 5.0 6.0 7.0 8.0; do
    for st in cold hot; do
      for rep in 1 2; do
        echo "$PY k04_embedded_sweep.py 4 $r 1.0 $T $st 20000 $rep" >> "$JOBS"
        echo "$PY k04_embedded_sweep.py 6 $r 1.0 $T $st 8000  $rep" >> "$JOBS"
      done
    done
  done
done
echo "jobs: $(wc -l < "$JOBS")"
xargs -a "$JOBS" -L1 -P 22 -I CMD bash -c 'CMD' >> "$OUT" 2>k04_embedded_errors.log
echo "DONE $(wc -l < "$OUT")"
