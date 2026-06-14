#!/bin/bash
# Aging/coarsening + percolation sweep (user priorities #2 and #3).
# Aging: d_trapped vs HOLD after boot-window ramps (R = 50, 200) — does coarsening
#        reduce the trapped fraction, and do spectra stay extensive?
# Percolation: L = 8 ramp curve (default hold) — does the crystallised network span?
# (L = 4 percolation is finite-size-trivial; L = 8 is the working size.)
cd "$(dirname "$0")"
PY=/home/dave/tenpy-env/bin/python
OUT=k04_aging_results.jsonl
JOBS=k04_aging_jobs.txt
: > "$JOBS"
for R in 50 200; do
  for H in 20 250 1000 5000 25000; do
    for rep in 1 2 3 4; do
      echo "$PY k04_embedded_sweep.py 6 2.0 1.0 0.5 ramp $R $rep $H" >> "$JOBS"
    done
  done
  for H in 20 250 1000 5000; do
    for rep in 1 2 3; do
      echo "$PY k04_embedded_sweep.py 8 2.0 1.0 0.5 ramp $R $rep $H" >> "$JOBS"
    done
  done
done
for R in 25 100 400 1600 6400; do
  for rep in 1 2 3; do
    echo "$PY k04_embedded_sweep.py 8 2.0 1.0 0.5 ramp $R $rep" >> "$JOBS"
  done
done
echo "jobs: $(wc -l < "$JOBS")"
xargs -a "$JOBS" -P 22 -I CMD bash -c 'CMD' >> "$OUT" 2>k04_aging_errors.log
echo "DONE $(wc -l < "$OUT")"
