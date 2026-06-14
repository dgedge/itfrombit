#!/bin/bash
# KZ ramp sweep on the Z3-embedded ensemble: d_trapped(cooling rate) + spectra +
# excess energy. r = 2.0 fixed; T: 6.0 -> 0.5 geometric over R sweeps (R = the KZ
# control parameter; 1 sweep ~ 1 tick, so the eta-survival boot window is R <~ 204).
cd "$(dirname "$0")"
PY=/home/dave/tenpy-env/bin/python
OUT=k04_kz_results.jsonl
JOBS=k04_kz_jobs.txt
: > "$JOBS"
for R in 10 25 50 100 200 400 800 1600 3200 6400; do
  for rep in 1 2 3 4 5 6 7 8 9 10; do
    echo "$PY k04_embedded_sweep.py 4 2.0 1.0 0.5 ramp $R $rep" >> "$JOBS"
  done
  for rep in 1 2 3 4 5 6; do
    echo "$PY k04_embedded_sweep.py 6 2.0 1.0 0.5 ramp $R $rep" >> "$JOBS"
  done
done
for R in 50 200 800 3200; do
  for rep in 1 2 3 4; do
    echo "$PY k04_embedded_sweep.py 8 2.0 1.0 0.5 ramp $R $rep" >> "$JOBS"
  done
done
echo "jobs: $(wc -l < "$JOBS")"
xargs -a "$JOBS" -P 22 -I CMD bash -c 'CMD' >> "$OUT" 2>k04_kz_errors.log
echo "DONE $(wc -l < "$OUT")"
