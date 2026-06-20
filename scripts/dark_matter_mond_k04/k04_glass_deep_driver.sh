#!/bin/bash
# K04 constraint-glass DEEP scaling driver: chi4(L) cooperative length + jamming + RFOT.
# One rep per job, fanned across cores (matches k04_d1_deep_driver.sh).
# Launch on deep (24 cores), venv active:
#   nohup bash k04_glass_deep_driver.sh > k04_glass_deep.log 2>&1 &
set -u
cd "$(dirname "$0")"
OUT=k04_glass_deep_results.jsonl
PROCS=${PROCS:-22}
REPS_T=${REPS_T:-192}       # T (chi4 + jamming) reps at L=8,10,12
REPS_T16=${REPS_T16:-64}    # T reps at L=16 (more expensive)
REPS_R=${REPS_R:-48}        # R (replica overlap) reps at L=8,10,12
: > "$OUT"
echo "[start] $(date) PROCS=$PROCS REPS_T=$REPS_T REPS_T16=$REPS_T16 REPS_R=$REPS_R"
{
  for L in 8 10 12; do for r in $(seq 0 $((REPS_T - 1))); do echo "T $L $r"; done; done
  for r in $(seq 0 $((REPS_T16 - 1))); do echo "T 16 $r"; done
  for L in 8 10 12; do for r in $(seq 0 $((REPS_R - 1))); do echo "R $L $r"; done; done
} | xargs -P "$PROCS" -L 1 sh -c 'python k04_glass_deep_worker.py "$0" "$1" "$2"' >> "$OUT"
echo "[done] $(date)  lines=$(wc -l < "$OUT")"
