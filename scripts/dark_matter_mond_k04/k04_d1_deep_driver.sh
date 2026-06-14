#!/bin/bash
# D1' deep driver: L-escalation at canonical R + R-scan at L=10.
cd "$(dirname "$0")"
OUT=d1_deep_results.jsonl
: > $OUT
{
for s in $(seq 0 11); do
  for R in 248 350 496 700 992; do echo "8 $R $s"; done
done
for s in $(seq 0 9); do
  for R in 350 496 700; do echo "10 $R $s"; done
done
for s in $(seq 0 7); do echo "12 496 $s"; done
for s in $(seq 0 5); do echo "16 496 $s"; done
} | xargs -P 22 -L 1 sh -c 'python k04_d1_deep_worker.py $0 $1 $2' >> $OUT
wc -l $OUT
