#!/bin/bash
cd "$(dirname "$0")"
OUT=d1_winding_results.jsonl
: > $OUT
{
for s in $(seq 0 11); do echo "8 496 $s"; done
for s in $(seq 0 9);  do echo "10 496 $s"; done
for s in $(seq 0 7);  do echo "12 496 $s"; done
for s in $(seq 0 5);  do echo "16 496 $s"; done
for s in $(seq 0 5);  do for R in 248 992; do echo "10 $R $s"; done; done
} | xargs -P 22 -L 1 sh -c 'python k04_d1_winding_worker.py $0 $1 $2' >> $OUT
wc -l $OUT
