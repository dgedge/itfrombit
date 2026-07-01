#!/usr/bin/env bash
set -euo pipefail

# Deep/overnight launcher for the full-SU(3) bond-bipyramid scale-setting run.
#
# Examples:
#   PROFILE=deep-preflight MAX_JOBS=2 ./python_code/record_grammar_tch_su3_bulk_scale_setting_deep_driver.sh
#   PROFILE=deep MAX_JOBS=3 NTS="2 3 4 5 6" ./python_code/record_grammar_tch_su3_bulk_scale_setting_deep_driver.sh
#
# Each Nt writes independent JSON, JSONL checkpoints, and a log file under
# python_code/su3_bulk_scale_setting_deep/.  This makes the run resumable by
# rerunning only the failed Nt value.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ -z "${PYTHON:-}" ]]; then
  for candidate in \
    /Users/davidelliman/bin/py13_7/bin/python \
    /home/dave/tenpy-env/bin/python \
    /home/dave/py311/bin/python \
    python3 \
    python; do
    if command -v "$candidate" >/dev/null 2>&1; then
      PYTHON="$(command -v "$candidate")"
      break
    fi
  done
fi
if [[ -z "${PYTHON:-}" ]]; then
  echo "No Python interpreter found. Set PYTHON=/path/to/python." >&2
  exit 2
fi
PROFILE="${PROFILE:-deep-preflight}"
MAX_JOBS="${MAX_JOBS:-2}"
SEED="${SEED:-20260630}"
OUTDIR="${OUTDIR:-python_code/su3_bulk_scale_setting_deep}"

if [[ -z "${NTS:-}" ]]; then
  if [[ "$PROFILE" == "deep" ]]; then
    NTS="2 3 4 5 6"
  else
    NTS="2 3 4"
  fi
fi

mkdir -p "$OUTDIR"
echo "profile=$PROFILE max_jobs=$MAX_JOBS seed=$SEED nts=$NTS outdir=$OUTDIR"
echo "python=$PYTHON"

running_jobs() {
  jobs -rp | wc -l | tr -d ' '
}

for nt in $NTS; do
  while [[ "$(running_jobs)" -ge "$MAX_JOBS" ]]; do
    sleep 10
  done

  out="$OUTDIR/${PROFILE}_Nt${nt}.json"
  jsonl="$OUTDIR/${PROFILE}_Nt${nt}.jsonl"
  log="$OUTDIR/${PROFILE}_Nt${nt}.log"
  echo "launch Nt=$nt -> $log"
  "$PYTHON" -u python_code/record_grammar_tch_su3_bulk_scale_setting.py \
    --profile "$PROFILE" \
    --nt "$nt" \
    --seed "$((SEED + 100 * nt))" \
    --output "$out" \
    --jsonl "$jsonl" \
    > "$log" 2>&1 &
done

wait
echo "all Nt jobs complete"
echo "logs/results in $OUTDIR"
