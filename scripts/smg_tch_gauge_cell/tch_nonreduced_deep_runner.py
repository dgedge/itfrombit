#!/usr/bin/env python3
"""Deep-box runner for heavy non-reduced TCH spin-network rows."""

from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "python_code" / "smg_dmrg_runs"
MODEL_PATH = ROOT / "python_code" / "tch_nonreduced_spin_network_scaling.py"


def load_model():
    spec = importlib.util.spec_from_file_location("tch_nonreduced", MODEL_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


CASES = {
    "minimal_2x1_cutoff4_t02": {
        "case": {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 0.2,
            "n_lanczos": 60,
            "progress_every": 10000,
        },
        "out": RUN_DIR / "tch_nonreduced_deep_minimal_cutoff4_t02.json",
    },
    "minimal_2x1_cutoff4_t4": {
        "case": {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "minimal",
            "beta": 0.5,
            "hopping": 4.0,
            "n_lanczos": 60,
            "progress_every": 10000,
        },
        "out": RUN_DIR / "tch_nonreduced_deep_minimal_cutoff4_t4.json",
    },
    "extended_1x1_cutoff4_t02": {
        "case": {
            "nx": 1,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 0.2,
            "n_lanczos": 60,
            "progress_every": 10000,
        },
        "out": RUN_DIR / "tch_nonreduced_deep_extended_1x1_cutoff4_t02.json",
    },
    "extended_2x1_cutoff4_t02": {
        "case": {
            "nx": 2,
            "ny": 1,
            "states_per_charge": 3,
            "energy_cutoff": 4.0,
            "rep_set": "extended",
            "beta": 0.5,
            "hopping": 0.2,
            "n_lanczos": 50,
            "progress_every": 10000,
        },
        "out": RUN_DIR / "tch_nonreduced_deep_extended_2x1_cutoff4_t02.json",
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cases", nargs="+", default=list(CASES))
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--summary",
        type=Path,
        default=RUN_DIR / "tch_nonreduced_deep_summary.json",
    )
    args = parser.parse_args()

    model = load_model()
    RUN_DIR.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name in args.cases:
        spec = CASES[name]
        out = spec["out"]
        if out.exists() and not args.force:
            print(f"SKIP {name}: {out}", flush=True)
            summary[name] = json.loads(out.read_text())
            continue
        case = dict(spec["case"])
        print(f"START {name}: {case}", flush=True)
        start = time.time()
        row = model.run_case(case)
        row["wall_seconds"] = time.time() - start
        row["case_name"] = name
        out.write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
        summary[name] = row
        print(
            "DONE "
            f"{name}: dim={row['dimension']} "
            f"avgT={row['avg_transitions_per_basis_state']:.6g} "
            f"maxT={row['max_transitions_per_basis_state']} "
            f"gap={row['full_gap']:.12g} "
            f"matter_gap={row['matter_dominated_gap']:.12g} "
            f"wall={row['wall_seconds']:.1f}s",
            flush=True,
        )
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"WROTE {args.summary}", flush=True)


if __name__ == "__main__":
    main()
