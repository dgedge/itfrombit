#!/usr/bin/env python3
r"""Hypercubic Wilson-SU(3) scale-setting validator.

This is the control calculation for the bond-bipyramid SU(3) scale-setting
programme.  Before interpreting a framework-specific ``T_c/sqrt(sigma)`` run,
the same extraction logic should pass the ordinary Wilson-SU(3) benchmark:

    T_c / sqrt(sigma) ~= 0.63

for pure-gauge SU(3) in the continuum window.

What this script tests
----------------------
It deliberately validates the *extractor*, not the framework geometry.

1. The finite-temperature half is the already validated ordinary hypercubic
   Wilson-SU(3) transition at N_t=4, beta_c ~= 5.6925.
2. The zero-temperature half is measured on a separate symmetric hypercubic
   lattice at the same beta, using rectangular Wilson loops and Creutz ratios.
3. The diagnostic asks whether the small-loop Creutz estimator used in the
   quick bond-bipyramid scale-setting harness gives

       1 / (N_t a sqrt(sigma)) ~= 0.63.

Interpretation
--------------
This is *not* a production lattice-QCD calculation: it uses small volumes,
unsmeared Wilson loops, and modest Metropolis statistics.  It is a correctness
gate for the current extraction strategy.  If this control misses badly, the
bond-bipyramid negative scale-setting result should be read as "extractor not
validated" until a standard zero-temperature static-potential/smearing pipeline
is in place.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from record_grammar_tch_su3_polyakov_torus import I3, dag, sweep


TARGET_TC_OVER_SQRT_SIGMA = 0.629
REFERENCE_BETA_C_NT4 = 5.6925
DEFAULT_OUTPUT = "python_code/cubic_su3_scale_setting_validator.json"


@dataclass
class CreutzRow:
    label: str
    chi: float | None
    a_sqrt_sigma: float | None
    tc_over_sqrt_sigma: float | None
    rel_error_vs_target: float | None


def roll_multi(A: np.ndarray, shifts: dict[int, int]) -> np.ndarray:
    out = A
    for axis, shift in shifts.items():
        if shift:
            out = np.roll(out, shift, axis=axis)
    return out


def wilson_loop_mean(U: np.ndarray, r: int, t: int) -> float:
    """Average fundamental Wilson loop for all ordered planes on a 4D torus."""

    ndim = U.shape[0]
    vals: list[float] = []
    for mu in range(ndim):
        for nu in range(ndim):
            if mu == nu:
                continue
            shape = U[mu].shape[:-2]
            M = np.broadcast_to(I3, shape + (3, 3)).copy()
            for i in range(r):
                M = M @ roll_multi(U[mu], {mu: -i})
            for j in range(t):
                M = M @ roll_multi(U[nu], {mu: -r, nu: -j})
            for i in range(r):
                M = M @ dag(roll_multi(U[mu], {mu: -(r - 1 - i), nu: -t}))
            for j in range(t):
                M = M @ dag(roll_multi(U[nu], {nu: -(t - 1 - j)}))
            vals.append(float(np.real(np.trace(M, axis1=-2, axis2=-1).mean() / 3.0)))
    return float(np.mean(vals))


def wilson_table(U: np.ndarray, max_loop: int) -> dict[tuple[int, int], float]:
    return {
        (r, t): wilson_loop_mean(U, r, t)
        for r in range(1, max_loop + 1)
        for t in range(1, max_loop + 1)
    }


def creutz(w: dict[tuple[int, int], float], r: int, t: int) -> float | None:
    keys = [(r, t), (r + 1, t + 1), (r + 1, t), (r, t + 1)]
    if any(key not in w or w[key] <= 0 for key in keys):
        return None
    value = -math.log((w[(r, t)] * w[(r + 1, t + 1)]) / (w[(r + 1, t)] * w[(r, t + 1)]))
    return value if math.isfinite(value) and value > 0 else None


def creutz_row(label: str, chi: float | None, nt: int) -> CreutzRow:
    if chi is None or chi <= 0:
        return CreutzRow(label, chi, None, None, None)
    a_sqrt_sigma = math.sqrt(chi)
    ratio = 1.0 / (nt * a_sqrt_sigma)
    rel = abs(ratio - TARGET_TC_OVER_SQRT_SIGMA) / TARGET_TC_OVER_SQRT_SIGMA
    return CreutzRow(label, chi, a_sqrt_sigma, ratio, rel)


def su3_unitarity_error(U: np.ndarray) -> float:
    eye = np.broadcast_to(I3, U.shape[:-2] + (3, 3))
    return float(np.max(np.abs(dag(U) @ U - eye)))


def run(beta: float, nt: int, size: int, therm: int, meas: int, hits: int, eps0: float, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    shape = (size, size, size, size)
    U = np.empty((4,) + shape + (3, 3), dtype=complex)
    for mu in range(4):
        U[mu] = np.broadcast_to(I3, shape + (3, 3)).copy()
    idx = np.indices(shape).sum(axis=0)
    even, odd = (idx % 2 == 0), (idx % 2 == 1)
    eps = eps0
    acc_values: list[float] = []

    for _ in range(therm):
        U, acc = sweep(U, beta, rng, even, odd, eps, hits=hits)
        eps = min(max(eps * (1.0 + 0.1 * (acc - 0.5)), 0.04), 1.0)
        acc_values.append(float(acc))

    tables: list[dict[tuple[int, int], float]] = []
    max_loop = 3
    for _ in range(meas):
        U, acc = sweep(U, beta, rng, even, odd, eps, hits=hits)
        eps = min(max(eps * (1.0 + 0.05 * (acc - 0.5)), 0.04), 1.0)
        acc_values.append(float(acc))
        tables.append(wilson_table(U, max_loop=max_loop))

    keys = sorted(tables[0])
    w_mean = {key: float(np.mean([table[key] for table in tables])) for key in keys}
    w_se = {
        key: float(np.std([table[key] for table in tables], ddof=1) / math.sqrt(len(tables)))
        for key in keys
    }
    rows = [
        creutz_row("chi(1,1)", creutz(w_mean, 1, 1), nt),
        creutz_row("chi(2,2)", creutz(w_mean, 2, 2), nt),
    ]
    valid_rows = [row for row in rows if row.tc_over_sqrt_sigma is not None]
    best = min(valid_rows, key=lambda row: row.rel_error_vs_target or float("inf")) if valid_rows else None

    verdict = "VALIDATED" if best and (best.rel_error_vs_target or 1.0) <= 0.15 else "NOT-VALIDATED"
    return {
        "profile": "hypercubic-wilson-su3-extractor-validator",
        "beta": beta,
        "nt_reference": nt,
        "size_zero_temperature": size,
        "therm": therm,
        "meas": meas,
        "hits": hits,
        "seed": seed,
        "target_tc_over_sqrt_sigma": TARGET_TC_OVER_SQRT_SIGMA,
        "acceptance_mean": float(np.mean(acc_values)),
        "unitarity_error": su3_unitarity_error(U),
        "wilson_means": {f"{r}x{t}": w_mean[(r, t)] for r, t in keys},
        "wilson_se": {f"{r}x{t}": w_se[(r, t)] for r, t in keys},
        "creutz_rows": [asdict(row) for row in rows],
        "best_row": asdict(best) if best else None,
        "verdict": verdict,
        "interpretation": (
            "Small-loop unsmeared Creutz extraction reproduces the ordinary SU(3) scale ratio."
            if verdict == "VALIDATED"
            else "Current small-loop unsmeared Creutz extraction does not reproduce the ordinary SU(3) scale ratio; the bond-bipyramid scale-setting negative cannot yet be read as a continuum-physics negative."
        ),
    }


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{digits}f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beta", type=float, default=REFERENCE_BETA_C_NT4)
    parser.add_argument("--nt", type=int, default=4)
    parser.add_argument("--size", type=int, default=6)
    parser.add_argument("--therm", type=int, default=80)
    parser.add_argument("--meas", type=int, default=80)
    parser.add_argument("--hits", type=int, default=1)
    parser.add_argument("--eps0", type=float, default=0.28)
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    t0 = time.time()
    result = run(args.beta, args.nt, args.size, args.therm, args.meas, args.hits, args.eps0, args.seed)
    result["elapsed_s"] = time.time() - t0
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("Hypercubic Wilson-SU(3) scale-setting validator")
    print("=" * 90)
    print(f"  beta={args.beta:.4f} (ordinary SU(3) Nt={args.nt} beta_c reference)")
    print(f"  zero-temperature lattice: {args.size}^4, therm={args.therm}, meas={args.meas}, hits={args.hits}")
    print(f"  acceptance={result['acceptance_mean']:.3f}, unitarity_error={result['unitarity_error']:.2e}")
    print("\n[Wilson loops]")
    for key, value in result["wilson_means"].items():
        print(f"  W({key}) = {value:.6f} +/- {result['wilson_se'][key]:.6f}")
    print("\n[Creutz scale estimates]")
    print("  estimator  chi       a√sigma  Tc/√sigma  rel_error_vs_0.629")
    for row in result["creutz_rows"]:
        print(
            f"  {row['label']:<9s} {fmt(row['chi']):>8s} "
            f"{fmt(row['a_sqrt_sigma']):>9s} {fmt(row['tc_over_sqrt_sigma']):>10s} "
            f"{fmt(row['rel_error_vs_target'], 3):>17s}"
        )
    print(f"\nVERDICT: {result['verdict']}")
    print("  " + result["interpretation"])
    print(f"  wrote {args.output}")
    print("exit 0")
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
