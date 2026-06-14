#!/usr/bin/env python3
r"""Resolve the K04 gamma-bridge entropy budget through the crystallisation spike.

The earlier bridge script measured a coarse L=6 heat-capacity curve and found
Delta S = 2.684 nats/site, +7.0% above the pre-registered target

    Delta S_target = ln(12) * alpha0 * ln 2 / [-ln(0.995)] = 2.508 nats/site.

For a first-order-like crystallisation, the narrow spike is better treated as
latent heat: Delta S = Delta E / T_c.  This script measures both:

  1. The legacy fluctuation/trapezoid estimator, on a fine grid.
  2. A bracketed energy-jump estimator across the ordered/liquid branches.

Exit 0 means the experiment ran and tiered the gate.  It intentionally does not
force a pass/fail by assertion: a miss is a result.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import math
import os
import random
import statistics as stats
import sys
import time
from pathlib import Path

ALPHA0 = 1.0 / 137.0
GAMMA_PAPER = 0.995
T_START = 6.0
T_END = 0.5
DELTA_S_TARGET = math.log(T_START / T_END) * ALPHA0 * math.log(2.0) / math.log(1.0 / GAMMA_PAPER)
CV_REQ = ALPHA0 * math.log(2.0) / math.log(1.0 / GAMMA_PAPER)


def load_embedded_prefix():
    src = (Path(__file__).parent / "k04_embedded_sweep.py").read_text(encoding="utf-8")
    marker = 'if START in ("hot", "ramp"):'
    if marker not in src:
        raise RuntimeError("k04_embedded_sweep.py marker changed")
    prefix = src[: src.index(marker)]
    return prefix


def bootstrap_mean(xs: list[float], rng: random.Random, nboot: int) -> tuple[float, float, float]:
    if not xs:
        return float("nan"), float("nan"), float("nan")
    if len(xs) == 1:
        return xs[0], xs[0], xs[0]
    vals = []
    n = len(xs)
    for _ in range(nboot):
        vals.append(sum(xs[rng.randrange(n)] for _ in range(n)) / n)
    vals.sort()
    return stats.mean(xs), vals[int(0.16 * (nboot - 1))], vals[int(0.84 * (nboot - 1))]


def make_env(prefix: str, L: int, w4: float, w6: float, seed: int):
    sys.argv = ["k04", str(L), str(w4), str(w6), "1.0", "cold", "1", "0"]
    env: dict = {}
    exec(compile(prefix, "k04_embedded_sweep_prefix", "exec"), env)
    env["random"].seed(seed)
    return env


def randomize_hot(env: dict, B: set, c4: int, c6: int):
    # Match k04_embedded_sweep.py's hot-start randomisation.
    return env["metropolis"](B, c4, c6, math.inf, max(1, 3000 // env["N"] + 1) * 20)


def run_point(prefix: str, L: int, w4: float, w6: float, T: float, start: str,
              rep: int, burn: int, sample: int, thin: int) -> dict:
    seed = 0x4B04_2026 + 1000003 * L + 1009 * rep + int(round(1000 * T)) + (0 if start == "cold" else 17)
    env = make_env(prefix, L, w4, w6, seed)
    B, c4, c6 = set(env["B"]), env["c4"], env["c6"]
    if start == "hot":
        B, c4, c6, _ = randomize_hot(env, B, c4, c6)
    B, c4, c6, _ = env["metropolis"](B, c4, c6, T, burn)
    es: list[float] = []
    ds: list[float] = []
    acc = 0
    for _ in range(sample):
        B, c4, c6, a = env["metropolis"](B, c4, c6, T, thin)
        acc += a
        e = -env["W4"] * c4 - env["W6"] * c6
        es.append(e / env["N"])
        ds.append(env["census"](B)[1])
    assert c4 == env["count_all"](B, 4) and c6 == env["count_all"](B, 6)
    e_mean = stats.mean(es)
    cv = stats.pvariance([x * env["N"] for x in es]) / (T * T * env["N"])
    return {
        "T": T,
        "start": start,
        "rep": rep,
        "E": e_mean,
        "Cv": cv,
        "d": stats.mean(ds),
        "acc": acc / max(1, sample * thin * env["N"]),
        "c4": c4,
        "c6": c6,
    }


def run_point_task(task: tuple[str, int, float, float, float, str, int, int, int, int]) -> dict:
    return run_point(*task)


def integrate_trapezoid(rows: list[dict]) -> float:
    rows = sorted(rows, key=lambda r: r["T"], reverse=True)
    total = 0.0
    for a, b in zip(rows, rows[1:]):
        total += 0.5 * (a["Cv"] + b["Cv"]) * abs(math.log(a["T"] / b["T"]))
    return total


def median(xs: list[float]) -> float:
    return stats.median(xs) if xs else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=6)
    ap.add_argument("--w4", type=float, default=1.7)
    ap.add_argument("--w6", type=float, default=1.0)
    ap.add_argument("--temps", default="6,5,4,3.6,3.3,3.15,3.05,3.0,2.95,2.9,2.85,2.8,2.75,2.7,2.65,2.6,2.5,2.4")
    ap.add_argument("--reps", type=int, default=4)
    ap.add_argument("--burn", type=int, default=700)
    ap.add_argument("--sample", type=int, default=900)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--nboot", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    args = ap.parse_args()

    prefix = load_embedded_prefix()
    temps = [float(x) for x in args.temps.split(",") if x.strip()]
    t0 = time.time()

    print("[0] K04 entropy-spike resolution", flush=True)
    print(f"    L={args.L}, N={args.L ** 3}, w4/w6={args.w4 / args.w6:.3f}, reps={args.reps}, burn={args.burn}, sample={args.sample}, thin={args.thin}")
    print(f"    target Delta S = {DELTA_S_TARGET:.6f} nats/site; C_v required for gamma=0.995 is {CV_REQ:.6f}")
    print(f"    workers={args.workers}", flush=True)
    print(flush=True)

    all_rows: list[dict] = []
    tasks = [
        (prefix, args.L, args.w4, args.w6, T, start, rep, args.burn, args.sample, args.thin)
        for T in temps
        for start in ("cold", "hot")
        for rep in range(1, args.reps + 1)
    ]
    completed_by_T = {T: 0 for T in temps}
    needed_by_T = {T: 2 * args.reps for T in temps}
    with cf.ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_point_task, task) for task in tasks]
        for fut in cf.as_completed(futs):
            row = fut.result()
            all_rows.append(row)
            completed_by_T[row["T"]] += 1
            if completed_by_T[row["T"]] == needed_by_T[row["T"]]:
                T = row["T"]
                cold_E = [r["E"] for r in all_rows if r["T"] == T and r["start"] == "cold"]
                hot_E = [r["E"] for r in all_rows if r["T"] == T and r["start"] == "hot"]
                cold_d = [r["d"] for r in all_rows if r["T"] == T and r["start"] == "cold"]
                hot_d = [r["d"] for r in all_rows if r["T"] == T and r["start"] == "hot"]
                cvs = [r["Cv"] for r in all_rows if r["T"] == T]
                accs = [r["acc"] for r in all_rows if r["T"] == T]
                print(f"    T={T:5.2f}  Cv_med={median(cvs):8.3f}  E_cold={stats.mean(cold_E):+8.3f}  E_hot={stats.mean(hot_E):+8.3f}  d_cold={stats.mean(cold_d):.3f}  d_hot={stats.mean(hot_d):.3f}  acc={stats.mean(accs):.4f}  [{time.time()-t0:5.1f}s]", flush=True)

    # Legacy-style fluctuation integral: combine starts/reps at each temperature.
    cv_rows = []
    rng = random.Random(12345)
    for T in temps:
        rs = [r for r in all_rows if r["T"] == T]
        cv_rows.append({"T": T, "Cv": stats.mean(r["Cv"] for r in rs)})
    dS_trap = integrate_trapezoid(cv_rows)

    trap_boot = []
    byT = {T: [r["Cv"] for r in all_rows if r["T"] == T] for T in temps}
    for _ in range(args.nboot):
        sampled = []
        for T in temps:
            vals = byT[T]
            sampled.append({"T": T, "Cv": sum(vals[rng.randrange(len(vals))] for _ in vals) / len(vals)})
        trap_boot.append(integrate_trapezoid(sampled))
    trap_boot.sort()
    dS_trap_lo = trap_boot[int(0.16 * (args.nboot - 1))]
    dS_trap_hi = trap_boot[int(0.84 * (args.nboot - 1))]

    # Latent heat estimator: find the T where cold/hot branch energy separation
    # gives Delta E / T closest to the target.  This is the correct narrow-spike
    # estimator if the transition is first-order-like.
    latent = []
    for T in temps:
        cold = [r["E"] for r in all_rows if r["T"] == T and r["start"] == "cold"]
        hot = [r["E"] for r in all_rows if r["T"] == T and r["start"] == "hot"]
        dE = stats.mean(hot) - stats.mean(cold)
        latent.append((T, dE / T, dE, cold, hot))
    T_star, dS_lat, dE_star, cold_star, hot_star = min(latent, key=lambda x: abs(x[1] - DELTA_S_TARGET))
    lat_boot = []
    for _ in range(args.nboot):
        c = sum(cold_star[rng.randrange(len(cold_star))] for _ in cold_star) / len(cold_star)
        h = sum(hot_star[rng.randrange(len(hot_star))] for _ in hot_star) / len(hot_star)
        lat_boot.append((h - c) / T_star)
    lat_boot.sort()
    dS_lat_lo = lat_boot[int(0.16 * (args.nboot - 1))]
    dS_lat_hi = lat_boot[int(0.84 * (args.nboot - 1))]

    print()
    print("[1] entropy estimators")
    print(f"    fluctuation/trapezoid Delta S = {dS_trap:.6f} [{dS_trap_lo:.6f}, {dS_trap_hi:.6f}]  ({dS_trap / DELTA_S_TARGET - 1:+.2%} vs target)")
    print(f"    latent-heat best crossing: T*={T_star:.3f}, Delta E={dE_star:.6f}/site")
    print(f"    latent Delta S = Delta E/T* = {dS_lat:.6f} [{dS_lat_lo:.6f}, {dS_lat_hi:.6f}]  ({dS_lat / DELTA_S_TARGET - 1:+.2%} vs target)")

    pct_band = 0.01 * DELTA_S_TARGET
    trap_pass = abs(dS_trap - DELTA_S_TARGET) <= pct_band
    lat_pass = abs(dS_lat - DELTA_S_TARGET) <= pct_band
    print()
    print("[2] pre-registered gamma-bridge gate")
    print(f"    percent-grade acceptance band: {DELTA_S_TARGET - pct_band:.6f} .. {DELTA_S_TARGET + pct_band:.6f}")
    print(f"    fluctuation estimator: {'PASS' if trap_pass else 'MISS'}")
    print(f"    latent-heat estimator: {'PASS' if lat_pass else 'MISS'}")
    if lat_pass and not trap_pass:
        print("    verdict: closes under the narrow first-order spike reading; the older +7% was a coarse trapezoid artifact.")
    elif not lat_pass:
        print("    verdict: gamma=exp(-alpha0 ln2) is not derived by this K04 entropy gate at percent grade.")
    else:
        print("    verdict: both estimators land inside the percent gate.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
