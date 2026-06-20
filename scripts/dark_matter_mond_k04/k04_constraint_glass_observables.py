#!/usr/bin/env python3
r"""K04 constraint-glass observable pilot.

This is the next layer after k04_constraint_glass_kcm_map.py.  It measures the
standard KCM/glass observables on the faithful embedded K04 plaquette dynamics:

  * activity large-deviation pilot:
        K_t = accepted plaquette flips,
        psi(s) = t^{-1} log <exp(-s K_t)>;

  * persistence and four-point susceptibility:
        p_i(t) = 1 if site i's incident-bond neighbourhood has never changed,
        P(t) = <N^{-1} sum_i p_i(t)>,
        chi_4(t) = N Var[N^{-1} sum_i p_i(t)];

  * bootstrap/jamming closure:
        greedily apply all legal non-uphill plaquette moves until none remain;

  * coupled-replica overlap:
        evolve a second copy with an -epsilon Q(B,B0) bond-overlap bias.

The purpose is classification and experiment design, not cosmological
prediction.  Default parameters are small enough for a laptop smoke test.  Use
larger L/reps on deep for scaling.

Run:
    ~/bin/py13_7/bin/python python_code/k04_constraint_glass_observables.py

Optional:
    ~/bin/py13_7/bin/python python_code/k04_constraint_glass_observables.py \
        --Ls 4,6 --reps 12 --sweeps 18 --out k04_constraint_glass_observables.json
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import fmean, pstdev

import numpy as np


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_code"))
import k04_constraint_glass_kcm_map as k04  # noqa: E402


AX = k04.AX


def parse_ints(text: str) -> list[int]:
    return [int(x) for x in text.split(",") if x.strip()]


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in text.split(",") if x.strip()]


def all_bonds(L: int) -> set[tuple[tuple[int, int, int], int]]:
    return {(p, a) for p in k04.sites(L) for a in range(3)}


def endpoints(edge: tuple[tuple[int, int, int], int], L: int) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    p, a = edge
    return p, k04.wrap(k04.add(p, AX[a]), L)


def touched_sites(move, L: int) -> set[tuple[int, int, int]]:
    rm, ad = move
    out: set[tuple[int, int, int]] = set()
    for edge in (*rm, *ad):
        out.update(endpoints(edge, L))
    return out


def apply_move(B: set[tuple[tuple[int, int, int], int]], move) -> None:
    rm, ad = move
    B.discard(rm[0])
    B.discard(rm[1])
    B.add(ad[0])
    B.add(ad[1])


def incident_signature(B: set[tuple[tuple[int, int, int], int]], p: tuple[int, int, int], L: int) -> tuple[int, ...]:
    bits = []
    for a in range(3):
        bits.append(1 if (p, a) in B else 0)
        m = k04.wrap(k04.add(p, (-AX[a][0], -AX[a][1], -AX[a][2])), L)
        bits.append(1 if (m, a) in B else 0)
    return tuple(bits)


def site_signatures(B: set[tuple[tuple[int, int, int], int]], L: int) -> dict[tuple[int, int, int], tuple[int, ...]]:
    return {p: incident_signature(B, p, L) for p in k04.sites(L)}


def choose_move(B: set[tuple[tuple[int, int, int], int]], L: int, rng: random.Random):
    moves = list(k04.flippable_moves(B, L))
    if not moves:
        return None
    return rng.choice(moves)


def prepare_state(L: int, seed: int, random_sweeps: int, prep_T: float, prep_sweeps: int):
    B = k04.randomize(set(k04.tiling(L)), L, sweeps=random_sweeps, seed=seed)
    B, _, _ = k04.metropolis(B, L, T=prep_T, sweeps=prep_sweeps, seed=seed + 17)
    return B


@dataclass
class Trajectory:
    L: int
    rep: int
    activity: int
    attempted: int
    accepted_fraction: float
    final_defect: float
    samples: list[dict[str, float]]


def run_trajectory(
    L: int,
    T: float,
    sweeps: int,
    rep: int,
    random_sweeps: int,
    prep_T: float,
    prep_sweeps: int,
    sample_every: int,
) -> Trajectory:
    rng = random.Random(10_000_019 + 1009 * L + rep)
    B = prepare_state(L, seed=31_337 + 97 * rep + L, random_sweeps=random_sweeps, prep_T=prep_T, prep_sweeps=prep_sweeps)
    sig0 = site_signatures(B, L)
    persistent = set(k04.sites(L))
    N = L**3
    total_steps = sweeps * N
    sample_step = max(1, sample_every * N)
    samples: list[dict[str, float]] = []
    activity = attempted = 0

    for step in range(1, total_steps + 1):
        move = choose_move(B, L, rng)
        if move is None:
            if step % sample_step == 0:
                samples.append(
                    {
                        "sweep": step / N,
                        "persistence": len(persistent) / N,
                        "defect": k04.defect_fraction(B, L),
                        "activity": activity,
                    }
                )
            continue
        attempted += 1
        dE = k04.delta_energy(B, move, L)
        if dE <= 0.0 or rng.random() < math.exp(-dE / T):
            changed = touched_sites(move, L)
            apply_move(B, move)
            activity += 1
            for p in changed:
                if p in persistent and incident_signature(B, p, L) != sig0[p]:
                    persistent.remove(p)
        if step % sample_step == 0:
            samples.append(
                {
                    "sweep": step / N,
                    "persistence": len(persistent) / N,
                    "defect": k04.defect_fraction(B, L),
                    "activity": activity,
                }
            )

    return Trajectory(
        L=L,
        rep=rep,
        activity=activity,
        attempted=attempted,
        accepted_fraction=(activity / attempted if attempted else 0.0),
        final_defect=k04.defect_fraction(B, L),
        samples=samples,
    )


def activity_ld(trajectories: list[Trajectory], s_values: list[float], sweeps: int) -> list[dict[str, float]]:
    Ks = np.array([tr.activity for tr in trajectories], dtype=float)
    out = []
    for s in s_values:
        logs = -s * Ks
        shift = float(np.max(logs))
        weights = np.exp(logs - shift)
        mean_weight = float(np.mean(weights))
        psi = (shift + math.log(mean_weight)) / sweeps
        k_s = float(np.sum(weights * Ks) / np.sum(weights) / sweeps)
        out.append({"s": s, "psi": psi, "activity_density": k_s})
    return out


def persistence_chi4(trajectories: list[Trajectory]) -> list[dict[str, float]]:
    by_t: dict[float, list[float]] = defaultdict(list)
    L = trajectories[0].L
    N = L**3
    for tr in trajectories:
        for row in tr.samples:
            by_t[row["sweep"]].append(row["persistence"])
    out = []
    for t in sorted(by_t):
        vals = np.array(by_t[t], dtype=float)
        out.append(
            {
                "sweep": t,
                "P": float(np.mean(vals)),
                "chi4": float(N * np.var(vals)),
                "n": int(len(vals)),
            }
        )
    return out


@dataclass
class Bootstrap:
    L: int
    rep: int
    accepted_nonuphill: int
    final_defect: float
    final_flippable_density: float
    residual_nonuphill: int | None
    stopped_by: str


def bootstrap_close(
    L: int,
    rep: int,
    random_sweeps: int,
    prep_T: float,
    prep_sweeps: int,
    max_rounds: int,
    exact_residual: bool,
) -> Bootstrap:
    rng = random.Random(50_000_029 + 2003 * L + rep)
    B = prepare_state(L, seed=72_001 + 101 * rep + L, random_sweeps=random_sweeps, prep_T=prep_T, prep_sweeps=prep_sweeps)
    accepted = 0
    stopped_by = "cap"
    for _ in range(max_rounds):
        moves = list(k04.flippable_moves(B, L))
        rng.shuffle(moves)
        chosen = None
        for mv in moves:
            if k04.delta_energy(B, mv, L) <= 0.0:
                chosen = mv
                break
        if chosen is None:
            stopped_by = "jammed"
            break
        apply_move(B, chosen)
        accepted += 1
    residual = None
    if exact_residual or stopped_by == "jammed":
        residual = sum(1 for mv in k04.flippable_moves(B, L) if k04.delta_energy(B, mv, L) <= 0.0)
    return Bootstrap(
        L=L,
        rep=rep,
        accepted_nonuphill=accepted,
        final_defect=k04.defect_fraction(B, L),
        final_flippable_density=len(list(k04.flippable_moves(B, L))) / len(k04.squares(L)),
        residual_nonuphill=residual,
        stopped_by=stopped_by,
    )


def bond_same_count(B: set[tuple[tuple[int, int, int], int]], B0: set[tuple[tuple[int, int, int], int]], universe) -> int:
    return sum(1 for edge in universe if ((edge in B) == (edge in B0)))


@dataclass
class Replica:
    L: int
    epsilon: float
    rep: int
    overlap: float
    accepted_fraction: float
    final_defect: float


def run_coupled_replica(
    L: int,
    T: float,
    epsilon: float,
    sweeps: int,
    rep: int,
    random_sweeps: int,
    prep_T: float,
    prep_sweeps: int,
) -> Replica:
    rng = random.Random(90_000_031 + 3001 * L + 37 * rep + round(1000 * epsilon))
    B0 = prepare_state(L, seed=101_003 + 107 * rep + L, random_sweeps=random_sweeps, prep_T=prep_T, prep_sweeps=prep_sweeps)
    B = k04.randomize(set(B0), L, sweeps=random_sweeps, seed=202_001 + 109 * rep + L)
    universe = all_bonds(L)
    same = bond_same_count(B, B0, universe)
    attempted = accepted = 0
    for _ in range(sweeps * L**3):
        move = choose_move(B, L, rng)
        if move is None:
            continue
        attempted += 1
        dE = k04.delta_energy(B, move, L)
        before = same
        B2 = set(B)
        apply_move(B2, move)
        after = bond_same_count(B2, B0, universe)
        dE_eff = dE - epsilon * (after - before)
        if dE_eff <= 0.0 or rng.random() < math.exp(-dE_eff / T):
            B = B2
            same = after
            accepted += 1
    return Replica(
        L=L,
        epsilon=epsilon,
        rep=rep,
        overlap=same / len(universe),
        accepted_fraction=(accepted / attempted if attempted else 0.0),
        final_defect=k04.defect_fraction(B, L),
    )


def mean_sd(vals: list[float]) -> tuple[float, float]:
    if len(vals) < 2:
        return vals[0], 0.0
    return fmean(vals), pstdev(vals)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ls", default="4")
    ap.add_argument("--reps", type=int, default=6)
    ap.add_argument("--bootstrap-reps", type=int, default=3)
    ap.add_argument("--bootstrap-max-rounds", type=int, default=12)
    ap.add_argument("--exact-bootstrap-residual", action="store_true")
    ap.add_argument("--replica-reps", type=int, default=4)
    ap.add_argument("--T", type=float, default=0.7)
    ap.add_argument("--sweeps", type=int, default=10)
    ap.add_argument("--random-sweeps", type=int, default=10)
    ap.add_argument("--prep-T", type=float, default=0.7)
    ap.add_argument("--prep-sweeps", type=int, default=8)
    ap.add_argument("--sample-every", type=int, default=2)
    ap.add_argument("--s-values", default="-0.08,-0.04,0,0.04,0.08")
    ap.add_argument("--eps-values", default="0,0.5,2.0")
    ap.add_argument("--out")
    args = ap.parse_args()

    Ls = parse_ints(args.Ls)
    s_values = parse_floats(args.s_values)
    eps_values = parse_floats(args.eps_values)

    print("[0] K04 CONSTRAINT-GLASS OBSERVABLE PILOT")
    print(f"    Ls={Ls}, reps={args.reps}, T={args.T}, sweeps={args.sweeps}")
    print("    finite-size smoke test: classify signals; do not claim transitions.")

    results: dict[str, object] = {"parameters": vars(args), "by_L": {}}

    for L in Ls:
        print(f"\n[1] L={L}: activity large deviations, persistence P(t), chi4(t)")
        trajectories = [
            run_trajectory(
                L=L,
                T=args.T,
                sweeps=args.sweeps,
                rep=rep,
                random_sweeps=args.random_sweeps,
                prep_T=args.prep_T,
                prep_sweeps=args.prep_sweeps,
                sample_every=args.sample_every,
            )
            for rep in range(args.reps)
        ]
        acts = [tr.activity / args.sweeps for tr in trajectories]
        afs = [tr.accepted_fraction for tr in trajectories]
        ds = [tr.final_defect for tr in trajectories]
        act_mean, act_sd = mean_sd(acts)
        af_mean, af_sd = mean_sd(afs)
        d_mean, d_sd = mean_sd(ds)
        print(f"    activity density K_t/t = {act_mean:.3f} +/- {act_sd:.3f}")
        print(f"    accepted/attempted     = {af_mean:.3f} +/- {af_sd:.3f}")
        print(f"    final defect fraction  = {d_mean:.3f} +/- {d_sd:.3f}")

        ld = activity_ld(trajectories, s_values, args.sweeps)
        print("    s-ensemble pilot:")
        for row in ld:
            print(f"      s={row['s']:+.3f} psi={row['psi']:+.5f} k_s={row['activity_density']:.3f}")

        pc = persistence_chi4(trajectories)
        chi_peak = max(pc, key=lambda row: row["chi4"])
        print("    persistence/chi4:")
        for row in pc:
            print(f"      t={row['sweep']:5.1f} P={row['P']:.3f} chi4={row['chi4']:.3f}")
        print(f"      peak chi4={chi_peak['chi4']:.3f} at t={chi_peak['sweep']:.1f}")

        print(f"\n[2] L={L}: bootstrap/jamming closure")
        boots = [
            bootstrap_close(
                L=L,
                rep=rep,
                random_sweeps=args.random_sweeps,
                prep_T=args.prep_T,
                prep_sweeps=args.prep_sweeps,
                max_rounds=args.bootstrap_max_rounds,
                exact_residual=args.exact_bootstrap_residual,
            )
            for rep in range(args.reps)
        ]
        b_def, b_def_sd = mean_sd([b.final_defect for b in boots])
        b_flip, b_flip_sd = mean_sd([b.final_flippable_density for b in boots])
        b_k, b_k_sd = mean_sd([float(b.accepted_nonuphill) for b in boots])
        exact_residuals = [b.residual_nonuphill for b in boots if b.residual_nonuphill is not None]
        skipped_residuals = sum(1 for b in boots if b.residual_nonuphill is None)
        stopped = {key: sum(1 for b in boots if b.stopped_by == key) for key in ("jammed", "cap")}
        print(f"    non-uphill flips until jam = {b_k:.2f} +/- {b_k_sd:.2f}")
        print(f"    jammed defect fraction     = {b_def:.3f} +/- {b_def_sd:.3f}")
        print(f"    jammed flippable density   = {b_flip:.3f} +/- {b_flip_sd:.3f}")
        print(f"    stopped by                 = {stopped}")
        if exact_residuals:
            residual = sum(exact_residuals)
            suffix = "" if skipped_residuals == 0 else f"; {skipped_residuals} capped rows skipped"
            print(f"    exact residual non-uphill  = {residual} over {len(exact_residuals)} rows{suffix}")
            assert residual == 0
        else:
            print("    exact residual non-uphill  = skipped (use --exact-bootstrap-residual on deep)")

        print(f"\n[3] L={L}: coupled-replica overlap scan")
        replicas: list[Replica] = []
        for eps in eps_values:
            rows = [
                run_coupled_replica(
                    L=L,
                    T=args.T,
                    epsilon=eps,
                    sweeps=args.sweeps,
                    rep=rep,
                    random_sweeps=args.random_sweeps,
                    prep_T=args.prep_T,
                    prep_sweeps=args.prep_sweeps,
                )
                for rep in range(args.replica_reps)
            ]
            replicas.extend(rows)
            q_mean, q_sd = mean_sd([r.overlap for r in rows])
            acc_mean, acc_sd = mean_sd([r.accepted_fraction for r in rows])
            print(f"    eps={eps:4.2f}: Q={q_mean:.3f} +/- {q_sd:.3f}; acc={acc_mean:.3f} +/- {acc_sd:.3f}")

        q_by_eps = {
            str(eps): [asdict(r) for r in replicas if abs(r.epsilon - eps) < 1e-12]
            for eps in eps_values
        }
        results["by_L"][str(L)] = {
            "trajectories": [asdict(tr) for tr in trajectories],
            "activity_ld": ld,
            "persistence_chi4": pc,
            "bootstrap": [asdict(b) for b in boots],
            "coupled_replica": q_by_eps,
        }

    print(
        """
[verdict]
This pilot turns the K04 constraint-glass programme into measured observables:

  * activity K_t and an s-ensemble estimator for active/inactive bias;
  * persistence P(t) and chi4(t) for dynamic heterogeneity;
  * bootstrap closure to a zero-temperature jammed core;
  * coupled-replica overlap Q(epsilon) as the RFOT gate.

Current small-L output is diagnostic only.  A publishable condensed-matter
claim needs scaling in L, longer times, and preferably trajectory sampling for
the s-ensemble.  The computational object is now well-posed.
exit 0"""
    )

    if args.out:
        Path(args.out).write_text(json.dumps(results, indent=2, sort_keys=True))
        print(f"[wrote] {args.out}")


if __name__ == "__main__":
    main()
