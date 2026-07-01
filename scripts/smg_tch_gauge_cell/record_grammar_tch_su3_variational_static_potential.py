#!/usr/bin/env python3
r"""Variational TCH SU(3) static-potential test on the bond-bipyramid bulk.

This is the next method-development rung after
``record_grammar_tch_su3_static_potential_pipeline.py``.

The straight-path TCH Wilson source was negative on the current block even
after the ordinary cubic SU(3) Cornell extractor was validated.  This script
keeps the same Cornell/static-potential readout but changes the source operator:
for each source separation it builds a geometry-native variational basis from
actual bond-graph paths with lengths R, R+2, R+4 when available.  The principal
correlator of that path basis is then used to extract V_eff(R,T).

Boundary
--------
This is still a diagnostic calculation.  A positive result would license a
larger/deeper production run; a negative result says the current finite block
and simple path basis still do not define a stable TCH static potential.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import deque
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

import numpy as np

from record_grammar_cubic_su3_static_potential_pipeline import (
    PotentialPoint,
    fit_cornell,
    fit_linear_tail,
    project_su3,
)
from record_grammar_tch_su3_bulk_polyakov_torus import I3, dag, staple_terms, sweep
from record_grammar_tch_su3_static_potential_pipeline import (
    ape_smear_spatial_irregular,
    spatial_link,
    spatial_staple_terms,
    straight_source_paths,
    su3_unitarity_error,
    thermal_lattice,
    time_link,
)


DEFAULT_OUTPUT = "python_code/tch_su3_variational_static_potential_larger.json"


def parse_size(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace(",", "x").split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("size must be NxMxK, e.g. 5x5x3")
    values = tuple(int(part) for part in parts)
    if any(value < 2 for value in values):
        raise argparse.ArgumentTypeError("size components must be >= 2")
    return values  # type: ignore[return-value]


def parse_ints(text: str) -> tuple[int, ...]:
    values = tuple(int(item) for item in text.replace(",", " ").split() if item.strip())
    if not values:
        raise argparse.ArgumentTypeError("at least one integer is required")
    return tuple(sorted(dict.fromkeys(values)))


def build_adjacency(lattice: dict) -> list[list[tuple[int, int, bool]]]:
    adj: list[list[tuple[int, int, bool]]] = [[] for _ in lattice["vertices"]]
    for edge_i, (left, right) in enumerate(lattice["edges"]):
        li = lattice["vidx"][left]
        ri = lattice["vidx"][right]
        adj[li].append((ri, edge_i, True))
        adj[ri].append((li, edge_i, False))
    for row in adj:
        row.sort(key=lambda item: (item[0], item[1], item[2]))
    return adj


def evenly_sample(items: list, limit: int) -> list:
    if limit <= 0 or len(items) <= limit:
        return list(items)
    if limit == 1:
        return [items[len(items) // 2]]
    return [items[round(i * (len(items) - 1) / (limit - 1))] for i in range(limit)]


def enumerate_path_groups(
    adj: list[list[tuple[int, int, bool]]],
    start: int,
    end: int,
    base_len: int,
    offsets: tuple[int, ...],
    paths_per_group: int,
    visit_cap: int,
) -> dict[int, list[list[tuple[int, bool]]]]:
    """Enumerate graph paths by length offset from the straight separation."""

    wanted = {base_len + offset: offset for offset in offsets}
    max_len = max(wanted)
    groups: dict[int, list[list[tuple[int, bool]]]] = {offset: [] for offset in offsets}
    queue = deque([(start, frozenset({start}), [])])
    visits = 0
    while queue and visits < visit_cap:
        node, seen, path = queue.popleft()
        visits += 1
        length = len(path)
        if length > max_len:
            continue
        if node == end and length in wanted:
            offset = wanted[length]
            if len(groups[offset]) < paths_per_group:
                groups[offset].append(path)
            if all(len(groups[offset]) >= paths_per_group for offset in offsets):
                break
            continue
        if length == max_len:
            continue
        for nxt, edge_i, forward in adj[node]:
            if nxt in seen:
                continue
            queue.append((nxt, seen | {nxt}, path + [(edge_i, forward)]))
    return groups


def build_basis(
    lattice: dict,
    max_r: int,
    offsets: tuple[int, ...],
    paths_per_group: int,
    max_pairs_per_r: int,
    visit_cap: int,
) -> dict[int, dict]:
    adj = build_adjacency(lattice)
    straight = straight_source_paths(lattice, max_r)
    basis: dict[int, dict] = {}
    for r in range(1, max_r + 1):
        items = evenly_sample(straight.get(r, []), max_pairs_per_r)
        raw_pairs = []
        offset_counts = {offset: 0 for offset in offsets}
        for start, end, _straight_path in items:
            groups = enumerate_path_groups(adj, start, end, r, offsets, paths_per_group, visit_cap)
            for offset in offsets:
                if groups[offset]:
                    offset_counts[offset] += 1
            raw_pairs.append({"start": start, "end": end, "groups": groups})
        # Use an R-specific basis.  The straight source is always required;
        # detour offsets are included only when they exist for at least two
        # endpoint pairs.  This prevents the small-block basis from vanishing at
        # larger R while still letting real detours enter where the geometry
        # supports them.
        active_offsets = tuple(offset for offset in offsets if offset == 0 or offset_counts[offset] >= 2)
        pairs = [
            item
            for item in raw_pairs
            if all(item["groups"].get(offset) for offset in active_offsets)
        ]
        basis[r] = {"offsets": active_offsets, "pairs": pairs, "offset_counts": offset_counts}
    return basis


def path_transport(U: np.ndarray, lattice: dict, path: list[tuple[int, bool]], t: int) -> np.ndarray:
    M = I3.copy()
    for edge_i, forward in path:
        link = U[spatial_link(lattice, edge_i, t)]
        M = M @ (link if forward else dag(link))
    return M


def group_transport(
    U: np.ndarray,
    lattice: dict,
    paths: list[list[tuple[int, bool]]],
    t: int,
) -> np.ndarray:
    mats = [path_transport(U, lattice, path, t) for path in paths]
    return sum(mats, np.zeros((3, 3), dtype=complex)) / len(mats)


def temporal_transport(U: np.ndarray, lattice: dict, vertex_i: int, t0: int, extent: int) -> np.ndarray:
    M = I3.copy()
    for dt in range(extent):
        M = M @ U[time_link(lattice, vertex_i, t0 + dt)]
    return M


def variational_correlator_table(
    U: np.ndarray,
    lattice: dict,
    basis: dict[int, list[dict]],
    offsets: tuple[int, ...],
    max_r: int,
    max_t: int,
) -> dict[int, dict[int, np.ndarray]]:
    out: dict[int, dict[int, np.ndarray]] = {}
    for r in range(1, max_r + 1):
        entry = basis.get(r, {"offsets": (), "pairs": []})
        active_offsets = tuple(entry["offsets"])
        dim = len(active_offsets)
        pairs = entry["pairs"]
        if not pairs:
            continue
        out[r] = {}
        for extent_t in range(1, max_t + 1):
            C = np.zeros((dim, dim), dtype=float)
            count = 0
            for t0 in range(lattice["Nt"]):
                top_t = t0 + extent_t
                for item in pairs:
                    start = int(item["start"])
                    end = int(item["end"])
                    tend = temporal_transport(U, lattice, end, t0, extent_t)
                    tstart = temporal_transport(U, lattice, start, t0, extent_t)
                    bottom = [
                        group_transport(U, lattice, item["groups"][offset], t0)
                        for offset in active_offsets
                    ]
                    top = [
                        group_transport(U, lattice, item["groups"][offset], top_t)
                        for offset in active_offsets
                    ]
                    for i in range(dim):
                        for j in range(dim):
                            loop = bottom[i] @ tend @ dag(top[j]) @ dag(tstart)
                            C[i, j] += float(np.real(np.trace(loop) / 3.0))
                    count += 1
            out[r][extent_t] = 0.5 * (C / count + (C / count).T)
    return out


def principal_correlators(
    table: dict[int, np.ndarray],
    t0: int,
    ridge: float,
) -> dict[int, float | None]:
    if t0 not in table:
        return {}
    C0 = 0.5 * (table[t0] + table[t0].T)
    evals, evecs = np.linalg.eigh(C0)
    keep = evals > max(ridge, ridge * abs(float(evals[-1])) if evals.size else ridge)
    if not np.any(keep):
        return {t: None for t in table}
    inv_sqrt = evecs[:, keep] @ np.diag(1.0 / np.sqrt(evals[keep])) @ evecs[:, keep].T
    out: dict[int, float | None] = {}
    for t, C in table.items():
        X = inv_sqrt @ (0.5 * (C + C.T)) @ inv_sqrt
        vals = np.linalg.eigvalsh(0.5 * (X + X.T))
        lam = float(vals[-1])
        out[t] = lam if math.isfinite(lam) and lam > 0.0 else None
    return out


def potential_points_from_samples(
    samples: list[dict[int, dict[int, np.ndarray]]],
    max_r: int,
    max_t: int,
    t0: int,
    t_min: int,
    ridge: float,
) -> list[PotentialPoint]:
    points: list[PotentialPoint] = []
    by_sample = {
        r: [principal_correlators(sample.get(r, {}), t0, ridge) for sample in samples]
        for r in range(1, max_r + 1)
    }
    for r in range(1, max_r + 1):
        estimates: list[float] = []
        used_t: list[int] = []
        for t in range(t_min, max_t):
            vals = []
            for pcs in by_sample[r]:
                left = pcs.get(t)
                right = pcs.get(t + 1)
                if left and right and left > 0 and right > 0:
                    vals.append(math.log(left / right))
            if vals:
                estimates.append(float(np.mean(vals)))
                used_t.append(t)
        if not estimates:
            points.append(PotentialPoint(r, None, None, 0, []))
            continue
        value = float(np.mean(estimates))
        se = float(np.std(estimates, ddof=1) / math.sqrt(len(estimates))) if len(estimates) > 1 else None
        n = sum(1 for pcs in by_sample[r] for t in used_t if pcs.get(t) and pcs.get(t + 1))
        points.append(PotentialPoint(r, value, se, n, used_t))
    return points


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{digits}f}"


def run(args: argparse.Namespace) -> dict:
    lattice = thermal_lattice(*args.size, args.time_extent)
    all_terms = staple_terms(lattice)
    smear_terms = spatial_staple_terms(lattice)
    basis = build_basis(
        lattice,
        args.max_r,
        args.basis_offsets,
        args.paths_per_group,
        args.max_pairs_per_r,
        args.path_visit_cap,
    )
    rng = np.random.default_rng(args.seed)
    U = np.broadcast_to(I3, (lattice["Nlinks"], 3, 3)).copy()
    eps = args.eps0
    acc_values: list[float] = []

    for _ in range(args.therm):
        U, acc = sweep(U, args.beta, rng, all_terms, lattice["Nlinks"], eps, hits=args.hits)
        eps = min(max(eps * (1.0 + 0.1 * (acc - 0.5)), 0.04), 1.0)
        acc_values.append(float(acc))

    samples: list[dict[int, dict[int, np.ndarray]]] = []
    for _ in range(args.meas):
        for _ in range(args.skip):
            U, acc = sweep(U, args.beta, rng, all_terms, lattice["Nlinks"], eps, hits=args.hits)
            eps = min(max(eps * (1.0 + 0.05 * (acc - 0.5)), 0.04), 1.0)
            acc_values.append(float(acc))
        Us = ape_smear_spatial_irregular(U, lattice, smear_terms, args.smear_steps, args.ape_alpha)
        samples.append(variational_correlator_table(Us, lattice, basis, args.basis_offsets, args.max_r, args.max_t))

    points = potential_points_from_samples(samples, args.max_r, args.max_t, args.gevp_t0, args.t_min, args.ridge)
    fit_r_max = args.fit_r_max if args.fit_r_max else args.max_r
    cornell = fit_cornell(points, args.smear_steps, args.nt_reference, args.fit_r_min, fit_r_max)
    linear = fit_linear_tail(points, args.smear_steps, args.nt_reference, args.fit_r_min, fit_r_max)
    verdict = "VARIATIONAL-POSITIVE" if cornell.sigma and cornell.sigma > 0 else "VARIATIONAL-NOT-STABLE"
    if cornell.sigma is None:
        verdict = "VARIATIONAL-NO-FIT"

    return {
        "profile": "bond-bipyramid-su3-variational-static-potential",
        "beta": args.beta,
        "nt_reference": args.nt_reference,
        "size": args.size,
        "time_extent": args.time_extent,
        "nv": int(lattice["Nv"]),
        "ne": int(lattice["Ne"]),
        "nlinks": int(lattice["Nlinks"]),
        "nplaquettes": int(lattice["Np"]),
        "basis_offsets": args.basis_offsets,
        "paths_per_group": args.paths_per_group,
        "max_pairs_per_r": args.max_pairs_per_r,
        "basis_pair_counts": {str(r): len(basis.get(r, {}).get("pairs", [])) for r in range(1, args.max_r + 1)},
        "active_offsets_by_r": {str(r): list(basis.get(r, {}).get("offsets", ())) for r in range(1, args.max_r + 1)},
        "offset_counts_by_r": {str(r): basis.get(r, {}).get("offset_counts", {}) for r in range(1, args.max_r + 1)},
        "therm": args.therm,
        "meas": args.meas,
        "skip": args.skip,
        "hits": args.hits,
        "seed": args.seed,
        "ape_alpha": args.ape_alpha,
        "smear_steps": args.smear_steps,
        "max_r": args.max_r,
        "max_t": args.max_t,
        "t_min": args.t_min,
        "gevp_t0": args.gevp_t0,
        "fit_r_min": args.fit_r_min,
        "fit_r_max": fit_r_max,
        "acceptance_mean": float(np.mean(acc_values)),
        "unitarity_error": su3_unitarity_error(U),
        "potential_points": [asdict(point) for point in points],
        "fit": asdict(cornell),
        "linear_tail_fit": asdict(linear),
        "verdict": verdict,
        "interpretation": (
            "Geometry-native variational path basis applied to the bond-bipyramid TCH static potential. "
            "Positive sigma would license a deeper production run; negative/no-fit keeps the TCH scale gate open."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beta", type=float, default=5.5)
    parser.add_argument("--nt-reference", type=int, default=4)
    parser.add_argument("--size", type=parse_size, default=parse_size("5x5x3"))
    parser.add_argument("--time-extent", type=int, default=8)
    parser.add_argument("--therm", type=int, default=80)
    parser.add_argument("--meas", type=int, default=80)
    parser.add_argument("--skip", type=int, default=2)
    parser.add_argument("--hits", type=int, default=1)
    parser.add_argument("--eps0", type=float, default=0.30)
    parser.add_argument("--ape-alpha", type=float, default=0.45)
    parser.add_argument("--smear-steps", type=int, default=4)
    parser.add_argument("--max-r", type=int, default=4)
    parser.add_argument("--max-t", type=int, default=6)
    parser.add_argument("--t-min", type=int, default=1)
    parser.add_argument("--gevp-t0", type=int, default=1)
    parser.add_argument("--ridge", type=float, default=1.0e-8)
    parser.add_argument("--fit-r-min", type=int, default=1)
    parser.add_argument("--fit-r-max", type=int, default=4)
    parser.add_argument("--basis-offsets", type=parse_ints, default=parse_ints("0,2,4"))
    parser.add_argument("--paths-per-group", type=int, default=3)
    parser.add_argument("--max-pairs-per-r", type=int, default=80)
    parser.add_argument("--path-visit-cap", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260702)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.max_t >= args.time_extent:
        raise SystemExit("--max-t must be smaller than --time-extent")
    if args.gevp_t0 < 1 or args.gevp_t0 >= args.max_t:
        raise SystemExit("--gevp-t0 must satisfy 1 <= t0 < max_t")

    t0 = time.time()
    result = run(args)
    result["elapsed_s"] = time.time() - t0
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("Bond-bipyramid SU(3) variational static-potential pipeline")
    print("=" * 96)
    print(f"  beta={args.beta:.3f}; Nt reference={args.nt_reference}; size={args.size}; time extent={args.time_extent}")
    print(f"  therm={args.therm}; meas={args.meas}; skip={args.skip}; max R,T={args.max_r},{args.max_t}")
    print(f"  basis offsets={args.basis_offsets}; paths/group={args.paths_per_group}; max pairs/R={args.max_pairs_per_r}")
    print(f"  basis pair counts: {result['basis_pair_counts']}")
    print(f"  acceptance={result['acceptance_mean']:.3f}; unitarity_error={result['unitarity_error']:.2e}")
    print("\n[variational potential points]")
    for point in result["potential_points"]:
        print(
            f"  R={point['r']}: V={fmt(point['v'])} se={fmt(point['se'])} "
            f"n={point['n']} T={point['t_values']}"
        )
    fit = result["fit"]
    linear = result["linear_tail_fit"]
    print("\n[Cornell fit]")
    print(
        f"  sigma={fmt(fit['sigma'])} se={fmt(fit['sigma_se'])} "
        f"a√sigma={fmt(fit['a_sqrt_sigma'])} Tc/√sigma={fmt(fit['tc_over_sqrt_sigma'])} "
        f"alpha={fmt(fit['alpha'])} V0={fmt(fit['v0'])} chi2/dof={fmt(fit['chi2_dof'])} "
        f"verdict={fit['verdict']}"
    )
    print("[linear tail control]")
    print(
        f"  sigma={fmt(linear['sigma'])} se={fmt(linear['sigma_se'])} "
        f"a√sigma={fmt(linear['a_sqrt_sigma'])} Tc/√sigma={fmt(linear['tc_over_sqrt_sigma'])} "
        f"V0={fmt(linear['v0'])} chi2/dof={fmt(linear['chi2_dof'])} "
        f"verdict={linear['verdict']}"
    )
    print(f"\nVERDICT: {result['verdict']}")
    print("  " + result["interpretation"])
    print(f"  wrote {args.output}")
    print("exit 0")
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
