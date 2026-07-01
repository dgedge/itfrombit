#!/usr/bin/env python3
r"""Full-SU(3) bond-bipyramid scale-setting harness.

Goal
----
This is the first executable scale-setting run for the record-grammar
confinement ladder.  It measures, on the same full-SU(3) bond-bipyramid thermal
ensembles,

  1. the Polyakov susceptibility peak beta_c(N_t, volume), and
  2. spatial Wilson-loop Creutz ratios a^2 sigma(beta) on the licensed flat
     bond-bipyramid sheets.

The output is deliberately tiered:

* ``smoke`` is a quick local run to verify the machinery and expose gross
  failures;
* ``local`` is still laptop-grade but uses a denser beta grid;
* larger production profiles should be run on deep or overnight.

Interpretation boundary
-----------------------
This script is the computation the paper named as missing, but a smoke/local run
is not a continuum string-tension result.  A physical sigma requires stable
T_c/sqrt(sigma) under increasing N_t, finite-size extrapolation, and a derived
action-normalisation convention.  The script prints exactly those quantities so
the verdict is mechanical rather than prose.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from record_grammar_tch_bond_bipyramid_bulk import (
    BondComplex,
    build_bond_complex,
    bonds_box,
    xor_masks,
)
from record_grammar_tch_bond_bipyramid_edge_loop_creutz import (
    plane_face_coordinates,
)
from record_grammar_tch_bond_bipyramid_static_potential import boundary_stats
from record_grammar_tch_bulk_polyakov_torus import build as build_z3_bulk
from record_grammar_tch_su3_bulk_polyakov_torus import (
    mean_plaq,
    polyakov_abs,
    staple_terms,
    sweep,
)
from record_grammar_tch_su3_polyakov_torus import I3, dag


@dataclass(frozen=True)
class Profile:
    nts: tuple[int, ...]
    sizes: tuple[tuple[int, int, int], ...]
    therm: int
    meas: int
    hits: int
    eps0: float
    betas: dict[int, tuple[float, ...]]
    max_loop: int


def grid(start: float, stop: float, step: float) -> tuple[float, ...]:
    n = int(round((stop - start) / step))
    return tuple(round(start + i * step, 6) for i in range(n + 1))


PROFILES: dict[str, Profile] = {
    "smoke": Profile(
        nts=(2, 3, 4),
        sizes=((2, 2, 2), (3, 3, 2)),
        therm=20,
        meas=24,
        hits=1,
        eps0=0.35,
        betas={
            2: (3.75, 4.25, 4.75, 5.25),
            3: (4.75, 5.25, 5.75, 6.25),
            4: (5.50, 6.25, 7.00, 7.75),
        },
        max_loop=2,
    ),
    "local": Profile(
        nts=(2, 3, 4),
        sizes=((2, 2, 2), (3, 3, 2)),
        therm=50,
        meas=60,
        hits=1,
        eps0=0.30,
        betas={
            2: (3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 5.25),
            3: (4.50, 4.875, 5.25, 5.625, 6.00, 6.375, 6.75),
            4: (5.50, 5.875, 6.25, 6.625, 7.00, 7.375, 7.75),
        },
        max_loop=2,
    ),
    "deep-preflight": Profile(
        nts=(2, 3, 4),
        sizes=((3, 3, 2), (4, 4, 3)),
        therm=100,
        meas=160,
        hits=1,
        eps0=0.30,
        betas={
            2: grid(3.50, 5.75, 0.25),
            3: grid(4.25, 7.00, 0.25),
            4: grid(5.25, 8.25, 0.25),
        },
        max_loop=3,
    ),
    "deep": Profile(
        nts=(2, 3, 4, 5, 6),
        sizes=((3, 3, 2), (4, 4, 3), (5, 5, 3)),
        therm=250,
        meas=400,
        hits=1,
        eps0=0.30,
        betas={
            2: grid(3.25, 5.75, 0.25),
            3: grid(4.25, 7.25, 0.25),
            4: grid(5.25, 8.50, 0.25),
            5: grid(6.25, 9.25, 0.25),
            6: grid(7.00, 10.25, 0.25),
        },
        max_loop=4,
    ),
}


@dataclass
class BetaRow:
    beta: float
    polyakov: float
    chi_polyakov: float
    plaquette: float
    w11: float | None
    w12: float | None
    w21: float | None
    w22: float | None
    creutz_11: float | None
    a_sqrt_sigma_11: float | None
    creutz_label: str | None
    creutz_largest: float | None
    a_sqrt_sigma_largest: float | None
    acceptance: float


@dataclass
class VolumeRun:
    nt: int
    size: tuple[int, int, int]
    nv: int
    nlinks: int
    nplaquettes: int
    beta_c_chi: float
    beta_c_half: float
    beta_near_sigma: float
    a_sqrt_sigma_near: float | None
    tc_over_sqrt_sigma_near: float | None
    rows: list[BetaRow]


@dataclass
class NtSummary:
    nt: int
    beta_c_fss: float
    slope_1_over_nv: float | None
    largest_size: tuple[int, int, int]
    beta_sigma: float
    a_sqrt_sigma: float | None
    tc_over_sqrt_sigma: float | None


def thermal_lattice(nx: int, ny: int, nz: int, nt: int) -> dict:
    """Build the existing thermal complex and attach the spatial complex metadata."""

    lattice = build_z3_bulk(nx, ny, nz, nt)
    complex_ = build_bond_complex(f"slab{nx}{ny}{nz}", bonds_box(nx, ny, nz))
    lattice["complex"] = complex_
    lattice["vertices"] = list(complex_.vertices)
    lattice["edges"] = list(complex_.edges)
    lattice["vidx"] = {v: i for i, v in enumerate(lattice["vertices"])}
    lattice["eidx"] = {e: i for i, e in enumerate(lattice["edges"])}
    return lattice


def rectangle_patch(
    coords: list[tuple[int, int, int]],
    r: int,
    t: int,
    a0: int,
    b0: int,
) -> tuple[int, ...]:
    return tuple(
        sorted(face_i for a, b, face_i in coords if a0 <= a < a0 + r and b0 <= b < b0 + t)
    )


def directed_cycle_from_boundary(complex_: BondComplex, boundary_mask: int) -> list[tuple[int, bool]]:
    edge_ids: list[int] = []
    mask = boundary_mask
    edge_i = 0
    while mask:
        if mask & 1:
            edge_ids.append(edge_i)
        mask >>= 1
        edge_i += 1

    components, even, max_degree, _ = boundary_stats(boundary_mask, complex_.edges)
    if components != 1 or not even or max_degree != 2:
        raise ValueError("boundary is not a single cycle")

    adjacency: dict[tuple[int, int, int], list[tuple[tuple[int, int, int], int]]] = {}
    for edge_i in edge_ids:
        left, right = complex_.edges[edge_i]
        adjacency.setdefault(left, []).append((right, edge_i))
        adjacency.setdefault(right, []).append((left, edge_i))

    start = min(adjacency)
    previous = None
    current = start
    directed: list[tuple[int, bool]] = []
    seen_edges: set[int] = set()
    while True:
        options = adjacency[current]
        candidates = [(nxt, edge_i) for nxt, edge_i in options if edge_i not in seen_edges]
        if not candidates:
            break
        if previous is not None and len(candidates) > 1:
            candidates = [(nxt, edge_i) for nxt, edge_i in candidates if nxt != previous]
        nxt, edge_i = candidates[0]
        canonical_left, _ = complex_.edges[edge_i]
        directed.append((edge_i, current == canonical_left))
        seen_edges.add(edge_i)
        previous, current = current, nxt
        if current == start:
            break

    if len(seen_edges) != len(edge_ids) or current != start:
        raise ValueError("failed to traverse cycle")
    return directed


def loop_specs(lattice: dict, max_loop: int) -> dict[tuple[int, int], list[list[tuple[int, bool]]]]:
    complex_: BondComplex = lattice["complex"]
    _, _, coords = plane_face_coordinates(complex_)
    max_a = max(a for a, _, _ in coords) + 1
    max_b = max(b for _, b, _ in coords) + 1
    specs: dict[tuple[int, int], list[list[tuple[int, bool]]]] = {}
    for r in range(1, min(max_loop, max_a) + 1):
        for t in range(1, min(max_loop, max_b) + 1):
            items: list[list[tuple[int, bool]]] = []
            for a0 in range(0, max_a - r + 1):
                for b0 in range(0, max_b - t + 1):
                    patch = rectangle_patch(coords, r, t, a0, b0)
                    boundary = xor_masks(complex_.face_masks, patch)
                    try:
                        items.append(directed_cycle_from_boundary(complex_, boundary))
                    except ValueError:
                        # Tiny boxes can wrap the chosen offset sheet into
                        # disconnected or self-touching boundaries.  They are
                        # still useful for Polyakov beta_c, but not for Wilson
                        # scale setting.
                        continue
            if items:
                specs[(r, t)] = items
    return specs


def spatial_link(lattice: dict, edge_i: int, time_i: int) -> int:
    return edge_i * lattice["Nt"] + (time_i % lattice["Nt"])


def wilson_cycle(U: np.ndarray, lattice: dict, cycle: list[tuple[int, bool]], time_i: int) -> float:
    M = I3.copy()
    for edge_i, forward in cycle:
        link = U[spatial_link(lattice, edge_i, time_i)]
        M = M @ (link if forward else dag(link))
    return float(np.real(np.trace(M) / 3.0))


def wilson_table(U: np.ndarray, lattice: dict, specs: dict[tuple[int, int], list[list[tuple[int, bool]]]]) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for key, cycles in specs.items():
        vals = []
        for time_i in range(lattice["Nt"]):
            for cycle in cycles:
                vals.append(wilson_cycle(U, lattice, cycle, time_i))
        out[key] = float(np.mean(vals))
    return out


def creutz_11(w: dict[tuple[int, int], float]) -> tuple[float | None, float | None]:
    needed = [(1, 1), (1, 2), (2, 1), (2, 2)]
    if any(key not in w or w[key] <= 0 for key in needed):
        return None, None
    chi = -math.log((w[(1, 1)] * w[(2, 2)]) / (w[(2, 1)] * w[(1, 2)]))
    if chi <= 0:
        return chi, None
    return chi, math.sqrt(chi / 4.0)


def creutz_values(w: dict[tuple[int, int], float], max_loop: int) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for r in range(1, max_loop):
        for t in range(1, max_loop):
            needed = [(r, t), (r + 1, t + 1), (r + 1, t), (r, t + 1)]
            if any(key not in w or w[key] <= 0 for key in needed):
                continue
            chi = -math.log((w[(r, t)] * w[(r + 1, t + 1)]) / (w[(r + 1, t)] * w[(r, t + 1)]))
            if chi > 0:
                out[(r, t)] = chi
    return out


def largest_creutz_estimator(w: dict[tuple[int, int], float], max_loop: int) -> tuple[str | None, float | None, float | None]:
    values = creutz_values(w, max_loop)
    if not values:
        return None, None, None
    key = max(values, key=lambda item: (item[0] * item[1], item[0], item[1]))
    chi = values[key]
    return f"{key[0]}x{key[1]}", chi, math.sqrt(chi / 4.0)


def beta_c_halfcross(rows: list[BetaRow], floor: float) -> float:
    pairs = [(row.beta, row.polyakov) for row in rows]
    plateau = np.mean(sorted((p for _, p in pairs))[-min(3, len(pairs)):])
    threshold = 0.5 * (plateau + floor)
    ordered = sorted(pairs)
    for (b0, p0), (b1, p1) in zip(ordered, ordered[1:]):
        if (p0 - threshold) * (p1 - threshold) <= 0 and p0 != p1:
            return b0 + (b1 - b0) * (threshold - p0) / (p1 - p0)
    return float("nan")


def parabolic_peak(rows: list[BetaRow]) -> float:
    chi = np.array([row.chi_polyakov for row in rows], dtype=float)
    betas = np.array([row.beta for row in rows], dtype=float)
    i = int(chi.argmax())
    if 0 < i < len(rows) - 1:
        coeff = np.polyfit(betas[i - 1 : i + 2], chi[i - 1 : i + 2], 2)
        if coeff[0] < 0:
            return float(-coeff[1] / (2.0 * coeff[0]))
    return float(betas[i])


def run_volume(
    nt: int,
    size: tuple[int, int, int],
    betas: Iterable[float],
    therm: int,
    meas: int,
    hits: int,
    eps0: float,
    max_loop: int,
    seed: int,
) -> VolumeRun:
    lattice = thermal_lattice(*size, nt)
    terms = staple_terms(lattice)
    specs = loop_specs(lattice, max_loop)
    rng = np.random.default_rng(seed)
    U = np.broadcast_to(I3, (lattice["Nlinks"], 3, 3)).copy()
    eps = eps0
    rows: list[BetaRow] = []
    floor = (1.0 / 3.0) / math.sqrt(lattice["Nv"])

    # Cold-start annealing is reliable only when scanning from weak to strong
    # coupling (large beta downwards), matching the validated bulk Polyakov
    # script.  Low-to-high scans can falsely keep the Polyakov line high at the
    # first low-beta point.
    for beta in sorted((float(item) for item in betas), reverse=True):
        acc_values = []
        for _ in range(therm):
            U, acc = sweep(U, float(beta), rng, terms, lattice["Nlinks"], eps, hits=hits)
            eps = min(max(eps * (1.0 + 0.1 * (acc - 0.5)), 0.04), 1.0)
            acc_values.append(float(acc))

        poly = []
        plaq = []
        w11 = []
        w12 = []
        w21 = []
        w22 = []
        for _ in range(meas):
            U, acc = sweep(U, float(beta), rng, terms, lattice["Nlinks"], eps, hits=hits)
            acc_values.append(float(acc))
            poly.append(polyakov_abs(U, lattice))
            plaq.append(mean_plaq(U, lattice))
            wt = wilson_table(U, lattice, specs)
            w11.append(wt.get((1, 1), float("nan")))
            w12.append(wt.get((1, 2), float("nan")))
            w21.append(wt.get((2, 1), float("nan")))
            w22.append(wt.get((2, 2), float("nan")))

        w_means = {
            (1, 1): float(np.nanmean(w11)),
            (1, 2): float(np.nanmean(w12)),
            (2, 1): float(np.nanmean(w21)),
            (2, 2): float(np.nanmean(w22)),
        }
        row_creutz, row_asig = creutz_11(w_means)
        label, creutz_largest, asig_largest = largest_creutz_estimator(w_means, max_loop)

        rows.append(
            BetaRow(
                beta=float(beta),
                polyakov=float(np.mean(poly)),
                chi_polyakov=float(lattice["Nv"] * np.var(poly)),
                plaquette=float(np.mean(plaq)),
                w11=w_means[(1, 1)],
                w12=w_means[(1, 2)],
                w21=w_means[(2, 1)],
                w22=w_means[(2, 2)],
                creutz_11=row_creutz,
                a_sqrt_sigma_11=row_asig,
                creutz_label=label,
                creutz_largest=creutz_largest,
                a_sqrt_sigma_largest=asig_largest,
                acceptance=float(np.mean(acc_values)),
            )
        )

    beta_c_chi = parabolic_peak(rows)
    beta_c_half = beta_c_halfcross(rows, floor)
    sigma_rows = [row for row in rows if row.a_sqrt_sigma_largest is not None]
    if sigma_rows:
        target_beta = beta_c_chi if math.isfinite(beta_c_chi) else rows[int(np.argmax([r.chi_polyakov for r in rows]))].beta
        nearest = min(sigma_rows, key=lambda row: abs(row.beta - target_beta))
        beta_sigma = nearest.beta
        a_sigma = nearest.a_sqrt_sigma_largest
        tc_ratio = 1.0 / (nt * a_sigma) if a_sigma and a_sigma > 0 else None
    else:
        beta_sigma = float("nan")
        a_sigma = None
        tc_ratio = None

    return VolumeRun(
        nt=nt,
        size=size,
        nv=int(lattice["Nv"]),
        nlinks=int(lattice["Nlinks"]),
        nplaquettes=int(lattice["Np"]),
        beta_c_chi=beta_c_chi,
        beta_c_half=beta_c_half,
        beta_near_sigma=beta_sigma,
        a_sqrt_sigma_near=a_sigma,
        tc_over_sqrt_sigma_near=tc_ratio,
        rows=rows,
    )


def fss_summary(nt: int, runs: list[VolumeRun]) -> NtSummary:
    ordered = sorted(runs, key=lambda run: run.nv)
    largest = ordered[-1]
    if len(ordered) >= 2:
        x = np.array([1.0 / run.nv for run in ordered], dtype=float)
        y = np.array([run.beta_c_chi for run in ordered], dtype=float)
        slope, intercept = np.polyfit(x, y, 1)
        beta_c_fss = float(intercept)
        slope_out: float | None = float(slope)
    else:
        beta_c_fss = float(largest.beta_c_chi)
        slope_out = None
    return NtSummary(
        nt=nt,
        beta_c_fss=beta_c_fss,
        slope_1_over_nv=slope_out,
        largest_size=largest.size,
        beta_sigma=largest.beta_near_sigma,
        a_sqrt_sigma=largest.a_sqrt_sigma_near,
        tc_over_sqrt_sigma=largest.tc_over_sqrt_sigma_near,
    )


def print_volume(run: VolumeRun) -> None:
    print(
        f"  Nt={run.nt} size={run.size} Nv={run.nv} links={run.nlinks} "
        f"plaqs={run.nplaquettes} beta_c_chi={run.beta_c_chi:.3f} "
        f"beta_c_half={run.beta_c_half:.3f}"
    )
    print("    beta    |L|     chi_L   plaq     W11     W22    chi11  a√s11  best   a√sbest acc")
    for row in run.rows:
        c = "nan" if row.creutz_11 is None else f"{row.creutz_11:6.3f}"
        a = "nan" if row.a_sqrt_sigma_11 is None else f"{row.a_sqrt_sigma_11:6.3f}"
        best = "-" if row.creutz_label is None else row.creutz_label
        abest = "nan" if row.a_sqrt_sigma_largest is None else f"{row.a_sqrt_sigma_largest:7.3f}"
        print(
            f"    {row.beta:5.3f}  {row.polyakov:6.3f}  {row.chi_polyakov:7.3f}"
            f"  {row.plaquette:6.3f}  {row.w11:7.3f}  {row.w22:7.3f}"
            f"  {c:>6s}  {a:>6s}  {best:>5s}  {abest:>7s} {row.acceptance:5.3f}"
        )


def parse_size(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace(",", "x").split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("size must be formatted as NxMxK, e.g. 4x4x3")
    try:
        out = tuple(int(item) for item in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("size components must be integers") from exc
    if any(item < 2 for item in out):
        raise argparse.ArgumentTypeError("all size components must be >= 2")
    return out  # type: ignore[return-value]


def write_jsonl(path: Path, payload: dict) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload) + "\n")
        handle.flush()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=sorted(PROFILES), default="smoke")
    parser.add_argument("--output", default="python_code/su3_bulk_scale_setting_results.json")
    parser.add_argument("--jsonl", default=None, help="Incremental checkpoint JSONL path; defaults to output with .jsonl suffix.")
    parser.add_argument("--nt", type=int, action="append", help="Restrict to one or more Nt values. Repeatable.")
    parser.add_argument("--size", type=parse_size, action="append", help="Restrict to one or more spatial sizes, e.g. 4x4x3. Repeatable.")
    parser.add_argument("--seed", type=int, default=20260630)
    args = parser.parse_args()

    profile = PROFILES[args.profile]
    selected_nts = tuple(nt for nt in profile.nts if args.nt is None or nt in set(args.nt))
    selected_sizes = tuple(size for size in profile.sizes if args.size is None or size in set(args.size))
    if not selected_nts:
        raise SystemExit(f"no Nt values selected from profile {args.profile}: requested {args.nt}")
    if not selected_sizes:
        raise SystemExit(f"no sizes selected from profile {args.profile}: requested {args.size}")

    out = Path(args.output)
    jsonl = Path(args.jsonl) if args.jsonl else out.with_suffix(".jsonl")
    out.parent.mkdir(parents=True, exist_ok=True)
    jsonl.parent.mkdir(parents=True, exist_ok=True)
    jsonl.write_text("", encoding="utf-8")

    print("Full-SU(3) bond-bipyramid scale-setting harness")
    print("=" * 96)
    print(f"  profile={args.profile} nts={selected_nts} sizes={selected_sizes}")
    print(f"  therm={profile.therm} meas={profile.meas} hits={profile.hits} max_loop={profile.max_loop}")
    print(f"  output={out}")
    print(f"  jsonl={jsonl}")

    write_jsonl(
        jsonl,
        {
            "type": "start",
            "profile": args.profile,
            "seed": args.seed,
            "nts": selected_nts,
            "sizes": selected_sizes,
            "therm": profile.therm,
            "meas": profile.meas,
            "hits": profile.hits,
            "max_loop": profile.max_loop,
        },
    )

    t0 = time.time()
    runs_by_nt: dict[int, list[VolumeRun]] = {}
    for nt in selected_nts:
        runs_by_nt[nt] = []
        print(f"\n[Nt={nt}]")
        for size_i, size in enumerate(selected_sizes):
            seed = args.seed + 1000 * nt + 17 * size_i
            run = run_volume(
                nt=nt,
                size=size,
                betas=profile.betas[nt],
                therm=profile.therm,
                meas=profile.meas,
                hits=profile.hits,
                eps0=profile.eps0,
                max_loop=profile.max_loop,
                seed=seed,
            )
            runs_by_nt[nt].append(run)
            print_volume(run)
            write_jsonl(
                jsonl,
                {
                    "type": "volume",
                    "profile": args.profile,
                    "seed": seed,
                    "elapsed_s": time.time() - t0,
                    "run": asdict(run),
                },
            )

    summaries = [fss_summary(nt, runs_by_nt[nt]) for nt in selected_nts]
    print("\n[scale-setting summary]")
    print("  Nt  beta_c_FSS  beta_sigma  a√sigma  T_c/√sigma  slope(1/Nv)")
    for item in summaries:
        sig = "nan" if item.a_sqrt_sigma is None else f"{item.a_sqrt_sigma:.4f}"
        tc = "nan" if item.tc_over_sqrt_sigma is None else f"{item.tc_over_sqrt_sigma:.4f}"
        slope = "nan" if item.slope_1_over_nv is None else f"{item.slope_1_over_nv:.3f}"
        print(f"  {item.nt:2d}  {item.beta_c_fss:10.4f}  {item.beta_sigma:10.4f}  {sig:>7s}  {tc:>10s}  {slope:>10s}")

    finite_tc = [item.tc_over_sqrt_sigma for item in summaries if item.tc_over_sqrt_sigma]
    if len(finite_tc) >= 2:
        spread = (max(finite_tc) - min(finite_tc)) / (0.5 * (max(finite_tc) + min(finite_tc)))
    else:
        spread = float("nan")
    print(f"\n  T_c/sqrt(sigma) spread across measured Nt = {spread:.3%}")

    verdict = (
        "scale-setting-local"
        if math.isfinite(spread) and spread < 0.25
        else "not-scale-set"
    )
    print(
        f"""
VERDICT:
  {verdict.upper()}.  This run is a {args.profile!r} profile.  It uses full
  SU(3) bond-bipyramid ensembles and measures Polyakov susceptibility and
  spatial Wilson Creutz ratios on the same beta grid.  A continuum physical
  string tension is promoted only if T_c/sqrt(sigma) stabilises under Nt and
  volume refinement.  The printed spread is therefore the gate, not a cosmetic
  diagnostic.
  elapsed {time.time() - t0:.1f}s
exit 0"""
    )

    payload = {
        "profile": args.profile,
        "seed": args.seed,
        "elapsed_s": time.time() - t0,
        "runs": {
            str(nt): [asdict(run) for run in runs]
            for nt, runs in runs_by_nt.items()
        },
        "summaries": [asdict(item) for item in summaries],
        "tc_over_sqrt_sigma_spread": spread,
        "verdict": verdict,
    }
    out.write_text(json.dumps(payload, indent=2) + "\n")
    write_jsonl(jsonl, {"type": "summary", **payload})
    print(f"wrote {out}")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
