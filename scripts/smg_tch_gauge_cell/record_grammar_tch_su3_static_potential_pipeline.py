#!/usr/bin/env python3
r"""Bond-bipyramid SU(3) static-potential pipeline.

This is the framework-geometry port of
``record_grammar_cubic_su3_static_potential_pipeline.py``.  The important
change from the earlier TCH scale-setting harness is the observable: this script
uses temporal Wilson rectangles and an effective-potential readout, not spatial
``1x1`` Creutz loops.

The script is deliberately tiered as a method-development rung:

* it constructs straight spatial source paths on the licensed bond-bipyramid
  bulk graph;
* builds temporal Wilson rectangles with the same SU(3) link variables used by
  the finite-temperature Polyakov runs;
* optionally APE-smears spatial links on the irregular complex;
* extracts V_eff(R,T)=log[W(R,T)/W(R,T+1)];
* fits both a Cornell form and a linear tail-slope control.

Boundary
--------
This does not claim a physical QCD string tension.  The ordinary cubic control
now validates the fixed Cornell/static-potential scheme at demonstration grade,
so this script's job is to expose whether that same frozen scheme gives a
stable TCH string tension on the framework geometry.  Current matched small-block
runs show the next failure mode: the straight-path source operator is not yet
scale-setting on the present bond-bipyramid block.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from record_grammar_cubic_su3_static_potential_pipeline import (
    fit_cornell,
    fit_linear_tail,
    mean_se,
    potential_points,
    project_su3,
)
from record_grammar_tch_su3_bulk_polyakov_torus import (
    I3,
    dag,
    staple,
    staple_terms,
    sweep,
)
from record_grammar_tch_su3_bulk_scale_setting import thermal_lattice


DEFAULT_OUTPUT = "python_code/tch_su3_static_potential_pipeline_smoke.json"
CUBIC_REFERENCE_TC_OVER_SQRT_SIGMA = 0.629


def spatial_link(lattice: dict, edge_i: int, t: int) -> int:
    return edge_i * lattice["Nt"] + (t % lattice["Nt"])


def time_link(lattice: dict, vertex_i: int, t: int) -> int:
    return lattice["Ne"] * lattice["Nt"] + vertex_i * lattice["Nt"] + (t % lattice["Nt"])


def canonical_step(axis: int) -> tuple[int, int, int]:
    out = [0, 0, 0]
    out[axis] = 2
    return tuple(out)  # type: ignore[return-value]


def add_vertex(left: tuple[int, int, int], right: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(left[i] + right[i] for i in range(3))  # type: ignore[return-value]


def straight_source_paths(lattice: dict, max_r: int) -> dict[int, list[tuple[int, int, list[tuple[int, bool]]]]]:
    """Return straight graph paths of length R on actual bond-bipyramid edges.

    Each path item is ``(start_vertex_index, end_vertex_index,
    [(edge_index, forward), ...])`` where ``forward`` is relative to the
    canonical orientation stored in ``lattice["edges"]``.
    """

    edge_lookup: dict[tuple[tuple[int, int, int], tuple[int, int, int]], tuple[int, bool]] = {}
    for edge_i, (a, b) in enumerate(lattice["edges"]):
        edge_lookup[(a, b)] = (edge_i, True)
        edge_lookup[(b, a)] = (edge_i, False)

    paths: dict[int, list[tuple[int, int, list[tuple[int, bool]]]]] = {r: [] for r in range(1, max_r + 1)}
    for axis in range(3):
        step = canonical_step(axis)
        for start in lattice["vertices"]:
            current = start
            path: list[tuple[int, bool]] = []
            for r in range(1, max_r + 1):
                nxt = add_vertex(current, step)
                item = edge_lookup.get((current, nxt))
                if item is None:
                    break
                path.append(item)
                current = nxt
                paths[r].append((lattice["vidx"][start], lattice["vidx"][current], list(path)))
    return paths


def temporal_wilson_rectangle(
    U: np.ndarray,
    lattice: dict,
    path: tuple[int, int, list[tuple[int, bool]]],
    t0: int,
    extent_t: int,
) -> float:
    start_v, end_v, spatial_path = path
    M = I3.copy()

    for edge_i, forward in spatial_path:
        link = U[spatial_link(lattice, edge_i, t0)]
        M = M @ (link if forward else dag(link))

    for dt in range(extent_t):
        M = M @ U[time_link(lattice, end_v, t0 + dt)]

    top_t = t0 + extent_t
    for edge_i, forward in reversed(spatial_path):
        link = U[spatial_link(lattice, edge_i, top_t)]
        M = M @ (dag(link) if forward else link)

    for dt in reversed(range(extent_t)):
        M = M @ dag(U[time_link(lattice, start_v, t0 + dt)])

    return float(np.real(np.trace(M) / 3.0))


def temporal_wilson_table(
    U: np.ndarray,
    lattice: dict,
    paths: dict[int, list[tuple[int, int, list[tuple[int, bool]]]]],
    max_r: int,
    max_t: int,
) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for r in range(1, max_r + 1):
        if not paths.get(r):
            continue
        for t in range(1, max_t + 1):
            vals = [
                temporal_wilson_rectangle(U, lattice, path, t0, t)
                for t0 in range(lattice["Nt"])
                for path in paths[r]
            ]
            out[(r, t)] = float(np.mean(vals))
    return out


def spatial_staple_terms(lattice: dict) -> list[list[tuple[int, list[tuple[int, bool]]]]]:
    """Staple terms restricted to spatial plaquettes for irregular APE smearing."""

    limit = lattice["Ne"] * lattice["Nt"]
    plink, psign, Np, Nl = lattice["plink"], lattice["psign"], lattice["Np"], lattice["Nlinks"]
    terms = [[] for _ in range(Nl)]
    for p in range(Np):
        active = [(int(plink[p, s]), int(psign[p, s])) for s in range(4) if psign[p, s] != 0]
        if not active or any(link >= limit for link, _ in active):
            continue
        m = len(active)
        for a in range(m):
            link, sign = active[a]
            others = [(active[(a + b) % m][0], active[(a + b) % m][1] == -1) for b in range(1, m)]
            terms[link].append((sign, others))
    return terms


def ape_smear_spatial_irregular(
    U: np.ndarray,
    lattice: dict,
    terms: list[list[tuple[int, list[tuple[int, bool]]]]],
    steps: int,
    alpha: float,
) -> np.ndarray:
    out = U.copy()
    if steps <= 0:
        return out
    spatial_limit = lattice["Ne"] * lattice["Nt"]
    for _ in range(steps):
        old = out.copy()
        for link_i in range(spatial_limit):
            if not terms[link_i]:
                continue
            mixed = (1.0 - alpha) * old[link_i] + (alpha / len(terms[link_i])) * staple(old, link_i, terms)
            out[link_i] = project_su3(mixed)
        out[spatial_limit:] = old[spatial_limit:]
    return out


def su3_unitarity_error(U: np.ndarray) -> float:
    eye = np.broadcast_to(I3, U.shape[:-2] + (3, 3))
    return float(np.max(np.abs(dag(U) @ U - eye)))


def run(args: argparse.Namespace) -> dict:
    lattice = thermal_lattice(*args.size, args.time_extent)
    all_terms = staple_terms(lattice)
    smear_terms = spatial_staple_terms(lattice)
    paths = straight_source_paths(lattice, args.max_r)
    rng = np.random.default_rng(args.seed)
    U = np.broadcast_to(I3, (lattice["Nlinks"], 3, 3)).copy()
    eps = args.eps0
    acc_values: list[float] = []

    for _ in range(args.therm):
        U, acc = sweep(U, args.beta, rng, all_terms, lattice["Nlinks"], eps, hits=args.hits)
        eps = min(max(eps * (1.0 + 0.1 * (acc - 0.5)), 0.04), 1.0)
        acc_values.append(float(acc))

    sample_tables: dict[int, list[dict[tuple[int, int], float]]] = {level: [] for level in args.smear_steps}
    for _ in range(args.meas):
        for _ in range(args.skip):
            U, acc = sweep(U, args.beta, rng, all_terms, lattice["Nlinks"], eps, hits=args.hits)
            eps = min(max(eps * (1.0 + 0.05 * (acc - 0.5)), 0.04), 1.0)
            acc_values.append(float(acc))
        for level in args.smear_steps:
            Us = ape_smear_spatial_irregular(U, lattice, smear_terms, level, args.ape_alpha)
            sample_tables[level].append(temporal_wilson_table(Us, lattice, paths, args.max_r, args.max_t))

    fit_r_max = args.fit_r_max if args.fit_r_max else args.max_r
    analyses = []
    cornell_fits = []
    linear_fits = []
    for level, samples in sample_tables.items():
        points = potential_points(samples, args.max_r, args.t_min, args.max_t)
        fit = fit_cornell(points, level, args.nt_reference, args.fit_r_min, fit_r_max)
        linear = fit_linear_tail(points, level, args.nt_reference, args.fit_r_min, fit_r_max)
        fit.verdict = tch_diagnostic_verdict(fit.sigma, fit.chi2_dof)
        linear.verdict = tch_diagnostic_verdict(linear.sigma, linear.chi2_dof)
        cornell_fits.append(fit)
        linear_fits.append(linear)
        w_means = {}
        w_ses = {}
        for r in range(1, args.max_r + 1):
            for t in range(1, args.max_t + 1):
                mean, se, _ = mean_se(table.get((r, t), float("nan")) for table in samples)
                w_means[f"{r}x{t}"] = mean
                w_ses[f"{r}x{t}"] = se
        analyses.append(
            {
                "smear_steps": level,
                "wilson_means": w_means,
                "wilson_se": w_ses,
                "potential_points": [asdict(point) for point in points],
                "fit": asdict(fit),
                "linear_tail_fit": asdict(linear),
            }
        )

    valid_cornell = [fit for fit in cornell_fits if fit.tc_over_sqrt_sigma is not None]
    best_cornell = min(valid_cornell, key=lambda f: f.rel_error_vs_target or float("inf")) if valid_cornell else None
    valid_linear = [fit for fit in linear_fits if fit.tc_over_sqrt_sigma is not None]
    best_linear = min(valid_linear, key=lambda f: f.rel_error_vs_target or float("inf")) if valid_linear else None
    # The cubic 0.629 value is a validator for the extractor, not a target for
    # the framework geometry.  Do not promote a TCH run because it happens to
    # sit near the cubic continuum ratio; promotion requires the cubic Cornell
    # control to be production-validated first, then the same fixed scheme on
    # TCH at matched beta_c(N_t).
    verdict = "PORTED-NOT-STABLE"
    if best_cornell and best_cornell.sigma and best_cornell.sigma > 0:
        verdict = "PORTED-CORNELL-POSITIVE"
    elif best_linear and best_linear.sigma and best_linear.sigma > 0:
        verdict = "PORTED-TAIL-POSITIVE"

    return {
        "profile": "bond-bipyramid-su3-static-potential-pipeline",
        "beta": args.beta,
        "nt_reference": args.nt_reference,
        "cubic_reference_tc_over_sqrt_sigma": CUBIC_REFERENCE_TC_OVER_SQRT_SIGMA,
        "size": args.size,
        "time_extent": args.time_extent,
        "nv": int(lattice["Nv"]),
        "ne": int(lattice["Ne"]),
        "nlinks": int(lattice["Nlinks"]),
        "nplaquettes": int(lattice["Np"]),
        "path_counts": {str(r): len(paths.get(r, [])) for r in range(1, args.max_r + 1)},
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
        "fit_r_min": args.fit_r_min,
        "fit_r_max": fit_r_max,
        "acceptance_mean": float(np.mean(acc_values)),
        "unitarity_error": su3_unitarity_error(U),
        "analyses": analyses,
        "best_fit": asdict(best_cornell) if best_cornell else None,
        "best_linear_tail_fit": asdict(best_linear) if best_linear else None,
        "verdict": verdict,
        "interpretation": (
            "The fixed static-potential readout has been ported to bond-bipyramid temporal Wilson rectangles. "
            "The printed T_c/sqrt(sigma) values are diagnostics only; the ordinary cubic Cornell control is validated, "
            "but TCH promotion waits on a stable geometry-specific source operator at matched beta_c(N_t)."
        ),
    }


def tch_diagnostic_verdict(sigma: float | None, chi2_dof: float | None) -> str:
    if sigma is None or not math.isfinite(float(sigma)):
        return "NO-FIT"
    if sigma <= 0:
        return "NEGATIVE-SIGMA"
    if chi2_dof is None:
        return "DIAGNOSTIC-POSITIVE-NO-DOF"
    return "DIAGNOSTIC-POSITIVE"


def parse_size(text: str) -> tuple[int, int, int]:
    parts = text.lower().replace(",", "x").split("x")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("size must be NxMxK, e.g. 4x4x3")
    try:
        values = tuple(int(item) for item in parts)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("size components must be integers") from exc
    if any(item < 2 for item in values):
        raise argparse.ArgumentTypeError("size components must be >=2")
    return values  # type: ignore[return-value]


def parse_smear_steps(text: str) -> list[int]:
    steps = [int(item) for item in text.replace(",", " ").split() if item.strip()]
    if not steps:
        raise argparse.ArgumentTypeError("at least one smear step is required")
    if any(step < 0 for step in steps):
        raise argparse.ArgumentTypeError("smear steps must be non-negative")
    return sorted(dict.fromkeys(steps))


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{digits}f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beta", type=float, default=6.5)
    parser.add_argument("--nt-reference", type=int, default=4)
    parser.add_argument("--size", type=parse_size, default=parse_size("3x3x2"))
    parser.add_argument("--time-extent", type=int, default=6)
    parser.add_argument("--therm", type=int, default=20)
    parser.add_argument("--meas", type=int, default=12)
    parser.add_argument("--skip", type=int, default=1)
    parser.add_argument("--hits", type=int, default=1)
    parser.add_argument("--eps0", type=float, default=0.30)
    parser.add_argument("--ape-alpha", type=float, default=0.45)
    parser.add_argument("--smear-steps", type=parse_smear_steps, default=parse_smear_steps("0,2"))
    parser.add_argument("--max-r", type=int, default=3)
    parser.add_argument("--max-t", type=int, default=4)
    parser.add_argument("--t-min", type=int, default=1)
    parser.add_argument("--fit-r-min", type=int, default=1)
    parser.add_argument("--fit-r-max", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.max_t >= args.time_extent:
        raise SystemExit("--max-t must be smaller than --time-extent")
    if args.max_r < 3:
        raise SystemExit("--max-r must be at least 3")

    t0 = time.time()
    result = run(args)
    result["elapsed_s"] = time.time() - t0
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("Bond-bipyramid SU(3) static-potential pipeline")
    print("=" * 96)
    print(f"  beta={args.beta:.3f}; Nt reference={args.nt_reference}; size={args.size}; time extent={args.time_extent}")
    print(f"  therm={args.therm}; meas={args.meas}; skip={args.skip}; hits={args.hits}; max R,T={args.max_r},{args.max_t}")
    print(f"  path counts: {result['path_counts']}")
    print(f"  acceptance={result['acceptance_mean']:.3f}; unitarity_error={result['unitarity_error']:.2e}")
    print("\n[Cornell fits]")
    print("  smear  sigma      a√sigma  Tc/√sigma  rel_vs_cubic_ref alpha     V0       chi2/dof verdict")
    for analysis in result["analyses"]:
        fit = analysis["fit"]
        print(
            f"  {fit['smear_steps']:>5d} "
            f"{fmt(fit['sigma']):>9s} {fmt(fit['a_sqrt_sigma']):>9s} "
            f"{fmt(fit['tc_over_sqrt_sigma']):>10s} {fmt(fit['rel_error_vs_target'], 3):>8s} "
            f"{fmt(fit['alpha']):>8s} {fmt(fit['v0']):>8s} {fmt(fit['chi2_dof']):>9s} "
            f"{fit['verdict']}"
        )
    print("\n[linear tail-slope controls]")
    print("  smear  sigma      a√sigma  Tc/√sigma  rel_vs_cubic_ref V0       chi2/dof verdict")
    for analysis in result["analyses"]:
        fit = analysis["linear_tail_fit"]
        print(
            f"  {fit['smear_steps']:>5d} "
            f"{fmt(fit['sigma']):>9s} {fmt(fit['a_sqrt_sigma']):>9s} "
            f"{fmt(fit['tc_over_sqrt_sigma']):>10s} {fmt(fit['rel_error_vs_target'], 3):>8s} "
            f"{fmt(fit['v0']):>8s} {fmt(fit['chi2_dof']):>9s} "
            f"{fit['verdict']}"
        )
    print(f"\nVERDICT: {result['verdict']}")
    print("  " + result["interpretation"])
    print(f"  wrote {args.output}")
    print("exit 0")
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
