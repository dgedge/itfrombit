"""Command-line guards for the SMG DMRG prototype."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .model_3450 import (
    MASS_PAIRS,
    Model3450Spec,
    centered_left_position,
    centered_left_unit_cell,
    build_3450_terms,
    free_spectrum_3450,
    mass_correlation_term,
    paper_mass_correlation_term,
    paper_single_fermion_correlation_term,
    single_fermion_correlation_term,
)
from .ed import solve_exact
from .models import ModelSpec, build_model_terms
from .tenpy_backend import require_tenpy, run_ground_dmrg, run_ground_dmrg_3450
from .terms import charge_violations


def parse_u_values(raw: str) -> list[float]:
    return [float(value.strip()) for value in raw.split(",") if value.strip()]


def parse_optional_float(raw: str) -> float | None:
    if raw.lower() in {"none", "off", "false", "no"}:
        return None
    return float(raw)


def spec_from_args(args: argparse.Namespace, *, u: float | None = None) -> ModelSpec:
    return ModelSpec(
        cells=args.cells,
        flavors=args.flavors,
        t=args.t,
        r=args.r,
        mass=args.mass,
        u=args.u if u is None else u,
        mirror_charge_sign=args.mirror_charge_sign,
        opposite_wilson_for_mirror=not args.same_wilson,
        quartic_sign=args.quartic_sign,
    )


def spec_3450_from_args(args: argparse.Namespace) -> Model3450Spec:
    g = getattr(args, "g", None)
    g1 = args.g1 if args.g1 is not None else (0.0 if g is None else g)
    g2 = args.g2 if args.g2 is not None else (0.0 if g is None else g)
    return Model3450Spec(
        unit_cells=args.unit_cells,
        t1=args.t1,
        t2=args.t2,
        g1=g1,
        g2=g2,
        include_free=not args.no_free,
        include_interaction=not args.no_interaction,
    )


def print_model_summary(spec: ModelSpec) -> tuple[int, int]:
    modes, terms = build_model_terms(spec)
    violations = charge_violations(terms, modes)
    print(f"cells                 = {spec.cells}")
    print(f"flavors               = {spec.flavors}")
    print(f"total modes           = {spec.total_modes}")
    print(f"term count            = {len(terms)}")
    print(f"charge violations     = {len(violations)}")
    if violations:
        for term in violations[:5]:
            print(f"  violation: {term.label} {term}")
    return len(terms), len(violations)


def cmd_ed_sanity(args: argparse.Namespace) -> None:
    spec = spec_from_args(args)
    modes, terms = build_model_terms(spec)
    violations = charge_violations(terms, modes)
    print_model_summary(spec)
    if violations:
        raise SystemExit("charge-conservation gate failed")
    result = solve_exact(spec, terms)
    print(f"hermiticity_residual  = {result.hermiticity_residual:.3e}")
    print(f"many_body_gap         = {result.gap:.12g}")
    print(f"mirror_bilinear_m2    = {result.mirror_bilinear_m2:.12g}")
    print("lowest_levels         = " + ", ".join(f"{x:.8g}" for x in result.eigenvalues[:8]))
    if result.hermiticity_residual > 1e-9:
        raise SystemExit("Hermiticity gate failed")
    print("ALL ED GUARDS PASSED.")


def cmd_scan_ed(args: argparse.Namespace) -> None:
    print(f"{'U':>10s} {'gap':>16s} {'m2':>16s} {'herm':>12s}")
    for value in parse_u_values(args.u_values):
        spec = spec_from_args(args, u=value)
        modes, terms = build_model_terms(spec)
        violations = charge_violations(terms, modes)
        if violations:
            raise SystemExit(f"charge-conservation failed at U={value}")
        result = solve_exact(spec, terms)
        print(
            f"{value:10.4g} {result.gap:16.8g} "
            f"{result.mirror_bilinear_m2:16.8g} {result.hermiticity_residual:12.3e}"
        )


def cmd_write_config(args: argparse.Namespace) -> None:
    spec = spec_from_args(args)
    path = Path(args.output)
    path.write_text(json.dumps(asdict(spec), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def cmd_export_terms(args: argparse.Namespace) -> None:
    spec = spec_from_args(args)
    modes, terms = build_model_terms(spec)
    payload = {
        "spec": asdict(spec),
        "modes": [asdict(mode) for mode in modes],
        "terms": [
            {
                "coeff": [term.coeff.real, term.coeff.imag],
                "ops": [asdict(op) for op in term.ops],
                "label": term.label,
            }
            for term in terms
        ],
        "charge_violations": len(charge_violations(terms, modes)),
    }
    path = Path(args.output)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def cmd_dmrg_ground(args: argparse.Namespace) -> None:
    spec = spec_from_args(args)
    modes, terms = build_model_terms(spec)
    violations = charge_violations(terms, modes)
    if violations:
        raise SystemExit(f"charge-conservation gate failed: {len(violations)} violations")
    try:
        result = run_ground_dmrg(
            spec,
            chi_max=args.chi_max,
            max_sweeps=args.max_sweeps,
            svd_min=args.svd_min,
            mixer=not args.no_mixer,
        )
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"energy                = {result['energy']}")
    print(f"mirror_bilinear_m2    = {result['mirror_bilinear_m2']:.12g}")
    print("DMRG ground-state run finished.")


def cmd_tenpy_check(_: argparse.Namespace) -> None:
    try:
        tenpy = require_tenpy()
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"TeNPy available: {getattr(tenpy, '__version__', 'unknown')}")


def print_3450_summary(spec: Model3450Spec) -> tuple[int, int]:
    modes, terms = build_3450_terms(spec)
    violations = charge_violations(terms, modes)
    print(f"unit_cells            = {spec.unit_cells}")
    print(f"lattice orbitals      = {spec.lattice_orbitals}")
    print(f"edge sites            = {spec.edge_sites}")
    print(f"total modes           = {spec.total_modes}")
    print(f"t1, t2                = {spec.t1}, {spec.t2}")
    print(f"g1, g2                = {spec.g1}, {spec.g2}")
    print(f"term count            = {len(terms)}")
    print(f"charge violations     = {len(violations)}")
    if violations:
        for term in violations[:8]:
            print(f"  violation: {term.label} {term}")
    return len(terms), len(violations)


def cmd_3450_check(args: argparse.Namespace) -> None:
    spec = spec_3450_from_args(args)
    print_3450_summary(spec)
    modes, terms = build_3450_terms(spec)
    violations = charge_violations(terms, modes)
    if violations:
        raise SystemExit("3-4-5-0 charge-conservation gate failed")
    print("ALL 3-4-5-0 TERM-LIST GUARDS PASSED.")


def cmd_3450_free_spectrum(args: argparse.Namespace) -> None:
    spec = spec_3450_from_args(args)
    evals, hermiticity, gap = free_spectrum_3450(spec)
    half = len(evals) // 2
    window = args.window
    print(f"unit_cells            = {spec.unit_cells}")
    print(f"total modes           = {spec.total_modes}")
    print(f"free hermiticity      = {hermiticity:.3e}")
    print(f"half-fill gap         = {gap:.12g}")
    lo = max(0, half - window)
    hi = min(len(evals), half + window)
    print("levels near E=0       = " + ", ".join(f"{x:.8g}" for x in evals[lo:hi]))
    if hermiticity > 1e-9:
        raise SystemExit("free single-particle Hermiticity gate failed")


def cmd_3450_export_terms(args: argparse.Namespace) -> None:
    spec = spec_3450_from_args(args)
    modes, terms = build_3450_terms(spec)
    payload = {
        "model": "3450",
        "spec": asdict(spec),
        "modes": [asdict(mode) for mode in modes],
        "terms": [
            {
                "coeff": [term.coeff.real, term.coeff.imag],
                "ops": [asdict(op) for op in term.ops],
                "label": term.label,
            }
            for term in terms
        ],
        "charge_violations": len(charge_violations(terms, modes)),
    }
    path = Path(args.output)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def cmd_3450_observable_terms(args: argparse.Namespace) -> None:
    distance = args.distance
    if args.position_grid == "paper":
        left = centered_left_unit_cell(args.unit_cells, distance) if args.centered else args.left_position
        fermion_builder = paper_single_fermion_correlation_term
        mass_builder = paper_mass_correlation_term
    else:
        left = centered_left_position(2 * args.unit_cells, distance) if args.centered else args.left_position
        fermion_builder = single_fermion_correlation_term
        mass_builder = mass_correlation_term
    payload = {
        "model": "3450",
        "edge": args.edge,
        "position_grid": args.position_grid,
        "left_position": left,
        "distance": distance,
        "fermion_correlators": {
            flavor: fermion_builder(args.edge, flavor, left, distance)
            for flavor in ("3", "4", "5", "0")
        },
        "mass_correlators": {
            f"{kind}_{a}{b}": mass_builder(args.edge, kind, (a, b), left, distance)
            for kind in ("dirac", "majorana")
            for a, b in MASS_PAIRS
        },
    }
    path = Path(args.output)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


def cmd_3450_dmrg_ground(args: argparse.Namespace) -> None:
    spec = spec_3450_from_args(args)
    modes, terms = build_3450_terms(spec)
    violations = charge_violations(terms, modes)
    if violations:
        raise SystemExit(f"3-4-5-0 charge-conservation gate failed: {len(violations)} violations")
    try:
        result = run_ground_dmrg_3450(
            spec,
            chi_max=args.chi_max,
            max_sweeps=args.max_sweeps,
            svd_min=args.svd_min,
            mixer=not args.no_mixer,
            conserve=None if args.conserve == "none" else args.conserve,
            max_trunc_err=args.max_trunc_err,
            combine=not args.no_combine,
        )
    except (RuntimeError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    print(f"energy                = {result['energy']}")
    if args.measure_correlations:
        correlations = measure_3450_correlations(
            result["psi"],
            spec,
            edge=args.measure_edge,
            max_distance=args.max_distance,
            position_grid=args.position_grid,
        )
        path = Path(args.correlation_output)
        path.write_text(json.dumps(correlations, indent=2) + "\n", encoding="utf-8")
        print(f"correlations          = {path}")
    print("3-4-5-0 DMRG ground-state run finished.")


def complex_payload(value: complex) -> list[float]:
    return [float(value.real), float(value.imag)]


def measure_3450_correlations(
    psi,
    spec: Model3450Spec,
    *,
    edge: str,
    max_distance: int,
    position_grid: str,
) -> dict:
    """Measure paper-style correlators on an MPS returned by TeNPy DMRG."""

    if max_distance <= 0:
        raise ValueError("max_distance must be positive")
    payload = {
        "model": "3450",
        "edge": edge,
        "position_grid": position_grid,
        "unit_cells": spec.unit_cells,
        "g1": spec.g1,
        "g2": spec.g2,
    }
    edges = ("A", "B") if edge == "both" else (edge,)
    if len(edges) == 1:
        payload["distances"] = measure_3450_correlations_one_edge(
            psi,
            spec,
            edge=edges[0],
            max_distance=max_distance,
            position_grid=position_grid,
        )
    else:
        payload["edges"] = {
            one_edge: measure_3450_correlations_one_edge(
                psi,
                spec,
                edge=one_edge,
                max_distance=max_distance,
                position_grid=position_grid,
            )
            for one_edge in edges
        }
    return payload


def measure_3450_correlations_one_edge(
    psi,
    spec: Model3450Spec,
    *,
    edge: str,
    max_distance: int,
    position_grid: str,
) -> list[dict]:
    if position_grid == "paper":
        max_distance = min(max_distance, spec.unit_cells - 1)
        left_builder = lambda distance: centered_left_unit_cell(spec.unit_cells, distance)
        fermion_builder = paper_single_fermion_correlation_term
        mass_builder = paper_mass_correlation_term
    elif position_grid == "edge":
        max_distance = min(max_distance, spec.edge_sites - 1)
        left_builder = lambda distance: centered_left_position(spec.edge_sites, distance)
        fermion_builder = single_fermion_correlation_term
        mass_builder = mass_correlation_term
    else:
        raise ValueError(f"unknown position_grid {position_grid!r}")
    distances = []
    for distance in range(1, max_distance + 1):
        left = left_builder(distance)
        entry = {
            "distance": distance,
            "left_position": left,
            "fermion": {},
            "mass": {},
        }
        for flavor in ("3", "4", "5", "0"):
            term = fermion_builder(edge, flavor, left, distance)
            entry["fermion"][flavor] = complex_payload(
                psi.expectation_value_term(term, autoJW=True)
            )
        for kind in ("dirac", "majorana"):
            for pair in MASS_PAIRS:
                key = f"{kind}_{pair[0]}{pair[1]}"
                term = mass_builder(edge, kind, pair, left, distance)
                entry["mass"][key] = complex_payload(
                    psi.expectation_value_term(term, autoJW=True)
                )
        distances.append(entry)
    return distances


def add_model_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--cells", type=int, default=2)
    parser.add_argument("--flavors", type=int, default=4)
    parser.add_argument("--t", type=float, default=1.0)
    parser.add_argument("--r", type=float, default=1.0)
    parser.add_argument("--mass", type=float, default=0.0)
    parser.add_argument("--u", type=float, default=1.0)
    parser.add_argument("--mirror-charge-sign", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--same-wilson", action="store_true")
    parser.add_argument("--quartic-sign", type=float, choices=(-1.0, 1.0), default=1.0)


def add_3450_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--unit-cells", type=int, default=20)
    parser.add_argument("--t1", type=float, default=1.0)
    parser.add_argument("--t2", type=float, default=0.5)
    parser.add_argument("--g", type=float, default=None, help="set g1=g2=g")
    parser.add_argument("--g1", type=float, default=None)
    parser.add_argument("--g2", type=float, default=None)
    parser.add_argument("--no-free", action="store_true")
    parser.add_argument("--no-interaction", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    ed = sub.add_parser("ed-sanity", help="run small dense exact-diagonalisation guards")
    add_model_args(ed)
    ed.set_defaults(func=cmd_ed_sanity)

    scan = sub.add_parser("scan-ed", help="scan U with small dense exact diagonalisation")
    add_model_args(scan)
    scan.add_argument("--u-values", default="0,0.5,1,2")
    scan.set_defaults(func=cmd_scan_ed)

    cfg = sub.add_parser("write-config", help="write a JSON model config")
    add_model_args(cfg)
    cfg.add_argument("--output", default="smg_dmrg_config.json")
    cfg.set_defaults(func=cmd_write_config)

    export = sub.add_parser("export-terms", help="export generated modes and terms as JSON")
    add_model_args(export)
    export.add_argument("--output", default="smg_dmrg_terms.json")
    export.set_defaults(func=cmd_export_terms)

    dmrg = sub.add_parser("dmrg-ground", help="run finite TeNPy DMRG ground-state calculation")
    add_model_args(dmrg)
    dmrg.add_argument("--chi-max", type=int, default=500)
    dmrg.add_argument("--max-sweeps", type=int, default=20)
    dmrg.add_argument("--svd-min", type=float, default=1e-10)
    dmrg.add_argument("--no-mixer", action="store_true")
    dmrg.set_defaults(func=cmd_dmrg_ground)

    check = sub.add_parser("tenpy-check", help="check whether TeNPy is importable")
    check.set_defaults(func=cmd_tenpy_check)

    check3450 = sub.add_parser("3450-check", help="check the 3-4-5-0 term list and U(1)xU(1) charges")
    add_3450_args(check3450)
    check3450.set_defaults(func=cmd_3450_check)

    spectrum3450 = sub.add_parser("3450-free-spectrum", help="single-particle free-ladder spectrum check")
    add_3450_args(spectrum3450)
    spectrum3450.add_argument("--window", type=int, default=8)
    spectrum3450.set_defaults(func=cmd_3450_free_spectrum)

    export3450 = sub.add_parser("3450-export-terms", help="export 3-4-5-0 modes and terms as JSON")
    add_3450_args(export3450)
    export3450.add_argument("--output", default="smg_dmrg_3450_terms.json")
    export3450.set_defaults(func=cmd_3450_export_terms)

    obs3450 = sub.add_parser("3450-observable-terms", help="export paper-style fermion and mass correlator terms")
    obs3450.add_argument("--unit-cells", type=int, default=20)
    obs3450.add_argument("--edge", choices=("A", "B"), default="B")
    obs3450.add_argument("--position-grid", choices=("paper", "edge"), default="paper")
    obs3450.add_argument("--distance", type=int, default=2)
    obs3450.add_argument("--left-position", type=int, default=0)
    obs3450.add_argument("--centered", action="store_true")
    obs3450.add_argument("--output", default="smg_dmrg_3450_observables.json")
    obs3450.set_defaults(func=cmd_3450_observable_terms)

    dmrg3450 = sub.add_parser("3450-dmrg-ground", help="run finite TeNPy DMRG for the 3-4-5-0 model")
    add_3450_args(dmrg3450)
    dmrg3450.add_argument("--chi-max", type=int, default=800)
    dmrg3450.add_argument("--max-sweeps", type=int, default=20)
    dmrg3450.add_argument("--svd-min", type=float, default=1e-10)
    dmrg3450.add_argument("--no-mixer", action="store_true")
    dmrg3450.add_argument(
        "--conserve",
        default="q3450",
        help="q3450 conserves Q1,Q2,parity; use none only for debugging the unblocked path",
    )
    dmrg3450.add_argument(
        "--max-trunc-err",
        type=parse_optional_float,
        default=None,
        help="TeNPy consistency threshold; default none downgrades high truncation error to a warning",
    )
    dmrg3450.add_argument("--no-combine", action="store_true")
    dmrg3450.add_argument("--measure-correlations", action="store_true")
    dmrg3450.add_argument("--measure-edge", choices=("A", "B", "both"), default="B")
    dmrg3450.add_argument("--position-grid", choices=("paper", "edge"), default="paper")
    dmrg3450.add_argument("--max-distance", type=int, default=8)
    dmrg3450.add_argument("--correlation-output", default="smg_dmrg_3450_correlations.json")
    dmrg3450.set_defaults(func=cmd_3450_dmrg_ground)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
