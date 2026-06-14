#!/usr/bin/env python3
"""Frontier ledger for the actual TCH continuum/gauged mirror-gap problem.

This is a triage script, not another expensive DMRG driver. It reduces the
already-computed local/proxy artifacts into the one remaining TCH-specific gate:

    construct the actual gauge-cell dressed mirror interaction, then test its
    symmetric gap under Gauss projection / Wilson fluctuations.

The script deliberately separates four things that were easy to conflate:

* closed preconditions: SM anomalies, Witten, trivial cobordism;
* validated mechanism: the 3-4-5-0 reference SMG harness;
* demoted/negative proxies: the Z3 colour strip and walk-sourced overlap kernel;
* still-open TCH object: the non-abelian SU(2)_L chiral measure or, equivalently
  the actual gauge-dressed mirror-gap operator on the TCH gauge web.

Standard library only. Writes a compact JSON verdict.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PYTHON_CODE = ROOT / "python_code"
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_gauged_mirror_frontier_verdict.json"


FERMIONS = [
    ("Q_L", F(1, 6), 2, +3),
    ("u_R^c", F(-2, 3), 1, -3),
    ("d_R^c", F(1, 3), 1, -3),
    ("L_L", F(-1, 2), 2, +1),
    ("e_R^c", F(1), 1, +1),
    ("nu_R^c", F(0), 1, +1),
]


def mult(field: tuple[str, F, int, int]) -> int:
    return field[2] * abs(field[3])


def anomaly_preconditions() -> dict[str, Any]:
    anomalies = {
        "n_weyl": sum(mult(f) for f in FERMIONS),
        "grav_u1": sum(mult(f) * f[1] for f in FERMIONS),
        "u1_cubic": sum(mult(f) * f[1] ** 3 for f in FERMIONS),
        "u1_su2_squared": sum(abs(f[3]) * f[1] for f in FERMIONS if f[2] == 2),
        "u1_su3_squared": sum(f[2] * f[1] for f in FERMIONS if abs(f[3]) == 3),
        "su3_cubic": sum(
            f[2] * (1 if f[3] > 0 else -1) for f in FERMIONS if abs(f[3]) == 3
        ),
        "su2_doublets_mod2": sum(abs(f[3]) for f in FERMIONS if f[2] == 2) % 2,
        "z16_weyl_mod16": sum(mult(f) for f in FERMIONS) % 16,
    }
    assert anomalies == {
        "n_weyl": 16,
        "grav_u1": 0,
        "u1_cubic": 0,
        "u1_su2_squared": 0,
        "u1_su3_squared": 0,
        "su3_cubic": 0,
        "su2_doublets_mod2": 0,
        "z16_weyl_mod16": 0,
    }
    return {
        "status": "closed_precondition",
        "meaning": "SM gauge/gravity/global anomalies vanish; U(1)_Y Luescher precondition is met.",
        "anomalies": {k: str(v) for k, v in anomalies.items()},
    }


def route2_gap_minresid(traj: list[list[float | None]]) -> tuple[float, float, int]:
    usable = [(int(i), float(g), float(r)) for i, g, r in traj if r is not None and g > 1e-3]
    i, gap, resid = min(usable, key=lambda x: x[2])
    return gap, resid, i


def route2_colour_proxy() -> dict[str, Any]:
    data_path = PYTHON_CODE / "route2_strip_tridiagonals.json"
    data = json.loads(data_path.read_text())
    gaps: dict[str, float] = {}
    details: dict[str, dict[str, Any]] = {}
    for key in ["3_0.2", "3_4.0", "4_0.2", "4_4.0", "6_0.2", "6_4.0"]:
        row = data[key]
        gap, resid, iteration = route2_gap_minresid(row["gap_traj"])
        gaps[key] = gap
        details[key] = {
            "states_per_charge": row["spc"],
            "t": row["t"],
            "gap_min_resid": gap,
            "resid": resid,
            "iteration": iteration,
        }
    collapse = (gaps["6_0.2"] - gaps["6_4.0"]) / gaps["6_0.2"]
    assert gaps["6_4.0"] < 0.5
    assert collapse > 0.8
    return {
        "status": "computed_negative_demoted_proxy",
        "meaning": (
            "The Z3/colour strip gap collapses at n/q=6, so the old CSS colour proxy "
            "does not supply continuum protection. It is demoted because colour is vector-like; "
            "it is not the actual SU(2)_L chiral obstruction."
        ),
        "nq6_t4_gap": gaps["6_4.0"],
        "nq6_t_response_collapse": collapse,
        "details": details,
    }


def walk_kernel_obstruction() -> dict[str, Any]:
    even = []
    odd = []
    for bits in [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]:
        parity = sum(bits) % 2
        kx, ky, kz = (math.pi * bit for bit in bits)
        v = [
            8 * math.sin(kx) * math.cos(ky) * math.cos(kz),
            8 * math.cos(kx) * math.sin(ky) * math.cos(kz),
            8 * math.cos(kx) * math.cos(ky) * math.sin(kz),
        ]
        w_body = 8 * (1 - math.cos(kx) * math.cos(ky) * math.cos(kz))
        norm = math.sqrt(sum(x * x for x in v) + 4 * w_body * w_body)
        (even if parity == 0 else odd).append({"corner": bits, "norm": norm})
    assert max(row["norm"] for row in even) < 1e-9
    assert min(row["norm"] for row in odd) > 1.0
    return {
        "status": "native_walk_kernel_fails",
        "meaning": (
            "The walk-sourced Wilson/overlap kernel vanishes at all four even Brillouin-zone "
            "corners, leaving unlifted momentum tastes. A single chiral species needs an imported "
            "exponentially-local overlap kernel, not the native ultralocal walk kernel."
        ),
        "even_corner_norms": even,
        "odd_corner_norms": odd,
    }


def reference_smg_harness() -> dict[str, Any]:
    verdict_path = ROOT / "python_code" / "smg_dmrg_runs" / "recovered_verdict.json"
    if not verdict_path.exists():
        return {
            "status": "missing_local_verdict",
            "meaning": "Run python_code/smg_dmrg/verdict_from_runs.py to regenerate the reference-harness verdict.",
        }
    data = json.loads(verdict_path.read_text())
    verdict = data.get("verdict", data)
    status = verdict.get("status", "")
    assert "reference 3-4-5-0 SMG benchmark confirmed" in status
    return {
        "status": "validated_reference_mechanism",
        "meaning": (
            "The TeNPy/DMRG pipeline reproduces the published 3-4-5-0 SMG benchmark. "
            "This validates the mechanism and harness, not the TCH gauge-coupled construction."
        ),
        "source_status": status,
        "scope": verdict.get("phase1_question") or verdict.get("scope"),
    }


def gauge_dressing_gate() -> dict[str, Any]:
    required = [
        PYTHON_CODE / "css_gauge_compatibility.py",
        PYTHON_CODE / "css_tch_operator_candidates.py",
        PYTHON_CODE / "css_ks_wilson_proxy_gap.py",
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    assert not missing, missing
    return {
        "status": "finite_next_gate_identified",
        "meaning": (
            "Bare CSS X-stabilizers are not gauge invariant, but css_gauge_compatibility.py "
            "constructs finite dual-orbit matter-plus-link dressings. css_tch_operator_candidates.py "
            "selects the closed-flux/Bianchi X-plaquette plus Gauss-orbit quotient as the concrete "
            "operator pattern to seek on the actual TCH gauge cell."
        ),
        "known_computed_numbers": {
            "dual_orbit_dimensions": [65, 104],
            "dressed_commutator_residual_scale": "1e-14",
            "candidate_internal_gap_scan": "positive at finite beta after closed-flux projection; trend not a continuum proof",
        },
        "next_test": [
            "identify the dual-orbit compensator/link spaces with TCH truncated-cube gauge-cell Wilson/plaquette DOF",
            "construct a Hermitian gauge-cell X operator whose boundary product is prod_l X_l",
            "include vertex Gauss generators that move local signs/frames while preserving the boundary product",
            "run the dressed mirror-gap survival test under background Wilson twists or dynamical gauge fluctuations",
        ],
    }


def build_verdict() -> dict[str, Any]:
    checks = {
        "anomaly_preconditions": anomaly_preconditions(),
        "reference_smg_harness": reference_smg_harness(),
        "route2_colour_proxy": route2_colour_proxy(),
        "walk_kernel_obstruction": walk_kernel_obstruction(),
        "gauge_dressing_gate": gauge_dressing_gate(),
    }
    return {
        "verdict": "progress_without_closure",
        "summary": (
            "The actual TCH continuum/gauged mirror-gap problem is not closed. Progress is a "
            "sharper reduction: the old reference/proxy computations are reconciled, and the "
            "remaining finite TCH-specific task is the explicit gauge-cell dressed mirror operator "
            "plus its gap test. Full continuum closure still rests on the field-open non-abelian "
            "SU(2)_L chiral measure or an equivalent successful symmetric mirror-gapping construction."
        ),
        "checks": checks,
        "do_not_spend_compute_on": [
            "rerunning the already-validated 3-4-5-0 reference harness",
            "escalating the demoted Z3/vector-like colour strip as if it were the SU(2)_L chiral problem",
            "trying to bypass the Luescher cohomology obstruction with finite-codebook/topology language",
        ],
        "highest_value_next_artifact": (
            "a TCH gauge-cell operator implementation that realizes the finite dual-orbit dressing and "
            "closed-flux/Bianchi boundary product, followed by a dressed-gap survival audit"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    print("TCH gauged mirror-gap frontier verdict")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    for name, check in verdict["checks"].items():
        print(f"  {name}: {check['status']}")
    print("\nNext finite gate:")
    print(f"  {verdict['highest_value_next_artifact']}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
