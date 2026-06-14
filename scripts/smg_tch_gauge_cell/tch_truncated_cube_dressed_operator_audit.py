#!/usr/bin/env python3
"""Audit the TCH truncated-cube dressed-operator route.

This script investigates the next finite target after the frontier triage:

    actual gauge-cell dressed operator  ->  dressed mirror-gap survival audit

The honest status is layered:

* The TCH architecture requires bridge/gauge-cell degrees of freedom, not
  boundary-sharing between matter cells.
* The SU(3)-center piece of the dressed CSS X operator is constructible:
  split each X stabilizer by triality, attach a qutrit clock compensator, and
  impose a Z3 Bianchi/closed-flux projector on the four oriented boundary links.
* The finite open-link Z3/qutrit strong-coupling gap survives, and the charged
  mirror-Fock chain keeps the local SMG matter gap open in a small Gauss-law
  audit.
* The first non-center Peter-Weyl SU(3) truncation does NOT survive with a
  Wilson plaquette alone; a gauge-invariant TCH/intertwiner lock is required.

Therefore this is progress, not closure: the center-layer dressed operator is
real, but the full truncated-cube SU(3)/SU(2)_L gauge-cell operator is still
missing unless the TCH lock term and its coefficient are derived from geometry.

Standard library + numpy. Writes a JSON verdict.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PYTHON_CODE = ROOT / "python_code"
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_truncated_cube_dressed_operator_verdict.json"


def load_module(name: str, relative_path: str):
    path = PYTHON_CODE / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def center_dressed_operator_audit() -> dict[str, Any]:
    comp = load_module("css_su3_z3_compensator_mod", "css_su3_z3_compensator.py")

    matter_center = comp.matter_center_unitary()
    clock, shift = comp.qutrit_ops()
    total_center = np.kron(matter_center, shift)
    identity = np.eye(256 * 3, dtype=complex)

    rows = []
    max_bare_residual = 0.0
    max_dressed_residual = 0.0
    max_hermiticity = 0.0
    max_involution = 0.0
    for index, codeword in enumerate(comp.GENERATORS):
        bare = sum(comp.x_components_by_triality(codeword))
        dressed, components = comp.dressed_x_on_one_center_link(codeword)
        bare_residual = float(np.linalg.norm(matter_center @ bare @ matter_center.conj().T - bare))
        dressed_residual = float(np.linalg.norm(total_center @ dressed @ total_center.conj().T - dressed))
        hermiticity = float(np.linalg.norm(dressed - dressed.conj().T))
        involution = float(np.linalg.norm(dressed @ dressed - identity))
        component_counts = [
            int(np.count_nonzero(np.abs(component) > 1e-12))
            for component in components
        ]
        max_bare_residual = max(max_bare_residual, bare_residual)
        max_dressed_residual = max(max_dressed_residual, dressed_residual)
        max_hermiticity = max(max_hermiticity, hermiticity)
        max_involution = max(max_involution, involution)
        rows.append(
            {
                "x_stabilizer": index,
                "triality_transition_entries_q012": component_counts,
                "bare_center_residual": bare_residual,
                "dressed_center_residual": dressed_residual,
                "hermiticity_error": hermiticity,
                "involution_error": involution,
            }
        )

    assert max_bare_residual > 1.0
    assert max_dressed_residual < 1e-10
    assert max_hermiticity < 1e-10
    assert max_involution < 1e-10

    clocks = comp.link_ops(clock)
    shifts = {power: np.linalg.matrix_power(shift, power) for power in range(3)}
    link_identity = np.eye(3**4, dtype=complex)
    bianchi = clocks[0] @ clocks[1] @ clocks[2].conj().T @ clocks[3].conj().T
    projector = (link_identity + bianchi + bianchi.conj().T) / 3

    gauss_projector = np.zeros_like(link_identity)
    for vertex_shifts in __import__("itertools").product(range(3), repeat=4):
        gauss_projector += comp.gauss_transform(vertex_shifts, shifts)
    gauss_projector /= 3**4

    max_gauss_comm = 0.0
    for vertex in range(4):
        generator = comp.gauss_generator(vertex, shifts)
        max_gauss_comm = max(
            max_gauss_comm,
            float(np.linalg.norm(generator @ bianchi - bianchi @ generator)),
            float(np.linalg.norm(generator @ projector - projector @ generator)),
        )

    projector_error = float(np.linalg.norm(projector @ projector - projector))
    hermiticity_error = float(np.linalg.norm(projector - projector.conj().T))
    gauss_projector_error = float(np.linalg.norm(gauss_projector @ gauss_projector - gauss_projector))
    gauss_rank = round(float(np.real(np.trace(gauss_projector))))
    bianchi_rank = round(float(np.real(np.trace(projector))))
    physical_zero_flux_rank = round(float(np.real(np.trace(gauss_projector @ projector))))

    assert projector_error < 1e-10
    assert hermiticity_error < 1e-10
    assert max_gauss_comm < 1e-10
    assert gauss_projector_error < 1e-10
    assert gauss_rank == 3
    assert bianchi_rank == 27
    assert physical_zero_flux_rank == 1

    return {
        "status": "center_layer_constructed",
        "x_stabilizers": rows,
        "max_bare_center_residual": max_bare_residual,
        "max_dressed_center_residual": max_dressed_residual,
        "max_hermiticity_error": max_hermiticity,
        "max_involution_error": max_involution,
        "bianchi_projector": {
            "projector_error": projector_error,
            "hermiticity_error": hermiticity_error,
            "max_gauss_commutator": max_gauss_comm,
            "gauss_projector_rank": gauss_rank,
            "zero_flux_rank_before_gauss": bianchi_rank,
            "zero_flux_rank_after_gauss": physical_zero_flux_rank,
        },
        "meaning": (
            "The SU(3)-center/qutrit part of the TCH dressed X operator is explicit: "
            "dressed X operators are Hermitian involutions and the oriented four-link "
            "Bianchi projector is Gauss-invariant."
        ),
    }


def open_link_center_gap_audit() -> dict[str, Any]:
    qutrit = load_module("css_openlink_qutrit_mod", "css_openlink_qutrit_su3_tch.py")
    rows = []
    min_gap = float("inf")
    for eta in [0.0, 0.25, 0.5, 1.0]:
        for beta in [0.25, 0.5, 0.75, 1.0]:
            eigvals = qutrit.spectrum(beta=beta, kappa=1.0, eta=eta)
            levels, degeneracies = qutrit.unique_levels(eigvals)
            gap = levels[1] - levels[0]
            min_gap = min(min_gap, gap)
            rows.append(
                {
                    "beta": beta,
                    "eta": eta,
                    "gap": gap,
                    "ground_degeneracy": degeneracies[0],
                    "ground_energy": levels[0],
                }
            )
    assert min_gap > 0.1
    return {
        "status": "center_open_link_gap_survives",
        "minimum_gap_beta_le_1": min_gap,
        "rows": rows,
        "meaning": (
            "After explicit vertex Gauss projection, the Z3/qutrit open-link "
            "center sector has a finite physical gap throughout the beta<=1 scan."
        ),
    }


def peter_weyl_lock_threshold_audit() -> dict[str, Any]:
    pw = load_module("css_peter_weyl_mod", "css_peter_weyl_su3_truncation.py")
    beta_values = [0.25, 0.5, 1.0]
    lock_values = [0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0]
    rows = []
    wilson_only_min = float("inf")
    threshold_for_gap_1 = None
    threshold_for_gap_2 = None

    for strength in lock_values:
        min_gap = float("inf")
        beta_gaps = {}
        for beta in beta_values:
            hamiltonian, allowed, _, _ = pw.hamiltonian(beta, tch_lock=strength)
            eigvals = np.linalg.eigvalsh(hamiltonian)
            gap = float(eigvals[1] - eigvals[0])
            beta_gaps[str(beta)] = gap
            min_gap = min(min_gap, gap)
        if strength == 0.0:
            wilson_only_min = min_gap
        if threshold_for_gap_1 is None and min_gap >= 1.0:
            threshold_for_gap_1 = strength
        if threshold_for_gap_2 is None and min_gap >= 2.0:
            threshold_for_gap_2 = strength
        rows.append({"lock_strength": strength, "min_gap": min_gap, "gaps_by_beta": beta_gaps})

    assert wilson_only_min < 0.1
    assert threshold_for_gap_1 is not None
    assert threshold_for_gap_2 is not None
    return {
        "status": "full_su3_minimal_truncation_needs_lock",
        "wilson_only_min_gap_beta_le_1": wilson_only_min,
        "threshold_lock_for_gap_ge_1": threshold_for_gap_1,
        "threshold_lock_for_gap_ge_2": threshold_for_gap_2,
        "rows": rows,
        "meaning": (
            "The first non-center Peter-Weyl SU(3) truncation is nearly degenerate "
            "with Wilson plaquette alone. A gauge-invariant TCH/intertwiner lock of "
            "O(1) strength is needed to restore a gap comparable to the local mirror gap."
        ),
    }


def charged_mirror_survival_audit() -> dict[str, Any]:
    chain = load_module("css_spatial_mirror_chain_mod", "css_spatial_mirror_fock_gauge_gap.py")
    cases = [
        {"length": 2, "states_per_charge": 4, "beta": 0.5, "hopping": 0.2},
        {"length": 3, "states_per_charge": 4, "beta": 0.5, "hopping": 0.2},
        {"length": 4, "states_per_charge": 2, "beta": 0.5, "hopping": 0.2},
        {"length": 3, "states_per_charge": 4, "beta": 1.0, "hopping": 0.2},
    ]
    rows = []
    min_full_gap = float("inf")
    min_matter_gap = float("inf")
    for case in cases:
        hamiltonian, charge_counts, components = chain.build_chain_hamiltonian(**case)
        eigvals = np.linalg.eigvalsh(hamiltonian)
        levels, degeneracies = chain.unique_levels(eigvals)
        full_gap = levels[1] - levels[0]
        matter_gap, matter_delta = chain.matter_dominated_gap(hamiltonian, components["matter"])
        min_full_gap = min(min_full_gap, full_gap)
        min_matter_gap = min(min_matter_gap, matter_gap)
        rows.append(
            {
                **case,
                "dimension": int(hamiltonian.shape[0]),
                "charge_counts": {str(k): int(v) for k, v in charge_counts.items()},
                "full_gap": full_gap,
                "matter_dominated_gap": matter_gap,
                "matter_expectation_delta": matter_delta,
                "ground_degeneracy": degeneracies[0],
            }
        )
    assert min_full_gap > 0.1
    assert min_matter_gap > 1.0
    return {
        "status": "finite_charged_mirror_gap_survives_small_audit",
        "minimum_full_gap": min_full_gap,
        "minimum_matter_gap": min_matter_gap,
        "rows": rows,
        "meaning": (
            "In the finite charge-resolved Z3 open-link chain, the local mirror SMG "
            "gap survives small gauge-covariant hopping under exact Gauss-law constraints. "
            "This is a low-energy 1D audit, not a 3+1D continuum result."
        ),
    }


def build_verdict() -> dict[str, Any]:
    checks = {
        "center_dressed_operator": center_dressed_operator_audit(),
        "open_link_center_gap": open_link_center_gap_audit(),
        "peter_weyl_lock_threshold": peter_weyl_lock_threshold_audit(),
        "charged_mirror_survival": charged_mirror_survival_audit(),
    }
    return {
        "verdict": "center_operator_closed_full_operator_conditional",
        "summary": (
            "The center-layer TCH dressed operator can be written and passes Gauss/Bianchi "
            "checks; small Z3 open-link and charged-mirror survival audits stay gapped. "
            "The first non-center Peter-Weyl SU(3) truncation remains nearly degenerate "
            "unless an O(1) gauge-invariant TCH lock is supplied. Therefore the actual "
            "truncated-cube dressed operator is not fully closed: deriving that lock from "
            "the TCH gauge-cell geometry is the remaining finite target before any stronger "
            "mirror-gap claim."
        ),
        "checks": checks,
        "next_required_artifact": [
            "derive the Peter-Weyl/intertwiner lock from the truncated-cube gauge-cell X-plaquette operator, not by adding it by hand",
            "include at least 8 and sextet channels in the SU(3) truncation and check convergence of the strong-coupling gap",
            "couple the derived lock to the charged mirror-Fock/SMG block and repeat the Gauss-law survival audit beyond the Z3 center truncation",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    print("TCH truncated-cube dressed-operator audit")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    for name, check in verdict["checks"].items():
        print(f"  {name}: {check['status']}")
    print("\nKey obstruction:")
    print(f"  {verdict['checks']['peter_weyl_lock_threshold']['meaning']}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
