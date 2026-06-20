#!/usr/bin/env python3
r"""Intrinsic-gravity gate: RT + first law + service-current T_{\mu\nu}.

Question
--------
Can the framework's verified RT/min-cut entanglement, the entanglement first
law, and an explicit QEC service-current stress tensor produce the linearized
Einstein equation?

This script is deliberately conservative.  It does not claim to derive the
observed Planck hierarchy, and it does not claim an exact all-scale discrete
Einstein equation.  It checks the finite ingredients needed for the standard
Jacobson/Faulkner linearized-gravity route:

  1. RT/min-cut first variations span the substrate edge perturbations.
  2. The service ledger supplies the boost modular Hamiltonian.
  3. The service ledger supplies a symmetric, positive, conserved T^{mu nu}.
  4. All-ball/all-interval first-law constraints localize the residual.
  5. The coefficient identity 2 pi / eta = 8 pi G follows from eta=1/(4G).

Exit 0 means:
  * the Einstein/source FORM is intrinsic at leading-order continuum grade;
  * the observed M_Pl hierarchy remains open and nonlocal/horizon-scale.
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PYTHON_CODE = ROOT / "python_code"

HBARC = 0.197327
A0_FM = 0.594
LAMBDA_QCD = HBARC / A0_FM
S_CELL = 55.0 / 8.0
M_PL_OBS = 1.2209e19


def run_owner(script: str, required_phrases: list[str]) -> str:
    """Run an owner audit and require its key verdict phrases."""

    proc = subprocess.run(
        [sys.executable, str(PYTHON_CODE / script)],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    out = proc.stdout
    missing = [phrase for phrase in required_phrases if phrase not in out]
    if missing:
        print(f"FAIL: {script} did not print required phrases:")
        for phrase in missing:
            print(f"  - {phrase}")
        raise SystemExit(1)
    print(f"PASS: {script}")
    for phrase in required_phrases:
        print(f"      found: {phrase}")
    return out


def interval_modular_kernel(n: int = 25) -> np.ndarray:
    """Build a 1D all-interval modular-kernel matrix.

    In a CFT vacuum interval, the modular Hamiltonian kernel is proportional
    to (R^2-(x-x0)^2)/(2R) inside the interval.  This toy does not attempt to
    prove the CFT theorem; it checks the finite linear-algebra fact needed by
    the all-balls argument: the family of local interval kernels separates all
    pointwise residuals on the lattice.
    """

    xs = np.arange(n, dtype=float)
    rows: list[np.ndarray] = []
    # Put interval endpoints on the half-grid.  If endpoints are placed exactly
    # on lattice sites, the two outer sites are always endpoint-zeroes and the
    # finite matrix artificially loses rank.
    endpoints = np.arange(-0.5, n + 0.5, 1.0)
    for left in endpoints:
        for right in endpoints:
            if right - left < 2.0:
                continue
            center = 0.5 * (left + right)
            radius = 0.5 * (right - left)
            row = np.zeros(n, dtype=float)
            for k in range(n):
                if not (left < xs[k] < right):
                    continue
                value = (radius * radius - (xs[k] - center) ** 2) / (2.0 * radius)
                if value > 1.0e-14:
                    row[k] = value
            if np.linalg.norm(row) > 0.0:
                rows.append(row)
    return np.vstack(rows)


def jacobson_coefficient() -> float:
    """Return the Einstein prefactor in units where G=1."""

    eta = 1.0 / 4.0
    return 2.0 * math.pi / eta


def emergent_bare_planck_mass() -> float:
    """M_Pl,bare from 1/(4G)=s_cell/a0^2."""

    return 2.0 * math.sqrt(S_CELL) * LAMBDA_QCD


def main() -> None:
    print("INTRINSIC GRAVITY / LINEARIZED EINSTEIN GATE")
    print("=" * 88)

    print("\n[1] Owner audits")
    run_owner(
        "rt_first_law_linearized_gravity_audit.py",
        [
            "RT subgradient-family combined rank: 80/80",
            "first-law constraint operator verified",
            "ALL ASSERTIONS PASSED",
        ],
    )
    run_owner(
        "modular_hamiltonian_service_ledger.py",
        [
            "STATIC-CODE OBSTRUCTION",
            "H_mod = 2pi",
            "ALL ASSERTIONS PASSED",
        ],
    )
    run_owner(
        "gravity_service_current_stress_tensor_gate.py",
        [
            "SERVICE-CURRENT STRESS-TENSOR GATE",
            "T^{mu nu}_svc",
            "ALL ASSERTIONS PASSED",
        ],
    )

    print("\n[2] All-interval localization toy")
    kernel = interval_modular_kernel(n=25)
    rank = int(np.linalg.matrix_rank(kernel, tol=1.0e-12))
    print(f"    interval rows = {kernel.shape[0]}, lattice points = {kernel.shape[1]}")
    print(f"    rank(all interval kernels) = {rank}/{kernel.shape[1]}")
    assert rank == kernel.shape[1]
    print("    -> if the first-law residual vanishes for all local intervals,")
    print("       the pointwise residual is zero in this finite model.")

    print("\n[3] Jacobson/Faulkner coefficient")
    prefactor = jacobson_coefficient()
    print(f"    2 pi / eta = {prefactor:.12f} G, eta=1/(4G)")
    print(f"    8 pi G     = {8.0 * math.pi:.12f} G")
    assert abs(prefactor - 8.0 * math.pi) < 1.0e-12
    print("    -> the entanglement first law with area density 1/(4G)")
    print("       gives the standard Einstein coupling.")

    print("\n[4] Hierarchy boundary")
    m_pl_bare = emergent_bare_planck_mass()
    hierarchy = M_PL_OBS / m_pl_bare
    enhancement = hierarchy * hierarchy
    print(f"    M_Pl,bare = 2 sqrt(55/8) Lambda_QCD = {m_pl_bare:.6g} GeV")
    print(f"    M_Pl,obs / M_Pl,bare = {hierarchy:.6e}")
    print(f"    required 1/G enhancement = {enhancement:.6e}")
    assert 1.0e18 < hierarchy < 1.0e20
    assert enhancement > 1.0e37
    print("    -> no local RT/service-current ingredient supplies this factor.")

    print(
        r"""
[5] VERDICT
    CLOSED AT LEADING-ORDER / CONTINUUM-JACOBSON GRADE:
      RT first variations + boost modular Hamiltonian + service-current
      T_{\mu\nu} give the linearized Einstein/source structure.  The framework
      no longer has to import a generic matter stress tensor; it has one from
      the service ledger, and the all-region first-law constraint has the right
      localizing rank.

    NOT CLOSED:
      the observed Planck hierarchy.  The intrinsic area density gives a
      lattice/bare M_Pl of order GeV, while the measured Planck mass requires
      a nonlocal/horizon-scale enhancement of order 10^38 in 1/G.  Exact
      discrete modular flow, gravitational dressing, and the HBC/dS continuum
      lift also remain theorem targets.

    Therefore the precise answer is:

      Einstein equation FORM:     intrinsic, conditional-closed.
      Observed G / M_Pl hierarchy: still horizon/nonlocal-theorem open.
"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
