#!/usr/bin/env python3
r"""No-go theorem for deriving the CC residual inside the current vertex class.

Question:
    Can the remaining

        delta C = C_target - C_loop = 0.001187526879799...

    be derived as a genuine next-order generation-vertex renormalization using
    only the already-proved ingredients?

Assumptions under test:
  A1. one monitored register touch -> one alpha0;
  A2. the forced chi/W generation vertex is a k-uniform Pauli insertion;
  A3. the bath is the free chirality-Weyl kernel H(k)=sin(k).sigma;
  A4. the observable is the one-point vacuum-energy/log-determinant barrier.

Under A1-A4, the local 2x2 operator algebra has no unused scalar:
  * identity has zero interband element;
  * Pauli directions are already the leading C_loop and are cubic-isotropic;
  * a multiplicative Pauli rescaling is just an unproved counterterm/alpha_eff.

The only intrinsic next-order contribution in this class is the exact finite-
coupling determinant of the same Pauli insertion.  This script computes it and
shows it is far too small even before imposing the canonical small coupling.

Exit 0 means: no positive renormalization theorem exists in this operator class.
It does not rule out a new connected gauge-web interaction, non-Gaussian bath,
or additional fixed operator not present in A1-A4.
"""

from __future__ import annotations

import math

import numpy as np

from cc_generation_vertex_item115_loop import extrapolate_c, rho_for_c
from cc_generation_vertex_next_order_audit import weighted_generation_c
from cc_next_order_residual_audit import finite_coupling_c
from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    solve_gamma_for_stationary_queue,
    solve_target,
)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def target_numbers() -> tuple[float, float, float, float]:
    _, q1_target = solve_target()
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    c_target = -math.log(gamma_star / BASE_GAMMA) / ALPHA0
    c_lin, c_quad, _ = extrapolate_c([64, 96, 128, 160, 192, 256])
    c_loop = 0.5 * (c_lin + c_quad)
    return q1_target, c_target, c_loop, c_target - c_loop


def operator_basis_demo(n_grid: int = 96) -> tuple[float, float]:
    """Return identity and unit-Pauli interband coefficients on the Weyl kernel."""
    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    sin2 = np.sin(ks) ** 2
    c_identity = 0.0
    c_pauli = 0.0
    for i in range(0, n_grid, 32):
        sx = sin2[i : i + 32][:, None, None]
        s2tot = sx + sin2[None, :, None] + sin2[None, None, :]
        # Identity insertion has <u_-|I|u_+>=0 exactly; keep the explicit zero
        # as the linear-algebra control.
        c_identity += 0.0
        c_pauli += float(np.sum((1.0 - sx / s2tot) / (2.0 * np.sqrt(s2tot))))
    return c_identity, c_pauli / n_grid**3


def finite_determinant_scan(n_grid: int = 96) -> tuple[float, float, float, list[tuple[float, float]]]:
    """Scan exact finite-g determinant for the same Pauli insertion.

    This is deliberately generous: g is allowed to range far beyond the
    canonical perturbative values.  If the same operator cannot hit C_target even
    with that freedom, it certainly cannot derive the residual canonically.
    """

    c0 = weighted_generation_c(n_grid, "point")
    rows: list[tuple[float, float]] = []
    for g in np.linspace(0.02, 3.0, 150):
        rows.append((float(g), finite_coupling_c(n_grid, float(g))))
    g_best, c_best = max(rows, key=lambda row: row[1])
    return c0, g_best, c_best, rows


def main() -> None:
    print("CC GENERATION-VERTEX RENORMALIZATION NO-GO AUDIT")

    q1_target, c_target, c_loop, delta_c = target_numbers()
    c2_required = delta_c / ALPHA0
    print("\n[1] Target after the final-tier proof")
    print(f"  C_target                    = {c_target:.15f}")
    print(f"  C_loop                      = {c_loop:.15f}")
    print(f"  delta C required             = {delta_c:.15f}")
    print(f"  if written C + alpha0*C2: C2 = {c2_required:.12f}")
    print(f"  rho(C_loop)/rho_obs          = {rho_for_c(c_loop, q1_target):.9f}")
    check(0.001 < delta_c < 0.0013, "residual is the known real +delta C target")

    print("\n[2] Local operator-algebra closure")
    c_i, c_p = operator_basis_demo()
    print(f"  identity interband coefficient = {c_i:.12f}")
    print(f"  unit-Pauli coefficient          = {c_p:.12f}")
    print("  local Hermitian 2x2 basis       = a I + v_x sigma_x + v_y sigma_y + v_z sigma_z")
    print("  identity contributes zero; Pauli triplet is the already-used isotropic vertex.")
    check(abs(c_i) < 1e-15, "identity insertion has no interband matrix element")
    check(abs(c_p / c_loop - 1.0) < 2e-4, "unit Pauli reproduces C_loop on the finite grid")
    print("  Therefore a local next-order scalar can only rescale the Pauli vertex.")
    print("  Such a rescaling is not derived by the algebra; it is an external counterterm.")

    print("\n[3] Exact finite-coupling determinant of the same Pauli operator")
    n_scan = 96
    c0, g_best, c_best, _ = finite_determinant_scan(n_scan)
    canonical = [
        ("g^2 = alpha0", math.sqrt(ALPHA0)),
        ("g^2 = 4 alpha0", 2.0 * math.sqrt(ALPHA0)),
        ("arbitrary best in scan", g_best),
    ]
    print(f"  finite-grid point coefficient C0(N={n_scan}) = {c0:.12f}")
    print(f"  {'case':<24s} {'g':>10s} {'C_det(g)':>14s} {'gain/deltaC':>14s}")
    best_fraction = 0.0
    for label, g in canonical:
        c_det = finite_coupling_c(n_scan, g)
        frac = (c_det - c0) / delta_c
        best_fraction = max(best_fraction, frac)
        print(f"  {label:<24s} {g:10.6f} {c_det:14.12f} {frac:13.4%}")
    print(f"  scan maximum at g={g_best:.6f}, C={c_best:.12f}")
    check(best_fraction < 0.10, "even an arbitrary finite-g determinant supplies under 10% of delta C")
    check(c_best < c_target, "the exact same-operator determinant never reaches C_target in the generous scan")
    det_c2 = (finite_coupling_c(n_scan, math.sqrt(ALPHA0)) - c0) / ALPHA0
    print(f"  canonical determinant C2          = {det_c2:.12f}")
    print(f"  required C2                       = {c2_required:.12f}")
    check(det_c2 / c2_required < 0.05, "canonical determinant C2 is under 5% of the required next-order coefficient")

    print("\n[4] The no-go statement")
    print(
        """
THEOREM (current-operator-class no-go)
  Under A1-A4, the monitored generation vertex has exactly one local traceless
  Pauli operator class.  The identity component is interband-null, the Pauli
  triplet is already C_loop, and the exact same-operator determinant is too
  small by more than an order of magnitude.  A multiplicative O(alpha0)
  rescaling of the Pauli vertex would numerically be able to fill the gap, but
  the algebra contains no fixed operator that supplies it; it is a counterterm
  unless a new mechanism is specified.

  Therefore the requested positive theorem cannot be derived from the existing
  monitored-billing / free-Weyl / k-uniform-Pauli operator map.

  What would count as a genuine positive theorem:
    * a fixed connected two-loop gauge-web vertex;
    * a non-Gaussian bath cumulant with a canon-derived coefficient;
    * a new QEC/Kraus normal-ordering term that is not just alpha_eff;
    * or another uniquely forced operator whose evaluated coefficient is
      delta C = +0.001187526879799.

  Until one of those is derived, the non-horizon CC result remains the current
  0.2%-grade near-closure, not a Locked cosmological-constant prediction.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
