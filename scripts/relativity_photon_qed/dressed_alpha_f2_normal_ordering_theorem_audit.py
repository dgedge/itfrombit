#!/usr/bin/env python3
r"""Dressed alpha: Ward-preserving finite F^2 normal-ordering audit.

Question
--------
Can the remaining dressed-alpha residual be closed by a finite normal-ordering
or matching theorem for the substrate Peierls Maxwell term?

The desired theorem would keep the Peierls/Ward identity exact while adding a
finite local counterterm to the electromagnetic Hessian.  In momentum space
the only local gauge-invariant quadratic term has the Maxwell tensor

    T_ij(q) = |q|^2 delta_ij - q_i q_j,

so the most general finite local repair is

    K_ij(q) -> K_ij(q) + C_fin T_ij(q).

Result
------
This script proves the theorem-shaped negative under the current canon:
Ward transversality and pure-gauge no-response are true for every C_fin.  They
license a finite F^2 term but do not select its value.  The value needed to
turn the current web-dressed Ward/Kubo Peierls slot into the observed
alpha^{-1} residual is large and specific:

    Z_fin = 2.537609,  q_eff = 1.592988,
    additive C_fin = 0.605928 of the observed residual.

No existing Peierls, service-projector, or unit-charge condition chooses those
numbers.  Therefore a Ward-preserving finite F^2 theorem can be formulated
only as a conditional theorem with a new independent substrate selector; it is
not derivable from the present axioms.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_ward_kubo_peierls_observable_audit as wk


ALPHA0_INV = 137.0
ALPHA_PHYS_INV = 137.035999084
DELTA_TARGET = ALPHA_PHYS_INV - ALPHA0_INV


@dataclass(frozen=True)
class SelectorTest:
    name: str
    fixes_c: bool
    verdict: str


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def transverse_f2_tensor(q: np.ndarray) -> np.ndarray:
    """The unique isotropic local Maxwell Hessian at quadratic order."""

    q = np.asarray(q, dtype=float)
    return float(q @ q) * np.eye(q.size) - np.outer(q, q)


def ward_defect(kernel: np.ndarray, q: np.ndarray) -> float:
    """Dimensionless two-sided Ward defect for a symmetric kernel."""

    denom = max(float(np.linalg.norm(kernel) * np.linalg.norm(q)), 1.0e-30)
    return float((np.linalg.norm(kernel @ q) + np.linalg.norm(q @ kernel)) / denom)


def pure_gauge_response(kernel: np.ndarray, q: np.ndarray) -> float:
    """Quadratic response to a pure gauge A_i proportional to q_i."""

    return float(q @ kernel @ q)


def finite_f2_family(q: np.ndarray, ward_ratio: float, c_fin: float) -> np.ndarray:
    """Toy local Hessian family with the same Ward tensor for all finite parts."""

    tensor = transverse_f2_tensor(q)
    return (ward_ratio + c_fin) * tensor


def selector_tests() -> list[SelectorTest]:
    return [
        SelectorTest(
            "Ward transversality",
            False,
            "requires q_i K_ij = 0; every local F^2 finite term already obeys it",
        ),
        SelectorTest(
            "pure-gauge no-response",
            False,
            "requires K_ij q_j = 0; again true for every finite F^2 coefficient",
        ),
        SelectorTest(
            "normal-order vacuum subtraction",
            False,
            "removes tadpoles/constant pieces but a local transverse F^2 term vanishes on constants",
        ),
        SelectorTest(
            "unit Peierls charge",
            True,
            "fixes the gauge derivative to q_eff = 1, which is not the required q_eff",
        ),
        SelectorTest(
            "minimal finite subtraction",
            True,
            "sets C_fin = 0 by convention and leaves the Ward/Kubo slot at the undershooting value",
        ),
        SelectorTest(
            "match observed alpha",
            True,
            "sets the needed coefficient but is the fit we are trying to avoid",
        ),
    ]


def main() -> None:
    print("DRESSED ALPHA: WARD-PRESERVING FINITE F^2 NORMAL-ORDERING AUDIT")
    print("=" * 100)

    ward = wk.ward_kubo_dispersion()
    ward_ratio = float(ward["ratio_web"])
    missing_fraction = 1.0 - ward_ratio
    z_fin = 1.0 / ward_ratio
    q_eff = math.sqrt(z_fin)
    missing_delta = DELTA_TARGET * missing_fraction

    print("[1] Current Ward/Kubo residual")
    print(f"  target delta(alpha^-1)             = {DELTA_TARGET:.9f}")
    print(f"  web-dressed Ward/Kubo delta/target = {ward_ratio:.6f}")
    print(f"  missing finite F^2 fraction        = {missing_fraction:.6f}")
    print(f"  missing finite delta(alpha^-1)     = {missing_delta:.9f}")
    print(f"  multiplicative Z_fin equivalent    = {z_fin:.6f}")
    print(f"  Peierls charge equivalent q_eff    = {q_eff:.6f}")
    check(0.35 < ward_ratio < 0.45, "physical Ward/Kubo slot undershoots the observed residual")
    check(0.55 < missing_fraction < 0.65, "additive F^2 repair is a large finite coefficient")
    check(abs(q_eff - 1.0) > 0.5, "multiplicative repair is not the unit Peierls charge")

    print("\n[2] Ward-preserving F^2 family")
    q_vectors = (
        np.array([1.0, 0.0, 0.0]),
        np.array([1.0, 2.0, -1.0]),
        np.array([0.25, -0.5, 0.75]),
    )
    c_values = (-2.0, -0.25, 0.0, missing_fraction, 1.0, 2.537609)
    max_ward = 0.0
    max_pure = 0.0
    for q in q_vectors:
        for c_fin in c_values:
            kernel = finite_f2_family(q, ward_ratio, c_fin)
            max_ward = max(max_ward, ward_defect(kernel, q))
            max_pure = max(max_pure, abs(pure_gauge_response(kernel, q)))
    print(f"  max Ward defect over finite family       = {max_ward:.3e}")
    print(f"  max pure-gauge quadratic response        = {max_pure:.3e}")
    print("  The finite local F^2 coefficient is invisible to these gauge tests.")
    check(max_ward < 1e-12, "every finite local F^2 term preserves the Ward identity")
    check(max_pure < 1e-12, "every finite local F^2 term gives zero pure-gauge response")

    print("\n[3] Candidate selector inventory")
    selected = [row for row in selector_tests() if row.fixes_c]
    noncircular = [
        row
        for row in selected
        if row.name not in {"unit Peierls charge", "minimal finite subtraction", "match observed alpha"}
    ]
    for row in selector_tests():
        status = "fixes a convention" if row.fixes_c else "leaves C_fin free"
        print(f"  {row.name:<34} {status:<22} {row.verdict}")
    check(not noncircular, "no current non-circular substrate selector fixes the finite coefficient")

    print("\nTHEOREM ATTEMPT")
    print(
        r"""
  Conditional theorem:

    Let K_ij(q) be the Peierls electromagnetic Hessian produced by the bridge
    and Wilson web.  Ward exactness requires q_i K_ij(q)=0.  At quadratic
    order in small q, the only local isotropic gauge-invariant finite
    counterterm is C_fin (|q|^2 delta_ij - q_i q_j), i.e. a finite F^2
    normal-ordering term.  Therefore Ward exactness fixes the tensor structure
    but leaves C_fin arbitrary.

    If a new substrate rule independently selected

        C_fin / delta_obs = 0.605928

    or equivalently Z_fin = 2.537609 in the Ward/Kubo slot, then the
    dressed-alpha residual would close without the retired N1=31 count.

  No-selector corollary:

    Under the present canon this theorem is not predictive.  Ward identity,
    pure-gauge decoupling, normal-ordering of constants, Peierls unit charge,
    and the service-projector algebra do not select the required finite
    coefficient.  Matching alpha would select it only circularly.

  Verdict:

    A Ward-preserving finite F^2 theorem is structurally allowed but not
    derivable yet.  The dressed 137.035999... value therefore remains a bounded
    electromagnetic residual unless a genuinely new EM normal-ordering or
    sector-billing primitive is added and justified independently of alpha.
exit 0 -- finite F^2 matching reduces to one unselected Ward-preserving scalar.
"""
    )


if __name__ == "__main__":
    main()
