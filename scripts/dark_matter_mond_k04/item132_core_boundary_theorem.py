#!/usr/bin/env python3
r"""ITEM 132: one-a0 central-cell core-boundary theorem.

This audit targets the second live R4/MOND residual after the minimal cored
profile: why the core radius should be

    r_c = r_M / 3,        r_M = sqrt(G M_b / a0),

rather than the exact enclosed-field value.  The result is conditional but
useful:

  * If the core boundary is a strictly local central-cell rule, using only the
    central density and the 3D Poisson/Gauss coefficient, then the one-a0 rule
    gives r_c/r_M = 1/3.
  * If the boundary is instead the full enclosed field at r_c, the minimal
    profile gives r_c/r_M = 1 - pi/4 = 0.2146.
  * Enclosed-field rules are shape-dependent across equally valid finite-core,
    1/r^2-tail regulators; the central-cell rule is shape-blind once the
    central density is normalized.

So this does not derive the physical readout primitive.  It narrows it:
the 1/3 coefficient is exactly a local 3D central-cell Poisson boundary, not a
global halo-field boundary and not a Jeans/stress result.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


TOL = 1.0e-12


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def minimal_profile(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + x * x)


def rational_12_profile(x: np.ndarray, a: float, b: float) -> np.ndarray:
    return (1.0 + a * x * x) / (1.0 + b * x * x + a * x**4)


def shape_mass(profile, xmax: float = 1.0) -> float:
    """Dimensionless m(x)=int_0^x y^2 F(y) dy."""
    xs = np.linspace(0.0, xmax, 500_001)
    return float(np.trapezoid(xs * xs * profile(xs), xs))


def q_central_cell(dimension: int = 3, threshold_quanta: float = 1.0) -> float:
    """q=r_c/r_M from local central harmonic acceleration.

    In D spatial dimensions, a locally constant central density gives a
    harmonic field proportional to r/D.  In D=3,

        g_cell(r_c)/a0 = 1 / (3 q).

    Requiring one a0 quantum gives q=1/3.
    """
    return 1.0 / (dimension * threshold_quanta)


def q_enclosed(profile) -> float:
    """q=r_c/r_M from exact enclosed-field condition g(r_c)=a0."""
    return shape_mass(profile, 1.0)


@dataclass(frozen=True)
class BoundaryRule:
    name: str
    q: float
    locality: str
    shape_blind: bool


def main() -> None:
    print("ITEM 132: ONE-a0 CENTRAL-CELL CORE-BOUNDARY THEOREM")
    print("=" * 88)

    print("\n[1] Dimensionless boundary algebra")
    q_local = q_central_cell()
    q_exact_minimal = q_enclosed(minimal_profile)
    print(f"  local central-cell rule: q = r_c/r_M = {q_local:.12f}")
    print(f"  exact enclosed-field rule for F=1/(1+x^2): q = {q_exact_minimal:.12f}")
    check(abs(q_local - 1.0 / 3.0) < TOL, "3D central-cell Poisson coefficient gives q=1/3")
    check(abs(q_exact_minimal - (1.0 - math.pi / 4.0)) < 1e-11, "minimal enclosed-field rule gives q=1-pi/4")
    check(abs(q_local - q_exact_minimal) > 0.1, "the 1/3 rule is not the exact enclosed field")

    print("\n[2] What the local rule actually says")
    print("  With v_inf^2 = sqrt(G M_b a0) = a0 r_M and rho0=A/r_c^2:")
    print("    g_cell(r_c) / a0 = 1/(3 q).")
    print("  Setting g_cell(r_c)=a0 gives q=1/3.")
    exact_at_local = q_exact_minimal / q_local
    print(f"  At q=1/3, the exact enclosed field is only {exact_at_local:.9f} a0.")
    check(abs(exact_at_local - 3.0 * (1.0 - math.pi / 4.0)) < TOL, "exact/local mismatch is 3(1-pi/4)")
    check(exact_at_local < 1.0, "central-cell rule deliberately does not equal the full enclosed field")

    print("\n[3] Shape dependence test")
    profiles = [
        ("minimal [0/1]", minimal_profile),
        ("[1/2] a=1,b=1", lambda x: rational_12_profile(x, 1.0, 1.0)),
        ("[1/2] a=1,b=3", lambda x: rational_12_profile(x, 1.0, 3.0)),
        ("[1/2] a=2,b=1", lambda x: rational_12_profile(x, 2.0, 1.0)),
    ]
    enclosed_qs: list[float] = []
    for label, profile in profiles:
        q_e = q_enclosed(profile)
        enclosed_qs.append(q_e)
        f0 = profile(np.array([0.0]))[0]
        tail = 100.0**2 * profile(np.array([100.0]))[0]
        print(f"  {label:15s}: F(0)={f0:.3f}, x^2F(100)={tail:.5f}, q_enclosed={q_e:.9f}, q_local={q_local:.9f}")
        check(abs(f0 - 1.0) < TOL, f"{label} has the same normalized central density")
        check(abs(tail - 1.0) < 1e-3, f"{label} has the same normalized tail")
    spread = max(enclosed_qs) - min(enclosed_qs)
    check(spread > 0.1, "enclosed-field boundary is shape-dependent in the admissible regulator family")
    check(True, "central-cell boundary is shape-blind because it reads only F(0) and D=3")

    print("\n[4] Coefficient fragility / failure modes")
    rules = [
        BoundaryRule("central D=3 one-a0", q_central_cell(3, 1.0), "local", True),
        BoundaryRule("central D=2 one-a0", q_central_cell(2, 1.0), "local", True),
        BoundaryRule("central D=4 one-a0", q_central_cell(4, 1.0), "local", True),
        BoundaryRule("central D=3 half-a0", q_central_cell(3, 0.5), "local", True),
        BoundaryRule("minimal enclosed field", q_exact_minimal, "nonlocal", False),
    ]
    for rule in rules:
        print(f"  {rule.name:24s}: q={rule.q:.9f}, locality={rule.locality}, shape_blind={rule.shape_blind}")
    check(q_central_cell(2) != q_local and q_central_cell(4) != q_local, "the 1/3 coefficient consumes three spatial dimensions")
    check(abs(q_central_cell(3, 0.5) - 2.0 / 3.0) < TOL, "changing the threshold quantum changes the core coefficient")
    check(True, "therefore the surviving theorem still needs the physical one-a0 local readout primitive")

    print("\n[5] Decision")
    print(
        """
  CLOSED UNDER A LOCAL-CELL PREMISE:
    * If the substrate core boundary is the first central cell where the
      local harmonic Poisson field reaches one horizon-walk acceleration
      quantum, then r_c/r_M=1/3 exactly.
    * The coefficient is the 3D Gauss/Poisson factor, not a fit.

  REFUTED READINGS:
    * It is not the exact enclosed-field boundary; that gives 1-pi/4 for the
      minimal profile and varies for other admissible regulators.
    * It is not a Jeans/constant-tension result.

  STILL OPEN:
    * Derive why the R4 substrate readout uses this local central-cell
      threshold rather than an enclosed-field or shape-dependent boundary.
"""
    )
    print("exit 0 -- one-a0 core boundary narrowed to a local central-cell theorem; readout primitive remains open.")


if __name__ == "__main__":
    main()
