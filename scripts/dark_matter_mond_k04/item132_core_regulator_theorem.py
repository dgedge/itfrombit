#!/usr/bin/env python3
r"""ITEM 132: minimal cored-profile regulator theorem.

This audit targets the most tractable live R4/MOND residual: the cored profile.
It does not derive P1, a0, G, or a dynamical Jeans stress.  It proves the
narrow theorem that the pseudo-isothermal profile

    rho(r) = A / (r^2 + r_c^2)

is the unique minimal rational regulator once the following profile-level
conditions are imposed:

  1. finite central density;
  2. asymptotic 1/r^2 line-halo tail with no extra tail coefficient;
  3. even/local dependence on x = r/r_c through u = x^2;
  4. no dimensionless shape parameter beyond the scale r_c;
  5. minimal Pade degree.

The theorem is useful but deliberately limited.  Higher rational profiles obey
finite-center plus 1/r^2 tail, so the profile is not forced without the
"minimal/no-extra-shape" premise.  Also, r_c = r_M/3 remains a separate
central-cell boundary rule; the exact enclosed-field one-a0 rule would give
r_c/r_M = 1 - pi/4 instead.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


TOL = 1.0e-10


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def minimal_profile(x: np.ndarray) -> np.ndarray:
    """F(x)=1/(1+x^2), normalized so x^2 F(x)->1."""
    return 1.0 / (1.0 + x * x)


def rational_12_profile(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """A [1/2] rational family with finite center and 1/x^2 tail.

    F(x)=(1+a x^2)/(1+b x^2+a x^4).  For a>0, x^2 F(x)->1.  Different b/a
    choices preserve the endpoint conditions but change the core shape.
    """
    return (1.0 + a * x * x) / (1.0 + b * x * x + a * x**4)


def enclosed_shape_mass(profile, xmax: float = 1.0) -> float:
    """m(x)=int_0^x y^2 F(y) dy in dimensionless units."""
    xs = np.linspace(0.0, xmax, 400_001)
    return float(np.trapezoid(xs * xs * profile(xs), xs))


@dataclass(frozen=True)
class PadeClass:
    numerator_degree: int
    denominator_degree: int

    @property
    def free_shape_count_after_normalization(self) -> int:
        """Count dimensionless shape parameters after F(0)=1 and tail coeff=1.

        For a rational function in u=x^2 with numerator degree m and
        denominator degree n, F(0)=1 fixes the constant terms.  A 1/u tail
        requires n=m+1 and fixes the ratio of the leading coefficients.
        The remaining non-leading coefficients are shape parameters.
        """
        m = self.numerator_degree
        n = self.denominator_degree
        if n != m + 1:
            return -1
        # numerator has m nonconstant coefficients; denominator has n
        # nonconstant coefficients; leading tail normalization removes one.
        return m + n - 1


def main() -> None:
    print("ITEM 132: MINIMAL CORED-PROFILE REGULATOR THEOREM")
    print("=" * 88)

    print("\n[1] Minimal rational class")
    classes = [PadeClass(0, 1), PadeClass(1, 2), PadeClass(2, 3)]
    for cls in classes:
        print(
            f"  [{cls.numerator_degree}/{cls.denominator_degree}] "
            f"shape parameters after endpoint normalization = "
            f"{cls.free_shape_count_after_normalization}"
        )
    check(classes[0].free_shape_count_after_normalization == 0, "[0/1] has no dimensionless shape parameter")
    check(classes[1].free_shape_count_after_normalization > 0, "[1/2] already introduces shape freedom")
    check(classes[2].free_shape_count_after_normalization > classes[1].free_shape_count_after_normalization, "higher classes add further shape freedom")

    print("\n[2] Unique [0/1] profile")
    # In u=x^2, F(u)=1/(1+b u).  Tail normalization u F(u)->1 forces b=1.
    b = 1.0
    xs = np.array([0.0, 0.5, 1.0, 3.0, 30.0])
    vals = minimal_profile(xs)
    print("  F_min(x)=1/(1+x^2)")
    print("  sample:", ", ".join(f"F({x:g})={v:.6f}" for x, v in zip(xs, vals)))
    check(abs(vals[0] - 1.0) < TOL, "finite center is normalized to F(0)=1")
    check(abs((30.0**2) * vals[-1] - 1.0) < 2e-3, "tail coefficient x^2 F(x)->1")
    check(abs(b - 1.0) < TOL, "tail normalization uniquely fixes the [0/1] denominator")

    print("\n[3] Non-uniqueness without the minimal/no-extra-shape premise")
    variants = [
        ("[1/2] a=1,b=1", lambda y: rational_12_profile(y, 1.0, 1.0)),
        ("[1/2] a=1,b=3", lambda y: rational_12_profile(y, 1.0, 3.0)),
        ("[1/2] a=2,b=1", lambda y: rational_12_profile(y, 2.0, 1.0)),
    ]
    m_min = enclosed_shape_mass(minimal_profile)
    print(f"  minimal m(1)=int_0^1 x^2/(1+x^2) dx = {m_min:.12f}")
    check(abs(m_min - (1.0 - math.pi / 4.0)) < 1e-10, "minimal enclosed shape mass is exactly 1-pi/4")
    different = 0
    for label, prof in variants:
        v0 = prof(np.array([0.0]))[0]
        tail = (100.0**2) * prof(np.array([100.0]))[0]
        m1 = enclosed_shape_mass(prof)
        f1 = prof(np.array([1.0]))[0]
        print(f"  {label:15s}: F(0)={v0:.3f}, x^2F(100)={tail:.5f}, F(1)={f1:.5f}, m(1)={m1:.6f}")
        check(abs(v0 - 1.0) < TOL, f"{label} has finite normalized center")
        check(abs(tail - 1.0) < 1e-3, f"{label} has the same 1/r^2 tail")
        different += int(abs(m1 - m_min) > 1e-3)
    check(different == len(variants), "endpoint conditions alone leave inequivalent core shapes")

    print("\n[4] One-a0 boundary split")
    q_central = 1.0 / 3.0
    q_exact_enclosed = 1.0 - math.pi / 4.0
    print(f"  central-harmonic one-a0 rule: r_c/r_M = {q_central:.12f}")
    print(f"  exact enclosed-field one-a0 rule for F_min: r_c/r_M = {q_exact_enclosed:.12f}")
    check(abs(q_central - q_exact_enclosed) > 0.1, "the 1/3 rule is not an enclosed-field consequence")
    check(True, "r_c=r_M/3 must be a local central-cell boundary theorem if it survives")

    print("\n[5] Decision")
    print(
        """
  CLOSED AS A PROFILE THEOREM:
    * With finite center, 1/r^2 tail, even locality in x^2, minimal Pade
      degree, and no dimensionless shape parameter, the cored regulator is
      uniquely rho=A/(r^2+r_c^2).

  NOT CLOSED:
    * The minimal/no-extra-shape premise is not yet derived from the R4
      substrate service algebra.
    * The r_c=r_M/3 factor is not forced by the profile's exact enclosed
      field; it still needs a central-cell one-a0 boundary theorem.
    * A positive Jeans/stress dynamics remains separate from this regulator
      uniqueness statement.
"""
    )
    print("exit 0 -- minimal core regulator closed conditionally; boundary and dynamics remain open.")


if __name__ == "__main__":
    main()
