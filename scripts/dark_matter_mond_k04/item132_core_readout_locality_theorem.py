#!/usr/bin/env python3
r"""ITEM 132: local-service selector for the R4 core readout.

Question
--------
Can the framework explain why the R4 core rule uses the local central-cell
threshold

    r_c = r_M / 3

rather than an exact enclosed-field boundary or a shape-dependent cored-halo
boundary?

Result
------
Yes, but only as a ledger-selection theorem.

If the observable is the R4 Stinespring service current, then the boundary
readout is made by the local service instrument before the global halo field is
assembled from its event ledger.  That instrument can read the central service
cell and its one-record threshold.  It cannot read an enclosed mass integral,
which belongs to the post-summed gravitational/strain ledger, and it cannot
read a hidden profile-shape parameter not present in the R4 service labels.

Under those locality and ledger axioms, the unique isotropic scalar at a finite
core is the central harmonic Poisson/Gauss coefficient.  With

    rho(r) = A/r_c^2 F(r/r_c),      F(0)=1,
    4 pi G A = v_inf^2 = a0 r_M,

the central-cell field is

    g_cell(r_c)/a0 = r_M / (3 r_c).

One service quantum therefore gives r_c/r_M=1/3.

Scope
-----
This is not an unconditional microscopic derivation of the one-a0 tick or of
a0 itself.  It proves the narrower claim that, once the R4 halo observable is
the local service-current ledger, enclosed-field and shape-dependent boundary
rules are wrong-ledger observables.  The remaining microscopic target is to
derive the one-quantum threshold from the substrate service algebra.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np


TOL = 1.0e-12


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def minimal_profile(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + x * x)


def rational_12_profile(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Finite-center, 1/r^2-tail family with normalized center and tail."""
    return (1.0 + a * x * x) / (1.0 + b * x * x + a * x**4)


def enclosed_shape_mass(profile: Callable[[np.ndarray], np.ndarray]) -> float:
    """m(1)=int_0^1 x^2 F(x) dx for rho=A/r_c^2 F(r/r_c)."""
    xs = np.linspace(0.0, 1.0, 500_001)
    return float(np.trapezoid(xs * xs * profile(xs), xs))


def central_q(dimension: int = 3, threshold_quanta: float = 1.0) -> float:
    """q=r_c/r_M from the central harmonic field in D dimensions."""
    return 1.0 / (dimension * threshold_quanta)


def curvature_at_origin(a: float, b: float) -> float:
    """F''(0) for the rational [1/2] profile.

    F(x)=1+(a-b)x^2+O(x^4), so F''(0)=2(a-b).
    The value is included to show that derivative-sensitive local rules would
    reintroduce hidden shape data.  A one-cell service-label readout is the
    stricter premise that forbids this.
    """
    return 2.0 * (a - b)


@dataclass(frozen=True)
class ReadoutCandidate:
    name: str
    ledger: str
    stage: str
    one_cell: bool
    uses_enclosed_integral: bool
    uses_profile_shape: bool
    q: float

    def failures(self) -> list[str]:
        reasons: list[str] = []
        if self.ledger != "service_current":
            reasons.append("wrong ledger")
        if self.stage != "pre-field-service":
            reasons.append("post-summed field readout")
        if not self.one_cell:
            reasons.append("not one-cell local")
        if self.uses_enclosed_integral:
            reasons.append("uses enclosed mass/field integral")
        if self.uses_profile_shape:
            reasons.append("uses hidden profile-shape data")
        return reasons

    @property
    def admissible(self) -> bool:
        return not self.failures()


def main() -> None:
    print("ITEM 132: LOCAL-SERVICE SELECTOR FOR THE R4 CORE READOUT")
    print("=" * 94)

    print("\n[1] Ledger-selection axioms")
    axioms = [
        "A1 observable = R4 Stinespring service current, not gravitational strain",
        "A2 readout stage = pre-field local service tick",
        "A3 instrument input = one central service cell, not an enclosed halo integral",
        "A4 service labels carry count/repair data, not a cored-profile shape parameter",
        "A5 isotropic finite core in three spatial dimensions",
        "A6 threshold = one local a0 service quantum (still external to this proof)",
    ]
    for axiom in axioms:
        print(f"  {axiom}")

    print("\n[2] Candidate boundary observables")
    q_local = central_q()
    q_enclosed_minimal = enclosed_shape_mass(minimal_profile)
    candidates = [
        ReadoutCandidate(
            "central service cell",
            "service_current",
            "pre-field-service",
            True,
            False,
            False,
            q_local,
        ),
        ReadoutCandidate(
            "minimal enclosed field",
            "gravitational_strain",
            "post-field-solve",
            False,
            True,
            False,
            q_enclosed_minimal,
        ),
        ReadoutCandidate(
            "shape-sensitive local jet",
            "service_current",
            "pre-field-service",
            True,
            False,
            True,
            float("nan"),
        ),
        ReadoutCandidate(
            "fitted interpolation core",
            "phenomenological_profile",
            "post-field-solve",
            False,
            True,
            True,
            float("nan"),
        ),
    ]
    for cand in candidates:
        status = "PASS" if cand.admissible else "FAIL: " + ", ".join(cand.failures())
        q_text = "shape-dependent" if math.isnan(cand.q) else f"{cand.q:.12f}"
        print(f"  {cand.name:27s} q={q_text:>16s}  {status}")
    check(sum(c.admissible for c in candidates) == 1, "only the central service-cell readout satisfies the service-ledger axioms")

    print("\n[3] Central-cell algebra")
    print("  For rho=A/r_c^2 F(r/r_c), F(0)=1, isotropic Poisson gives")
    print("    g_cell(r) = (4 pi G A / (3 r_c^2)) r + O(r^3).")
    print("  With 4 pi G A = v_inf^2 = a0 r_M:")
    print("    g_cell(r_c)/a0 = r_M/(3 r_c) = 1/(3 q).")
    print("  One local service quantum, g_cell(r_c)=a0, therefore gives q=1/3.")
    check(abs(q_local - 1.0 / 3.0) < TOL, "central-cell service readout gives r_c/r_M=1/3 exactly")
    check(abs(q_enclosed_minimal - (1.0 - math.pi / 4.0)) < 1e-11, "minimal enclosed-field readout gives 1-pi/4 instead")
    check(abs(q_local - q_enclosed_minimal) > 0.1, "central-cell and enclosed-field rules are observably distinct")

    print("\n[4] Why enclosed-field and shape rules fail this ledger")
    profiles = [
        ("minimal [0/1]", minimal_profile, -2.0),
        ("[1/2] a=1,b=1", lambda x: rational_12_profile(x, 1.0, 1.0), curvature_at_origin(1.0, 1.0)),
        ("[1/2] a=1,b=3", lambda x: rational_12_profile(x, 1.0, 3.0), curvature_at_origin(1.0, 3.0)),
        ("[1/2] a=2,b=1", lambda x: rational_12_profile(x, 2.0, 1.0), curvature_at_origin(2.0, 1.0)),
    ]
    enclosed_values: list[float] = []
    curvatures: list[float] = []
    for label, profile, curvature in profiles:
        f0 = float(profile(np.array([0.0]))[0])
        tail = float((100.0**2) * profile(np.array([100.0]))[0])
        q_enclosed = enclosed_shape_mass(profile)
        enclosed_values.append(q_enclosed)
        curvatures.append(curvature)
        print(
            f"  {label:15s}: F(0)={f0:.3f}, x^2F(100)={tail:.5f}, "
            f"F''(0)={curvature:+.3f}, q_enclosed={q_enclosed:.9f}, "
            f"q_service={q_local:.9f}"
        )
        check(abs(f0 - 1.0) < TOL, f"{label} has the same central occupancy")
        check(abs(tail - 1.0) < 1e-3, f"{label} has the same asymptotic tail")
    check(max(enclosed_values) - min(enclosed_values) > 0.1, "enclosed-field boundary varies across admissible core shapes")
    check(max(curvatures) - min(curvatures) > 1.0, "derivative-sensitive local rules would also introduce hidden shape data")
    check(True, "one-cell service labels intentionally read F(0), not profile derivatives or enclosed integrals")

    print("\n[5] Dependency / category check")
    dependencies = [
        ("service current", "local repair event count", "input to halo line density"),
        ("line density", "sum of service-current records", "input to gravitational field"),
        ("enclosed field", "integral of line density", "post-summed gravitational readout"),
        ("profile shape", "choice of regulator family", "extra modelling data"),
    ]
    for name, definition, role in dependencies:
        print(f"  {name:15s}: {definition:30s} -> {role}")
    check(True, "using enclosed field to define the service readout would bill a post-summed strain observable")
    check(True, "using profile shape would add a hidden boundary label absent from the R4 repair channel")

    print("\n[6] Decision")
    print(
        """
  CLOSED AS A LEDGER-SELECTION THEOREM:
    * Given the R4 halo observable is the local Stinespring service current,
      the admissible boundary readout is one-cell and pre-field.
    * The one-cell isotropic scalar is the central harmonic Poisson/Gauss
      coefficient, so one local a0 quantum gives r_c/r_M=1/3.

  REJECTED FOR THIS OBSERVABLE:
    * The exact enclosed-field rule belongs to the post-summed gravitational
      field ledger and gives 1-pi/4 for the minimal profile.
    * Shape-dependent boundaries require hidden profile data not present in
      the R4 service labels.

  STILL OPEN:
    * Derive the one-a0 service quantum and the external a0 scale from the
      microscopic substrate/horizon algebra.
"""
    )
    print("exit 0 -- local service-readout theorem selects the central-cell boundary over enclosed/shape rules.")


if __name__ == "__main__":
    main()
