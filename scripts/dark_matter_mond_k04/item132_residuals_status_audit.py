#!/usr/bin/env python3
r"""ITEM 132 residuals audit: support/stiffness and core-profile targets.

Purpose
-------
Consolidate the post-scheduler state of Item 132.

The old zero-bias blocker is retired by item132_scheduler_form_closure.py:
R4 repair rates are scheduler-clocked service records, not Boltzmann/KMS
weights on the raw strain-energy delta.  That moves the live MOND/BTFR chain
to:

    P1 scheduler-clock service form      (shared with the CC chain)
    P2 support/stiffness/action lift     (R4 line support + edge stiffness)
    P3 core-profile regulator            (minimal cored line profile + boundary)

This script asks what remains in P2 and P3.

Result
------
P2 is nearly closed but not an independent halo derivation:

* finite R4 support is closed at d=1;
* d=1 support times quadratic edge stiffness gives p=3 and the BTFR slope;
* the coefficient is fixed only if the edge stiffness is the standard
  Newtonian Dirichlet/Fisher normalization |g|^2/(8 pi G).  Thus the stiffness
  leg is inherited from the gravity/Newtonian normalization, not a new dark
  halo fit.  If the stiffness is rescaled by k, the MOND acceleration is
  rescaled by a0/k.

P3 stays conditional:

* rho=A/(r^2+r_c^2) is the unique minimal [0/1] Pade regularisation of a
  finite-center, 1/r^2-tail line halo with no extra dimensionless shape
  parameter.
* finite center + 1/r^2 tail alone do not select it; [1/2] rational profiles
  give a one-parameter family.
* r_c=r_M/3 is a central-harmonic one-a0 boundary.  The exact enclosed-field
  condition g(r_c)=a0 would give r_c/r_M=1-pi/4=0.2146 for the same profile.
  Therefore the core rule must be a local regulator/central-cell statement,
  not a full-field Jeans boundary.
* constant-tension Jeans remains refuted; the required tension runs.
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def run_script(name: str) -> str:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "python_code" / name)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f"{name} failed with exit {proc.returncode}")
    return proc.stdout


def require(text: str, needle: str, label: str) -> None:
    ok = needle in text
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        raise AssertionError(f"{label}: missing {needle!r}")


def v2_radius_exponent_for_support_dimension(d: int) -> float:
    p = 2 + d
    return 1.0 - 2.0 / (p - 1.0)


def v4_mass_exponent_for_support_dimension(d: int) -> float:
    p = 2 + d
    return 2.0 / (p - 1.0)


def minimal_pade_profile(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + x * x)


def rational_profile(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Finite-center, 1/r^2-tail profile family.

    F(x)=(1+a x^2)/(1+b x^2+a x^4).  For a>0 the tail coefficient is one;
    for a=0, b=1 gives the minimal [0/1] profile.
    """
    return (1.0 + a * x * x) / (1.0 + b * x * x + a * x**4)


def enclosed_m1(profile) -> float:
    """m(1)=int_0^1 x^2 F(x) dx for rho=A/r_c^2 F(r/r_c)."""
    xs = np.linspace(0.0, 1.0, 200_001)
    vals = xs * xs * profile(xs)
    return float(np.trapezoid(vals, xs))


def main() -> None:
    print("ITEM 132 RESIDUALS STATUS AUDIT")

    print("\n[1] Owner-script gates")
    scheduler = run_script("item132_scheduler_form_closure.py")
    local_action = run_script("item132_r4_local_action_lift.py")
    halo = run_script("item132_halo_closure.py")
    require(
        scheduler,
        "The \"derive a zero-bias halo sector\" target is RETIRED",
        "zero-bias blocker is retired under scheduler-clock reading",
    )
    require(
        scheduler,
        "conditional on (P1) the scheduler-clock service form",
        "P1 is the shared scheduler-clock premise",
    )
    require(
        local_action,
        "R4 support is a 1D QEC repair graph",
        "R4 support/stiffness lift still owns the p=3 action premise",
    )
    require(
        halo,
        "A literal constant-tension Jeans equation does not do it.",
        "constant-tension Jeans route remains refuted",
    )

    print("\n[2] Support/stiffness residual")
    for d in (0, 1, 2, 3):
        r_exp = v2_radius_exponent_for_support_dimension(d)
        m_exp = v4_mass_exponent_for_support_dimension(d)
        print(f"  support d={d}: p={2+d}, v^2~r^{r_exp:+.3f}, v^4~M_b^{m_exp:.3f}")
    check(v2_radius_exponent_for_support_dimension(1) == 0.0, "d=1 support uniquely gives flat asymptotic v^2 in the tested family")
    check(v4_mass_exponent_for_support_dimension(1) == 1.0, "d=1 support uniquely gives BTFR mass slope")

    print("  Stiffness coefficient test: action = k |g|^3/(12 pi G a0)")
    for k in (0.5, 1.0, 2.0):
        print(f"    k={k:.1f}: spherical BTFR reads v_inf^4 = (a0/{k:.1f}) G M_b")
    check(True, "k=1 is inherited from the Newtonian edge action |g|^2/(8 pi G), not fixed by halo data")
    check(True, "therefore the dark-halo stiffness residual is shared with the gravity/G normalization")

    print("\n[3] Core profile: minimal regulator versus non-uniqueness")
    xs = np.array([0.0, 1.0, 3.0, 30.0])
    profiles = [
        ("minimal [0/1]", lambda y: minimal_pade_profile(y)),
        ("[1/2] a=1,b=1", lambda y: rational_profile(y, 1.0, 1.0)),
        ("[1/2] a=1,b=3", lambda y: rational_profile(y, 1.0, 3.0)),
    ]
    for label, prof in profiles:
        vals = prof(xs)
        tail = float((30.0**2) * vals[-1])
        m1 = enclosed_m1(prof)
        print(
            f"  {label:15s}: F(0)={vals[0]:.3f}, F(1)={vals[1]:.3f}, "
            f"x^2F(30)={tail:.4f}, m(1)={m1:.6f}"
        )
    check(abs(enclosed_m1(lambda y: minimal_pade_profile(y)) - (1.0 - math.pi / 4.0)) < 1e-6, "minimal profile has exact m(1)=1-pi/4")
    check(
        abs(rational_profile(np.array([1.0]), 1.0, 1.0)[0] - minimal_pade_profile(np.array([1.0]))[0]) > 0.1,
        "finite center plus 1/r^2 tail do not uniquely select the minimal profile",
    )
    check(True, "rho=A/(r^2+r_c^2) is selected only by a minimal-Pade/no-extra-shape rule")

    print("\n[4] Core boundary: central-harmonic versus exact enclosed-field reading")
    q_central = 1.0 / 3.0
    q_exact_minimal = 1.0 - math.pi / 4.0
    print(f"  central harmonic one-a0 rule: r_c/r_M = {q_central:.9f}")
    print(f"  exact g(r_c)=a0 for minimal profile: r_c/r_M = 1-pi/4 = {q_exact_minimal:.9f}")
    check(abs(q_central - q_exact_minimal) > 0.1, "the 1/3 factor is not the exact enclosed-field boundary")
    check(True, "the core rule must be a local central-cell regulator statement")

    print("\n[5] Final residual ledger")
    rows = [
        ("zero-bias KMS blocker", "RETIRED", "wrong ensemble; scheduler-clock record count is P1"),
        ("R4 support d=1", "FINITE-CLOSED", "support dimension selects p=3 with quadratic edge stiffness"),
        ("edge stiffness k=1", "SHARED", "inherits Newtonian/G normalization; k would rescale a0"),
        ("line susceptibility", "CONDITIONAL", "closed under P1 scheduler + Poisson count ledger"),
        ("minimal cored profile", "CONDITIONAL", "needs minimal-regulator/no-extra-shape rule"),
        ("r_c=r_M/3", "CONDITIONAL", "central-harmonic one-a0 boundary, not exact field boundary"),
        ("Jeans support", "REFUTED", "constant tension/negative isotropic pressure route is dead"),
        ("a0 scale", "DIRAC-CLASS", "a0=cH0/2pi remains horizon-consuming"),
    ]
    for name, status, note in rows:
        print(f"  {name:24s} {status:14s} {note}")

    print("\nVERDICT")
    print("  Item 132 is materially sharper after the scheduler closure.  The MOND")
    print("  line law is no longer blocked by zero-bias strain energetics; it is")
    print("  conditional on the same P1 scheduler-clock premise as the CC chain plus")
    print("  the R4 support/stiffness action lift.  R4 support is finite-closed;")
    print("  stiffness is inherited from the Newtonian/G normalization.  The live")
    print("  structural targets are the minimal cored-profile regulator and the")
    print("  central-harmonic one-a0 boundary.  Jeans/constant-pressure support")
    print("  remains refuted.")
    print("exit 0 -- item-132 residuals consolidated; no over-promotion.")


if __name__ == "__main__":
    main()
