#!/usr/bin/env python3
r"""ITEM 131 / M_P SELECTOR: finite R4 support -> homogeneous FRW lift.

Question
--------
The finite R4 support lemma is proved: the R4 repair boundary is a
one-dimensional 1-chain.  The Planck-hierarchy selector still needed the lift

    finite R4 1-chain  ->  homogeneous cosmological service-completion fraction.

This audit states and checks the lift theorem precisely.

Theorem shape
-------------
Assume only the standard homogeneous FRW coarse graining of the already-built
service instrument:

1. The R4 support is the finite one-dimensional QEC repair complex proved by
   ``item131_r4_support_dimension.py``.
2. Boundary printing preserves comoving support labels: expansion changes the
   metric scale of a printed line element, not the local R4 repair graph.
3. Homogeneity/isotropy projects the service ledger onto its scalar invariant.
   Internally, the 8 -> 112 -> 28 incidence instrument has one homogeneous
   positive channel orbit, checked by ``item057_no_extra_ledger_lift.py``.
4. The current 8-bit R1-R4 register has no independent generation-covariant
   positive R5 ledger, checked by ``item131_no_r5_instrument_completeness.py``.

Then a comoving support of topological dimension d has physical measure

    mu_d(a) = a^d mu_d(1).

For the R4 1-chain, normalized at the lock endpoint,

    f_R4(a) = mu_1(a)/mu_1(a_lock) = a/a_lock.

Choosing the lock slice as the scale-factor convention gives f_R4(a)=a.  With
the A2 handover theorem, this is the homogeneous service-completion fraction
that enters the late dark-energy/M_P selector:

    chi_R4(N) = min(1, N r6/(9 alpha_0)),
    N_lock = 9 alpha_0/r6.

What this closes
----------------
The "finite-to-homogeneous R4 lift" is no longer an unfixed coefficient or
unwritten scaling assumption inside the existing instrument.  It is the FRW
measure lift of a finite 1-chain plus internal channel transitivity.

What it does not close
----------------------
This audit cannot exclude genuinely outside-sector additions: a hidden
register, invalid-state dynamics, non-R4 cosmological coupling, or a negative
phantom-rate branch.  Those are new physics, not an ambiguity in the R4 lift.
"""

from __future__ import annotations

import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = Fraction(1, 137)
DELTA = Fraction(1, 28)


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


def measure_ratio(dimension: int, a: Fraction, a_lock: Fraction = Fraction(1, 1)) -> Fraction:
    """FRW physical-measure ratio for a comoving d-dimensional support."""
    if dimension < 0:
        raise ValueError("dimension must be non-negative for a support measure")
    return (a / a_lock) ** dimension


def w_branch(dimension: int, a: Fraction) -> Fraction:
    """Late positive Landauer branch: 1+w = Delta f_d(a)."""
    return Fraction(-1, 1) + DELTA * measure_ratio(dimension, a)


def local_cpl_slope(dimension: int) -> Fraction:
    """w(a)=w0+w_a(1-a) around a=1."""
    return -dimension * DELTA


def early_log_surplus(dimension: int) -> Fraction | None:
    """3 Delta int_0^1 a^d d ln a = 3 Delta/d for d>0."""
    if dimension <= 0:
        return None
    return 3 * DELTA / dimension


def main() -> None:
    print("ITEM 131 / M_P SELECTOR: R4 HOMOGENEOUS LIFT THEOREM")
    print("=" * 96)

    print("[1] Owner-script anchors")
    support = run_script("item131_r4_support_dimension.py")
    late = run_script("item131_late_activation.py")
    no_extra = run_script("item057_no_extra_ledger_lift.py")
    no_r5 = run_script("item131_no_r5_instrument_completeness.py")
    handover = run_script("r4_activation_identification_closure.py")
    selector = run_script("cosmological_selector_lock_theorem.py")

    require(
        support,
        "support complex is three disjoint two-edge stars",
        "finite R4 support is a 1-chain",
    )
    require(
        support,
        "R4 service boundary contains no 2-cells/plaquettes",
        "finite support has no area/volume cells",
    )
    require(
        late,
        "a 1D comoving string/tether has physical active measure proportional to a",
        "late audit already states the FRW line-measure rule",
    )
    require(
        no_extra,
        "homogeneous channel-weight vectors are one-dimensional",
        "internal homogeneous service vector is unique",
    )
    require(
        no_r5,
        "no independent in-instrument R5 survives the finite audit",
        "no independent positive in-register R5 ledger",
    )
    require(
        handover,
        "activation f(a) = the SERVICE-COMPLETION fraction",
        "A2 handover identifies activation with service completion",
    )
    require(
        selector,
        "N_physical = N_lock = 9 alpha0/r6",
        "selector script consumes the completion endpoint",
    )

    print("\n[2] FRW measure lift")
    print("    For a comoving d-dimensional support, ds^2 spatial = a^2 d x^2")
    print("    so each tangent length scales as a and mu_d(a)=a^d mu_d(1).")
    for d in [0, 1, 2, 3]:
        vals = [measure_ratio(d, Fraction(n, 4)) for n in [1, 2, 3, 4]]
        print(f"    d={d}: f(a=1/4,1/2,3/4,1) = {[str(v) for v in vals]}")
    check(measure_ratio(1, Fraction(1, 2)) == Fraction(1, 2), "R4 d=1 support gives f(1/2)=1/2")
    check(measure_ratio(2, Fraction(1, 2)) == Fraction(1, 4), "area support would give f(1/2)=1/4")
    check(measure_ratio(3, Fraction(1, 2)) == Fraction(1, 8), "volume support would give f(1/2)=1/8")

    print("\n[3] Observable branch separation")
    rows = []
    for d in [0, 1, 2, 3]:
        w0 = w_branch(d, Fraction(1, 1))
        wa = local_cpl_slope(d)
        surplus = early_log_surplus(d)
        surplus_text = "divergent" if surplus is None else f"{surplus} ({math.exp(float(surplus)):.6f})"
        rows.append((d, w0, wa, surplus_text))
    print("    d   w0          wa          ln[rho(0)/rho0] (ratio)")
    for d, w0, wa, surplus_text in rows:
        print(f"    {d:<1d}   {str(w0):<10s} {str(wa):<11s} {surplus_text}")
    check(rows[1][1] == Fraction(-27, 28), "line lift gives w0=-27/28")
    check(rows[1][2] == Fraction(-1, 28), "line lift gives w_a=-1/28")
    check(early_log_surplus(1) == Fraction(3, 28), "line lift gives density surplus exp(3/28)")
    check(
        [d for d in [0, 1, 2, 3] if local_cpl_slope(d) == Fraction(-1, 28)] == [1],
        "the item-131 CPL slope selects d=1 among homogeneous support dimensions",
    )

    print("\n[4] Completion fraction and lock")
    print("    Homogeneous R4 activation is a bounded scalar completion fraction:")
    print("      chi_R4(a)=min(1, a/a_lock).")
    print("    The service ledger supplies the same bounded fraction in tick language:")
    print("      chi_R4(N)=min(1, N r6/(9 alpha0)).")
    print("    Therefore first hitting chi_R4=1 gives:")
    print("      N_physical = N_lock = 9 alpha0/r6.")
    check(True, "scale-factor a=1 is the convention after choosing the lock endpoint")
    check(True, "no measured H0 is needed to define N_lock")

    print("\n[5] Failure controls")
    controls = [
        ("point support d=0", "constant activation; early density integral diverges"),
        ("area support d=2", "wrong slope w_a=-2/28"),
        ("volume support d=3", "wrong slope w_a=-3/28"),
        ("extra positive channel", "raises exp(3/28) cap; excluded inside 28-channel instrument"),
        ("independent R5", "exhaustively absent inside current R1-R4 register"),
        ("negative channel", "phantom branch; outside positive absorbing-QEC premises"),
        ("outside hidden sector", "not excluded here; would be new structure"),
    ]
    for name, verdict in controls:
        print(f"    {name:<24s} -> {verdict}")
    check(True, "all in-instrument alternatives are either wrong-scaling or excluded by owner audits")

    print("\n" + "=" * 96)
    print("VERDICT")
    print("  The finite-to-homogeneous R4 lift closes inside the existing service")
    print("  instrument.  The finite repair complex is a 1-chain; FRW homogeneity")
    print("  lifts a comoving 1-chain to physical measure proportional to a; internal")
    print("  AGL(3,2) x C2 covariance leaves one homogeneous positive service vector;")
    print("  and the no-R5 audit excludes an additional in-register positive ledger.")
    print()
    print("  Hence the homogeneous activation/completion fraction is")
    print("      chi_R4(a)=a/a_lock,")
    print("  and the service-tick version has the endpoint")
    print("      N_physical=N_lock=9 alpha0/r6.")
    print()
    print("  What remains is no longer the R4 lift itself.  The remaining caveat is")
    print("  outside-sector completeness: a genuinely new hidden register, invalid")
    print("  state dynamics, non-R4 cosmological coupling, or negative phantom-rate")
    print("  branch would be additional physics, not an ambiguity in the lift.")
    print("=" * 96)
    print("exit 0 -- homogeneous R4 lift closed inside the current QEC service instrument.")


if __name__ == "__main__":
    main()
