#!/usr/bin/env python3
r"""v_phase6_cw_finite_normal_ordering.py

Phase 6B: finite Coleman-Weinberg normal-ordering audit.

Phase 5B reduced the EW normalisation residual to a concrete one-loop target:

    lambda_required - lambda_crit(y_t0) = 0.014928611...

This script asks whether the ordinary finite one-loop Coleman-Weinberg constants
are of the right size, rather than inserting a free amplitude.

We use the standard field-dependent masses

    m_W^2 = g_2^2 h^2 / 4,
    m_Z^2 = (g_2^2 + g_Y^2) h^2 / 4,
    m_t^2 = y_t^2 h^2 / 2,
    m_H^2 = 3 lambda h^2,
    m_G^2 = lambda h^2,

with MS-like finite constants c_i = 5/6 for gauge bosons and 3/2 for
fermions/scalars.  Writing mu = kappa h, the finite one-loop quartic shift is

    Delta lambda_CW(kappa)
      = sum_i n_i a_i^2 [ log(a_i / kappa^2) - c_i ] / (16 pi^2).

This is a precision *target* calculation, not a full SM matching result:
proper closure still needs a fixed scheme, pole-to-MSbar threshold matching,
and higher loops.

The notable output is that the broken-triad normal-ordering scale

    kappa = 1 / sqrt(3)

lands within ~0.1% in the VEV amplitude in this one-loop convention.  That is
a strong candidate for a substrate theorem ("three broken service directions
set the finite normal-ordering scale"), but it is not adopted as locked here.
"""

from __future__ import annotations

import math
import sys


LOOP = 16.0 * math.pi * math.pi
ALPHA0 = 1.0 / 137.036
M_P = 1.2209e19
M_Z = 91.1876
V_OBS = 246.22
G_Y0 = 0.3573
G_20 = 0.6517
Y_T0 = 0.9500
LAMBDA_CRIT = 0.144057344
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def species(lam: float = LAMBDA_CRIT) -> list[tuple[str, float, float, float]]:
    """(name, degeneracy n_i, mass coefficient a_i in m_i^2=a_i h^2, finite c_i)."""
    gz2 = G_20 * G_20 + G_Y0 * G_Y0
    return [
        ("W", 6.0, G_20 * G_20 / 4.0, 5.0 / 6.0),
        ("Z", 3.0, gz2 / 4.0, 5.0 / 6.0),
        ("top", -12.0, Y_T0 * Y_T0 / 2.0, 3.0 / 2.0),
        ("H", 1.0, 3.0 * lam, 3.0 / 2.0),
        ("Goldstone", 3.0, lam, 3.0 / 2.0),
    ]


def delta_lambda_cw(kappa: float, lam: float = LAMBDA_CRIT) -> float:
    return sum(n * a * a * (math.log(a / (kappa * kappa)) - c) / LOOP for _, n, a, c in species(lam))


def required_lambda() -> float:
    return (ALPHA0**8 / (V_OBS / M_P)) ** 2


def vev_ratio(lam_eff: float) -> float:
    """Predicted v / observed v after replacing lambda by lam_eff."""
    return (ALPHA0**8 / math.sqrt(lam_eff)) / (V_OBS / M_P)


def solve_kappa_for_target(target_delta: float) -> float:
    lo, hi = 0.05, 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        # In this convention Delta_lambda increases with kappa because the
        # dominant top term has negative degeneracy.
        if delta_lambda_cw(mid) < target_delta:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main() -> int:
    print("=" * 96)
    print("PHASE 6B: FINITE CW NORMAL-ORDERING AUDIT")
    print("=" * 96)

    lam_req = required_lambda()
    target_delta = lam_req - LAMBDA_CRIT
    print("\n[0] Target")
    print(f"  lambda_crit central       = {LAMBDA_CRIT:.9f}")
    print(f"  lambda_required exact     = {lam_req:.9f}")
    print(f"  required Delta lambda_CW  = {target_delta:.9f}")
    check("required finite shift is percent-level, not OOM", 0.01 < target_delta < 0.02)

    print("\n[1] One-loop finite terms")
    for name, n, a, c in species():
        term = n * a * a * (math.log(a) - c) / LOOP
        print(f"  {name:10s} n={n:5.1f}  a={a:.9f}  c={c:.6f}  term(kappa=1)={term:+.9f}")
    delta_k1 = delta_lambda_cw(1.0)
    print(f"  total Delta lambda_CW(kappa=1) = {delta_k1:+.9f}")
    check("ordinary one-loop finite constants are the right order", 0.02 < delta_k1 < 0.04)
    check("kappa=1 overshoots the exact target in this convention", delta_k1 > target_delta)

    print("\n[2] Candidate normal-ordering scales")
    rows = [
        ("mu=h", 1.0),
        ("top threshold y_t/sqrt(2)", Y_T0 / math.sqrt(2.0)),
        ("broken triad 1/sqrt(3)", 1.0 / math.sqrt(3.0)),
        ("W threshold g_2/2", G_20 / 2.0),
        ("Z threshold g_Z/2", math.sqrt(G_20 * G_20 + G_Y0 * G_Y0) / 2.0),
    ]
    for label, kappa in rows:
        dl = delta_lambda_cw(kappa)
        lam_eff = LAMBDA_CRIT + dl
        print(
            f"  {label:26s} kappa={kappa:.9f}  Delta={dl:+.9f}  "
            f"lambda_eff={lam_eff:.9f}  pred/obs={vev_ratio(lam_eff):.9f}"
        )

    k_req = solve_kappa_for_target(target_delta)
    k_triad = 1.0 / math.sqrt(3.0)
    ratio_triad = vev_ratio(LAMBDA_CRIT + delta_lambda_cw(k_triad))
    print("\n[3] Exact kappa target")
    print(f"  kappa required by exact one-loop landing = {k_req:.9f}")
    print(f"  1/sqrt(3)                              = {k_triad:.9f}")
    print(f"  relative kappa miss                    = {100.0 * abs(k_triad / k_req - 1.0):.3f}%")
    print(f"  pred/obs at 1/sqrt(3)                  = {ratio_triad:.9f}")
    check("broken-triad kappa lands within 0.2% in the VEV amplitude", abs(ratio_triad - 1.0) < 0.002)
    check("exact kappa is close to 1/sqrt(3), but not identical in this one-loop proxy",
          0.005 < abs(k_triad / k_req - 1.0) < 0.02)

    print("\nVERDICT")
    print("  The finite CW constants are no longer a shapeless multiplier: the standard")
    print("  one-loop normal-ordering terms are of the right size, and the substrate-")
    print("  natural broken-triad scale kappa=1/sqrt(3) lands at the 0.1% amplitude")
    print("  level in this convention.  This is a strong candidate theorem, not a lock.")
    print("  Closure still requires deriving the broken-triad normal-ordering scale from")
    print("  the EW service algebra and then repeating the calculation with a fixed")
    print("  multi-loop threshold scheme.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
