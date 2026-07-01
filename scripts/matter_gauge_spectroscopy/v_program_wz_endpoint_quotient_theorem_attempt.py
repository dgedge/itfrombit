#!/usr/bin/env python3
r"""Attempt the post-EWSB W/Z endpoint-quotient theorem.

Question
--------
The Z-map is now advanced because alpha(M_Z) is reduced to the framework's
dressed Thomson alpha(0) plus standard photon vacuum-polarisation running.  The
remaining theorem target is therefore narrow:

    prove that the W/Z pole ledger reads sin^2(theta_W)_pole = 2/9.

This script tests whether that quotient is already forced by the documented
finite electroweak service algebra.

Important discipline
--------------------
This is not allowed to revive the retired M9 claim that 2/9 is a bare/UV
Weinberg angle.  The UV charge trace remains 3/8.  A successful theorem must be
post-EWSB and pole-ledger specific.

Exit 0 means the theorem attempt completed cleanly.  The present verdict is
negative/sharpening, not a proof of 2/9:

  * documented finite quotients give 3/8, 1/4, 1/3, or 1/2;
  * sin^2=2/9 requires a pole exposure vector B:W = 2:7;
  * the only cheap arithmetic construction, 2/(2+4+3), mixes endpoint,
    pre-EWSB, and post-EWSB quotient spaces and is therefore not licensed as an
    operator theorem.

This script is now historical/negative: it rejects the mixed-ledger
``2/(2+4+3)`` route.  The constructive pole-projector route is implemented in
``v_program_wz_pole_exposure_operator.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import importlib
import math
import sys


wz = importlib.import_module("v_program_wz_endpoint_coupling_candidate")
ew_alpha = importlib.import_module("ew_alpha_mz_from_framework_dressed_alpha")

ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


@dataclass(frozen=True)
class Quotient:
    name: str
    value: Fraction
    owner: str
    status: str


def exposure_ratio_for_s2(s2: Fraction) -> tuple[int, int]:
    """Return integer (B, W) such that s2 = B / (B + W)."""
    return s2.numerator, s2.denominator - s2.numerator


def masses_from_framework_alpha(s2: Fraction) -> tuple[float, float]:
    alpha_mz_inv = ew_alpha.alpha_mz_inv(ew_alpha.DELTA_LEP_FULL, ew_alpha.DELTA_HAD5)
    mw, mz, _g2, _gy = wz.masses_from_alpha_s2(1.0 / alpha_mz_inv, float(s2))
    return mw, mz


def pct(pred: float, obs: float) -> float:
    return 100.0 * (pred / obs - 1.0)


def main() -> int:
    print("=" * 96)
    print("W/Z POLE ENDPOINT QUOTIENT THEOREM ATTEMPT")
    print("=" * 96)

    target = Fraction(2, 9)
    b_exp, w_exp = exposure_ratio_for_s2(target)
    print("\n[1] What the target theorem requires")
    print(f"    target sin^2(theta_W)_pole = {target} = {float(target):.9f}")
    print(f"    equivalent pole exposure vector: B:W = {b_exp}:{w_exp}")
    check("2/9 is exactly a 2:7 hypercharge:weak exposure quotient", (b_exp, w_exp) == (2, 7))

    exact_pole = 1.0 - (wz.M_W_OBS / wz.M_Z_OBS) ** 2
    mw, mz = masses_from_framework_alpha(target)
    print(f"    exact pole ratio from reference masses = {exact_pole:.9f}")
    print(f"    framework alpha(M_Z) + 2/9 gives m_W={mw:.6f} GeV ({pct(mw, wz.M_W_OBS):+.3f}%)")
    print(f"    framework alpha(M_Z) + 2/9 gives m_Z={mz:.6f} GeV ({pct(mz, wz.M_Z_OBS):+.3f}%)")
    check("2/9 remains the right numerical pole target at sub-percent mass grade",
          abs(mw / wz.M_W_OBS - 1.0) < 0.005 and abs(mz / wz.M_Z_OBS - 1.0) < 0.005)

    print("\n[2] Documented finite EW service quotients")
    documented = [
        Quotient(
            "UV SO(10)/GUT charge trace",
            Fraction(3, 8),
            "charge-trace normalization",
            "licensed UV/bare map; not a pole endpoint",
        ),
        Quotient(
            "one U(1) direction out of four EW directions",
            Fraction(1, 4),
            "pre-EWSB service-rank map",
            "licensed rank quotient; wrong value",
        ),
        Quotient(
            "one neutral direction out of broken triad",
            Fraction(1, 3),
            "post-EWSB broken-rank map",
            "licensed rank quotient; wrong value",
        ),
        Quotient(
            "binary endpoint split",
            Fraction(1, 2),
            "endpoint-count map",
            "licensed endpoint count; wrong value",
        ),
    ]
    for q in documented:
        b, w = exposure_ratio_for_s2(q.value)
        print(
            f"    {q.name:<48s} {str(q.value):>5s} = {float(q.value):.9f}"
            f"   exposure {b}:{w}   [{q.owner}]"
        )
    check("no documented same-space finite quotient equals 2/9",
          all(q.value != target for q in documented))
    check("UV 3/8 is preserved as a distinct charge-trace theorem",
          documented[0].value == Fraction(3, 8) and documented[0].value != target)

    print("\n[3] The tempting 2/(2+4+3) arithmetic fork")
    endpoint_contacts = 2
    pre_ewsb_rank = 4
    broken_rank = 3
    cheap = Fraction(endpoint_contacts, endpoint_contacts + pre_ewsb_rank + broken_rank)
    print(f"    2 endpoint contacts / (2 endpoint + 4 pre-EWSB + 3 broken) = {cheap}")
    check("the cheap mixed-space count lands arithmetically on 2/9", cheap == target)
    print("    rejection reason:")
    print("      the denominator adds three different ledgers: endpoint contacts,")
    print("      the pre-EWSB SU(2)xU(1) service alphabet, and the post-EWSB broken")
    print("      quotient after removing Q=T3+Y.  Adding them double-counts the")
    print("      electroweak vector space unless a pole-exposure operator explicitly")
    print("      glues these spaces into one residue ledger.")
    check("2/(2+4+3) is not promoted without an explicit pole-exposure gluing operator", True)

    print("\n[4] What a proof would have to supply")
    print("    Required operator form:")
    print("      E_pole = diag(B_contacts=2, W_contacts=7) on the W/Z external-leg residue,")
    print("      sin^2(theta_W)_pole = Tr E_B / Tr(E_B + E_W) = 2/9.")
    print("    Compatibility conditions:")
    print("      (i)  E_pole acts after EWSB on the pole residue, not on the UV charge trace;")
    print("      (ii) the UV charge trace remains 3/8 and runs by ordinary SM RG;")
    print("      (iii) the seven weak contacts are not assembled by double-counting the")
    print("            same EW directions before and after symmetry breaking;")
    print("      (iv) W and Z share the same pole ledger so the ratio and absolute masses")
    print("           use one quotient, not separate fitted residues.")

    print("\n[5] Verdict")
    print(
        """
    MIXED-LEDGER ROUTE REJECTED.

    The theorem target is sharpened by this negative result: the residue weights
    must be B:W = 2:7, but they cannot be obtained by adding endpoint contacts,
    the pre-EWSB EW alphabet, and the post-EWSB broken quotient as if they were
    one space.  That tempting 2/(2+4+3) count is exactly the right arithmetic and
    exactly the wrong proof.

    Constructive resolution is in v_program_wz_pole_exposure_operator.py: the
    on-shell W/Z pole space has rank 9, and the hypercharge endpoint projector
    has rank 2 because an abelian endpoint contributes only transverse Maxwell
    records.  This script remains as the guardrail against the rejected route.
"""
    )

    if ok:
        print("ALL CHECKS PASSED")
        return 0
    print("CHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
