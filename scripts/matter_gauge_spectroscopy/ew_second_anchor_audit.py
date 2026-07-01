#!/usr/bin/env python3
r"""The EW/nuclear SECOND ABSOLUTE SCALE — a pre-registered second-anchor audit (E6 #2, ANCHOR §16.2).

The framework has ONE saturated absolute anchor: Lambda_QCD (confinement / K04 crystallisation). Canon
HOMES the electroweak scale structurally -- v = 246 GeV is the R4 Feshbach-resonance pole / R4
constraint-violation energy (ANCHOR §5.8, §2.13) -- but ASSIGNS its magnitude by hand (§15 item 53:
"the Higgs VEV is phenomenologically assigned ... without a rigorous substrate-level injection
theorem"). E6 #2 lists this as the one genuinely-open NATIVE coefficient.

Question (the user's): a SECOND saturated anchor, or a mechanism that derives v from {Lambda, M_P, alpha0}.

This audit does three disciplined things instead of number-matching:
  [A] PRECISION ASYMMETRY -- the structural reason v cannot be a clean consequence of Lambda:
      Lambda_QCD is SOFT (scheme/scale dependent, ~few %); v and m_top are SHARP (0.01-0.4%). A sharp
      quantity cannot be DERIVED from a soft one below the soft one's own precision -> any v=f(Lambda)
      is untestable below ~few %, and at ~few % the form space is dense (numerology floor, shown in [B]).
  [B] ACCIDENT GRADE -- enumerate low-complexity native forms v/Lambda = 2^p alpha0^q c and count how
      many land within each tolerance of the measured 741.6. If many land at the Lambda-floor (~few %)
      and none lands at v's own precision (0.01%), then NO number-match is informative.
  [C] THE CANDIDATE SECOND ANCHOR -- top-Yukawa saturation. y_t = sqrt2 m_top/v is the UNIQUE O(1)
      Yukawa (an RG IR quasi-fixed point); y_t ~= 1 => v = sqrt2 m_top. The heaviest fermion's
      saturated Yukawa IS a genuine second SATURATED anchor (parallel to Lambda = QCD saturation) --
      but it anchors v to m_top (one irreducible EW input), it does NOT reduce the EW scale to Lambda.

Verdict: the second saturated anchor is identified (top-Yukawa quasi-fixed-point; m_top == v), and it is
IRREDUCIBLE to the QCD anchor -- the framework genuinely needs TWO mass anchors. Confirms + sharpens E6 #2.
"""
from __future__ import annotations

import itertools
import math

# --- measured inputs (PDG) ---
V_EW = 246.22          # GeV, Higgs VEV (from G_F; ~0.01%)
M_TOP_POLE = 172.76    # GeV pole (~0.4%); MSbar m_t(m_t) ~= 163 -> y_t scheme spread
M_TOP_MSB = 163.0
LAMBDA = 0.332         # GeV, the framework's QCD anchor
LAMBDA_REL_UNC = 0.05  # ~few % : Lambda_QCD is scheme/lattice/scale dependent (the SOFT anchor)
ALPHA0 = 1.0 / 137.0
M_N = 2.0 * math.sqrt(2.0) * LAMBDA   # framework proton anchor M_N = 2 sqrt2 Lambda


def bar(c="="):
    print(c * 100)


def main():
    bar()
    print("EW SECOND ABSOLUTE SCALE — pre-registered second-anchor audit (E6 #2, ANCHOR §16.2)")
    bar()

    # ---------- [A] precision asymmetry: the structural no-derive-from-Lambda argument ----------
    print("[A] PRECISION ASYMMETRY (the structural obstruction)")
    print(f"    Lambda_QCD = {LAMBDA} GeV, SOFT: scheme/scale dependent, rel. uncertainty ~{LAMBDA_REL_UNC:.0%}")
    print(f"    v          = {V_EW} GeV, SHARP: ~0.01% (from G_F)")
    print(f"    m_top      = {M_TOP_POLE} GeV (pole, ~0.4%)")
    r = V_EW / LAMBDA
    print(f"    => v/Lambda = {r:.1f}, but inherits Lambda's softness: {r:.0f} +/- {r*LAMBDA_REL_UNC:.0f} (~{LAMBDA_REL_UNC:.0%})")
    print("    A SHARP quantity (v, 0.01%) cannot be DERIVED from a SOFT one (Lambda, few %) below")
    print("    Lambda's own precision -> no v=f(Lambda) is testable better than ~few %.")
    assert LAMBDA_REL_UNC > 50 * 1e-4, "Lambda softness (~5%) dwarfs v's ~0.01% precision -> v/Lambda is a few-% statement"

    # ---------- [B] accident grade: how dense is the low-complexity form space near 741.6? ----------
    print("\n[B] ACCIDENT GRADE — low-complexity native forms v/Lambda = 2^p * alpha0^q * c")
    target = V_EW / LAMBDA
    consts = {"1": 1.0, "sqrt2": math.sqrt(2), "sqrt3": math.sqrt(3), "2": 2.0, "3": 3.0,
              "pi": math.pi, "phi": (1 + 5 ** 0.5) / 2, "5": 5.0, "7": 7.0}
    powers2 = [p / 2 for p in range(0, 25)]        # 2^0 .. 2^12 in half-steps
    powersA = [-2, -1, 0, 1, 2]                     # alpha0^q
    forms = []
    for p, q, (cname, c) in itertools.product(powers2, powersA, consts.items()):
        val = (2.0 ** p) * (ALPHA0 ** q) * c
        forms.append((abs(val / target - 1.0), f"2^{p:g} a0^{q} {cname}", val))
    forms.sort()
    for tol in (0.05, 0.02, 0.005, 0.001):
        n = sum(1 for d, _, _ in forms if d <= tol)
        print(f"    within {tol*100:>4.1f}% of v/Lambda={target:.1f}:  {n:>4d} of {len(forms)} forms")
    print("    best 4 low-complexity forms (all incompatible, none at v's 0.01% precision):")
    for d, name, val in forms[:4]:
        print(f"      {name:<16s} = {val:7.1f}  ({(val/target-1)*100:+.2f}%)")
    n_floor = sum(1 for d, _, _ in forms if d <= LAMBDA_REL_UNC)
    n_sharp = sum(1 for d, _, _ in forms if d <= 0.001)
    assert n_floor >= 5, "at the Lambda softness floor (~few %) the form space is dense -> numerology"
    assert n_sharp == 0 or forms[0][0] > 0.001, "no low-complexity form reaches v's own 0.01-0.1% precision"
    print(f"    => nothing within ~2.4%; the {n_floor} forms within Lambda's own ~{LAMBDA_REL_UNC:.0%} softness are "
          f"scattered/incompatible; at 0.1%: {n_sharp}. No v=f(Lambda) survives Lambda's softness.")

    # ---------- [C] the candidate SECOND SATURATED anchor: top-Yukawa saturation ----------
    print("\n[C] CANDIDATE SECOND SATURATED ANCHOR — top-Yukawa quasi-fixed-point saturation")
    yt_pole = math.sqrt(2) * M_TOP_POLE / V_EW
    yt_msb = math.sqrt(2) * M_TOP_MSB / V_EW
    print(f"    y_top = sqrt2 m_top / v = {yt_pole:.3f} (pole) .. {yt_msb:.3f} (MSbar)  -- the UNIQUE O(1) Yukawa")
    print("    y_t ~= 1 is an RG IR quasi-fixed point (Pendleton-Ross/Hill): the heaviest fermion's")
    print("    Yukawa SATURATES at the natural O(1) value => v = sqrt2 m_top / y_t ~= sqrt2 m_top.")
    print(f"    So the EW scale is anchored by the TOP MASS via saturation: v = sqrt2 * {M_TOP_POLE} / {yt_pole:.3f} = {math.sqrt(2)*M_TOP_POLE/yt_pole:.1f} GeV.")
    assert 0.90 < yt_pole < 1.05, "top Yukawa is O(1) (the saturation), within scheme spread of 1"
    # this is a genuine SECOND saturated anchor, parallel to Lambda = QCD-confinement saturation,
    # but it anchors v to m_top -- it does NOT reduce the EW scale to Lambda.
    mtop_over_lambda = M_TOP_POLE / LAMBDA
    print(f"    Residual: m_top/Lambda = {mtop_over_lambda:.0f} -- and by [A]/[B] this has NO clean native")
    print(f"    form testable below Lambda's ~{LAMBDA_REL_UNC:.0%} (e.g. 2^9 a0^0 1 = {2**9} is {(2**9/mtop_over_lambda-1)*100:+.1f}%, numerology-grade).")

    print(f"""
{"=" * 100}
VERDICT (pre-registered, exit 0):  a second saturated anchor IS identified -- the top-Yukawa quasi-
fixed point (y_t ~= 1 => v = sqrt2 m_top) -- and it is IRREDUCIBLE to the QCD anchor.

  FOUND (the constructive answer): the natural second SATURATED anchor is top-Yukawa saturation.
  y_t = sqrt2 m_top/v ~= 1 is the unique O(1) Yukawa and an RG IR quasi-fixed point, so the heaviest
  fermion's coupling saturates and v = sqrt2 m_top. This is a genuine saturation anchor, parallel to
  Lambda being the QCD-confinement saturation -- the framework's TWO mass anchors are then {{Lambda (QCD,
  light hadrons), m_top == v (EW, via y_t->1)}}, and all EW masses (W, Z, h, the Yukawa hierarchy ratios
  already probed in tch_yukawa_hierarchy_probe.py) follow from v + standard EW.

  IRREDUCIBLE (the honest limit): the EW scale does NOT reduce to Lambda. [A] Lambda is SOFT (~few %)
  while v/m_top are SHARP (0.01-0.4%): a sharp scale cannot be derived from a soft one below the soft
  one's precision. [B] NO low-complexity native form lands within ~2.4% of v/Lambda, and the handful
  within Lambda's own ~5% softness are scattered and mutually incompatible (untestable below Lambda's
  precision); none reaches v's 0.01%. So no v=f(Lambda) number-match is informative. The framework
  genuinely needs a SECOND independent mass input; it cannot be economised away.

  This CONFIRMS and SHARPENS E6 #2 / §16.2: not "find a number for v from Lambda" (shown impossible/
  numerology), but "accept m_top==v as the irreducible second saturated anchor (top quasi-fixed point)."
  The only deeper reduction that could exist is the canon-named Feshbach INJECTION theorem (§15 item 53):
  derive the energy UNIT that dimensionalises the R4 pole at E=1 -- and that unit, to beat numerology,
  must be pinned to a SHARP scale (M_P, not the soft Lambda), with kill criterion: a unique
  <0.5%-accident-grade form, not a few-% near-hit.
{"=" * 100}""")
    print(f"exit 0 -- EW second anchor = top-Yukawa saturation (y_t={yt_pole:.2f}~=1 => v=sqrt2 m_top); "
          f"irreducible to soft Lambda (numerology floor {n_floor} forms at {LAMBDA_REL_UNC:.0%}, 0 at 0.1%); confirms+sharpens E6 #2.")


if __name__ == "__main__":
    main()
