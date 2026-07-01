#!/usr/bin/env python3
r"""fermion_q4_dynamical_coupling.py -- the dynamical-coupling attempt to PROMOTE the Q4 +1 selector.

Condition (2) (the §16.3 count) could not force the |Q|-neutrality selector: the outcome space is
four-fold (25% baseline) and neutrality vs isospin are degenerate on the derived lepton case. The only
route left is a DYNAMICAL mechanism -- a coupling argument that *derives* which doublet member carries
the +1 frozen defect. This script tries it and follows it honestly.

THE DERIVED LEPTON MECHANISM (canon). The base d=2 channels into the sterile RH reference nu_R are:
  nu_L -> nu_R  via chi-flip  (I3 PRESERVED, Delta Q = 0: both neutral),
  e_R  -> nu_R  via I3-flip   (I3 ACTIVE,    Delta Q = 1: e_R charged, nu_R neutral).
Canon derives the +1 as the CNOT-FROZEN I3 bit: the neutrino's I3 is frozen (its channel preserves
I3), the electron's is active (its channel flips I3). So the +1 attaches to the member whose channel
to the reference PRESERVES I3 = the member with the SAME isospin as the reference. nu_R is up-type, so
the +1 goes to the up-type lepton (nu). This is mechanism A = ISOSPIN-ALIGNMENT (frozen-I3), and it
is the *derived* lepton rule.

EXTENDING A TO QUARKS. The structural analogue of nu_R (the up-type RH singlet) is u_R. The
I3-preserving (chi-flip) channel into u_R comes from u_L (the up-type quark). So mechanism A puts the
quark +1 on the UP quark -> the set {nu, u}.

THE COMPETITOR (mechanism B = NEUTRALITY). The +1 = extra sterility/decoupling, scaling with
neutrality (less gauge charge -> more decoupled -> more defect): +1 on the smaller-|Q| member -> {nu,d}.
B matches the empirical {nu,d} (this is the candidate). A and B AGREE on leptons (nu is both up-type
AND neutral) and SPLIT on quarks.

THE TEST. Predict all four twists delta = (base 2 + the +1)/N_eff under A and B, compare to data.
"""
from fractions import Fraction as F

MEM = {  # (doublet, isospin up=+1/down=-1, |Q|, coloured)
    "nu": ("L", +1, F(0),    False),
    "e":  ("L", -1, F(1),    False),
    "u":  ("Q", +1, F(2, 3), True),
    "d":  ("Q", -1, F(1, 3), True),
}
DBL = {"L": ["nu", "e"], "Q": ["u", "d"]}
N_EFF = {m: (27 if MEM[m][3] else 9) for m in MEM}       # quarks colour-share (N*N_c), leptons N
EMPIRICAL = {"e": F(2, 9), "nu": F(3, 9), "u": F(2, 27), "d": F(1, 9)}
REF_ISOSPIN = +1   # the RH reference is the UP-TYPE singlet (nu_R for leptons, u_R for quarks)


def plus_one_isospin(m):   # mechanism A: +1 on the member whose isospin matches the up-type reference
    return 1 if MEM[m][1] == REF_ISOSPIN else 0


def plus_one_neutrality(m):  # mechanism B: +1 on the smaller-|Q| member of the doublet
    dbl = MEM[m][0]
    return 1 if m == min(DBL[dbl], key=lambda x: abs(MEM[x][2])) else 0


def delta(m, plus_one_fn):
    return F(2 + plus_one_fn(m), N_EFF[m])


def show(label, fn):
    print(f"  {label}")
    print(f"    {'sector':6s} {'+1?':>4s} {'delta_pred':>10s} {'delta_emp':>10s} {'match':>6s}")
    oks = {}
    for m in ("e", "nu", "u", "d"):
        dp = delta(m, fn)
        ok = (dp == EMPIRICAL[m])
        oks[m] = ok
        print(f"    {m:6s} {plus_one(fn, m):>4d} {str(dp):>10s} {str(EMPIRICAL[m]):>10s} {str(ok):>6s}")
    lept = oks["e"] and oks["nu"]
    quark = oks["u"] and oks["d"]
    print(f"    -> leptons {'OK' if lept else 'FAIL'}; quarks {'OK' if quark else 'FAIL'}")
    return oks


def plus_one(fn, m):
    return fn(m)


def main():
    print("=== Q4 dynamical-coupling attempt: isospin-alignment (derived) vs neutrality (candidate) ===\n")

    okA = show("Mechanism A  (frozen-I3 / ISOSPIN-alignment -- the DERIVED lepton rule):", plus_one_isospin)
    print()
    okB = show("Mechanism B  (NEUTRALITY / sterility ~ 1-|Q| -- the candidate):", plus_one_neutrality)

    print("\n  [split] A and B agree on leptons (nu is BOTH up-type and neutral), split on quarks:")
    print(f"    A -> +1 on UP quark  => delta_u=1/9, delta_d=2/27  (SWAPPED vs data: emp delta_u=2/27, delta_d=1/9)")
    print(f"    B -> +1 on DOWN quark => delta_u=2/27, delta_d=1/9  (matches data)")

    print("\n[verdict] DYNAMICAL-COUPLING ATTEMPT: FAILS to promote the |Q| selector -- and informatively.")
    print("  The ONLY *derived* dynamical mechanism is A (frozen-I3): the +1 attaches to the member whose")
    print("  channel to the RH reference preserves I3 = the same-isospin (up-type) member. A is FORCED for")
    print("  leptons (+1 on nu, matching canon) -- but extended to quarks it puts the +1 on the UP quark,")
    print("  giving the quark twists SWAPPED relative to data (delta_u<->delta_d). So the derived dynamics")
    print("  predicts {nu,u}, while the data require {nu,d}.")
    print("  The candidate B (neutrality) matches data but is NOT independently derived -- and it CONFLICTS")
    print("  with the derived mechanism A. The lepton case cannot adjudicate (nu is up-type AND neutral);")
    print("  the quark datum picks B, against the derived dynamical grain.")
    print("  NET: the +1 selector cannot be promoted dynamically. Worse than 'unforced': the substrate's")
    print("  derived coupling (frozen-I3) actively predicts the WRONG quark member. So the up/down delta")
    print("  asymmetry is a PHENOMENOLOGICAL neutrality pattern that the derived dynamics does not produce")
    print("  (and mildly contradicts). The +1 stays Proposition; both promotion routes (count + dynamics)")
    print("  are now exhausted negative. The base d=2 (condition 1) remains derived.")

    # gates: A correct on leptons, WRONG (swapped) on quarks; B correct everywhere
    assert okA["e"] and okA["nu"], "derived mechanism A must reproduce the lepton twists"
    assert not (okA["u"] or okA["d"]), "mechanism A must FAIL on BOTH quark twists (the swap)"
    assert all(okB.values()), "the neutrality candidate B must reproduce all four twists"
    assert delta("u", plus_one_isospin) == F(1, 9) and delta("d", plus_one_isospin) == F(2, 27), \
        "mechanism A gives the quark twists swapped (delta_u=1/9, delta_d=2/27)"
    print("\nGATES PASSED -- derived isospin mechanism A reproduces leptons but SWAPS the quark twists")
    print("(predicts {nu,u}); neutrality B matches data but is underived and conflicts with A. Dynamical")
    print("promotion of the +1 selector: NEGATIVE. Q4 +1 remains a phenomenological (data-only) rule.")
    print("exit 0")


if __name__ == "__main__":
    main()
