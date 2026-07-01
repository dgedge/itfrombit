#!/usr/bin/env python3
r"""fermion_q4_selector_nullcheck.py -- §16.3 null-model check on the Q4 |Q|-neutrality selector.

Promotion condition (2) for the Q4 charge-neutrality candidate (DRIFT M9 cluster): is the
"+1 defect on the smaller-|Q| (more neutral) doublet member" selector FORCED, or is it one alphabet
choice among many that happen to reproduce the empirical set {nu, d}? The §16.3 discipline: a match is
evidence only if the controlled hypothesis space makes it rare.

THE HARD STRUCTURAL FACT. A "selector" picks one member from each of the two weak doublets
{nu,e} and {u,d}. That is exactly 2^2 = 4 possible outcomes: {nu,u},{nu,d},{e,u},{e,d}. So a RANDOM
selector reproduces the empirical {nu,d} with probability 1/4 = 25%. Reproducing {nu,d} therefore
CANNOT be "rare" -- the outcome space is irreducibly four-fold. This already caps how forced any
selector can be.

THE CONTROLLED ALPHABET. We enumerate the natural per-member scalar quantum numbers and use each as a
min/max selector, tallying which of the 4 outcomes each gives, and in particular how MANY natural
selectors reproduce {nu,d}. If {nu,d} is hit by a whole class (not a lone selector), |Q| is not
distinguished. We also do the sharper, fairer test: CONDITION on the rule reproducing the DERIVED
lepton assignment (+1 on the neutral nu) -- since the lepton +1 is canon-derived -- and ask how the
surviving selectors split on the quark doublet.

EXPECTED HONEST OUTCOME. (i) {nu,d} is reproduced by the whole charge-magnitude class, not a lone
selector -> not rare. (ii) Crucially nu is BOTH the up-type (I3=+1/2) AND the neutral (Q=0) lepton, so
the derived lepton case CANNOT distinguish a neutrality rule (|Q|) from an isospin rule (I3); they are
DEGENERATE on leptons and split only on the quark datum (|Q|->d, correct; I3->u, wrong). So the
selector is data-selected, not a-priori forced. (iii) Canon's *derived* lepton mechanism is the
dynamics-based frozen-I3, which is firewalled on quarks (LQ=1 fires the CNOT on both -> no asymmetry,
DRIFT M9), so |Q| is a NEW value-based rule, not that mechanism extended.

VERDICT (honest): condition (2) does NOT force |Q|. |Q|-neutrality is the physically-favoured choice
(the +1 tracks neutrality/sterility; the isospin alternative would put the quark +1 on the *charged*
up quark, against that intuition), but it is NOT count-forced. The +1 selector stays Proposition; the
base d=2 (condition 1) remains derived. The candidate is half-promoted, honestly.
"""
from fractions import Fraction as F
from itertools import product

# LH weak doublets; hypercharge Y_LH, B-L
M = {
    "nu": dict(dbl="L", Q=F(0), I3=F(1, 2),  Y=F(-1), BL=F(-1)),
    "e":  dict(dbl="L", Q=F(-1), I3=F(-1, 2), Y=F(-1), BL=F(-1)),
    "u":  dict(dbl="Q", Q=F(2, 3), I3=F(1, 2),  Y=F(1, 3), BL=F(1, 3)),
    "d":  dict(dbl="Q", Q=F(-1, 3), I3=F(-1, 2), Y=F(1, 3), BL=F(1, 3)),
}
DBL = {"L": ["nu", "e"], "Q": ["u", "d"]}
TARGET = frozenset({"nu", "d"})

# per-member scalar quantum numbers to try as selectors
SCALARS = {
    "Q (signed)": lambda m: M[m]["Q"],
    "|Q|":        lambda m: abs(M[m]["Q"]),
    "Q^2":        lambda m: M[m]["Q"] ** 2,
    "I3 (signed)": lambda m: M[m]["I3"],
    "|I3|":       lambda m: abs(M[m]["I3"]),
    "Y":          lambda m: M[m]["Y"],
    "|Y|":        lambda m: abs(M[m]["Y"]),
    "B-L":        lambda m: M[m]["BL"],
}


def discriminates(scalar, dbl):
    a, b = DBL[dbl]
    return SCALARS[scalar](a) != SCALARS[scalar](b)


def pick(scalar, dbl, mode):
    a, b = DBL[dbl]
    f = SCALARS[scalar]
    return (min if mode == "min" else max)(DBL[dbl], key=f)


def selector_outcome(scalar, mode):
    # ill-defined if it fails to discriminate in either doublet
    if not (discriminates(scalar, "L") and discriminates(scalar, "Q")):
        return None
    return frozenset({pick(scalar, "L", mode), pick(scalar, "Q", mode)})


def main():
    print("=== §16.3 null-model check on the Q4 |Q|-neutrality selector ===\n")
    print("  Outcome space (one member per doublet): {nu,u},{nu,d},{e,u},{e,d}  -> 4 outcomes.")
    print("  RANDOM selector reproduces the empirical {nu,d} with probability 1/4 = 25%.")
    print("  => reproducing {nu,d} is structurally NOT rare; the space is irreducibly four-fold.\n")

    print("  Controlled alphabet of natural scalar selectors (min/max per doublet):")
    tally = {frozenset(s): [] for s in [("nu", "u"), ("nu", "d"), ("e", "u"), ("e", "d")]}
    ndefined = 0
    for scalar in SCALARS:
        for mode in ("min", "max"):
            out = selector_outcome(scalar, mode)
            if out is None:
                continue
            ndefined += 1
            tally[out].append(f"{mode} {scalar}")
    for s, sels in tally.items():
        label = "{" + ",".join(sorted(s)) + "}"
        star = "   <-- empirical" if s == TARGET else ""
        print(f"    {label:9s}: {len(sels)} selectors  {sels}{star}")
    hits = len(tally[TARGET])
    print(f"\n  natural selectors reproducing {{nu,d}}: {hits} of {ndefined} well-defined "
          f"({100*hits/ndefined:.0f}%).  Each of the 4 outcomes is hit by a whole CLASS, not a lone rule.")

    print("\n  [degeneracy] nu is BOTH up-type (I3=+1/2) AND neutral (Q=0). So among selectors that")
    print("  correctly give nu in the LEPTON doublet (the canon-DERIVED case), how do they split on quarks?")
    for scalar in SCALARS:
        for mode in ("min", "max"):
            out = selector_outcome(scalar, mode)
            if out is None or "nu" not in out:
                continue
            q = next(iter(out - {"nu"}))
            verdict = "-> {nu,d} (neutrality, CORRECT)" if q == "d" else "-> {nu,u} (isospin/charge-sign, wrong on quarks)"
            print(f"    {mode} {scalar:11s} selects nu  {verdict}")
    print("  => neutrality (|Q|/Q^2) and isospin (I3) are DEGENERATE on leptons; only the QUARK datum")
    print("     splits them ({nu,d} vs {nu,u}). The lepton derivation cannot fix the selector.")

    print("\n  [mechanism] canon's *derived* lepton +1 is the dynamics-based frozen-I3 rule, which is")
    print("  FIREWALLED on quarks (LQ=1 fires the CNOT on both up and down -> no asymmetry; DRIFT M9).")
    print("  So |Q|-neutrality is a NEW value-based rule, not the derived mechanism extended.")

    print("\n[verdict] CONDITION (2) NOT PASSED -- the |Q| selector is NOT forced by §16.3:")
    print("  (i) the outcome space is four-fold (25% random baseline) -- {nu,d} cannot be rare;")
    print("  (ii) it is reproduced by the whole charge-magnitude CLASS, not a lone selector;")
    print("  (iii) neutrality and isospin are lepton-degenerate (nu is neutral AND up-type), split only")
    print("        by the quark datum -- so |Q| is DATA-SELECTED, not a-priori derived;")
    print("  (iv) the derived lepton mechanism (frozen-I3) is firewalled on quarks, so |Q| is not a")
    print("       derivation from it. |Q| remains the physically-FAVOURED choice (the +1 tracks")
    print("       neutrality/sterility; isospin would put the quark +1 on the charged up quark), but")
    print("       FAVOURED != FORCED. The +1 selector stays Proposition.")
    print("  NET Q4: base d=2 DERIVED (condition 1); +1 selector PROPOSITION (condition 2 fails to")
    print("  force). The candidate is half-promoted -- the numerator's base is a theorem, its +1 is a")
    print("  physically-motivated, data-selected hypothesis. Honest ceiling.")

    # gates -- these assert the HONEST (negative) structure, not a forced selector
    assert hits >= 2, "the charge-magnitude class (>=2 natural selectors) reproduces {nu,d} -> not rare"
    # neutrality and isospin both select nu on leptons but differ on quarks (the degeneracy)
    assert selector_outcome("|Q|", "min") == TARGET, "min|Q| must give {nu,d}"
    assert selector_outcome("I3 (signed)", "max") == frozenset({"nu", "u"}), "max I3 must give {nu,u} (the lepton-degenerate competitor)"
    assert "nu" in selector_outcome("|Q|", "min") and "nu" in selector_outcome("I3 (signed)", "max"), \
        "both |Q|-min and I3-max select nu on leptons (the degeneracy that the lepton case cannot break)"
    print("\nGATES PASSED -- {nu,d} reproduced by a class (not rare); |Q|-neutrality and I3-isospin are")
    print("lepton-degenerate competitors split only by the quark datum. Condition (2) honestly NOT forced.")
    print("exit 0")


if __name__ == "__main__":
    main()
