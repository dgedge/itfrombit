#!/usr/bin/env python3
r"""Item 87 -- run at the canon's standing hypothesis (ANCHOR line 3460): is the
Koide phase delta the chiral Q^3 ANOMALY coefficient (pi-free rational), derived
via Fujikawa rather than geometry?  The lepton coincidence delta_l = 2/9 = |Q^3|
(total, = -2/9 exactly, line 351) looks promising.  The honest test of a mechanism
vs a coincidence is CROSS-SECTOR: a real anomaly origin must reproduce delta in
ALL four fermion sectors, not just charged leptons.

delta values (canon d/N, line 756):  lepton 2/9, neutrino 3/9=1/3, down 1/9, up 2/27.
We compute each sector's chiral Q^3 from the SM charges (exact, via Fraction) and
ask whether delta_sector equals that sector's OWN |Q^3|, or the per-generation
TOTAL |Q^3|, consistently.
"""
from fractions import Fraction as F


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


# SM electric charges, one generation (color multiplicity folded in for quarks)
Q3 = {
    "charged lepton (Q=-1)":   F(-1) ** 3,                # -1
    "neutrino      (Q=0)":     F(0) ** 3,                 # 0
    "up   x3 color (Q=+2/3)":  3 * F(2, 3) ** 3,          # 8/9
    "down x3 color (Q=-1/3)":  3 * F(-1, 3) ** 3,         # -1/9
}
TOTAL = sum(Q3.values())

# canon delta = d/N (defect-to-plaquette), per sector
DELTA = {
    "charged lepton (Q=-1)": F(2, 9),
    "neutrino      (Q=0)":   F(1, 3),
    "down x3 color (Q=-1/3)": F(1, 9),
    "up   x3 color (Q=+2/3)": F(2, 27),
}


def main():
    print("ITEM 87 -- is delta the chiral Q^3 anomaly coefficient? (cross-sector test)")
    print("=" * 82)

    print("\n[1] The exact per-generation total Q^3 anomaly (line 351)")
    print(f"    sum Q^3 = (-1)^3 + 3(2/3)^3 + 3(-1/3)^3 = {TOTAL}")
    check(TOTAL == F(-2, 9), "total Q^3 = -2/9 exactly")
    check(abs(F(2, 9)) == abs(TOTAL), "lepton delta=2/9 numerically equals |total Q^3| (the promising coincidence)")

    print("\n[2] Per-sector: does delta match the sector's OWN |Q^3| or the TOTAL?")
    print(f"    |total Q^3| = {abs(TOTAL)}")
    rule = {}
    for sec in DELTA:
        d = DELTA[sec]
        own = abs(Q3[sec])
        match_own = (d == own)
        match_total = (d == abs(TOTAL))
        rule[sec] = (match_own, match_total)
        tag = "OWN" if match_own else ("TOTAL" if match_total else "NEITHER")
        print(f"    {sec:<24s} delta={str(d):>5s}  |own Q^3|={str(own):>4s}  -> matches {tag}")

    print("\n[3] The decisive falsifier: a chargeless sector has zero anomaly")
    check(Q3["neutrino      (Q=0)"] == 0 and DELTA["neutrino      (Q=0)"] != 0,
          "neutrino: Q^3 = 0 but delta_nu = 1/3 != 0  -> delta CANNOT be the Q^3 anomaly")

    print("\n[4] No consistent rule across sectors")
    check(rule["charged lepton (Q=-1)"] == (False, True), "lepton matches TOTAL only (not its own -1)")
    check(rule["down x3 color (Q=-1/3)"] == (True, False), "down matches its OWN (1/9), not the total")
    check(rule["up   x3 color (Q=+2/3)"] == (False, False), "up (2/27) matches NEITHER (own 8/9, total 2/9)")
    check(rule["neutrino      (Q=0)"] == (False, False), "neutrino (1/3) matches NEITHER (own 0, total 2/9)")

    print(
        """
[5] VERDICT (the Fujikawa-anomaly route FAILS -- honest negative)
    Tested as a mechanism, "delta = chiral Q^3 anomaly coefficient" does not
    survive the cross-sector test:
      * neutrino is the clean falsifier -- it is electrically neutral, so its Q^3
        anomaly is EXACTLY 0, yet delta_nu = 1/3.  A charge-cube anomaly cannot
        produce a nonzero phase for a chargeless sector.
      * up (2/27) matches neither its own |Q^3|=8/9 nor the total 2/9.
      * the two sectors that DO match use INCONSISTENT references: the charged
        lepton matches the per-generation TOTAL (2/9), while the down quark
        matches its OWN contribution (1/9).  No single "delta = (some) Q^3" rule
        fits all four.
    So the promising delta_l = 2/9 = |Q^3_total| is an isolated COINCIDENCE of the
    charge-cube arithmetic aligning with the lepton defect count; it does not
    generalise.  The chiral anomaly (canon line 3460 hypothesis) is NOT the
    mechanism for delta.

    NET (the honest floor for delta): all three candidate mechanisms now FAIL --
    geometry/arc (documented coupling is phaseless), holonomy/winding (excluded by
    the fit; would carry 2pi/pi), and chiral anomaly (refuted cross-sector here).
    delta = d/N (defect counting) stands as a phenomenological law with NO derived
    substrate mechanism; its value is fit/count-selected.  (This is delta only --
    the CP phase Phi remains a genuine Z3 holonomy 2pi/3, a separate object.)
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
