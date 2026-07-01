#!/usr/bin/env python3
r"""fermion_q4_charge_selector.py -- a candidate resolution of Q4 (the up/down delta asymmetry).

THE OPEN QUESTION (Q4, ANCHOR item 35 / DRIFT M9). The four-sector Koide twists are
   lepton e: 2/9    neutrino nu: 3/9    up u: 2/27    down d: 1/9.
Canon writes delta = d/N_eff with d the nu_R-defect count and N_eff the (colour-enlarged) plaquette.
It reproduces e, nu (d=2, d_eff=3 via the frozen I_3) and asserts up (d=2, N_eff=N*N_c=27). The DOWN
value delta_d=1/9 is the named-open Q4: canon reads it as "single-path" (d=1, N_eff=9), with the d=1
numerator "unnamed", and the nu_R mechanism is FIREWALLED on quarks (LQ=1), so the lepton d cannot be
borrowed. Prior routes (cubic-charge-trace, AB-holonomy, Atiyah-Singer, Heron) are all RETIRED.

THE CANDIDATE (this script). Read delta_d = 1/9 = 3/27, i.e. read BOTH quark triplets as uniformly
colour-shared (N_eff = N*N_c = 27 -- more natural than singling down out as single-path). Then the
four-sector DEFECT pattern is
   d(e)=2   d(nu)=3   d(u)=2   d(d)=3,
i.e. a BASE d=2 in every sector plus a +1 that lands on {nu, d}. The proposal: the +1 frozen defect
attaches to the SMALLER-|Q| (more electrically neutral) member of each weak doublet -- generalising
the lepton mechanism (the neutral nu, not e, carries the frozen extra defect) to quarks (the more
neutral d, |Q|=1/3, not u, |Q|=2/3, carries it). One named referent: electric-charge neutrality.

   delta_sector = [2 + (1 if member is the smaller-|Q| of its doublet else 0)] / N_eff,
   N_eff = 9 * (N_c=3 if coloured else 1).

THE TEST (the honest part). A rule that reproduces 4 numbers is only evidence if the SELECTOR is not
one of many. So we enumerate every natural within-doublet selector (I_3 sign, |Q| magnitude, |Y|,
Q^2, up/down-type) and report which ones pick exactly the empirical set {nu, d}. Within a weak doublet
hypercharge Y is EQUAL (doublet structure), so |Y| cannot discriminate at all; I_3 discriminates by
SIGN (giving up-type {nu,u} or down-type {e,d}); only |Q| / Q^2 discriminate by MAGNITUDE and give
{nu, d}. If |Q|-magnitude is the UNIQUE clean discriminator yielding {nu,d}, the asymmetry has a
single named referent.

HONEST SCOPE: this is a candidate (Proposition tier), NOT a closure. It (a) re-reads delta_d with
N_eff=27 (a choice, defended by triplet uniformity); (b) needs the BASE d=2 on quarks to be a
quark-native referent (2 single-flip channels to the RH quark singlet, the firewall-free analogue of
the lepton's 2 nu_R channels) -- asserted here, not proven; (c) does NOT address the separate down-R
anomaly (R_d = sqrt(12/5) ~ 1.549, "no colour story") or the absolute scale (No-Go-1).
"""
import numpy as np
from fractions import Fraction as F

# weak-doublet members: (sector, doublet, is_coloured, I3 [+/-1/2], Q, Y_LH)
MEM = {
    "e":  dict(doublet="lepton", coloured=False, I3=F(-1, 2), Q=F(-1, 1), Y=F(-1, 1)),
    "nu": dict(doublet="lepton", coloured=False, I3=F(+1, 2), Q=F(0, 1),  Y=F(-1, 1)),
    "u":  dict(doublet="quark",  coloured=True,  I3=F(+1, 2), Q=F(+2, 3), Y=F(+1, 3)),
    "d":  dict(doublet="quark",  coloured=True,  I3=F(-1, 2), Q=F(-1, 3), Y=F(+1, 3)),
}
DOUBLETS = {"lepton": ["nu", "e"], "quark": ["u", "d"]}
TARGET = {"e": F(2, 9), "nu": F(3, 9), "u": F(2, 27), "d": F(1, 9)}  # the empirical twists
PLUS_ONE_SET = {"nu", "d"}                                            # who carries the +1 (d=3)


def n_eff(m):
    return 9 * (3 if MEM[m]["coloured"] else 1)


def smaller_absQ_member(doublet):
    a, b = DOUBLETS[doublet]
    return a if abs(MEM[a]["Q"]) < abs(MEM[b]["Q"]) else b


def predict_delta(m):
    plus = 1 if m == smaller_absQ_member(MEM[m]["doublet"]) else 0
    return F(2 + plus, n_eff(m))


def main():
    print("=== Q4 candidate: +1 frozen defect on the smaller-|Q| doublet member; quarks colour-share ===\n")
    print("  [1] Does the rule reproduce all four twists exactly?")
    ok = True
    print(f"      {'sector':6s} {'N_eff':>5s} {'d':>2s} {'delta_pred':>10s} {'delta_obs':>10s} {'match':>6s}")
    for m in ("e", "nu", "u", "d"):
        dp = predict_delta(m)
        d_num = 2 + (1 if m == smaller_absQ_member(MEM[m]["doublet"]) else 0)
        match = (dp == TARGET[m])
        ok &= match
        print(f"      {m:6s} {n_eff(m):5d} {d_num:2d} {str(dp):>10s} {str(TARGET[m]):>10s} {str(match):>6s}")
    print(f"      -> all four reproduced: {ok}\n")

    print("  [2] Is the |Q|-magnitude selector UNIQUE? (which within-doublet selectors pick {nu,d}?)")
    selectors = {
        "smaller |Q|":      lambda dbl: min(DOUBLETS[dbl], key=lambda x: abs(MEM[x]["Q"])),
        "larger  |Q|":      lambda dbl: max(DOUBLETS[dbl], key=lambda x: abs(MEM[x]["Q"])),
        "smaller Q^2":      lambda dbl: min(DOUBLETS[dbl], key=lambda x: MEM[x]["Q"] ** 2),
        "I3 = +1/2 (up-type)":   lambda dbl: max(DOUBLETS[dbl], key=lambda x: MEM[x]["I3"]),
        "I3 = -1/2 (down-type)": lambda dbl: min(DOUBLETS[dbl], key=lambda x: MEM[x]["I3"]),
        "smaller |Y|":      lambda dbl: min(DOUBLETS[dbl], key=lambda x: abs(MEM[x]["Y"])),
        "larger  |Y|":      lambda dbl: max(DOUBLETS[dbl], key=lambda x: abs(MEM[x]["Y"])),
    }
    hits = []
    for name, sel in selectors.items():
        chosen = {sel("lepton"), sel("quark")}
        # |Y| is degenerate within a doublet -> selection ill-defined; flag it
        ydeg = abs(MEM["nu"]["Y"]) == abs(MEM["e"]["Y"]) and abs(MEM["u"]["Y"]) == abs(MEM["d"]["Y"])
        note = "  (|Y| DEGENERATE within doublet -> cannot discriminate)" if "|Y|" in name and ydeg else ""
        match = (chosen == PLUS_ONE_SET)
        if match:
            hits.append(name)
        print(f"      {name:22s} -> {{{', '.join(sorted(chosen))}}}  {'== {nu,d} HIT' if match else ''}{note}")
    print(f"\n      selectors reproducing the empirical {{nu,d}}: {hits}")

    print("\n  [3] Why I_3 and Y FAIL (the structural point):")
    print("      - within a weak doublet, hypercharge Y is EQUAL by construction (lepton Y=-1 both;")
    print("        quark Y=+1/3 both) -> |Y| cannot discriminate doublet members at all.")
    print("      - I_3 discriminates by SIGN -> up-type {nu,u} or down-type {e,d}, never {nu,d}.")
    print("      - only |Q| (= electric-charge MAGNITUDE) discriminates by size and yields {nu,d}.")
    print("      So the asymmetry is genuinely a NEUTRALITY rule, NOT the obvious isospin/type split.")

    print("\n[verdict] CANDIDATE Q4 resolution (Proposition tier):")
    print("  Reading both quark triplets as colour-shared (N_eff=27) makes the four-sector defect")
    print("  pattern (e,nu,u,d)=(2,3,2,3): a BASE d=2 + a +1 on the SMALLER-|Q| (more neutral) doublet")
    print("  member. This reproduces all four twists EXACTLY with ONE named referent (electric-charge")
    print("  neutrality), generalising the lepton frozen-neutrino mechanism to quarks (down, not up,")
    print("  carries the +1). |Q|-magnitude is the UNIQUE clean within-doublet selector (I_3 gives the")
    print("  wrong sign-split; Y is doublet-degenerate). This is exactly the promotion target: a single")
    print("  rule, named referent, no borrowed lepton d (the +1 is |Q|-selected; base d=2 is the quark's")
    print("  OWN RH-singlet flip count).")
    print("  NOT a closure: (a) base d=2 on quarks (2 flips to the RH quark singlet) is asserted, needs")
    print("  the explicit firewall-free channel count; (b) N_eff=27 down-reading is a choice (defended");
    print("  by triplet uniformity); (c) the down-R anomaly (sqrt(12/5)) and the absolute scale are")
    print("  untouched. Tier: Proposition -- a candidate referent for Q4, to be promoted by the channel")
    print("  count + a §16.3 null-model check on the |Q| selector across a wider hypothesis set.")

    assert ok, "the rule must reproduce all four twists exactly"
    assert hits == ["smaller |Q|", "smaller Q^2"], "the |Q|-magnitude selector must be the unique clean discriminator giving {nu,d}"
    assert smaller_absQ_member("lepton") == "nu" and smaller_absQ_member("quark") == "d"
    print("\nGATES PASSED -- rule reproduces (2/9,3/9,2/27,1/9) exactly; |Q|-magnitude is the unique clean")
    print("within-doublet selector for {nu,d}; I_3/Y excluded. Q4 reframed to a single named referent.")
    print("exit 0")


if __name__ == "__main__":
    main()
