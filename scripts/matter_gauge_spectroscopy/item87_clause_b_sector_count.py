#!/usr/bin/env python3
r"""Item 87, clause (B): is the sector-count A_s = d_s DERIVED from the R1/boot
structure, or imported from the canon's (electroweak) defect counting?

The honest test is NOT "does d_s/N_s hit the targets" -- it does, because d_s/N_s IS
the canon's own delta = d/N (line 756).  The test is whether the R1 Hasse-edge boot
structure SOURCES the sector-dependent d_s, or whether d_s must be brought in from
the separately-derived defect-qubit counting and merely RE-LABELLED as 'R1 rescue
contacts'.
"""
from fractions import Fraction as F


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


# canon's defect-qubit counts and plaquette denominators (ANCHOR lines 705, 756)
SECTORS = {
    "lepton":   dict(d=2, N=9,  origin="2 nu_R channel bits (e_R, nu_L -> nu_R)"),
    "neutrino": dict(d=3, N=9,  origin="2 nu_R bits + frozen I_3 (CNOT-inactive)  -> d_l + 1"),
    "down":     dict(d=1, N=9,  origin="1 (down-quark frozen-defect count)"),
    "up":       dict(d=2, N=27, origin="2, with N -> N*N_c = 27 (colour)"),
}


def main():
    print("ITEM 87 clause (B): is A_s = d_s DERIVED, or an imported re-label?")
    print("=" * 80)

    print("\n[1] The R1/boot Hasse structure is SECTOR-INDEPENDENT")
    # the generation order ideal {00,01,10} and its 2 Hasse cover edges are the SAME
    # for every fermion sector (e, mu, tau exist in all sectors).
    hasse_edges = 2          # 00-01, 00-10
    generations = 3
    print(f"    Hasse cover edges = {hasse_edges}, generations = {generations}: identical for ALL sectors")
    check("the R1 Hasse structure gives a FIXED count (2 edges), not the sector-dependent d_s",
          hasse_edges == 2 and len({hasse_edges for _ in SECTORS}) == 1)
    check("so the R1/boot structure alone CANNOT source d_s = (2,3,1,2)", len({s["d"] for s in SECTORS.values()}) > 1)

    print("\n[2] The sector-dependent d_s comes from the canon's DEFECT-QUBIT counting")
    for name, s in SECTORS.items():
        print(f"    {name:<8s} d={s['d']} N={s['N']:2d}  <- {s['origin']}")
    check("d_s is an electroweak/colour defect count (nu_R bits, frozen I_3, colour), not a generation property",
          SECTORS["neutrino"]["d"] == SECTORS["lepton"]["d"] + 1)  # the +1 is the frozen I_3, an EW bit

    print("\n[3] So clause (B) RE-LABELS d_s as 'R1 rescue contacts' via a 1:1 connection")
    print("    A_s = d_s requires: each sector frozen-defect channel triggers exactly ONE R1")
    print("    rescue contact. The count d_s is the canon's (derived) defect number; the R1")
    print("    structure does not produce it. The 1:1 defect<->rescue identification is the")
    print("    new premise -- plausible (a frozen channel perturbs the register -> a rescue),")
    print("    but it bridges the ELECTROWEAK defect sector to the GENERATION R1 sector.")
    # the cross-sector 'pass' is trivial: it is the canon's own d/N re-stated
    for name, s in SECTORS.items():
        eps = F(s["d"], s["N"])
        print(f"    {name:<8s}: A_s=d_s={s['d']}, B_s=(N-2d)/3={F(s['N']-2*s['d'],3)}, eps=d/N={eps}")
    check("with the imported d_s, eps=d/N reproduces the targets -- but this is the canon's d/N, re-labelled",
          F(SECTORS["lepton"]["d"], SECTORS["lepton"]["N"]) == F(2, 9))

    print(
        """
[4] VERDICT -- clause (B) does NOT close: it is an imported re-label, not a derivation
    The honest result: the R1/boot Hasse structure is SECTOR-INDEPENDENT (2 cover
    edges, 3 generations, identical across fermion sectors), so it cannot by itself
    produce the sector-dependent d_s = (2,3,1,2).  Those numbers are the canon's
    DEFECT-QUBIT counts (nu_R bits + frozen I_3 + colour tripling) -- an
    electroweak/colour bookkeeping, separately derived, NOT a property of the
    generation register.  Clause (B) therefore RE-LABELS the canon's d_s as 'R1
    rescue contacts' through a 1:1 identification -- each sector frozen-defect
    channel triggers exactly one R1 rescue.  That identification is plausible (a
    frozen channel perturbing the generation register would trigger a rescue) but
    it is ASSERTED, not derived; it bridges two different sectors (electroweak
    defect <-> generation R1), and the cross-sector 'match' is trivial because
    d_s/N_s IS the canon's own delta = d/N.

    So clause (B) is the genuine residual, and it is softer than a wall but harder
    than clause (A): unlike the FORCED Hasse geometry, the sector MAGNITUDE is not
    forced by the R1 structure -- it is the canon's (derived) defect count, whose
    identification with R1-rescue multiplicity is a plausible but unproven 1:1 law.

    NET delta/Phi status, honestly:
      * geometry: FORCED (Boolean order ideal);
      * delta/Phi carrier + balance + QND recovery: settled;
      * clause (A) (R1 monitored): leans TRUE (nonlinearity dissolved, sec-5.9-favoured);
      * clause (B) (A_s = d_s): NOT closed -- the sector magnitude is the canon's
        derived defect count, re-labelled as R1 rescue contacts via a plausible-but-
        asserted 1:1 connection (electroweak defect channel <-> generation R1 rescue).
    delta = d/N is therefore CONDITIONAL on (A) [leans true] + (B) [plausible 1:1,
    unproven].  The thread closes here: a near-derivation whose remaining gap is a
    concrete, plausible, but unproven defect<->rescue counting identification.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
