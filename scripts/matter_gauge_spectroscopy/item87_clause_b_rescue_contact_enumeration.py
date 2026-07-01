#!/usr/bin/env python3
r"""Item 87, clause (B), concrete enumeration: count the sector-visible R1 rescue
contacts directly on the register, and test whether they equal d_s = (2,3,1,2).

The R1 rule lives on the generation register (G0,G1): forbidden corner 11, valid
ideal {00,01,10}.  An R1 rescue fires when a state drifts toward 11 and is corrected
to a legal neighbour.  The 'rescue contacts' are therefore a property of the
GENERATION register -- and that register is the SAME in every fermion sector (every
sector has the three generations).  So this enumeration tests, faithfully, whether
the R1 rescue-contact count can even be sector-dependent.
"""
import itertools


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


GEN = ["00", "01", "10", "11"]
FORBIDDEN = "11"
VALID = ["00", "01", "10"]
d_canon = {"lepton": 2, "neutrino": 3, "down": 1, "up": 2}   # canon defect-qubit counts
N_canon = {"lepton": 9, "neutrino": 9, "down": 9, "up": 27}


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    print("ITEM 87 clause (B): enumerate R1 rescue contacts on the register")
    print("=" * 80)

    print("\n[1] R1 rescue structure on the generation register (same object for all sectors)")
    # rescue edges: legal one-bit moves that rescue a drift toward the forbidden 11.
    # equivalently the valid generations one flip from 11 (vulnerable), or the Hasse
    # cover edges of the order ideal -- both give the same count of 2.
    vulnerable = [g for g in VALID if hamming(g, FORBIDDEN) == 1]      # one flip from 11
    cover_edges = [(a, b) for a, b in itertools.combinations(VALID, 2) if hamming(a, b) == 1]
    print(f"    valid generations one flip from forbidden 11: {vulnerable}  (count {len(vulnerable)})")
    print(f"    Hasse cover edges of the order ideal: {cover_edges}  (count {len(cover_edges)})")
    check("the generation R1 rescue structure has exactly 2 contacts", len(vulnerable) == 2 and len(cover_edges) == 2)

    print("\n[2] This count is SECTOR-INDEPENDENT (the generation register is shared)")
    # every fermion sector contains the same three generations -> same R1 rescue structure
    rescue_count = {s: len(cover_edges) for s in d_canon}             # 2 for every sector
    for s in d_canon:
        print(f"    {s:<8s}: R1 generation rescue contacts = {rescue_count[s]}   (canon d_s = {d_canon[s]})")
    check("R1 generation rescue contacts = 2 for every sector (sector-independent)", set(rescue_count.values()) == {2})

    print("\n[3] So the R1 rescue count does NOT reproduce d_s = (2,3,1,2)")
    matches = {s: (rescue_count[s] == d_canon[s]) for s in d_canon}
    for s in d_canon:
        print(f"    {s:<8s}: R1 rescue {rescue_count[s]} vs d_s {d_canon[s]}  -> {'match' if matches[s] else 'MISMATCH'}")
    check("lepton (2) and up (2) coincide; neutrino (3) and down (1) MISMATCH the rescue count",
          matches["lepton"] and matches["up"] and not matches["neutrino"] and not matches["down"])

    print("\n[4] The (2,3,1,2) is the canon's ELECTROWEAK/colour defect count, a different object")
    print("    lepton 2 = nu_R coupling channels; neutrino 3 = 2 nu_R + frozen I_3 (an EW bit);")
    print("    down 1; up 2 with N->N*N_c. None of these is a generation-register property.")
    check("d_s is sector-dependent; the R1 generation rescue count cannot be (it is shared)",
          len(set(d_canon.values())) > 1 and len(set(rescue_count.values())) == 1)

    print(
        """
[5] VERDICT -- clause (B) REFUTED by enumeration: R1 rescue contacts are not d_s
    Counted directly on the register, the R1 rescue structure has exactly 2 contacts
    (the two generations one flip from the forbidden corner / the two Hasse cover
    edges), and this is SECTOR-INDEPENDENT because every fermion sector shares the
    same generation register.  It therefore cannot equal the sector-dependent
    d_s = (2,3,1,2): it coincides with the lepton (2) and up (2) by accident but
    MISMATCHES the neutrino (3) and down (1).

    The numbers (2,3,1,2) are the canon's ELECTROWEAK/colour defect-qubit counts
    (nu_R channels, frozen I_3, colour tripling) -- a different, sector-dependent
    structure that is NOT a property of the generation register where R1 lives.  So
    the boot-R1 rescue mechanism, counted faithfully, supplies a sector-independent
    "2", not the magnitude law.  Clause (B) as stated -- 'sector s exposes d_s
    R1-polarised rescue contacts' -- is refuted: the R1 rescue contacts are 2 for
    every sector.

    NET (final, honest): the boot-R1 mechanism delivers delta's TYPE (CP-even
    ellipticity), SHAPE (forced Hasse geometry), and MONITORING (clause A, leans
    true) -- but NOT the sector MAGNITUDE.  The magnitude d_s is the canon's
    electroweak/colour defect count; it is NOT the R1 rescue-contact count (which is
    sector-independent), so identifying the two fails the enumeration.  delta = d/N
    therefore remains: d, N derived separately (defect/plaquette counts), the
    CP-even packaging from boot-R1, but the magnitude NOT mechanistically the R1
    rescue count.  Clause (B) does not close; the boot-R1 carrier explains the form,
    not the sector numbers.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
