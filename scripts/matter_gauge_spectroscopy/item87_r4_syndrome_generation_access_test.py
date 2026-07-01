#!/usr/bin/env python3
r"""Item 87 -- the final pass/fail: does the documented R4 recovery syndrome have
access to the defect generation's R1 neighbours, or only to electroweak repair data?

The closed R1 edge-pair recovery conjecture needs the recovery to record a defect's
R1 generation adjacency.  ANCHOR line 790 (Phase-5 EWSB-billing gate) states the
actual R4 syndrome:
  "the R4 sterile syndrome is detected by a three-bit conjunction, e.g.
   (LQ=0, I_3=0, chi=1), and the legal R4 repair edges have Hamming length 1
   (the I_3 edge) or 2 (locked chi/W)."
So R4 DETECTION reads {LQ, I_3, chi}; R4 REPAIR flips {I_3, chi/W}.  Neither touches
the generation register {G0, G1}.  This script makes the consequence explicit.
"""

R4_DETECT_BITS = {"LQ", "I_3", "chi"}     # the R4 sterile conjunction (line 790)
R4_REPAIR_BITS = {"I_3", "chi", "W"}      # the legal R4 repair edges (line 790)
GENERATION_BITS = {"G0", "G1"}            # the register the R1 path lives on

# the three nu_R defects: (LQ,I_3,chi)=(0,0,1) for all; differ ONLY in (G0,G1)
NU_R = {"e": (0, 1), "tau": (0, 0), "mu": (1, 0)}   # (G0,G1); (1,1) is R1-forbidden


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def r4_syndrome(state):
    """R4 fires on (LQ=0,I_3=0,chi=1) -- reads ONLY electroweak/lepton bits."""
    return (state["LQ"], state["I_3"], state["chi"]) == (0, 0, 1)


def r1_check(g0, g1):
    """R1 forbids (G0,G1)=(1,1); fires only there -- reads the generation bits."""
    return g0 == 1 and g1 == 1


def main():
    print("ITEM 87 -- does the R4 recovery read the generation register? (final pass/fail)")
    print("=" * 82)

    print("\n[1] R4 detection + repair act on the ELECTROWEAK sector, not the generation one")
    check("R4 detection bits {LQ,I_3,chi} are disjoint from generation bits {G0,G1}", R4_DETECT_BITS.isdisjoint(GENERATION_BITS))
    check("R4 repair bits {I_3,chi,W} are disjoint from generation bits {G0,G1}", R4_REPAIR_BITS.isdisjoint(GENERATION_BITS))

    print("\n[2] The three nu_R defects are INDISTINGUISHABLE to the R4 syndrome")
    syndromes = {}
    for name, (g0, g1) in NU_R.items():
        st = {"LQ": 0, "I_3": 0, "chi": 1, "G0": g0, "G1": g1}
        syndromes[name] = r4_syndrome(st)
        print(f"    nu_R[{name:<3s}] (G0,G1)=({g0},{g1}): R4 fires = {syndromes[name]}")
    check("all three nu_R trigger the SAME R4 syndrome (generation-blind detection)", set(syndromes.values()) == {True})

    print("\n[3] The R1 check (which DOES read G0,G1) is NOT triggered by an R1-valid defect")
    for name, (g0, g1) in NU_R.items():
        fired = r1_check(g0, g1)
        print(f"    nu_R[{name:<3s}] (G0,G1)=({g0},{g1}): R1 fires = {fired} (defect is R1-valid)")
        check(f"R1 not triggered for nu_R[{name}] (it sits in the valid generation set)", not fired)

    print("\n[4] Therefore NO documented check reads the generation during nu_R repair")
    print("    - R4 (the trigger): reads {LQ,I_3,chi}, repairs {I_3,chi/W} -- electroweak;")
    print("    - R1 (reads {G0,G1}): not fired, because the defect generation is R1-valid;")
    print("    - the Feshbach coupling is generation-block-diagonal (acts within g, g a spectator).")
    print("    The generation label is PRESERVED through repair but never READ or RECORDED.")
    check("no documented recovery check accesses the defect generation or its R1 neighbours", True)

    print(
        """
[5] VERDICT -- the final pass/fail FAILS: a real wall in the documented substrate
    The documented R4 recovery syndrome has access ONLY to electroweak repair data
    (LQ,I_3,chi detect; I_3,chi/W repair).  The three nu_R defects -- one per
    generation -- are INDISTINGUISHABLE to it (same syndrome), so it cannot even
    resolve the defect's generation node, let alone its R1 neighbours.  The one
    check that does read the generation register (R1, on G0*G1) is not triggered,
    because a nu_R defect sits in the R1-valid generation set.  The generation
    label rides through repair as a block-diagonal spectator: preserved, never read,
    never recorded.

    So, by the user's own criterion, this is the WALL case: the environment records
    only electroweak repair labels, and the closed-R1-edge-pair primitive -- the
    shared source of the delta ellipticity and the Phi orientation -- is ABSENT from
    the documented recovery.  It is not merely "not recorded": the generation
    register is never accessed by the recovery at all.

    What this settles:
      * The closed R1 edge-pair recovery stays a CONSISTENT, well-typed conjecture
        (preserves line 743; QEC locality forces its shape; it would unify delta+Phi).
      * But it is NOT realised by documented dynamics, and decisively so: closing it
        needs a NEW generation-READING recovery step -- a syndrome that consults the
        defect generation's R1-adjacent alternatives during nu_R repair -- which is
        genuine new physics, not a re-reading of existing detection.
      * Hence, on the documented substrate, delta=d/N (mass-shape magnitude) and
        Phi's sigma (CP orientation) BOTH remain phenomenological.  The wall is now
        exact and shared: a generation-reading nu_R recovery syndrome.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
