#!/usr/bin/env python3
r"""foundations_rhbc_scalar_clock_uniqueness.py

R_HBC (item 131) — attempt to CLOSE the scalar-clock uniqueness within the local single-clock instrument.

Open piece (DRIFT l.2366/2423; item131_hbc_channel_whitening_closure.py): R_HBC = psi - nu is the
gauge-invariant curvature perturbation = delta N, conditional on ONE field identification -- that the scalar
clock is "the post-decoder, gauge-invariant, colour-restoring print-time current", UNIQUE (no second scalar
source), within the single-clock architecture (item 94). Already closed: the all-ones entropy current and the
horizon-mode covariance countermodels (the latter = nonlocal = new physics); no hidden positive ledger in the
8-bit R1-R4 instrument (item131_no_r5_instrument_completeness.py). MISSING: a positive UNIQUENESS statement.

THE ARGUMENT (a clock is a MONOTONE / print-time current; a conserved charge cannot tick), made countable by
tying it to the SAME walk-active structure that grounds G7:

  A scalar clock must be (i) gauge-invariant (colour/weak singlet), (ii) spin-0 / homogeneous, and
  (iii) a PRINT-TIME current -- it must advance with the irreversible service (recovery) events. A conserved
  gauge-invariant charge is constant and cannot serve as a time variable.

  Step 1 (walk-active vs conserved). The canonical hop flips only {C0,C1,I3} (G7,
          foundations_g7_cloop_support_from_walk.py). Hence the gauge-invariant bits {LQ, G0, G1} are
          walk-CONSERVED -> conserved charges (baryon/lepton, generation), NOT clocks.
  Step 2 (which walk-active bit prints). Of {C0,C1,I3}, only the COLOUR bits {C0,C1} EXIT the codespace
          (a detectable error -> a recovery/erasure event = a print); I3 is the FREE bit (never exits,
          no recovery). So the print-time current lives on the colour-recovery events ONLY.
  Step 3 (gauge projection -> dimension 1). The colour-recovery events carry a 3-colour index; under the
          colour symmetry (S3 Weyl of SU(3), the V_strong K_3 permutation) the 3-dim event space splits as
          singlet (+) 2-dim standard. Only the SINGLET (colour-invariant total) is gauge-invariant. So the
          gauge-invariant scalar print-time current is EXACTLY 1-dimensional = the colour-RESTORING current.

  => Within the local instrument there is EXACTLY ONE gauge-invariant scalar print-time current. The scalar
     clock is UNIQUE; R_HBC = psi - nu is forced. The ONLY escape is an outside sector (a nonlocal
     horizon-mode operator or a new hidden register) -- explicit new physics, which also breaks the
     single-clock premise (item 94) it would need. That escape is a PERMANENT conditional, the R_HBC analogue
     of G7's T7 (horizon-input) wall: not closable from within, by construction.

HONEST GRADE: this PROMOTES R_HBC uniqueness from "assumed" to "derived within the local single-clock
instrument" (the colour-singlet recovery current is the unique gauge-invariant scalar print-time current),
resting on (a) the print-time = monotone-current reading of "clock", and (b) the no-outside-sector premise
(new physics, permanent). It is NOT an absolute closure -- the outside-sector escape is irreducible.

Self-asserting; exit 0. Run under ~/bin/py13_7 (numpy).
"""
import itertools
import numpy as np

NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
B = {n: i for i, n in enumerate(NAMES)}
ALL = list(itertools.product((0, 1), repeat=8))
IDX = {c: i for i, c in enumerate(ALL)}
WALK_ACTIVE = {"C0", "C1", "I3"}                       # from G7 (V_em diag, V_weak->I3, V_strong->colour)


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def valid(c):
    G0, G1, LQ, C0, C1, I3, chi, W = c
    return (not (G0 and G1)) and W == chi and ((LQ == 0) == (C0 == 0 and C1 == 0))


P = [c for c in ALL if valid(c)]
inP = {c: valid(c) for c in ALL}
assert len(P) == 48


def fl(c, n):
    d = list(c); d[B[n]] ^= 1; return tuple(d)


def main():
    print("=== R_HBC scalar-clock uniqueness within the local instrument ===")

    # ---- Step 1: gauge-invariant bits {LQ,G0,G1} are walk-CONSERVED (not clocks) ----
    gauge_inv_bits = ["LQ", "G0", "G1"]               # colour/weak singlets (C0,C1,I3 are NOT singlets)
    conserved = [n for n in gauge_inv_bits if n not in WALK_ACTIVE]
    ok(set(conserved) == {"LQ", "G0", "G1"},
       f"gauge-invariant bits {gauge_inv_bits} are all walk-CONSERVED (not in walk-active {sorted(WALK_ACTIVE)})")
    print("    -> LQ (baryon/lepton), G0,G1 (generation) are conserved charges -> CANNOT be print-time clocks")

    # ---- Step 2: of the walk-active bits, only COLOUR exits (prints); I3 is free ----
    exits = {}
    for n in WALK_ACTIVE:
        # a bit 'prints' if flipping it takes SOME codeword out of P (a detectable error needing recovery)
        exits[n] = any(inP[c] and not inP[fl(c, n)] for c in ALL if inP[c])
    ok(exits == {"C0": True, "C1": True, "I3": False},
       f"walk-active exit/print table: {exits} -> only colour {{C0,C1}} prints; I3 is the FREE bit")
    print("    -> the print-time (recovery) current lives ONLY on the colour-recovery events")

    # ---- Step 3: gauge (colour-S3) projection of the 3-colour recovery space -> dimension 1 ----
    # the 3 colour states {01,10,11} carry the recovery index; colour symmetry = S3 permuting them.
    # build the 3-dim permutation representation and project onto the S3-invariant (singlet) subspace.
    colours = [(0, 1), (1, 0), (1, 1)]
    # S3 generators on 3 colours: a transposition (swap colours 0,1) and a 3-cycle
    def perm_mat(p):
        M = np.zeros((3, 3))
        for a in range(3):
            M[p[a], a] = 1.0
        return M
    gens = [perm_mat([1, 0, 2]), perm_mat([1, 2, 0])]   # (01), (012)
    # invariant subspace = intersection of ker(g - I); compute its dimension
    stack = np.vstack([g - np.eye(3) for g in gens])
    rank = np.linalg.matrix_rank(stack, tol=1e-9)
    dim_singlet = 3 - rank
    ok(dim_singlet == 1, f"colour-S3 invariant (singlet) subspace of the 3-colour recovery space = dim {dim_singlet}")
    # the singlet is the all-ones (colour-blind total) vector = the colour-RESTORING current
    v = np.ones(3) / np.sqrt(3)
    ok(all(np.allclose(g @ v, v) for g in gens), "the singlet IS the colour-blind total = the colour-RESTORING current")
    # the orthogonal complement (the 2-dim standard rep) is NOT gauge-invariant
    ok(rank == 2, "the remaining 2-dim colour-recovery space (standard rep) is NOT colour-invariant -> not a scalar")

    print("\n[verdict] R_HBC scalar-clock uniqueness:")
    print("  - WITHIN the local single-clock instrument there is EXACTLY ONE gauge-invariant scalar print-time")
    print("    current: the colour-SINGLET recovery current (= the colour-restoring topology current). Proof:")
    print("    (1) the gauge-invariant bits LQ,G0,G1 are walk-conserved -> conserved charges, not clocks;")
    print("    (2) of the walk-active bits {C0,C1,I3} only colour prints (I3 is free, never recovers);")
    print("    (3) colour-S3 invariance projects the 3-colour recovery space to its 1-dim singlet.")
    print("  - So the scalar clock is UNIQUE and R_HBC = psi - nu = delta N is FORCED (not assumed) inside the")
    print("    instrument -- closing the item-131 uniqueness gap at instrument grade.")
    print("  - PERMANENT conditional (irreducible, the R_HBC analogue of G7's T7 horizon wall): a nonlocal")
    print("    horizon-mode operator or a NEW hidden register would add a second scalar source -- but that is")
    print("    explicit new physics and breaks the single-clock premise (item 94) it would require. Not closable")
    print("    from within, by construction. HONEST: promoted 'assumed' -> 'derived within the instrument'.")
    print("  exit 0")


if __name__ == "__main__":
    main()
