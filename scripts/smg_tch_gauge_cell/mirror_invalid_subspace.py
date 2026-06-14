#!/usr/bin/env python3
"""
The item-98(-adjacent) test: do the lattice MIRROR doublers map into the 208-dim INVALID
subspace that the code projection P annihilates?  If yes (and if P acts as a symmetric
interaction, not a hard projection), the code gaps the mirrors -- the SMG closure.

HONESTY / SCOPE CORRECTION (from reading ANCHOR §15 item 98):
  item 98 is the C_8-vertex-ring <-> Q_3-face INTERNAL dual-locus reconciliation (two
  presentations of where the 8 bits live on the cell). It is NOT the external-BZ <->
  internal-codeword map. That map (lattice doubler -> codeword) is NOT closed anywhere in
  the corpus (items 77 HDR / 97 K6-Bloch / 102 k.p are the relevant OPEN pieces). So this
  is a CANDIDATE-MAP test, grounded in (a) a standard lattice fact and (b) the framework's
  own R2 -- NOT a framework-closed computation. Stated plainly so it can't be misread.

The candidate map (minimal, standard, framework-consistent):
  * Standard lattice doubling fact: a doubler at a BZ corner with n components = pi has
    Dirac chirality (-1)^n relative to the origin. So doubling FLIPS the spacetime Dirac
    chirality. In the register that bit is chi (face 110; ANCHOR §2.1, g5 = sigma_y(x)I,
    = the Q_3 bipartite grading, item 116).
  * Internal gauge quantum numbers are NOT spatial: the weak charge W (face 111) and all
    flavour/colour bits are carried by the field independent of which BZ corner it sits at,
    so the spatial doubling cannot flip them.
  => the spatial mirror of a codeword = flip chi, keep W (and everything else).
     R2 is exactly "W = chi" ("chirality locks to weak charge -- left-handed weak
     interactions", §2.2). So a mirror has W != chi -> R2-VIOLATED -> invalid.

This script:
  (A) state the scope correction (no faking the open map);
  (B) compute: the 16 valid codewords (W=chi) and their 16 spatial mirrors (flip chi);
      show ALL mirrors land in the 208-dim invalid subspace, specifically R2-violating;
      and that the 48 valid + 48 mirror = 96 states are exactly the vector-like doubled
      content the code splits chiral-vs-mirror;
  (C) the LEGITIMACY caveat -- the crux: a HARD projection deleting the mirrors is itself
      Nielsen-Ninomiya-forbidden; SMG legitimacy needs R2 realised as a SYMMETRIC
      INTERACTION (the code's stabiliser dynamics), legal iff anomaly-free (the 16,
      smg_code_projection.py). Mechanism identified; dynamical realisation still open;
  (D) the SAME-chirality residual (walk_kernel_overlap's 4 even-corner species): R2 does
      NOT distinguish same-chirality species -> not addressed here; and that residual may
      be an artifact of the simplified body-diagonal kernel, not the real TCH band
      structure (item 97). Doubly open.

numpy-free; pure enumeration; self-asserting.
"""
import itertools

def Hd(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# register (octant-indexed, ANCHOR §2.1):  c = (G0,G1,LQ,C0,C1,I3,chi,W)
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
def R1(c): return not (c[G0] == 1 and c[G1] == 1)
def R2(c): return c[W] == c[CHI]
def R3(c):
    return (c[C0], c[C1]) == (0, 0) if c[LQ] == 0 else (c[C0], c[C1]) != (0, 0)
def valid(c): return R1(c) and R2(c) and R3(c)          # the 48-codeword physical subspace

allwords = [tuple(c) for c in itertools.product([0, 1], repeat=8)]
valid48  = [c for c in allwords if valid(c)]
invalid208 = [c for c in allwords if not valid(c)]
assert len(valid48) == 48 and len(invalid208) == 208

# ----------------------------------------------------------------------------------
Hd("(A) SCOPE: item 98 is the INTERNAL C_8/Q_3 dual; the external map is NOT closed")
print("  item 98  = C_8-vertex-ring <-> Q_3-face reconciliation (where the 8 bits live).")
print("  NEEDED   = BZ-corner-doubler -> codeword map (items 77/97/102 territory, OPEN).")
print("  => this is a CANDIDATE-MAP test (standard doubler chirality + R2), not framework-closed.")

# ----------------------------------------------------------------------------------
Hd("(B) Do the spatial mirrors land in the 208-dim invalid subspace?")
def mirror(c):                       # spatial doubling: flip Dirac chirality chi, keep all gauge bits
    m = list(c); m[CHI] ^= 1; return tuple(m)

mirrors = [mirror(c) for c in valid48]
assert len(set(mirrors)) == 48                                   # chi-flip is a bijection: 48 distinct
in_invalid = [m for m in mirrors if not valid(m)]
r2_violating = [m for m in mirrors if not R2(m)]
r1r3_ok      = [m for m in mirrors if R1(m) and R3(m)]
print(f"  valid codewords (W=chi)                  : {len(valid48)}")
print(f"  their spatial mirrors (flip chi, keep W) : {len(mirrors)} distinct")
print(f"  mirrors landing in the 208 invalid space : {len(in_invalid)} / 48")
print(f"  mirrors that are R2-violating (W != chi) : {len(r2_violating)} / 48")
print(f"  mirrors still R1&R3-OK (only R2 broken)   : {len(r1r3_ok)} / 48")
assert len(in_invalid) == 48 and len(r2_violating) == 48 and len(r1r3_ok) == 48
print("  -> EVERY spatial mirror is invalid, and invalid PRECISELY by R2 (W != chi).")
print("     R2 ('chirality locks to weak charge') IS the anti-mirror constraint, exactly.")

# the chiral generation + its mirror = the vector-like doubled content, split by R2
union = set(valid48) | set(mirrors)
print(f"\n  48 valid + 48 mirror = {len(union)} states = the W-doubling of the 48 (chi free, W free")
print(f"     with R1,R3 fixed). The code keeps the {len(valid48)} chiral (W=chi), discards the")
print(f"     {len(mirrors)} mirror (W!=chi). Remaining invalid: {208 - len(mirrors)} states violate R1/R3")
print("     (wrong generation / wrong colour) -- a different kind of non-state, not mirrors.")
assert len(union) == 96 and (208 - len(mirrors)) == 160

# the projection P annihilates the mirror sector (P = diag indicator of valid48 in codeword basis)
def P(c): return 1 if valid(c) else 0
assert all(P(m) == 0 for m in mirrors) and all(P(c) == 1 for c in valid48)
print("  code projection P (256->48) annihilates all 48 mirror states; P is ULTRALOCAL")
print("  (R2 is a 2-bit check on faces chi=110, W=111 -- finite range on one cell).")

# ----------------------------------------------------------------------------------
Hd("(C) THE CRUX: hard projection is N-N-forbidden; SMG needs R2 as a SYMMETRIC INTERACTION")
print("""  A naive reading -- 'P just deletes the mirrors' -- is NOT a legitimate closure:
  deleting half the doublers of a FREE lattice fermion by a kinematic projection is
  exactly what Nielsen-Ninomiya forbids (it cannot be done locally without breaking the
  chiral symmetry). The mirror cannot be PROJECTED out; it must be GAPPED.

  SMG legitimacy: R2 must be realised as a SYMMETRIC INTERACTION (the code's Z-stabiliser
  dynamics on the matter cell, §2.x; the 'internal walk' on Q_3) that dynamically gaps the
  W!=chi mirror sector while preserving the chiral gauge symmetry. That is legal IFF the
  matter content is anomaly-free incl. the global Z_16 -- which holds (the 16,
  smg_code_projection.py). And the mirror sector here is the chirality-conjugate of the
  physical 16, hence the anomaly-free SMG partner by construction.""")
# sanity: the mirror sector is the chirality-conjugate of the physical sector (1-1, chi-flipped)
assert set(mirror(m) for m in mirrors) == set(valid48)          # mirror of mirror = original
print("  CHECK: mirror(mirror(.)) = identity on the 48 -> mirror sector is the exact")
print("         chirality-conjugate of the physical generation (the SMG-able partner). OK.")
print("\n  STATUS of (C): mechanism IDENTIFIED (R2 = anti-mirror; mirrors are the SMG partner),")
print("  precondition MET (anomaly-free 16); but R2-as-symmetric-gapping (not hard projection)")
print("  is the DYNAMICAL step still to be constructed. This is the real open piece.")

# ----------------------------------------------------------------------------------
Hd("(D) The SAME-chirality residual (walk_kernel_overlap's 4 even-corner species)")
print("""  walk_kernel_overlap.py's surviving species were 4 SAME-chirality (even-corner) modes.
  R2 distinguishes opposite chiralities (W vs chi) -- it does NOT separate same-chirality
  species (they share chi). So R2/the code does NOT address that multiplicity.
  Two honest caveats:
   * that 4-fold multiplicity is a feature of the SIMPLIFIED body-diagonal kinetic kernel
     (gw_nogo's 'enhanced doubling'), NOT necessarily the real TCH band structure, which is
     the K6/O_h dispersion of item 97 (spectrum {+5,-1^5} at Gamma; unevaluated doublers).
   * resolving it would need the open corner->codeword map (items 77/97/102), the same map
     this script had to SUBSTITUTE a candidate for.
  So the same-chirality residual is doubly open: possibly not a real problem, and not
  decidable without the unowned external map.""")

# ----------------------------------------------------------------------------------
Hd("VERDICT — what the item-98 test shows")
print("""POSITIVE, GROUNDED (computed above, modulo the candidate map):
  * Under the minimal standard map (doubling flips chi, leaves W), the OPPOSITE-CHIRALITY
    mirror doublers are EXACTLY the W!=chi states, and ALL 48 land in the 208-dim invalid
    subspace -- invalid PRECISELY by R2. R2 ('left-handed weak interactions') IS the
    anti-mirror constraint; the mirror sector is the chirality-conjugate SMG partner of the
    anomaly-free 16. The projection that targets them is ULTRALOCAL (a 2-bit check).

STILL OPEN (not faked):
  1. R2-as-symmetric-gapping. A hard projection is N-N-forbidden; closure needs R2 realised
     as the code's symmetric INTERACTION dynamically gapping the mirror. Precondition met,
     construction not done. THIS is the crux.
  2. The same-chirality residual (walk_kernel_overlap's 4 species) -- not touched by R2;
     possibly an artifact of the simplified kernel; needs item 97's real band structure.
  3. The external corner->codeword map itself (items 77/97/102) -- substituted by a
     candidate here, not closed.

NET: the test moves the SMG route from 'precondition met' to 'precondition met AND the
mirror-removing operator IDENTIFIED (R2) with the mirrors shown to sit in the invalid
subspace it targets ultralocally'. It does NOT reach continuum-Locked: items (1)-(3)
remain, and (1) is exactly the dynamical SMG step. Honest direction: encouraging on the
chirality-mirror half, silent on the same-chirality half, and pinned on the dynamics.

CORRECTION to record in canon: item 98 is the INTERNAL C_8/Q_3 dual, not the external map;
the external corner->codeword map is unowned (items 77/97/102 are its nearest open pieces).""")
print("\nALL ASSERTS PASSED (computed claims hold; open pieces stated, not papered over).")
