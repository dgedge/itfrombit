#!/usr/bin/env python3
r"""FINISH the E2 foundational-circularity audit (sec 2.8 charges / item 99 / item 133).

E2 (2026-05-30) left: GREEN-AMBER on sec 2.8 (the colour one-hot normalization rests on
a flagged-open step), GREEN+reframe on item 99, AMBER on item 133 (numerator-48
un-derived; 3-generations conditional on 3D). The stake, per the outside-auditor pass:
if the foundations were reverse-engineered to their answers, every downstream
"derivation" is a consistency check. This script finishes all three legs with zero-knob
computations and renders the root-level verdict.

LEG A — sec 2.8: the one-hot normalization is ANOMALY-FORCED (closes the soft spot).
LEG B — item 133: the covering map is REFUTED in its natural formalization (O_h has no
        (Z2)^4 subgroup to carry the code's logical group); the quotient does not track
        the register content under the R1 counterfactual; the operative origin of "3"
        is the disclosed axiom R1 (= 2^2 - 1). Reflections = alpha/beta swap (computed)
        grounds 48 = 24 x 2 only as cardinality, not as a covering.
LEG C — item 99: harness re-verified (4 flat trivalent vertices; {4,f,f} flatness gives
        f = 8 uniquely); E2's re-credit standing checked textually (reported, not edited).
LEG D — DISCLOSED tension later resolved by axis_booking_resolution.py: sec 2.4 used
        axis wording for colour while item 133/E2 used axis wording for generations.
        The resolution is category separation, not a shared axis theorem.

Self-asserting; exit 0 = every number in the prose verified."""
import itertools as it
from fractions import Fraction as Fr
import numpy as np

# =================== LEG A: the colour one-hot is anomaly-forced ===================
# Framework inputs (E2-verified, upstream): Z_f (isospin sign), Z_p (CPT tag) anchor the
# LEPTON charges nu = 0, e = -1, i.e. Y_lep = -1/2 in the GN split Q = T3 + Y with
# T3 = Z_f/2. ALSO framework content: 16 Weyl per generation incl. nu_R (item 116).
# Anomaly conditions then determine the quark hypercharges with NO colour convention:
Yl = Fr(-1, 2)                                        # lepton doublet (nu=0, e=-1 anchors)
Ye, Ynu = Fr(-1), Fr(0)                               # e_R, nu_R singlets
Yq = -Yl / 3                                          # [SU(2)]^2 U(1): 3 Yq + Yl = 0
assert Yq == Fr(1, 6)
# [SU(3)]^2 U(1): 2 Yq - Yu - Yd = 0
S = 2 * Yq                                            # Yu + Yd
# [grav]^2 U(1): 6Yq - 3Yu - 3Yd + 2Yl - Ye - Ynu = 0 (check it is automatic)
assert 6 * Yq - 3 * S + 2 * Yl - Ye - Ynu == 0
# [U(1)]^3: 6Yq^3 - 3(Yu^3+Yd^3) + 2Yl^3 - Ye^3 - Ynu^3 = 0  ->  Yu^3+Yd^3:
C = (6 * Yq**3 + 2 * Yl**3 - Ye**3 - Ynu**3) / 3
Pprod = (S**3 - C) / (3 * S)                          # from u^3+v^3 = S^3 - 3PS
disc = S * S - 4 * Pprod
import math
r = Fr(math.isqrt(disc.numerator), math.isqrt(disc.denominator))
assert r * r == disc                                  # exact square -> rational roots
Yu, Yd = (S + r) / 2, (S - r) / 2
assert {Yu, Yd} == {Fr(2, 3), Fr(-1, 3)}
print("LEG A — the one-hot normalization is ANOMALY-FORCED, closing sec 2.8's soft spot:")
print(f"   inputs: lepton anchors Y_lep = {Yl} (from Z_f/Z_p — E2-verified upstream), nu_R in content.")
print(f"   [SU2]^2 Y  =>  Y_q = {Yq} (unique);  [SU3]^2 Y => Yu+Yd = {S};  [U(1)]^3 => {{Yu,Yd}} = {{{Yu},{Yd}}}")
print(f"   (exact rational solve; [grav]^2 Y holds automatically). In sec 2.8's variables")
print(f"   Y = (1/3) Sum Z_c + Z_p/2, so the anomaly-forced Y values demand Sum Z_c = -3")
print(f"   (leptons) and -1 (quarks) — EXACTLY the one-hot eigenvalue pattern. The flagged-")
print(f"   open F2<->one-hot dual mapping is BYPASSED for charge bookkeeping: any convention")
print(f"   must reproduce the anomaly-forced values, and one-hot is that unique convention")
print(f"   (canon's 2026-06-01 note had separately solved the local reconstruction; the")
print(f"   residual openness concerns SMG colour GAUGING only). sec 2.8: GREEN-AMBER -> GREEN.")
for s_lep, s_q in [(-3, -1)]:
    assert Fr(1,3)*s_lep + Fr(1,2) == Yl and Fr(1,3)*s_q + Fr(1,2) == Yq

# =================== LEG B: item 133's covering map, decisively tested ===================
# O_h as the 48 signed permutation matrices.
Oh = []
for perm in it.permutations(range(3)):
    for signs in it.product((1, -1), repeat=3):
        M = np.zeros((3, 3), dtype=int)
        for j in range(3):
            M[j, perm[j]] = signs[j]
        Oh.append(M)
assert len(Oh) == 48
key = lambda M: M.tobytes()
# (i) which O_h elements swap the alpha/beta sublattices (Hamming-parity octant classes)?
def parity_class_swap(M):
    swaps = set()
    for v in range(8):
        c = np.array([2*((v >> 2) & 1)-1, 2*((v >> 1) & 1)-1, 2*(v & 1)-1])
        cp = M @ c
        w = sum(1 for x in cp if x == 1)              # Hamming weight of image octant
        w0 = sum(1 for x in c if x == 1)
        swaps.add((w0 % 2) != (w % 2))
    assert len(swaps) == 1                            # uniform over octants
    return swaps.pop()
def prod_signs(M):
    return int(round(np.prod([M[i, np.argmax(np.abs(M[i]))] for i in range(3)])))
preserve = [M for M in Oh if not parity_class_swap(M)]
for M in Oh:                                          # swap <=> product of signs = -1 (NOT det!)
    assert parity_class_swap(M) == (prod_signs(M) == -1)
assert len(preserve) == 24
assert any(round(np.linalg.det(M)) == -1 for M in preserve)        # T_d: contains mirrors
assert any(round(np.linalg.det(M)) == +1 and parity_class_swap(M) for M in Oh)  # C4 swaps
print("\nLEG B — item 133 (N_g = |O_h|/2^k):")
print("   (i) computed: the alpha/beta swap criterion is prod(signs) = -1, NOT det = -1.")
print("       The sublattice-preserving subgroup has order 24 and CONTAINS improper")
print("       elements (transposition mirrors) while EXCLUDING proper C4 rotations —")
print("       it is T_d, not the rotation group O. (Corrects the 'reflections distinguish")
print("       alpha/beta' gloss: the swap coset is O_h \\ T_d.) So the bipartite reading")
print("       grounds the numerator only as the CARDINALITY 48 = |T_d| x 2, and:")
# (ii) the natural covering formalization fails: no (Z2)^4 in O_h
def subgroup_closure(gens):
    G = {key(np.eye(3, dtype=int))}
    mats = {key(np.eye(3, dtype=int)): np.eye(3, dtype=int)}
    frontier = list(gens)
    while frontier:
        g = frontier.pop()
        for hk in list(mats):
            for prod in (g @ mats[hk], mats[hk] @ g):
                if key(prod) not in G:
                    G.add(key(prod)); mats[key(prod)] = prod; frontier.append(prod)
    return list(mats.values())
# Sylow-2 = stabilizer of the z-axis line: signed perms fixing {+-e_z}
syl = [M for M in Oh if abs(M[2, 2]) == 1 and M[2, 0] == M[2, 1] == 0 and M[0, 2] == M[1, 2] == 0]
assert len(syl) == 16                                 # a Sylow-2 subgroup (index 3)
invol = sum(1 for M in syl if np.array_equal(M @ M, np.eye(3, dtype=int)) and not np.array_equal(M, np.eye(3, dtype=int)))
has_order4 = any(not np.array_equal(np.linalg.matrix_power(M, 2), np.eye(3, dtype=int))
                 and np.array_equal(np.linalg.matrix_power(M, 4), np.eye(3, dtype=int)) for M in syl)
nonabelian = any(not np.array_equal(A @ B, B @ A) for A in syl for B in syl)
print(f"   (ii) Sylow-2 of O_h (order 16, all conjugate by Sylow II): involutions = {invol}")
print(f"        (a (Z2)^4 needs 15), contains order-4 elements: {has_order4}, non-abelian:")
print(f"        {nonabelian}. The code's logical group is (Z2)^4 — abelian, exponent 2 —")
print(f"        so O_h contains NO subgroup isomorphic to it: the 'codebook faithfully")
print(f"        tiles O_h' covering CANNOT be a subgroup/coset structure carrying the 16")
print(f"        logical states. The covering map is REFUTED in its natural formalization;")
print(f"        |O_h|/2^k is cardinality division, not group theory.")
assert invol == 11 and has_order4 and nonabelian
# (iii) the R1 counterfactual: the quotient does not track the register content
def count_states(r1_on):
    n = 0
    for c in it.product((0, 1), repeat=8):
        G0, G1, LQ, C0, C1, I3, chi, W = c
        if r1_on and (G0 and G1):
            continue
        if W != chi or ((LQ == 0) != ((C0, C1) == (0, 0))):
            continue
        n += 1
    return n
n_with, n_without = count_states(True), count_states(False)
assert (n_with, n_without) == (48, 64)
print(f"   (iii) R1 counterfactual: with R1 the register holds {n_with} = 3 x 16 states (3")
print(f"        G-configs = 2^2 - 1, canon's OWN sec 2.2 gloss: R1 'limits to three")
print(f"        generations'); without R1, {n_without} = 4 x 16 — yet |O_h|/16 = 3 regardless.")
print(f"        The quotient does NOT track the register content: the operative origin of")
print(f"        '3' is the AXIOM R1, and the quotient is a numerical coincidence")
print(f"        (|O_h| = 48 = #codewords is order-equality of two structurally unrelated")
print(f"        quantities: 6x8 vs 3x16).")
# (iv) E2's Platonic leg, reproduced: only the cube gives an integer
plat = {"T_d (tetra)": 24, "O_h (cube)": 48, "I_h (icosa)": 120}
for nm, o in plat.items():
    print(f"        {nm:14s} |G|/16 = {Fr(o,16)}")
assert [o % 16 == 0 for o in plat.values()] == [False, True, False]

# =================== LEG C: item 99 harness re-verified ===================
flat = [v for v in [(4,8,8),(6,6,6),(3,12,12),(4,6,12)] if Fr(1) - Fr(3,2) + sum(Fr(1,f) for f in v) == 0]
assert len(flat) == 4                                 # flatness alone does NOT select 4.8.8
f48 = [f for f in range(3, 100) if Fr(1) - Fr(3,2) + Fr(1,4) + Fr(2,f) == 0]
assert f48 == [8]                                     # within {4,f,f}, flatness -> f=8 unique
print("\nLEG C — item 99 harness re-verified: 4 flat trivalent vertex figures exist")
print("   (4.8.8, 6.6.6, 3.12.12, 4.6.12) — flatness is necessary-not-sufficient (E2's")
print("   correction stands); within the C4-forced {4,f,f} family, flatness gives f = 8")
print("   UNIQUELY. The load-bearing selector is Lemma 3.1's code-structural C4 plaquette")
print("   (from the 4+4 bit partition — a disclosed INPUT), then flatness filters. Item 99")
print("   stays GREEN as a derivation-from-disclosed-postulates.")

# =================== LEG D: historical axis-booking tension ===================
print("""
LEG D — historical disclosed tension (axis-booking): sec 2.4 states 'the colour-axis <->
   Q3-orientation identification (c in {0,1,2} -> x, y, z bipyramid axis) is a structural
   prediction'; item 133/E2 reads the three generations as the three spatial axes
   (N_g = |O_h : D_4h| = #axes). The SAME three axes are booked for COLOUR and for
   GENERATION with no stated mechanism distinguishing the two triplications. Until one
   leg is re-grounded, at most one of the two axis-identifications can be structural.

   Superseding note (2026-06-10): python_code/axis_booking_resolution.py closes this
   as a category split. Colour's 3 is the internal one-hot orbit; generation's 3 is
   R1; spatial rank 3 is TCH/Z^3. No shared axis theorem is adopted.

=================================================================================
E2 VERDICT — FINISHED (the root-level reverse-engineering question, answered):
  sec 2.8 charges   GREEN (upgraded): weights forced by 3 disclosed inputs (E2) AND the
                    one-hot normalization now ANOMALY-FORCED (this audit) — the charge
                    formula is a genuine derivation; reverse-engineering DISCONFIRMED.
  item 99 (4.8.8)   GREEN (E2 re-credit stands): derived from disclosed postulates via
                    the C4 lemma; flatness is the auxiliary filter. NOT reverse-engineered
                    (the bit partition is an upstream input, used openly).
  item 133 (3 gen)  DEMOTED to consistency observation: the covering map is refuted in
                    its natural formalization (no (Z2)^4 in O_h — computed); the quotient
                    does not track the register content (R1 counterfactual); the bipartite
                    split is T_d-vs-coset (computed; corrects 'reflections = swap') and
                    grounds only the cardinality 48 = |T_d| x 2. The
                    operative origin of '3' is the DISCLOSED AXIOM R1 (3 = 2^2 - 1),
                    formerly geometrically motivated by the conditional axis picture,
                    now superseded by the category split. This is the ONE confirmed instance of the post-hoc pattern
                    at the foundations — localized to item 133's framing, NOT the axioms.
  ROOT AXIOMS       The register, octant map, R1-R4, and bit partition are DISCLOSED
                    POSTULATES (inputs allowed to be inputs); the damaging scenario
                    ('axioms secretly tuned, all derivations consistency checks') is
                    DISCONFIRMED for charges and 4.8.8, CONFIRMED-localized for item
                    133's quotient framing only. Generation count: INPUT, not derived.
=================================================================================""")
print("ALL ASSERTIONS PASSED — every number above is verified. exit 0")
