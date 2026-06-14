#!/usr/bin/env python3
"""
CIRCULARITY AUDIT of ANCHOR §15 item 99 (4.8.8 vertex-figure uniqueness) and
item 133 (3 generations = |O_h|/2^k = 48/16). EXPLORATORY (uncommitted), 2026-05-30.

Same input-vs-derived method as audit_charge_circularity.py: separate what is FORCED
(theorem) from what is CHOSEN (input that could have been otherwise). Compute the
counterfactuals. Self-asserting where exact. Needs only fractions.
"""
from fractions import Fraction as F
def H(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)

H("ITEM 99 (4.8.8 uniqueness) -- (A) what does Gauss-Bonnet flatness ACTUALLY select?")
# combinatorial (Knill) curvature of a trivalent vertex with faces [f1,f2,f3]:
# K_v = 1 - deg/2 + sum(1/f_i) = -1/4 + sum(1/f_i) for deg=3. Flat <=> K=0.
def Kvertex(faces): return F(1) - F(3,2) + sum(F(1,f) for f in faces)
print("  trivalent-vertex curvature K = 1 - 3/2 + sum 1/f_i  (flat <=> 0):")
fam={'4.8.8  (the claim)':[4,8,8],'4.10.10 (adjacent x-wire)':[4,10,10],
     '4.6.6':[4,6,6],'6.6.6  (honeycomb)':[6,6,6],'3.12.12':[3,12,12],'4.6.12':[4,6,12]}
flat=[]
for name,f in fam.items():
    k=Kvertex(f); tag='FLAT' if k==0 else ('hyperbolic' if k<0 else 'spherical')
    if k==0: flat.append(name)
    print(f"    {name:26s} faces={str(f):12s} K={str(k):6s} {tag}")
print(f"  FLAT trivalent vertices found: {flat}")
print("  *** HONEST FINDING (corrects a draft overclaim): flatness ALONE does NOT select 4.8.8.")
print("      6.6.6 (honeycomb), 4.8.8, 3.12.12, 4.6.12 are ALL flat trivalent Archimedean")
print("      vertices. So canon's 'Gauss-Bonnet flatness selects 4.8.8' is an OVERSIMPLIFICATION.")
# what flatness DOES do: among {4,f,f} (one square + two EQUAL faces), it picks f=8 uniquely
ff=[f for f in range(3,13) if Kvertex([4,f,f])==0]
print(f"  Within the {{4,f,f}} family (one C4 + two EQUAL faces): flat <=> f in {ff} -> f=8 UNIQUE.")
print("  So the real selection = (forced C4 plaquette, Lemma 3.1, from the CODE) + (two equal")
print("  faces) + (flatness). Flatness only finishes a job the C4-plaquette structure starts.")
assert Kvertex([4,8,8])==0 and Kvertex([6,6,6])==0 and ff==[8]
print("  AUDIT: flatness is INDEPENDENTLY motivated (observed Omega_k~0) -- not reverse-")
print("  engineered. But it is NECESSARY-not-sufficient; the LOAD-BEARING selector is the C4")
print("  plaquette (Lemma 3.1), which depends on the 4+4 bit partition + 2-regularity (upstream).")

H("ITEM 99 -- (B) what is the load-bearing CHOICE? the 4+4 internal/external bit partition")
print("  The proof chain (Lemmas 3.1-3.4, Thms 3.5-3.6) is FORCED given: (i) 8 bits split as")
print("  V_int={I3,chi,W,LQ} + E_ext={G0,G1,C0,C1}; (ii) trivalent vertex (degree-3 budget);")
print("  (iii) flatness. (i) is the upstream physical assignment (which bits are matter-")
print("  registers vs routing-ports). It is PRIOR to the geometry, like §2.8's Z_f/Z_p.")
print("  AUDIT verdict item 99: GREEN. Theorem chain solid + flatness independently motivated;")
print("  the one input (bit partition) is upstream and disclosed, not tuned to the 4.8.8 answer.")

H("ITEM 133 (3 generations = |O_h|/2^k) -- (A) Platonic-cell counterfactual")
# does the integer-3 only come out for the cube, or is the cell cherry-picked?
groups={'tetrahedron T_d':24,'cube/octahedron O_h':48,'icosa/dodeca I_h':120,
        'cube ROTATIONS O':24,'tetra ROT T':12,'icosa ROT I':60}
print("  N_g = |G| / 2^k with 2^k = 16 (code logical dim):")
for name,order in groups.items():
    ng=F(order,16); print(f"    {name:22s} |G|={order:3d}  |G|/16 = {ng}  {'INTEGER' if ng.denominator==1 else 'non-integer'}{' = 3 ✓' if ng==3 else ''}")
print("  => Among Platonic FULL point groups only the cube's O_h=48 gives an integer, and it")
print("     gives exactly 3. NON-trivial: not 'any cell gives 3'. BUT note the choice below.")
assert F(48,16)==3 and F(24,16).denominator!=1 and F(120,16).denominator!=1

H("ITEM 133 -- (B) the 48-vs-24 freedom: reflections are LOAD-BEARING")
print(f"  full O_h (with reflections) = 48 -> 48/16 = {F(48,16)} = 3  (the claim)")
print(f"  rotations-only O           = 24 -> 24/16 = {F(24,16)} = 3/2 (NOT integer)")
print("  => the result REQUIRES using the full point group incl. improper operations. Using")
print("     rotations-only fails. Canon does not justify why fermion-copy-count = |O_h| rather")
print("     than |O|. This is a genuine UNMOTIVATED CHOICE (the rep-theory 'faithful covering'")
print("     gap the framework itself flags). LOAD-BEARING and currently un-derived.")

H("ITEM 133 -- (C) the denominator-choice freedom: why 2^k=16 and not 48, 256, 45?")
for d,lab in [(16,'2^k logical dim = 1 generation Weyl content'),(48,'codewords'),
              (256,'full register 2^8'),(45,'valid SM fermions')]:
    ng=F(48,d); print(f"    |O_h|/{d:3d} ({lab:42s}) = {ng}  {'INTEGER' if ng.denominator==1 else ''}")
print("  => 16 is the choice that gives 3. BUT 16 = one generation's Weyl content is the")
print("     physically RIGHT denominator for a 'capacity / content = count' argument, and is")
print("     independently anchored (§15 item 116). So the denominator is defensible -- the")
print("     numerator (|O_h| as total capacity) is the weak leg, not the denominator.")

H("ITEM 133 -- (D) the deepest point: '3 generations' is IDENTIFIED WITH '3 spatial axes'")
print("  Canon's physical interpretation: gen 1,2,3 = X,Y,Z orbits of the Q_3 cell. So the")
print("  '3' of generations is the SAME '3' as spatial dimension -- not an independent 3.")
print("  CONSEQUENCE for the Occam comparison: the framework does NOT derive 3 generations")
print("  ex nihilo. It TIES generation-count to spatial-dimension-count (3D). Given 3D space")
print("  (input, or item-99-derived-as-flat-but-still-3D-by-assumption), 3 generations follows.")
print("  This is a genuine ECONOMY (one '3' doing two jobs vs SO(10)+SU(3)_F's separate family")
print("  symmetry) -- but it is 'gen-count <=> dim-count', contingent on 3D, not '3 from nothing'.")

H("VERDICT -- items 99 & 133")
for ln in [
 "ITEM 99 (4.8.8 uniqueness): GREEN. The flatness criterion that selects 4.8.8 is verified",
 "  (4.8.8 uniquely flat among trivalent {4,f,f}; 4.10.10 hyperbolic) and is INDEPENDENTLY",
 "  motivated (observed flat space), not reverse-engineered. The proof chain is forced given",
 "  the 4+4 bit partition, which is an upstream physical assignment (disclosed, like §2.8's",
 "  Z_f/Z_p). Not circular. The honest residual is the usual one: '4.8.8 vertex figure' is the",
 "  2D local figure; the 3D-bulk uniqueness is a separate open item (canon already flags this).",
 "",
 "ITEM 133 (3 generations): AMBER. The formula has a sensible COUNTING FORM (total symmetry",
 "  capacity / per-generation content). Two legs check out: (i) among Platonic cells only the",
 "  cube gives an integer N_g, and it gives exactly 3 -- non-trivial; (ii) the denominator 16",
 "  is the physically-right, independently-anchored per-generation content. BUT two real soft",
 "  spots: (a) the numerator MUST be the full O_h=48 (with reflections) -- rotations-only O=24",
 "  gives 3/2; canon does not justify fermion-copy-count = |O_h|, and itself flags the rep-",
 "  theory 'faithful covering' step as heuristic. This is load-bearing and un-derived. (b) the",
 "  physical interpretation IDENTIFIES the 3 generations WITH the 3 spatial axes (X,Y,Z), so",
 "  the result is 'generation-count = spatial-dimension-count', contingent on 3D space being",
 "  input -- NOT a parameter-free derivation of 3 from nothing.",
 "",
 "  IMPACT ON THE OCCAM NUANCE: the earlier claim 'framework wins on inter-generation count'",
 "  must be QUALIFIED. The framework does not derive 3 ex nihilo; it ties generation-count to",
 "  spatial-dimension-count. That is a genuine economy (one '3' vs SO(10)+SU(3)_F's separate",
 "  family symmetry) IF 3D is granted -- but the 48-not-24 'faithful covering' gap means even",
 "  that economy is currently heuristic, not proven. Net: the inter-generation-count Occam",
 "  credit is REAL but CONDITIONAL (on 3D input + the open faithful-covering rep-theory), not",
 "  bankable as stated. Tier: Proposition, not Theorem.",
]: print(ln)
print("\nALL ASSERTS PASSED.")
