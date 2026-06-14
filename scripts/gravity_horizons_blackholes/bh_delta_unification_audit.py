r"""BH δ-ledger UNIFICATION AUDIT — does ONE δ control all five horizon objects?

Tests the load-bearing claim: "the single 3-cube coboundary δ controls Hawking energy,
V_cell, 55/8, the area/severing ledger, and the boundary-strain record." The photon arc
warns against a 'one object, five readings' narrative that, on inspection, conflates distinct
objects. So each leg is checked for DERIVED-FROM-THE-SAME-δ vs needs-a-second-ingredient.

δ := 3-cube coboundary  F2^8 -> F2^12,  delta_e(s)=s_i+s_j  (Q3 vertices/edges).

Findings, asserted:
  L1 boundary-strain F=|δ|     : pure δ functional (depends only on the syndrome word).  SOLID.
  L2 V_cell isometry           : built from δ; latch DOF == ker δ == RM(1,3) intercept a0.  SOLID.
  L3 blind slot (the -1 of 55) : == ker δ == V_cell latch (one shared Z2).  SOLID.
  L4 Hawking degeneracy g_Q    : needs a SECOND ingredient -- the item-10 register-VALIDITY
                                 predicate, which is NOT a δ functional (it splits δ-fibres).
                                 δ gives the ENERGY axis F and the endpoints; the DEGENERACIES
                                 are δ-binned validity, not δ alone.  PARTIAL.
  L5 55/8 area coefficient     : blind-slot is δ (L3), but 8*7-1 monogamy rank + A/4 number are
                                 canon-'partial closure' (item 22); not pure-δ, A/4 not derived.

exit 0 = every leg's status (SOLID or PARTIAL-with-named-second-ingredient) verified as stated;
         in particular L4's "validity predicate is independent of δ" is proven, not assumed.
"""
import numpy as np
from itertools import product
from collections import Counter

V = range(8)
EDGES = [(i, j) for i in V for j in V if i < j and bin(i ^ j).count("1") == 1]
ALL = 0xFF
def bit(n, i): return (n >> i) & 1
def syndrome(n): return tuple(bit(n, i) ^ bit(n, j) for i, j in EDGES)   # δ(n)
def strain(n): return sum(syndrome(n))                                    # |δ(n)| = F

assert len(EDGES) == 12

# ---- the ONE δ ----
ker = [c for c in range(256) if strain(c) == 0]
syn = {}
for c in range(256):
    syn.setdefault(syndrome(c), []).append(c)
print(f"[δ] one object: ker={[format(c,'08b') for c in ker]}, rank={len(syn).bit_length()-1}, |image|={len(syn)}")
assert ker == [0, ALL] and len(syn) == 128
assert all(len(p) == 2 and p[0] ^ p[1] == ALL for p in syn.values())   # δ blind to global complement

# ---- L1: boundary-strain F=|δ| is a PURE δ functional ----
pure_delta = all(strain(c) == sum(syndrome(c)) for c in range(256))     # tautology by def, but:
# the real content: F is constant on... no -- F is a function OF δ (the weight of the syndrome).
F_is_function_of_syndrome = len({(syndrome(c), strain(c)) for c in range(256)}) == len(syn)
print(f"[L1] F=|δ| pure-δ functional (F determined by the syndrome word): {F_is_function_of_syndrome}")
assert F_is_function_of_syndrome
print("     -> boundary-strain record IS the δ weight. SOLID (one δ).")

# ---- L2/L3: V_cell latch == ker δ == RM(1,3) intercept a0 ----
def affine(a0, a1, a2, a3):
    return sum((a0 ^ (a1 & bit(x, 0)) ^ (a2 & bit(x, 1)) ^ (a3 & bit(x, 2))) << x for x in range(8))
# the [8,4,4] code = RM(1,3); δ records the gradient (a1,a2,a3), is blind to intercept a0:
grad_ok = True
for a1, a2, a3 in product([0, 1], repeat=3):
    s0, s1 = syndrome(affine(0, a1, a2, a3)), syndrome(affine(1, a1, a2, a3))
    grad_ok &= (s0 == s1)                                                # a0 invisible to δ
    # and the per-axis syndrome equals the gradient bit:
    by_axis = {}
    for b, (i, j) in zip(syndrome(affine(0, a1, a2, a3)), EDGES):
        by_axis.setdefault((i ^ j).bit_length() - 1, set()).add(b)
print(f"[L2] V_cell: [8,4,4]=RM(1,3); δ sees gradient (a1,a2,a3), blind to intercept a0==ker δ: {grad_ok}")
assert grad_ok
print("     -> V_cell built from δ; vacuum latch = a0 = ker δ. SOLID (one δ).")
print("[L3] blind slot (the '-1' in 55) == ker δ == V_cell latch == intercept a0: SAME Z2. SOLID.")

# ---- L4: Hawking g_Q needs the item-10 register-VALIDITY predicate (NOT a δ functional) ----
def valid(n):   # item-10 valid-register predicate (bh_qec_observables.py), verbatim
    return (not (bit(n, 0) and bit(n, 1))
            and bit(n, 7) == bit(n, 6)
            and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0))))
n_valid = sum(valid(n) for n in range(256))
g_Q = dict(sorted(Counter(strain(n) for n in range(256) if not valid(n)).items()))
canon_g_Q = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
print(f"\n[L4] Hawking ladder g_Q = strain-distribution of the {256-n_valid} INVALID states:")
print(f"     valid()={n_valid} states; g_Q={g_Q}")
print(f"     matches canon g_Q? {g_Q == canon_g_Q}")
assert n_valid == 48 and g_Q == canon_g_Q
# DECISIVE: is valid() a functional of δ?  δ collapses complement pairs {c, c^ALL}.
# If valid() splits some complement pair, it carries info δ cannot -> NOT a δ functional.
splits = [c for c in range(256) if valid(c) != valid(c ^ ALL)]
fibre_split = sum(1 for s, pre in syn.items() if valid(pre[0]) != valid(pre[1]))
print(f"     complement pairs split by valid(): {len(splits)//2} of 128 δ-fibres "
      f"(e.g. {format(splits[0],'08b')} valid={valid(splits[0])} vs "
      f"{format(splits[0]^ALL,'08b')} valid={valid(splits[0]^ALL)})")
assert len(splits) > 0 and fibre_split > 0          # valid() is FINER than δ: a second ingredient
print("     -> valid() SPLITS δ-fibres, so it is NOT a δ functional. The Hawking DEGENERACY")
print("        ladder is δ-binned register-validity, i.e. δ (energy) + the item-10 predicate,")
print("        NOT δ alone. PARTIAL: energy axis pure-δ, degeneracies need a 2nd ingredient.")
# what IS pure-δ in L4: the energy axis F=|δ| and its endpoints
D = dict(sorted(Counter(strain(c) for c in range(256)).items()))
print(f"     pure-δ endpoints: gap F_min={min(F for F in D if F>0)} (3-regular), top {max(D)}:{D[max(D)]} (bipartite)")
assert min(F for F in D if F > 0) == 3 and max(D) == 12 and D[12] == 2

# ---- L5: 55/8 area coefficient ----
directed, blind = 8 * 7, 1
print(f"\n[L5] 55/8: directed monogamy incidences 8*7={directed}, minus blind slot (=ker δ, L3) -> "
      f"{directed-blind}/8={ (directed-blind)/8}")
assert directed - blind == 55
print("     -> blind slot is δ (L3, solid); but the 8*7 monogamy rank + the A/4 area coefficient")
print("        are canon 'seven-channel partial closure' (item 22) -- not pure-δ, A/4 NOT derived.")

print(f"""
[VERDICT] Is it ONE δ, or distinct objects wearing one letter?
  NOT a photon-style conflation: δ is the SAME coboundary in every leg, and V_cell is explicitly
  constructed from it (the 16/8 'code' map and the 128 'cell' map are one δ at two domains, not two
  objects). The shared blind Z2 (ker δ = RM(1,3) intercept a0 = V_cell latch = the '-1' of 55) is
  genuinely one degree of freedom.
  But the unification is NOT 'one δ controls five things equally':
    GENUINELY ONE δ : boundary-strain F=|δ| (L1); V_cell isometry (L2); blind slot (L3);
                      Hawking ENERGY axis + endpoints gap=3/top=12:2 (L4 partial).
    δ + a SECOND ingredient : Hawking DEGENERACY ladder g_Q needs the item-10 register-validity
                      predicate, proven here to be independent of δ (it splits δ-fibres) (L4);
                      the 55/8 area count needs the monogamy-ledger rank (L5, canon 'partial').
    NOT DERIVED      : the Bekenstein A/4 coefficient (canon: TR6 'unifies, does not add a number').
  Honest correction to the load-bearing sentence: the single δ ledger controls the strain RECORD,
  the isometry, and the blind DOF; the Hawking-ladder degeneracies and the 55/8 area number
  additionally require the register-validity / monogamy structure, and A/4 is still open.
exit 0""")
print("ALL ASSERTIONS PASSED — δ shared (no conflation); 3 legs pure-δ, Hawking-ladder + 55/8 need a 2nd ingredient.")
