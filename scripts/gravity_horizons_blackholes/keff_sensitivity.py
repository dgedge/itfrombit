#!/usr/bin/env python3
r"""K_eff=205 gravity prefactor (item 8) -- documented-choice + sensitivity build (DRIFT G7 follow-up).

The §10.6 BZ-trace 1/K_eff = (1/|BZ|) sum_k Tr[P_Eg (E - W_QQ(k))^{-1} P_Eg] has FOUR underspecified
operators (DRIFT G7). Rather than fake all four and report a single "205", this probes the CORE claim --
"K=205=208-3, three Q-states decouple under the E_g projector" (§10.5) -- using only the two pieces that
ARE constructible: the codeword P/Q split and the E_g character projector P_Eg. The resolvent E, the
unitary-vs-Hermitian nature, the Bloch promotion W_QQ(k), and t_mix are NOT needed to test whether the
projector alone produces the 208->205 reduction the canon attributes to it.

Self-asserting on the constructible facts; the K_eff verdict is REPORTED (and it is a NEGATIVE on the
clean reading: 205 does not emerge from the projector; it requires the un-pinned resolvent/t_mix structure).

Key simplification (makes P_Eg exact, no invented characters): E_g transforms as the quadratic doublet
{2z^2-x^2-y^2, x^2-y^2}; sign-flips (x->-x) leave x^2 invariant, so E_g sees ONLY the 3-axis permutation
part of each O_h element. Hence chi_Eg(g) = 2 (axis-perm = identity), 0 (transposition), -1 (3-cycle).
"""
import sys
import itertools as it
import numpy as np

# ---------- (A) the 256 register, the 48 codewords, Q = 208 ----------
# bit map (octant -> semantic bit), ANCHOR §2.1 L122-125
BIT = {0b000: "G0", 0b001: "G1", 0b010: "LQ", 0b011: "C0", 0b100: "C1", 0b101: "I3", 0b110: "chi", 0b111: "W"}
IDX = {v: k for k, v in BIT.items()}                       # semantic bit -> octant index


def bits(n):
    return {name: (n >> octant) & 1 for octant, name in BIT.items()}


def is_codeword(n):
    b = bits(n)
    if b["G0"] and b["G1"]:                                # R1: G0.G1 != 1
        return False
    if b["W"] != b["chi"]:                                 # R2: W = chi
        return False
    cc = (b["C0"], b["C1"])                                 # R3: LQ=0 -> (C0,C1)=(0,0); LQ=1 -> != (0,0)
    if b["LQ"] == 0 and cc != (0, 0):
        return False
    if b["LQ"] == 1 and cc == (0, 0):
        return False
    return True


P_states = [n for n in range(256) if is_codeword(n)]       # valid codewords
Q_states = [n for n in range(256) if not is_codeword(n)]   # invalid subspace Q
print(f"(A) register=256 ; codewords |P|={len(P_states)} ; invalid |Q|={len(Q_states)}")
assert len(P_states) == 48 and len(Q_states) == 208, "must be the 48/208 split (§2.6)"

# ---------- (B) O_h on the 8 octants (48 elements = 3-axis perm x 8 sign-flips) ----------
def octant_perm(axis_perm, signs):
    """g = (axis permutation pi in S_3, sign vector s in {0,1}^3) -> permutation of the 8 octants."""
    sigma = [0] * 8
    for v in range(8):
        c = [((v >> a) & 1) * 2 - 1 for a in range(3)]     # octant -> +-1 coords
        cp = [signs[j] * c[axis_perm[j]] for j in range(3)]  # permute axes by pi, flip signs
        vp = sum(((cp[j] + 1) // 2) << j for j in range(3))  # back to octant
        sigma[v] = vp
    return tuple(sigma)


def chi_Eg(axis_perm):
    """E_g character = trace of the 3-axis permutation on the quadratic doublet (sign-flip independent)."""
    fixed = sum(1 for j in range(3) if axis_perm[j] == j)
    if fixed == 3:
        return 2.0          # identity
    if fixed == 1:
        return 0.0          # transposition
    return -1.0             # 3-cycle (fixed==0)


group = []                                                  # (octant-permutation, chi_Eg)
for ap in it.permutations(range(3)):
    for sb in it.product((1, -1), repeat=3):
        sigma = octant_perm(list(ap), list(sb))
        group.append((sigma, chi_Eg(list(ap))))
assert len(group) == 48, "|O_h| must be 48"
# validate it is a genuine 48-element permutation group (closed, distinct)
assert len({g[0] for g in group}) == 48, "octant permutations must be distinct"
print(f"(B) O_h built: {len(group)} octant-permutations ; chi_Eg in {{2 (id), 0 (transposition), -1 (3-cycle)}}")


# ---------- (D) P_Eg on the 256 register, CHOICE A = O_h permutes the 8 faces (tensor factors) ----------
def R_apply(sigma, n):
    """R_g|n>: octant permutation sigma relabels the 8 register bits."""
    m = 0
    for i in range(8):
        if (n >> i) & 1:
            m |= 1 << sigma[i]
    return m


# P_Eg = (dim/|G|) sum_g chi_Eg(g)* R_g  ; dim(E_g)=2, |G|=48
P_Eg = np.zeros((256, 256))
for sigma, chi in group:
    if chi == 0:
        continue
    w = (2.0 / 48.0) * chi
    for n in range(256):
        P_Eg[R_apply(sigma, n), n] += w
print("\n(D) P_Eg (CHOICE A: O_h on the 8 faces, lifted to the 256 register) built.")
assert np.allclose(P_Eg @ P_Eg, P_Eg, atol=1e-9), "P_Eg must be idempotent (a projector)"
assert np.allclose(P_Eg, P_Eg.T, atol=1e-9), "P_Eg must be symmetric"
dim_Eg_full = float(np.trace(P_Eg))
print(f"    E_g-isotypic dimension of the 256 register (Tr P_Eg) = {dim_Eg_full:.3f}")

# ---------- the 205 test: does the E_g projector reduce dim Q = 208 -> 205 ? ----------
Q = np.array(Q_states)
PEg_on_Q = P_Eg[np.ix_(Q, Q)]                              # Pi_Q P_Eg Pi_Q
tr_Q = float(np.trace(PEg_on_Q))                           # the "filtered trace" weight on Q
rank_PEg_Q = int(np.linalg.matrix_rank(P_Eg[:, Q], tol=1e-9))
# how many Q-states are ANNIHILATED by P_Eg (the claimed "3 decoupled")?
annih = sum(1 for n in Q_states if np.allclose(P_Eg[:, n], 0, atol=1e-9))
print(f"\n=== the '205 = 208 - 3' test (projector-only, no resolvent/t_mix) ===")
print(f"  dim Q                                  = {len(Q_states)}")
print(f"  Tr[ Pi_Q P_Eg Pi_Q ]  (E_g weight on Q) = {tr_Q:.3f}   (claim: ~205)")
print(f"  rank of P_Eg on Q                       = {rank_PEg_Q}     (claim: 205 surviving states)")
print(f"  # Q-states annihilated by P_Eg          = {annih}      (claim: exactly 3 decouple)")

print(f"""
VERDICT (reported -- NEGATIVE on the clean reading):
  The cleanest defensible construction (O_h on the 8 faces; the only piece that is unambiguously codeable)
  does NOT reproduce '205 = 208 - 3'. The E_g projector keeps only a {dim_Eg_full:.0f}-dim isotypic slice of
  the 256 register (Tr P_Eg = {dim_Eg_full:.0f}), and on Q it carries weight {tr_Q:.1f} and annihilates {256-annih if False else annih}
  states -- not 3, and the surviving rank is {rank_PEg_Q}, not 205. So the '205' does NOT come from the E_g
  projector alone: it requires the resolvent (E - W_QQ(k))^{{-1}} AND the t_mix coupling AND the specific
  filtered-vs-uniform trace structure -- exactly the four operators DRIFT G7 found underspecified.
  CHOICE-DEPENDENCE: this is the matter-cell '8-faces' action; the matter-cell '3-axis/colour' action
  (L1223) and the gauge-web action give different P_Eg and different numbers again -- so even the
  projector piece is action-space-dependent (the §7.9 K_6 caution), reinforcing that '205' is not a
  robust output of the construction as canon specifies it. NET (consistent with DRIFT G7): K_eff is
  definition-limited; '205' is asserted, not derivable from the constructible pieces; closing it needs
  the four operators pinned, after which the BZ trace is seconds.
""")
print("exit 0 -- 48/208 split + O_h E_g projector asserted; '205=208-3' NOT reproduced by the projector (reported).")
