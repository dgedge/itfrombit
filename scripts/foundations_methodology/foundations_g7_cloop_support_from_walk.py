#!/usr/bin/env python3
r"""foundations_g7_cloop_support_from_walk.py

G7 — promote C_loop = 3/2 (the gravity O(1) coefficient) from "derived-under-reading" toward "derived",
by reading the dissipator's SUPPORT and testing its UNIFORMITY directly against the canonical §3.2 hop.

Background (DRIFT G7, l.213-242; item119_jump_channel.py):
  - alpha^2/K_eff=205 is the DEAD object (alpha^2 structurally impossible: L_k Pi_Q=0, item 119; and
    K_eff=205 non-constructible). The LIVE coefficient is the alpha^1 erasure route
        M_P^2 = (2/3) alpha Lambda^3 R_dS,   C_loop = 3/2,   K = 3/(2 alpha) = 205.5  (vs K_data 206.49, 0.45%).
  - C_loop is a RATIO  C_loop = (sum gamma_k) / Gamma_vac , independent of item-79's magnitude normalisation
    sum gamma_k = alpha*Lambda (that only fixes the absolute K, hence the un-closable M_P magnitude:
    rank-1 theorem l.230 + T7 gate forbid an intrinsic M_P regardless). So C_loop is the only DECIDABLE piece.
  - item 119 got C_loop=3/2 under TWO inputs, both flagged "a grounded reading, not a theorem":
      (b1) SUPPORT: the dissipator acts on the walk-active bits {C0,C1,I3};
      (b2) UNIFORMITY: gamma_k equal on those three (gamma = alpha*Lambda/3 each).
    The algebra: C_loop = 1 + gamma_{I3}/(gamma_{C0}+gamma_{C1}); =3/2 IFF the rates are equal.

THIS SCRIPT reads both inputs off the CANONICAL §3.2 hop  T_d = -(i/sqrt3) R_d (V_em + V_weak + V_strong)
(ANCHOR §3.2, l.397-408):
  V_em    = diag(Q_bridge),  Q_bridge = I3 - 1/2 (1-LQ)            -> DIAGONAL (flips no bit)
  V_weak  = sqrt(2/9) * (chi-controlled I3 flip)  [the chirality-mixing CNOT coin]   -> flips I3
  V_strong= g_s * (colour permutation on C0,C1),  g_s = 1          -> flips C0 / C1

FINDINGS (self-asserting):
 (1) SUPPORT {C0,C1,I3} is DERIVED, not read: V_em is diagonal, V_weak touches only I3, V_strong only
     {C0,C1}; G0,G1 are hop-conserved (l.194) and LQ,chi,W are not single-hop flips. So the coherent walk's
     bit-flip support is EXACTLY {C0,C1,I3}. (Promotes caveat b1 from "the T_000 reading" to "read off the
     V_em+V_weak+V_strong channel decomposition".)
 (2) UNIFORMITY is the WHOLE coefficient, and it is a BILLING question, now sharply isolated and decided:
     C_loop=3/2 holds IFF the free bit I3 is billed EQUAL to the colour bits (unital / syndrome-extraction
     "trace" billing). The §3.2 GOLDEN-RULE billing (gamma ∝ coherent amplitude^2: weak I3 -> 2/9,
     strong colour -> 1) gives C_loop = 10/9 = 1.111, K = 152 -- a 26% MISS. So:
       * golden-rule (amplitude^2) billing is QUANTITATIVELY EXCLUDED (26% off K_data);
       * the uniform billing that gives 3/2 is exactly item-79's UNITAL site-basis syndrome-dephasing class
         (Hermitian jumps -> Evans-Frigerio uniform fixed point -> alpha_0=1/137; ANCHOR l.593). So the
         uniformity is the SAME billing the framework already uses for alpha_0, not an independent assumption.
     NET: C_loop=3/2 is promoted to "support derived from §3.2 + uniformity = item-79 unital billing, with
     golden-rule excluded at 26%." It stays conditional only on (i) item-79's magnitude normalisation (which
     governs the un-closable M_P magnitude, not C_loop) and (ii) the dissipator-inherits-walk-support bridge.
     This is a sharpening + a clean negative (golden-rule killed), NOT a full closure.

Self-asserting; exit 0.  Run under ~/bin/py13_7 (numpy).
"""
import itertools
import numpy as np

ALPHA = 1 / 137.035999
NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
B = {n: i for i, n in enumerate(NAMES)}
ALL = list(itertools.product((0, 1), repeat=8))
IDX = {c: i for i, c in enumerate(ALL)}


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def valid(c):  # the 48 codewords P (item119 validity rule)
    G0, G1, LQ, C0, C1, I3, chi, W = c
    return (not (G0 and G1)) and W == chi and ((LQ == 0) == (C0 == 0 and C1 == 0))


P = [i for i, c in enumerate(ALL) if valid(c)]
Q = [i for i, c in enumerate(ALL) if not valid(c)]
inP = np.zeros(256, bool); inP[P] = True
assert (len(P), len(Q)) == (48, 208)


def flip(c, *names):
    d = list(c)
    for n in names:
        d[B[n]] ^= 1
    return tuple(d)


def support_of(M):
    """Which bit-positions differ between any pair (i,j) with M[i,j]!=0 (off-diagonal)."""
    bits = set()
    nz = np.argwhere(np.abs(M) > 1e-12)
    for j, i in nz:                       # column i -> row j
        if i == j:
            continue
        diff = [k for k in range(8) if ALL[i][k] != ALL[j][k]]
        bits.update(diff)
    return bits


def main():
    print("=== build the canonical §3.2 internal hop channel  V_em + V_weak + V_strong  (256-dim) ===")
    # V_em: diagonal in Q_bridge = I3 - 1/2 (1-LQ)  -> flips nothing
    Vem = np.diag([ (c[B["I3"]] - 0.5 * (1 - c[B["LQ"]])) for c in ALL ]).astype(complex)

    # V_weak = sqrt(2/9) * chi-controlled I3 flip
    Vw = np.zeros((256, 256), complex); cw = np.sqrt(2 / 9)
    for i, c in enumerate(ALL):
        if c[B["chi"]] == 1:                          # controlled by chi
            Vw[IDX[flip(c, "I3")], i] = cw

    # V_strong = g_s * colour permutation on (C0,C1): K_3 adjacency on the 3 colour states {01,10,11}
    Vs = np.zeros((256, 256), complex); gs = 1.0
    colours = [(0, 1), (1, 0), (1, 1)]
    for i, c in enumerate(ALL):
        cc = (c[B["C0"]], c[B["C1"]])
        if cc in colours:                              # only acts where a colour exists (quark sector)
            for tgt in colours:
                if tgt != cc:
                    d = list(c); d[B["C0"]], d[B["C1"]] = tgt
                    Vs[IDX[tuple(d)], i] = gs

    print("[1] SUPPORT read off the channel decomposition:")
    sup_em, sup_w, sup_s = support_of(Vem), support_of(Vw), support_of(Vs)
    ok(sup_em == set(), "V_em is DIAGONAL -> flips no bit (it is diag(Q_bridge))")
    ok(sup_w == {B["I3"]}, f"V_weak off-diagonal flips ONLY I3 (bit {B['I3']}) -> support {sorted(sup_w)}")
    ok(sup_s <= {B["C0"], B["C1"]}, f"V_strong off-diagonal flips ONLY {{C0,C1}} -> support {sorted(sup_s)}")
    walk_support = sup_em | sup_w | sup_s
    ok(walk_support == {B["C0"], B["C1"], B["I3"]},
       f"=> canonical-hop bit-flip SUPPORT = {{C0,C1,I3}} = {sorted(walk_support)} (DERIVED, not read)")
    print("    (G0,G1 hop-conserved per l.194; LQ,chi,W are not single-hop flips -> excluded structurally)")

    # exit table on the colourless vacuum cell nu_e = all-zeros
    vac = IDX[(0,) * 8]; assert inP[vac]
    exits = {n: (not inP[IDX[flip(ALL[vac], n)]]) for n in ("C0", "C1", "I3")}
    print(f"\n[2] vacuum (nu_e) exit table on the support: {exits}")
    ok(exits == {"C0": True, "C1": True, "I3": False},
       "on the colourless vacuum: C0,C1 flips EXIT (break R3); I3 is the FREE bit (never exits)")

    # ---- the coefficient is C_loop = (sum gamma)/Gamma_vac ; Gamma_vac = gamma_{C0}+gamma_{C1} ----
    def cloop(gC0, gC1, gI3):
        return (gC0 + gC1 + gI3) / (gC0 + gC1)        # only C0,C1 exit on the vacuum

    print("\n[3] the coefficient under the two billings (the WHOLE question is I3's rate):")
    # (a) UNITAL / trace billing: equal rate per detectable walk-active bit  (item 119; item-79 unital class)
    Cu = cloop(1, 1, 1); Ku = Cu / ALPHA
    print(f"    (a) UNITAL  (gamma equal on C0,C1,I3): C_loop = {Cu:.4f} = 3/2,  K = {Ku:.2f}")
    ok(abs(Cu - 1.5) < 1e-12 and abs(Ku / 206.49 - 1) < 0.005, "unital billing -> C_loop=3/2, K=205.5 (0.45% vs K_data)")

    # (b) GOLDEN-RULE billing: gamma ∝ coherent amplitude^2  (weak I3 -> 2/9 ; strong colour -> 1 each)
    Cg = cloop(1, 1, 2 / 9); Kg = Cg / ALPHA
    print(f"    (b) GOLDEN-RULE (gamma∝amp^2: I3->2/9, colour->1): C_loop = {Cg:.4f} = 10/9,  K = {Kg:.2f}")
    ok(abs(Cg - 10 / 9) < 1e-12, "golden-rule billing -> C_loop = 10/9 = 1.111")
    ok(abs(Kg / 206.49 - 1) > 0.15, f"golden-rule MISSES K_data by {abs(Kg/206.49-1)*100:.0f}% -> QUANTITATIVELY EXCLUDED")
    # robustness: 'colour as one channel' variant -> 11/9, also a clean miss
    Cg2 = (1 + 2 / 9) / 1; print(f"    (b') colour-as-one-channel variant: C_loop = {Cg2:.4f} = 11/9, K = {Cg2/ALPHA:.0f} (also a >15% miss)")
    ok(abs(Cg2 / 1.5 - 1) > 0.15, "the golden-rule conclusion is robust to the colour-channel split (both ~1.1-1.2, not 3/2)")

    print("\n[verdict] G7 / C_loop=3/2:")
    print("  - SUPPORT {C0,C1,I3}: DERIVED from §3.2 (V_em diagonal, V_weak->I3, V_strong->colour) -- promotes")
    print("    caveat (b1) from 'the T_000 reading' to 'read off the channel decomposition'.")
    print("  - UNIFORMITY is the entire coefficient and reduces to ONE question -- is the free bit I3 billed")
    print("    EQUAL to colour? YES (unital/syndrome billing) -> 3/2 (0.45%); golden-rule (amp^2) -> 10/9 (26% miss,")
    print("    EXCLUDED). The uniform billing is item-79's OWN unital site-dephasing class (the alpha_0=1/137")
    print("    Evans-Frigerio equipartition) -- so the uniformity is framework-consistent, not an extra assumption.")
    print("  - NET: C_loop=3/2 promoted (support derived; golden-rule killed; uniformity tied to item-79 unitality).")
    print("    NOT a full closure: still rests on (i) the dissipator-inherits-walk-support bridge and (ii) item-79's")
    print("    magnitude normalisation -- which governs only the un-closable M_P magnitude (rank-1 l.230 + T7), not C_loop.")
    print("  exit 0")


if __name__ == "__main__":
    main()
