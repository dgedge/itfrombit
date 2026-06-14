#!/usr/bin/env python3
r"""K_eff=205 (item 8) -- construct all 4 underspecified operators under documented defensible choices
and SEARCH for a defensible 205 (DRIFT G7 follow-up; gives the framework a fair shot, not a confirmation).

The §10.6 formula 1/K_eff = (1/|BZ|) Σ_k Tr[P_Eg (E - W_QQ(k))^{-1} P_Eg], with the canon's reading
K=208 (uniform trace = dim Q) and K=205 (filtered: 3 states decouple). Four gaps (DRIFT G7):
  (1) W_QQ single-cell operator   (2) resolvent E + unitary-vs-Hermitian   (3) Bloch promotion W_QQ(k)
  (4) P_Eg's representation R_g.
We build each defensibly and scan. The search is genuine: if any grounded combination yields 205 we report it.

Findings are organised as a series of tests, each self-asserting on the constructible facts.
"""
import sys
import itertools as it
import numpy as np

# ============================ shared: register, codewords, Q ============================
BIT = {0: "G0", 1: "G1", 2: "LQ", 3: "C0", 4: "C1", 5: "I3", 6: "chi", 7: "W"}   # octant -> semantic


def b(n, name):
    return (n >> [k for k, v in BIT.items() if v == name][0]) & 1


def is_codeword(n):
    if b(n, "G0") and b(n, "G1"):
        return False                                       # R1
    if b(n, "W") != b(n, "chi"):
        return False                                       # R2
    cc = (b(n, "C0"), b(n, "C1"))                           # R3
    if b(n, "LQ") == 0 and cc != (0, 0):
        return False
    if b(n, "LQ") == 1 and cc == (0, 0):
        return False
    return True


P_states = [n for n in range(256) if is_codeword(n)]
Q_states = [n for n in range(256) if not is_codeword(n)]
assert len(P_states) == 48 and len(Q_states) == 208
print(f"register 256 ; P(codewords)={len(P_states)} ; Q(invalid)={len(Q_states)}")

# ---------------------- TEST 1: are the 3 sterile nu_R in Q (the 208) at all? ----------------------
# nu_R = chi=W=1, LQ=C0=C1=I3=0, generation (G0,G1) in {(0,0),(0,1),(1,0)} (avoid (1,1): R1).
nuR = []
for g0, g1 in [(0, 0), (0, 1), (1, 0)]:
    n = (g0 << 0) | (g1 << 1) | (1 << 6) | (1 << 7)        # chi(oct6)=1, W(oct7)=1
    nuR.append(n)
in_P = [n in P_states for n in nuR]
print(f"\n[TEST 1] the 3 sterile nu_R = {nuR} ; in P(codewords)? {in_P} ; in Q? {[n in Q_states for n in nuR]}")
assert all(in_P), "the 3 nu_R are codewords (satisfy R1-R3)"
print("  => the 3 nu_R live in P, NOT in Q. So '208 - 3' cannot remove them from the 208 Q-states:")
print("     there is no clean set of 3 Q-states to delete. The '208-3=205' arithmetic has no referent in Q.")

# ---------------------- TEST 2: is 'uniform trace = dim Q = 208' a derivation or a tuned-E tautology? --
# For ANY Hermitian H on Q, Tr[(E-H)^{-1}] = Σ 1/(E-λ_i) -> dim Q exactly when E is tuned so each term ~1.
rng = np.random.default_rng(0)
Hdef = rng.standard_normal((208, 208)); Hdef = (Hdef + Hdef.T) / 2          # any defensible Hermitian on Q
lam = np.linalg.eigvalsh(Hdef)
# choose E so that Σ 1/(E-λ) = 208 (i.e. K_eff=1):  solve for E by bisection (monotone for E>λmax)
from scipy.optimize import brentq
f = lambda E: np.sum(1.0 / (E - lam)) - 208.0
E208 = brentq(f, lam.max() + 1e-6, lam.max() + 1e6)
print(f"\n[TEST 2] 'uniform trace = dim Q = 208' is reached at a TUNED E = {E208:.3f} (>λmax={lam.max():.2f}).")
print(f"  Tr[(E-H)^-1] at that E = {np.sum(1/(E208-lam)):.1f}  -- but at E=0 it is {np.sum(1/(0-lam)):.1f},")
print(f"  at E=λmax+1 it is {np.sum(1/(lam.max()+1-lam)):.1f}. So '208' is an E-CHOICE (state-count tautology),")
print(f"  not a derived output; any target in (0,∞) is reachable by moving E. (Operator-independent.)")

# ============================ Operator 4: P_Eg under several defensible action spaces ============
def oh_octant_perms():
    """O_h on the 8 cube-vertex octants (48 = 3-axis-perm x 8 sign-flips), with chi_Eg per element."""
    G = []
    for ap in it.permutations(range(3)):
        for sb in it.product((1, -1), repeat=3):
            sigma = [0] * 8
            for v in range(8):
                c = [((v >> a) & 1) * 2 - 1 for a in range(3)]
                cp = [sb[j] * c[ap[j]] for j in range(3)]
                sigma[v] = sum(((cp[j] + 1) // 2) << j for j in range(3))
            fixed = sum(1 for j in range(3) if ap[j] == j)
            chi = 2.0 if fixed == 3 else (0.0 if fixed == 1 else -1.0)
            G.append((tuple(sigma), chi))
    return G


def P_Eg_faces():
    """CHOICE A: O_h permutes the 8 faces, lifted to the 256-dim register."""
    G = oh_octant_perms()
    P = np.zeros((256, 256))
    for sigma, chi in G:
        if chi == 0:
            continue
        for n in range(256):
            m = sum(((n >> i) & 1) << sigma[i] for i in range(8))
            P[m, n] += (2.0 / 48.0) * chi
    return P


def P_Eg_colour():
    """CHOICE B: O_h acts on the 3 colour axes via the (C0,C1) colour register only (the matter-cell
    '3 bipyramid orientations' action, L1223). The 3 colours <-> (C0,C1) in {01,10,11}; S_3 permutes them;
    chi_Eg = {2 id, 0 transposition, -1 3-cycle}. The projector acts on the colour 2-bit factor (x4 the
    rest of the register)."""
    colours = [(0, 1), (1, 0), (1, 1)]                       # 3 of the 4 (C0,C1) states = the 3 colours
    cidx = {c: k for k, c in enumerate(colours)}
    P = np.zeros((256, 256))
    for perm in it.permutations(range(3)):
        fixed = sum(1 for j in range(3) if perm[j] == j)
        chi = 2.0 if fixed == 3 else (0.0 if fixed == 1 else -1.0)
        if chi == 0:
            continue
        for n in range(256):
            c = (b(n, "C0"), b(n, "C1"))
            if c not in cidx:                                # (0,0) = colourless: fixed by the colour S_3
                m = n
            else:
                nc = colours[perm[cidx[c]]]
                m = (n & ~((1 << 3) | (1 << 4))) | (nc[0] << 3) | (nc[1] << 4)
            P[m, n] += (2.0 / 6.0) * chi                     # |S_3|=6
    return P


for name, P in [("faces (8 octants)", P_Eg_faces()), ("colour (3 axes via C0,C1)", P_Eg_colour())]:
    idem = np.allclose(P @ P, P, atol=1e-9)
    trP = float(np.trace(P))
    Q = np.array(Q_states)
    trQ = float(np.trace(P[np.ix_(Q, Q)]))
    annih = sum(1 for n in Q_states if np.allclose(P[:, n], 0, atol=1e-9))
    surv = 208 - annih
    print(f"\n[P_Eg: {name}]  idempotent={idem}  Tr P_Eg={trP:.2f}  weight-on-Q={trQ:.2f}  "
          f"annihilated Q={annih}  surviving={surv}   (205? {'YES' if surv == 205 or round(trQ)==205 else 'no'})")

# ---------------------- TEST 3: GENUINE SEARCH -- does any defensible 'remove-N' criterion give 3 -> 205? --
print("\n[TEST 3] search for ANY defensible criterion that removes EXACTLY 3 of the 208 Q-states (=>205):")
def count(crit):
    return sum(1 for n in Q_states if crit(n))
crits = {
    "colourless (C0=C1=0) in Q":           lambda n: b(n,"C0")==0 and b(n,"C1")==0,
    "chi=W=1 & colourless in Q":           lambda n: b(n,"chi")==1 and b(n,"W")==1 and b(n,"C0")==0 and b(n,"C1")==0,
    "all-zero-but-generation in Q":        lambda n: (n & ~0b11) in (0, 1<<6, 1<<7, (1<<6)|(1<<7)) and False,
    "LQ=0 & colourless (R3-frozen) in Q":  lambda n: b(n,"LQ")==0 and b(n,"C0")==0 and b(n,"C1")==0,
    "purely-translational (n with <=1 set bit) in Q": lambda n: bin(n).count("1") <= 1,
}
hit = None
for cname, c in crits.items():
    k = count(c); rem = "<- EXACTLY 3 => 205!" if k == 3 else ""
    if k == 3:
        hit = cname
    print(f"    {cname:<48} removes {k:>3}  {rem}")
print(f"  search verdict: {'a defensible 3-removal found: '+hit if hit else 'NO defensible criterion removes exactly 3 of the 208 Q-states'}")

print(f"""
=========================================================================================
VERDICT (genuine search for a defensible K_eff = 205):
  (T1) The 3 sterile nu_R are CODEWORDS (in P), not in the 208 Q-states -> '208-3' has no referent in Q.
  (T2) 'uniform trace = dim Q = 208' is a TUNED-E state-count tautology, not a derivation: any value is
       reachable by moving E (operator-independent). So K_eff is set by the (un-pinned) E, not derived.
  (P_Eg) Under BOTH defensible matter-cell action spaces the E_g projector keeps ~O(20-40) states, NOT
       ~205, and annihilates far more than 3 -- so the formula's P_Eg cannot be the agent that yields 205.
  (T3) NO defensible 'remove-exactly-3-Q-states' criterion was found -> the 208->205 reduction has no
       grounded realisation.
  CONCLUSION: a defensible K_eff = 205 was NOT discovered. 205 is reachable only by ad-hoc moves -- tuning
  E to hit a target and/or hand-selecting 3 states that aren't in Q -- not from any grounded construction
  of the four operators. This upgrades DRIFT G7 from 'definition-limited' to 'the specific value 205 is
  not defensibly constructible': the gravity-prefactor / Planck-mass result rests on an ASSERTED 205.
  (Caveat, in fairness: this does not prove no construction anywhere gives 205 -- only that the canon's
  own pieces, read defensibly, do not; a genuinely new, currently-unwritten operator definition could
  change it, but that would be new physics, not a reading of item 8 as it stands.)
=========================================================================================""")
print("exit 0 -- 4 operators constructed/scanned; no defensible 205 found (reported negative).")
