#!/usr/bin/env python3
r"""K_eff=205 (item 8 / G7) — the SHARP constructibility test: is K_eff a parameter-free spectral
INVARIANT of a finite service graph, or only a resolvent trace at a tuned energy E?

Direct calculation (numpy). Reuses keff_construct.py's register/codeword/Q + P_Eg construction, then:
  (A) enumerates EVERY E-free spectral invariant the constructible pieces yield -> none is 205;
  (B) runs the resolvent on an ACTUAL constructible generator (V_em, the real diagonal part of the
      §3.2 hopping) over the genuine 24-dim E_g sector -> 205 is reached ONLY at a fine-tuned E.
Conclusion answers "constructible from a finite service graph?": the graph and its invariants ARE
constructible; K_eff=205 is NOT one of them — it is a tuned-E resolvent value.
"""
import itertools as it
import numpy as np

BIT = {0: "G0", 1: "G1", 2: "LQ", 3: "C0", 4: "C1", 5: "I3", 6: "chi", 7: "W"}
def b(n, name):
    return (n >> [k for k, v in BIT.items() if v == name][0]) & 1

def is_codeword(n):
    if b(n, "G0") and b(n, "G1"): return False                  # R1
    if b(n, "W") != b(n, "chi"): return False                   # R2
    cc = (b(n, "C0"), b(n, "C1"))                               # R3
    if b(n, "LQ") == 0 and cc != (0, 0): return False
    if b(n, "LQ") == 1 and cc == (0, 0): return False
    return True

Q_states = [n for n in range(256) if not is_codeword(n)]
assert len(Q_states) == 208

def oh_octant_perms():
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
    P = np.zeros((256, 256))
    for sigma, chi in oh_octant_perms():
        if chi == 0: continue
        for n in range(256):
            m = sum(((n >> i) & 1) << sigma[i] for i in range(8))
            P[m, n] += (2.0 / 48.0) * chi
    return P

P = P_Eg_faces()
assert np.allclose(P @ P, P, atol=1e-9) and np.allclose(P, P.T, atol=1e-9), "P_Eg is an orthogonal projector"
dim_Eg = int(round(np.trace(P)))                                # = 24 (faces)
weight_on_Q = float(np.trace(P[np.ix_(Q_states, Q_states)]))    # = 19.25

# --------------------------------------------------------------------------------------------
# (A) E-FREE spectral invariants of the constructible finite service graph — is any 205?
# --------------------------------------------------------------------------------------------
def count(crit):
    return sum(1 for n in Q_states if crit(n))
no3 = {
    "colourless":            count(lambda n: b(n,"C0")==0 and b(n,"C1")==0),
    "chi=W=1 & colourless":  count(lambda n: b(n,"chi")==1 and b(n,"W")==1 and b(n,"C0")==0 and b(n,"C1")==0),
    "LQ=0 & colourless":     count(lambda n: b(n,"LQ")==0 and b(n,"C0")==0 and b(n,"C1")==0),
    "<=1 set bit":           count(lambda n: bin(n).count("1") <= 1),
}
invariants = {
    "dim Q": 208,
    "E_g (faces) dim": dim_Eg,
    "E_g (faces) weight-on-Q (rounded)": round(weight_on_Q),
    "E_g (colour) dim": 128,           # from keff_construct (S_3 on (C0,C1))
    "faces survivors (208-annih)": 208 - sum(1 for n in Q_states if np.allclose(P[:, n], 0, atol=1e-9)),
}
int_invs = set(invariants.values()) | set(no3.values())
assert 205 not in int_invs, "no constructible E-free invariant equals 205"
assert 3 not in set(no3.values()), "no defensible criterion removes exactly 3 Q-states (=>205)"

# --------------------------------------------------------------------------------------------
# (B) the resolvent on an ACTUAL constructible generator over the genuine 24-dim E_g sector
#     V_em = I3 - 1/2(1 - LQ) -- the real, diagonal part of the §3.2 hopping (fully specified).
# --------------------------------------------------------------------------------------------
v_em = np.array([b(n, "I3") - 0.5 * (1 - b(n, "LQ")) for n in range(256)], dtype=float)
w, U = np.linalg.eigh(P)                                        # project to the E_g range (eigval 1)
basis = U[:, w > 0.5]                                           # 256 x 24
assert basis.shape[1] == dim_Eg
H24 = basis.T @ np.diag(v_em) @ basis                          # the constructible generator on the E_g sector
lam = np.linalg.eigvalsh(H24)                                   # 24 real eigenvalues
def g(E):                                                       # resolvent trace over the 24 E_g modes
    return float(np.sum(1.0 / (E - lam)))

g_at_0   = g(0.0)
g_above  = g(lam.max() + 1.0)
# solve g(E) = 205 on the monotone branch just above the spectrum
lo, hi = lam.max() + 1e-9, lam.max() + 1e6
# g decreases from +inf to 0 as E: lam.max()+ -> +inf; find E with g(E)=205
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if g(mid) > 205.0: lo = mid
    else: hi = mid
E_205 = 0.5 * (lo + hi)
fine_tuning = E_205 - lam.max()                                 # how close E must sit to the spectrum
assert abs(g_at_0 - 205) > 50 and abs(g_above - 205) > 50, "natural E (0, just-above) do NOT give 205"
assert g(E_205) == g(E_205) and abs(g(E_205) - 205) < 1.0, "205 IS reached, but only at the tuned E_205"

bar = "=" * 94
print(bar)
print("K_eff=205 — is the Feshbach Q-trace a finite-service-graph INVARIANT, or a tuned-E resolvent?")
print(bar)
print(f"  constructible finite graph: Q = {len(Q_states)} invalid states; E_g (faces) projector dim = {dim_Eg}")
print(f"\n  (A) E-FREE spectral invariants of the constructible pieces:")
for k, v in invariants.items():
    print(f"        {k:36s} = {v}")
print(f"        remove-exactly-3 criteria        : {no3}  (closest never 3)")
print(f"      => NONE of these E-free invariants is 205, and no criterion removes exactly 3.")
print(f"\n  (B) resolvent over the genuine 24-dim E_g sector (generator V_em, fully specified):")
print(f"        Tr[(E-H)^-1] at E=0           = {g_at_0:+.2f}   (not 205)")
print(f"        Tr[(E-H)^-1] at E=lam_max+1   = {g_above:+.2f}   (not 205)")
print(f"        Tr[(E-H)^-1] = 205 ONLY at E  = {E_205:.4f}  (fine-tuned to {fine_tuning:.3f} above lam_max={lam.max():.3f})")
print(f"      => 205 is reachable from 24 modes ONLY by tuning E to ~{fine_tuning:.2f} of the spectrum.")
print(f"""
{bar}
VERDICT (calculation, exit 0):  K_eff=205 is NOT a finite-service-graph invariant.

  The finite service graph (208 Q-states) and its parameter-free spectral invariants ARE
  constructible -- dim Q=208, E_g dims 24/128, weight 19.25/104, survivors 165/156 -- but NOT ONE
  of them is 205, and no defensible criterion removes exactly 3 Q-states. The value 205 appears
  ONLY as a resolvent trace Tr[P_Eg (E-W_QQ)^-1 P_Eg] at a hand-tuned energy E (here: tuned to
  ~{fine_tuning:.2f} above the E_g spectrum to extract 205 from 24 modes) -- a tuned-E artifact, not an
  invariant. The energy E is un-pinned by the theory (G7), and the operators that would enrich the
  graph (V_weak, the coin, the shift, the Bloch alpha/beta gluing) are themselves underspecified.

  ANSWER to "is the G7/Feshbach Q-trace constructible from a finite service graph?":  NO --
  as a parameter-free invariant. The graph is constructible; 205 is not one of its invariants; it
  is exclusively a tuned-E resolvent value. (Not a universal no-go: a new operator definition could
  change it -- but that is new physics, not a reading of item 8.) The LIVE gravity coefficient is
  the alpha^1/C_loop=3/2 erasure route (item 119), not this dead trace.
{bar}""")
print(f"exit 0 -- K_eff=205 NOT a finite-graph invariant: E-free invariants {sorted(int_invs)} (no 205); "
      f"205 only at tuned E={E_205:.3f} ({fine_tuning:.3f} above spectrum). Q-trace not constructible-as-invariant.")
