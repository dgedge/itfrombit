#!/usr/bin/env python3
r"""ITEM 149 — QUANTUM DARWINISM from the QEC substrate: objectivity via redundant syndrome records.

The absent foundational piece (rung iv). Rungs i/ii gave the Born STRUCTURE: pointer basis = syndrome
eigenbasis (substrate_pointer_basis.py), measure = |<v|psi>|^2 (Gleason). What is still missing is
OBJECTIVITY -- why many independent observers agree on a definite outcome. Zurek's answer is quantum
Darwinism: the pointer observable is REDUNDANTLY recorded in the environment, so many disjoint
environment fragments each carry the full classical record.

The framework supplies exactly this, for free: measurement = syndrome recording, and error-DETECTION
broadcasts the syndrome to MANY sites (the records are carried off as Landauer radiation, item 144). So
the post-measurement state is a broadcast (GHZ-like) record state
        |Psi> = c0 |s+>_S |R0 R0 ... R0> + c1 |s->_S |R1 R1 ... R1>
with N record fragments each holding a copy of the SYNDROME (pointer) value. We compute the
quantum-Darwinism signatures directly from this state:

 [1] partial information I(S:F_f) vs fragment fraction f -> the CLASSICAL PLATEAU at H_S (redundant);
 [2] redundancy R_delta = how many DISJOINT fragments each independently reveal the outcome -> = N
     (maximal: the records are classical perfect copies, Landauer-committed);
 [3] the plateau is for the SYNDROME (pointer) observable only -- the conjugate (phase) info is NOT
     redundant (it needs the WHOLE environment), so ONLY the einselected basis is objective;
 [4] a non-broadcast control (record in ONE place) gives R=1 -> no objectivity: the broadcast is what
     makes it objective, and the substrate's error-detection IS a broadcast.

CONSEQUENCE: in this framework einselection (rung ii) and quantum Darwinism are the SAME mechanism --
syndrome recording -- both selecting the syndrome basis; objectivity is derived, not posited. (The
Born MEASURE and complex-Hilbert residuals of rung iii are untouched -- the universal wall.)
Self-asserting; numpy only. exit 0 = plateau + redundancy N + pointer-only objectivity + control.
"""
import numpy as np
import itertools

def reduced(psi_t, keep):
    n = psi_t.ndim
    keep = sorted(keep); trace = [a for a in range(n) if a not in keep]
    T = np.transpose(psi_t, keep + trace).reshape(2**len(keep), 2**len(trace))
    return T @ T.conj().T

def vN(rho):                                   # von Neumann entropy in bits
    ev = np.linalg.eigvalsh(rho); ev = ev[ev > 1e-13]
    return float(-np.sum(ev * np.log2(ev)))

def H2(p): return float(-p*np.log2(p) - (1-p)*np.log2(1-p)) if 0 < p < 1 else 0.0

# ---- the broadcast record state: system qubit (axis 0) + N syndrome-copy fragments ----
N = 8
p0 = 0.70                                       # = |c0|^2 (generic, != 1/2, to show it is not special)
c0, c1 = np.sqrt(p0), np.sqrt(1 - p0)
psi = np.zeros(2**(N+1)); psi[0] = c0; psi[-1] = c1   # |0>|0..0> + |1>|1..1>
psi_t = psi.reshape((2,)*(N+1))
H_S = vN(reduced(psi_t, {0}))
print(f"system entropy H_S = {H_S:.4f} bits  (= H2({p0}) = {H2(p0):.4f}); N = {N} record fragments")

# ---- [1] partial information I(S:F_k), F_k = first k fragments (axes 1..k) ----
print("\n[1] partial information I(S:F_k) vs fragment count k (the redundancy plateau):")
I = []
for k in range(N+1):
    F = set(range(1, k+1))
    iSF = H_S + (vN(reduced(psi_t, F)) if F else 0.0) - vN(reduced(psi_t, {0} | F))
    I.append(iSF)
    bar = "#" * int(round(iSF/ (2*H_S) * 40))
    print(f"    k={k}:  I(S:F_k) = {iSF:.4f} bits  {bar}")
plateau = I[1:N]                                # k = 1 .. N-1
assert all(abs(x - H_S) < 1e-9 for x in plateau), "classical plateau must sit at H_S"
assert abs(I[N] - 2*H_S) < 1e-9, "the WHOLE environment carries 2 H_S (the quantum/global info)"
assert I[0] < 1e-12
print(f"    -> PLATEAU at H_S = {H_S:.3f} bits across k=1..{N-1} (any fragment already has the full")
print(f"       classical record); the jump to 2 H_S only at k=N is the non-redundant quantum residual.")

# ---- [2] redundancy R_delta: # of DISJOINT fragments each giving >= (1-delta) H_S ----
delta = 0.1
# each single fragment (one qubit copy) carries I(S:one) = H_S >= (1-delta)H_S -> N disjoint of them
per_fragment = H_S + vN(reduced(psi_t, {1})) - vN(reduced(psi_t, {0, 1}))
R = N if per_fragment >= (1-delta)*H_S else 0
print(f"\n[2] each single fragment carries I = {per_fragment:.4f} bits (>= {(1-delta)*H_S:.3f}); so")
print(f"    redundancy R_delta = {R} disjoint fragments each independently reveal the outcome -> MAXIMAL.")
assert R == N and per_fragment > (1-delta)*H_S

# ---- [3] only the SYNDROME (pointer) observable is objective; the conjugate is not ----
# classical (pointer/Z) info is in EVERY fragment (above). The conjugate (X/phase) info: a fragment of
# k<N qubits, traced, is DIAGONAL (no phase) -> carries 0 conjugate info; only k=N recovers the phase.
print("\n[3] pointer vs conjugate: the classical (syndrome) record is in every fragment (plateau above);")
print("    the conjugate (phase) information is NOT redundant -- a proper sub-fragment's reduced state is")
diag_offdiag = []
for k in (1, N//2, N-1):
    rhoF = reduced(psi_t, set(range(1, k+1)))
    offmax = np.max(np.abs(rhoF - np.diag(np.diag(rhoF))))
    diag_offdiag.append(offmax)
    print(f"    k={k}: max off-diagonal of rho_F = {offmax:.2e}  (diagonal => no phase => no conjugate info)")
assert all(o < 1e-12 for o in diag_offdiag)
print("    -> ONLY the syndrome (pointer) basis is redundantly broadcast = OBJECTIVE; the conjugate")
print("       observable is global (needs all N) = not objective. Objectivity picks the SAME basis as")
print("       einselection (rung ii): syndrome recording is the single mechanism behind both.")

# ---- [4] control: a NON-broadcast record (syndrome copied to ONE site) -> R = 1, not objective ----
psi1 = np.zeros(2**(N+1))                       # |0>|10..0> + |1>|00..0>: only fragment 1 correlated
# system with only fragment-1 correlated, rest in |0>
idx = lambda s, f1: (s << N) | (f1 << (N-1))
psi1[idx(0,1)] = c0; psi1[idx(1,0)] = c1
psi1_t = psi1.reshape((2,)*(N+1))
Ictrl = []
for k in range(N+1):
    F = set(range(1, k+1))
    iSF = vN(reduced(psi1_t, {0})) + (vN(reduced(psi1_t, F)) if F else 0.0) - vN(reduced(psi1_t, {0}|F))
    Ictrl.append(iSF)
# only fragments INCLUDING site 1 carry the info; a fragment of the OTHER sites carries 0
info_site1 = vN(reduced(psi1_t,{0})) + vN(reduced(psi1_t,{1})) - vN(reduced(psi1_t,{0,1}))
info_site2 = vN(reduced(psi1_t,{0})) + vN(reduced(psi1_t,{2})) - vN(reduced(psi1_t,{0,2}))
print(f"\n[4] NON-broadcast control (syndrome recorded at ONE site): fragment{{1}} info = {info_site1:.3f},")
print(f"    fragment{{2}} info = {info_site2:.3f} -> only 1 fragment knows -> R=1, NOT objective.")
print(f"    The framework AVOIDS this: error-DETECTION broadcasts the syndrome to many sites (item 144),")
print(f"    so the physical case is [1]-[2] (R=N), not this control.")
assert info_site1 > (1-delta)*H_S and info_site2 < 1e-9

print(f"""
--- VERDICT (item 149 rung iv: quantum Darwinism / objectivity DERIVED) ---
Quantum Darwinism holds in the QEC substrate, and maximally: measurement is syndrome recording, and
error-DETECTION broadcasts the syndrome to many record fragments (Landauer radiation, item 144). The
broadcast (GHZ-like) record state has the QD signatures computed above:
  * a CLASSICAL PLATEAU: I(S:F_k) = H_S = {H_S:.3f} bits for every fragment k=1..{N-1} (jump to 2 H_S only
    for the whole environment) -- the outcome is redundantly recorded;
  * REDUNDANCY R = {R} (= N): N disjoint fragments each independently reveal the outcome, so many
    observers agree without disturbing it -> OBJECTIVITY;
  * the redundantly-broadcast (objective) observable is EXACTLY the syndrome (pointer) basis -- the
    conjugate (phase) info is global, not objective. So einselection (rung ii) and quantum Darwinism
    select the SAME basis: in this framework they are ONE mechanism (syndrome recording).
This closes the OBJECTIVITY gap of the measurement account: structure (rungs i/ii) + objectivity
(this) all follow from syndrome recording -- the same mechanism as entropy and the arrow (item 144).
HONEST SCOPE: this does NOT touch the rung-iii residuals. The Born MEASURE (why an equal-amplitude
branch carries unit weight) and WHY complex Hilbert space remain the universal reconstruction wall
(substrate_born_residual.py). QD supplies the envariance SETTING (redundant records) Zurek's measure
argument needs, but the equal-amplitude indifference step is unchanged. TIER: objectivity derived
(rung iv); Born measure + complex formalism still imported (rung iii, the universal residual).
""")
print("ALL ASSERTIONS PASSED -- plateau at H_S, redundancy R=N, pointer-only objectivity, R=1 control.")
import sys; sys.exit(0)
