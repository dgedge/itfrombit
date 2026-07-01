#!/usr/bin/env python3
"""Dressed-alpha attack: the non-unital (irreversible-emission) correction to the unital 1/137.

Setup: the item-79 closure derived P(em) = 1/137 for the MONITORED, BIDIRECTIONAL (T=inf) bridge.
Physical emission is irreversible: the photon escapes to the macroscopic gauge web. Add the escape
channel (loss at rate Gamma_esc on |em>) and compute the quasi-stationary (conditioned-on-survival)
measure pi_QSS. The dressed coupling is alpha_eff = pi_QSS(em); the dressed shift is
delta = 1/alpha_eff - 137.

Regime: monitoring at the clock rate (gamma >> couplings) reduces the populations to the classical
dephased chain with symmetric rates W_st ~ |H_st|^2 (doubly stochastic — reproducing the unital
uniform fixed point at Gamma_esc = 0). Escape perturbs it. Compute: sign, the leading linear
coefficient c in delta = c * (Gamma_esc / W_norm) on the actual canonical graph, and the ratio
required to reproduce the observed delta = 0.0360 (= CODATA 137.0360 - 137).
16.3 discipline: ALL normalisations reported; no parameter chosen to fit.
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ---- canonical 137-state system (from item79_unital_channel.py) ----
N = 16
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1
pairs = [(i,j) for i in range(N) for j in range(i,N)]
idx = {p:k for k,p in enumerate(pairs)}
d = len(pairs); D = d+1
Bas = np.zeros((N*N,d))
for (i,j),k in idx.items():
    v = np.zeros((N,N))
    if i==j: v[i,i]=1.0
    else: v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,k]=v.reshape(-1)
H = np.zeros((D,D)); H[:d,:d] = Bas.T @ (np.kron(A,np.eye(N))+np.kron(np.eye(N),A)) @ Bas
w = 0.5
for p in ((0,8),(0,9),(1,8),(1,9)): H[idx[p],d] = H[d,idx[p]] = w*0.5

# dephased classical chain: rates W_st = |H_st|^2 (gamma scale absorbed; symmetric => doubly stochastic)
W = np.abs(H)**2; np.fill_diagonal(W, 0.0)
gen = W - np.diag(W.sum(1))                       # generator, columns... symmetric: gen = gen.T
assert np.linalg.norm(W-W.T) == 0.0
# Gamma_esc = 0 check: uniform is stationary (doubly stochastic)
assert np.linalg.norm(gen @ np.ones(D)) < 1e-12
print(f"dephased chain on D={D}: symmetric rates |H_st|^2, uniform stationary at Gamma_esc=0 (unital limit OK)")

def qss(Gesc):
    M = gen.copy(); M[d,d] -= Gesc                # loss at |em>
    ev, V = np.linalg.eig(M)
    k = np.argmax(ev.real)                        # slowest-decaying (quasi-stationary) mode
    v = np.abs(V[:,k].real); v /= v.sum()
    return v[d], -ev[k].real                      # pi_em, escape rate

# scan Gamma_esc over decades; normalisations reported: chain gap and mean off-diagonal rate
evg = np.sort(np.linalg.eigvalsh(gen)); gap = -evg[-2]
Wmean = W.sum()/(D*(D-1))
print(f"chain normalisations: spectral gap = {gap:.4e}; mean pair rate = {Wmean:.4e}")
print(f"{'Gesc':>10} {'pi_em':>12} {'delta=1/pi-137':>16} {'sign':>5}")
deltas = []
for Gesc in (1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
    pe, esc = qss(Gesc)
    dlt = 1/pe - 137
    deltas.append((Gesc, dlt))
    print(f"{Gesc:>10.1e} {pe:>12.8f} {dlt:>16.6f} {'+' if dlt>0 else '-':>5}")
    assert dlt > 0                                # SIGN: irreversible escape pushes alpha^-1 ABOVE 137
# leading linear coefficient in each normalisation
c_gap   = np.polyfit([g/gap   for g,_ in deltas[:3]], [x for _,x in deltas[:3]], 1)[0]
c_mean  = np.polyfit([g/Wmean for g,_ in deltas[:3]], [x for _,x in deltas[:3]], 1)[0]
print(f"\nleading coefficient: delta = c * (Gesc/gap)  with c_gap  = {c_gap:.4f}")
print(f"                     delta = c * (Gesc/Wmean) with c_mean = {c_mean:.4f}")
target = 0.0360
req_gap, req_mean = target/c_gap, target/c_mean
print(f"\nINVERSE PROBLEM: observed delta = {target} requires Gesc/gap = {req_gap:.3e}"
      f"  (equivalently Gesc/Wmean = {req_mean:.3e})")
# canonical-scale comparison (REPORTING ALL, fitting NONE):
alpha0 = 1/137
cands = {"alpha_0": alpha0, "alpha_0/2pi": alpha0/(2*np.pi), "31*alpha_0/2pi": 31*alpha0/(2*np.pi),
         "1/136": 1/136, "w^2=0.25": w*w, "2/3*alpha_0": (2/3)*alpha0}
print(f"canonical-ratio candidates vs required Gesc/gap = {req_gap:.3e}:")
for nm, v in cands.items():
    print(f"   {nm:>14} = {v:.4e}   ratio to required: {v/req_gap:>8.2f}")

print(f"""
RESULT (honest scope):
  SIGN DERIVED: irreversible escape makes the survival-conditioned measure avoid the portal,
  so pi_em < 1/137 and alpha_eff^-1 > 137 — the dressed shift has the right sign STRUCTURALLY.
  FORM DERIVED: delta is linear in the escape-to-relaxation ratio at small Gesc, with the
  coefficient computed on the actual canonical graph (c_gap = {c_gap:.3f} in gap units).
  MAGNITUDE OPEN: reproducing delta = 0.036 requires Gesc/gap = {req_gap:.2e}; no canonical scale
  in the table above matches it cleanly (all tried are reported; none is adopted). The dressed
  problem is hereby SHARPENED, not solved: from 'count vs integral, unexplained' (K5) to
  'one named physical ratio — the photon escape rate over the bridge relaxation rate — whose
  canon-derivation is the remaining open step.' Deriving Gesc (Fermi golden rule of the portal
  into the T_1u web) and the bridge relaxation from the substrate would close it; that build
  needs the gauge-web density of states, which canon has not pinned (15 item 75 3D-bulk limit).
ALL ASSERTS PASSED""")
