#!/usr/bin/env python3
"""Item 119 jump-operator channel, RUN — answering DRIFT G7's posed question:
'does a virtual graviton P->Q->P loop contain exactly one non-unitary syndrome projection,
and what is C_loop?'  Self-asserting; exit 0 = every quoted number verified.

Construction (ANCHOR 15.119): L_k = sqrt(gamma_k) Pi_Q X_k Pi_P, sum gamma_k = alpha*Lambda.
Two defensible gamma supports compared:
  (8-bit)   uniform over all k=0..7 (item 119's literal template);
  (walk-3)  uniform over the bits the canonical walk actually flips, {C0,C1,I3} (2.5: the
            reference hop T_000 flips ONLY C0,C1,I3; the coin flips I3) - the dissipator as
            'the non-unitary shadow of W' (item 119's own characterisation).
Units: Lambda=1, rates per tick. alpha = 1/137.035999.
"""
import itertools, math, numpy as np

alpha = 1/137.035999
NAMES = ["G0","G1","LQ","C0","C1","I3","chi","W"]
def valid(c):
    G0,G1,LQ,C0,C1,I3,chi,W = c
    return not(G0 and G1) and W==chi and ((LQ==0) == (C0==0 and C1==0))
ALL = list(itertools.product((0,1),repeat=8))
P   = [i for i,c in enumerate(ALL) if valid(c)]
Q   = [i for i,c in enumerate(ALL) if not valid(c)]
assert (len(P),len(Q)) == (48,208)
inP = np.zeros(256,bool); inP[P]=True

def flip(i,k):  # X_k as index map on the computational basis
    return i ^ (1 << (7-k))   # bit k in NAMES order = MSB-first
# sanity: flipping W (k=7) toggles the last bit
assert ALL[flip(0,7)][7] == 1 and ALL[flip(0,0)][0] == 1

# exit table: does flipping bit k take codeword i out of P?
exit_k = {i: [not inP[flip(i,k)] for k in range(8)] for i in P}

# ---------- structural fact 1: I3 is the unique free bit; chi/W/LQ always exit ----------
for i in P:
    assert exit_k[i][5] == False                      # I3 flip never exits
    assert exit_k[i][6] and exit_k[i][7] and exit_k[i][2]   # chi, W, LQ always exit
print("fact 1: I3 never exits; chi, W, LQ always exit (verified on all 48)")

# ---------- structural fact 2: one-jump-per-loop is algebraic ----------
# L_k = Pi_Q X_k Pi_P annihilates Q (Pi_P|q> = 0): no Q->P jumps exist in the channel.
# Build dense L_k once (256x256 permutation-projector products) and verify.
def Lmat(k):
    L = np.zeros((256,256))
    for i in P:
        j = flip(i,k)
        if not inP[j]: L[j,i] = 1.0
    return L
Ls = [Lmat(k) for k in range(8)]
PiQ = np.diag((~inP).astype(float))
for k in range(8):
    assert np.linalg.norm(Ls[k] @ PiQ) == 0.0         # L_k Pi_Q = 0 exactly
print("fact 2: L_k Pi_Q = 0 for all k -> NO dissipative Q->P channel exists;")
print("        a P->Q->P virtual loop contains EXACTLY ONE jump event (the return is coherent).")
print("        => the erasure count is alpha^1, not alpha^2.   [G7 question: answered YES]")

# ---------- C_loop under the two supports ----------
def rates(support):
    g = {k: (alpha/len(support) if k in support else 0.0) for k in range(8)}
    return {i: sum(g[k] for k in range(8) if exit_k[i][k]) for i in P}
sup8, sup3 = list(range(8)), [3,4,5]                  # all bits vs {C0,C1,I3}
r8, r3 = rates(sup8), rates(sup3)

vac = ALL.index((0,)*8)                               # the all-zeros codeword (nu_e / vacuum cell)
assert vac in P
print("\nC_loop = alpha / Gamma_vac  (vacuum-cell leakage; K = C_loop/alpha):")
for name, r in (("8-bit uniform", r8), ("walk-3 {C0,C1,I3}", r3)):
    G = r[vac]; C = alpha/G; K = C/alpha
    print(f"  {name:20s} Gamma_vac = {G/alpha:.4f}*alpha   C_loop = {C:.4f}   K = {K:.2f}")
# the walk-3 result is exact 3/2: vacuum is colourless -> C0,C1 flips exit (2 of 3), I3 free
assert abs(alpha/r3[vac] - 1.5) < 1e-12
K3 = 1.5/alpha
Kdata = 206.49                                         # DRIFT G7 (measured M_P in dS variables)
print(f"  walk-3: C_loop = 3/2 EXACTLY (2 exiting channels of 3 walk-active);")
print(f"  K = 3/(2 alpha) = {K3:.2f} vs K_data = {Kdata} -> {abs(K3/Kdata-1)*100:.2f}% (G7's 0.45%)")
assert abs(K3 - 205.55) < 0.01 and abs(K3/Kdata-1) < 0.005

# sector table (the 3/2 is colourless-sector-specific, not universal):
from collections import Counter
sect = Counter()
for i in P:
    c = ALL[i]
    s = ("lepton" if c[2]==0 else f"quark-{('B' if (c[3],c[4])==(1,1) else 'R/G')}")
    sect[(s, round(r3[i]/alpha,4))] += 1
print("  sector breakdown (walk-3, Gamma/alpha):", dict(sect))

# ---------- CPTP + stroboscopic envelope, quantified ----------
rng = np.random.default_rng(3)
tau = 0.05
def channel(rho, tau, Ls, g):
    Lsum = sum(gk * (L.T @ L) for L,gk in zip(Ls,g))
    M0 = np.eye(256) - 0.5*tau*Lsum                  # H_eff = 0 for the trace test
    out = M0 @ rho @ M0.T
    for L,gk in zip(Ls,g):
        out += tau*gk * (L @ rho @ L.T)
    return out, Lsum
g3 = [ (alpha/3 if k in sup3 else 0.0) for k in range(8) ]
# trace preservation order: residual operator = (tau^2/4) Lsum^2
_, Lsum = channel(np.zeros((256,256)), tau, Ls, g3)
res = np.linalg.norm( sum(tau*gk*(L.T@L) for L,gk in zip(Ls,g3)) + ( (np.eye(256)-0.5*tau*Lsum)@(np.eye(256)-0.5*tau*Lsum) ) - np.eye(256) )
res_half = np.linalg.norm( sum((tau/2)*gk*(L.T@L) for L,gk in zip(Ls,g3)) + ( (np.eye(256)-0.25*tau*Lsum)@(np.eye(256)-0.25*tau*Lsum) ) - np.eye(256) )
print(f"\nCPTP: trace-preservation residual ||sum M^t M - I|| = {res:.3e} at tau={tau}, "
      f"{res_half:.3e} at tau/2 (ratio {res/res_half:.2f} ~ 4 => O(tau^2), as canon states)")
assert 3.6 < res/res_half < 4.4

# stroboscopic claim E = exp(tau L): compare on a random density matrix (Taylor exp to converge)
rho0 = rng.standard_normal((256,256)); rho0 = rho0@rho0.T; rho0/=np.trace(rho0)
def lind(rho):
    out = np.zeros_like(rho)
    for L,gk in zip(Ls,g3):
        out += gk*(L@rho@L.T - 0.5*(L.T@L@rho + rho@L.T@L))
    return out
E_rho,_ = channel(rho0, tau, Ls, g3)
exp_rho = rho0.copy(); term = rho0.copy()
for n in range(1,30):
    term = lind(term)*tau/n; exp_rho = exp_rho + term
d1 = np.linalg.norm(E_rho-exp_rho)
E_rho2,_ = channel(rho0, tau/2, Ls, g3)
exp_rho2 = rho0.copy(); term = rho0.copy()
for n in range(1,30):
    term = lind(term)*(tau/2)/n; exp_rho2 = exp_rho2 + term
d2 = np.linalg.norm(E_rho2-exp_rho2)
print(f"stroboscopic: ||E[rho] - exp(tau Lind)[rho]|| = {d1:.3e} at tau={tau}, {d2:.3e} at tau/2 "
      f"(ratio {d1/d2:.2f} ~ 4): the boxed 'E = exp(tau L)' holds to O(tau^2) per tick, not exactly")
assert 3.5 < d1/d2 < 4.5
print("\nALL ASSERTS PASSED")
