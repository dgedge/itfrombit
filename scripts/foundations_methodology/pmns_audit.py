#!/usr/bin/env python3
"""
pmns_audit.py
=============
Reproducibility artifact for the PMNS-sector audit (DRIFT M12 addendum (c) +
refinement; ANCHOR §5.10.1 Q4 falsification, §14 row, item 87). Self-asserting.

Part A -- the §5.10.1 Q4 "multiplicative cross-talk" is falsified: U_nu^dag U_L
  with U_nu=TBM and U_L=R23(38)R13(6) gives (th12,th23,th13)=(4.6,82,29), NOT the
  claimed (33,42,8.5). Structural reason: TBM's exact th23=45 + U_L's 38 compose
  toward ~83, never 42. A phase optimiser in the claimed order cannot place th12
  and th23 both physical (so "no phase rescues it" is supported -- but as finite-
  search evidence, NOT a theorem: E4 scope bound).

Part B -- the end-to-end walk-derived U_L gives th23_L ~ 46 (near-maximal, a
  mu-tau M^2 near-degeneracy artifact) and th13_L = 0, reproducing NEITHER the
  claimed 38 nor 6; the TBM product then collapses th23 -> ~1-3 deg.

numpy + scipy.
"""
import numpy as np
import scipy.linalg as la
from scipy.optimize import minimize

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

TBM = np.array([[np.sqrt(2/3), 1/np.sqrt(3), 0],
                [-1/np.sqrt(6), 1/np.sqrt(3), -1/np.sqrt(2)],
                [-1/np.sqrt(6), 1/np.sqrt(3), 1/np.sqrt(2)]], complex)
def R23(t): return np.array([[1, 0, 0], [0, np.cos(t), np.sin(t)], [0, -np.sin(t), np.cos(t)]], complex)
def R13(t, d=0.0):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s*np.exp(-1j*d)], [0, 1, 0], [-s*np.exp(1j*d), 0, c]], complex)
def ang(U):
    s13 = abs(U[0, 2])**2; t13 = np.degrees(np.arcsin(np.sqrt(min(1, s13))))
    t23 = np.degrees(np.arcsin(np.sqrt(min(1, abs(U[1, 2])**2/(1-s13)))))
    t12 = np.degrees(np.arcsin(np.sqrt(min(1, abs(U[0, 1])**2/(1-s13)))))
    return np.array([t12, t23, t13])

print("=" * 72)
print("PART A -- §5.10.1 Q4 cross-talk falsified")
print("=" * 72)
UL = R23(np.radians(38)) @ R13(np.radians(6))
a = ang(TBM.conj().T @ UL)
print(f"U_nu^dag U_L angles (th12,th23,th13) = {np.round(a,2)}   claimed (33,42,8.5)")
check(a[1] > 70 and a[2] > 20, "cross-talk gives th23~82, th13~29 -- violently unphysical, NOT (42,8.5)")
print(f"structural: TBM th23=45 (exact) + U_L th23=38 -> 45+38={45+38} or 45-38={45-38}, never 42")
# phase optimiser in the CLAIMED order cannot place th12 AND th23 both physical
np.random.seed(0)
def Mof(x):
    d, p1, p2 = x; Pd = np.diag([1, np.exp(1j*p1), np.exp(1j*p2)])
    return TBM.conj().T @ Pd @ (R23(np.radians(38)) @ R13(np.radians(6), d))
best = min((minimize(lambda x: np.linalg.norm(ang(Mof(x))-np.array([33., 45., 8.5])),
                     np.random.uniform(0, 2*np.pi, 3), method='Nelder-Mead').fun
            for _ in range(150)))
print(f"best phased fit (claimed order, 150 multistart) to (33,45,8.5): {best:.1f} deg off")
check(best > 3, "phases do NOT rescue the claimed order (supports failure; finite-search, not a theorem)")

print("\n" + "=" * 72)
print("PART B -- bare end-to-end walk-derived U_L")
print("=" * 72)
dlt = 2/9
Ak = np.zeros(8, complex); Ak[0] = np.sqrt(1-dlt)
for k in range(1, 8): Ak[k] = np.sqrt(dlt/7)*np.exp(1j*k*np.pi/4)
W = np.zeros((256, 256), complex)
for k in range(8):
    ctrl, targ = (2-k) % 8, (5-k) % 8; Uk = np.zeros((256, 256))
    for i in range(256):
        if (i >> ctrl) & 1: Uk[i ^ (1 << targ), i] = 1.0
        else: Uk[i, i] = 1.0
    W += Ak[k]*Uk
check(np.linalg.norm(W.conj().T@W - np.eye(256)) > 1, "same non-unitary walk as ckm_audit Part B")
M2 = (W.conj().T@W) @ (W.conj().T@W)
idx = [2+32, 1+32, 0+32]  # e, mu, tau  (LQ=0,I3=1 charged leptons)
Hl = M2[np.ix_(idx, idx)].real
print(f"mu-tau mass block diag: {Hl[1,1]:.4f} vs {Hl[2,2]:.4f} (near-degenerate, diff {abs(Hl[1,1]-Hl[2,2]):.4f})")
check(abs(Hl[1, 1]-Hl[2, 2]) < 0.05, "mu,tau near-degenerate in M^2 (forces ~maximal mixing)")
ev, U23 = la.eigh(Hl[1:, 1:]); UL2 = np.eye(3, dtype=complex); UL2[1:, 1:] = U23
th23L = np.degrees(np.arcsin(min(1, abs(UL2[1, 2])))); th13L = np.degrees(np.arcsin(min(1, abs(UL2[0, 2]))))
print(f"walk-derived U_L: th23_L = {th23L:.1f} deg, th13_L = {th13L:.1f} deg  (framework claims 38, 6)")
check(abs(th23L-38) > 5 and th13L < 1, "walk U_L = (~46, 0), reproducing NEITHER claimed 38 nor 6")
th23_prod = ang(TBM.conj().T @ UL2)[1]
print(f"TBM product collapses th23 -> {th23_prod:.1f} deg (exact value eigenvector-convention noise; collapse robust)")
check(th23_prod < 10, "atmospheric angle collapses in the TBM product")

print("\n" + "=" * 72)
import sys
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- PMNS cross-talk falsified; walk-U_L matches neither claimed angle.")
print("Scope (E4): the simple product fails robustly; 'impossible under all phases' is")
print("finite-search evidence, not a theorem; Koide/TBM tension sharpens open item 87.")
print("=" * 72)
