#!/usr/bin/env python3
"""
ckm_audit.py
============
Reproducibility artifact for the CKM-sector audit findings (DRIFT M12 addendum +
addendum (b); ANCHOR §15 item 134 demotion). Self-asserting: exit 0 == every
recorded number reproduces.

Part A -- the "unitarity completion" is Anti-pattern 2 (Heron/Jarlskog tautology):
  the standard PDG parametrization is unitary for ANY 4 inputs, so "the other 5
  elements + J match" carries zero independent information; the content is the 4
  inputs only, which miss PDG by up to ~1.6%. Includes the Jarlskog formula check
  (c13^2, not c23^2) and the §16.3 competitor count for the Cabibbo angle.

Part B -- the bare end-to-end walk engine: an honest NEGATIVE result. The 256-dim
  operator is a weighted SUM of CNOTs (NON-unitary, ||W^dag W - I|| ~ 8.3), the
  "mass matrix" M^2=(W^dag W)^2 is nontrivial only because of that, and the bare
  CKM comes out qualitatively wrong (huge mixing in the wrong generation pair).
  The "topological GIM" zero is Gen1-Gen2 (Hamming-2), NOT Gen1-Gen3.

numpy + scipy.
"""
import numpy as np
import scipy.linalg as la

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

# ---------------------------------------------------------------- Part A
print("=" * 72)
print("PART A -- CKM unitarity-completion is Anti-pattern 2 (Heron tautology)")
print("=" * 72)
s13 = (1/3)*np.sin(np.pi/14)**3
c13 = np.sqrt(1 - s13**2)
s12 = np.sin(np.pi/14)/c13
s23 = (np.sqrt(1/28)*np.sin(np.pi/14))/c13
delta = 3*np.pi/8

def ckm(s12, s13, s23, d):
    c12, c13, c23 = [np.sqrt(1-s**2) for s in (s12, s13, s23)]
    e = np.exp(-1j*d)
    return np.array([
        [ c12*c13, s12*c13, s13*e],
        [-s12*c23 - c12*s23*s13*np.conj(e),  c12*c23 - s12*s23*s13*np.conj(e), s23*c13],
        [ s12*s23 - c12*c23*s13*np.conj(e), -c12*s23 - s12*c23*s13*np.conj(e), c23*c13]])

V = ckm(s12, s13, s23, delta); A = np.abs(V)
print("framework |V_ij|:\n", np.round(A, 6))
check(np.linalg.norm(V.conj().T@V - np.eye(3)) < 1e-9, "framework inputs give a unitary matrix")
# the tautology: random/wrong inputs are ALSO unitary -> unitarity tests nothing
randomly_unitary = all(
    np.linalg.norm((Vr := ckm(a, b, c, d)).conj().T@Vr - np.eye(3)) < 1e-9
    for (a, b, c, d) in [(0.3, 0.05, 0.1, 1.0), (0.5, 0.2, 0.3, 2.0), (0.1, 0.01, 0.4, 0.5)])
check(randomly_unitary, "RANDOM/wrong 4 inputs are ALSO unitary -> 'unitarity exact' tests nothing")

# the 9 magnitudes + J are a deterministic function of the 4 inputs (4 numbers, not 9)
PDG = np.array([[0.97435, 0.22501, 0.003732],
                [0.22487, 0.97349, 0.04183],
                [0.00858, 0.04111, 0.999118]])
pct = 100*(A - PDG)/PDG
print("inputs vs PDG-2024 (the ONLY content):"
      f"  |Vus| {pct[0,1]:+.1f}%   |Vcb| {pct[1,2]:+.1f}%   |Vub| {pct[0,2]:+.1f}%")
check(abs(pct[0,1]) > 1 and abs(pct[0,2]) > 1,
      "the 4 inputs (esp. |Vus|, |Vub|) are the WORST-matching elements (>1%)")

# Jarlskog: matrix value matches the correct formula (c13^2 only), NOT the c23^2 typo
J = np.imag(V[0,1]*V[1,2]*np.conj(V[0,2])*np.conj(V[1,1]))
c12, c23 = np.sqrt(1-s12**2), np.sqrt(1-s23**2)
J_correct = s12*c12*s23*c23*s13*c13**2*np.sin(delta)
J_typo    = s12*c12*s23*c23**2*s13*c13**2*np.sin(delta)
check(abs(J - J_correct) < 1e-12, f"Jarlskog matches c13^2 formula (J={J:.4e})")
check(abs(J - J_typo) > 1e-9, "Jarlskog does NOT match the c23^2 typo (ANCHOR item 135 corrected)")

# §16.3: the Cabibbo angle has multiple alphabet expressions within 1%
import math
T = 0.2243
comp = [f"{p}/{q}" for q in range(2, 30) for p in range(1, q) if abs(p/q - T)/T < 0.01]
comp += [f"sin(pi/{n})" for n in range(3, 40) if abs(math.sin(math.pi/n) - T)/T < 0.01]
comp += [f"sqrt(1/{n})" for n in range(3, 40) if abs(math.sqrt(1/n) - T)/T < 0.01]
print(f"§16.3 competitors within 1% of Cabibbo 0.2243: {sorted(set(comp))}")
check(any('2/9' in x for x in comp) and any('14' in x for x in comp),
      "BOTH 2/9 and sin(pi/14) hit |Vus| within 1% (two framework expressions, §16.3)")

# ---------------------------------------------------------------- Part B
print("\n" + "=" * 72)
print("PART B -- bare end-to-end walk CKM engine (honest NEGATIVE result)")
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
unit_err = np.linalg.norm(W.conj().T@W - np.eye(256))
print(f"||W^dag W - I|| = {unit_err:.3f}")
check(unit_err > 1, "the walk operator is NON-unitary (weighted SUM of CNOTs, not a product)")
M2 = (W.conj().T@W) @ (W.conj().T@W)
gens = [(0, 1), (1, 0), (0, 0)]; cols = [(0, 1), (1, 0), (1, 1)]
def gi(g, c, i3): return g[0] | (g[1] << 1) | (1 << 2) | (c[0] << 3) | (c[1] << 4) | (i3 << 5)
up = [gi(g, c, 0) for g in gens for c in cols]
w = np.array([np.exp(-dlt*sum(c)) for c in cols]); w /= np.linalg.norm(w)
P = np.zeros((9, 3));  [P.__setitem__((slice(g*3, g*3+3), g), w) for g in range(3)]
Hu = (P.T @ M2[np.ix_(up, up)] @ P)
print("bare up mass matrix |H_up|:\n", np.round(np.abs(Hu), 4))
hd = lambda a, b: sum(x != y for x, y in zip(a, b))
check(hd(gens[0], gens[1]) == 2 and abs(Hu[0, 1]) < 1e-9,
      "GIM zero is at Gen1-Gen2 (Hamming-2 pair): H_up[0,1]=0")
check(abs(Hu[0, 2]) > 1e-4,
      f"Gen1-Gen3 (Hamming-1) is NONZERO ({abs(Hu[0,2]):.4f}) -- not the 'forbidden' pair")
print("=> bare CKM is qualitatively wrong (suppresses the Cabibbo pair); a NEGATIVE result.")

print("\n" + "=" * 72)
import sys
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- CKM audit reproduces. The 'unitarity completion' is the")
print("Anti-pattern 2 tautology; the bare walk engine is a non-unitary NEGATIVE result.")
print("=" * 72)
