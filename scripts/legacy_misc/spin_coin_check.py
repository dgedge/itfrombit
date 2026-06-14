#!/usr/bin/env python3
"""
Verification for §15 item 138 spin-coin assessment.

Checks (all reproducible, no fitting):
  1. E_1/2 (C^2 SU(2)) coin: R(2pi) = -1, R(4pi) = +1  [the core double-cover claim]
  2. The 8 body-diagonal DIRECTION orbit carries +1 under 2pi (single-valued)
  3. The 8-direction perm rep of O_h = A1g + A2u + T1u + T2g (all integer-spin)
  4. §3.5's (alpha_i, beta) are a genuine Cl(3,1) Dirac set  => standard
     Foldy-Wouthuysen reduction gives tree-level g=2 IF the spinor is the
     rotation-locked coin (not internal chi (x) I3).
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

print("CHECK 1 — E_1/2 coin double cover")
def R(theta, n):
    n = np.array(n, float); n = n / np.linalg.norm(n)
    return np.cos(theta/2)*I2 - 1j*np.sin(theta/2)*(n[0]*sx + n[1]*sy + n[2]*sz)
print("  R(2pi) = -I :", np.allclose(R(2*np.pi, [0, 0, 1]), -I2))
print("  R(4pi) = +I :", np.allclose(R(4*np.pi, [0, 0, 1]),  I2))

print("\nCHECK 2 — 8 body-diagonal direction orbit under 2pi")
verts = [np.array([(1 if b & 4 else -1), (1 if b & 2 else -1), (1 if b & 1 else -1)], float)
         for b in range(8)]
def SO3(theta, axis):
    n = np.array(axis, float); n /= np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(theta)*K + (1 - np.cos(theta))*K @ K
M = SO3(2*np.pi, [0, 0, 1])
P = np.zeros((8, 8))
for i, v in enumerate(verts):
    w = M @ v
    j = min(range(8), key=lambda k: np.linalg.norm(verts[k] - w))
    P[j, i] = 1
print("  2pi permutation = identity (value +1) :", np.allclose(P, np.eye(8)),
      " trace =", int(round(np.trace(P))))

print("\nCHECK 3 — perm-rep decomposition (cube vertices under O_h)")
# class order: E 8C3 6C2' 6C4 3C2(=C4^2) i 6S4 8S6 3sigma_h 6sigma_d
order = {'E':1,'8C3':8,'6C2p':6,'6C4':6,'3C2':3,'i':1,'6S4':6,'8S6':8,'3sh':3,'6sd':6}
perm  = {'E':8,'8C3':2,'6C2p':0,'6C4':0,'3C2':0,'i':0,'6S4':0,'8S6':0,'3sh':0,'6sd':4}
# O_h single-valued character table (rows = irreps), same class order
irr = {
 'A1g':{'E':1,'8C3':1,'6C2p':1,'6C4':1,'3C2':1,'i':1,'6S4':1,'8S6':1,'3sh':1,'6sd':1},
 'A2g':{'E':1,'8C3':1,'6C2p':-1,'6C4':-1,'3C2':1,'i':1,'6S4':-1,'8S6':1,'3sh':1,'6sd':-1},
 'Eg' :{'E':2,'8C3':-1,'6C2p':0,'6C4':0,'3C2':2,'i':2,'6S4':0,'8S6':-1,'3sh':2,'6sd':0},
 'T1g':{'E':3,'8C3':0,'6C2p':-1,'6C4':1,'3C2':-1,'i':3,'6S4':1,'8S6':0,'3sh':-1,'6sd':-1},
 'T2g':{'E':3,'8C3':0,'6C2p':1,'6C4':-1,'3C2':-1,'i':3,'6S4':-1,'8S6':0,'3sh':-1,'6sd':1},
 'A1u':{'E':1,'8C3':1,'6C2p':1,'6C4':1,'3C2':1,'i':-1,'6S4':-1,'8S6':-1,'3sh':-1,'6sd':-1},
 'A2u':{'E':1,'8C3':1,'6C2p':-1,'6C4':-1,'3C2':1,'i':-1,'6S4':1,'8S6':-1,'3sh':-1,'6sd':1},
 'Eu' :{'E':2,'8C3':-1,'6C2p':0,'6C4':0,'3C2':2,'i':-2,'6S4':0,'8S6':1,'3sh':-2,'6sd':0},
 'T1u':{'E':3,'8C3':0,'6C2p':-1,'6C4':1,'3C2':-1,'i':-3,'6S4':-1,'8S6':0,'3sh':1,'6sd':1},
 'T2u':{'E':3,'8C3':0,'6C2p':1,'6C4':-1,'3C2':-1,'i':-3,'6S4':1,'8S6':0,'3sh':1,'6sd':-1},
}
g = sum(order.values())  # 48
for name, ch in irr.items():
    n = sum(order[c]*perm[c]*ch[c] for c in order) / g
    if abs(n) > 1e-9:
        print(f"  {name}: multiplicity {n:.0f}")
print("  => A1g + A2u + T1u + T2g (1+1+3+3=8); all single-valued. No 2O spinor irrep present.")

print("\nCHECK 4 — §3.5 matrices are a genuine Dirac set")
ax = [np.kron(sx, sx), np.kron(sx, sy), np.kron(sx, sz)]
beta = np.kron(sz, I2)
anti = all(np.allclose(ax[i]@ax[j] + ax[j]@ax[i], 2*(i == j)*np.eye(4))
           for i in range(3) for j in range(3))
antib = all(np.allclose(ax[i]@beta + beta@ax[i], np.zeros((4, 4))) for i in range(3))
print("  {alpha_i,alpha_j}=2 delta_ij, {alpha_i,beta}=0, beta^2=I :",
      anti and antib and np.allclose(beta@beta, np.eye(4)))
print("  => genuine Dirac set; FW reduction => tree-level g=2 PROVIDED the spinor")
print("     is the rotation-locked coin. Internal chi(x)I3 is not rotation-locked,")
print("     so g=2 is not formulable there. (g=2 is the Dirac baseline, not new.)")
