#!/usr/bin/env python3
"""
hgp_gauss_structure.py
======================
Verifies the *correct* structural part of the HGP-coarse-graining proposal
(DRIFT TQO Target-3 step 4): the HGP([8,4,4]) checks ARE discrete lattice-gauge
operators -- X-checks = Gauss's-law (coboundary) form, Z-checks = Wilson-plaquette
(dual coboundary) form. This is the standard fact that homological CSS codes are
discrete Z_2 gauge theories. Self-asserting: exit 0 == verified.

It ALSO records the decisive negative: the HGP code is a commuting stabilizer
code => gapped => NO massless photon, independent of any geometric embedding /
coarse-graining. So the Gauss/Wilson structure being real does NOT rescue the
photon (a Z_2 Gauss law is gapped; the U(1) Maxwell Gauss law that carries a
photon is gapless). This is why the proposed "Item 144 RESOLVED" is declined.

numpy only.
"""
import numpy as np
import sys

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

H = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
              [0, 0, 0, 0, 1, 1, 1, 1],
              [0, 0, 1, 1, 0, 0, 1, 1],
              [0, 1, 0, 1, 0, 1, 0, 1]])
m, n = H.shape                                       # m=4 checks, n=8 bits
Hx = np.hstack([np.kron(H, np.eye(n, dtype=int)), np.kron(np.eye(m, dtype=int), H.T)]) % 2
Hz = np.hstack([np.kron(np.eye(n, dtype=int), H), np.kron(H.T, np.eye(m, dtype=int))]) % 2

def vv(i, j): return i * n + j                       # VV qubit: internal bit i in spatial cell j
def cc(a, b): return n * n + a * m + b               # CC qubit: internal check a x spatial check b

# X-check (a,j): {VV(i,j): H[a,i]=1}  U  {CC(a,b): H[b,j]=1}   -- Gauss's-law coboundary form
ok_gauss = all(
    set(np.where(Hx[a * n + j])[0]) ==
    (set(vv(i, j) for i in range(n) if H[a, i]) | set(cc(a, b) for b in range(m) if H[b, j]))
    for a in range(m) for j in range(n))
check(ok_gauss, "HGP X-checks = discrete Gauss's-law (coboundary) operators on the 2x2x2 cube")

# Z-check (i,b): {VV(i,j): H[b,j]=1}  U  {CC(a,b): H[a,i]=1}   -- Wilson-plaquette dual form
ok_wilson = all(
    set(np.where(Hz[i * m + b])[0]) ==
    (set(vv(i, j) for j in range(n) if H[b, j]) | set(cc(a, b) for a in range(m) if H[a, i]))
    for i in range(n) for b in range(m))
check(ok_wilson, "HGP Z-checks = macroscopic Wilson-plaquette (dual coboundary) operators")

# decisive negative: commuting stabilizer code => gapped => no massless photon (embedding-independent)
commute = (Hx @ Hz.T % 2 == 0).all()
check(commute, "all stabilizers commute => commuting-projector Hamiltonian => GAPPED (no massless photon)")

print("\n" + "=" * 70)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- the Gauss's-law / Wilson-plaquette structure is REAL (homological")
print("codes ARE discrete Z_2 gauge theories). But the code is a gapped Z_2 stabilizer code,")
print("so it has NO massless photon -- embedding-independent. The structure does not rescue the")
print("photon; a Z_2 Gauss law is gapped, the U(1) Maxwell Gauss law (with a photon) is gapless.")
print("=" * 70)
