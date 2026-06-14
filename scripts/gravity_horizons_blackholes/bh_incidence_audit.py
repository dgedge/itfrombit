#!/usr/bin/env python3
"""
bh_incidence_audit.py
=====================
Reproducibility artifact for the ANCHOR §15 item 122 correction (2026-06-05):
the incidence-homomorphism mode IDENTITY was backwards. Self-asserting.

Item 122 pushes the four [8,4,4] CSS X-stabilizers (R1,R2,R3,W) through the 6x8
face->vertex incidence matrix K and claims the 2 survivors are "the two transverse
polarisations of the photon." The linear algebra says the opposite:
  - Phi(R1) = (4,-4,0,0,0,0)  =>  J_x+ = J_x-  (SYMMETRIC), not J_x+ + J_x- = 0;
  - so each R_i removes the ANTISYMMETRIC dipole and keeps the symmetric mode;
  - Phi(W) = (4,...,4)  =>  s_x+s_y+s_z = 0  (TRACELESS), not "div J = 0";
  - the 2 survivors are the symmetric traceless diagonal quadrupole = E_g, which
    by the framework's own assignments (§10.1; lines 1169/1211: E_g=graviton,
    T_1u=photon) is the GRAVITON. The homomorphism projects OUT the T_1u photon
    (antisymmetric dipoles) and KEEPS the E_g graviton.
The 6->2 count stands. The 2/8 = 1/4 is a §16.3 consistency-check, not a
derivation (trivial arithmetic; "transmission ratio = entropy/area" is asserted).

numpy only.
"""
import numpy as np
import itertools

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

oct8 = list(itertools.product([1, -1], repeat=3))           # 8 octant faces
verts = ['X+', 'X-', 'Y+', 'Y-', 'Z+', 'Z-']                # 6 axis-pole vertices
def onside(v, o):
    ax = {'X': 0, 'Y': 1, 'Z': 2}[v[0]]; sg = 1 if v[1] == '+' else -1
    return 1 if o[ax] == sg else 0
K = np.array([[onside(v, o) for o in oct8] for v in verts])  # 6x8
check(all(K[:, c].sum() == 3 for c in range(8)) and all(K[r].sum() == 4 for r in range(6)),
      "incidence K: each face touches 3 vertices, each vertex 4 faces (octahedral)")

R1 = np.array([o[0] for o in oct8]); R2 = np.array([o[1] for o in oct8])
R3 = np.array([o[2] for o in oct8]); Wst = np.ones(8)
phi = {'R1': K@R1, 'R2': K@R2, 'R3': K@R3, 'W': K@Wst}
for nm, p in phi.items(): print(f"  Phi({nm}) = {p.astype(int)}")
check(list((K@R1).astype(int)) == [4, -4, 0, 0, 0, 0],
      "Phi(R1) = (4,-4,0,0,0,0) -> 4 J_x+ - 4 J_x- = 0 -> J_x+ = J_x- (SYMMETRIC, not antisym)")

C = np.array([K@R1, K@R2, K@R3, K@Wst], float)             # 4x6 constraint
u, s, vh = np.linalg.svd(C); rank = int((s > 1e-9).sum())
ns = vh[rank:]
check(6 - rank == 2, f"surviving DOF = 6 - {rank} = 2 (the count is correct)")
# each survivor symmetric (J_i+ = J_i-) and traceless (sum over + poles = 0) -> E_g
sym = all(abs(r[0]-r[1]) < 1e-9 and abs(r[2]-r[3]) < 1e-9 and abs(r[4]-r[5]) < 1e-9 for r in ns)
trc = all(abs(r[0]+r[2]+r[4]) < 1e-9 for r in ns)
check(sym and trc, "survivors are symmetric (J_i+=J_i-) AND traceless = the E_g quadrupole")
# the photon (antisymmetric T_1u dipole) is PROJECTED OUT
dip = [np.array([1, -1, 0, 0, 0, 0.]), np.array([0, 0, 1, -1, 0, 0.]), np.array([1, -1, 1, -1, 1, -1.])]
check(all(np.linalg.norm(C@d) > 1e-9 for d in dip),
      "antisymmetric dipole (T_1u photon) modes are KILLED by the projection")

print("\n  => survivors = E_g (graviton, per §10.1/1169/1211), NOT the T_1u photon.")
print("     Item 122's '2 transverse photon polarisations' is backwards; corrected to E_g graviton.")
print("     2/8 = 1/4 is a §16.3 consistency-check (trivial arithmetic), not a derivation.")

print("\n" + "=" * 64)
import sys
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- the incidence homomorphism keeps the E_g graviton and")
print("projects out the T_1u photon; Item 122's photon identification was backwards.")
print("=" * 64)
