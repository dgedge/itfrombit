#!/usr/bin/env python3
"""Item 48 -- first explicit 3D construction: does the rho's L(P5)=P4 phi-eigenvalue survive
embedding in the 3D TCH gauge web? Self-asserting on the ROBUST facts; the Delta-scan is REPORTED
as model-dependent (see CAVEATS).

Geometry from ANCHOR: gauge on edges / meson = C8 octagon flux, L(P5)=P4 leading eigval phi (L1765);
z=5 TCH bulk coordination (L51/107); SC bulk gauge-web adjacency band [-6,6] (L1084, K=6-2 sum cos);
8 body-diagonal bridges (L358); the framework's claimed survival mechanism is HDR confinement (L2173).

RESULT
  (1) the phi-eigenvector does NOT structurally decouple from the bulk (win-state 1 dead);
  (2) phi=1.618 lies INSIDE the SC band [-6,6] -> a resonance, not a bound state;
  (3) under the bare line-graph embedding the phi-mode DISSOLVES (low path-weight), recovering phi
      only in the large-gap HDR limit -- and even then with a residual downward shift.
=> Item 48's topological / symmetry-selection shortcut is RETIRED (the T_1u bulk leak is symmetry-
   ALLOWED); survival is a DYNAMICAL condition on the confinement gap, NOT a geometric guarantee.

CAVEATS (the Delta-scan threshold + residual shift are MODEL-DEPENDENT, not substrate values):
  - bulk = SC-web proxy (the true line-graph leak is to other edges of L(TCH), not SC centres);
  - t_leak = t_path = 1 (unweighted line graph); the real leak amplitude (from the §3.1 bridges) is
    UNPINNED and is the single biggest lever -- a weaker leak makes survival easier;
  - the confinement gap is modelled as a bulk on-site offset Delta (a simplification of HDR string tension);
  - the 3D path is ASSUMED to be the octagon PERIMETER; Item 48 leaves open whether the 3D flux instead
    traces the §3.1 body-diagonal bridges (a different graph whose leading eigenvalue need not be phi).
  Robust (low-arbitrariness): (1) decoupling-fails and (2) phi-in-band. The Delta-scan (3) is illustrative.
"""
import numpy as np

phi = (1 + 5 ** 0.5) / 2

# (1) harness: isolated P4 = L(P5) leading eigenvalue is phi; Perron (all-positive) eigenvector
A_P4 = np.array([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]], float)
w, V = np.linalg.eigh(A_P4)
u = V[:, -1] * np.sign(V[0, -1])
assert abs(w[-1] - phi) < 1e-9, w[-1]
print(f"[harness] P4 leading eigenvalue = {w[-1]:.4f} = phi ; u_phi = {np.round(u, 4)}")

# (2) structural decoupling test: a bulk channel at octagon vertex v_i couples to (u_{i-1}+u_i).
#     The Perron eigenvector is all-positive, so no adjacent pair can cancel -> no decoupling.
ch = {"v0": u[0], "v1": u[0] + u[1], "v2": u[1] + u[2], "v3": u[2] + u[3], "v4": u[3]}
maxov = max(abs(x) for x in ch.values())
ov = ", ".join(f"{k}:{v:+.3f}" for k, v in ch.items())
print(f"[decouple] overlaps {{ {ov} }} ; max={maxov:.3f}")
assert maxov > 0.1, "exact decoupling would close Item 48; it does NOT hold"
print("[decouple] FAILS -> phi-mode couples to the bulk at every vertex (win-state 1 dead).")

# (3) SC gauge-web band and phi's position in it
L = 6
xs = range(L)
C = [(x, y, z) for x in xs for y in xs for z in xs]
idx = {c: i for i, c in enumerate(C)}
Nb = len(C)
A_sc = np.zeros((Nb, Nb))
for c in C:
    for d in range(3):
        for s in (1, -1):
            cc = list(c); cc[d] = (cc[d] + s) % L
            A_sc[idx[c], idx[tuple(cc)]] = 1.0
bw = np.linalg.eigvalsh(A_sc)
assert bw[0] < phi < bw[-1], "phi must sit inside the SC band for the resonance argument"
print(f"[band] SC web band [{bw[0]:.2f},{bw[-1]:.2f}] ; phi={phi:.3f} INSIDE -> resonance, not bound state.")

# (4) Delta-scan (MODEL-DEPENDENT -- see CAVEATS). Track the eigenvalue nearest +phi and its path-weight.
grp = [[idx[(1, 1, 1)], idx[(1, 1, 2)], idx[(1, 2, 1)]], [idx[(2, 2, 2)], idx[(2, 2, 3)], idx[(2, 3, 2)]],
       [idx[(3, 3, 3)], idx[(3, 3, 4)], idx[(3, 4, 3)]], [idx[(4, 4, 4)], idx[(4, 4, 5)], idx[(4, 5, 4)]]]
t = 1.0
print("[scan] (illustrative)  Delta :  lam(+phi)   dev    path-weight")
for D in [0, 2, 6, 10, 20, 40]:
    H = np.zeros((4 + Nb, 4 + Nb)); H[:4, :4] = A_P4; H[4:, 4:] = A_sc + D * np.eye(Nb)
    for p, g in enumerate(grp):
        for sN in g:
            H[p, 4 + sN] = t; H[4 + sN, p] = t
    ew, EV = np.linalg.eigh(H)
    pw = (EV[:4, :] ** 2).sum(0)
    j = int(np.argmin(np.abs(ew - phi)))
    print(f"               {D:5.1f} : {ew[j]:8.4f}  {abs(ew[j] - phi) / phi:6.3f}    {pw[j]:.3f}")
print("[scan] small gap -> phi dissolves (path-weight ~0.25); re-localizes only at large Delta, shifted below phi.")
print("\nITEM 48: symmetry/topological shortcut RETIRED; survival is a dynamical gap condition, not geometric.")
