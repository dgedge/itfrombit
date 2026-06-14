#!/usr/bin/env python3
"""Item 48 — Feshbach self-energy + substrate-constant evaluation (companion to
item48_rho_3d_embedding.py). Self-asserting on the arithmetic; the verdict is REPORTED with its
two input-identifications flagged (t~=Lambda, Delta=sigma*a0) -- those are physical estimates, not
first-principles derivations, and the body-diagonal §3.1 path is unexamined.

Chain: integrating out the SC bulk gives a self-energy Sigma_ii ~= z_leak * t^2 * G_SC(E,0); for a
state above the band, G_SC(E,0) ~ 1/E, so the dressed shift is |dlam| ~= z_leak * t^2 / Delta. A 2%
match (the 760 MeV margin) then needs Delta/t >= z_leak/(0.02*phi). The framework's anchored constants
(sqrt(sigma)=4/3 Lambda §7.17 ; a0 = 1/Lambda) give Delta = sigma*a0 = 16/9 Lambda and t ~ Lambda, i.e.
Delta/t ~ 1.78 -- a ~52x shortfall => the rho mode dissolves under the substrate's own scales.
"""
import numpy as np

phi = (1 + 5 ** 0.5) / 2

# (A) 2%-survival threshold (equal hoppings t): |dlam|/phi = z*t/(Delta*phi) <= 0.02
thr = lambda z: z / (0.02 * phi)
print(f"(A) 2%-survival threshold Delta/t :  z_leak=3 (TCH z=5) -> {thr(3):.1f} ;  z_leak=4 (SC-cube z=6) -> {thr(4):.1f}")
assert abs(thr(3) - 92.7) < 0.6

# (B) substrate evaluation: t ~ Lambda ; Delta = sigma*a0 = (16/9 Lambda^2)(1/Lambda) = 16/9 Lambda
ratio = 16 / 9
print(f"(B) substrate Delta/t = (16/9 Lambda)/(Lambda) = {ratio:.3f}   [sqrt(sigma)=4/3 Lambda §7.17 ; a0=1/Lambda]")
print(f"    shortfall vs requirement:  z=3 -> {thr(3) / ratio:.0f}x ;  z=4 -> {thr(4) / ratio:.0f}x")
assert abs(ratio - 1.778) < 0.01

# (C) the two escape routes (each needs an implausible ~50x effect)
print("(C) escapes to survival (both implausible):")
print(f"    leak suppressed: t_leak/Lambda <= {ratio / thr(3):.3f}  (~{thr(3) / ratio:.0f}x below Lambda; boundary coupling is O(1)/silver-ratio)")
print(f"    gap enhanced:    Delta/Lambda  >= {thr(3):.1f}     (~{thr(3) / ratio:.0f}x above sigma*a0; a single-link gap cannot)")

# (D) confirm: at the substrate Delta/t=1.78 the rho mode is dissolved (dense SC, z_leak=4)
L = 12
xs = range(L)
C = [(x, y, z) for x in xs for y in xs for z in xs]
idx = {c: i for i, c in enumerate(C)}
N = len(C)
A = np.zeros((N, N))
for c in C:
    for d in range(3):
        for s in (1, -1):
            cc = list(c); cc[d] = (cc[d] + s) % L
            A[idx[c], idx[tuple(cc)]] = 1.0
cen = L // 2
path = [idx[(cen + i, cen, cen)] for i in range(4)]
print(f"\n(D) dense SC L={L} (z_leak=4), embedded P4, t=1:")
print("     Delta    lam_rel    dev     path-weight")
for D in [0.0, 1.78, 10.0, 92.7, 123.6, 300.0]:
    H = A.copy()
    for p in path:
        H[p, p] += D
    ew, EV = np.linalg.eigh(H)
    pw = (EV ** 2)[path, :].sum(0)
    j = int(np.argmin(np.abs(ew - (D + phi))))
    print(f"     {D:6.1f}   {ew[j] - D:7.4f}  {abs(ew[j] - D - phi) / phi:6.3f}    {pw[j]:.3f}")

print("\nVERDICT: robustly negative for the octagon-perimeter route under the framework's anchored constants")
print("  (Delta/t~=1.78 vs required ~93). NOT a closed theorem: (i) t~=Lambda and Delta=sigma*a0 are physical")
print("  identifications, not derivations; (ii) the §3.1 body-diagonal path is unexamined. Blast radius is the")
print("  2D-line-graph mass cluster (rho §9.1, N §11.2) -- already §16.3-weak; the structural core is untouched.")
