#!/usr/bin/env python3
"""Item 48 (DRIFT H6) -- BODY-DIAGONAL route (s3.1) survival test. Self-asserting on the COMPUTATION.
Under a GAPLESS bulk web, the meson phi-eigenvalue dissolves on the body-diagonal route too -- three
route-independent reasons + a concrete embedding. CALIBRATION (item 77): HDR exempts the *dynamic*
meson (no rescue) but PROTECTS the *static* baryon, so M_N=2sqrt2*Lambda survives and is NOT flagged.
And two residuals keep the meson at LEANS-negative, not closed: the unpinned leak amplitude t_leak,
and whether the meson's coloured flux leaks into a gapless bulk (s7.3 [-6,6] -> dissolves) or a
confinement-gapped bulk above phi*Lambda (-> survives). The firm result: survival is DYNAMICAL.
"""
import numpy as np

phi = (1 + 5 ** 0.5) / 2

# bulk = body-diagonal web: each cell -> 8 (+-1,+-1,+-1) neighbours (s3.1 shift operator S)
L = 8; xs = range(L); C = [(x, y, z) for x in xs for y in xs for z in xs]; idx = {c: i for i, c in enumerate(C)}; N = len(C)
A = np.zeros((N, N))
for c in C:
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                A[idx[c], idx[((c[0] + sx) % L, (c[1] + sy) % L, (c[2] + sz) % L)]] = 1.0

# (1) wider bulk band (deg 8, dispersion 8 cos cos cos) -> phi deeper inside -> resonance
bw = np.linalg.eigvalsh(A)
print(f"(1) body-diagonal band [{bw[0]:.2f},{bw[-1]:.2f}] (deg 8) vs SC octagon-web [-6,6]: phi={phi:.3f} deeper inside")
assert bw[-1] > 7.8 and bw[0] < -7.8 and bw[0] < phi < bw[-1]

# (2) higher leak coordination (z_leak=6) -> higher Feshbach 2%-survival threshold
thr = lambda z: z / (0.02 * phi); sub = 16 / 9
print(f"(2) z_leak=6 -> 2%-survival needs Delta/t>={thr(6):.0f}; substrate Delta/t={sub:.2f} -> {thr(6) / sub:.0f}x short (octagon: z=3->{thr(3):.0f}, 52x)")
assert abs(thr(6) - 185) < 1 and round(thr(6) / sub) == 104

# (3) phi is not even native to the girth-4 body-diagonal web (the pentagon is octagon-specific)
P3 = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], float)
print(f"(3) girth=4 (no triangles: three +-1 never sum to 0); half-4cycle P3 leading eig={np.linalg.eigvalsh(P3)[-1]:.4f}=sqrt2, NOT phi")
assert abs(np.linalg.eigvalsh(P3)[-1] - 2 ** 0.5) < 1e-9

# (4) even a body-diagonal geodesic that DOES carry phi dissolves at the substrate gap
pj = [idx[(i, i, i)] for i in range(4)]
assert abs(np.linalg.eigvalsh(A[np.ix_(pj, pj)])[-1] - phi) < 1e-6   # isolated geodesic = P4 = phi
print("(4) P4 geodesic on the web (isolated eig=phi, leak z=6):  Delta/t  dev-from-phi  path-weight")
pw = {}
for D in [1.78, 50, 185, 400]:
    H = A.copy()
    for p in pj: H[p, p] += D
    ew, EV = np.linalg.eigh(H); w = (EV ** 2)[pj, :].sum(0)
    j = int(np.argmin(np.abs(ew - (D + phi)))); pw[D] = w[j]
    print(f"       {D:7.1f}    {abs(ew[j] - D - phi) / phi:6.3f}      {w[j]:.3f}")
assert pw[1.78] < 0.20, "at substrate Delta/t=1.78 the phi-mode is delocalized (not a clean phi-state)"
assert pw[400] > 0.90,  "only a >100x gap localizes it -> threshold effect consistent with Delta/t~185"

print("\nBODY-DIAGONAL: under a GAPLESS bulk the meson dissolves here too (104x vs 52x), with no native phi.")
print("Calibration (item 77): HDR exempts the dynamic meson [no rescue] but protects the static baryon [M_N survives].")
print("Residuals (t_leak unpinned; gapless-vs-confining-gapped bulk) keep the meson Item 48 LEANS-negative, NOT closed.")
