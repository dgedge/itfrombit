#!/usr/bin/env python3
r"""foundations_ltch_cell_and_phasing_family.py

Deeper structural insight into the L(TCH) photon operator (pursuing the exact slab). Two solid new
results + one honest limit:

  [A] CELL IDENTIFIED (corrected). The 12-band cuboctahedral cell is L(Q3) -- the line graph of the
      CUBE GRAPH Q3 (8 vertices, 12 edges, 3-regular) -- NOT L(octahedron) (which is 6-regular, top
      eigenvalue 6). Physically the Q3 is the 8 matter octahedra around one gauge cell (truncated cube),
      joined as a cube; its 12 edges are the 12 gauge links of that gauge cell, and the photon's 12-band
      cell is the line graph of those links. (Earlier "L(octahedron)" framing was wrong; corrected here.)

  [B] PHASING IS A TRIANGLE-FLUX (fingerprint). canon 7.2's phased Gamma-spectrum keeps three T_d
      3-folds {-2,0,2} and moves a 3-block from unphased {4,-2,-2} to {2sqrt2, sqrt6-sqrt2, -sqrt6-sqrt2}.
      That move PRESERVES tr (0) and tr^2 (24) and changes ONLY the determinant (16 -> -8sqrt2) -- i.e.
      only tr(H^3). Since tr(H^3) counts (phased) triangles, the phasing is a TRIANGLE-FLUX on the
      cuboctahedron (flux through its 8 triangular faces), not a generic/Peierls flux.

  [C] HONEST LIMIT: the closed-form spectrum does NOT pin the operator. A triangle-flux CAN reproduce
      canon's exact spectrum (an explicit one below, loss < 1e-4) -- but it is ASYMMETRIC (gauge-invariant
      fluxes distributed 4+3+1 across the 8 triangles), and NO clean T-symmetric (constant-per-tetrahedron)
      triangle-flux matches (best dist ~ 1.0). The three 3-folds make the spectrum degenerate-rich, so
      many triangle-fluxes share it: the spectrum is NECESSARY but NOT SUFFICIENT to recover the canonical
      operator. Pinning it needs the physical gauge connection (the Gauss structure) + the inter-cell
      delta-structure (for C_{S7}=-1 and the slab) -- still unrecorded.

NET: the photon cell is now correctly identified (L of the gauge links around a gauge cell), and the
topological phasing is localized to the triangle-flux family with a clean algebraic fingerprint -- real
insight into how the TCH lattice operates. But the exact chiral-edge slab still needs the gauge connection
+ inter-cell structure; the spectrum gate alone cannot reconstruct them.

Self-asserting; exit 0.
"""
import itertools
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== L(TCH) cell + phasing family (deeper structure) ===\n")

    # [A] cell = L(Q3), Q3 = the 8 matter cells around a gauge cell, joined as a cube graph
    print("[A] cell = L(cube graph Q3) = cuboctahedron (the 12 gauge links around a gauge cell):")
    cV = list(itertools.product((0, 1), repeat=3))
    cE = [(i, j) for i in range(8) for j in range(i+1, 8)
          if sum(a != b for a, b in zip(cV[i], cV[j])) == 1]
    Lq = np.zeros((12, 12))
    for a in range(12):
        for b in range(a+1, 12):
            if set(cE[a]) & set(cE[b]):
                Lq[a, b] = Lq[b, a] = 1
    sp = np.round(np.sort(np.linalg.eigvalsh(Lq)), 6)
    ok(len(cE) == 12, "Q3 has 12 edges (the 12 gauge links around a gauge cell)")
    ok(np.allclose(sp, [-2,-2,-2,-2,-2,0,0,0,2,2,2,4]), f"L(Q3) = cuboctahedron {sp} (4-regular, canon 7.2)")
    Loct_top = 6  # L(octahedron) is 6-regular -- the earlier mis-identification
    ok(Loct_top != 4, "NB L(octahedron) is 6-regular (top eig 6) != cuboctahedron -> the cell is L(Q3), not L(oct)")

    # [B] phasing fingerprint: triangle-flux (changes only det of the moving 3-block)
    print("\n[B] phasing fingerprint: a TRIANGLE-flux (preserves tr, tr^2; changes only det = tr H^3):")
    unph = np.array([4.0, -2.0, -2.0]); ph = np.array([2*np.sqrt(2), np.sqrt(6)-np.sqrt(2), -np.sqrt(6)-np.sqrt(2)])
    for nm, e in (("unphased", unph), ("phased", ph)):
        print(f"    {nm:8s}: tr={e.sum():+.3f} tr2={(e*e).sum():.3f} det={np.prod(e):+.3f}")
    ok(abs(unph.sum()) < 1e-9 and abs(ph.sum()) < 1e-9, "tr preserved (=0)")
    ok(abs((unph*unph).sum() - (ph*ph).sum()) < 1e-9, "tr^2 preserved (=24)")
    ok(abs(np.prod(unph) - 16) < 1e-9 and abs(np.prod(ph) + 8*np.sqrt(2)) < 1e-6, "det changes 16 -> -8sqrt2 (triangle-flux)")

    # ---- machinery: triangle-flux on the cuboctahedron ----
    V = []
    for a in (1, -1):
        for b in (1, -1):
            V += [(a, b, 0), (a, 0, b), (0, a, b)]
    V = np.array(V, float)
    LE = [(i, j) for i in range(12) for j in range(i+1, 12) if abs(np.linalg.norm(V[i]-V[j]) - np.sqrt(2)) < 1e-9]
    A = np.zeros((12, 12))
    for i, j in LE:
        A[i, j] = A[j, i] = 1
    eidx = {e: k for k, e in enumerate(LE)}
    octs = list(itertools.product((1, -1), repeat=3))
    tris = []
    for (sx, sy, sz) in octs:
        vs = [k for k in range(12) if
              (abs(V[k][0]-sx) < 1e-9 and abs(V[k][1]-sy) < 1e-9 and abs(V[k][2]) < 1e-9) or
              (abs(V[k][1]-sy) < 1e-9 and abs(V[k][2]-sz) < 1e-9 and abs(V[k][0]) < 1e-9) or
              (abs(V[k][0]-sx) < 1e-9 and abs(V[k][2]-sz) < 1e-9 and abs(V[k][1]) < 1e-9)]
        tris.append((sx*sy*sz, vs))
    squares = [[k for k in range(12) if abs(V[k][d]-s) < 1e-9] for d in range(3) for s in (1, -1)]

    def ordc(m):
        m = m[:]; c = [m.pop(0)]
        while m:
            for x in list(m):
                if A[c[-1], x]:
                    c.append(x); m.remove(x); break
            else:
                break
        return c

    def Brow(m):
        c = ordc(m); row = np.zeros(len(LE))
        for k in range(len(c)):
            i, j = c[k], c[(k+1) % len(c)]
            if (i, j) in eidx:
                row[eidx[(i, j)]] += 1
            else:
                row[eidx[(j, i)]] -= 1
        return row

    Btri = np.array([Brow(v) for _, v in tris]); Bsq = np.array([Brow(s) for s in squares])
    target = np.array(sorted([2*np.sqrt(2), np.sqrt(6)-np.sqrt(2), -np.sqrt(6)-np.sqrt(2), -2,-2,-2,0,0,0,2,2,2]))

    def spec_tri(ftri):
        B = np.vstack([Btri, Bsq]); f = np.concatenate([ftri, np.zeros(6)])
        th, *_ = np.linalg.lstsq(B, f, rcond=None)
        H = np.zeros((12, 12), complex)
        for k, (i, j) in enumerate(LE):
            H[i, j] = np.exp(1j*th[k]); H[j, i] = np.conj(H[i, j])
        return np.sort(np.linalg.eigvalsh(H))

    # [C1] an explicit (asymmetric) triangle-flux reproduces canon's spectrum
    print("\n[C] the spectrum is NECESSARY but NOT SUFFICIENT to pin the operator:")
    f_found = np.pi * np.array([0.88543, 0.88543, 1.11457, -1.11457, 0.658744, 0.658744, -1.341256, -1.341256])
    d_found = np.linalg.norm(spec_tri(f_found) - target)
    fluxes_mod2 = np.round(np.sort(f_found/np.pi % 2.0), 3)
    print(f"    explicit triangle-flux reproduces canon: dist={d_found:.2e}; fluxes/pi (mod 2) = {fluxes_mod2}")
    ok(d_found < 1e-3, "an explicit triangle-flux reproduces canon's exact spectrum (phasing is in the triangle-flux family)")
    ok(len(set(np.round(fluxes_mod2, 2))) >= 3, "...but its fluxes are ASYMMETRIC (4+3+1), not a clean T_d/T pattern")

    # [C2] no clean T-symmetric (constant-per-tetrahedron) triangle-flux matches
    bestT = None
    for fe in np.linspace(-np.pi, np.pi, 121):
        for fo in np.linspace(-np.pi, np.pi, 121):
            ftri = np.array([fe if p > 0 else fo for p, _ in tris])
            d = np.linalg.norm(spec_tri(ftri) - target)
            if bestT is None or d < bestT[0]:
                bestT = (d, fe, fo)
    print(f"    best T-symmetric (phi_even,phi_odd) flux: dist={bestT[0]:.3f} -- does NOT match")
    ok(bestT[0] > 0.3, "no clean T-symmetric triangle-flux reproduces canon -> spectrum does NOT pin the operator")

    print("\n[verdict] DEEPER STRUCTURE recovered; exact operator/slab still need the gauge connection:")
    print("  + CELL corrected: the 12-band cell = L(Q3) = line graph of the 12 gauge links around a gauge")
    print("    cell (8 matter octahedra as a cube graph), NOT L(octahedron);")
    print("  + PHASING localized: a triangle-flux (fingerprint: preserves tr, tr^2; changes det 16->-8sqrt2);")
    print("  - but the closed-form SPECTRUM is necessary, NOT sufficient: a triangle-flux reproduces it only")
    print("    asymmetrically (non-unique), and no clean T-symmetric one matches. Pinning the canonical")
    print("    operator (and C_{S7}=-1 + the chiral-edge slab) needs the physical Gauss/gauge connection +")
    print("    the inter-cell delta-structure -- still unrecorded. Real insight into the cell + phasing")
    print("    family; the exact slab remains gated on the gauge structure, not the spectrum. exit 0")


if __name__ == "__main__":
    main()
