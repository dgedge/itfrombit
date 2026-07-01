#!/usr/bin/env python3
r"""foundations_intercell_assembly_layers.py

The inter-cell assembly of the gauge-link web -- its STRUCTURE, clarified into two layers, with an honest
bound on the exact micro crystal.

GEOMETRY (item113): matter octahedra on the SC lattice (Z+1/2)^3; gauge cells (truncated cubes) on Z^3
(the dual). Each gauge cell is surrounded by 8 matter cells (a cube Q3); its 12 gauge links -> the
cuboctahedron cell L(Q3) (the photon's 12-band cell, prior result). KEY sharing fact: each gauge link
(a matter-octahedron bond) is shared by 4 gauge cells -- so the matter-bond web has only 3 links per
PRIMITIVE cell. This forces the two-layer structure:

  MACRO (canon 7.3): the photon we propagate is the scalar Wilson-Maxwell field on the SC matter lattice,
    K(k) = 6 - 2 sum cos(k_i): MASSLESS (K(0)=0), isotropic omega ~ |k| with the (a0 k)^2 correction.
    PLAIN Maxwell -- it carries NO Berry/pi-4 flux. (Adding the Berry flux to a vector link model GAPS
    and MASSES the photon -- see [3] -- so the flux cannot live on the macro photon.)
  MICRO (canon 7.2): the topological 12-band cuboctahedron carrying the derived Berry/pi-4 flux and the
    Chern C_{S7}=-1. This is the line graph of the full truncated-cube edge set; its exact inter-cell
    delta-structure (the conventional-cell bookkeeping of the shared truncated-cube edges) is the genuine
    unrecorded geometric input -- the wall for the exact chiral-edge slab.

So the Berry/pi-4 connection (derived last turn) is a MICRO effect; the chiral topology lives in the
12-band micro, not the macro photon. The inter-cell assembly is thus UNDERSTOOD structurally; the exact
micro crystal needs the truncated-cube full-edge geometry.

  [1] sharing: each gauge link belongs to 4 gauge cells -> 3 links/primitive cell (the macro web).
  [2] macro photon = SC scalar Maxwell K(k): massless + isotropic to (a0 k)^2 (canon 7.3), no Berry flux.
  [3] a matter-bond VECTOR link model carries an intrinsic Weyl monopole (charge +-1), which the Berry
      flux GAPS/masses -> the flux is micro-only (it would unphysically mass the macro photon).

Self-asserting; exit 0. Tier: the two-layer assembly structure + the macro=plain-Maxwell verification +
the Berry-flux-is-micro-only demonstration are SOLID; the exact 12-band micro inter-cell delta-structure
(truncated-cube full edges) remains the unrecorded geometric input for the chiral-edge slab.
"""
import itertools
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== Inter-cell assembly of the gauge-link web: two layers ===\n")

    # [1] each gauge link (matter bond) is shared by 4 gauge cells
    print("[1] link sharing (gauge cells on Z^3, octahedra on (Z+1/2)^3):")
    # an x-link connects octahedra s=(0.5,0.5,0.5) and s+x=(1.5,0.5,0.5); which integer cubes have both as corners?
    s, sx = np.array([0.5, 0.5, 0.5]), np.array([1.5, 0.5, 0.5])
    cubes = []
    for g in itertools.product(range(-1, 3), repeat=3):
        g = np.array(g)
        corners = [g + np.array(o)*0.5 for o in itertools.product((1, -1), repeat=3)]
        if any(np.allclose(s, c) for c in corners) and any(np.allclose(sx, c) for c in corners):
            cubes.append(tuple(g))
    print(f"    the x-link (0.5,0.5,0.5)-(1.5,0.5,0.5) is owned by gauge cells: {cubes}")
    ok(len(cubes) == 4, "each gauge link is shared by 4 gauge cells -> 3 links/primitive cell (the macro web)")

    # [2] macro photon = SC scalar Wilson-Maxwell K(k): massless + isotropic to (a0 k)^2 (canon 7.3)
    print("\n[2] macro photon (canon 7.3) = SC scalar Maxwell K(k)=6-2 sum cos(k_i):")
    K = lambda k: 6 - 2*sum(np.cos(np.array(k)))
    print(f"    K(0,0,0) = {K([0,0,0]):.3f} (massless); K small-k ~ |k|^2:")
    for kk in (0.2, 0.1, 0.05):
        print(f"      |k|={kk}: K={K([kk,0,0]):.5f}, K/|k|^2={K([kk,0,0])/kk**2:.4f} (->1, massless linear photon)")
    ok(abs(K([0, 0, 0])) < 1e-12, "K(0)=0: the macro photon is massless (no Berry-flux gap)")
    # isotropy to (a0 k)^2: K([k,0,0]) vs K(k/sqrt3 each) at fixed |k|
    kmag = 0.3
    k100 = K([kmag, 0, 0]); k111 = K([kmag/np.sqrt(3)]*3)
    ok(abs(k100 - k111)/k100 < 0.02, "macro photon isotropic to (a0 k)^2 (K[100] ~ K[111] at small |k|)")

    # [3] a matter-bond VECTOR link model has an intrinsic Weyl monopole, GAPPED by the Berry flux
    print("\n[3] vector link model: intrinsic Weyl monopole (phi=0), GAPPED by the Berry flux:")
    def Hk(k, phi, t=1.0, tp=0.6):
        kx, ky, kz = k; w = np.exp(1j*phi/3)
        hxy = tp*(1+np.exp(1j*kx))*(1+np.exp(-1j*ky))*w
        hyz = tp*(1+np.exp(1j*ky))*(1+np.exp(-1j*kz))*w
        hzx = tp*(1+np.exp(1j*kz))*(1+np.exp(-1j*kx))*w
        return np.array([[2*t*np.cos(kx), hxy, np.conj(hzx)],
                         [np.conj(hxy), 2*t*np.cos(ky), hyz],
                         [hzx, np.conj(hyz), 2*t*np.cos(kz)]], complex)

    def monopole(band, phi, r=0.15, N=28):
        th = np.linspace(0.02, np.pi-0.02, N); ph = np.linspace(0, 2*np.pi, N, endpoint=False)
        U = np.empty((N, N, 3), complex)
        for i, t_ in enumerate(th):
            for j, p_ in enumerate(ph):
                k = r*np.array([np.sin(t_)*np.cos(p_), np.sin(t_)*np.sin(p_), np.cos(t_)])
                _, v = np.linalg.eigh(Hk(k, phi)); U[i, j] = v[:, band]
        F = 0.0
        for i in range(N-1):
            for j in range(N):
                jp = (j+1) % N
                F += np.angle(np.vdot(U[i, j], U[i, jp])*np.vdot(U[i, jp], U[i+1, jp]) *
                              np.vdot(U[i+1, jp], U[i+1, j])*np.vdot(U[i+1, j], U[i, j]))
        return F/(2*np.pi)
    m0 = [round(monopole(b, 0.0)) for b in range(3)]
    mphi = [round(monopole(b, np.pi/4)) for b in range(3)]
    print(f"    phi=0   (no Berry): band monopole charges = {m0}  (an intrinsic Weyl point)")
    print(f"    phi=pi/4 (Berry)  : band monopole charges = {mphi}  (gapped -> trivial)")
    ok(sorted(m0) == [-1, 0, 1], "bare vector link web carries a Weyl monopole (+-1): intrinsic geometric topology")
    ok(mphi == [0, 0, 0], "the Berry flux GAPS it -> would mass the photon -> the flux is MICRO-only, not on the macro photon")

    print("\n[verdict] INTER-CELL ASSEMBLY = TWO LAYERS (structure understood):")
    print("  - gauge cells on Z^3; each gauge link shared by 4 cells -> 3 links/primitive cell;")
    print("  - MACRO (canon 7.3): the propagating photon = SC scalar Wilson-Maxwell K(k), massless +")
    print("    isotropic to (a0 k)^2, PLAIN (no Berry flux) -- the (a0 k)^2 photon of the whole arc;")
    print("  - MICRO (canon 7.2): the 12-band cuboctahedron carrying the derived Berry/pi-4 flux + the")
    print("    Chern C_{S7}=-1 -- the chiral topology lives HERE, not on the macro photon (adding the flux")
    print("    to a vector link model gaps/masses it, [3]);")
    print("  - the Berry connection (prior turn) is thus a MICRO effect; the exact 12-band micro CRYSTAL")
    print("    (the conventional-cell delta-bookkeeping of the shared truncated-cube edges) is the genuine")
    print("    unrecorded geometric input -- the wall for the exact chiral-edge slab. Assembly understood;")
    print("    macro = plain Maxwell (computed); micro exact-delta remains the gap. exit 0")


if __name__ == "__main__":
    main()
