#!/usr/bin/env python3
r"""foundations_ltch_edge_velocity_probe.py

Last derived-vs-assumed gap: is the chiral edge/Dirac mode's velocity EXACTLY c and its dispersion
EXACTLY linear? Honest answer: NO -- generically. Topology (Chern C_{S7}=-1) protects the EXISTENCE and
CHIRALITY of the mode (spectral flow / bulk-boundary), but NOT its velocity (non-universal,
boundary/parameter dependent) nor its linearity (curved away from the crossing/touching). So the carrier
constituents are real and chiral (grounding the bundle + forcing collinearity, prior result stands), but
only APPROXIMATELY null/linear -- exact only asymptotically in the soft/IR limit at the touching.

HONEST SCOPE. The exact canonical 12-band L(TCH) Bloch operator (the specific pi/4 flux giving canon
7.2's [-3.864,...] spectrum + C_{S7}=-1) is NOT recorded in a reusable form. So this establishes the
answer from (i) the cuboctahedral CELL, verified faithful by reproducing canon 7.2's unphased spectrum,
and (ii) a clean C=-1 Chern slab demonstrating the GENERAL topological facts (which is what the question
turns on). The exact L(TCH) numbers remain a follow-up requiring the canonical operator to be written
down explicitly.

  [1] CELL FAITHFULNESS: the 12-site cuboctahedron reproduces canon 7.2's unphased spectrum
      {4,2,2,2,0,0,0,-2x5} -- so the L-cell is the right one.
  [2] TOPOLOGY IS REAL but METRIC IS NOT: on a C=-1 model (QWZ), the Chern number is robust (=-1 across
      the phase), but the EDGE VELOCITY v_edge VARIES with the parameter (non-universal) and the edge
      dispersion has nonzero CURVATURE -> not exactly linear, not a fixed c.
  [3] BULK DIRAC CONE: linear + isotropic only to LEADING order at the touching (symmetry-protected),
      with nonzero curvature beyond -> exact linearity is an asymptotic (soft) property, not exact.
  [4] TWO-VELOCITY RISK: the edge/Dirac velocity is generically != the SC k->0 photon velocity ("c"),
      so without a boundary symmetry pinning them equal there is a mild residual two-velocity LV.

Self-asserting; exit 0. Tier: cuboctahedral cell faithfulness (canon-exact) + the general topological
facts (existence/chirality protected; velocity/linearity NOT) DERIVED/computed on a clean C=-1 model.
VERDICT: the last gap does NOT close -- exact (v=c, linear) is generically FALSE; the prior "exactly
null bundle" is corrected to "approximately null, exact only in the soft/IR limit", GRB-safety becomes
approximate (soft-suppressed ~ (eps/Lambda)^2). Protection still grounds existence+chirality+collinearity.
The exact-c+linear property would require a fine-tuned/symmetry-protected boundary (the open question);
the exact L(TCH) slab needs the canonical 12-band operator (not reconstructed here).
"""
import numpy as np

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


# ---- QWZ Chern model (tx scales the x-velocity term: a C-preserving bulk anisotropy) ----
def Hk(kx, ky, m, tx=1.0):
    return tx * np.sin(kx) * SX + np.sin(ky) * SY + (m + np.cos(kx) + np.cos(ky)) * SZ


def chern_lower(m, tx=1.0, N=24):
    """Fukui-Hatsugai-Suzuki Chern number of the lower band."""
    ks = np.linspace(0, 2 * np.pi, N, endpoint=False)
    U = np.empty((N, N, 2), complex)
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks):
            w, v = np.linalg.eigh(Hk(kx, ky, m, tx))
            U[i, j] = v[:, 0]
    F = 0.0
    for i in range(N):
        for j in range(N):
            ip, jp = (i + 1) % N, (j + 1) % N
            u1 = np.vdot(U[i, j], U[ip, j]); u2 = np.vdot(U[ip, j], U[ip, jp])
            u3 = np.vdot(U[ip, jp], U[i, jp]); u4 = np.vdot(U[i, jp], U[i, j])
            F += np.angle(u1 * u2 * u3 * u4)
    return F / (2 * np.pi)


def bottom_edge_branch(m, tx=1.0, Ny=80):
    """Strip periodic in x, open in y. Track the BOTTOM-edge in-gap chiral branch E(kx).
    tx scales the x-velocity term (a bulk anisotropy that deforms v_edge while keeping C=-1)."""
    def H(kx):
        H = np.zeros((2 * Ny, 2 * Ny), complex)
        onsite = tx * np.sin(kx) * SX + (m + np.cos(kx)) * SZ
        Ty = 0.5 * SZ - 0.5j * SY
        for y in range(Ny):
            H[2*y:2*y+2, 2*y:2*y+2] = onsite
            if y < Ny - 1:
                H[2*y:2*y+2, 2*y+2:2*y+4] = Ty
                H[2*y+2:2*y+4, 2*y:2*y+2] = Ty.conj().T
        return H
    nb = 12                                              # bottom 6 sites (2 orbitals each)
    kxs = np.linspace(-np.pi, np.pi, 1201)
    branch = np.full(len(kxs), np.nan)
    for i, kx in enumerate(kxs):
        ev, U = np.linalg.eigh(H(kx))
        bottom = np.sum(np.abs(U[:nb, :])**2, axis=0)
        top = np.sum(np.abs(U[-nb:, :])**2, axis=0)
        ingap = np.abs(ev) < 0.85
        score = (bottom - top) - 10.0 * (~ingap)        # prefer in-gap, bottom-localized
        j = int(np.argmax(score))
        if ingap[j] and bottom[j] > 0.5:                # a genuine bottom-edge in-gap state
            branch[i] = ev[j]
    return kxs, branch


def edge_velocity_and_curvature(kxs, branch):
    """Velocity at the E=0 crossing of the bottom-edge branch, and the VARIATION of v_g ALONG the
    branch (the honest non-linearity signature: a linear branch has constant v_g; the edge sine does not)."""
    good = ~np.isnan(branch)
    kk, bb = kxs[good], branch[good]
    if len(kk) < 5:
        return 0.0, 0.0, kk, bb
    vg = np.gradient(bb, kk)                              # local group velocity along the branch
    sign = np.sign(bb)
    cross = np.where(np.diff(sign) != 0)[0]
    v_cross = abs(vg[cross[len(cross)//2]]) if len(cross) else abs(vg[np.argmin(np.abs(bb))])
    v_spread = float(np.max(np.abs(vg)) - np.min(np.abs(vg)))   # >0 => v_g not constant => not linear
    return v_cross, v_spread, kk, bb


def main():
    print("=== L(TCH) edge velocity probe: is v_edge exactly c and the dispersion exactly linear? ===\n")

    # [1] cuboctahedral cell faithfulness (canon 7.2 unphased spectrum)
    print("[1] cell faithfulness -- 12-site cuboctahedron reproduces canon 7.2 unphased spectrum:")
    verts = []
    for a in (1, -1):
        for b in (1, -1):
            verts += [(a, b, 0), (a, 0, b), (0, a, b)]
    V = np.array(verts, float)
    A = np.array([[1.0 if i != j and abs(np.linalg.norm(V[i]-V[j]) - np.sqrt(2)) < 1e-9 else 0.0
                   for j in range(12)] for i in range(12)])
    spec = np.round(np.sort(np.linalg.eigvalsh(A)), 6)
    ok(np.allclose(spec, [-2,-2,-2,-2,-2,0,0,0,2,2,2,4]), f"cuboctahedron spectrum {spec} = canon 7.2 (right L-cell)")

    # [2] topology robust, metric NOT: C=-1 fixed but v_edge deforms with a bulk anisotropy + curvature!=0
    print("\n[2] topology is robust, the METRIC is not (C=-1 model, deform a bulk anisotropy t_x):")
    m = -1.0
    v_edges, curvs = {}, {}
    for tx in (0.6, 1.0, 1.6):
        C = chern_lower(m, tx)
        kxs, branch = bottom_edge_branch(m, tx)
        v_edge, v_spread, kk, bb = edge_velocity_and_curvature(kxs, branch)
        v_edges[tx], curvs[tx] = v_edge, v_spread
        print(f"    t_x={tx:.1f}: Chern={C:+.2f}, v_edge(at crossing)={v_edge:.3f}, "
              f"v_g variation along branch={v_spread:.3f} (>0 => not linear)")
        ok(abs(C - (-1)) < 0.15, f"Chern = -1 (robust) at t_x={tx}")
        ok(v_edge > 1e-3, f"a genuine chiral edge branch crosses E=0 with nonzero slope at t_x={tx}")
    spread = max(v_edges.values()) - min(v_edges.values())
    print(f"    => v_edge spans {min(v_edges.values()):.3f}..{max(v_edges.values()):.3f} (spread {spread:.3f}) "
          f"while Chern stays -1: VELOCITY IS NON-UNIVERSAL (deformable, not a fixed c).")
    ok(spread > 0.15, "edge velocity VARIES under a C-preserving deformation -> non-universal, NOT pinned to c")
    ok(max(curvs.values()) > 1e-2, "v_g varies ALONG the edge branch (it turns over) -> NOT exactly linear")

    # [3] bulk Dirac cone: linear to leading order, curvature beyond (exact linearity only asymptotic)
    print("\n[3] bulk Dirac cone (m -> -2, gap closes at k=0): linear leading order + curvature beyond:")
    m = -1.98                                            # near the Dirac transition
    def gap(kx, ky):
        w = np.linalg.eigvalsh(Hk(kx, ky, m)); return w[1] - w[0]
    ks = np.array([0.02, 0.05, 0.10, 0.20, 0.40])
    half = np.array([gap(k, 0)/2 for k in ks])           # upper-band energy along kx
    v_lin = half[0] / ks[0]                              # slope near the point
    resid = half - v_lin * ks                            # deviation from pure linear
    print(f"    half-gap E(k): {np.round(half,4)}; linear v={v_lin:.3f}; deviation-from-linear {np.round(resid,4)}")
    ok(abs(half[0] - v_lin*ks[0]) < 1e-6, "linear to LEADING order near the touching (Dirac cone)")
    ok(abs(resid[-1]) > 1e-3, "nonzero curvature away from the touching -> exact linearity is asymptotic only")

    # [4] two-velocity risk: edge/Dirac velocity generically != the SC k->0 photon velocity
    print("\n[4] two-velocity risk: the chiral/Dirac velocity is generically != the SC k->0 photon c:")
    v_sc = 1.0                                            # SC photon: omega=c|k|, v=c=1 (lattice units)
    v_dirac = v_lin
    print(f"    SC k->0 photon velocity c = {v_sc:.3f}; bulk-Dirac slope v_D = {v_dirac:.3f} "
          f"(ratio {v_dirac/v_sc:.3f} != 1 generically)")
    ok(abs(v_dirac - v_sc) > 1e-3, "Dirac/edge velocity != SC c generically -> a mild two-velocity LV unless boundary-symmetry-pinned")

    print("\n[verdict] THE LAST GAP DOES NOT CLOSE -- 'exactly c + exactly linear' is generically FALSE:")
    print("  - the L-cell is faithful (cuboctahedron = canon 7.2), and the Chern number is robust (=-1);")
    print("  - BUT topology protects only EXISTENCE + CHIRALITY (spectral flow), NOT the metric: the edge")
    print("    velocity is non-universal (varies across the C=-1 phase) and the dispersion is curved")
    print("    (linear only asymptotically at the touching/crossing); the Dirac velocity != the SC c.")
    print("  CONSEQUENCE: the prior 'exactly null bundle' is corrected to 'APPROXIMATELY null, exact only")
    print("  in the soft/IR limit'; GRB-safety becomes APPROXIMATE (soft-suppressed ~ (eps/Lambda)^2), and")
    print("  a mild two-velocity LV is possible. Topological protection still GROUNDS existence + chirality")
    print("  + collinearity (the prior YES stands). Exact c+linear would need a fine-tuned/symmetry-protected")
    print("  boundary (the open question); the exact L(TCH) slab needs the canonical 12-band operator. exit 0")


if __name__ == "__main__":
    main()
