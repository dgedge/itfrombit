#!/usr/bin/env python3
r"""foundations_causet_dirac_slash.py

WAYS FORWARD with the causal-set Dirac construction (the trans-Lambda_QCD frontier "Dirac spin needs a spin
frame"). Refines the 2026-06-25 "order-alone no-go bypassed" note: the Dirac operator is
D = gamma^a e_a^mu d_mu + m -- THREE pieces, and the code supplies only the first:
  gamma^a  : the Clifford algebra            -- code-supplied (§3.5), verified here.
  e_a^mu   : the vierbein / frame field       -- trivial (global) in flat space; the deep piece in curved.
  d_mu     : a Lorentz-invariant directional derivative on the sprinkling -- the genuine open content.

PROBE of Route A (embedded Feynman-slash): D = sum_{xy} slash(Delta_xy) f(s^2_xy), Delta_xy = y - x the
embedding separation, slash(Delta) = gamma_mu Delta^mu, f a weight depending only on the invariant interval
s^2 = Delta^2.

  RESULT 1 (positive, decisive): this operator is MANIFESTLY LORENTZ-COVARIANT. Because slash(Lambda Delta)
  = S(Lambda) slash(Delta) S(Lambda)^-1 (the code gamma's realise S^-1 gamma^mu S = Lambda^mu_nu gamma^nu)
  and s^2 is boost-invariant, D -> (I_N (x) S) D (I_N (x) S^-1) under a boost: machine-exact (~1e-15). So the
  code frame RESOLVES the spin-frame/Lorentz concern -- not merely "bypasses" it; it gives a constructively
  Lorentz-covariant Dirac operator. The order-alone no-go is about a BARE causal set; the coded one has the
  frame.

  RESULT 2 (the real open piece): LOCALITY. A Lorentz-invariant weight f(s^2) CANNOT localise the operator,
  because the light cone is s^2=0, so f(s^2)~f(0) along it: arbitrarily-distant near-null pairs couple. The
  max spatial separation of pairs with 0<s^2<ell^2 grows without bound with the box size. This is the
  STANDARD causal-set non-locality (the light-cone "tail"), NOT spinor-specific -- the SCALAR d'Alembertian
  has the identical tail, cured by the Benincasa-Dowker LAYERED alternating sum.

WAYS FORWARD:
  Route A (here): Feynman-slash -- Lorentz-covariant (done); locality open.
  Route B (recommended, tractable): port the BD layered alternating sum to the slash (= Johnston's causal-set
    Dirac) -> inherits Lorentz invariance (A) AND locality (BD). Next computation: build the BD-layered slash,
    verify D^2 -> box + m^2 and the relativistic dispersion E^2 = p^2 + m^2.
  Route C (physical/deep): the §3.5 LATTICE Dirac (walk W = S.C) already works at the lattice scale with the
    vierbein = the lattice geometry (3 colour-bipyramid axes + time); the trans-Lambda lift asks whether that
    vierbein survives the continuum limit as a smooth Lorentz-invariant frame field (ties to the SMG/TCH
    emergent-LI continuum extrapolation). Curved-space vierbein (gravity coupling) is the deepest part.

Self-asserting; exit 0.
"""
import numpy as np
from scipy.linalg import expm


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    I2 = np.eye(2); sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
    g0 = np.kron(sz, I2); g = [g0] + [g0 @ np.kron(sx, s) for s in (sx, sy, sz)]   # g^0=beta, g^i=beta alpha_i
    eta = np.diag([1, -1, -1, -1.])
    slash = lambda D: g[0] * D[0] - g[1] * D[1] - g[2] * D[2] - g[3] * D[3]        # gamma_mu Delta^mu

    print("=== code spin frame + Feynman-slash causal-set Dirac (Route A probe) ===")
    ok(all(np.allclose(g[m] @ g[n] + g[n] @ g[m], 2 * eta[m, n] * np.eye(4)) for m in range(4) for n in range(4)),
       "[1] code gamma matrices Clifford Cl(3,1) exact")

    rng = np.random.default_rng(1); D = rng.normal(size=4)
    ok(np.allclose(slash(D) @ slash(D), (D[0]**2 - D[1]**2 - D[2]**2 - D[3]**2) * np.eye(4)),
       "[2] slash identity (gamma_mu Delta^mu)^2 = s^2 I (the invariant interval)")

    phi = 0.7
    Lam = np.array([[np.cosh(phi), np.sinh(phi), 0, 0], [np.sinh(phi), np.cosh(phi), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1.]])
    S = expm(0.5 * phi * g[0] @ g[1])
    ok(np.allclose(np.linalg.inv(S) @ g[1] @ S, np.cosh(phi) * g[1] + np.sinh(phi) * g[0]) and
       np.allclose(np.linalg.inv(S) @ g[0] @ S, np.cosh(phi) * g[0] + np.sinh(phi) * g[1]),
       "[3] spinor boost S=exp(phi/2 g0 g1): S^-1 gamma^mu S = Lambda^mu_nu gamma^nu")

    worst = 0.0
    for _ in range(300):
        D = rng.normal(size=4); D[0] = abs(D[0]) + np.linalg.norm(D[1:]) + 0.1   # future-timelike
        worst = max(worst, np.max(np.abs(slash(Lam @ D) - S @ slash(D) @ np.linalg.inv(S))))
    ok(worst < 1e-12, f"[4] operator covariance slash(Lambda Delta) = S slash(Delta) S^-1 over 300 pairs (err {worst:.1e}) -> D MANIFESTLY Lorentz-covariant")

    # [5] locality: invariant weight cannot localise (light-cone tail). 2D slice, fixed density, ref at origin.
    print("  [5] LOCALITY (the real open piece): near-null reach vs box size (fixed density), ref at origin")
    reach = {}
    for L in (10, 20, 40):
        N = int(20 * (2 * L) * L)                       # fixed density in t in [0,2L], x in [-L,L]
        t = rng.uniform(0, 2 * L, N); x = rng.uniform(-L, L, N)
        s2 = t * t - x * x; nn = (s2 > 0) & (s2 < 1.0)   # 0 < s^2 < 1 : the invariant 'near' set
        reach[L] = float(np.max(np.abs(x[nn]))) if nn.any() else 0.0
        print(f"      box L={L:3d}: max |x-separation| of pairs with 0<s^2<1 = {reach[L]:5.1f}")
    ok(reach[40] > reach[20] > reach[10] > 1.0,
       "near-null spatial reach GROWS with box size -> operator is NON-LOCAL (the light-cone tail; s^2=0 on the cone)")

    # [6] the bridge to Route B: the BD locality cure (layered sum over causal-interval cardinality) is
    # Lorentz-INVARIANT, so it ports to the covariant slash without spoiling covariance.
    pts = rng.uniform(-6, 6, size=(400, 2)); pts[:, 0] = rng.uniform(0, 12, 400)
    causal = lambda a, b: (b[0] - a[0] > abs(b[1] - a[1]))           # 2D timelike future
    x_, y_ = None, None
    for i in range(400):
        for j in range(400):
            if causal(pts[i], pts[j]) and sum(causal(pts[i], pts[k]) and causal(pts[k], pts[j]) for k in range(400)) >= 3:
                x_, y_ = i, j; break
        if x_ is not None:
            break
    card = lambda P: sum(causal(P[x_], P[k]) and causal(P[k], P[y_]) for k in range(400))
    B = np.array([[np.cosh(0.9), np.sinh(0.9)], [np.sinh(0.9), np.cosh(0.9)]])   # 2D boost
    boosted = pts @ B.T
    print(f"  [6] causal-interval cardinality |{{z: x<z<y}}| = {card(pts)} (rest) = {card(boosted)} (boosted)")
    ok(card(pts) == card(boosted) and card(pts) >= 3,
       "BD layer index (interval cardinality) is Lorentz-INVARIANT -> the locality cure PORTS to the covariant slash (Route B viable)")

    print("\n[verdict] ways forward with the causal-set Dirac construction:")
    print("  - the code frame gives a MANIFESTLY LORENTZ-COVARIANT Dirac operator (Feynman-slash + invariant")
    print("    weight): spin-frame/Lorentz concern RESOLVED (machine-exact covariance). Route A done on covariance.")
    print("  - the remaining open piece is LOCALITY (the light-cone tail) -- the SCALAR sector's already-solved")
    print("    problem, NOT spinor-specific. RECOMMENDED next step (Route B): port the Benincasa-Dowker layered")
    print("    alternating sum to the slash (= Johnston's causal-set Dirac); verify D^2 -> box + m^2 and the")
    print("    relativistic dispersion. Route C (deep): lattice Dirac (§3.5) -> trans-Lambda vierbein continuum")
    print("    limit; curved-space vierbein (gravity) is the deepest part. exit 0")


if __name__ == "__main__":
    main()
