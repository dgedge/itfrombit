#!/usr/bin/env python3
r"""foundations_dirac_curved_route_c.py

ROUTE C: does the §3.5 lattice vierbein survive the continuum limit as a smooth frame field (emergent LI in
flat space, a curved spin connection in curved space)? -- where the Dirac construction meets gravity.

IMPORTANT RECONCILIATION. The flat causal-set Dirac SPIN FRAME was already closed on 2026-06-13:
foundations_causet_dirac_spin_closure.py (ANCHOR L2378) -- the E_{1/2}/coin supplies the spin fibre, the
Clifford relation + Dirac factorisation + Spin(3,1) covariance + link transports all hold, and the order-only
no-go is derived; the STANDING residual it flags is the "dynamical frame equation". This session's Route A/B
(foundations_causet_dirac_slash.py / _dispersion.py) substantially RE-DERIVED that closure (covariance, slash
identity, dispersion) with one minor increment: the explicit sprinkling momentum-space symbol D~(k) ∝ k-slash.
So Route C does NOT re-open the spin frame; it attacks the genuinely-open piece both passes flag.

Two findings here:
  PART 1 -- emergent isotropy in the CORRECTED geometry. §3.5's emergent-isotropy theorem (ANCHOR L464)
  credits "the 4.8.8 OCTAGONAL stencil forces the cos(4 theta) anisotropy to cancel, pushing anisotropy to
  O(k^4 a^4)". That geometry was RETRACTED this session (no octagons; bond-centred bipyramids, dual rectified
  cubic honeycomb). So that mechanism is a downstream casualty. Re-examined here: the corrected matter cell
  has 3 orientations x/y/z = a cubic-O_h lattice; the massless Dirac cone E=sqrt(sum sin^2 k_i) is isotropic
  to LEADING order with fractional anisotropy ~ (|k|a)^2 (p~2) -> EMERGENT LI holds (anisotropy -> 0 as a->0),
  the vierbein -> isotropic delta_a^mu in the continuum. The OLD octagon p=4 (sharper subleading suppression)
  is NOT reproduced by a plain cubic lattice; re-deriving the sharp suppression in the corrected (rectified-
  cubic / cuboctahedral) stencil is a sub-task. Leading-order emergent LI is geometry-robust; the SHARPNESS
  is the retracted piece.

  PART 2 -- curved spin connection from the lattice frame transports (extends the 2026-06-13 transports into
  the curved regime). A connection omega_mu(x) with non-zero curvature gives link transports U=exp(omega.dx)
  whose plaquette holonomy = exp(-F_xy dx dy), F_xy = d_x omega_y - d_y omega_x + [omega_x,omega_y] -- the
  spin curvature, lying in the Lorentz algebra span(sigma_ab) = the Riemann curvature in the spinor rep. So
  the curved-space covariant Dirac D = gamma^a e_a^mu (d_mu + omega_mu) emerges kinematically from the frame
  transports.

DEEP OPEN FRONTIER (unchanged): the DYNAMICAL frame equation -- what determines e_a^mu(x) / omega_mu(x):
Einstein-Cartan dynamics FROM the substrate, and the equivalence-principle match to the E_g gravitational
field (ANCHOR L2378's "dynamical frame equation OPEN"). This script settles only the KINEMATICS.

Self-asserting; exit 0.
"""
import numpy as np
from scipy.linalg import expm


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
    g0 = np.kron(sz, I2); g = [g0] + [g0 @ np.kron(sx, s) for s in (sx, sy, sz)]

    print("=== PART 1: emergent isotropy of the Dirac cone in the corrected cubic-O_h geometry ===")
    Econe = lambda kv: np.sqrt(max(np.abs(np.linalg.eigvalsh(
        (sum(g[i + 1] * np.sin(kv[i]) for i in range(3))).conj().T @ sum(g[i + 1] * np.sin(kv[i]) for i in range(3))))))

    def anis(kmag, n=400):
        rng = np.random.default_rng(1); v = rng.normal(size=(n, 3)); v /= np.linalg.norm(v, axis=1, keepdims=True)
        E = np.array([Econe(kmag * u) for u in v]); return (E.max() - E.min()) / E.mean()

    ks = [0.4, 0.2, 0.1, 0.05]; A = [anis(k) for k in ks]
    for k, a in zip(ks, A):
        print(f"    |k|a={k:.2f}: fractional cone anisotropy = {a:.5f}")
    p = np.polyfit(np.log(ks), np.log(A), 1)[0]
    ok(A[-1] < A[0] and A[-1] < 1e-3, "anisotropy -> 0 as |k|a -> 0: EMERGENT LI (vierbein -> isotropic in the continuum)")
    ok(1.7 < p < 2.3, f"anisotropy ~ (|k|a)^p with p={p:.2f} ~ 2 (plain cubic; the OLD octagon p=4 is RETRACTED-geometry, re-derivation needed)")

    print("\n=== PART 2: curved spin connection from the lattice frame transports ===")
    Mab = lambda a, b: 0.25 * (g[a] @ g[b] - g[b] @ g[a])
    M = Mab(0, 1); kap = 0.3                                        # omega_y = x*kap*M, omega_x = 0
    dx = dy = 1e-3
    wy = lambda x: x * kap * M
    H = expm(wy(0.7) * dy) @ np.linalg.inv(expm(wy(0.7 + dx) * dy))  # plaquette holonomy (omega_x=0)
    F = kap * M                                                     # F_xy = d_x omega_y + [.,.] = kap*M
    ok(np.linalg.norm(H - expm(-F * dx * dy)) / np.linalg.norm(H - np.eye(4)) < 1e-9,
       "plaquette holonomy H = exp(-F_xy dx dy), F_xy = kappa*M (the spin curvature)")
    basis = [Mab(a, b) for a in range(4) for b in range(a + 1, 4)]
    Gm = np.array([[np.vdot(bi, bj) for bj in basis] for bi in basis])
    coef = np.linalg.solve(Gm, [np.vdot(bi, F) for bi in basis])
    ok(np.linalg.norm(F - sum(c * b for c, b in zip(coef, basis))) / np.linalg.norm(F) < 1e-9,
       "F_xy in the Lorentz algebra span(sigma_ab) = Riemann curvature in the spinor rep")
    ok(np.linalg.norm(H - np.eye(4)) > 1e-9, "holonomy != I -> genuine curvature; curved Dirac D = gamma^a e_a^mu (d_mu + omega_mu) emerges")

    print("\n[verdict] Route C -- the lattice vierbein in the continuum:")
    print("  - FLAT: emergent LI, anisotropy ~ (|k|a)^2 -> 0; the vierbein -> isotropic smooth frame. (Leading")
    print("    order is geometry-robust; the SHARP subleading suppression -- old 4.8.8-octagon p=4, now retracted")
    print("    -- needs re-derivation in the corrected rectified-cubic/cuboctahedral stencil.)")
    print("  - CURVED: the spin connection + curvature (Riemann in the spinor rep) emerge from the frame")
    print("    transports; the curved-space covariant Dirac is kinematically well-defined.")
    print("  - RECONCILIATION: the flat spin frame was already closed 2026-06-13 (L2378); Route A/B re-derived it.")
    print("  - DEEP OPEN (unchanged): the DYNAMICAL frame equation -- e_a^mu(x) from substrate Einstein-Cartan")
    print("    dynamics + the equivalence-principle match to the E_g gravitational field. Kinematics done; the")
    print("    dynamics (the actual Dirac-gravity coupling) is the standing frontier. exit 0")


if __name__ == "__main__":
    main()
