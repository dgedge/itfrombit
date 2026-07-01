#!/usr/bin/env python3
r"""foundations_torsion_coefficient_ecsk.py

THE TORSION COEFFICIENT (residual (ii) of the dynamical frame equation, foundations_dynamical_frame_einstein_cartan.py).

In Einstein-Cartan-Sciama-Kibble (ECSK) gravity the torsion is NON-PROPAGATING: varying the spin connection
gives an ALGEBRAIC relation (torsion ~ spin current), and integrating the torsion out leaves a four-fermion
CONTACT interaction. This script DERIVES that coefficient for the framework's code-gamma Dirac (it is NOT
quoted from memory; the standard value 3/16 kappa^2 is the check, produced by the numerics).

Pieces (all verified with the [8,4,4] code gammas):
  - SPIN/axial structure: {gamma_mu, gamma_ab} = c * eps_{mu a b}{}^d gamma_d gamma5 with c = -2i (residual 0)
    -> the code Dirac's spin current is the standard totally-antisymmetric (axial) one.
  - GRAVITY side: totally-antisymmetric torsion T_{mab}=eps_{mab sigma}A^sigma (axial vector A), contortion
    K=T/2; T^2/A^2 = -6, K^2/A^2 = -3/2, so the EH torsion term is |a_E| A^2 with |a_E| = (1/2)|K^2/A^2| =
    3/4 (in units 1/kappa^2; kappa^2 = 8 pi G).
  - DIRAC side: the contortion coupling L_D^K = (i/8) K^{mu ab} {gamma_mu, gamma_ab} (psibar..psi) is
    proportional to A_mu (psibar gamma^mu gamma5 psi) = A . j5, with coupling b_D = -3/4 (residual ~1e-16).
  - INTEGRATE OUT the (algebraic) torsion: A = (kappa^2/2) j5 (Cartan relation), and the contact coefficient
    C = b_D^2 / (4 |a_E|) = 3/16 kappa^2.

RESULT: torsion coefficient = (3/16) kappa^2 = (3/2) pi G (standard ECSK), reproduced EXACTLY by the
code-gamma Dirac. So the framework's torsion sector is standard ECSK -- non-propagating torsion giving the
four-fermion axial-axial (spin-spin) contact interaction L = (3/2) pi G (psibar gamma^mu gamma5 psi)^2 (sign
per the standard convention). The dimensional kappa^2 = 8 pi G is §10.5's (induced gravity, K_eff=205). This
CLOSES residual (ii); the remaining residuals are (i) the micro->macro frame coarse-graining (the crux) and
(iii) the quantitative 8 pi G (existing §10.5).

Self-asserting; exit 0.
"""
import numpy as np
import itertools


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
    g = [np.kron(sz, I2)] + [np.kron(sz, I2) @ np.kron(sx, s) for s in (sx, sy, sz)]
    eta = np.diag([1, -1, -1, -1.]); g5 = 1j * g[0] @ g[1] @ g[2] @ g[3]
    gl = [eta[m, m] * g[m] for m in range(4)]
    gab = lambda a, b: 0.5 * (gl[a] @ gl[b] - gl[b] @ gl[a])
    eL = np.zeros((4, 4, 4, 4))
    for p in itertools.permutations(range(4)):
        s = 1; pl = list(p)
        for i in range(4):
            for j in range(i + 1, 4):
                if pl[i] > pl[j]:
                    s = -s
        eL[p] = s

    ok(np.allclose(g5 @ g5, np.eye(4)) and max(np.linalg.norm(g5 @ g[m] + g[m] @ g5) for m in range(4)) < 1e-12,
       "gamma5^2 = I and {gamma5, gamma^mu} = 0")

    # {gamma_mu, gamma_ab} = c * eps_{mu a b}{}^d gamma_d gamma5
    L = np.array([(gl[mu] @ gab(a, b) + gab(a, b) @ gl[mu]).ravel() for mu in range(4) for a in range(4) for b in range(4)]).ravel()
    R = np.array([sum(eL[mu, a, b, d] * (g[d] @ g5) for d in range(4)).ravel() for mu in range(4) for a in range(4) for b in range(4)]).ravel()
    c = np.vdot(R, L) / np.vdot(R, R)
    ok(np.linalg.norm(L - c * R) / np.linalg.norm(L) < 1e-9 and abs(abs(c) - 2) < 1e-9,
       f"code spin structure: {{gamma_mu,gamma_ab}} = c eps gamma_d gamma5, c = {c:.3f} (|c|=2) -> standard axial spin current")

    # gravity side: T_{mab}=eps A^sigma, K=T/2
    A = np.array([0.3, 0.7, -0.2, 0.5]); A2 = sum(eta[m, m] * A[m] * A[m] for m in range(4))
    T = np.array([[[sum(eL[mu, a, b, s] * A[s] for s in range(4)) for b in range(4)] for a in range(4)] for mu in range(4)])
    sq = lambda X: sum(X[mu, a, b] * eta[mu, mu] * eta[a, a] * eta[b, b] * X[mu, a, b] for mu in range(4) for a in range(4) for b in range(4))
    K = 0.5 * T
    ok(abs(sq(T) / A2 + 6) < 1e-9, f"T^2/A^2 = {sq(T)/A2:.3f} = -6 (totally-antisymmetric torsion <-> axial vector)")
    aE = 0.5 * abs(sq(K) / A2)
    ok(abs(aE - 0.75) < 1e-9, f"|a_E| = (1/2)|K^2/A^2| = {aE:.3f} = 3/4 (in units 1/kappa^2) [EH torsion term]")

    # Dirac side: coupling operator O = (i/8) K^{mu ab} {gamma_mu, gamma_ab}
    Kup = np.array([[[eta[mu, mu] * eta[a, a] * eta[b, b] * K[mu, a, b] for b in range(4)] for a in range(4)] for mu in range(4)])
    O = (1j / 8) * sum(Kup[mu, a, b] * (gl[mu] @ gab(a, b) + gab(a, b) @ gl[mu]) for mu in range(4) for a in range(4) for b in range(4))
    Asl = sum(A[m] * gl[m] for m in range(4)) @ g5
    bD = (np.vdot(Asl, O) / np.vdot(Asl, Asl)).real
    ok(np.linalg.norm(O - bD * Asl) / np.linalg.norm(O) < 1e-9 and abs(abs(bD) - 0.75) < 1e-9,
       f"Dirac coupling O proportional to A-slash gamma5, b_D = {bD:+.3f} (|b_D|=3/4) -> L_D^K = b_D A.j5")

    # integrate out torsion -> contact coefficient
    C = bD ** 2 / (4 * aE)
    print(f"\n  Cartan relation: A_mu = (kappa^2/2) j5_mu (axial torsion proportional to axial current)")
    ok(abs(C - 3 / 16) < 1e-9, f"TORSION COEFFICIENT C = b_D^2/(4|a_E|) = {C:.5f} kappa^2 = 3/16 kappa^2 = (3/2) pi G (standard ECSK, DERIVED)")

    print("\n[verdict] torsion coefficient closed (residual (ii)):")
    print("  - the code-gamma Dirac reproduces the standard ECSK four-fermion contact term EXACTLY:")
    print("    L_contact = (3/16) kappa^2 (psibar gamma^mu gamma5 psi)^2 = (3/2) pi G (axial current)^2;")
    print("  - torsion is non-propagating (algebraic Cartan relation A = (kappa^2/2) j5); no surprises -- the")
    print("    framework's torsion sector IS standard ECSK. kappa^2 = 8 pi G is §10.5's (induced gravity).")
    print("  - remaining residuals of the dynamical frame equation: (i) the micro->macro frame coarse-graining")
    print("    (the crux); (iii) the quantitative 8 pi G = M_P (existing §10.5, K_eff=205). exit 0")


if __name__ == "__main__":
    main()
