#!/usr/bin/env python3
r"""foundations_dynamical_frame_einstein_cartan.py

ATTACK on the dynamical frame equation (the standing residual from ANCHOR L2378 / Route C): what determines
e_a^mu(x)? Claim: it is the EINSTEIN-CARTAN equation, and the framework already supplies the pieces -- this
script ASSEMBLES them, verifies the assembly is consistent, and PINS the genuine residual. It does NOT
derive Einstein's equation from scratch (that is §10.5's induced-gravity claim, with its own residuals).

The assembly:
  - FRAME = the substrate cell-strain (the Dirac vierbein, Route C). Its irrep decomposition is the exact
    Einstein-Cartan field content: 16 = 1 (trace/dilaton) + 9 (traceless-symmetric = spin-2 GRAVITON, the
    E_g of §10.1) + 6 (antisymmetric = local Lorentz / spin connection). So the fermion's frame CONTAINS the
    gravitating strain -- the equivalence-principle identification (the Dirac sees the E_g graviton).
  - DYNAMICS = the Sakharov/holographic induced Einstein-Hilbert action (§10.5). KEY: in induced gravity G is
    induced by ALL matter loops and ALL matter couples to the SAME induced metric, so the equivalence-
    principle UNIVERSALITY (same G for every species) is AUTOMATIC -- not a separate postulate.
  - SOURCE = the Dirac matter. Verified here: the Hilbert stress-energy T^munu = (1/2)(k^nu j^mu + k^mu j^nu)
    is SYMMETRIC and CONFORMAL (trace = m * ubar u -> 0 for m=0) -> a proper graviton/Einstein-RHS source;
    the spin current S^{mu ab} ~ psibar {gamma^mu, sigma^ab} psi is TOTALLY ANTISYMMETRIC -> the Cartan
    torsion source. So delta S / delta e -> G_munu = 8 pi G T_munu ; delta S / delta omega -> torsion ~ spin.

RESIDUAL (the genuine open part, honestly pinned):
  (i) the MICRO->MACRO coarse-graining -- the microscopic cell-vierbein (Route C) must coarse-grain to the
      §10.5 macroscopic induced metric (ties to the SMG/TCH continuum extrapolation). THIS is the crux and is
      not touched computationally here.
  (ii) the torsion coefficient / its (non-)propagation;  (iii) the quantitative 8 pi G = M_P (existing §10.5,
      K_eff = 205, with its own residuals).
So this REDUCES "dynamical frame equation OPEN" to a structured Einstein-Cartan problem with the crux pinned
to the micro->macro frame coarse-graining; it is a structural reduction + consistency check, NOT a closure.

Self-asserting; exit 0.
"""
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
    g0 = np.kron(sz, I2); g = [g0] + [g0 @ np.kron(sx, s) for s in (sx, sy, sz)]; eta = np.diag([1, -1, -1, -1.])

    print("=== (1) frame field = trace(1) + graviton(9) + local-Lorentz(6): the EC field content ===")
    H = np.random.default_rng(0).normal(size=(4, 4))
    S = (H + H.T) / 2; A = (H - H.T) / 2; tr = np.trace(S) / 4 * np.eye(4); STL = S - tr
    ok(np.allclose(tr + STL + A, H), "decomposition 16 = 1 + 9 + 6 reconstructs the frame perturbation")
    ok(abs(np.vdot(STL, tr)) + abs(np.vdot(STL, A)) + abs(np.vdot(tr, A)) < 1e-12,
       "trace / traceless-symmetric (=E_g graviton) / antisymmetric parts mutually orthogonal")

    print("\n=== (2) Dirac stress-energy sources the GRAVITON (symmetric + conformal) ===")
    m, p = 0.6, 0.8; E = np.sqrt(p * p + m * m); k = np.array([E, p, 0, 0])
    slk = g[0] * k[0] - g[1] * k[1]; w, V = np.linalg.eig(slk); u = V[:, np.argmin(abs(w - m))]
    ub = u.conj() @ g[0]; j = np.array([(ub @ g[mu] @ u) for mu in range(4)]).real
    T = np.array([[0.5 * (k[nu] * j[mu] + k[mu] * j[nu]) for nu in range(4)] for mu in range(4)])
    trT = sum(eta[a, a] * T[a, a] for a in range(4))
    ok(np.allclose(T, T.T), "T^munu symmetric -> couples to the symmetric metric / graviton")
    ok(abs(trT - (m * (ub @ u)).real) < 1e-9, f"trace T^mu_mu = m * ubar u = {trT:.4f} (Hilbert stress-energy)")
    m0 = 0.0; k0 = np.array([p, p, 0, 0]); slk0 = g[0] * k0[0] - g[1] * k0[1]; w0, V0 = np.linalg.eig(slk0)
    u0 = V0[:, np.argmin(abs(w0 - m0))]; ub0 = u0.conj() @ g[0]; j0 = np.array([(ub0 @ g[mu] @ u0) for mu in range(4)]).real
    T0 = np.array([[0.5 * (k0[nu] * j0[mu] + k0[mu] * j0[nu]) for nu in range(4)] for mu in range(4)])
    ok(abs(sum(eta[a, a] * T0[a, a] for a in range(4))) < 1e-9, "massless: trace T = 0 (CONFORMAL) -> couples to the traceless graviton, not the dilaton")

    print("\n=== (3) Dirac spin current sources TORSION (totally antisymmetric) ===")
    sig = lambda a, b: 0.5j * (g[a] @ g[b] - g[b] @ g[a]); N = lambda M: np.linalg.norm(M)
    swap = max(N(g[mu] @ sig(a, b) + sig(a, b) @ g[mu] + g[a] @ sig(mu, b) + sig(mu, b) @ g[a])
               for mu in range(4) for a in range(4) for b in range(4))
    diag0 = max(N(g[mu] @ sig(mu, b) + sig(mu, b) @ g[mu]) for mu in range(4) for b in range(4))
    ok(swap < 1e-9 and diag0 < 1e-9,
       "{gamma^mu, sigma^ab} totally antisymmetric -> spin current = Cartan torsion source (antisym under mu<->a; vanishes mu=a)")

    print("\n[verdict] the dynamical frame equation = EINSTEIN-CARTAN (structural reduction; residual pinned):")
    print("  - frame = substrate strain; spin-2 part = E_g graviton (§10.1) -> equivalence-principle identification;")
    print("  - dynamics = §10.5 Sakharov/holographic induced EH; G induced by ALL matter loops, ALL matter couples")
    print("    to the SAME metric -> equivalence-principle UNIVERSALITY automatic; ")
    print("  - Dirac is a proper source: symmetric+conformal T (graviton/Einstein RHS) + totally-antisymmetric spin")
    print("    current (torsion); delta/delta e -> G=8piG T, delta/delta omega -> torsion ~ spin.")
    print("  RESIDUAL (genuine, not closed): (i) the MICRO->MACRO coarse-graining of the cell-vierbein to the §10.5")
    print("  induced metric (the crux; ties to SMG/TCH); (ii) the torsion coefficient; (iii) quantitative 8piG=M_P")
    print("  (existing §10.5, K_eff=205). So 'dynamical frame equation OPEN' is REDUCED to (i), not closed. exit 0")


if __name__ == "__main__":
    main()
