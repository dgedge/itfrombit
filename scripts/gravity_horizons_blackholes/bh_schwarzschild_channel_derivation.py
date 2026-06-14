#!/usr/bin/env python3
r"""Schwarzschild shell channel V_Sch(M) from the finite V_cell map.

The finite stabilizer map V_cell is already constructed in
bh_isometry_v_construction.py.  The remaining black-hole question is whether
there is a canonical radial composition over an actual Schwarzschild horizon:

  * which cells enter the map,
  * what scheduler/rates act on the invalid Q subspace,
  * whether the KMS line weights follow from the localized-mass steady state.

This script closes the algebraic/radial part to the honest boundary.

Result:

  V_Sch(M) is the direct-sum isometry over the radial bonds/cells crossing the
  frozen-coin horizon shell:

      |x>_shell |[s]>_B |gamma>_vac
          -> |x>_R |delta(s)>_syndrome |gamma>_latch ,

  where x labels a horizon-shell cell, [s] is the Q3 state modulo global
  complement, and gamma is the vacuum/complement latch.  Orthogonal shell labels
  make the direct sum isometric.

  The radial scheduler is the local one-bit jump generator on Q with rates

      W_{c->c'}(x) = Gamma_x exp[- beta_H Delta E_infty(c,c') / 2]

  for one-bit neighbors c,c' in Q.  Equivalently, using local cell energy
  E_loc(c)=epsilon_F F(c) and Tolman temperature T_loc=T_H/sqrt(f), the exponent
  is beta_loc Delta E_loc / 2.  The Tolman redshift cancels:

      beta_loc E_loc = beta_infty E_infty.

  Therefore the stationary line weights are exactly

      P(F) proportional to g_Q(F) exp[- beta_eff F],

  with beta_eff = epsilon_F / T_H = 1/(2 phi T_H) in the canon's item-10 units.

Honest boundary:

  This derives KMS from the localized-mass horizon steady state only under the
  standard Schwarzschild regularity/Davies-scheduler premise: the shell dynamics
  is stationary with respect to the horizon Killing time and reversible against
  the local Tolman bath.  The finite V_cell algebra alone cannot force that
  premise.  Also still open: the exact dispersive freeze surfaces r_F(F;M), i.e.
  the item-9 rho_c(F) relation and absolute flux normalisation.
"""

from __future__ import annotations

from collections import Counter
import math

import numpy as np


PHI = (math.sqrt(5.0) - 1.0) / 2.0
LAMBDA_QCD_GEV = 0.332
M_PLANCK_GEV = 1.220890e19
GEV_PER_KG = 5.60959e26
M_SUN_GEV = 1.98892e30 * GEV_PER_KG
ALL = (1 << 8) - 1
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    return (
        not (bit(n, 0) and bit(n, 1))
        and bit(n, 7) == bit(n, 6)
        and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0)))
    )


def strain(n: int) -> int:
    return sum(1 for i, j in EDGES if bit(n, i) != bit(n, j))


def syndrome_word(s: int) -> int:
    out = 0
    for n, (i, j) in enumerate(EDGES):
        out |= (bit(s, i) ^ bit(s, j)) << n
    return out


def complement_representatives() -> list[int]:
    return [s for s in range(256) if s < (s ^ ALL)]


def invalid_subspace() -> list[int]:
    return [s for s in range(256) if not valid(s)]


def record_row(shell_cell: int, syn: int, latch: int) -> tuple[int, int, int]:
    return (shell_cell, syn, latch)


def verify_direct_sum_v_sch(n_shell_test: int = 5) -> None:
    """A finite proxy for the direct-sum shell isometry.

    The real shell has N_H(M) enormous.  Injectivity is independent of N_H:
    rows are (cell label, syndrome, latch), and the cell label is orthogonal.
    """

    reps = complement_representatives()
    rows = []
    for x in range(n_shell_test):
        for rep in reps:
            for gamma in (0, 1):
                s = rep ^ (ALL if gamma else 0)
                rows.append(record_row(x, syndrome_word(s), gamma))
    counts = Counter(rows)
    collisions = [row for row, count in counts.items() if count > 1]
    assert not collisions
    print(
        f"    V_Sch direct-sum injectivity: {n_shell_test} shell cells x "
        f"256 columns -> {len(rows)} unique rows"
    )


def horizon_radius(M_gev: float) -> float:
    # Natural units: G = 1/M_P^2, r_s = 2GM.
    return 2.0 * M_gev / M_PLANCK_GEV**2


def surface_gravity(M_gev: float) -> float:
    # kappa = 1/(2 r_s) = 1/(4GM).
    return 1.0 / (2.0 * horizon_radius(M_gev))


def hawking_temperature_infinity(M_gev: float) -> float:
    return surface_gravity(M_gev) / (2.0 * math.pi)


def schwarzschild_f(r: float, M_gev: float) -> float:
    return 1.0 - horizon_radius(M_gev) / r


def shell_cell_count(M_gev: float) -> float:
    # Canon node area A_node=1/(4 Lambda^2), so N=A/A_node=16 pi r_s^2 Lambda^2.
    rs = horizon_radius(M_gev)
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV * LAMBDA_QCD_GEV


def verify_tolman_cancellation() -> None:
    eps_f = 1.0 / (2.0 * PHI)
    for m_solar in (3.0, 30.0, 4.3e6):
        M = m_solar * M_SUN_GEV
        beta_inf = 1.0 / hawking_temperature_infinity(M)
        rs = horizon_radius(M)
        for eps in (1.0e-6, 1.0e-4, 1.0e-2):
            r = rs * (1.0 + eps)
            f = schwarzschild_f(r, M)
            beta_loc = beta_inf * math.sqrt(f)
            for F in (3, 4, 6, 9, 12):
                e_loc = eps_f * F
                e_inf = math.sqrt(f) * e_loc
                assert abs(beta_loc * e_loc / (beta_inf * e_inf) - 1.0) < 1e-13
    print("    Tolman/Killing exponent cancellation: PASS")


def local_kms_generator(states: list[int], beta_eff: float, gamma: float = 1.0) -> np.ndarray:
    index = {s: i for i, s in enumerate(states)}
    fvals = np.array([strain(s) for s in states], dtype=float)
    gen = np.zeros((len(states), len(states)), dtype=float)
    for i, s in enumerate(states):
        for k in range(8):
            t = s ^ (1 << k)
            j = index.get(t)
            if j is not None:
                gen[i, j] = gamma * math.exp(-0.5 * beta_eff * (fvals[j] - fvals[i]))
    gen[np.diag_indices_from(gen)] = -gen.sum(axis=1)
    return gen


def verify_kms_stationary(beta_eff: float) -> dict[int, float]:
    states = invalid_subspace()
    gq = dict(sorted(Counter(strain(s) for s in states).items()))
    assert gq == TARGET_GQ
    gen = local_kms_generator(states, beta_eff)

    fvals = np.array([strain(s) for s in states], dtype=float)
    pi = np.exp(-beta_eff * fvals)
    pi /= pi.sum()
    detailed_balance = np.max(np.abs(pi[:, None] * gen - pi[None, :] * gen.T))
    assert detailed_balance < 2e-14
    stationary_residual = np.max(np.abs(pi @ gen))
    assert stationary_residual < 2e-13

    line = {F: 0.0 for F in TARGET_GQ}
    for s, p in zip(states, pi):
        line[strain(s)] += float(p)
    raw = {F: g * math.exp(-beta_eff * F) for F, g in TARGET_GQ.items()}
    z = sum(raw.values())
    target = {F: v / z for F, v in raw.items()}
    assert max(abs(line[F] - target[F]) for F in TARGET_GQ) < 1e-14
    return dict(sorted(line.items()))


def main() -> None:
    print("[0] V_Sch(M) finite-cell direct sum")
    verify_direct_sum_v_sch()
    for m_solar in (3.0, 30.0, 4.3e6):
        n_shell = shell_cell_count(m_solar * M_SUN_GEV)
        temp = hawking_temperature_infinity(m_solar * M_SUN_GEV)
        print(
            f"    M={m_solar:g} Msun: N_H=16 pi r_s^2 Lambda^2 = {n_shell:.3e}, "
            f"T_H(infinity)={temp:.3e} GeV"
        )

    print("\n[1] Which cells enter")
    print("    The channel acts on radial bonds/cells straddling the frozen-coin")
    print("    surface: exterior endpoint still has a functioning coin, interior")
    print("    endpoint is in the obligatory-inward/frozen-coin phase.")
    print("    To leading thin-shell order this is the horizon shell with N_H cells.")
    print("    The exact dispersive surfaces r_F(F;M) require the still-open item-9")
    print("    rho_c(F) / radial-frustration-gradient relation.")

    print("\n[2] Schwarzschild/Tolman rate exponent")
    verify_tolman_cancellation()
    print("    A local jump with E_loc=F/(2 phi) read at T_loc=T_H/sqrt(f)")
    print("    has the same exponent as the redshifted energy E_inf=sqrt(f)E_loc")
    print("    read at the asymptotic Hawking temperature T_H.")

    print("\n[3] KMS scheduler on the 208-state Q subspace")
    for beta_eff in (0.25, 0.70, 1.00, 1.60):
        line = verify_kms_stationary(beta_eff)
        ratio = line[4] / line[3]
        target_ratio = (TARGET_GQ[4] / TARGET_GQ[3]) * math.exp(-beta_eff)
        print(
            f"    beta_eff={beta_eff:4.2f}: I4/I3={ratio:.8f}; "
            f"target={(target_ratio):.8f}; line P(F)={line}"
        )

    print(
        "\n[4] Verdict\n"
        "    The Schwarzschild shell channel is fixed in form:\n"
        "      V_Sch(M) = direct-sum_x V_cell,x over frozen-coin shell cells,\n"
        "      with x carrying the orthogonal horizon-cell address.\n"
        "    The scheduler is the local one-bit reversible/Davies generator on Q,\n"
        "      W_{c->c'} = Gamma_x exp[-beta_eff (F(c')-F(c))/2].\n"
        "    Under Schwarzschild regularity/Tolman local equilibrium, beta_eff is\n"
        "      the item-10 beta = 1/(2 phi T_H), and KMS follows: the stationary\n"
        "      radiation lines are exactly g_Q(F) exp[-beta_eff F].\n"
        "    Not closed by this algebra: the exact r_F(F;M) shell dispersion and\n"
        "      absolute flux normalisation.  Those remain item-9 / freeze-surface\n"
        "      data, not V_cell data.\n"
        "ALL ASSERTIONS PASSED"
    )


if __name__ == "__main__":
    main()
