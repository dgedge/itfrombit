#!/usr/bin/env python3
"""
Gauss-projected U(1), SU(2), and SU(3) one-plaquette proxy tests.

Purpose:
  css_tch_operator_candidates.py found an abelian center/Bianchi candidate in a
  finite compact-link proxy. This script tests the same idea in the physical
  Gauss-projected one-plaquette Hilbert space:

    U(1): electric/rotor basis |n>, n in Z.
    SU(2): gauge-invariant character basis chi_j(U).
    SU(3): gauge-invariant character basis chi_(p,q)(U).

  In the non-abelian cases, Gauss projection reduces pure one-plaquette states
  to class functions of the plaquette holonomy. The Wilson term acts by
  multiplication with the fundamental character:

    SU(2): chi_1/2 chi_j = chi_{j+1/2} + chi_{j-1/2}
    SU(3): chi_3 chi_(p,q) =
             chi_(p+1,q) + chi_(p-1,q+1) + chi_(p,q-1)

What it checks:
  1. Center sectors are exactly isospectral after Gauss projection:
       SU(2): z = +/-1
       SU(3): z in Z_3
     This is the non-abelian analogue of "sector crossings are gauge/center
     copies", not independent physical mirror states.

  2. The lowest internal gap in the projected plaquette problem remains finite
     along the Kogut-Susskind weak-coupling trajectory beta = 1/g^2.

Limitations:
  This is still a pure-gauge one-plaquette proxy. It does not construct the TCH
  truncated-cube X-stabilizer, does not include mirror fermions, and does not
  prove the interacting chiral gauge theory. It only says that the simple
  non-abelian twist-closing can be lifted by proper Gauss projection in this
  character-basis proxy.

numpy; self-asserting.
"""

import numpy as np


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def gap(eigvals):
    return float(eigvals[1] - eigvals[0])


def u1_rotor_eigs(n_max, beta, g2):
    ns = np.arange(-n_max, n_max + 1)
    hamiltonian = np.diag((g2 / 2) * ns**2).astype(complex)
    for index in range(2 * n_max):
        hamiltonian[index, index + 1] += -beta / 2
        hamiltonian[index + 1, index] += -beta / 2
    return np.linalg.eigvalsh(hamiltonian)


def su2_character_eigs(n_max, beta, g2, center=1):
    # n = 2j, so j = n/2.
    hamiltonian = np.diag(
        [g2 * (n / 2) * (n / 2 + 1) for n in range(n_max + 1)]
    ).astype(complex)
    for n in range(n_max):
        hamiltonian[n, n + 1] += -beta * center
        hamiltonian[n + 1, n] += -beta * center
    return np.linalg.eigvalsh(hamiltonian)


def su3_reps(level):
    return [(p, q) for p in range(level + 1) for q in range(level + 1 - p)]


def su3_casimir(p, q):
    return (p * p + q * q + p * q + 3 * p + 3 * q) / 3


def su3_character_eigs(level, beta, g2, center=1):
    reps = su3_reps(level)
    index = {rep: i for i, rep in enumerate(reps)}
    dim = len(reps)

    mult_fund = np.zeros((dim, dim), complex)
    for col, (p, q) in enumerate(reps):
        for target in [(p + 1, q), (p - 1, q + 1), (p, q - 1)]:
            if target in index:
                mult_fund[index[target], col] += 1

    electric = np.diag([g2 * su3_casimir(*rep) for rep in reps]).astype(complex)
    magnetic = -(beta / 2) * (
        center * mult_fund + np.conj(center) * mult_fund.conj().T
    )
    hamiltonian = electric + magnetic
    hermiticity_error = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
    assert hermiticity_error < 1e-10
    return np.linalg.eigvalsh(hamiltonian)


def report_u1():
    hd("U(1) Projected Rotor: beta-limit check")
    fixed_gaps = []
    ks_gaps = []
    for beta in [0.5, 1, 2, 5, 10, 20, 50, 100]:
        fixed = gap(u1_rotor_eigs(n_max=80, beta=beta, g2=1.0))
        ks = gap(u1_rotor_eigs(n_max=80, beta=beta, g2=1 / beta))
        fixed_gaps.append(fixed)
        ks_gaps.append(ks)
        print(f"  beta={beta:<5}: fixed-g2 gap={fixed:.6g}; KS g2=1/beta gap={ks:.6g}")
    assert ks_gaps[-1] > 0.9
    print(
        "  -> Along the KS trajectory beta=1/g^2, the projected U(1) rotor gap "
        "approaches a finite value near 1. The previous finite-Z4 decreasing trend "
        "was not the continuum trajectory."
    )


def report_su2():
    hd("SU(2) Gauss-Projected Character Proxy")
    max_center_spread = 0.0
    ks_gaps = []
    for beta in [0.5, 1, 2, 5, 10, 20, 50, 100]:
        g2 = 1 / beta
        spectra = [
            su2_character_eigs(n_max=180, beta=beta, g2=g2, center=center)[:8]
            for center in [+1, -1]
        ]
        center_spread = float(np.max(np.ptp(np.array(spectra), axis=0)))
        max_center_spread = max(max_center_spread, center_spread)
        current_gap = gap(spectra[0])
        ks_gaps.append(current_gap)
        print(
            f"  beta={beta:<5}: KS gap={current_gap:.6g}; "
            f"center-sector spread={center_spread:.3e}"
        )
    assert max_center_spread < 1e-9
    assert min(ks_gaps) > 1.5
    print(
        "  -> SU(2) center sectors are exactly isospectral after Gauss projection, "
        "and the projected character gap stays finite along beta=1/g^2."
    )


def report_su3():
    hd("SU(3) Gauss-Projected Character Proxy")
    centers = [np.exp(2j * np.pi * k / 3) for k in range(3)]
    max_center_spread = 0.0
    ks_gaps = []
    for beta, level in [(0.5, 18), (1, 18), (2, 20), (5, 24), (10, 28), (20, 32), (30, 36)]:
        g2 = 1 / beta
        spectra = [
            su3_character_eigs(level=level, beta=beta, g2=g2, center=center)[:8]
            for center in centers
        ]
        center_spread = float(np.max(np.ptp(np.array(spectra), axis=0)))
        max_center_spread = max(max_center_spread, center_spread)
        current_gap = gap(spectra[0])
        ks_gaps.append(current_gap)
        print(
            f"  beta={beta:<5}, trunc={level:<2}: KS gap={current_gap:.6g}; "
            f"Z3-sector spread={center_spread:.3e}"
        )
    assert max_center_spread < 1e-8
    assert min(ks_gaps) > 1.2
    print(
        "  -> SU(3) center sectors are exactly isospectral after Gauss projection, "
        "and the projected character gap stays finite in this truncated scan. "
        "The high-beta rows use larger irrep cutoffs because the wavefunction spreads "
        "to higher representations."
    )


def main():
    hd("Scope")
    print(
        "This tests non-abelian Gauss projection in the one-plaquette character "
        "basis. It is a proxy for the missing TCH X-plaquette operator, not that "
        "operator itself."
    )
    report_u1()
    report_su2()
    report_su3()
    hd("Verdict")
    print(
        "The abelian and non-abelian one-plaquette proxies all support the same "
        "candidate pattern: raw sector crossings are center/Gauss copies after "
        "projection, and the projected internal gap is finite along the KS "
        "beta=1/g^2 trajectory in these truncations."
    )
    print(
        "\nThis is encouraging but still below the real wall. The actual next object "
        "needed is the TCH truncated-cube X-plaquette/Gauss operator acting on the "
        "CSS mirror interaction, followed by the same gap scan with mirror fermions "
        "and dynamical non-abelian gauge fields."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
