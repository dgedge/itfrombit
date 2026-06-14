#!/usr/bin/env python3
"""
Matter-coupled Gauss-projected plaquette proxies for the CSS SMG gap.

Why this exists:
  css_nonabelian_gauss_proxy.py was deliberately pure gauge. It showed that
  one-plaquette U(1), SU(2), and SU(3) character-basis proxies have finite
  internal gauge gaps after Gauss projection, but it did not test whether the
  CSS mirror/matter stabilizer sectors remain separated.

This script adds the next minimal layer:

  * Keep the exact [[8,0,4]] CSS stabilizer-sector reduction. The eight
    stabilizers commute, so each matter sector is labeled by

        x_i = +/-1  for the four X stabilizers,
        z_i = +/-1  for the four Z stabilizers.

  * Couple the X-sector signs to a normalized Wilson/class-function operator in
    a Gauss-projected one-plaquette gauge basis:

        U(1):  cos(theta)
        SU(2): chi_{1/2}(U) / 2
        SU(3): Re chi_3(U) / 3

    In an x-sector the gauge Hamiltonian is

        H_gauge(x) = H_E - (beta + kappa * sum_i x_i) W.

    This is the one-plaquette character-basis shadow of replacing bare X_i by
    gauge-dressed X_i W_i. It is not the missing open-link TCH construction.

What is measured:
  full/gauge gap:
    The first excitation in the fixed CSS ground sector. This can be just a
    gauge oscillator and must not be confused with the mirror gap.

  raw X-matter gap:
    Lowest sector with any changed X stabilizer. This includes sector crossings
    that may be Gauss/center copies.

  closed-flux X gap:
    Same, but restricted to prod_i x_i = +1, the Z2 Bianchi/closed-flux
    candidate selected by css_tch_operator_candidates.py.

  closed-flux + center-quotient X gap:
    For Z2-center shadows (U(1) even-link and SU(2)), additionally identify the
    global sign flip x_i -> -x_i as the center partner of the all-plus sector.
    This is the direct matter-coupled version of "remove sector crossings or
    show they are pure gauge/Gauss-equivalent."

Important limitation:
  The SU(3) center is Z3, while the CSS stabilizer eigenvalues here are Pauli
  signs (Z2). The SU(3) rows are therefore a positive-source Wilson-character
  shadow, not a real SU(3) center quotient. A real TCH/SU(3) build must explain
  how the Pauli CSS X-half embeds into, or is replaced by, the Z3 center data.
  Numerically this is not a wholesale SU(3) collapse: the proxy recovers the
  CSS matter gap in the weak-coupling direction. The failure is concentrated at
  small beta / large g^2, i.e. the confining color regime. This is the executable
  shadow of ANCHOR §2.8's open "F2-closed colour algebra" to one-hot projector
  dual mapping.

numpy; self-asserting.
"""

import itertools

import numpy as np


GROUND_X = (1, 1, 1, 1)
GLOBAL_FLIP_X = (-1, -1, -1, -1)
BETA_VALUES = [0.5, 1, 2, 5, 10, 20, 50, 100]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def parity(xs):
    out = 1
    for value in xs:
        out *= value
    return out


def x_sectors(closed_flux=None):
    sectors = list(itertools.product([1, -1], repeat=4))
    if closed_flux is not None:
        sectors = [xs for xs in sectors if parity(xs) == closed_flux]
    return sectors


def su3_reps(level):
    return [(p, q) for p in range(level + 1) for q in range(level + 1 - p)]


def su3_casimir(p, q):
    return (p * p + q * q + p * q + 3 * p + 3 * q) / 3


class U1Rotor:
    name = "U(1) rotor"
    center_quotient = True

    def __init__(self, n_max=120):
        self.n_max = n_max
        ns = np.arange(-n_max, n_max + 1)
        self.electric_diag = ns**2 / 2
        dim = len(ns)
        self.wilson = np.zeros((dim, dim), complex)
        for index in range(dim - 1):
            self.wilson[index, index + 1] = 0.5
            self.wilson[index + 1, index] = 0.5

    def eigs(self, rho, g2):
        hamiltonian = np.diag(g2 * self.electric_diag).astype(complex)
        hamiltonian -= rho * self.wilson
        return np.linalg.eigvalsh(hamiltonian)[:2]


class SU2Character:
    name = "SU(2) character"
    center_quotient = True

    def __init__(self, n_max=240):
        self.n_max = n_max
        self.electric_diag = np.array(
            [(n / 2) * (n / 2 + 1) for n in range(n_max + 1)]
        )
        self.wilson = np.zeros((n_max + 1, n_max + 1), complex)
        for n in range(n_max):
            # Multiplication by chi_{1/2}/2.
            self.wilson[n, n + 1] = 0.5
            self.wilson[n + 1, n] = 0.5

    def eigs(self, rho, g2):
        hamiltonian = np.diag(g2 * self.electric_diag).astype(complex)
        hamiltonian -= rho * self.wilson
        return np.linalg.eigvalsh(hamiltonian)[:2]


class SU3Character:
    name = "SU(3) character shadow"
    center_quotient = False

    def __init__(self, level):
        self.level = level
        reps = su3_reps(level)
        index = {rep: pos for pos, rep in enumerate(reps)}
        dim = len(reps)
        mult_fund = np.zeros((dim, dim), complex)
        for col, (p, q) in enumerate(reps):
            for target in [(p + 1, q), (p - 1, q + 1), (p, q - 1)]:
                if target in index:
                    mult_fund[index[target], col] += 1

        self.electric_diag = np.array([su3_casimir(*rep) for rep in reps])
        # Re chi_3 / 3, normalized to one at the group identity.
        self.wilson = (mult_fund + mult_fund.conj().T) / 6
        self.dim = dim

    def eigs(self, rho, g2):
        hamiltonian = np.diag(g2 * self.electric_diag).astype(complex)
        hamiltonian -= rho * self.wilson
        hermiticity_error = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
        assert hermiticity_error < 1e-10
        return np.linalg.eigvalsh(hamiltonian)[:2]


def gauge_spectra_by_source(model, beta, kappa):
    g2 = 1 / beta
    return {
        source: model.eigs(beta + kappa * source, g2)
        for source in [-4, -2, 0, 2, 4]
    }


def min_x_gap(spectra, candidates):
    ground_energy = float(spectra[4][0])
    gaps = []
    for xs in candidates:
        if xs == GROUND_X:
            continue
        gaps.append(float(spectra[sum(xs)][0] - ground_energy))
    return min(gaps)


def diagnostics(model, beta, kappa=1.0):
    spectra = gauge_spectra_by_source(model, beta, kappa)
    ground_spectrum = spectra[4]

    gauge_gap = float(ground_spectrum[1] - ground_spectrum[0])
    z_matter_gap = 2.0
    raw_x_gap = min_x_gap(spectra, x_sectors())
    closed_x_gap = min_x_gap(spectra, x_sectors(closed_flux=+1))

    quotient_candidates = [
        xs
        for xs in x_sectors(closed_flux=+1)
        if xs not in {GROUND_X, GLOBAL_FLIP_X}
    ]
    quotient_x_gap = min_x_gap(spectra, quotient_candidates)

    return {
        "gauge_gap": gauge_gap,
        "raw_x_gap": raw_x_gap,
        "closed_x_gap": closed_x_gap,
        "quotient_x_gap": quotient_x_gap,
        "raw_matter_gap": min(z_matter_gap, raw_x_gap),
        "closed_matter_gap": min(z_matter_gap, closed_x_gap),
        "quotient_matter_gap": min(z_matter_gap, quotient_x_gap),
        "raw_full_gap": min(gauge_gap, z_matter_gap, raw_x_gap),
        "quotient_full_gap": min(gauge_gap, z_matter_gap, quotient_x_gap),
    }


def report_model(model, beta_values=BETA_VALUES):
    hd(model.name)
    print(
        "  beta       gauge      raw-X   P_X=+ X  quotient-X   "
        "matter(q)   full(q)"
    )
    mins = {
        "gauge_gap": float("inf"),
        "raw_x_gap": float("inf"),
        "quotient_x_gap": float("inf"),
        "quotient_matter_gap": float("inf"),
    }

    for beta in beta_values:
        row = diagnostics(model, beta)
        for key in mins:
            mins[key] = min(mins[key], row[key])
        print(
            f"  {beta:<8g} "
            f"{row['gauge_gap']:<10.6g} "
            f"{row['raw_x_gap']:<8.6g} "
            f"{row['closed_x_gap']:<8.6g} "
            f"{row['quotient_x_gap']:<11.6g} "
            f"{row['quotient_matter_gap']:<10.6g} "
            f"{row['quotient_full_gap']:<8.6g}"
        )

    print(
        f"  minima: gauge={mins['gauge_gap']:.6g}, raw-X={mins['raw_x_gap']:.6g}, "
        f"quotient-X={mins['quotient_x_gap']:.6g}, "
        f"quotient matter={mins['quotient_matter_gap']:.6g}"
    )
    return mins


def report_su3():
    hd("SU(3) character shadow")
    print(
        "  The CSS X sectors are Z2 signs, while the SU(3) center is Z3. "
        "The quotient columns below only remove the CSS global sign partner; "
        "they are not a real SU(3) center quotient."
    )
    print(
        "  Read this column by coupling regime: the matter gap recovers at "
        "beta>=2, but collapses at beta<=1, the strong-coupling/color-confining "
        "side of the scan."
    )
    print(
        "  beta       trunc  dim    gauge      raw-X   P_X=+ X  no-flip-X  "
        "matter(no-flip)"
    )
    beta_levels = [
        (0.5, 16),
        (1, 16),
        (2, 18),
        (5, 22),
        (10, 26),
        (20, 30),
        (30, 34),
    ]
    mins = {
        "gauge_gap": float("inf"),
        "raw_x_gap": float("inf"),
        "quotient_x_gap": float("inf"),
        "quotient_matter_gap": float("inf"),
    }
    for beta, level in beta_levels:
        model = SU3Character(level=level)
        row = diagnostics(model, beta)
        for key in mins:
            mins[key] = min(mins[key], row[key])
        print(
            f"  {beta:<8g} "
            f"{level:<6d} "
            f"{model.dim:<6d} "
            f"{row['gauge_gap']:<10.6g} "
            f"{row['raw_x_gap']:<8.6g} "
            f"{row['closed_x_gap']:<8.6g} "
            f"{row['quotient_x_gap']:<10.6g} "
            f"{row['quotient_matter_gap']:<10.6g}"
        )
    print(
        f"  minima: gauge={mins['gauge_gap']:.6g}, raw-X={mins['raw_x_gap']:.6g}, "
        f"no-flip-X={mins['quotient_x_gap']:.6g}, "
        f"no-flip matter={mins['quotient_matter_gap']:.6g}"
    )
    return mins


def main():
    hd("Scope")
    print(
        "This is the first matter-coupled Gauss-projected proxy: exact CSS "
        "stabilizer-sector bookkeeping plus one-plaquette Wilson-character gauge "
        "dynamics along beta=1/g^2. It separates gauge oscillator gaps from "
        "CSS matter-sector gaps."
    )

    u1_mins = report_model(U1Rotor())
    su2_mins = report_model(SU2Character())
    su3_mins = report_su3()

    hd("Checks")
    print(
        "  Raw X-sector gaps can be small at coarse beta, so raw sector crossings "
        "remain a bad diagnostic."
    )
    print(
        "  After imposing P_X=+ and quotienting the Z2 global center partner, the "
        "U(1)/SU(2) matter gap stays close to the original CSS Z-stabilizer gap "
        "2 in this scan rather than drifting toward zero."
    )
    print(
        "  The SU(3) Pauli-sign shadow does not pass the same matter-gap test at "
        "coarse beta. It recovers at weak coupling, but the strong-coupling "
        "failure is treated as a real warning about the Z2 CSS X signs versus "
        "the Z3 SU(3) center, not as a positive SU(3) result."
    )
    print(
        "  The first full gap is often a gauge oscillator gap; the matter(q) "
        "column is the relevant mirror-gap proxy."
    )

    assert u1_mins["raw_x_gap"] < 1.0
    assert su2_mins["raw_x_gap"] < 1.0
    assert u1_mins["quotient_matter_gap"] > 1.9
    assert su2_mins["quotient_matter_gap"] > 1.9
    assert su3_mins["quotient_matter_gap"] < 1.0

    hd("Verdict")
    print(
        "The matter-coupled proxy supports the candidate pattern only in the Z2 "
        "center cases: the raw crossings are not the right physical question, "
        "while closed-flux plus Z2 center/Gauss quotient leaves a finite CSS "
        "matter gap through the scanned beta=1/g^2 trajectory for U(1)/SU(2)."
    )
    print(
        "\nThis is still not the real TCH construction. The single character-basis "
        "plaquette has already integrated out the open-link Gauss structure, and "
        "the SU(3) center mismatch is now an actual strong-coupling diagnostic: "
        "Pauli CSS X signs do not by themselves realize a Z3 center quotient or "
        "protect the color matter gap in the confining regime. This sharpens "
        "ANCHOR §2.8's open F2-colour to one-hot colour mapping into the concrete "
        "item-13 requirement: restore open links/Gauss generators and identify a "
        "Z3-aware TCH gauge-cell operator for the non-abelian Bianchi constraint "
        "on the mirror interaction."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
