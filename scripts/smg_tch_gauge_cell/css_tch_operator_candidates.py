#!/usr/bin/env python3
"""
Candidate TCH gauge-cell operators for rescuing the CSS SMG gap.

Starting point:
  css_ks_wilson_proxy_gap.py found that a raw compact-U(1) Wilson proxy has
  zero first gap because different X-stabilizer sectors cross. This script asks
  which plausible TCH operator structures could remove those crossings or make
  them gauge copies.

Candidates tested in the same one-plaquette Wilson proxy:

  A. Gauss-orbit quotient:
     local vertex gauge transformations move minus signs around the plaquette.
     The only one-plaquette invariant is P_X = prod_l x_l. If the same-P_X
     sectors are exactly isospectral, those degeneracies are pure gauge copies.

  B. Closed-flux / Bianchi plaquette projector:
     P_B = (1 + prod_l X_l)/2 keeps the zero discrete-flux sector. In a TCH
     implementation this would be the product of the four X-type gauge-cell
     stabilizers around the elementary closed plaquette, not a matter-cell
     hard projection.

  C. Soft plaquette-parity lock:
     H_lock = mu (1 - prod_l X_l)/2 energetically enforces the same sector
     instead of projecting it. This is the Hamiltonian version of B.

Important limitation:
  These are plausible algebraic candidates, not the missing explicit TCH
  non-abelian operator. The real item-13 operator still has to be identified on
  truncated-cube gauge cells and then tested in the full mirror/dynamical gauge
  problem. The positive finite-beta projected gaps below are not a continuum
  extrapolation; for the beta-limit and non-abelian Gauss-projected character
  test, see css_nonabelian_gauss_proxy.py.

numpy; self-asserting.
"""

import importlib.util
import itertools
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PROXY_PATH = HERE / "css_ks_wilson_proxy_gap.py"


def load_proxy():
    spec = importlib.util.spec_from_file_location("css_ks_wilson_proxy_gap", PROXY_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


proxy = load_proxy()


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def representatives():
    return {
        +1: (1, 1, 1, 1),
        -1: (-1, 1, 1, 1),
    }


def sector_spectrum(model, x_sector, flux, n_eigs=10):
    eigvals, hermiticity_error = model.sector_eigs(x_sector, flux, n_eigs=n_eigs)
    assert hermiticity_error < 1e-10
    return eigvals


def unique_sorted(values, tol=1e-8):
    out = []
    for value in sorted(values):
        if all(abs(value - seen) > tol for seen in out):
            out.append(float(value))
    return out


def same_parity_spread(model, flux):
    max_spread = 0.0
    for parity in [+1, -1]:
        spectra = []
        for x_sector in itertools.product([1, -1], repeat=4):
            if np.prod(x_sector) == parity:
                spectra.append(sector_spectrum(model, x_sector, flux, n_eigs=5))
        spectra = np.array(spectra)
        max_spread = max(max_spread, float(np.max(np.ptp(spectra, axis=0))))
    return max_spread


def projected_gap(model, parity, flux):
    x_sector = representatives()[parity]
    eigvals = sector_spectrum(model, x_sector, flux, n_eigs=10)

    # CSS Z half: z=++++ contributes -4. One Z stabilizer excitation costs +2.
    levels = [float(eigval - 4) for eigval in eigvals]
    levels.append(float(eigvals[0] - 2))
    levels = unique_sorted(levels)
    return levels[1] - levels[0], levels[:5]


def penalty_gap(model, flux, mu):
    levels = []
    for parity, x_sector in representatives().items():
        eigvals = sector_spectrum(model, x_sector, flux, n_eigs=10)
        penalty = mu * (1 - parity) / 2
        for eigval in eigvals:
            levels.append(float(eigval - 4 + penalty))
        levels.append(float(eigvals[0] - 2 + penalty))
    levels = unique_sorted(levels)
    return levels[1] - levels[0], levels[:5]


def scan(model, fn, n_flux=33):
    best = (float("inf"), None, None)
    for flux in np.linspace(0, 2 * np.pi, n_flux):
        gap, levels = fn(flux)
        if gap < best[0]:
            best = (gap, flux, levels)
    return best


def main():
    hd("Scope")
    print(
        "Using the smallest even compact-link truncation N_link=4 so a pi link-sign "
        "shift is represented exactly. Odd truncations are bad for this question "
        "because they cannot represent U_l -> -U_l exactly."
    )

    beta_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    n_flux = 33

    hd("A. Gauss-Orbit Quotient Candidate")
    max_spread = 0.0
    for beta in beta_values:
        model = proxy.FourLinkPlaquette(n_link=4, beta=beta, g2=1.0, kappa=1.0)
        beta_spread = 0.0
        for flux in np.linspace(0, 2 * np.pi, n_flux):
            beta_spread = max(beta_spread, same_parity_spread(model, flux))
        max_spread = max(max_spread, beta_spread)
        print(f"  beta={beta:<4}: max same-P_X spectral spread = {beta_spread:.3e}")
    assert max_spread < 1e-10
    print(
        "  -> Sectors with the same plaquette parity P_X=prod X_l are exactly "
        "isospectral in the even-link Wilson proxy. These crossings are credible "
        "Gauss copies, not physical levels."
    )

    hd("B. Closed-Flux / Bianchi Projector Candidate")
    projected_min = float("inf")
    for beta in beta_values:
        model = proxy.FourLinkPlaquette(n_link=4, beta=beta, g2=1.0, kappa=1.0)
        for parity in [+1, -1]:
            best_gap, best_flux, levels = scan(
                model, lambda flux, parity=parity: projected_gap(model, parity, flux),
                n_flux=n_flux,
            )
            projected_min = min(projected_min, best_gap)
            print(
                f"  beta={beta:<4}, P_X={parity:+}: min internal gap={best_gap:.6g} "
                f"at Phi={best_flux:.6f}; levels={np.round(levels[:4], 6)}"
            )
    assert projected_min > 1e-3
    print(
        "  -> Fixing one closed-flux sector gives a positive internal gap throughout "
        "this finite-beta scan. This is the strongest concrete candidate, but the "
        "trend is not itself a continuum-limit statement: the TCH X-plaquette "
        "operator should implement a Bianchi/Gauss closed-flux sector, and the "
        "beta-limit must be tested separately."
    )

    hd("C. Soft Plaquette-Parity Lock Candidate")
    for mu in [0.5, 1.0, 2.0, 5.0]:
        min_gap = float("inf")
        for beta in beta_values[1:]:
            model = proxy.FourLinkPlaquette(n_link=4, beta=beta, g2=1.0, kappa=1.0)
            best_gap, best_flux, levels = scan(
                model, lambda flux, mu=mu: penalty_gap(model, flux, mu), n_flux=n_flux
            )
            min_gap = min(min_gap, best_gap)
            print(
                f"  mu={mu:<3}, beta={beta:<4}: min full gap={best_gap:.6g} "
                f"at Phi={best_flux:.6f}; levels={np.round(levels[:4], 6)}"
            )
        print(f"  mu={mu:<3}: worst-case scanned gap = {min_gap:.6g}\n")

    hd("Candidate Ranking")
    print(
        "1. Best candidate: closed-flux/Bianchi X-plaquette constraint "
        "P_B=(1+prod_l X_l)/2, interpreted as a gauge-cell plaquette/Gauss "
        "condition rather than a matter-cell projection."
    )
    print(
        "2. Necessary companion: local Gauss-orbit quotient/symmetrizer. It makes "
        "same-P_X sector degeneracies pure gauge and is exactly supported by the "
        "even-link Wilson proxy."
    )
    print(
        "3. Softer Hamiltonian version: H_lock=mu(1-prod_l X_l)/2. It can bias the "
        "physical sector, but its required scale and compatibility with the CSS "
        "SMG gap still need a real TCH derivation."
    )
    print(
        "\nWhat to look for in the actual TCH construction: an X-type stabilizer on "
        "the truncated-cube gauge cell whose boundary product is prod_l X_l, with "
        "vertex Gauss generators that move signs along the plaquette but leave the "
        "boundary product fixed. That is the operator pattern this proxy selects."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
