#!/usr/bin/env python3
"""
Kogut-Susskind/Wilson proxy gap scan for the CSS SMG interaction.

Why this exists:
  css_wilson_twist_gap.py showed that the undynamical cell CSS gap closes
  under simple nonuniform background twists. The next requested test is to
  replace those ad hoc diagonal twists with actual Wilson/plaquette gauge
  degrees of freedom.

What the corpus currently supplies:
  photon_paper2_final.tex gives the standard compact U(1) Kogut-Susskind form

      T_gauged = T x U_l
      H_gauge = g^2/2 sum_l E_l^2
              + 1/g^2 sum_p (1 - Re prod_{l in dp} U_l e^{i Phi_p})

  It does NOT yet supply the explicit TCH non-abelian gauge-cell plaquette
  operators that item 13 asks for. So this script is a canon-grounded proxy:
  four compact U(1) links around one plaquette are coupled to the four CSS
  X-stabilizers, and the low spectrum is scanned versus plaquette flux Phi.

Model:
  Work in the simultaneous eigenbasis of the 8 commuting CSS stabilizers.
  The four Z stabilizers contribute E_Z = -sum_i z_i. The four X stabilizers
  label sectors x_i = +/-1. In a sector x, the finite-link gauge Hamiltonian is

      H_x(Phi) = g^2/2 sum_l E_l^2
               - kappa/2 sum_l x_l (U_l + U_l^dagger)
               - beta/2 (e^{i Phi} U0 U1 U2^dag U3^dag + h.c.).

  This is not the full SM chiral gauge theory. It is the smallest executable
  Wilson/plaquette replacement for the background-twist scan.

Interpretation:
  If the gap closes here, the Wilson proxy does not rescue the cell SMG gap.
  If it stayed open, that would only be evidence for this compact U(1) proxy,
  not for the missing TCH non-abelian item-13 construction.

numpy; self-asserting.
"""

import itertools
import numpy as np


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


class FourLinkPlaquette:
    def __init__(self, n_link=3, kappa=1.0, beta=1.0, g2=1.0):
        self.n_link = n_link
        self.kappa = kappa
        self.beta = beta
        self.g2 = g2
        self.identity = np.eye(n_link, dtype=complex)
        self.u_link, self.e_link = self._single_link_ops(n_link)
        self.u_ops = [
            self._kron4(self.u_link, self.identity, self.identity, self.identity),
            self._kron4(self.identity, self.u_link, self.identity, self.identity),
            self._kron4(self.identity, self.identity, self.u_link, self.identity),
            self._kron4(self.identity, self.identity, self.identity, self.u_link),
        ]
        self.e_ops = [
            self._kron4(self.e_link, self.identity, self.identity, self.identity),
            self._kron4(self.identity, self.e_link, self.identity, self.identity),
            self._kron4(self.identity, self.identity, self.e_link, self.identity),
            self._kron4(self.identity, self.identity, self.identity, self.e_link),
        ]
        self.electric = sum((g2 / 2) * (e_op @ e_op) for e_op in self.e_ops)
        self.link_cos = [(u_op + u_op.conj().T) / 2 for u_op in self.u_ops]
        self.plaquette = (
            self.u_ops[0]
            @ self.u_ops[1]
            @ self.u_ops[2].conj().T
            @ self.u_ops[3].conj().T
        )

    @staticmethod
    def _single_link_ops(n_link):
        shift = np.zeros((n_link, n_link), complex)
        for state in range(n_link):
            shift[(state + 1) % n_link, state] = 1
        electric_values = np.arange(n_link) - n_link // 2
        electric = np.diag(electric_values.astype(float))
        return shift, electric

    @staticmethod
    def _kron4(a, b, c, d):
        return np.kron(np.kron(np.kron(a, b), c), d)

    def hamiltonian(self, x_sector, flux):
        h = self.electric.copy()
        for sign, cos_link in zip(x_sector, self.link_cos):
            h -= self.kappa * sign * cos_link
        h -= (self.beta / 2) * (
            np.exp(1j * flux) * self.plaquette
            + np.exp(-1j * flux) * self.plaquette.conj().T
        )
        return h

    def sector_eigs(self, x_sector, flux, n_eigs=4):
        h = self.hamiltonian(x_sector, flux)
        hermiticity_error = float(np.linalg.norm(h - h.conj().T))
        eigvals = np.linalg.eigvalsh(h)
        return eigvals[:n_eigs], hermiticity_error


def low_levels(model, flux):
    levels = []
    for x_sector in itertools.product([1, -1], repeat=4):
        eigvals, hermiticity_error = model.sector_eigs(x_sector, flux)
        assert hermiticity_error < 1e-10

        # Z-stabilizer energy: all z=+ gives -4. A single Z excitation costs +2.
        for eig in eigvals:
            levels.append((float(eig - 4), x_sector, "gauge, z=++++"))
        levels.append((float(eigvals[0] - 2), x_sector, "first Z-stab excitation"))

    return sorted(levels, key=lambda item: item[0])


def gap_at(model, flux):
    levels = low_levels(model, flux)
    return levels[1][0] - levels[0][0], levels[:8]


def flux_scan(n_link=3, beta=1.0, g2=1.0, kappa=1.0, n_flux=33):
    model = FourLinkPlaquette(n_link=n_link, kappa=kappa, beta=beta, g2=g2)
    best = (float("inf"), None, None)
    values = []
    for flux in np.linspace(0, 2 * np.pi, n_flux):
        gap, levels = gap_at(model, flux)
        values.append(gap)
        if gap < best[0]:
            best = (gap, flux, levels)
    return best, values


def report_scan(n_link, beta, g2=1.0, kappa=1.0, n_flux=33):
    best, values = flux_scan(
        n_link=n_link, beta=beta, g2=g2, kappa=kappa, n_flux=n_flux
    )
    gap, flux, levels = best
    print(
        f"  N_link={n_link}, beta={beta:<4}, g2={g2:<4}, kappa={kappa:<4} "
        f"-> min gap={gap:.12g} at Phi={flux:.6f}"
    )
    print(
        "     lowest levels: "
        + ", ".join(
            f"{energy:.6f} x={sector} {tag}" for energy, sector, tag in levels[:3]
        )
    )
    return gap


def main():
    hd("Scope")
    print(
        "This is a compact-U(1) Kogut-Susskind proxy using the Wilson plaquette "
        "form written in photon_paper2_final.tex. It is not the missing explicit "
        "TCH non-abelian item-13 plaquette construction."
    )

    hd("Flux Scan")
    gaps = []
    for beta in [0.0, 0.5, 1.0, 2.0, 5.0]:
        gaps.append(report_scan(n_link=3, beta=beta, n_flux=33))

    hd("Stronger Link Truncation Spot Check")
    for beta in [0.5, 1.0, 2.0]:
        gaps.append(report_scan(n_link=5, beta=beta, n_flux=9))

    hd("Verdict")
    min_gap = min(gaps)
    print(f"  minimum gap observed across proxy scans: {min_gap:.12g}")
    assert min_gap < 1e-8

    print(
        "\nIn this compact-U(1) Wilson proxy, the dynamical plaquette/link degrees "
        "do not rescue the cell CSS gap. With beta > 0, finite plaquette fluxes "
        "produce exact or numerical-zero degeneracies between distinct X-stabilizer "
        "sectors. The unique CSS ground state is therefore not robust in this "
        "minimal Wilson replacement either."
    )
    print(
        "\nThis is useful but narrow: it tests a standard Wilson U(1) plaquette "
        "proxy coupled to CSS X sectors. Because the actual TCH gauge-cell "
        "X-stabilizers and Gauss-law physical subspace are still unidentified, "
        "this is not a no-go theorem for item 13. It is a concrete warning that "
        "generic Wilson dressing is not enough; the real TCH operator must remove "
        "or gauge-identify these sector crossings and then pass the gap scan."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
