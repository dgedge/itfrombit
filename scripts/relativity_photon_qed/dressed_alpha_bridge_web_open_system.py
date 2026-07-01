#!/usr/bin/env python3
"""Dressed-alpha bridge-to-web open-system calculation.

This is the next-step version of the non-unital escape route.  It does not use
the old N1=31 count as a loop integral.  Instead it computes the photon escape
rate from the bridge into the 3D Wilson photon web using the Peierls current
portal, then feeds that rate into the non-unital response coefficient.

Pinned here, from substrate arguments rather than alpha matching:
  * portal: Peierls current on the two bridge links, not the older flat portal;
  * eta=-1: the pi/4-per-edge bridge flux gives a pi bridge-cycle phase;
  * E_ref=band bottom: Keldysh/Fermi-golden-rule energies are excitation
    energies above the confined pair ground state;
  * g=1/3: Grover gauge-matter portal coupling.

Still not pinned by the present canon:
  * t_m/t_web, the matter-pair hopping scale relative to the Wilson photon web;
  * gamma_mon, the monitoring/dephasing clock that sets the bridge relaxation
    gap when t_m != t_web.

The script therefore reports only the strict pinned cases and the remaining
one-parameter scaling law.  It is a progress audit, not a closure.
"""
from __future__ import annotations

import numpy as np


ALPHA0 = 1.0 / 137.0
DELTA_TARGET = 0.035999084
ETA_PIN = -1
G_PORTAL = 1.0 / 3.0


def build_pair_system() -> tuple[np.ndarray, list[tuple[int, int]], dict[tuple[int, int], int], np.ndarray]:
    n = 16
    edges = [(i, (i + 1) % 8) for i in range(8)]
    edges += [(8 + i, 8 + (i + 1) % 8) for i in range(8)]
    edges += [(0, 8), (1, 9)]

    a = np.zeros((n, n))
    for i, j in edges:
        a[i, j] = a[j, i] = 1.0

    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    idx = {p: k for k, p in enumerate(pairs)}
    bas = np.zeros((n * n, len(pairs)))
    for (i, j), k in idx.items():
        v = np.zeros((n, n))
        if i == j:
            v[i, i] = 1.0
        else:
            v[i, j] = v[j, i] = 1.0 / np.sqrt(2.0)
        bas[:, k] = v.reshape(-1)

    h2 = bas.T @ (np.kron(a, np.eye(n)) + np.kron(np.eye(n), a)) @ bas
    return h2, pairs, idx, bas


def current_portal(idx: dict[tuple[int, int], int], eta: int) -> np.ndarray:
    """Peierls-current portal on the two physical bridge links."""
    v = np.zeros(len(idx))
    v[idx[(0, 8)]] = 1.0 / np.sqrt(2.0)
    v[idx[(1, 9)]] = eta / np.sqrt(2.0)
    return v


def photon_form_factor(n_grid: int = 160):
    """Wilson photon DOS and two-link interference factor on the SC gauge web."""
    kk = (np.arange(n_grid) + 0.5) * 2.0 * np.pi / n_grid
    kx, ky, kz = np.meshgrid(kk, kk, kk, indexing="ij")
    omega = np.sqrt(4.0 * (np.sin(kx / 2.0) ** 2 + np.sin(ky / 2.0) ** 2 + np.sin(kz / 2.0) ** 2)).ravel()
    cosx = np.cos(kx).ravel()

    omega_max = float(omega.max())
    bins = 400
    hist, edges = np.histogram(omega, bins=bins, range=(0.0, omega_max + 1e-9), density=True)
    hist_c, _ = np.histogram(omega, bins=bins, range=(0.0, omega_max + 1e-9), weights=cosx)
    hist_n, _ = np.histogram(omega, bins=bins, range=(0.0, omega_max + 1e-9))
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean_c = np.where(hist_n > 0, hist_c / np.maximum(hist_n, 1), 0.0)

    def rho(w: float) -> float:
        if not (0.0 < w < omega_max):
            return 0.0
        return float(np.interp(w, centers, hist))

    def s_eta(w: float, eta: int) -> float:
        if not (0.0 < w < omega_max):
            return 0.0
        r = float(np.interp(w, centers, hist))
        c = r * float(np.interp(w, centers, mean_c))
        return max(r + eta * c, 0.0)

    # Low-energy pi-flux branch suppresses soft emission, as expected.
    assert s_eta(0.15, -1) < 0.15 * rho(0.15)
    assert s_eta(0.15, +1) > 1.7 * rho(0.15)
    return s_eta, omega_max


def response_coefficient(h2: np.ndarray, idx: dict[tuple[int, int], int]) -> tuple[float, float]:
    """Non-unital response coefficient for loss through the two current links."""
    i08, i19 = idx[(0, 8)], idx[(1, 9)]
    w = np.abs(h2) ** 2
    np.fill_diagonal(w, 0.0)
    gen = w - np.diag(w.sum(1))
    gap_unit = 2.0 * (-np.sort(np.linalg.eigvalsh(gen))[-2])

    def qss_2site(g_esc: float) -> float:
        m = gen.copy()
        m[i08, i08] -= g_esc / 2.0
        m[i19, i19] -= g_esc / 2.0
        ev, vec = np.linalg.eig(m)
        k = np.argmax(ev.real)
        v = np.abs(vec[:, k].real)
        v /= v.sum()
        return (v[i08] + v[i19]) / 2.0

    d = len(idx)
    deltas = []
    for g_esc in (1e-4, 3e-4, 1e-3):
        p_occ = qss_2site(g_esc)
        alpha_eff = p_occ / (1.0 / d) * ALPHA0
        deltas.append((g_esc / gap_unit, 1.0 / alpha_eff - 137.0))

    c_j = np.polyfit([x for x, _ in deltas], [y for _, y in deltas], 1)[0]
    assert 4.5 < c_j < 5.5
    return float(c_j), float(gap_unit)


def gamma_escape(
    h2: np.ndarray,
    idx: dict[tuple[int, int], int],
    s_eta,
    *,
    t_m: float,
    eta: int,
    g_portal: float,
) -> tuple[float, float, float]:
    """Keldysh/FGR escape rate from pair eigenstates into the photon web."""
    evals_unit, evecs = np.linalg.eigh(h2)
    j = current_portal(idx, eta)
    overlaps = (evecs.T @ j) ** 2
    evals = t_m * evals_unit
    e_ref = float(evals.min())  # pinned band-bottom excitation reference
    omega = evals - e_ref
    spectral = np.array([overlaps[k] * s_eta(float(omega[k]), eta) for k in range(len(evals))])
    total = float(spectral.sum())
    gamma = 2.0 * np.pi * ALPHA0 * g_portal * g_portal * total
    mean_omega = float((spectral @ omega) / total) if total > 0.0 else float("nan")
    return gamma, mean_omega, total


def main() -> None:
    h2, _pairs, idx, _bas = build_pair_system()
    s_eta, omega_max = photon_form_factor()
    c_j, gap_unit = response_coefficient(h2, idx)

    j = current_portal(idx, ETA_PIN)
    portal_diag = float(j @ h2 @ j)
    assert abs(portal_diag) < 1e-12

    print("DRESSED ALPHA BRIDGE-WEB OPEN SYSTEM")
    print("=" * 78)
    print("Pinned ingredients:")
    print("  portal        = Peierls current on bridge links (0,8) and (1,9)")
    print("  eta           = -1 (pi bridge-cycle phase)")
    print("  E_ref         = pair-band bottom (excitation-energy reference)")
    print("  g_portal      = 1/3 (Grover gauge-matter portal)")
    print(f"  photon band   = Wilson web, omega_max = {omega_max:.4f} Lambda")
    print(f"  response      = c_J {c_j:.3f}; unit bridge gap {gap_unit:.4f} Lambda")
    print("  check         = <J|H_pair|J> = 0 exactly for the Peierls current portal")

    print("\nDirect Keldysh/FGR escape rates with consistent bridge-gap scaling")
    print("  gap(t_m,gamma) = gap_unit * t_m^2 / gamma_mon")
    print(f"  {'case':<26} {'Gamma_esc':>11} {'<omega>':>9} {'gap':>11} {'delta':>11} {'/target':>9}")

    rows = []
    for label, t_m, gamma_mon in (
        ("shared clock t_m=1", 1.0, 1.0),
        ("Grover, global monitor", 1.0 / 3.0, 1.0),
        ("Grover, sector monitor", 1.0 / 3.0, 1.0 / 3.0),
        ("Grover, legacy scaling", 1.0 / 3.0, 1.0 / 9.0),
    ):
        g_esc, mean_omega, _weight = gamma_escape(h2, idx, s_eta, t_m=t_m, eta=ETA_PIN, g_portal=G_PORTAL)
        gap = gap_unit * t_m * t_m / gamma_mon
        delta = c_j * g_esc / gap
        rows.append((label, g_esc, mean_omega, gap, delta))
        print(f"  {label:<26} {g_esc:>11.3e} {mean_omega:>9.3f} {gap:>11.3e} {delta:>11.3e} {delta / DELTA_TARGET:>9.2f}")

    print("\nRemaining single scaling target")
    g_esc, mean_omega, _weight = gamma_escape(h2, idx, s_eta, t_m=1.0 / 3.0, eta=ETA_PIN, g_portal=G_PORTAL)
    target_gap = c_j * g_esc / DELTA_TARGET
    gamma_required = gap_unit * (1.0 / 3.0) ** 2 / target_gap
    print(f"  With t_m=1/3, eta=-1, E_ref=band-bottom, g=1/3:")
    print(f"    Gamma_esc       = {g_esc:.3e} Lambda")
    print(f"    <omega>_emit    = {mean_omega:.3f} Lambda")
    print(f"    gap required    = {target_gap:.3e} Lambda")
    print(f"    gamma_mon needed= {gamma_required:.3f} Lambda")
    print("    This is the only remaining continuous scale in this pinned calculation.")

    assert rows[0][4] < DELTA_TARGET
    assert rows[1][4] > DELTA_TARGET
    assert 0.05 < gamma_required < 1.0
    print("\nVERDICT")
    print("  The count-to-integral route is gone.  The Peierls-current bridge/web")
    print("  calculation gives a direct Gamma_esc, pins eta and E_ref, and exposes the")
    print("  real unresolved object: the monitored bridge relaxation clock when")
    print("  t_m != t_web.  Deriving gamma_mon from the substrate would make this a")
    print("  prediction; choosing it from the target would be another fit.")
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
