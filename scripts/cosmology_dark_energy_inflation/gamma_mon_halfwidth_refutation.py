#!/usr/bin/env python3
"""Test the proposed half-linewidth rule for the dressed-alpha monitor clock.

Claim under test:

    gamma_mon = gap_unit / 2

for the Peierls-current bridge readout.

This script separates three different objects:

  1. gap_unit/2, half of the classical population-generator gap used in the
     Zeno/dephased bridge reduction;
  2. the Peierls-current spectral width of the Hamiltonian current vertex;
  3. the critical-damping scale of the minimal current/Krylov two-level
     reduction.

Only (2) or (3) would make the rule a physical half-linewidth theorem.  If the
number comes only from (1), it is a graph-internal clock postulate, not a
derivation of gamma_mon.
"""
from __future__ import annotations

import numpy as np

import dressed_alpha_bridge_web_open_system as bw


T_M = 1.0 / 3.0


def weighted_quantile(values: np.ndarray, weights: np.ndarray, p: float) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cdf = np.cumsum(w)
    cdf /= cdf[-1]
    return float(np.interp(p, cdf, v))


def peierls_current_spectrum(h2: np.ndarray, idx: dict[tuple[int, int], int], *, t_m: float):
    j = bw.current_portal(idx, bw.ETA_PIN)
    evals, evecs = np.linalg.eigh(t_m * h2)
    weights = (evecs.T @ j) ** 2
    weights /= weights.sum()
    mean = float(weights @ evals)
    variance = float(weights @ ((evals - mean) ** 2))
    q05 = weighted_quantile(evals, weights, 0.05)
    q25 = weighted_quantile(evals, weights, 0.25)
    q75 = weighted_quantile(evals, weights, 0.75)
    q95 = weighted_quantile(evals, weights, 0.95)
    return {
        "mean": mean,
        "rms_width": variance**0.5,
        "q05": q05,
        "q25": q25,
        "q75": q75,
        "q95": q95,
        "central_90_width": q95 - q05,
        "iqr": q75 - q25,
    }


def current_krylov_scale(h2: np.ndarray, idx: dict[tuple[int, int], int], *, t_m: float) -> float:
    """Return Omega = ||(H - <H>)|J>|| for the Peierls-current portal."""
    j = bw.current_portal(idx, bw.ETA_PIN)
    h = t_m * h2
    mean = float(j @ h @ j)
    variance = float(j @ (h @ h) @ j - mean * mean)
    return variance**0.5


def main() -> None:
    h2, _pairs, idx, _bas = bw.build_pair_system()
    c_j, gap_unit = bw.response_coefficient(h2, idx)
    half_gap = 0.5 * gap_unit

    unit_spectrum = peierls_current_spectrum(h2, idx, t_m=1.0)
    matter_spectrum = peierls_current_spectrum(h2, idx, t_m=T_M)
    omega_unit = current_krylov_scale(h2, idx, t_m=1.0)
    omega_matter = current_krylov_scale(h2, idx, t_m=T_M)

    # Minimal two-level current/Krylov reduction:
    # z'' + gamma z' + 4 Omega^2 z = 0 if gamma denotes coherence decay.
    # Thus gamma_crit = 4 Omega.  If the line is quoted by HWHM = gamma/2,
    # the critical HWHM is 2 Omega.  Either convention is far from gap_unit/2.
    gamma_crit_coherence = 4.0 * omega_matter
    gamma_crit_hwhm = 2.0 * omega_matter

    print("GAMMA_MON HALFWIDTH / CRITICAL-DAMPING TEST")
    print("=" * 78)
    print("Population-generator object:")
    print(f"  c_J                         = {c_j:.6f}")
    print(f"  gap_unit                    = {gap_unit:.6f} Lambda")
    print(f"  gap_unit / 2                = {half_gap:.6f} Lambda")
    print()
    print("Peierls-current spectral object:")
    print(f"  unit-current rms width       = {unit_spectrum['rms_width']:.6f} Lambda")
    print(f"  unit-current central 90%     = {unit_spectrum['central_90_width']:.6f} Lambda")
    print(f"  matter-current rms width     = {matter_spectrum['rms_width']:.6f} Lambda")
    print(f"  matter-current central 90%   = {matter_spectrum['central_90_width']:.6f} Lambda")
    print(f"  (gap_unit/2) / rms_matter    = {half_gap / matter_spectrum['rms_width']:.6f}")
    print()
    print("Minimal current/Krylov critical damping:")
    print(f"  Omega_unit = sqrt(<J|H^2|J>) = {omega_unit:.6f} Lambda")
    print(f"  Omega_matter                = {omega_matter:.6f} Lambda")
    print(f"  gamma_crit, coherence rate  = {gamma_crit_coherence:.6f} Lambda")
    print(f"  gamma_crit, HWHM convention = {gamma_crit_hwhm:.6f} Lambda")
    print(f"  (gap_unit/2) / gamma_crit   = {half_gap / gamma_crit_coherence:.6f}")
    print(f"  (gap_unit/2) / HWHM_crit    = {half_gap / gamma_crit_hwhm:.6f}")
    print()
    print("VERDICT")
    print("  Refuted as stated.  gap_unit/2 is half of the dephased population")
    print("  generator's unit-clock gap, not the Peierls-current spectral half-width")
    print("  and not the critical-damping scale of the current readout.  Adopting")
    print("  gamma_mon = gap_unit/2 would therefore be an extra bridge-clock postulate,")
    print("  not a derived linewidth theorem in the current service algebra.")

    assert abs(omega_unit**2 - 6.0) < 1e-12
    assert matter_spectrum["rms_width"] > 5.0 * half_gap
    assert gamma_crit_hwhm > 10.0 * half_gap
    assert gamma_crit_coherence > 20.0 * half_gap
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
