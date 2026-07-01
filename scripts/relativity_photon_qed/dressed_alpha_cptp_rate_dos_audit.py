#!/usr/bin/env python3
"""Rate and DOS-convergence audit for the dressed-alpha CPTP near-closure.

The combined monitor-plus-web CPTP audit found a near-closure:

    native item-79 pair-basis monitor, rate = 1
    + full Peierls/Wilson web width operator
    -> delta/target ~= 0.99.

This script checks the two load-bearing quantitative questions:

  1. At monitor rate = 1, does the Wilson DOS continuum extrapolate toward
     the observed shift or plateau below it?
  2. How special is rate = 1?  Where does delta/target cross 1?

It uses the py13_7 SciPy sparse eigensolver through
dressed_alpha_monitor_web_cptp_audit.py.  The Wilson DOS is accumulated
slab-by-slab to avoid the full 3D mesh allocation and to make larger n_grid
checks cheap enough to run locally.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_cptp_audit as cptp


NATIVE_RATE = 1.0


@dataclass(frozen=True)
class NativeResult:
    rate: float
    gamma_web: float
    pole: float
    delta: float
    ratio: float
    residual: float
    mode: str
    solver: str


def streamed_photon_form_factor(n_grid: int, bins: int = 400):
    """Wilson photon DOS/form-factor estimate without materialising the 3D mesh."""
    kk = (np.arange(n_grid) + 0.5) * 2.0 * np.pi / n_grid
    sy2 = np.sin(kk / 2.0) ** 2
    yz = sy2[:, None] + sy2[None, :]
    if n_grid % 2 == 0:
        omega_max = 2.0 * math.sqrt(3.0) * math.cos(math.pi / (2.0 * n_grid))
    else:
        omega_max = 2.0 * math.sqrt(3.0)
    hist_n = np.zeros(bins, dtype=float)
    hist_c = np.zeros(bins, dtype=float)
    hist_range = (0.0, omega_max + 1e-9)

    for sx2, kx in zip(sy2, kk):
        omega = np.sqrt(4.0 * (sx2 + yz)).ravel()
        counts, edges = np.histogram(omega, bins=bins, range=hist_range)
        hist_n += counts
        hist_c += math.cos(float(kx)) * counts

    bin_width = (hist_range[1] - hist_range[0]) / bins
    density = hist_n / (float(n_grid) ** 3 * bin_width)
    mean_c = np.divide(hist_c, hist_n, out=np.zeros_like(hist_c), where=hist_n > 0)
    centers = 0.5 * (edges[:-1] + edges[1:])

    def rho(w: float) -> float:
        if not (0.0 < w < omega_max):
            return 0.0
        return float(np.interp(w, centers, density))

    def s_eta(w: float, eta: int) -> float:
        if not (0.0 < w < omega_max):
            return 0.0
        r = float(np.interp(w, centers, density))
        c = r * float(np.interp(w, centers, mean_c))
        return max(r + eta * c, 0.0)

    soft_ratio = s_eta(0.15, -1) / max(rho(0.15), 1e-300)
    return s_eta, omega_max, soft_ratio


def native_result(
    h: np.ndarray,
    gamma_web: np.ndarray,
    pairs: list[tuple[int, int]],
    idx: dict[tuple[int, int], int],
    rate: float,
) -> NativeResult:
    labels = cptp.pair_basis_labels(pairs, idx)
    current_minus = bw.current_portal(idx, -1)
    current_plus = bw.current_portal(idx, +1)
    candidate = cptp.MonitorCandidate(
        "native pair-basis monitor",
        "pair_site",
        rate,
        "rate scan",
        "item-79 site/pair dephasing at scanned rate",
    )
    matvec = cptp.make_matvec(h, gamma_web, candidate, labels, current_minus, current_plus)
    eigenvalue, eigenvector, residual, solver = cptp.scipy_dominant(matvec, len(idx) * len(idx))
    obs = cptp.qss_observables(eigenvalue, eigenvector, idx, len(idx))
    mode = cptp.row_mode(obs, residual)
    return NativeResult(
        rate=rate,
        gamma_web=float(np.trace(gamma_web)),
        pole=obs["pole"],
        delta=obs["delta"],
        ratio=obs["delta"] / bw.DELTA_TARGET,
        residual=residual,
        mode=mode,
        solver=solver,
    )


def width_from_streamed_dos(h2: np.ndarray, idx: dict[tuple[int, int], int], n_grid: int, bins: int):
    s_eta, _omega_max, soft_ratio = streamed_photon_form_factor(n_grid=n_grid, bins=bins)
    gamma_web, gamma_escape, bright_count, dark_dim = cptp.grouped_web_width_operator(h2, idx, s_eta)
    assert bright_count == 32
    assert dark_dim == 104
    return gamma_web, gamma_escape, soft_ratio


def fit_continuum(grid_rows: list[tuple[int, int, NativeResult]]) -> tuple[float, float]:
    """Fit ratio(n) = ratio_inf + a/n^2.  This is a diagnostic, not a proof."""
    x = np.array([1.0 / (n * n) for n, _bins, _res in grid_rows])
    y = np.array([res.ratio for _n, _bins, res in grid_rows])
    slope, intercept = np.polyfit(x, y, 1)
    return float(intercept), float(slope)


def find_crossing(
    h: np.ndarray,
    gamma_web: np.ndarray,
    pairs: list[tuple[int, int]],
    idx: dict[tuple[int, int], int],
    low: float,
    high: float,
    *,
    steps: int = 8,
) -> tuple[float, NativeResult]:
    lo_res = native_result(h, gamma_web, pairs, idx, low)
    hi_res = native_result(h, gamma_web, pairs, idx, high)
    assert lo_res.mode == "qss" and hi_res.mode == "qss"
    assert (lo_res.ratio - 1.0) * (hi_res.ratio - 1.0) <= 0.0
    best = lo_res
    for _ in range(steps):
        mid = 0.5 * (low + high)
        mid_res = native_result(h, gamma_web, pairs, idx, mid)
        assert mid_res.mode == "qss"
        best = mid_res
        if (lo_res.ratio - 1.0) * (mid_res.ratio - 1.0) <= 0.0:
            high = mid
            hi_res = mid_res
        else:
            low = mid
            lo_res = mid_res
    return best.rate, best


def main() -> None:
    h2, pairs, idx, _bas = bw.build_pair_system()
    h = cptp.T_M * h2

    print("DRESSED-ALPHA CPTP RATE + DOS CONVERGENCE AUDIT")
    print("=" * 104)
    print("Native monitor under test:")
    print("  instrument = item-79 pair-basis dephasing")
    print("  rate       = scanned; rate=1 is native substrate tick")
    print("  web loss   = full Peierls/Wilson width operator, not scalar gamma_mon")
    print("  solver     = py13_7 SciPy sparse eigs via matrix-vector products")

    print("\n[1] DOS/grid convergence at native rate = 1")
    print(f"  {'n_grid':>7} {'bins':>5} {'soft S-/rho':>11} {'Gamma_web':>12} {'delta':>12} {'/target':>10} {'pole':>12} {'resid':>9} mode")
    grid_rows: list[tuple[int, int, NativeResult]] = []
    for n_grid, bins in ((160, 400), (200, 400), (240, 400), (280, 400), (320, 400)):
        gamma_web, gamma_escape, soft_ratio = width_from_streamed_dos(h2, idx, n_grid, bins)
        result = native_result(h, gamma_web, pairs, idx, NATIVE_RATE)
        grid_rows.append((n_grid, bins, result))
        print(
            f"  {n_grid:>7d} {bins:>5d} {soft_ratio:>11.4f} {gamma_escape:>12.6e}"
            f" {result.delta:>12.6e} {result.ratio:>10.6f} {result.pole:>12.6e}"
            f" {result.residual:>9.1e} {result.mode}"
        )
        assert result.mode == "qss"
        assert result.solver == "scipy-eigs"

    ratio_inf, slope = fit_continuum(grid_rows[-4:])
    print(f"\n  linear fit over n=200..320 in 1/n^2: ratio_inf = {ratio_inf:.6f}, slope = {slope:.3e}")
    print("  Treat this as a convergence diagnostic only: histogram binning remains visible.")

    print("\n[2] Histogram-bin sensitivity at n_grid = 240, rate = 1")
    print(f"  {'bins':>5} {'Gamma_web':>12} {'delta':>12} {'/target':>10} {'resid':>9} mode")
    bin_rows: list[tuple[int, NativeResult]] = []
    for bins in (300, 400, 600, 800):
        gamma_web, gamma_escape, _soft_ratio = width_from_streamed_dos(h2, idx, 240, bins)
        result = native_result(h, gamma_web, pairs, idx, NATIVE_RATE)
        bin_rows.append((bins, result))
        print(
            f"  {bins:>5d} {gamma_escape:>12.6e} {result.delta:>12.6e}"
            f" {result.ratio:>10.6f} {result.residual:>9.1e} {result.mode}"
        )
        assert result.mode == "qss"

    print("\n[3] Native-monitor rate scan at n_grid = 240, bins = 400")
    gamma_web_240, gamma_escape_240, _soft_ratio_240 = width_from_streamed_dos(h2, idx, 240, 400)
    print(f"  fixed Gamma_web = {gamma_escape_240:.6e}")
    print(f"  {'rate':>9} {'delta':>12} {'/target':>10} {'pole':>12} {'resid':>9} mode")
    scan_rates = (0.25, 1.0 / 3.0, 0.45, 0.60, 0.75, 0.90, 1.00, 1.10, 1.25, 1.50, 2.00)
    scan_rows: list[NativeResult] = []
    for rate in scan_rates:
        result = native_result(h, gamma_web_240, pairs, idx, rate)
        scan_rows.append(result)
        print(
            f"  {rate:>9.6f} {result.delta:>12.6e} {result.ratio:>10.6f}"
            f" {result.pole:>12.6e} {result.residual:>9.1e} {result.mode}"
        )
        assert result.mode == "qss"

    brackets = []
    for left, right in zip(scan_rows, scan_rows[1:]):
        if (left.ratio - 1.0) * (right.ratio - 1.0) <= 0.0:
            brackets.append((left.rate, right.rate))
    assert brackets
    crossings = [find_crossing(h, gamma_web_240, pairs, idx, left, right) for left, right in brackets]

    print("\n[4] Rate crossing")
    print(f"  {'bracket':<23} {'crossing':>10} {'delta/target':>13} {'native/cross':>13} {'cross-native':>13}")
    for (left, right), (crossing_rate, crossing_result) in zip(brackets, crossings):
        print(
            f"  [{left:.6f}, {right:.6f}]"
            f" {crossing_rate:>10.6f} {crossing_result.ratio:>13.6f}"
            f" {NATIVE_RATE / crossing_rate:>13.6f}"
            f" {crossing_rate - NATIVE_RATE:>+13.6f}"
        )

    print("\nVERDICT")
    print("  The native-rate result is genuinely rate-sensitive.  At fixed n=240,")
    print("  the rate response is U-shaped and crosses the observed target twice;")
    print("  rate=1 lies close to the upper crossing but is not itself forced by the")
    print("  CPTP calculation.")
    print("  The DOS sequence and bin sweep both move the native-rate answer through")
    print("  the one-percent band, so the current histogram DOS cannot decide whether")
    print("  the continuum limit lands exactly on 1.0 or plateaus slightly low/high.")
    print("  Honest next step: replace histogram DOS with an adaptive/smoothed")
    print("  continuum evaluation of S_eta at the 32 bright-shell transition energies,")
    print("  then rerun the rate-crossing audit.")

    assert all(0.97 < res.ratio < 1.02 for _n, _bins, res in grid_rows)
    assert all(abs(result.ratio - 1.0) < 0.003 for _rate, result in crossings)
    assert any(0.5 < rate < 0.7 for rate, _result in crossings)
    assert any(1.0 < rate < 1.2 for rate, _result in crossings)
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
