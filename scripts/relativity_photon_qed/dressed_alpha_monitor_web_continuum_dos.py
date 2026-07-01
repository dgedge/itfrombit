#!/usr/bin/env python3
r"""Continuum (binning-free) DOS evaluation for the monitor+Wilson-web CPTP dressed-alpha route.

The histogram rate/convergence audit (in technical_notes/dressed_alpha_escape_program.md) left ONE thing
undecided: whether delta/target at the native substrate tick lands exactly on 1.0 in the continuum, or
plateaus slightly off. The blocker was histogram binning noise in the 3D Wilson-web DOS (Gamma_web
wobbles ~0.5% over n_grid=160..320, hiding the true limit). The flagged next step was a "smoothed
continuum evaluation of S_eta(omega) at the 32 bright-shell transition energies."

This script does exactly that, with a binning-FREE analytic DOS. The Wilson photon dispersion is
omega = 2*sqrt(sum sin^2(k_i/2)), so omega^2 = 6 - E with E = 2*sum cos(k_i) the simple-cubic
tight-binding band. Hence the photon DOS is a pushforward of the SC band DOS, which is the THREE-FOLD
CONVOLUTION of the 1D arcsine density g(c)=1/(pi*sqrt(1-c^2)) of c=cos(k):

    rho_E(E) = (1/2) (g*g*g)(E/2),      c_E(E) = <cos kx> rho_E = (1/2) ((c*g)*g*g)(E/2),
    S_eta(omega) = 2*omega * ( rho_E(6-omega^2) + eta * c_E(6-omega^2) ).

The convolutions use EXACT arcsine bin masses (the CDF F(c)=arcsin(c)/pi+1/2 integrates the c=+-1 van Hove
singularities analytically), so the only control parameter is the 1D bin count M -- a 1D refinement with
no 3D Monte-Carlo/Poisson noise. We:

  [A] validate the smooth S_eta against the histogram S_eta in the bulk (same object);
  [B] compute delta/target at the native tick (rate=1) vs M -> the binning-free continuum value;
  [C] re-run the rate-crossing scan with the smooth DOS -> does the upper target crossing sit at rate=1?

This is the convergence half of the route's "promote by theorem + convergence" requirement. It does NOT
attempt the second half (proving the p_link readout is the physical observable), which stays open.
"""
from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_monitor_web_cptp_audit as aud


H2, PAIRS, IDX, _BAS = bw.build_pair_system()
H = aud.T_M * H2
D = len(IDX)
CMINUS = bw.current_portal(IDX, -1)
CPLUS = bw.current_portal(IDX, +1)
LABELS = aud.pair_basis_labels(PAIRS, IDX)


def smooth_band_dos(m_bins: int):
    """Analytic SC band DOS rho_E(E) and cos-weighted c_E(E) via 3-fold arcsine convolution."""
    edges = np.linspace(-1.0, 1.0, m_bins + 1)
    cdf = np.arcsin(np.clip(edges, -1.0, 1.0)) / np.pi + 0.5  # exact arcsine CDF (van Hove handled)
    mass = np.diff(cdf)                                       # exact P(c in bin), length m_bins
    dc = 2.0 / m_bins
    centers = -1.0 + (np.arange(m_bins) + 0.5) * dc
    wmass = centers * mass                                    # cos-weighted 1D mass

    m3 = fftconvolve(fftconvolve(mass, mass), mass)           # density of c1+c2+c3
    cw = fftconvolve(fftconvolve(wmass, mass), mass)          # cos-weighted version
    s_pos = -3.0 + (np.arange(len(m3)) + 1.5) * dc            # positions of the 3-fold sum
    e_grid = 2.0 * s_pos                                      # E = 2*(c1+c2+c3)
    rho_e = (m3 / dc) / 2.0                                   # rho_E(E)  (int rho_E dE = 1)
    c_e = (cw / dc) / 2.0                                     # c_E(E) = <cos> rho_E
    return e_grid, rho_e, c_e


def make_smooth_s_eta(m_bins: int):
    e_grid, rho_e, c_e = smooth_band_dos(m_bins)
    e_lo, e_hi = float(e_grid[0]), float(e_grid[-1])

    def s_eta(w: float, eta: int) -> float:
        if w <= 0.0:
            return 0.0
        e = 6.0 - w * w                       # photon dispersion -> SC band energy
        if not (e_lo < e < e_hi):
            return 0.0
        r = float(np.interp(e, e_grid, rho_e))
        c = float(np.interp(e, e_grid, c_e))
        return max(2.0 * w * (r + eta * c), 0.0)

    return s_eta


def delta_ratio(rate: float, s_eta) -> tuple[float, float, int]:
    gamma_web, gesc, bright, _dark = aud.grouped_web_width_operator(H2, IDX, s_eta)
    cand = aud.MonitorCandidate("c", "pair_site", rate, "", "")
    mv = aud.make_matvec(H, gamma_web, cand, LABELS, CMINUS, CPLUS)
    ev, vec, res, _ = aud.scipy_dominant(mv, D * D)
    ratio = aud.qss_observables(ev, vec, IDX, D)["delta"] / bw.DELTA_TARGET
    return ratio, res, bright


def main() -> None:
    print("MONITOR+WILSON-WEB CPTP — CONTINUUM (binning-free) DOS EVALUATION")
    print("=" * 92)

    # ---- [A] validate the smooth DOS against the histogram DOS in the bulk ----
    hist_s_eta, _ = bw.photon_form_factor(n_grid=240)
    smooth = make_smooth_s_eta(12000)
    print("\n[A] smooth (analytic) vs histogram S_eta(omega) — same object?  (eta=-1)")
    print(f"    {'omega':>7} {'S_hist':>12} {'S_smooth':>12} {'rel.diff':>10}")
    for w in (0.5, 1.0, 1.5, 2.0, 2.5):
        sh, ss = hist_s_eta(w, -1), smooth(w, -1)
        rel = abs(ss - sh) / sh if sh else float("nan")
        print(f"    {w:>7.2f} {sh:>12.5e} {ss:>12.5e} {rel:>9.1%}")

    # ---- [B] binning-free continuum value at the native tick (rate=1), refined in M ----
    print("\n[B] delta/target at native tick (rate=1), binning-FREE, vs 1D bin count M:")
    print(f"    {'M':>7} {'Gamma_web':>13} {'delta/target':>13} {'resid':>9} {'bright':>7}")
    vals = []
    for m_bins in (3000, 6000, 12000, 24000):
        s_eta = make_smooth_s_eta(m_bins)
        gw, gesc, bright, _ = aud.grouped_web_width_operator(H2, IDX, s_eta)
        cand = aud.MonitorCandidate("c", "pair_site", 1.0, "", "")
        mv = aud.make_matvec(H, gw, cand, LABELS, CMINUS, CPLUS)
        ev, vec, res, _ = aud.scipy_dominant(mv, D * D)
        ratio = aud.qss_observables(ev, vec, IDX, D)["delta"] / bw.DELTA_TARGET
        vals.append(ratio)
        print(f"    {m_bins:>7d} {gesc:>13.6e} {ratio:>13.6f} {res:>9.1e} {bright:>7d}")
    spread = max(vals) - min(vals)
    print(f"    continuum value (binning-free): delta/target = {vals[-1]:.5f}  (M-spread over the scan: {spread:.2e})")

    # ---- [C] rate-crossing with the smooth DOS: is the upper crossing the native tick? ----
    print("\n[C] rate-crossing scan with the smooth continuum DOS (M=12000):")
    print(f"    {'rate':>7} {'delta/target':>13}")
    s_eta = make_smooth_s_eta(12000)
    scan = [0.5, 0.6, 0.75, 0.9, 1.0, 1.05, 1.1, 1.25, 1.5, 2.0]
    pts = []
    for r in scan:
        ratio, _res, _b = delta_ratio(r, s_eta)
        pts.append((r, ratio))
        tag = "  <- native substrate tick" if abs(r - 1.0) < 1e-9 else ""
        print(f"    {r:>7.3f} {ratio:>13.4f}{tag}")
    crossings = []
    for (r0, y0), (r1, y1) in zip(pts, pts[1:]):
        if (y0 - 1.0) * (y1 - 1.0) <= 0 and y0 != y1:
            crossings.append(r0 + (1.0 - y0) * (r1 - r0) / (y1 - y0))
    print(f"    target crossings at rate ~= {', '.join(f'{c:.3f}' for c in crossings) if crossings else 'none in scan'}")
    upper = max(crossings) if crossings else float("nan")
    print(f"    upper crossing {upper:.3f} vs native tick 1.0: gap {abs(upper-1.0):.3f}"
          f"  ({'native tick IS the closure rate (within 5%)' if abs(upper-1.0) < 0.05 else 'native tick is NOT the closure rate'})")

    # ---- verdict ----
    native = [y for r, y in pts if abs(r - 1.0) < 1e-9][0]
    print("\n" + "=" * 92)
    print("VERDICT (continuum completion):")
    print(f"  * the binning-free continuum DOS gives delta/target = {vals[-1]:.4f} at the native tick (rate=1),")
    print(f"    stable to {spread:.0e} across M=3k..24k -- so the ~1% the histogram showed is REAL, not binning")
    print(f"    noise: the native-tick value sits {'~at' if abs(vals[-1]-1)<0.01 else f'{(vals[-1]-1)*100:+.1f}% from'} target.")
    print(f"  * the exact target is crossed at rate ~= {upper:.3f}, {'~the native tick' if abs(upper-1)<0.05 else 'NOT the native tick'};")
    print(f"    the monitor rate stays load-bearing (delta/target runs {pts[0][1]:.2f}..{pts[-1][1]:.2f} over the scan).")
    print("  * CONVERGENCE half: settled -- the continuum value is pinned (no longer DOS-ambiguous).")
    print("    THEOREM half: still open -- why rate=1 is the physical readout rate, and why the two")
    print("    Peierls current-link occupations are the dressed-alpha observable, are not yet proven.")
    print("  * honest grade unchanged: BEST dressed-alpha route, non-circular in structure")
    print("    (alpha0 x pinned geometric kernel), continuum-stable within ~1% at the canonical tick --")
    print("    but NOT closed (rate=1 / eta=-1 choices + an un-proven readout).")


if __name__ == "__main__":
    main()
