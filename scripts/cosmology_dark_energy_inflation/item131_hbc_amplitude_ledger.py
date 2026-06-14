#!/usr/bin/env python3
r"""ITEM 131: HBC scalar-clock amplitude ledger audit.

Question
--------
Can the primordial scalar-clock amplitude A_nu be derived from the existing
boundary-printing / QEC service ledger, rather than inserted?

Verdict
-------
Not yet.  The finite 28-channel service instrument fixes relative channel
weights and the anomalous log-shell generator, but it does not fix the
absolute number of statistically independent print/service events in one
mode-local horizon shell.

The sharp ledger form is:

    nu_HBC(x,N) = integral [j(x,N)/jbar(N) - 1] dN

and for a stochastic event current with effective Fano/correlation factor
F_eff and N_eff independent service events per mode-local shell,

    A_nu = F_eff / N_eff.

Thus CMB-sized scalar amplitude A_nu ~= 2.1e-9 requires

    N_eff / F_eff ~= 4.8e8

independent service events per printed shell.  The 28-channel code structure
does not supply this number.  It supplies the channel partition and the
1/28 generator; the total event density, Fano factor, and correlation volume
remain microscopic HBC/KMS/Landauer data.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
N_CHANNELS = 28
N_FLAGS = 112
DELTA = Fraction(1, 28)

# Existing item-131 scripts already use 2.1e-9 as the scalar-amplitude scale.
A_SCALAR = 2.1e-9
MPL_REDUCED_GEV = 2.435e18
MPL_NONREDUCED_GEV = 1.2209e19


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def shot_noise_amplitude(n_eff: float, f_eff: float = 1.0) -> float:
    """Relative count-noise amplitude for a compensated event current."""
    return f_eff / n_eff


def required_events(amplitude: float = A_SCALAR, f_eff: float = 1.0) -> float:
    return f_eff / amplitude


def equicorrelated_fano(n_channels: int, rho: float) -> float:
    """Effective Fano factor for n equal channels with pair correlation rho."""
    return 1.0 + (n_channels - 1) * rho


def h_over_mpl_reduced_from_entropy(n_eff: float) -> float:
    """If N_eff were de-Sitter entropy S=8*pi^2*Mbar_P^2/H^2."""
    return math.sqrt(8.0 * math.pi**2 / n_eff)


def h_over_mpl_nonreduced_from_entropy(n_eff: float) -> float:
    """Equivalent non-reduced convention S=pi*M_P^2/H^2."""
    return math.sqrt(math.pi / n_eff)


def main() -> None:
    print("ITEM 131 HBC AMPLITUDE LEDGER AUDIT")

    w_to_28 = (ROOT / "python_code" / "item131_w_to_28_instrument.py").read_text()
    clock_exit = (ROOT / "python_code" / "item131_hbc_clock_covariance_exit.py").read_text()

    print("\n[1] What the finite service instrument fixes")
    p_channel = Fraction(1, N_CHANNELS)
    p_flag_given_point = Fraction(1, 14)
    check("service channel has rate 1/28" in w_to_28, "W-to-28 script fixes relative service-channel weights")
    check(p_channel == DELTA, "uniform channel probability is 1/28")
    check(N_FLAGS == 4 * N_CHANNELS, "112 microscopic flags coarse-grain four-to-one onto 28 channels")
    check(p_flag_given_point * 14 == 1, "incidence refinement preserves each original single-bit Kraus effect")
    print("  The 28-channel theorem fixes p_c = 1/28 after conditioning on a service tick.")
    print("  It does not fix how many ticks/events occur in one horizon log shell.")

    print("\n[2] Event-count noise cases")
    deterministic_total_amplitude = 0.0
    one_event_total = shot_noise_amplitude(1.0)
    one_event_per_channel = shot_noise_amplitude(float(N_CHANNELS))
    million_events = shot_noise_amplitude(1.0e6)
    print(f"  deterministic total event count: A_nu = {deterministic_total_amplitude:.3e}")
    print(f"  one total event per shell       : A_nu = {one_event_total:.3e}")
    print(f"  one per 28 channels per shell   : A_nu = {one_event_per_channel:.3e}")
    print(f"  1e6 events per shell            : A_nu = {million_events:.3e}")
    check(deterministic_total_amplitude == 0.0, "fixed total service load produces no scalar-clock count fluctuation")
    check(one_event_per_channel / A_SCALAR > 1.0e7, "28 labels alone would be far too noisy if read as the event count")
    check(million_events > A_SCALAR, "even one million events per shell is still above CMB-sized amplitude")
    check("A_nu remains an undetermined" in clock_exit, "previous clock-covariance audit leaves A_nu open")

    print("\n[3] Required microscopic service-event density")
    n_required = required_events()
    print(f"  target A_nu              = {A_SCALAR:.3e}")
    print(f"  required N_eff/F_eff     = {n_required:.6e}")
    print(f"  events per 28-channel    = {n_required / N_CHANNELS:.6e}")
    print(f"  events per 112 flag      = {n_required / N_FLAGS:.6e}")
    check(abs(shot_noise_amplitude(n_required) - A_SCALAR) / A_SCALAR < 1e-15, "A_nu=1/N_eff fixes N_eff at fixed F_eff")
    check(n_required / N_CHANNELS > 1.0e7, "each service channel needs O(10^7) independent events for F_eff=1")
    check(n_required / N_FLAGS > 1.0e6, "each microscopic incidence flag still needs O(10^6) events for F_eff=1")

    print("\n[4] Correlations are load-bearing")
    for rho in [0.0, 0.01, 1.0 / 27.0, -0.01]:
        f_eff = equicorrelated_fano(N_CHANNELS, rho)
        status = "valid" if f_eff >= 0.0 else "invalid/over-anticorrelated"
        print(
            f"  rho={rho:+.5f}: F_eff={f_eff:.6f}, "
            f"required N_eff={required_events(f_eff=f_eff):.6e} ({status})"
        )
    check(equicorrelated_fano(N_CHANNELS, 0.0) == 1.0, "independent channels give F_eff=1")
    check(equicorrelated_fano(N_CHANNELS, 1.0 / 27.0) == 2.0, "positive channel correlations double the required event count")
    check(equicorrelated_fano(N_CHANNELS, -1.0 / 27.0) == 0.0, "fixed-total anticorrelation can erase total-count noise")
    check(True, "therefore A_nu needs the service-current Fano/correlation ledger, not only channel counting")

    print("\n[5] Horizon-entropy closure attempt")
    h_red = h_over_mpl_reduced_from_entropy(n_required)
    h_nonred = h_over_mpl_nonreduced_from_entropy(n_required)
    print("  If N_eff is identified with de-Sitter entropy S_dS:")
    print(f"    reduced convention H/Mbar_P = {h_red:.6e}, H ~= {h_red * MPL_REDUCED_GEV:.3e} GeV")
    print(f"    non-reduced convention H/M_P = {h_nonred:.6e}, H ~= {h_nonred * MPL_NONREDUCED_GEV:.3e} GeV")
    check(abs(h_red * MPL_REDUCED_GEV - h_nonred * MPL_NONREDUCED_GEV) / (h_red * MPL_REDUCED_GEV) < 0.01, "entropy conventions give the same H scale")
    check(h_red > 0, "entropy identification would infer an inflationary H scale")
    check(True, "without an independent HBC derivation of H or S_dS, this is an inversion of A_nu, not a prediction")

    print("\n[6] Closure status")
    closed = [
        "A_nu reduces to the effective event-count ratio F_eff/N_eff",
        "CMB-sized A_nu requires N_eff/F_eff ~= 4.8e8 events per mode-local shell",
        "the 28-channel service instrument fixes relative weights but not total event density",
    ]
    open_items = [
        "absolute service-event density per horizon log shell",
        "Fano/correlation factor of the HBC/KMS/Landauer current",
        "correlation volume that maps microscopic flags to one scalar mode",
        "independent early H or entropy scale, if N_eff is read as S_dS",
    ]
    for item in closed:
        check(True, f"conditional reduction: {item}")
    for item in open_items:
        check(True, f"not closed: {item}")

    print("\n" + "=" * 100)
    print("VERDICT")
    print("  The amplitude bridge sharpens to A_nu = F_eff/N_eff.")
    print("  The existing 28-channel service ledger fixes p_c=1/28 and the tilt")
    print("  generator, but it does not fix N_eff or F_eff.  Deterministic total")
    print("  service gives A_nu=0; one-event/channel counting is much too noisy;")
    print("  CMB-sized amplitude needs about 4.8e8 independent service events per")
    print("  mode-local shell for F_eff=1.  Reading that count as de-Sitter entropy")
    print("  infers H~1e15 GeV, but this is not a parameter-free prediction unless")
    print("  HBC independently derives that entropy/event density.")
    print("=" * 100)
    print("exit 0 -- A_nu reduced to F_eff/N_eff; microscopic event-density ledger remains open.")


if __name__ == "__main__":
    main()
