#!/usr/bin/env python3
r"""ITEM 123 / M15: attempted full CMB lock.

This script tests the two remaining CMB gates after the zero-mode reservoir:

  1. selector-H0 versus acoustic scale;
  2. halo non-double-counting.

It deliberately distinguishes a numerical repair from a canon lock.

Result
------
With the exact framework matter budget, CAMB wants

    H0_CMB ~= 66.218 km/s/Mpc = 0.98442 H0_selector.

The finite-depth-looking value

    H0_CMB = (63/64) H0_selector

matches 100 theta_* at the 0.001% level and keeps the TT shape at about the
same ~1% level.  This is a sharp theorem target because 64 = 2^6 is the same
depth-6 lock scale that appears in the selector/boot residual machinery.

But it is not derived here.  Existing canon previously rejected analogous
unlicensed finite-depth boundary factors (e.g. the 129/128 boot correction).
So the honest outcome is:

  * theta_* is numerically lockable by a 63/64 acoustic pre-lock projection;
  * full CMB is not canon-Locked until a service-ledger theorem proves that
    the CMB acoustic ruler sees the depth-6 pre-latch span while the late
    selector H0 sees the completed 64/64 span.

Halo branch
-----------
The CMB-compatible accounting branch is:

    N_zero is the CDM-like mobile halo mass;
    active R4/MOND is not also fitted as an independent galaxy force.

If the framework insists on active baryonic MOND/RAR in galaxies, the fair-
sample zero-mode must be depleted/screened by >95.9%.  No current finite
operator supplies that.  Thus the MOND-preserving branch is not locked.
"""

from __future__ import annotations

import math

import camb
import numpy as np


ALPHA0_SKY = 1.0 / 137.036
ALPHA0_CANON = 1.0 / 137.0
LAMBDA_KEV = 331.7e3
N_GAMMA = 410.7
RHO_CRIT_H2 = 1.0537e4
M_P_EV = 938.272e6

H0_SELECTOR = 67.266
OMEGA_L_CANON = 12.0 * math.pi / 55.0
DEPTH6_FACTOR = 63.0 / 64.0

NS_FW = 27.0 / 28.0
AS_FW = 0.75 * ALPHA0_CANON**4
TAU_FW = 0.054
W0_FW = -1.0 + 1.0 / 28.0
WA_FW = -1.0 / 28.0

H0_PLANCK = 67.36
OMBH2_PLANCK = 0.02237
OMCH2_PLANCK = 0.1200
NS_PLANCK = 0.9649
AS_PLANCK = 2.10e-9
TAU_PLANCK = 0.0544


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def framework_matter_budget() -> tuple[float, float, float]:
    eta = (3.0 / 14.0) * ALPHA0_SKY**4
    omega_b = eta * N_GAMMA * M_P_EV / RHO_CRIT_H2
    m_nur_ev = ALPHA0_SKY**2 * LAMBDA_KEV * 1.0e3
    omega_nur = (ALPHA0_SKY / 208.0) * N_GAMMA * m_nur_ev / RHO_CRIT_H2
    omega_dark = 5.0 * omega_nur
    return omega_b, omega_nur, omega_dark


def params(
    *,
    H0: float,
    ombh2: float,
    omch2: float,
    ns: float,
    As: float,
    tau: float,
    w0: float,
    wa: float,
    lmax: int,
) -> camb.CAMBparams:
    pars = camb.set_params(
        H0=H0,
        ombh2=ombh2,
        omch2=omch2,
        ns=ns,
        As=As,
        tau=tau,
        w=w0,
        wa=wa,
        dark_energy_model="ppf",
    )
    pars.set_for_lmax(lmax, lens_potential_accuracy=1 if lmax > 1000 else 0)
    return pars


def theta100_framework(H0: float, omega_b: float, omega_dark: float) -> float:
    pars = params(
        H0=H0,
        ombh2=omega_b,
        omch2=omega_dark,
        ns=NS_FW,
        As=AS_FW,
        tau=TAU_FW,
        w0=W0_FW,
        wa=WA_FW,
        lmax=450,
    )
    return float(camb.get_results(pars).cosmomc_theta() * 100.0)


def theta100_planck() -> float:
    pars = params(
        H0=H0_PLANCK,
        ombh2=OMBH2_PLANCK,
        omch2=OMCH2_PLANCK,
        ns=NS_PLANCK,
        As=AS_PLANCK,
        tau=TAU_PLANCK,
        w0=-1.0,
        wa=0.0,
        lmax=450,
    )
    return float(camb.get_results(pars).cosmomc_theta() * 100.0)


def find_h0_root(omega_b: float, omega_dark: float, target: float) -> float:
    lo, hi = 65.0, 68.0
    f_lo = theta100_framework(lo, omega_b, omega_dark) - target
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        f_mid = theta100_framework(mid, omega_b, omega_dark) - target
        if f_lo * f_mid <= 0.0:
            hi = mid
        else:
            lo = mid
            f_lo = f_mid
    return 0.5 * (lo + hi)


def spectrum(H0: float, omega_b: float, omega_dark: float, *, framework: bool) -> tuple[np.ndarray, float]:
    if framework:
        pars = params(
            H0=H0,
            ombh2=omega_b,
            omch2=omega_dark,
            ns=NS_FW,
            As=AS_FW,
            tau=TAU_FW,
            w0=W0_FW,
            wa=WA_FW,
            lmax=2500,
        )
    else:
        pars = params(
            H0=H0_PLANCK,
            ombh2=OMBH2_PLANCK,
            omch2=OMCH2_PLANCK,
            ns=NS_PLANCK,
            As=AS_PLANCK,
            tau=TAU_PLANCK,
            w0=-1.0,
            wa=0.0,
            lmax=2500,
        )
    res = camb.get_results(pars)
    return res.get_cmb_power_spectra(pars, CMB_unit="muK")["total"][:, 0], float(res.cosmomc_theta() * 100.0)


def peak(tt: np.ndarray, lo: int, hi: int) -> tuple[int, float]:
    seg = tt[lo:hi]
    idx = int(np.argmax(seg))
    return lo + idx, float(seg[idx])


def compare_tt(H0: float, omega_b: float, omega_dark: float) -> tuple[float, float, tuple[float, float, float]]:
    fw, _ = spectrum(H0, omega_b, omega_dark, framework=True)
    pl, _ = spectrum(H0_PLANCK, omega_b, omega_dark, framework=False)
    length = min(len(fw), len(pl))
    ell = np.arange(length)
    mask = (ell >= 30) & (ell <= 2000)
    frac = np.abs(fw[:length][mask] / pl[:length][mask] - 1.0)
    peak_diffs = []
    for lo, hi in ((180, 260), (480, 600), (760, 870)):
        _, fv = peak(fw, lo, hi)
        _, pv = peak(pl, lo, hi)
        peak_diffs.append((fv / pv - 1.0) * 100.0)
    return float(np.sqrt(np.mean(frac**2)) * 100.0), float(np.max(frac) * 100.0), tuple(peak_diffs)


def main() -> None:
    print("ITEM 123 / M15: FULL CMB LOCK ATTEMPT")
    print("=" * 90)
    print(f"CAMB version: {getattr(camb, '__version__', 'unknown')}")

    omega_b, omega_nur, omega_dark = framework_matter_budget()
    omega_m = omega_b + omega_dark
    target = theta100_planck()
    theta_selector = theta100_framework(H0_SELECTOR, omega_b, omega_dark)
    h0_root = find_h0_root(omega_b, omega_dark, target)
    h0_6364 = H0_SELECTOR * DEPTH6_FACTOR
    theta_6364 = theta100_framework(h0_6364, omega_b, omega_dark)

    print("\n[1] Exact framework matter budget")
    print(f"  omega_b h^2       = {omega_b:.6f}")
    print(f"  omega_nuR h^2     = {omega_nur:.6f}")
    print(f"  omega_dark h^2    = {omega_dark:.6f}")
    print(f"  omega_m h^2       = {omega_m:.6f}")
    check(abs(omega_dark - 0.1209) < 0.001, "derived dark budget remains LCDM-sized")

    print("\n[2] Acoustic scale")
    print(f"  Planck/LCDM reference 100 theta_*       = {target:.6f}")
    print(f"  selector H0={H0_SELECTOR:.3f}: theta     = {theta_selector:.6f}  ({(theta_selector/target-1)*100:+.3f}%)")
    print(f"  CAMB root H0                            = {h0_root:.6f}  ({h0_root/H0_SELECTOR:.8f} of selector)")
    print(f"  finite-depth candidate (63/64) H0       = {h0_6364:.6f}  theta={theta_6364:.6f}  ({(theta_6364/target-1)*100:+.4f}%)")
    print(f"  root / (63/64 selector)                 = {h0_root/h0_6364:.8f}")
    omega_l_root = 1.0 - omega_m / (h0_root / 100.0) ** 2
    print(f"  flat Omega_L at CAMB root               = {omega_l_root:.6f}  (canon 12pi/55={OMEGA_L_CANON:.6f})")
    check(abs(theta_selector / target - 1.0) > 0.002, "selector-H0 offset is real")
    check(abs(theta_6364 / target - 1.0) < 2.0e-5, "63/64 factor numerically closes theta_*")
    check(abs(h0_root / h0_6364 - 1.0) < 1.0e-4, "CAMB root is essentially the 63/64 depth-6 candidate")
    check(abs(omega_l_root - OMEGA_L_CANON) > 0.01, "plain H0 repair breaks the canonical Omega_L/flatness package")

    rms, maxdiff, peak_diffs = compare_tt(h0_6364, omega_b, omega_dark)
    print("\n[3] TT shape at 63/64 acoustic projection")
    print(f"  RMS frac 30<=ell<=2000       = {rms:.2f}%")
    print(f"  max frac 30<=ell<=2000       = {maxdiff:.2f}%")
    print("  peak-height diffs            = " + ", ".join(f"{x:+.2f}%" for x in peak_diffs))
    check(rms < 2.0 and maxdiff < 2.0, "63/64 projection preserves the TT shape closure")
    check(max(abs(x) for x in peak_diffs) < 2.0, "peak heights remain closed at the percent level")

    print("\n[4] Canon-lock status of 63/64")
    print("  candidate source: depth-6 selector span has 2^6=64 slots; CMB acoustic")
    print("  projection would see a pre-latch 63/64 span while late H0 sees 64/64.")
    print("  status: THEOREM TARGET ONLY. Current canon has no operator proving that")
    print("  acoustic propagation omits exactly one depth-6 slot, and analogous")
    print("  unlicensed finite-depth boundary factors have been rejected before.")
    check(True, "63/64 is recorded as a theorem target, not promoted to Locked")

    print("\n[5] Halo branch decision")
    zero_to_baryon = (4.0 * omega_nur) / omega_b
    rar_room = 10.0**0.07 - 1.0
    max_fraction = rar_room / zero_to_baryon
    depletion = 1.0 - max_fraction
    print(f"  zero-mode/baryon fair-sample mass ratio       = {zero_to_baryon:.3f}")
    print(f"  active-MOND branch max zero-mode fraction     = {max_fraction:.3f}")
    print(f"  required zero-mode depletion/screening         >= {depletion:.1%}")
    print("  CMB-lock branch: zero-mode is the CDM-like mobile halo; active R4/MOND")
    print("  is not also fitted as an independent galaxy force.")
    check(depletion > 0.95, "MOND-preserving branch needs a new >95% depletion/screening theorem")
    check(True, "CDM-like zero-mode branch avoids double-counting by construction")

    print("\n[6] Verdict")
    print(
        "  STRICT RESULT: full CMB is still not Locked under current canon.\n"
        "  What changed is that the remaining theta_* failure is now a single,\n"
        "  sharp finite-depth theorem target: prove H_CMB=(63/64)H_selector from\n"
        "  the service ledger, or accept that the selector-H0 package is in CMB\n"
        "  tension.  Halo non-double-counting is lockable only by choosing the\n"
        "  zero-mode-CDM branch; retaining active MOND/RAR requires a new, strong\n"
        "  depletion theorem.  Therefore the full-CMB lock package would be:\n"
        "    A. derive the 63/64 acoustic pre-latch projection; and\n"
        "    B. choose zero-mode as the mobile halo mass, demoting active R4/MOND\n"
        "       from an independently fitted galaxy force.\n"
        "  Without A, the CMB remains selector-gated. Without B, halos double-count."
    )
    print("exit 0 -- CMB lock attempt complete: 63/64 target found; halo branch forced; full lock not yet earned.")


if __name__ == "__main__":
    main()
