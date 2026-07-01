#!/usr/bin/env python3
r"""Item 123 / M15: acoustic-clock theorem lifted to a native Boltzmann spectrum.

Question
--------
The local theorem now gives

    P_acoustic = P_busy = I - |111111><111111|,
    H_CMB / H_selector = Tr(P_busy)/64 = 63/64.

Does that stay merely a one-number theta_* repair, or does it lift to the full
linear Boltzmann spectrum as the framework's native CMB branch?

Result
------
At forward Boltzmann-code grade, it lifts.  The 63/64 acoustic-lapse branch is
numerically identical, to the relevant accuracy, to the direct CAMB theta-root
branch, and the same branch passes the existing pre-likelihood TT/EE/TE/lensing
and matter-power gate.  The sterile-only control still fails, so the zero-mode
reservoir remains load-bearing.

This is not a Planck plik/CamSpec/cobaya likelihood and not a nonlinear halo
catalogue.  It closes the "conditional acoustic clock -> native spectrum"
step: once the QND busy-lapse theorem is accepted, the native spectrum is the
63/64 Boltzmann branch, not a fitted H0 shift.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock
import item123_external_cmb_lensing_nonlinear_gate as external_gate


@dataclass(frozen=True)
class Anchor:
    path: Path
    patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text(encoding="utf-8")
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def anchors() -> list[Anchor]:
    return [
        Anchor(
            ROOT / "python_code" / "item123_substrate_busy_flag_law.py",
            (
                r"A_busy = I - \|111111><111111\|",
                r"every incomplete selector state exposes one",
                r"Tr\(P_busy\)/64",
            ),
            "Finite selector theorem: active-address demux gives the depth-6 busy projector.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_acoustic_clock_qnd_lapse_readout_audit.py",
            (
                r"acoustic phase is a clock",
                r"P_acoustic = P_busy",
                r"H_CMB / H_selector = Tr\(P_acoustic rho_selector\) = 63/64",
            ),
            "Sector theorem: the acoustic ruler is the QND service-lapse readout.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_external_cmb_lensing_nonlinear_gate.py",
            (
                r"TT/EE/TE/lensing spectra",
                r"framework 63/64",
                r"sterile-only control fails",
            ),
            "External gate: the 63/64 branch is checked against full CAMB spectra and controls.",
        ),
    ]


def fractional_rms(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a[mask] / b[mask] - 1.0) ** 2)))


def normalized_rms(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2)) / np.sqrt(np.mean(b[mask] ** 2)))


def branch_to_branch_metrics(data, ref) -> dict[str, float]:
    length = min(len(data.ell), len(ref.ell))
    ell = np.arange(length)
    acoustic = (ell >= 30) & (ell <= 2000)
    lens_len = min(len(data.lens_potential), len(ref.lens_potential))
    lens_ell = np.arange(lens_len)
    lens_mask = (lens_ell >= 10) & (lens_ell <= 1000) & (ref.lens_potential[:lens_len, 0] > 0)

    return {
        "theta_delta": data.theta100 / ref.theta100 - 1.0,
        "TT_rms": fractional_rms(data.total[:length, 0], ref.total[:length, 0], acoustic),
        "EE_rms": fractional_rms(data.total[:length, 1], ref.total[:length, 1], acoustic & (ref.total[:length, 1] > 1.0e-6)),
        "TE_norm": normalized_rms(data.total[:length, 3], ref.total[:length, 3], acoustic),
        "lens_rms": fractional_rms(data.lens_potential[:lens_len, 0], ref.lens_potential[:lens_len, 0], lens_mask),
    }


def pk_branch_metrics(data, ref, *, nonlinear: bool) -> dict[str, float]:
    kh = data.kh_nl if nonlinear else data.kh_lin
    pk = data.pk_nl if nonlinear else data.pk_lin
    rpk = ref.pk_nl if nonlinear else ref.pk_lin
    bands = {
        "lin_0.01_0.3": (kh >= 0.01) & (kh <= 0.3),
        "quasi_0.3_1": (kh > 0.3) & (kh <= 1.0),
        "nonlin_1_5": (kh > 1.0) & (kh <= 5.0),
    }
    out: dict[str, float] = {}
    for name, mask in bands.items():
        diffs = [fractional_rms(pk[iz], rpk[iz], mask) for iz in range(len(external_gate.REDSHIFTS))]
        out[name] = max(diffs)
    return out


def percent(value: float) -> str:
    return f"{100.0 * value:.4f}%"


def main() -> None:
    print("ITEM 123 / M15: CMB NATIVE BOLTZMANN-SPECTRUM LIFT")
    print("=" * 96)
    print(f"CAMB version: {getattr(lock.camb, '__version__', 'unknown')}")

    print("\n[1] Required theorem anchors")
    for anchor in anchors():
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(source_contains(anchor.path, anchor.patterns), f"{rel} contains the required clause")

    omega_b, omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_native = lock.H0_SELECTOR * lock.DEPTH6_FACTOR
    theta_native = lock.theta100_framework(h0_native, omega_b, omega_dark)
    theta_root = lock.theta100_framework(h0_root, omega_b, omega_dark)

    print("\n[2] Native acoustic clock versus direct CAMB theta root")
    print(f"  selector H0                     = {lock.H0_SELECTOR:.6f}")
    print(f"  native 63/64 H0                 = {h0_native:.6f}")
    print(f"  direct CAMB theta-root H0       = {h0_root:.6f}")
    print(f"  root/native ratio               = {h0_root / h0_native:.8f}")
    print(f"  native 100 theta_*              = {theta_native:.9f}")
    print(f"  root 100 theta_*                = {theta_root:.9f}")
    print(f"  Planck-comparator 100 theta_*   = {target:.9f}")
    check(abs(h0_root / h0_native - 1.0) < 1.0e-4, "63/64 native clock is the CAMB theta root")
    check(abs(theta_native / target - 1.0) < 2.0e-5, "native clock closes the acoustic scale")

    print("\n[3] Full Boltzmann spectra")
    specs = external_gate.models()
    results = {spec.name: external_gate.run(spec) for spec in specs}
    planck = results["Planck LCDM comparator"]
    native = results["framework 63/64"]
    root = results["framework theta-root"]
    selector = results["framework selector-H0"]
    sterile = results["sterile-only control"]

    native_vs_root = branch_to_branch_metrics(native, root)
    native_vs_planck = external_gate.spectrum_metrics(native, planck)
    selector_vs_planck = external_gate.spectrum_metrics(selector, planck)
    sterile_vs_planck = external_gate.spectrum_metrics(sterile, planck)

    print("  native-vs-root CMB metric       value")
    for key, value in native_vs_root.items():
        print(f"  {key:30s} {percent(value)}")
    check(abs(native_vs_root["theta_delta"]) < 1.0e-5, "native 63/64 and CAMB root have the same acoustic scale")
    check(native_vs_root["TT_rms"] < 5.0e-4, "native 63/64 and CAMB root are TT-identical at forward-spectrum grade")
    check(native_vs_root["EE_rms"] < 5.0e-4, "native 63/64 and CAMB root are EE-identical at forward-spectrum grade")
    check(native_vs_root["TE_norm"] < 5.0e-4, "native 63/64 and CAMB root are TE-identical at forward-spectrum grade")
    check(native_vs_root["lens_rms"] < 5.0e-4, "native 63/64 and CAMB root are lensing-identical at forward-spectrum grade")

    print("\n  native-vs-Planck pre-likelihood metric")
    for key in ("TT_rms", "EE_rms", "TE_norm", "lens_rms", "theta_delta", "sigma8_delta", "s8_delta"):
        print(f"  {key:30s} {percent(native_vs_planck[key])}")
    check(native_vs_planck["TT_rms"] < 0.025, "native branch passes the lensed TT pre-likelihood gate")
    check(native_vs_planck["EE_rms"] < 0.04, "native branch passes the lensed EE pre-likelihood gate")
    check(native_vs_planck["TE_norm"] < 0.06, "native branch passes the TE pre-likelihood gate")
    check(native_vs_planck["lens_rms"] < 0.08, "native branch passes the lensing-potential pre-likelihood gate")
    check(abs(native_vs_planck["theta_delta"]) < 5.0e-5, "native branch removes the selector theta_* failure")
    check(abs(selector_vs_planck["theta_delta"]) > 2.0e-3, "late selector branch still fails theta_*")
    check(sterile_vs_planck["TT_rms"] > 0.10, "sterile-only control still fails the CMB spectra")

    print("\n[4] Matter-power lift")
    native_root_lin = pk_branch_metrics(native, root, nonlinear=False)
    native_root_nl = pk_branch_metrics(native, root, nonlinear=True)
    native_planck_lin = external_gate.pk_metrics(native, planck, nonlinear=False)
    native_planck_nl = external_gate.pk_metrics(native, planck, nonlinear=True)
    print("  native-vs-root max RMS over redshift bands")
    for label, metrics in (("linear", native_root_lin), ("halofit", native_root_nl)):
        for band, value in metrics.items():
            print(f"  {label:8s} {band:14s} {percent(value)}")
            check(value < 5.0e-4, f"native/root {label} {band} is identical at forward-spectrum grade")

    print("\n  native-vs-Planck matter-power gate")
    for band, values in native_planck_lin.items():
        print(f"  linear  {band:14s} max={percent(max(values))}")
    for band, values in native_planck_nl.items():
        print(f"  halofit {band:14s} max={percent(max(values))}")
    check(max(native_planck_lin["lin_0.01_0.3"]) < 0.03, "native branch keeps linear P(k) close on structure scales")
    check(max(native_planck_nl["lin_0.01_0.3"]) < 0.04, "native branch keeps Halofit close on linear/quasi-linear scales")
    check(max(native_planck_nl["quasi_0.3_1"]) < 0.06, "native branch keeps Halofit controlled at k=0.3..1 h/Mpc")

    print("\n[5] Lift theorem status")
    print(
        """
  The local acoustic-clock theorem now composes with the Boltzmann solver:

      finite active-demux busy projector
        -> QND acoustic lapse
        -> H_CMB = (63/64) H_selector
        -> the full native CAMB branch.

  The 63/64 branch is not just a theta_* numerator trick.  Across TT, EE, TE,
  lensing potential, and linear/Halofit matter power it is numerically the
  same branch as the direct CAMB theta root, while the sterile-only control
  fails badly.  Therefore the CMB acoustic-clock conditional theorem is lifted
  to a native forward Boltzmann spectrum.

  Remaining boundaries:
    1. This is still a forward-spectrum / pre-likelihood result, not a full
       Planck TTTEEE+lensing likelihood with nuisance parameters.
    2. It assumes the zero-mode reservoir is the CDM-like mobile halo branch;
       nonlinear galaxy/halo catalogues remain external work.
    3. An exact finite Boltzmann-Keldysh derivation would be a deeper
       microscopic proof, but it is no longer the blocker for the native
       linear spectrum.
exit 0 -- acoustic-clock theorem lifted to native full Boltzmann-spectrum grade; full likelihood/halo modelling remain.
"""
    )


if __name__ == "__main__":
    main()
