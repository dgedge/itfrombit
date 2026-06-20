#!/usr/bin/env python3
r"""ITEM 123: external CAMB check of the CMB pressureless-slot target.

This is not a new derivation.  It is a Boltzmann-code target extraction for the
current Item-123/R4 CMB gate.

Canon already pins the analytic target:

    omega_b h^2       = 0.0224
    omega_nuR h^2     = 0.024   (20% sterile-nu_R cold anchor)
    missing omega_x   = 0.096   (the declared 80% R4 dark share)

The question checked here is whether adding omega_x as a genuinely cold,
pressureless, low-sound-speed Boltzmann component restores the acoustic third
peak, rather than merely fixing a background equality formula.

Implementation:
  * CAMB is optional and external.  If it is not importable, this script exits
    with a clear instruction instead of silently replacing the Boltzmann run.
  * The diagnostic sweep holds the observed acoustic scale fixed via
    cosmomc_theta and lets CAMB solve H0.  This isolates the peak-height issue
    from pure geometry.
  * omega_x is added as CDM-like matter.  That is exactly the phenomenological
    target the framework must derive: w=0, c_s^2=0, conserved rho~a^-3.

Expected durable result (CAMB 1.6.x, Planck-ish inputs):
  * omega_x=0 forces H0 ~ 147 to keep the acoustic scale and leaves D3/D2
    around 0.65, far below the restored case.
  * omega_x=0.096 gives omega_cold h^2=0.120, z_eq ~= 3409, H0 ~= 67.4,
    and D3/D2 ~= 0.98.

Verdict:
  The external Boltzmann gate confirms the analytic pressureless-slot target.
  Current active R4 line-current/exhaust branches cannot supply this component.
  A live framework route must derive an R4 service-phase / Stueckelberg /
  Brown-Kuchar dust variable with omega_x h^2 = 0.096, w=c_s^2=0, and no
  double-counting of late MOND/R4 or dark-energy ledgers.

Example invocation when CAMB is not installed globally:

    uv pip install --target /tmp/cmb_pkgs camb
    PYTHONPATH=/tmp/cmb_pkgs python3 python_code/item123_cmb_boltzmann_sweep.py
"""

from __future__ import annotations

from dataclasses import dataclass
import sys

try:
    import camb
    import numpy as np
    from scipy.signal import find_peaks
except Exception as exc:  # pragma: no cover - dependency gate
    raise SystemExit(
        "CAMB/SciPy dependency missing. Install externally, for example:\n"
        "  uv pip install --target /tmp/cmb_pkgs camb\n"
        "  PYTHONPATH=/tmp/cmb_pkgs python3 python_code/item123_cmb_boltzmann_sweep.py\n"
        f"Import failure: {type(exc).__name__}: {exc}"
    ) from exc


OMEGA_B_H2 = 0.0224
OMEGA_NUR_H2 = 0.024
OMEGA_X_TARGET_H2 = 0.096
OMEGA_COLD_TARGET_H2 = OMEGA_NUR_H2 + OMEGA_X_TARGET_H2

AS = 2.10e-9
NS = 27.0 / 28.0
N_EFF = 3.044
TAU = 0.054
Y_HE = 0.245
COSMOMC_THETA = 0.010411
LMAX = 2200

OMEGA_GAMMA_H2 = 2.469e-5
NEUTRINO_FACTOR = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
OMEGA_R_H2 = OMEGA_GAMMA_H2 * NEUTRINO_FACTOR


@dataclass(frozen=True)
class SweepRow:
    omega_x_h2: float
    omega_cold_h2: float
    h0_solved: float
    z_eq: float
    ell1: int
    d1: float
    ell2: int
    d2: float
    ell3: int
    d3: float

    @property
    def d3_d1(self) -> float:
        return self.d3 / self.d1

    @property
    def d3_d2(self) -> float:
        return self.d3 / self.d2


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def z_eq(omega_m_h2: float) -> float:
    return omega_m_h2 / OMEGA_R_H2 - 1.0


def camb_tt_spectrum(omega_cold_h2: float) -> tuple[np.ndarray, float]:
    pars = camb.CAMBparams()
    pars.set_cosmology(
        cosmomc_theta=COSMOMC_THETA,
        theta_H0_range=(1.0, 300.0),
        ombh2=OMEGA_B_H2,
        omch2=omega_cold_h2,
        omk=0.0,
        tau=TAU,
        YHe=Y_HE,
        nnu=N_EFF,
        mnu=0.06,
        num_massive_neutrinos=1,
    )
    pars.InitPower.set_params(As=AS, ns=NS)
    pars.set_for_lmax(LMAX, lens_potential_accuracy=0)
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK")
    return powers["total"][:, 0], float(results.Params.H0)


def first_three_acoustic_peaks(tt: np.ndarray) -> list[tuple[int, float]]:
    ell = np.arange(len(tt))
    mask = (ell >= 50) & (ell <= 1200)
    peak_indexes, _ = find_peaks(tt[mask], prominence=50.0, distance=120)
    peak_ells = ell[mask][peak_indexes]
    peak_heights = tt[mask][peak_indexes]
    order = np.argsort(peak_ells)
    peaks = [(int(l), float(h)) for l, h in zip(peak_ells[order], peak_heights[order])]
    if len(peaks) < 3:
        raise RuntimeError(f"CAMB spectrum exposed fewer than three acoustic peaks: {peaks}")
    return peaks[:3]


def run_row(omega_x_h2: float) -> SweepRow:
    omega_cold_h2 = OMEGA_NUR_H2 + omega_x_h2
    tt, h0 = camb_tt_spectrum(omega_cold_h2)
    (ell1, d1), (ell2, d2), (ell3, d3) = first_three_acoustic_peaks(tt)
    return SweepRow(
        omega_x_h2=omega_x_h2,
        omega_cold_h2=omega_cold_h2,
        h0_solved=h0,
        z_eq=z_eq(OMEGA_B_H2 + omega_cold_h2),
        ell1=ell1,
        d1=d1,
        ell2=ell2,
        d2=d2,
        ell3=ell3,
        d3=d3,
    )


def main() -> None:
    print("ITEM 123: CAMB BOLTZMANN SWEEP FOR THE PRESSURELESS CMB SLOT")
    print("=" * 88)
    print(f"CAMB version: {getattr(camb, '__version__', 'unknown')}")
    print(
        "Fixed inputs: "
        f"omega_b h^2={OMEGA_B_H2:.4f}, omega_nuR h^2={OMEGA_NUR_H2:.4f}, "
        f"n_s=27/28={NS:.6f}, N_eff={N_EFF}, tau={TAU}, theta={COSMOMC_THETA}"
    )
    print(
        "Sweep: omega_x is inserted as CDM-like matter.  This is the target "
        "a derived zero-mode must reproduce, not an assumed closure."
    )

    omega_x_values = [0.000, 0.024, 0.048, 0.072, 0.096, 0.120]
    rows = [run_row(x) for x in omega_x_values]
    target = next(r for r in rows if abs(r.omega_x_h2 - OMEGA_X_TARGET_H2) < 1.0e-12)
    base = next(r for r in rows if abs(r.omega_x_h2) < 1.0e-12)

    print("\n[1] Acoustic-scale-fixed CAMB sweep")
    print(
        "omega_x  om_cold  H0_solved  z_eq   "
        "ell1 D1      ell2 D2      ell3 D3      D3/D1 D3/D2"
    )
    for r in rows:
        print(
            f"{r.omega_x_h2:7.3f} {r.omega_cold_h2:8.3f} {r.h0_solved:9.2f} "
            f"{r.z_eq:6.0f} {r.ell1:5d} {r.d1:7.1f} {r.ell2:5d} {r.d2:7.1f} "
            f"{r.ell3:5d} {r.d3:7.1f} {r.d3_d1:6.3f} {r.d3_d2:6.3f}"
        )

    print("\n[2] Gate checks")
    check(abs(target.omega_cold_h2 - OMEGA_COLD_TARGET_H2) < 1.0e-12, "target row has omega_cold h^2 = 0.120")
    check(3400.0 < target.z_eq < 3420.0, "omega_x=0.096 restores z_eq ~= 3409")
    check(65.0 < target.h0_solved < 70.0, "omega_x=0.096 keeps a normal H0 at fixed acoustic scale")
    check(base.h0_solved > 120.0, "omega_x=0 requires pathological H0 at fixed acoustic scale")
    check(base.z_eq < 1200.0, "omega_x=0 leaves equality near recombination, not LCDM")
    check(base.d3_d2 / target.d3_d2 < 0.70, "omega_x=0 underproduces the relative third peak")
    check(0.90 < target.d3_d2 < 1.10, "omega_x=0.096 restores the third/second peak scale")

    print("\n[3] Verdict")
    print(
        "  CAMB confirms the budget-level gate: the cold sterile-nu_R anchor alone "
        "does not reproduce the acoustic pattern.  A pressureless, low-sound-speed "
        "component with omega_x h^2=0.096 restores equality and the relative third "
        "peak in the diagnostic Boltzmann system."
    )
    print(
        "  Therefore the framework must derive a genuine R4 service-phase / "
        "Stueckelberg / Brown-Kuchar dust variable with w=c_s^2=0 and the existing "
        "80% R4 abundance, or concede the CMB completion."
    )
    print("exit 0 -- external CAMB target extraction complete.")


if __name__ == "__main__":
    sys.exit(main())
