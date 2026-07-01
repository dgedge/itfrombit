#!/usr/bin/env python3
r"""ITEM 123: zero-mode Brown--Kuchar dust in a Boltzmann/structure run.

Purpose
-------
The dark source map is now conditional-theorem grade:

    n_nuR/n_gamma = alpha0/208,
    omega_zero    = 4 omega_nuR,
    omega_dark    = omega_nuR + omega_zero ~= 0.1209.

This script performs the remaining external implementation gate at linear
Boltzmann/structure level:

  * implement the zero-mode reservoir as Brown--Kuchar dust, i.e. ordinary CDM
    in CAMB: w = c_s^2 = 0 and rho proportional to a^-3;
  * compare the z=0 linear matter power and growth diagnostics with a
    Planck-2018 LCDM reference;
  * keep the halo branch discipline fixed: N_zero is the CDM-like mobile halo
    mass, while active R4/MOND is not also fitted as an independent force;
  * classify the acoustic-scale alternatives: selector H0, the 63/64 pre-latch
    H0, and the direct CAMB theta_* root.

This is not a nonlinear halo simulation and not a likelihood analysis.  It is
the durable "does the derived dust component behave like CDM in a Boltzmann
code?" gate.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import sys

try:
    import camb
    import numpy as np
except Exception as exc:  # pragma: no cover - dependency gate
    raise SystemExit(
        "CAMB/NumPy dependency missing. Use the project Python environment, e.g.\n"
        "  ~/bin/py13_7/bin/python python_code/item123_zero_mode_structure_implementation.py\n"
        f"Import failure: {type(exc).__name__}: {exc}"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock
from item123_halo_branch_fork_resolution import required_depletion_fraction


K_MIN = 1.0e-3
K_MAX = 1.0
K_POINTS = 160
K_LINEAR_MIN = 1.0e-2
K_LINEAR_MAX = 3.0e-1
K_NORM = 5.0e-2


@dataclass(frozen=True)
class Model:
    name: str
    H0: float
    ombh2: float
    omch2: float
    ns: float
    As: float
    tau: float
    w0: float
    wa: float
    note: str


@dataclass(frozen=True)
class ModelResult:
    model: Model
    theta100: float
    omega_m: float
    sigma8: float
    s8: float
    fsigma8: float
    kh: np.ndarray
    pk: np.ndarray


@dataclass(frozen=True)
class Comparison:
    rms_percent: float
    max_percent: float
    shape_rms_percent: float
    sigma8_delta_percent: float
    s8_delta_percent: float
    theta_delta_percent: float


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def camb_params(model: Model) -> camb.CAMBparams:
    pars = camb.set_params(
        H0=model.H0,
        ombh2=model.ombh2,
        omch2=model.omch2,
        ns=model.ns,
        As=model.As,
        tau=model.tau,
        w=model.w0,
        wa=model.wa,
        dark_energy_model="ppf",
    )
    pars.set_for_lmax(900, lens_potential_accuracy=0)
    pars.set_matter_power(redshifts=[0.0], kmax=2.0)
    return pars


def run_model(model: Model) -> ModelResult:
    pars = camb_params(model)
    results = camb.get_results(pars)
    kh, _z, pk = results.get_matter_power_spectrum(
        minkh=K_MIN,
        maxkh=K_MAX,
        npoints=K_POINTS,
    )
    sigma8 = float(results.get_sigma8()[0])
    omega_m = (model.ombh2 + model.omch2) / (model.H0 / 100.0) ** 2
    return ModelResult(
        model=model,
        theta100=float(results.cosmomc_theta() * 100.0),
        omega_m=omega_m,
        sigma8=sigma8,
        s8=sigma8 * math.sqrt(omega_m / 0.3),
        fsigma8=float(results.get_fsigma8()[0]),
        kh=kh,
        pk=pk[0],
    )


def compare_to_reference(result: ModelResult, reference: ModelResult) -> Comparison:
    mask = (result.kh >= K_LINEAR_MIN) & (result.kh <= K_LINEAR_MAX)
    frac = result.pk[mask] / reference.pk[mask] - 1.0
    norm_index = int(np.argmin(np.abs(result.kh - K_NORM)))
    shape = (result.pk / result.pk[norm_index]) / (reference.pk / reference.pk[norm_index]) - 1.0
    shape_frac = shape[mask]
    return Comparison(
        rms_percent=float(np.sqrt(np.mean(frac**2)) * 100.0),
        max_percent=float(np.max(np.abs(frac)) * 100.0),
        shape_rms_percent=float(np.sqrt(np.mean(shape_frac**2)) * 100.0),
        sigma8_delta_percent=(result.sigma8 / reference.sigma8 - 1.0) * 100.0,
        s8_delta_percent=(result.s8 / reference.s8 - 1.0) * 100.0,
        theta_delta_percent=(result.theta100 / reference.theta100 - 1.0) * 100.0,
    )


def models() -> tuple[list[Model], float]:
    omega_b, omega_nur, omega_dark = lock.framework_matter_budget()
    target_theta = lock.theta100_planck()
    h0_root = lock.find_h0_root(omega_b, omega_dark, target_theta)
    h0_6364 = lock.H0_SELECTOR * lock.DEPTH6_FACTOR

    return [
        Model(
            "Planck LCDM reference",
            lock.H0_PLANCK,
            lock.OMBH2_PLANCK,
            lock.OMCH2_PLANCK,
            lock.NS_PLANCK,
            lock.AS_PLANCK,
            lock.TAU_PLANCK,
            -1.0,
            0.0,
            "external Planck-2018 LCDM comparator",
        ),
        Model(
            "framework selector-H0",
            lock.H0_SELECTOR,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "zero-mode+nu_R as CDM; late selector H0",
        ),
        Model(
            "framework 63/64 acoustic",
            h0_6364,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "zero-mode+nu_R as CDM; pre-latch H_CMB=(63/64)H_selector",
        ),
        Model(
            "framework theta-root",
            h0_root,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "zero-mode+nu_R as CDM; direct CAMB theta_* root",
        ),
        Model(
            "sterile-only control",
            lock.H0_SELECTOR,
            omega_b,
            omega_nur,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "nu_R only; zero-mode omitted",
        ),
    ], h0_root


def main() -> None:
    print("ITEM 123: ZERO-MODE STRUCTURE IMPLEMENTATION")
    print("=" * 96)
    print(f"CAMB version: {getattr(camb, '__version__', 'unknown')}")
    print(
        "Implementation rule: Brown--Kuchar zero-mode dust has w=c_s^2=0, "
        "so in linear CAMB it is entered as CDM and added to the sterile nu_R cold share."
    )

    zero_to_baryon, scatter_room, max_fraction, depletion = required_depletion_fraction()
    print("\n[1] Halo branch discipline")
    print(f"  fair-sample zero-mode / baryon mass ratio = {zero_to_baryon:.3f}")
    print(f"  RAR scatter proxy room                    = {scatter_room:.3f}")
    print(f"  MOND branch max zero-mode fraction        = {max_fraction:.3f}")
    print(f"  required depletion/screening              = {depletion:.1%}")
    check(depletion > 0.95, "active-MOND branch still needs a >95% depletion theorem")
    print("  Adopted branch for this run: N_zero is the CDM-like mobile halo mass;")
    print("  active R4/MOND is not also fitted as an independent galaxy force.")

    model_list, h0_root = models()
    results = {model.name: run_model(model) for model in model_list}
    reference = results["Planck LCDM reference"]
    comparisons = {
        name: compare_to_reference(result, reference)
        for name, result in results.items()
        if name != "Planck LCDM reference"
    }

    print("\n[2] Linear structure diagnostics at z=0")
    print(
        "model                         H0       100theta*  Omega_m  sigma8   S8      fSigma8  "
        "PkRMS  PkMax  shapeRMS"
    )
    for model in model_list:
        result = results[model.name]
        comp = comparisons.get(model.name)
        if comp is None:
            rms = maxp = shape = 0.0
        else:
            rms, maxp, shape = comp.rms_percent, comp.max_percent, comp.shape_rms_percent
        print(
            f"{model.name:28s} {model.H0:8.3f} {result.theta100:9.6f} "
            f"{result.omega_m:8.4f} {result.sigma8:7.4f} {result.s8:7.4f} "
            f"{result.fsigma8:8.4f} {rms:6.2f}% {maxp:6.2f}% {shape:7.2f}%"
        )

    selector = comparisons["framework selector-H0"]
    acoustic = comparisons["framework 63/64 acoustic"]
    root = comparisons["framework theta-root"]
    sterile = comparisons["sterile-only control"]

    print("\n[3] Gate checks")
    check(selector.rms_percent < 1.0, "selector-H0 zero-mode run matches Planck linear P(k) at sub-percent RMS")
    check(selector.shape_rms_percent < 1.0, "selector-H0 zero-mode run has sub-percent linear-shape residual")
    check(abs(selector.s8_delta_percent) < 1.0, "selector-H0 S8 is within 1% of Planck comparator")
    check(abs(selector.theta_delta_percent) > 0.2, "selector-H0 keeps the known theta_* offset")
    check(acoustic.rms_percent < 2.5, "63/64 acoustic run keeps linear P(k) within a few percent")
    check(acoustic.shape_rms_percent < 2.0, "63/64 acoustic run preserves the linear shape to ~2%")
    check(abs(acoustic.theta_delta_percent) < 0.005, "63/64 acoustic run closes theta_* against Planck comparator")
    check(abs(root.theta_delta_percent) < 0.005, "direct theta-root also closes theta_*")
    check(abs(root.rms_percent - acoustic.rms_percent) < 0.05, "theta-root and 63/64 are numerically the same structure branch")
    check(sterile.rms_percent > 50.0, "sterile-only control catastrophically fails the structure amplitude/shape gate")
    check(sterile.sigma8_delta_percent < -60.0, "sterile-only control underproduces growth")

    print("\n[4] Acoustic-scale interpretation")
    print(f"  CAMB theta-root H0                  = {h0_root:.6f}")
    print(f"  finite-depth candidate (63/64) H0   = {lock.H0_SELECTOR * lock.DEPTH6_FACTOR:.6f}")
    print(f"  root / candidate                    = {h0_root / (lock.H0_SELECTOR * lock.DEPTH6_FACTOR):.8f}")
    print("  The 63/64 option fixes theta_* and barely changes linear structure,")
    print("  but it remains a sector-coupling theorem unless the acoustic ruler is")
    print("  proven to bill the substrate busy-flag exposure rather than a different")
    print("  perturbation variable.")

    print("\n[5] Verdict")
    print(
        """
  The Boltzmann/structure implementation is now explicit:
    Brown--Kuchar N_zero is entered as CDM.  With nu_R+N_zero in omch2,
    the framework is LCDM-like for linear matter growth.  The selector-H0
    run matches the Planck comparator's z=0 linear P(k) at sub-percent RMS;
    the 63/64 acoustic run fixes theta_* and keeps the linear P(k) within a
    few percent.  The sterile-only control fails badly.

  The halo fork is no longer a free modelling choice in this implementation:
    N_zero is the CDM-like mobile halo mass.  Active R4/MOND is excluded as an
    independently fitted galaxy response unless a future >95% depletion or
    screening theorem is supplied.

  Remaining work:
    (i) a full likelihood/MCMC for theta_* and background parameters;
    (ii) nonlinear halo/galaxy simulations for the zero-mode branch;
    (iii) a sector-coupling theorem if the 63/64 acoustic pre-latch readout is
    to be promoted from numerically excellent to canon-Locked.
exit 0 -- zero-mode Brown--Kuchar dust passes the linear Boltzmann/structure implementation gate.
"""
    )


if __name__ == "__main__":
    main()
