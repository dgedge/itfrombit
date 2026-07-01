#!/usr/bin/env python3
r"""ITEM 123 / M15: external CMB spectra + nonlinear structure gate.

Context
-------
The local acoustic theorem now identifies the CMB acoustic readout as

    P_acoustic = P_busy = I - |111111><111111|,
    H_CMB/H_selector = 63/64.

That removes the local projector as the blocker.  The remaining work is
phenomenological/external:

  1. full Planck TT/TE/EE/lensing likelihood, including nuisance parameters;
  2. nonlinear matter/halo modelling for the zero-mode-CDM branch.

This script is a durable *pre-likelihood* external gate.  It does not claim to
be Planck plik/CamSpec/cobaya.  Instead it runs CAMB for the exact framework
matter budget and compares against a Planck-2018 LCDM CAMB comparator in:

  * lensed TT, EE, TE and lensing-potential spectra;
  * linear and Halofit nonlinear matter power at z=0, 0.5, 1, 2;
  * sigma8/S8 and f sigma8 diagnostics;
  * the sterile-only control.

Verdict language:
  * if these spectra fail badly, the framework fails before a real likelihood;
  * if they pass at percent/few-percent level, the next object is genuinely a
    full likelihood and nonlinear halo catalogue, not another local theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import shutil
import sys

try:
    import camb
    from camb import model
    import numpy as np
except Exception as exc:  # pragma: no cover - dependency gate
    raise SystemExit(
        "CAMB/NumPy dependency missing. Use the project Python environment, e.g.\n"
        "  ~/bin/py13_7/bin/python python_code/item123_external_cmb_lensing_nonlinear_gate.py\n"
        f"Import failure: {type(exc).__name__}: {exc}"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock


LMAX = 2500
K_MIN = 1.0e-3
K_MAX = 5.0
K_POINTS = 220
REDSHIFTS = [0.0, 0.5, 1.0, 2.0]


@dataclass(frozen=True)
class RunSpec:
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
class RunData:
    spec: RunSpec
    theta100: float
    omega_m: float
    sigma8: float
    s8: float
    fsigma8_z0: float
    ell: np.ndarray
    total: np.ndarray
    lens_potential: np.ndarray
    kh_lin: np.ndarray
    pk_lin: np.ndarray
    kh_nl: np.ndarray
    pk_nl: np.ndarray


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def camb_params(spec: RunSpec, *, nonlinear: bool) -> camb.CAMBparams:
    pars = camb.set_params(
        H0=spec.H0,
        ombh2=spec.ombh2,
        omch2=spec.omch2,
        ns=spec.ns,
        As=spec.As,
        tau=spec.tau,
        w=spec.w0,
        wa=spec.wa,
        dark_energy_model="ppf",
    )
    pars.set_for_lmax(LMAX, lens_potential_accuracy=2)
    pars.set_matter_power(redshifts=REDSHIFTS, kmax=K_MAX)
    pars.NonLinear = model.NonLinear_both if nonlinear else model.NonLinear_none
    return pars


def get_pk(spec: RunSpec, *, nonlinear: bool) -> tuple[np.ndarray, np.ndarray]:
    pars = camb_params(spec, nonlinear=nonlinear)
    res = camb.get_results(pars)
    kh, redshifts, pk = res.get_matter_power_spectrum(
        minkh=K_MIN,
        maxkh=K_MAX,
        npoints=K_POINTS,
    )
    # CAMB returns redshifts sorted earliest first, so reorder back to REDSHIFTS.
    order = [list(np.round(redshifts, 8)).index(round(z, 8)) for z in REDSHIFTS]
    return kh, pk[order]


def run(spec: RunSpec) -> RunData:
    pars = camb_params(spec, nonlinear=True)
    res = camb.get_results(pars)
    powers = res.get_cmb_power_spectra(pars, CMB_unit="muK")
    kh_lin, pk_lin = get_pk(spec, nonlinear=False)
    kh_nl, pk_nl = get_pk(spec, nonlinear=True)
    omega_m = (spec.ombh2 + spec.omch2) / (spec.H0 / 100.0) ** 2
    # CAMB sorts transfer redshifts earliest first, so with REDSHIFTS
    # [0, 0.5, 1, 2] the z=0 growth entries are the final elements.
    sigma8_z0 = float(res.get_sigma8()[-1])
    fsigma8_z0 = float(res.get_fsigma8()[-1])
    return RunData(
        spec=spec,
        theta100=float(res.cosmomc_theta() * 100.0),
        omega_m=omega_m,
        sigma8=sigma8_z0,
        s8=sigma8_z0 * math.sqrt(omega_m / 0.3),
        fsigma8_z0=fsigma8_z0,
        ell=np.arange(powers["total"].shape[0]),
        total=powers["total"],
        lens_potential=powers["lens_potential"],
        kh_lin=kh_lin,
        pk_lin=pk_lin,
        kh_nl=kh_nl,
        pk_nl=pk_nl,
    )


def rms_frac(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    frac = a[mask] / b[mask] - 1.0
    return float(np.sqrt(np.mean(frac**2)))


def rms_norm(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a[mask] - b[mask]) ** 2)) / np.sqrt(np.mean(b[mask] ** 2)))


def spectrum_metrics(data: RunData, ref: RunData) -> dict[str, float]:
    length = min(len(data.ell), len(ref.ell))
    ell = np.arange(length)
    acoustic = (ell >= 30) & (ell <= 2000)
    high = (ell >= 1000) & (ell <= 2000)
    lens_len = min(len(data.lens_potential), len(ref.lens_potential))
    lens_ell = np.arange(lens_len)
    lens_mask = (lens_ell >= 10) & (lens_ell <= 1000) & (ref.lens_potential[:lens_len, 0] > 0)

    fw = data.total[:length]
    pl = ref.total[:length]
    out = {
        "TT_rms": rms_frac(fw[:, 0], pl[:, 0], acoustic),
        "EE_rms": rms_frac(fw[:, 1], pl[:, 1], acoustic & (pl[:, 1] > 1.0e-6)),
        "TE_norm": rms_norm(fw[:, 3], pl[:, 3], acoustic),
        "TT_high": rms_frac(fw[:, 0], pl[:, 0], high),
        "EE_high": rms_frac(fw[:, 1], pl[:, 1], high & (pl[:, 1] > 1.0e-6)),
        "lens_rms": rms_frac(
            data.lens_potential[:lens_len, 0],
            ref.lens_potential[:lens_len, 0],
            lens_mask,
        ),
        "theta_delta": data.theta100 / ref.theta100 - 1.0,
        "sigma8_delta": data.sigma8 / ref.sigma8 - 1.0,
        "s8_delta": data.s8 / ref.s8 - 1.0,
    }
    return out


def pk_metrics(data: RunData, ref: RunData, *, nonlinear: bool) -> dict[str, list[float]]:
    kh = data.kh_nl if nonlinear else data.kh_lin
    pk = data.pk_nl if nonlinear else data.pk_lin
    rpk = ref.pk_nl if nonlinear else ref.pk_lin
    bands = {
        "lin_0.01_0.3": (kh >= 0.01) & (kh <= 0.3),
        "quasi_0.3_1": (kh > 0.3) & (kh <= 1.0),
        "nonlin_1_5": (kh > 1.0) & (kh <= 5.0),
    }
    out: dict[str, list[float]] = {name: [] for name in bands}
    for iz in range(len(REDSHIFTS)):
        for name, mask in bands.items():
            out[name].append(rms_frac(pk[iz], rpk[iz], mask))
    return out


def models() -> list[RunSpec]:
    omega_b, omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_6364 = lock.H0_SELECTOR * lock.DEPTH6_FACTOR
    return [
        RunSpec(
            "Planck LCDM comparator",
            lock.H0_PLANCK,
            lock.OMBH2_PLANCK,
            lock.OMCH2_PLANCK,
            lock.NS_PLANCK,
            lock.AS_PLANCK,
            lock.TAU_PLANCK,
            -1.0,
            0.0,
            "external comparator, not a fit to framework data",
        ),
        RunSpec(
            "framework selector-H0",
            lock.H0_SELECTOR,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "late selector H0; no 63/64 acoustic readout",
        ),
        RunSpec(
            "framework 63/64",
            h0_6364,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "local busy-projector acoustic readout",
        ),
        RunSpec(
            "framework theta-root",
            h0_root,
            omega_b,
            omega_dark,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "direct CAMB theta root; should coincide with 63/64",
        ),
        RunSpec(
            "sterile-only control",
            lock.H0_SELECTOR,
            omega_b,
            omega_nur,
            lock.NS_FW,
            lock.AS_FW,
            lock.TAU_FW,
            lock.W0_FW,
            lock.WA_FW,
            "zero-mode omitted; known equality/structure failure control",
        ),
    ]


def planck_likelihood_probe() -> tuple[bool, str]:
    probes = ("clik", "cobaya", "planckpr4lensing")
    found = [name for name in probes if shutil.which(name)]
    if found:
        return True, ", ".join(found)
    try:
        import clik  # type: ignore  # noqa: F401
        return True, "python module clik"
    except Exception:
        return False, "no clik/cobaya/Planck likelihood command or module found"


def print_percent(value: float) -> str:
    return f"{100.0 * value:7.2f}%"


def main() -> None:
    print("ITEM 123 / M15: EXTERNAL CMB + LENSING + NONLINEAR STRUCTURE GATE")
    print("=" * 104)
    print(f"CAMB version: {getattr(camb, '__version__', 'unknown')}")
    available, detail = planck_likelihood_probe()
    print(f"Planck likelihood probe: {'available' if available else 'not available'} ({detail})")
    print("This script is a pre-likelihood CAMB comparator, not plik/CamSpec/cobaya.")

    specs = models()
    results = {spec.name: run(spec) for spec in specs}
    ref = results["Planck LCDM comparator"]

    print("\n[1] Background and growth diagnostics")
    print("model                     H0       100theta    Omega_m  sigma8   S8      fSigma8")
    for spec in specs:
        data = results[spec.name]
        print(
            f"{spec.name:25s} {spec.H0:8.3f} {data.theta100:10.6f} "
            f"{data.omega_m:8.4f} {data.sigma8:7.4f} {data.s8:7.4f} {data.fsigma8_z0:8.4f}"
        )

    print("\n[2] Lensed CMB spectra against the Planck LCDM CAMB comparator")
    print("model                     TT rms   EE rms   TE norm  lens rms  theta")
    spec_metrics: dict[str, dict[str, float]] = {}
    for name, data in results.items():
        if name == ref.spec.name:
            continue
        m = spectrum_metrics(data, ref)
        spec_metrics[name] = m
        print(
            f"{name:25s} {print_percent(m['TT_rms'])} {print_percent(m['EE_rms'])} "
            f"{print_percent(m['TE_norm'])} {print_percent(m['lens_rms'])} {print_percent(m['theta_delta'])}"
        )

    print("\n[3] Linear and Halofit nonlinear matter-power RMS differences")
    print("model                     type      z=0     z=0.5   z=1     z=2       band")
    pk_store: dict[tuple[str, str], dict[str, list[float]]] = {}
    for name, data in results.items():
        if name == ref.spec.name:
            continue
        for label, nonlinear in (("linear", False), ("halofit", True)):
            metrics = pk_metrics(data, ref, nonlinear=nonlinear)
            pk_store[(name, label)] = metrics
            for band, values in metrics.items():
                print(
                    f"{name:25s} {label:8s} "
                    + " ".join(print_percent(v) for v in values)
                    + f"  {band}"
                )

    acoustic = spec_metrics["framework 63/64"]
    selector = spec_metrics["framework selector-H0"]
    sterile = spec_metrics["sterile-only control"]
    acoustic_nl = pk_store[("framework 63/64", "halofit")]
    acoustic_lin = pk_store[("framework 63/64", "linear")]
    selector_lin = pk_store[("framework selector-H0", "linear")]

    print("\n[4] Gate checks")
    check(acoustic["TT_rms"] < 0.025, "63/64 branch keeps lensed TT within a few percent of LCDM comparator")
    check(acoustic["EE_rms"] < 0.04, "63/64 branch keeps lensed EE within a few percent")
    check(acoustic["TE_norm"] < 0.06, "63/64 branch keeps TE within a few percent in normalized RMS")
    check(acoustic["lens_rms"] < 0.08, "63/64 branch keeps lensing-potential spectrum within a few-to-10 percent")
    check(abs(acoustic["theta_delta"]) < 5.0e-5, "63/64 branch closes the acoustic scale")
    check(abs(selector["theta_delta"]) > 2.0e-3, "selector-H0 branch keeps the known theta-star failure")
    check(sterile["TT_rms"] > 0.10, "sterile-only control fails CMB spectra")
    check(max(acoustic_lin["lin_0.01_0.3"]) < 0.03, "63/64 linear P(k) remains close on structure scales")
    check(max(acoustic_nl["lin_0.01_0.3"]) < 0.04, "63/64 Halofit P(k) remains close on linear/quasi-linear scales")
    check(max(acoustic_nl["quasi_0.3_1"]) < 0.06, "63/64 Halofit P(k) remains controlled at k=0.3..1 h/Mpc")
    check(max(selector_lin["lin_0.01_0.3"]) < 0.02, "selector branch still has near-LCDM linear growth despite theta-star failure")

    print("\n[5] Interpretation")
    print(
        """
  What closes here:
    The zero-mode-CDM branch survives the next external CAMB gate.  With the
    local 63/64 acoustic readout, TT/EE/TE/lensing spectra and Halofit
    nonlinear matter power remain close to a Planck-2018 LCDM comparator.  The
    sterile-only control fails, so the zero-mode reservoir is still load-bearing.

  What does not close here:
    This is not a full Planck likelihood.  No plik/CamSpec/cobaya nuisance
    likelihood is available in this environment, and no nonlinear halo catalogue
    or galaxy-scale baryonic modelling is run.  A real publication-grade test
    still needs (i) a Planck TTTEEE+lensing likelihood with nuisance parameters
    and (ii) nonlinear halo modelling for the adopted zero-mode-CDM branch.

  Updated blocker:
    The blocker is no longer the local 63/64 projector.  It is the external
    likelihood/halo pipeline.  The first-pass spectra do not show an obvious
    kill-shot; they justify doing the expensive likelihood rather than arguing
    by prose.
exit 0 -- external CAMB comparator gate passed; full likelihood and halo catalogue remain.
"""
    )


if __name__ == "__main__":
    main()
