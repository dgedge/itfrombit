#!/usr/bin/env python3
r"""ITEM 123: halo-branch depletion audit.

Question
--------
The CMB branch wants the R4 zero-mode reservoir to be pressureless CDM-like
dust.  The galaxy branch already has an active R4/MOND line-current mechanism.
Those two uses double-count the galaxy discrepancy unless one of two things is
true:

  A. zero-mode is the CDM-like halo component and active R4/MOND is not fitted
     as an independent extra force in the same galaxy regime; or

  B. active R4/MOND is kept as the baryonic RAR law, and the zero-mode is
     depleted/screened from galaxies strongly enough to hide the fair-sampling
     CDM-like halo.

The global accounting gate says branch B needs about 96% zero-mode depletion.
This script adds a modest halo-model Monte Carlo: put the zero-mode into NFW
CDM-like halos, use simple disk radii, and ask how much local screening would
be needed at optical and outer-disk/RAR radii.

Verdict
-------
The halo model does not produce depletion.  Collisionless zero-mode-as-CDM is
present in galaxies and becomes more important at larger radii and in
low-baryon systems.  A MOND-separate branch would require order-90--97% local
screening, with >95% required for a majority of realistic abundance-matching
halos by 5--8 disk scale lengths.  No such screening operator exists in the
current zero-mode inventory.

Therefore the current non-double-counted branch is zero-mode-as-CDM halos plus
at most subleading active-R4 polarisation.  Keeping MOND as an independent
galaxy force requires a new depletion/screening theorem.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import random

import numpy as np

from item123_cmb_cold_matter_closure_status import (
    OMEGA_B_H2,
    OMEGA_C_OBS_H2,
    RAR_SCATTER_DEX,
    density_for_power,
)
from item123_r4_zero_mode_dust_hamiltonian import operator_inventory


SEED = 123132
N_HALOS = 50_000


@dataclass(frozen=True)
class HaloModel:
    name: str
    description: str


MODELS = (
    HaloModel(
        "cosmic-max-disk",
        "all cosmic baryons retained in the disk; deliberately minimises the zero-mode problem",
    ),
    HaloModel(
        "high-disk",
        "large retained disk, M_b=min(f_b M_200, 0.06 M_200)",
    ),
    HaloModel(
        "abundance-matching",
        "Moster-like stellar mass plus gas cap; closer to observed low baryon retention",
    ),
)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def f_nfw(x: float) -> float:
    return math.log1p(x) - x / (1.0 + x)


def disk_fraction(n_scale_lengths: float) -> float:
    """Mass fraction of an exponential disk inside n scale lengths."""

    return 1.0 - math.exp(-n_scale_lengths) * (1.0 + n_scale_lengths)


def concentration_median(m200: float) -> float:
    """Adequate z=0 concentration proxy for this fork audit."""

    return 10.0 * (m200 / 1.0e12) ** (-0.10)


def moster_stellar_mass(m200: float) -> float:
    """Simple Moster-style stellar-to-halo map; used only as a bracket."""

    n = 0.0351
    m1 = 10.0**11.59
    beta = 1.376
    gamma = 0.608
    return m200 * 2.0 * n / ((m200 / m1) ** (-beta) + (m200 / m1) ** gamma)


def gas_multiplier(mstar: float) -> float:
    """Crude gas-rich correction, capped so the bracket stays conservative."""

    return 1.0 + min(8.0, (1.0e10 / max(mstar, 1.0)) ** 0.45)


def baryon_mass(model: str, m200: float, f_baryon_total: float) -> float:
    if model == "cosmic-max-disk":
        return f_baryon_total * m200
    if model == "high-disk":
        return min(f_baryon_total * m200, 0.06 * m200)
    if model == "abundance-matching":
        mstar = moster_stellar_mass(m200)
        return min(f_baryon_total * m200, mstar * gas_multiplier(mstar))
    raise ValueError(model)


def halo_sample_rows(
    *,
    model: str,
    radius_scale_lengths: float,
    n_halos: int = N_HALOS,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Return zero/baryon acceleration ratios and required depletion fractions."""

    model_offsets = {
        "cosmic-max-disk": 101,
        "high-disk": 211,
        "abundance-matching": 307,
    }
    radius_offset = int(round(100.0 * radius_scale_lengths))
    rng = random.Random(seed + model_offsets[model] + 10_000 * radius_offset)
    _omega_nur, omega_zero, omega_dark, _zeq = density_for_power(1)
    omega_m_h2 = OMEGA_B_H2 + omega_dark
    f_baryon_total = OMEGA_B_H2 / omega_m_h2
    f_zero_total = omega_zero / omega_m_h2
    scatter_room = 10.0**RAR_SCATTER_DEX - 1.0

    ratios: list[float] = []
    depletions: list[float] = []
    for _ in range(n_halos):
        # Broad galaxy halo range.  Uniform-in-log is intentional: this is a
        # branch audit, not an observed halo catalogue likelihood.
        log_m = rng.uniform(10.0, 13.0)
        m200 = 10.0**log_m
        concentration = concentration_median(m200) * 10.0 ** rng.gauss(0.0, 0.11)
        spin = 0.035 * 10.0 ** rng.gauss(0.0, 0.25)
        r200_kpc = 206.0 * (m200 / 1.0e12) ** (1.0 / 3.0)
        rd_kpc = spin * r200_kpc / math.sqrt(2.0)
        radius_kpc = radius_scale_lengths * rd_kpc

        enclosed_zero = (
            f_zero_total
            * m200
            * f_nfw(concentration * radius_kpc / r200_kpc)
            / f_nfw(concentration)
        )
        enclosed_baryon = (
            baryon_mass(model, m200, f_baryon_total)
            * disk_fraction(radius_scale_lengths)
        )

        # Acceleration ratio at the same radius is just enclosed-mass ratio in
        # this sphericalised gate.  Disk geometry changes O(1) numbers, not the
        # branch verdict.
        ratio = enclosed_zero / max(enclosed_baryon, 1.0e-300)
        allowed_fraction = scatter_room / ratio if ratio > 0.0 else 1.0
        depletion = max(0.0, 1.0 - allowed_fraction)
        ratios.append(ratio)
        depletions.append(depletion)

    return np.array(ratios), np.array(depletions)


def summarize(arr: np.ndarray) -> tuple[float, float, float]:
    p16, p50, p84 = np.percentile(arr, [16.0, 50.0, 84.0])
    return float(p16), float(p50), float(p84)


def operator_has_screening() -> bool:
    for term in operator_inventory():
        if term.status.startswith("ADMITTED") and any(
            token in term.form.lower() or token in term.name.lower()
            for token in ("gradient", "hopping", "stiffness", "screen", "shear")
        ):
            return True
    return False


def main() -> None:
    print("ITEM 123: HALO-BRANCH ZERO-MODE DEPLETION HALO-MODEL AUDIT")
    print("=" * 104)

    omega_nur, omega_zero, omega_dark, zeq = density_for_power(1)
    zero_to_baryon = omega_zero / OMEGA_B_H2
    scatter_room = 10.0**RAR_SCATTER_DEX - 1.0
    max_fair_fraction = scatter_room / zero_to_baryon
    global_depletion = 1.0 - max_fair_fraction

    print("\n[1] Global non-double-counting gate")
    print(f"  omega_b h^2                         = {OMEGA_B_H2:.6f}")
    print(f"  omega_nuR h^2                       = {omega_nur:.6f}")
    print(f"  omega_zero h^2                      = {omega_zero:.6f}")
    print(f"  omega_dark h^2                      = {omega_dark:.6f}")
    print(f"  omega_c target h^2                  = {OMEGA_C_OBS_H2:.6f}")
    print(f"  z_eq                                = {zeq:.1f}")
    print(f"  zero-mode / baryon cosmic ratio     = {zero_to_baryon:.3f}")
    print(f"  RAR scatter room ({RAR_SCATTER_DEX:.2f} dex)           = {scatter_room:.3f}")
    print(f"  max fair-sample zero-mode in MOND branch = {max_fair_fraction:.3f}")
    print(f"  global depletion required           = {global_depletion:.1%}")
    check(abs(omega_dark / OMEGA_C_OBS_H2 - 1.0) < 0.01, "zero-mode + nu_R is a CDM-sized reservoir")
    check(global_depletion > 0.95, "independent-MOND branch requires >95% global zero-mode depletion")

    print("\n[2] Local halo-model depletion required if MOND is kept independent")
    print("  Quantity shown: zero-mode acceleration / baryon acceleration, then")
    print("  depletion needed so the residual zero-mode acceleration fits inside")
    print("  the RAR scatter room.  A true CDM-like zero-mode has no such depletion.")
    print(
        f"  {'model':<20} {'R/Rd':>5} {'g_z/g_b p16,p50,p84':>28} "
        f"{'depl p16,p50,p84':>27} {'P(depl>95%)':>13}"
    )
    summaries: dict[tuple[str, float], tuple[float, float, float, float, float, float, float]] = {}
    for radius in (2.2, 5.0, 8.0):
        for model in MODELS:
            ratios, depletions = halo_sample_rows(model=model.name, radius_scale_lengths=radius)
            r16, r50, r84 = summarize(ratios)
            d16, d50, d84 = summarize(depletions)
            p95 = float(np.mean(depletions > 0.95))
            summaries[(model.name, radius)] = (r16, r50, r84, d16, d50, d84, p95)
            print(
                f"  {model.name:<20} {radius:5.1f} "
                f"{r16:7.3f},{r50:7.3f},{r84:7.3f} "
                f"{d16:7.3f},{d50:7.3f},{d84:7.3f} "
                f"{p95:13.3f}"
            )

    moster_5 = summaries[("abundance-matching", 5.0)]
    moster_8 = summaries[("abundance-matching", 8.0)]
    cosmic_22 = summaries[("cosmic-max-disk", 2.2)]
    check(cosmic_22[4] < 0.75, "maximal baryon disk is the conservative lower-depletion bound")
    check(moster_5[4] > 0.94, "realistic outer-disk halo model needs about 95% median depletion")
    check(moster_8[4] > 0.96, "outer RAR radii need >95% median depletion")
    check(moster_8[6] > 0.65, "most realistic outer-disk samples require >95% depletion")

    print("\n[3] Does current zero-mode dynamics supply that depletion?")
    for term in operator_inventory():
        print(f"  {term.name:28s} {term.status:23s} {term.form}")
        print(f"    {term.reason}")
    has_screen = operator_has_screening()
    check(not has_screen, "current admitted zero-mode inventory has no spatial screening/depletion operator")
    print("  Collisionless CDM-like dust survives in the halo model with survival")
    print("  fraction 1.  Baryonic concentration would contract it rather than remove")
    print("  it.  Therefore the required 90--97% local suppression is not generated by")
    print("  the zero-mode-as-CDM dynamics.")

    print("\n[4] Fork verdict")
    print("  zero-mode-as-CDM branch:")
    print("    viable bookkeeping branch for CMB/mobile halo mass; active R4/MOND must")
    print("    be treated as subleading polarisation or not independently fitted.")
    print("  MOND-separate branch:")
    print("    still possible only after adding a new galaxy depletion/screening")
    print("    theorem.  The halo model shows the required screening is severe, and")
    print("    the current finite operator inventory does not supply it.")
    print("  forbidden branch:")
    print("    fair-sample zero-mode CDM halo plus independent MOND/R4 force; this")
    print("    double-counts the galaxy discrepancy.")
    print("exit 0 -- no >95% depletion theorem; adopt zero-mode-CDM bookkeeping or add new screening physics.")


if __name__ == "__main__":
    main()
