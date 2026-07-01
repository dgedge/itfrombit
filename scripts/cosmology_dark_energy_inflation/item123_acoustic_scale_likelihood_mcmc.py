#!/usr/bin/env python3
r"""ITEM 123 / M15: compressed acoustic-scale likelihood and MCMC.

Purpose
-------
The zero-mode / nu_R cold-matter ledger can restore the CMB equality scale, but
the acoustic scale has a separate gate:

    selector H0 = 67.266 km/s/Mpc

predicts 100 theta_* high unless the CMB acoustic ruler sees the depth-6
pre-latch exposure

    H_CMB = (63/64) H_selector .

Previous scripts solved this deterministically.  This script puts a genuine
likelihood around the same question:

  * build a CAMB grid for 100 theta_*(H0) with the framework matter budget;
  * use a compressed Planck acoustic-scale likelihood

        chi2 = ((100 theta_* - target) / sigma_theta)^2;

  * sample the one-parameter posterior with a random-walk Metropolis chain;
  * report the selector-H0 and 63/64 branches in sigma and Delta-chi2 units.

Scope
-----
This is deliberately a *compressed acoustic-scale* likelihood, not the full
Planck TT/TE/EE/lensing likelihood and not a nuisance-parameter analysis.  Its
job is narrower: quantify whether the selector-H0 acoustic residual is an
empirical pressure or a harmless bookkeeping offset.  It is enough to classify
the acoustic-scale gate, but it does not replace a full CLASS/CAMB likelihood
run.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import random
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock


# Planck 2018 reports 100 theta_MC to about 3e-4.  We use the same absolute
# scale as a compressed acoustic likelihood while taking the target value from
# the local CAMB Planck-reference call used elsewhere in canon.
THETA100_SIGMA = 3.0e-4

GRID_MIN_H0 = 65.60
GRID_MAX_H0 = 67.60
GRID_POINTS = 121

MCMC_STEPS = 60_000
BURN_IN = 10_000
PROPOSAL_SIGMA_H0 = 0.055
RNG_SEED = 20260621


@dataclass(frozen=True)
class Branch:
    name: str
    h0: float
    theta100: float
    z: float
    chi2: float
    note: str


@dataclass(frozen=True)
class PosteriorSummary:
    mean: float
    median: float
    q16: float
    q84: float
    q025: float
    q975: float
    accept: float


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def theta_grid(omega_b: float, omega_dark: float) -> tuple[np.ndarray, np.ndarray]:
    h0 = np.linspace(GRID_MIN_H0, GRID_MAX_H0, GRID_POINTS)
    theta = np.array([lock.theta100_framework(float(x), omega_b, omega_dark) for x in h0])
    return h0, theta


def interpolated_theta(h0: float, h0_grid: np.ndarray, theta_grid_values: np.ndarray) -> float:
    if h0 < h0_grid[0] or h0 > h0_grid[-1]:
        return math.nan
    return float(np.interp(h0, h0_grid, theta_grid_values))


def log_like(h0: float, h0_grid: np.ndarray, theta_values: np.ndarray, target: float) -> float:
    theta = interpolated_theta(h0, h0_grid, theta_values)
    if not math.isfinite(theta):
        return -math.inf
    z = (theta - target) / THETA100_SIGMA
    return -0.5 * z * z


def run_mcmc(
    h0_grid: np.ndarray,
    theta_values: np.ndarray,
    target: float,
    start: float,
) -> tuple[np.ndarray, float]:
    rng = random.Random(RNG_SEED)
    current = start
    current_ll = log_like(current, h0_grid, theta_values, target)
    samples: list[float] = []
    accepted = 0

    for _ in range(MCMC_STEPS):
        proposal = current + rng.gauss(0.0, PROPOSAL_SIGMA_H0)
        proposal_ll = log_like(proposal, h0_grid, theta_values, target)
        if math.isfinite(proposal_ll) and math.log(rng.random()) < proposal_ll - current_ll:
            current = proposal
            current_ll = proposal_ll
            accepted += 1
        samples.append(current)

    return np.array(samples[BURN_IN:], dtype=float), accepted / MCMC_STEPS


def summarize(samples: np.ndarray, accept: float) -> PosteriorSummary:
    q025, q16, median, q84, q975 = np.quantile(samples, [0.025, 0.16, 0.5, 0.84, 0.975])
    return PosteriorSummary(
        mean=float(np.mean(samples)),
        median=float(median),
        q16=float(q16),
        q84=float(q84),
        q025=float(q025),
        q975=float(q975),
        accept=accept,
    )


def branch(name: str, h0: float, theta100: float, target: float, note: str) -> Branch:
    z = (theta100 - target) / THETA100_SIGMA
    return Branch(name=name, h0=h0, theta100=theta100, z=z, chi2=z * z, note=note)


def gaussian_tail_log10(z: float) -> float:
    """Approximate one-sided Gaussian tail log10 for large positive |z|."""

    x = abs(z)
    if x < 8.0:
        tail = 0.5 * math.erfc(x / math.sqrt(2.0))
        return math.log10(tail)
    # Mills-ratio asymptotic: Q(x) ~= exp(-x^2/2)/(x sqrt(2pi)).
    return (-0.5 * x * x - math.log(x * math.sqrt(2.0 * math.pi))) / math.log(10.0)


def direct_slope_near_root(h0_root: float, omega_b: float, omega_dark: float) -> float:
    dh = 0.05
    up = lock.theta100_framework(h0_root + dh, omega_b, omega_dark)
    down = lock.theta100_framework(h0_root - dh, omega_b, omega_dark)
    return (up - down) / (2.0 * dh)


def main() -> None:
    print("ITEM 123 / M15: ACOUSTIC-SCALE COMPRESSED LIKELIHOOD + MCMC")
    print("=" * 96)
    print(f"CAMB version: {getattr(lock.camb, '__version__', 'unknown')}")
    print(f"Compressed likelihood sigma(100 theta_*) = {THETA100_SIGMA:.6f}")
    print(f"MCMC seed={RNG_SEED}, steps={MCMC_STEPS}, burn-in={BURN_IN}")

    omega_b, omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_6364 = lock.H0_SELECTOR * lock.DEPTH6_FACTOR
    theta_selector = lock.theta100_framework(lock.H0_SELECTOR, omega_b, omega_dark)
    theta_root = lock.theta100_framework(h0_root, omega_b, omega_dark)
    theta_6364 = lock.theta100_framework(h0_6364, omega_b, omega_dark)

    print("\n[1] Inputs and direct CAMB branch values")
    print(f"  omega_b h^2                 = {omega_b:.8f}")
    print(f"  omega_nuR h^2               = {omega_nur:.8f}")
    print(f"  omega_dark h^2              = {omega_dark:.8f}")
    print(f"  Planck/CAMB reference 100theta = {target:.9f}")
    branches = [
        branch(
            "selector H0",
            lock.H0_SELECTOR,
            theta_selector,
            target,
            "late selector H0, no acoustic pre-latch repair",
        ),
        branch(
            "CAMB theta root",
            h0_root,
            theta_root,
            target,
            "H0 freely shifted to match theta_*",
        ),
        branch(
            "63/64 pre-latch",
            h0_6364,
            theta_6364,
            target,
            "depth-6 busy/pre-latch acoustic readout candidate",
        ),
    ]
    print("  branch             H0          100theta       z       chi2        note")
    for b in branches:
        print(f"  {b.name:16s} {b.h0:10.6f}  {b.theta100:11.9f} {b.z:8.2f} {b.chi2:10.2f}  {b.note}")

    print("\n[2] CAMB-grid likelihood table")
    h0_values, theta_values = theta_grid(omega_b, omega_dark)
    grid_best_idx = int(np.argmax([log_like(float(x), h0_values, theta_values, target) for x in h0_values]))
    print(f"  grid range                    = [{h0_values[0]:.3f}, {h0_values[-1]:.3f}]")
    print(f"  grid points                   = {GRID_POINTS}")
    print(f"  best grid H0                  = {h0_values[grid_best_idx]:.6f}")
    print(f"  direct CAMB bisection root     = {h0_root:.6f}")
    check(abs(h0_values[grid_best_idx] - h0_root) < 0.02, "grid maximum resolves the CAMB theta root")
    check(np.all(np.diff(theta_values) > 0.0), "theta100(H0) is monotone over the likelihood grid")

    print("\n[3] One-parameter posterior")
    samples, accept = run_mcmc(h0_values, theta_values, target, h0_root)
    summary = summarize(samples, accept)
    factor_samples = samples / lock.H0_SELECTOR
    factor_summary = summarize(factor_samples, accept)
    slope = direct_slope_near_root(h0_root, omega_b, omega_dark)
    h0_sigma_linear = abs(THETA100_SIGMA / slope)
    print(f"  acceptance                    = {summary.accept:.3f}")
    print(f"  d(100theta)/dH0 at root        = {slope:.8f} per km/s/Mpc")
    print(f"  linear sigma(H0)              = {h0_sigma_linear:.5f}")
    print(
        f"  H0 posterior median           = {summary.median:.6f} "
        f"(-{summary.median-summary.q16:.6f}, +{summary.q84-summary.median:.6f})"
    )
    print(f"  H0 95% interval               = [{summary.q025:.6f}, {summary.q975:.6f}]")
    print(
        f"  readout factor median         = {factor_summary.median:.8f} "
        f"(-{factor_summary.median-factor_summary.q16:.8f}, +{factor_summary.q84-factor_summary.median:.8f})"
    )
    print(f"  readout factor 95% interval   = [{factor_summary.q025:.8f}, {factor_summary.q975:.8f}]")
    print(f"  63/64                         = {lock.DEPTH6_FACTOR:.8f}")
    check(summary.q025 < h0_root < summary.q975, "direct CAMB theta root lies inside the posterior")
    check(abs(factor_summary.median - lock.DEPTH6_FACTOR) < 5.0e-4, "posterior readout factor is centered on the 63/64 candidate")

    selector = branches[0]
    pre = branches[2]
    print("\n[4] Acoustic-scale gate verdict")
    print(f"  selector-H0 one-sided tail log10 p  = {gaussian_tail_log10(selector.z):.2f}")
    print(f"  selector-H0 Delta chi2              = {selector.chi2:.2f}")
    print(f"  63/64 Delta chi2                    = {pre.chi2:.5f}")
    print(f"  selector H0 / posterior sigma       = {(lock.H0_SELECTOR - summary.median) / h0_sigma_linear:.2f}")
    print(f"  (63/64 H0 - posterior median)/sigma = {(h0_6364 - summary.median) / h0_sigma_linear:.2f}")
    check(selector.chi2 > 100.0, "selector-H0 without acoustic repair is decisively rejected by the compressed theta likelihood")
    check(pre.chi2 < 1.0e-2, "63/64 pre-latch candidate is effectively at the acoustic maximum")

    print("\n[5] Interpretation")
    print(
        """
  The real-likelihood statement is now quantitative:
    With the framework matter budget held fixed, a compressed Planck acoustic
    likelihood places the CMB H0/readout factor at the CAMB theta root,
    H0 ~= 66.218 km/s/Mpc, i.e. almost exactly (63/64) H_selector.

    The unmodified selector-H0 branch is not a mild offset: in this compressed
    likelihood it costs Delta chi2 > 100.  Therefore the acoustic-scale issue
    is a genuine empirical gate.

    The 63/64 branch is statistically excellent in this one-number likelihood,
    but still not a full canon lock unless the CMB acoustic ruler is proved to
    bill the depth-6 busy/pre-latch service exposure rather than the completed
    late-time selector span.

  What this is not:
    It is not a full Planck likelihood, and it does not include TT/TE/EE
    nuisance parameters, lensing, BAO, supernovae, or nonlinear structure.
    It is the durable compressed likelihood for the specific theta_* gate.
exit 0 -- acoustic-scale compressed MCMC: selector-H0 rejected; 63/64 branch passes the theta gate conditionally.
"""
    )


if __name__ == "__main__":
    main()
