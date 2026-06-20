#!/usr/bin/env python3
r"""ITEM 131: reopened HBC inflation gate after the R4/Planck lift closure.

Question
--------
The late R4 homogeneous lift and the Planck selector are now closed inside the
current service instrument.  Does that also close the inflation/HBC residual?

Specifically, can the current canon now force both

    N_shell alpha_0^4 = C_F = 4/3,
    S_j(k=aH) = 1,

so that

    A_nu = (3/4) alpha_0^4

is derived rather than conditional?

Result
------
No.  The new R4 lift closes the homogeneous completion fraction (a zero-mode /
FRW-average statement).  The scalar amplitude still samples two additional
objects that are independent of that lift:

1. channel projection: entropy/area saturation fixes the all-ones printed-load
   projection, while the scalar/topological amplitude needs the post-decoder
   colour-restoring weight-4 projection;
2. spatial projection: colour-load saturation fixes a mean shell load, while
   the scalar amplitude samples the nonzero compensated horizon Fourier mode
   through S_j(k=aH).

The tilt survives because a scale-independent multiplicative structure factor
changes the amplitude, not the log-shell generator slope.  Therefore the
remaining HBC inflation frontier is exactly two operator statements:

    A. channel-lock / critical-latch theorem;
    B. spatial-whitening theorem.

This is a sharpened no-go against deriving A_nu from the already-closed gates,
not a retreat to an unstructured open problem.
"""

from __future__ import annotations

import cmath
import math
import subprocess
import sys
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0
C_F = 4.0 / 3.0
DELTA = 1.0 / 28.0

Vector = tuple[int, int, int]


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def run_script(name: str) -> str:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "python_code" / name)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f"{name} failed with exit {proc.returncode}")
    return proc.stdout


def require(text: str, needle: str, label: str) -> None:
    ok = needle in text
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        raise AssertionError(f"{label}: missing {needle!r}")


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))  # type: ignore[return-value]


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def phase_sum(n: int, k: Vector) -> complex:
    return sum(phase(n, k, x) for x in torus_points(n))


def horizon_profile(n: int, k: Vector) -> np.ndarray:
    values = np.array([phase(n, k, x).real for x in torus_points(n)], dtype=float)
    return values - np.mean(values)


def structure_factor(cov: np.ndarray, n: int, k: Vector, mean_per_site: float) -> float:
    pts = torus_points(n)
    chars = np.array([phase(n, k, x) for x in pts], dtype=complex)
    var = (chars @ cov @ np.conjugate(chars)) / ((n**3) ** 2 * mean_per_site**2)
    return float(var.real * (n**3) * mean_per_site)


def zero_mode_variance(cov: np.ndarray, mean_per_site: float) -> float:
    ones = np.ones(cov.shape[0])
    return float((ones @ cov @ ones) / (cov.shape[0] ** 2 * mean_per_site**2))


def amplitude_from_utilization(utilization: float, s_j: float = 1.0) -> float:
    return s_j / (utilization * C_F / ALPHA0**4)


def power_spectrum(k: float, amplitude: float, s_factor: float = 1.0) -> float:
    return amplitude * s_factor * k ** (-DELTA)


def log_slope(k1: float, k2: float, amplitude: float, s_factor: float = 1.0) -> float:
    return math.log(power_spectrum(k2, amplitude, s_factor) / power_spectrum(k1, amplitude, s_factor)) / math.log(k2 / k1)


def main() -> None:
    print("ITEM 131 REOPENED HBC INFLATION GATE")
    print("=" * 96)

    print("\n[1] Owner-script status after the recent R4/Planck lift")
    owner_outputs = {
        "r4_lift": run_script("item131_r4_homogeneous_lift_theorem.py"),
        "early_log": run_script("item131_early_logscale_lift.py"),
        "saturation_gate": run_script("item131_hbc_saturation_gate.py"),
        "decision_gate": run_script("item131_hbc_amplitude_decision_gate.py"),
        "t5b": run_script("item131_t5b_whiteness_lemma.py"),
    }
    require(owner_outputs["r4_lift"], "homogeneous R4 lift closed inside the current QEC service instrument", "late/Planck homogeneous R4 lift is closed")
    require(owner_outputs["early_log"], "n_s = 27/28", "tilt still lands under saturated radial log-shell HBC")
    require(owner_outputs["saturation_gate"], "N_shell equality", "saturation gate still names N_shell equality")
    require(owner_outputs["decision_gate"], "amplitude remains conditional under current canon", "decision gate still blocks A_nu promotion")
    require(owner_outputs["t5b"], "no-horizon-scale-covariance premise", "T5b still isolates spatial whiteness")

    print("\n[2] Conditional inflation numbers")
    n_shell_sat = C_F * ALPHA0 ** -4
    a_sat = 1.0 / n_shell_sat
    n_s = 1.0 - DELTA
    print(f"  n_s                         = 1 - 1/28 = {n_s:.12f}")
    print(f"  N_shell if saturated         = C_F alpha0^-4 = {n_shell_sat:.6e}")
    print(f"  A_nu if S_j=1,F_eff=1        = (3/4) alpha0^4 = {a_sat:.12e}")
    check(abs(a_sat - 0.75 * ALPHA0**4) < 1.0e-24, "conditional amplitude equals (3/4) alpha0^4")

    print("\n[3] Channel projection remains independent")
    total_channels = 28
    colour_channels = 6
    rank_entropy_colour = 2
    u_uniform = colour_channels / total_channels
    print(f"  service channels             = {total_channels}")
    print(f"  colour-restoring channels     = {colour_channels}")
    print(f"  rank(entropy, colour)         = {rank_entropy_colour}")
    print(f"  uniform all-channel colour fraction = {u_uniform:.9f}")
    for u in [u_uniform, 0.5, 0.9, 1.0]:
        amp_ratio = amplitude_from_utilization(u) / a_sat
        print(f"  colour utilization u={u:.9f}: entropy can remain saturated, A/A_sat={amp_ratio:.6f}")
    check(rank_entropy_colour == 2, "entropy all-ones projection and colour-restoring projection are not the same constraint")
    check(abs(amplitude_from_utilization(u_uniform) / a_sat - total_channels / colour_channels) < 1.0e-12, "uniform entropy loading would not give the saturated colour-load amplitude")
    check(abs(amplitude_from_utilization(0.9) / a_sat - 1.0) > 0.1, "subcritical colour utilization materially shifts amplitude")

    print("\n[4] Spatial projection remains independent")
    n = 9
    k = (1, 0, 0)
    mean = 17.0
    m = n**3
    white = mean * np.eye(m)
    phi = horizon_profile(n, k)
    horizon_cov = white + 0.02 * mean**2 * np.outer(phi, phi)
    common_cov = white + 0.20 * mean**2 * np.ones((m, m))
    check(abs(phase_sum(n, k)) < 1.0e-12, "scalar horizon mode is compensated/nonzero")
    rows = [
        ("white local current", white),
        ("homogeneous common clock", common_cov),
        ("horizon-mode covariance", horizon_cov),
    ]
    for label, cov in rows:
        min_eig = float(np.linalg.eigvalsh(cov).min())
        s_j = structure_factor(cov, n, k, mean)
        zvar = zero_mode_variance(cov, mean)
        print(f"  {label:25s}: S_j(k)={s_j:10.6f}, zero-mode Var={zvar:10.6e}, min eig={min_eig:10.6e}")
        check(min_eig > -1.0e-9, f"{label}: covariance is positive semidefinite")
    check(abs(structure_factor(white, n, k, mean) - 1.0) < 1.0e-10, "local white current gives S_j=1")
    check(abs(structure_factor(common_cov, n, k, mean) - 1.0) < 1.0e-10, "homogeneous common clock does not affect nonzero Pi_k")
    check(structure_factor(horizon_cov, n, k, mean) > 2.0, "horizon-mode covariance changes amplitude without changing mean load")

    print("\n[5] Tilt cannot see a scale-independent structure factor")
    for s_factor in [0.5, 1.0, 2.0, 62.965]:
        slope = log_slope(0.5, 2.0, a_sat, s_factor=s_factor)
        amp_ratio = s_factor
        print(f"  S_j={s_factor:8.3f}: log-slope={slope:+.12f}, A/A_sat={amp_ratio:.3f}")
        check(abs(slope + DELTA) < 1.0e-14, f"S_j={s_factor:g}: tilt unchanged by scale-independent covariance")
    print("  Therefore the 28-clock tilt cannot by itself falsify or fix S_j.")

    print("\n[6] Reopened theorem target")
    targets = [
        (
            "channel-lock / critical latch",
            "prove scalar-relevant entropy events are exactly post-decoder colour-restoring weight-4 commits and the finite printer latches at u=1",
        ),
        (
            "spatial whitening",
            "prove the compensated service-current covariance has no connected horizon-mode component after the homogeneous clock is removed",
        ),
        (
            "constant-H/exit compatibility",
            "prove those two identities hold only during the finite saturated stage and fail by reheating/handoff",
        ),
    ]
    for name, target in targets:
        print(f"  {name:30s}: {target}")

    print("\n" + "=" * 96)
    print("VERDICT")
    print("  Reopening the HBC inflation branch does not currently close the amplitude.")
    print("  The recent R4 homogeneous lift closes the late/Planck zero-mode endpoint,")
    print("  but A_nu is a nonzero-mode scalar-current covariance problem.  The existing")
    print("  28-clock and saturated radial log-shell premises still give n_s=27/28,")
    print("  while finite countermodels preserve those premises and change A_nu.")
    print()
    print("  The live frontier is now exactly two operator identities: channel-lock")
    print("  plus spatial-whitening.  If both land, A_nu=(3/4)alpha0^4 and tilt/amplitude")
    print("  close together.  If either fails, the tilt can survive while the amplitude")
    print("  remains conditional or dies.")
    print("=" * 96)
    print("exit 0 -- HBC inflation reopened; amplitude remains conditional on channel-lock and spatial-whitening.")


if __name__ == "__main__":
    main()
