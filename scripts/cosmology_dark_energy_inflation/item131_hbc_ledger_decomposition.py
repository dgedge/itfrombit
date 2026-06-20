#!/usr/bin/env python3
r"""ITEM 131: HBC entropy/colour/covariance ledger decomposition.

This is the deeper version of the HBC saturation gate.

Question
--------
Why exactly does

    HBC entropy/area saturation

not force

    colour-load saturation,  N_shell alpha_0^4 = C_F,

and why does colour-load saturation not force

    S_j(k=aH)=1 ?

Answer
------
Because they constrain different projections of the microscopic service
ledger.

1. Channel projection:

   The homogeneous entropy/area ledger fixes the all-ones projection over the
   28 service channels.  The colour-restoring topology ledger fixes a different
   projection: the post-decoder colour-silent weight-4 channels.  These vectors
   are not proportional.  Therefore fixing the entropy projection leaves a
   one-parameter family of colour utilizations.

2. Spatial projection:

   Colour-load saturation fixes the spatial zero mode of the topology current.
   The scalar amplitude samples a nonzero compensated Fourier mode through
   S_j(k=aH).  The zero-mode constraint does not fix nonzero-mode covariance.
   A positive semidefinite covariance component at the horizon wave number
   changes S_j while leaving the mean load and the entropy/area ledger intact.

Thus the theorem needed to close Item 131 is not "more saturation prose" but a
pair of operator statements:

   A. channel-lock theorem: the scalar-relevant entropy projection is identical
      to the post-decoder colour-restoring current, and the finite printer
      latches at its colour capacity;
   B. spatial-whitening theorem: after the homogeneous shell clock is removed,
      the HBC service current is local/exchangeable, with no connected
      covariance at k=aH.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0
C_F = 4.0 / 3.0


Point = tuple[int, int, int]
Hyperplane = tuple[int, int, int, int]
Channel = tuple[Hyperplane, int]
Vector = tuple[int, int, int]


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def dot_bits(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x & y for x, y in zip(a, b)) & 1


def points() -> list[Point]:
    return list(product([0, 1], repeat=3))  # type: ignore[return-value]


def hyperplanes() -> list[Hyperplane]:
    out: list[Hyperplane] = []
    for n in points():
        if n == (0, 0, 0):
            continue
        for b in (0, 1):
            out.append((*n, b))
    return out


def contains_point(h: Hyperplane, p: Point) -> bool:
    return dot_bits(h[:3], p) == h[3]


def support(h: Hyperplane) -> set[int]:
    return {
        (p[0] << 2) | (p[1] << 1) | p[2]
        for p in points()
        if contains_point(h, p)
    }


def all_channels() -> list[Channel]:
    return [(h, mode) for h in hyperplanes() for mode in (0, 1)]


def colour_silent_channel(ch: Channel) -> bool:
    """Colour-restoring channel criterion inherited from the C_F audit.

    In the 8-bit register convention used by the item-131 C_F script,
    C0,C1 are bit positions 3,4.  The admitted colour-silent hyperplanes are
    those whose support avoids those two colour bits.  The transverse mode does
    not change this channel-level support criterion, so it doubles the count.
    """
    h, _mode = ch
    return not bool(support(h) & {3, 4})


def rate_vector(colour_utilization: float) -> np.ndarray:
    """Entropy-normalized 28-channel rate vector.

    Sum_i r_i = 1 always.  The fraction on the colour-restoring subspace is
    `colour_utilization`.  This is a witness family, not a physical proposal.
    """
    channels = all_channels()
    colour = np.array([colour_silent_channel(ch) for ch in channels], dtype=bool)
    r = np.zeros(len(channels), dtype=float)
    r[colour] = colour_utilization / np.sum(colour)
    r[~colour] = (1.0 - colour_utilization) / np.sum(~colour)
    return r


def rank_of_constraints() -> int:
    channels = all_channels()
    entropy = np.ones(len(channels), dtype=float)
    colour = np.array([1.0 if colour_silent_channel(ch) else 0.0 for ch in channels], dtype=float)
    return int(np.linalg.matrix_rank(np.vstack([entropy, colour])))


def amplitude_from_utilization(colour_utilization: float, s_j: float = 1.0) -> float:
    return s_j / (colour_utilization * C_F / ALPHA0**4)


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))  # type: ignore[return-value]


def dot_int(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot_int(k, x) / n)


def phase_sum(n: int, k: Vector) -> complex:
    return sum(phase(n, k, x) for x in torus_points(n))


def horizon_profile(n: int, k: Vector) -> np.ndarray:
    return np.array([phase(n, k, x).real for x in torus_points(n)], dtype=float)


def structure_factor_from_covariance(cov: np.ndarray, n: int, k: Vector, mean_per_site: float) -> float:
    pts = torus_points(n)
    m = len(pts)
    chars = np.array([phase(n, k, x) for x in pts], dtype=complex)
    variance = (chars @ cov @ np.conjugate(chars)) / (m**2 * mean_per_site**2)
    n_shell = m * mean_per_site
    return float(variance.real * n_shell)


def zero_mode_variance(cov: np.ndarray, mean_per_site: float) -> float:
    m = cov.shape[0]
    ones = np.ones(m)
    return float((ones @ cov @ ones) / (m**2 * mean_per_site**2))


@dataclass(frozen=True)
class SpatialCase:
    name: str
    cov: np.ndarray
    note: str


def spatial_cases(n: int, mean_per_site: float, k: Vector) -> list[SpatialCase]:
    pts = torus_points(n)
    m = len(pts)
    identity = mean_per_site * np.eye(m)
    phi = horizon_profile(n, k)
    phi = phi - np.mean(phi)
    horizon_cov = 0.02 * mean_per_site**2 * np.outer(phi, phi)
    common = 0.20 * mean_per_site**2 * np.ones((m, m))
    return [
        SpatialCase("local Poisson", identity, "white local current"),
        SpatialCase("homogeneous common-rate", identity + common, "changes only the zero mode"),
        SpatialCase("horizon-mode covariance", identity + horizon_cov, "leaves mean load fixed but changes k=aH"),
    ]


def main() -> None:
    print("ITEM 131 HBC LEDGER DECOMPOSITION")

    print("\n[1] Source checks")
    checks = [
        ("28-channel instrument", "python_code/item131_w_to_28_instrument.py", "service channel has rate 1/28"),
        ("C_F colour support", "python_code/item131_cf_stop_rule_closure_attempt.py", "colour-silent logical channels"),
        ("HBC entropy saturation", "cosmological_qec_engine/cosmological_qec_engine.tex", "Bekenstein--Hawking entropy bound saturation"),
        ("T5b structure factor", "python_code/item131_t5b_correlation_volume_audit.py", "S_j(k)"),
        ("HBC saturation gate", "python_code/item131_hbc_saturation_gate.py", "independence witnesses"),
    ]
    for label, path, phrase in checks:
        check(contains(path, phrase), label)

    print("\n[2] Channel projection: entropy saturation is not colour-load saturation")
    channels = all_channels()
    colour = np.array([colour_silent_channel(ch) for ch in channels], dtype=bool)
    print(f"  service channels               = {len(channels)}")
    print(f"  colour-restoring channels       = {int(np.sum(colour))}")
    print(f"  non-colour channels             = {int(np.sum(~colour))}")
    print(f"  rank{{entropy vector, colour vector}} = {rank_of_constraints()}")
    check(len(channels) == 28, "finite service bridge has 28 channels")
    check(int(np.sum(colour)) == 6, "three colour-silent hyperplanes times two modes give six colour-restoring channels")
    check(rank_of_constraints() == 2, "entropy and colour projections are linearly independent")

    a_sat = amplitude_from_utilization(1.0)
    for u in (0.25, 0.50, 0.90, 1.00):
        r = rate_vector(u)
        entropy_projection = float(np.sum(r))
        colour_projection = float(np.sum(r[colour]))
        amp = amplitude_from_utilization(u)
        print(
            f"  utilisation={u:4.2f}: entropy={entropy_projection:.3f}, "
            f"colour={colour_projection:.3f}, A/A_sat={amp/a_sat:.3f}"
        )
        check(abs(entropy_projection - 1.0) < 1e-12, f"u={u}: entropy/area projection remains saturated")
        check(abs(colour_projection - u) < 1e-12, f"u={u}: colour projection varies independently")
    print("  Algebraic meaning: HBC area saturation fixes the all-ones projection.")
    print("  It does not identify that projection with the colour-restoring subspace.")

    print("\n[3] What extra channel-lock theorem must prove")
    print("  Needed map:")
    print("      entropy event -> post-decoder colour-restoring weight-4 commit")
    print("  plus a finite latch:")
    print("      coherent constant-H printing maximizes this current until")
    print("      lambda_shell = C_F, then exits for overload.")
    print("  Without that, a stable subcritical queue is algebraically allowed.")

    print("\n[4] Spatial projection: colour saturation is not S_j(k)=1")
    n = 9
    k = (1, 0, 0)
    mean = 17.0
    check(abs(phase_sum(n, k)) < 1e-12, "the scalar mode is compensated: sum_x e^{-ikx}=0")
    for case in spatial_cases(n, mean, k):
        eig_min = float(np.linalg.eigvalsh(case.cov).min())
        s_j = structure_factor_from_covariance(case.cov, n, k, mean)
        zvar = zero_mode_variance(case.cov, mean)
        print(
            f"  {case.name:25s}: S_j(k)={s_j:10.6f}, "
            f"zero-mode Var={zvar:10.6e}, min eig={eig_min:10.6e}  {case.note}"
        )
        check(eig_min > -1e-9, f"{case.name}: covariance is positive semidefinite")
    print("  The homogeneous ledger reads the mean and zero mode.  The scalar")
    print("  amplitude reads the nonzero horizon Fourier covariance.")
    print("  These are independent covariance components.")

    print("\n[5] Consequence for Item 131")
    conditions = [
        ("entropy-area saturation", "fixes total printed capacity / homogeneous shell clock", "not enough for A_nu"),
        ("colour-load saturation", "fixes N_shell alpha0^4=C_F", "not enough for S_j=1"),
        ("spatial whiteness", "fixes S_j(k=aH)=1", "not enough for N_shell"),
        ("all three", "fix N_shell and S_j with the 28-clock", "closes A_nu and n_s together"),
    ]
    for name, fixes, status in conditions:
        print(f"  {name:24s}: {fixes:56s} -> {status}")

    print("\nVERDICT")
    print("  The two caveats are not loose worries; they are independent projections.")
    print("  Entropy saturation is a channel zero-mode constraint.  Colour-load")
    print("  saturation is a different channel projection.  S_j(k=aH) is a spatial")
    print("  nonzero-mode covariance.  Existing HBC/QEC mechanics fix useful pieces")
    print("  of each, but not the operator identities that collapse them into one.")
    print("  The precise route forward is now:")
    print("    1. prove the channel-lock/latch theorem,")
    print("    2. prove the spatial-whitening theorem.")
    print("exit 0 -- entropy, colour load, and scalar covariance are independent ledger projections.")


if __name__ == "__main__":
    main()
