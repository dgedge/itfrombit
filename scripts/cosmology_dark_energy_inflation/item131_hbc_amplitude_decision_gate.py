#!/usr/bin/env python3
r"""ITEM 131: HBC scalar-amplitude decision gate.

Question
--------
Can the current HBC/QEC ledger prove both

    N_shell alpha_0^4 = C_F = 4/3
    S_j(k=aH) = 1

so that A_nu=(3/4) alpha_0^4 is Locked rather than conditional?

Answer
------
No, not from the current premises.  This script proves the obstruction by
constructing finite admissible ledgers that preserve every already-closed
ingredient but give different scalar amplitudes.

Closed ingredients preserved by the witnesses:
  * the 28-channel service instrument;
  * the post-decoder weight-4 / alpha_0^4 topology-current scale;
  * the C_F=4/3 colour-restoring load coefficient;
  * the compensated Fourier scalar projector Pi_k;
  * F_eff=1 for the CTMC count ledger;
  * homogeneous entropy/area saturation.

Two independent freedoms remain:
  1. channel freedom: the all-ones entropy projection is not the
     colour-restoring projection, so entropy saturation does not force
     colour-load utilization u=1;
  2. spatial freedom: the zero-mode load constraint does not fix connected
     covariance at the nonzero horizon mode k=aH, so colour saturation does
     not force S_j(k)=1.

Therefore the scalar amplitude remains a conditional recovery candidate until
two new theorems are supplied: a channel-lock/critical-latch theorem and a
spatial-whitening theorem.
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


def bitdot(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    return sum(x & y for x, y in zip(a, b)) & 1


def q3_points() -> list[Point]:
    return list(product((0, 1), repeat=3))  # type: ignore[return-value]


def hyperplanes() -> list[Hyperplane]:
    out: list[Hyperplane] = []
    for normal in q3_points():
        if normal == (0, 0, 0):
            continue
        for offset in (0, 1):
            out.append((*normal, offset))
    return out


def hyperplane_support(h: Hyperplane) -> set[int]:
    normal = h[:3]
    offset = h[3]
    return {
        (p[0] << 2) | (p[1] << 1) | p[2]
        for p in q3_points()
        if bitdot(normal, p) == offset
    }


def channels() -> list[Channel]:
    return [(h, mode) for h in hyperplanes() for mode in (0, 1)]


def colour_restoring(ch: Channel) -> bool:
    """Confinement-filtered colour-silent channel from the C_F audit.

    In the established register ordering, C0,C1 are vertex labels 3 and 4.
    The admitted channels are those whose hyperplane support avoids both.
    """

    h, _mode = ch
    return len(hyperplane_support(h) & {3, 4}) == 0


def channel_rate_vector(utilization: float) -> np.ndarray:
    """Entropy-normalized 28-channel rate witness.

    The all-ones entropy projection is fixed to one for every utilization.
    The colour-restoring projection equals utilization.
    """

    ch = channels()
    c = np.array([colour_restoring(x) for x in ch], dtype=bool)
    rates = np.zeros(len(ch), dtype=float)
    rates[c] = utilization / np.sum(c)
    rates[~c] = (1.0 - utilization) / np.sum(~c)
    return rates


def amplitude(utilization: float, s_j: float = 1.0) -> float:
    return s_j / (utilization * C_F / ALPHA0**4)


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))  # type: ignore[return-value]


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def phase_sum(n: int, k: Vector) -> complex:
    return sum(phase(n, k, x) for x in torus_points(n))


def profile(n: int, k: Vector) -> np.ndarray:
    values = np.array([phase(n, k, x).real for x in torus_points(n)], dtype=float)
    return values - np.mean(values)


def structure_factor(cov: np.ndarray, n: int, k: Vector, mean: float) -> float:
    pts = torus_points(n)
    chars = np.array([phase(n, k, x) for x in pts], dtype=complex)
    variance = (chars @ cov @ np.conjugate(chars)) / ((n**3) ** 2 * mean**2)
    return float((variance.real) * (n**3) * mean)


def zero_mode_load(cov: np.ndarray, mean: float) -> float:
    ones = np.ones(cov.shape[0])
    return float((ones @ cov @ ones) / (cov.shape[0] ** 2 * mean**2))


@dataclass(frozen=True)
class SpatialWitness:
    name: str
    cov: np.ndarray
    expected_status: str


def spatial_witnesses(n: int, k: Vector, mean: float) -> list[SpatialWitness]:
    m = n**3
    white = mean * np.eye(m)
    common = 0.1 * mean**2 * np.ones((m, m))
    phi = profile(n, k)
    horizon = 0.02 * mean**2 * np.outer(phi, phi)
    return [
        SpatialWitness("white local CTMC", white, "target S_j=1"),
        SpatialWitness("homogeneous common clock", white + common, "same nonzero S_j"),
        SpatialWitness("horizon-mode covariance", white + horizon, "allowed countermodel"),
    ]


def main() -> None:
    print("ITEM 131 HBC SCALAR-AMPLITUDE DECISION GATE")
    print("=" * 96)

    print("\n[1] Existing owner gates")
    owner_checks = [
        ("28-channel service bridge", "python_code/item131_w_to_28_instrument.py", "service channel has rate 1/28"),
        ("C_F load geometry", "python_code/item131_cf_stop_rule_closure_attempt.py", "2*(8/12)=4/3"),
        ("projector form", "python_code/item131_scalar_mode_projector.py", "Pi_k projector form derived"),
        ("T5b whiteness conditional", "python_code/item131_t5b_whiteness_lemma.py", "S_j(k=aH)=1"),
        ("ledger decomposition", "python_code/item131_hbc_ledger_decomposition.py", "independent ledger projections"),
        ("saturation gate", "python_code/item131_hbc_saturation_gate.py", "explicit countermodels"),
    ]
    for label, path, phrase in owner_checks:
        check(contains(path, phrase), label)

    print("\n[2] Channel projection witness")
    ch = channels()
    c = np.array([colour_restoring(x) for x in ch], dtype=bool)
    entropy_vector = np.ones(len(ch))
    colour_vector = c.astype(float)
    rank = int(np.linalg.matrix_rank(np.vstack([entropy_vector, colour_vector])))
    print(f"  channels                         = {len(ch)}")
    print(f"  colour-restoring channels         = {int(np.sum(c))}")
    print(f"  rank(all-ones, colour-projection) = {rank}")
    check(len(ch) == 28, "finite service instrument has 28 channels")
    check(int(np.sum(c)) == 6, "colour-restoring subspace has 6 channels")
    check(rank == 2, "entropy and colour projections are independent")

    a_target = amplitude(1.0, 1.0)
    for utilization in (0.25, 0.50, 0.90, 1.00):
        rates = channel_rate_vector(utilization)
        entropy = float(np.sum(rates))
        colour = float(np.sum(rates[c]))
        amp_ratio = amplitude(utilization, 1.0) / a_target
        print(
            f"  utilization={utilization:4.2f}: "
            f"entropy={entropy:.3f}, colour={colour:.3f}, A/A_sat={amp_ratio:.3f}"
        )
        check(abs(entropy - 1.0) < 1e-12, f"u={utilization}: entropy saturation preserved")
        check(abs(colour - utilization) < 1e-12, f"u={utilization}: colour projection freely set")
    check(abs(amplitude(0.5) / a_target - 2.0) < 1e-12, "half colour-load doubles the scalar amplitude")
    print("  Therefore HBC entropy/area saturation alone cannot prove")
    print("  N_shell alpha_0^4=C_F.")

    print("\n[3] Spatial covariance witness")
    n = 9
    k = (1, 0, 0)
    mean = 17.0
    check(abs(phase_sum(n, k)) < 1.0e-12, "k=aH toy mode is compensated/nonzero")
    for witness in spatial_witnesses(n, k, mean):
        eig_min = float(np.linalg.eigvalsh(witness.cov).min())
        sj = structure_factor(witness.cov, n, k, mean)
        zload = zero_mode_load(witness.cov, mean)
        print(
            f"  {witness.name:25s}: S_j={sj:10.6f}, "
            f"zero-load={zload:10.6e}, min-eig={eig_min:10.6e}  {witness.expected_status}"
        )
        check(eig_min > -1.0e-9, f"{witness.name}: covariance is positive semidefinite")
    sj_white = structure_factor(spatial_witnesses(n, k, mean)[0].cov, n, k, mean)
    sj_horizon = structure_factor(spatial_witnesses(n, k, mean)[2].cov, n, k, mean)
    check(abs(sj_white - 1.0) < 1.0e-10, "white local CTMC gives S_j=1")
    check(sj_horizon > 2.0, "horizon-mode connected covariance changes S_j materially")
    print("  Therefore colour-load saturation/zero-mode constraints cannot prove")
    print("  S_j(k=aH)=1.")

    print("\n[4] Decision table")
    rows = [
        ("finite 28-clock", "CLOSED", "relative shell generator"),
        ("Pi_k scalar projector", "CLOSED", "compensated Fourier shell"),
        ("F_eff=1 CTMC", "CLOSED", "count Fano factor"),
        ("C_F=4/3 coefficient", "SUBSTANTIALLY CLOSED", "post-decoder colour-restoring load"),
        ("N_shell alpha^4=C_F", "NOT PROVED", "channel utilization remains free"),
        ("S_j(k=aH)=1", "NOT PROVED", "nonzero-mode covariance remains free"),
        ("A_nu=(3/4)alpha^4", "CONDITIONAL", "requires the two unproved identities"),
    ]
    for item, status, note in rows:
        print(f"  {item:28s} {status:22s} {note}")

    print(
        """
[5] THEOREM TARGET LEFT
    To promote the scalar amplitude, canon must prove two operator statements:

      A. Channel-lock / critical-latch theorem:
         scalar-relevant printed entropy is exactly the post-decoder
         colour-restoring topology current, and finite coherent printing
         latches at utilization u=1 rather than any subcritical u<1.

      B. Spatial-whitening theorem:
         after subtracting the homogeneous shell clock, HBC service-current
         covariance has no connected component at the horizon wave number.

    Without A and B, the same closed finite ledger admits different A_nu.
"""
    )
    print("ALL ASSERTIONS PASSED -- amplitude remains conditional under current canon.")


if __name__ == "__main__":
    main()
