#!/usr/bin/env python3
r"""ITEM 131: HBC scalar-printer amplitude observable audit.

Question
--------
Can the scalar amplitude be promoted to

    A_nu = (3/4) alpha_0^4

or is there a real horizon-mode scalar-current operator that changes it?

Result
------
Both statements can be made sharply.

Within the local, single-clock, saturated HBC scalar-printer algebra:

  * the scalar readout is the post-decoder colour-restoring topology current,
    not all-channel entropy;
  * the six admitted scalar/readout labels all carry C_F=4/3 load;
  * saturated constant-H queue balance gives
        lambda_shell = N_shell alpha_0^4 = C_F;
  * product-local CTMC service, fixed-total exchangeable quota, and homogeneous
    clock noise give S_j(k!=0)=1 for compensated scalar modes.

Therefore the local branch gives

    N_shell = (4/3) alpha_0^-4,
    A_nu = 1/N_shell = (3/4) alpha_0^4.

The escape hatch is also real but is new physics: a nonlocal horizon-mode
current

    L_{k,c} propto sum_x cos(k.x) J_x,
    L_{k,s} propto sum_x sin(k.x) J_x

adds a rank-two Fourier covariance supported at +-k.  It preserves the
homogeneous shell count because the profiles are compensated, but it changes
S_j(k) and hence A_nu.  A positive CPTP/Lindblad noise source can only raise
the amplitude; lowering it requires nonlocal feedback/quota suppression, not
an added positive jump source.  Either branch is outside the current local
HBC/QEC service algebra and would reopen the tilt/amplitude audit.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0
C_F = 4.0 / 3.0
A_OBS = 2.10e-9
SIGMA_A_OBS = 0.03e-9

Point = tuple[int, int, int]
Hyperplane = tuple[int, int, int, int]
Vector = tuple[int, int, int]

NAMES = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")
COLOUR_BITS = {3, 4}
Q3_POINTS = tuple(range(8))
Q3_EDGES = tuple(
    (u, v)
    for u in Q3_POINTS
    for v in Q3_POINTS
    if u < v and int.bit_count(u ^ v) == 1
)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def dot_mod2(a: int, x: int) -> int:
    return int.bit_count(a & x) & 1


def hyperplanes() -> list[frozenset[int]]:
    out: set[frozenset[int]] = set()
    for normal in range(1, 8):
        for offset in (0, 1):
            out.add(frozenset(p for p in Q3_POINTS if dot_mod2(normal, p) == offset))
    return sorted(out, key=lambda h: tuple(sorted(h)))


def tripped_edges(support: frozenset[int]) -> int:
    return sum((u in support) ^ (v in support) for u, v in Q3_EDGES)


def support_name(support: frozenset[int]) -> str:
    return "{" + ",".join(NAMES[i] for i in sorted(support)) + "}"


@dataclass(frozen=True)
class ChannelRow:
    support: frozenset[int]
    scalar_readout: bool
    tripped: int
    load: Fraction


def finite_channel_rows() -> list[ChannelRow]:
    rows = []
    for support in hyperplanes():
        tripped = tripped_edges(support)
        rows.append(
            ChannelRow(
                support=support,
                scalar_readout=not bool(support & COLOUR_BITS),
                tripped=tripped,
                load=Fraction(2 * tripped, len(Q3_EDGES)),
            )
        )
    return rows


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))  # type: ignore[return-value]


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def phase_vector(n: int, k: Vector) -> np.ndarray:
    return np.array([phase(n, k, x) for x in torus_points(n)], dtype=complex)


def cos_profile(n: int, k: Vector) -> np.ndarray:
    return np.array([math.cos(2 * math.pi * dot(k, x) / n) for x in torus_points(n)], dtype=float)


def sin_profile(n: int, k: Vector) -> np.ndarray:
    return np.array([math.sin(2 * math.pi * dot(k, x) / n) for x in torus_points(n)], dtype=float)


def structure_factor(cov: np.ndarray, n: int, k: Vector, mean_per_site: float) -> float:
    chars = phase_vector(n, k)
    m = n**3
    variance = (chars @ cov @ np.conjugate(chars)) / (m**2 * mean_per_site**2)
    return float(variance.real * m * mean_per_site)


def zero_mode_variance(cov: np.ndarray, mean_per_site: float) -> float:
    ones = np.ones(cov.shape[0])
    return float((ones @ cov @ ones) / (cov.shape[0] ** 2 * mean_per_site**2))


def translation_permutation(n: int, shift: Vector) -> np.ndarray:
    pts = torus_points(n)
    index = {x: i for i, x in enumerate(pts)}
    perm = np.zeros((len(pts), len(pts)))
    for i, x in enumerate(pts):
        y = tuple((xi + si) % n for xi, si in zip(x, shift))
        perm[index[y], i] = 1.0
    return perm


def horizon_covariance(n: int, k: Vector, mean_per_site: float, kappa: float) -> np.ndarray:
    """Translation-covariant nonlocal +-k covariance with S_j(k)=1+kappa.

    The rank-two Fourier pair is normalized so kappa is directly the excess
    scalar structure factor at the selected horizon mode.
    """

    m = n**3
    c = cos_profile(n, k)
    s = sin_profile(n, k)
    white = mean_per_site * np.eye(m)
    strength = 2.0 * kappa / (m * mean_per_site)
    return white + strength * mean_per_site**2 * (np.outer(c, c) + np.outer(s, s))


def horizon_suppression_covariance(n: int, k: Vector, mean_per_site: float, kappa: float) -> np.ndarray:
    """Nonlocal feedback/quota covariance with S_j(k)=1-kappa, if PSD."""

    m = n**3
    c = cos_profile(n, k)
    s = sin_profile(n, k)
    white = mean_per_site * np.eye(m)
    strength = 2.0 * kappa / (m * mean_per_site)
    return white - strength * mean_per_site**2 * (np.outer(c, c) + np.outer(s, s))


def local_covariances(n: int, mean_per_site: float) -> dict[str, np.ndarray]:
    m = n**3
    white = mean_per_site * np.eye(m)
    fixed_total = mean_per_site * (np.eye(m) - np.ones((m, m)) / m)
    common_clock = white + 0.2 * mean_per_site**2 * np.ones((m, m))
    return {
        "product-local CTMC": white,
        "fixed-total shell quota": fixed_total,
        "homogeneous common clock": common_clock,
    }


def offdiag_density(mat: np.ndarray, tol: float = 1e-14) -> float:
    off = mat.copy()
    np.fill_diagonal(off, 0.0)
    return float(np.count_nonzero(np.abs(off) > tol) / off.size)


def amplitude_from_sj(s_j: float) -> float:
    n_shell = C_F * ALPHA0 ** -4
    return s_j / n_shell


def main() -> None:
    print("ITEM 131 HBC SCALAR-PRINTER OBSERVABLE AUDIT")
    print("=" * 104)

    require_text(
        "recent_papers/cosmology/cosmology.tex",
        (
            "The scalar readout is no longer arbitrary",
            "Fourier-pair horizon current",
            "Within the local single-clock scalar-printer algebra",
        ),
    )
    require_text(
        "python_code/item131_hbc_stop_rule_proof.py",
        (
            "lambda_shell = E[N_topological commits per printed shell] = C_F = 4/3",
            "the proof uses the HBC scalar-printer identification",
        ),
    )

    print("\n[1] Channel-lock / stop-rule calculation")
    rows = finite_channel_rows()
    scalar = [row for row in rows if row.scalar_readout]
    kernel = [row for row in rows if not row.scalar_readout]
    scalar_loads = {row.load for row in scalar}
    print(f"  hyperplane supports             = {len(rows)}")
    print(f"  transverse service labels        = {2 * len(rows)}")
    print(f"  scalar/readout labels            = {2 * len(scalar)}")
    print(f"  scalar-kernel labels             = {2 * len(kernel)}")
    for row in rows:
        tag = "scalar/readout" if row.scalar_readout else "readout-kernel"
        print(f"    {support_name(row.support):24s} {tag:14s} trips={row.tripped:2d}/12 load={row.load}")
    check(len(rows) == 14, "there are 14 weight-4 hyperplane supports")
    check(2 * len(rows) == 28, "the transverse service instrument has 28 labels")
    check(2 * len(scalar) == 6, "the gauge-projected scalar readout has six labels")
    check(scalar_loads == {Fraction(4, 3)}, "every scalar/readout label carries C_F=4/3")

    n_shell = C_F * ALPHA0 ** -4
    a_local = amplitude_from_sj(1.0)
    z_obs = (a_local - A_OBS) / SIGMA_A_OBS
    print("\n[2] Local saturated-printer amplitude")
    print("  queue balance: lambda_shell = N_shell alpha_0^4 = C_F")
    print(f"  N_shell                         = {n_shell:.9e}")
    print(f"  A_local = 1/N_shell             = {a_local:.12e}")
    print(f"          = (3/4) alpha_0^4")
    print(f"  comparison to A_s=2.10(3)e-9    = {z_obs:+.3f} sigma")
    check(abs(a_local - 0.75 * ALPHA0**4) < 1e-24, "local branch gives A_nu=(3/4)alpha_0^4")
    check(abs(z_obs) < 1.2, "local branch is inside the current one-sigma amplitude band")

    print("\n[3] Spatial-whitening calculation")
    n = 9
    k = (1, 0, 0)
    mean = 17.0
    m = n**3
    chars = phase_vector(n, k)
    print(f"  toy horizon lattice              = {n}^3 sites")
    print(f"  compensated mode sum             = {abs(np.sum(chars)):.3e}")
    for label, cov in local_covariances(n, mean).items():
        eig_min = float(np.linalg.eigvalsh(cov).min())
        s_j = structure_factor(cov, n, k, mean)
        zvar = zero_mode_variance(cov, mean)
        print(f"  {label:28s}: S_j={s_j:10.6f} zero-var={zvar:10.6e} min-eig={eig_min:10.6e}")
        check(eig_min > -1e-8, f"{label}: covariance is positive semidefinite")
        check(abs(s_j - 1.0) < 1e-9, f"{label}: compensated nonzero mode has S_j=1")

    print("\n[4] Nonlocal horizon-mode scalar-current operator")
    white_zero_var = zero_mode_variance(local_covariances(n, mean)["product-local CTMC"], mean)
    for kappa in (0.2, 1.0, 61.965):
        cov = horizon_covariance(n, k, mean, kappa)
        eig_min = float(np.linalg.eigvalsh(cov).min())
        s_j = structure_factor(cov, n, k, mean)
        zvar = zero_mode_variance(cov, mean)
        comm = np.linalg.norm(translation_permutation(n, (1, 0, 0)) @ cov @ translation_permutation(n, (1, 0, 0)).T - cov)
        print(
            f"  L_k Fourier pair kappa={kappa:7.3f}: "
            f"S_j={s_j:10.6f} zero-var={zvar:10.6e} "
            f"offdiag={offdiag_density(cov):.3f} T-comm={comm:.3e}"
        )
        check(eig_min > -1e-8, f"kappa={kappa:g}: positive covariance / CPTP noise source")
        check(abs(s_j - (1.0 + kappa)) < 1e-9, f"kappa={kappa:g}: S_j=1+kappa")
        check(abs(zvar - white_zero_var) < 1e-12, f"kappa={kappa:g}: zero-mode shell variance is unchanged from the white baseline")
        check(comm < 1e-9, f"kappa={kappa:g}: covariance is translation-covariant")

    print("\n[5] Can the nonlocal operator tune the observed central value?")
    target_ratio = A_OBS / a_local
    positive_kappa = target_ratio - 1.0
    suppression_kappa = 1.0 - target_ratio
    print(f"  observed/local amplitude ratio   = {target_ratio:.6f}")
    print(f"  positive Lindblad source kappa    = {positive_kappa:+.6f}")
    print(f"  suppression quota kappa           = {suppression_kappa:+.6f}")
    check(positive_kappa < 0.0, "the observed central value is slightly below the local branch")
    print("  A positive extra horizon-mode jump source raises A_nu, so it cannot tune")
    print("  the slightly lower central value.  A nonlocal feedback/quota suppressor can.")
    suppress = horizon_suppression_covariance(n, k, mean, suppression_kappa)
    s_suppress = structure_factor(suppress, n, k, mean)
    eig_suppress = float(np.linalg.eigvalsh(suppress).min())
    print(f"  nonlocal suppressor: S_j={s_suppress:.6f}, A={amplitude_from_sj(s_suppress):.12e}, min-eig={eig_suppress:.6e}")
    check(abs(s_suppress - target_ratio) < 1e-9, "mode suppressor can set the central amplitude")
    check(eig_suppress > -1e-8, "small horizon-mode suppressor remains PSD")

    print("\n[6] Decision")
    print(
        """
  The calculation supports the current Canon wording:

    * A_s=(3/4) alpha_0^4 is promoted as far as the local single-clock
      saturated scalar-printer algebra permits.
    * It is not unconditional, because a genuine horizon-mode scalar-current
      operator exists mathematically and changes S_j(k=aH) without changing
      the homogeneous shell count.
    * That operator is nonlocal over the horizon shell.  It is not generated
      by local fresh-ancilla service, fixed-total exchangeable allocation, or
      homogeneous clock noise.

  Therefore the remaining theorem is precise: either prove that the HBC/QEC
  scalar service generator contains only local plus homogeneous-clock current
  terms after compensation, or admit the Fourier-pair horizon current as a new
  scalar source and rerun both the tilt and amplitude audits.
"""
    )
    print("ALL ASSERTIONS PASSED -- local amplitude conditionally closes; nonlocal horizon-current escape is real new physics.")


if __name__ == "__main__":
    main()
