#!/usr/bin/env python3
r"""Mesoscopic Alexandrov-loop Maxwell action: candidate gates, not closure.

K23 isolated the viable manifestly gauge-invariant route:

    mesoscopic Alexandrov-loop holonomy^2,
    divided by a cardinality-estimated area,
    averaged over a covariant interval-size window.

This script pushes that candidate through the gates K23 named:

  1. non-arbitrary interval measure,
  2. 3+1 tensor contraction,
  3. locality expansion,
  4. overlap-correlation control,
  5. normalization.

The result is deliberately tiered.  The constant-field / Wick-rotated 3+1
tensor and normalization gates close as a computable scaffold; the full
Lorentzian causal-loop continuum theorem is still open.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


RNG = np.random.default_rng(20260613)
DIM = 4
N_BIVECTOR = DIM * (DIM - 1) // 2


def random_antisymmetric(dim: int = DIM) -> np.ndarray:
    a = RNG.normal(size=(dim, dim))
    return a - a.T


def upper_components(a: np.ndarray) -> np.ndarray:
    return np.array([a[i, j] for i in range(a.shape[0]) for j in range(i + 1, a.shape[1])])


def random_simple_bivector(dim: int = DIM) -> np.ndarray:
    """Unit simple bivector B=u wedge v, represented by upper-triangle components."""

    q, _ = np.linalg.qr(RNG.normal(size=(dim, 2)))
    u, v = q[:, 0], q[:, 1]
    b = np.outer(u, v) - np.outer(v, u)
    return upper_components(b)


def alexandrov_volume_4(tau: float) -> float:
    """Continuum 3+1 Alexandrov volume V_4(tau)=pi tau^4/24."""

    return math.pi * tau**4 / 24.0


def tau_from_cardinality(n_eps: float, rho: float = 1.0) -> float:
    """Order-only cardinality to proper-time scale for a 3+1 interval."""

    return (24.0 * n_eps / (math.pi * rho)) ** 0.25


@dataclass(frozen=True)
class IntervalWindow:
    n_target: float
    rho: float = 1.0
    rel_width: float = 0.15

    @property
    def tau(self) -> float:
        return tau_from_cardinality(self.n_target, self.rho)

    @property
    def lo(self) -> int:
        return math.ceil((1.0 - self.rel_width) * self.n_target)

    @property
    def hi(self) -> int:
        return math.floor((1.0 + self.rel_width) * self.n_target)


def interval_window_stats(window: IntervalWindow, trials: int = 100_000) -> tuple[float, float, float]:
    counts = RNG.poisson(window.n_target, trials)
    keep = counts[(counts >= window.lo) & (counts <= window.hi)]
    tau_hat = tau_from_cardinality(keep, window.rho)
    return float(len(keep) / trials), float(np.mean(tau_hat / window.tau)), float(np.std(tau_hat / window.tau))


def tensor_contraction_gate(samples: int = 120_000) -> tuple[float, float, float]:
    """Verify 3+1 isotropic bivector average.

    For a random unit simple bivector B in four Euclidean/Wick dimensions,

        E[(F.B)^2] = ||F||^2 / 6,

    where ||F||^2=sum_{mu<nu} F_{mu nu}^2.  Therefore the Euclidean Maxwell
    density 1/4 F_{mu nu}F_{mu nu}=1/2||F||^2 is

        3 E[(F.B)^2].
    """

    f = upper_components(random_antisymmetric())
    f_norm2 = float(np.dot(f, f))
    dots = np.array([np.dot(f, random_simple_bivector()) for _ in range(samples)])
    maxwell_true = 0.5 * f_norm2
    maxwell_est = 3.0 * float(np.mean(dots * dots))
    rel_err = (maxwell_est / maxwell_true) - 1.0
    return maxwell_true, maxwell_est, rel_err


def locality_gate() -> tuple[float, float, float]:
    """Symmetric mesoscopic loop cancels odd gradients; first error is O(tau^2).

    Model one projected field component on a loop plane:
        f(s,t)=f0+a_s s+a_t t+b(s^2+t^2).
    The symmetric surface average cancels the linear part exactly and differs
    from f0 by b h^2/6 for a square loop of side h.
    """

    f0, a_s, a_t, b = 1.7, 3.0, -2.0, 0.43
    hs = np.array([0.4, 0.2, 0.1, 0.05])
    errors = []
    for h in hs:
        grid = np.linspace(-h / 2.0, h / 2.0, 801)
        s, t = np.meshgrid(grid, grid, indexing="ij")
        values = f0 + a_s * s + a_t * t + b * (s * s + t * t)
        errors.append(float(np.mean(values) - f0))
    errors = np.array(errors)
    slope = np.polyfit(np.log(hs), np.log(np.abs(errors)), 1)[0]
    coeff = float(errors[-1] / (hs[-1] ** 2))
    expected = b / 6.0
    return float(slope), coeff, expected


def overlap_control_gate(cluster_size: int = 8, clusters: int = 1200, rho: float = 0.55) -> tuple[float, float, float]:
    """Block-estimate the effective sample size for overlapping interval families.

    Intervals with centers closer than O(tau) overlap and are not independent.
    This synthetic test models each local overlap cluster with equicorrelation
    rho and checks that block averaging recovers the inflated error bar:

        M_eff = M / (1 + (m-1)rho).
    """

    sigma = 1.0
    common = RNG.normal(scale=math.sqrt(rho) * sigma, size=(clusters, 1))
    private = RNG.normal(scale=math.sqrt(1.0 - rho) * sigma, size=(clusters, cluster_size))
    x = common + private
    naive_se = float(np.std(x.reshape(-1), ddof=1) / math.sqrt(x.size))
    block_means = np.mean(x, axis=1)
    block_se = float(np.std(block_means, ddof=1) / math.sqrt(clusters))
    inflation = block_se / naive_se
    expected = math.sqrt(1.0 + (cluster_size - 1) * rho)
    return inflation, expected, block_se


def main() -> None:
    print("[1] Non-arbitrary interval measure")
    window = IntervalWindow(n_target=4096.0, rho=1.0, rel_width=0.10)
    accept, mean_tau, std_tau = interval_window_stats(window)
    print("    measure: ordered pairs whose Alexandrov interval cardinality lies in a fixed")
    print("             covariant window N_epsilon +/- 10%; tau from V_4=tau^4*pi/24.")
    print(f"    N_epsilon={window.n_target:.0f}, tau={window.tau:.6f}, acceptance={accept:.3f}")
    print(f"    tau_hat/tau = {mean_tau:.6f} +/- {std_tau:.6f}")
    assert accept > 0.99 and abs(mean_tau - 1.0) < 0.003 and std_tau < 0.006

    print("\n[2] 3+1 tensor contraction and normalization")
    true, est, rel = tensor_contraction_gate()
    print("    Wick/Euclidean 4D isotropic simple-bivector average:")
    print("        E[(F.B)^2]=||F||^2/6, so 1/4 F_mn F_mn = 3 E[(F.B)^2].")
    print(f"    Maxwell density true={true:.6f}, estimator={est:.6f}, rel_err={rel:+.4%}")
    assert abs(rel) < 0.015

    print("\n[3] Locality expansion")
    slope, coeff, expected = locality_gate()
    print("    symmetric loop cancels the linear-gradient term; first correction is O(tau^2).")
    print(f"    fitted log-error slope={slope:.4f}, coefficient={coeff:.6f}, expected={expected:.6f}")
    assert abs(slope - 2.0) < 0.03 and abs(coeff - expected) < 0.002

    print("\n[4] Overlap-correlation control")
    infl, expected, block_se = overlap_control_gate()
    print("    overlapping intervals must be blocked by center-separation/tau cells.")
    print(f"    measured SE inflation={infl:.4f}, expected={expected:.4f}, block_SE={block_se:.5f}")
    assert abs(infl / expected - 1.0) < 0.05

    print("\n[5] Verdict")
    print("    The mesoscopic Alexandrov-loop route now has an executable scaffold:")
    print("      * interval-size measure: fixed cardinality window, order/covariant;")
    print("      * 3+1 tensor normalization: isotropic bivector average fixed;")
    print("      * locality: symmetric loop gives F(x)+O(tau^2 partial^2 F);")
    print("      * correlations: overlap handled by block/effective-sample counting.")
    print("    This is NOT a full Maxwell loop-action theorem: the remaining proof is the")
    print("    Lorentzian causal-loop construction and continuum limit for variable fields.")
    print("ALL ASSERTIONS PASSED -- mesoscopic loop action candidate sharpened, not locked.")


if __name__ == "__main__":
    main()
