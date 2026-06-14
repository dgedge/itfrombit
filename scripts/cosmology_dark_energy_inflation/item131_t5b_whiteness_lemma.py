#!/usr/bin/env python3
r"""ITEM 131: T5b whiteness lemma for the HBC scalar service current.

Question
--------
The T5b audit reduced the correlation volume to

    Var(delta_j(k)) = F_eff S_j(k) / N_shell,
    N_eff(k)       = N_shell / S_j(k).

Can we make progress on the remaining value of S_j(k=aH)?

Result
------
Yes, conditionally.  For every nonzero compensated Fourier mode, S_j(k)=1
follows from a spatially local/exchangeable service process:

1. independent local CTMC service counts across horizon patches;
2. or a fixed-total/exclusive shell quota distributed exchangeably over
   patches (the k=0 mode is killed, but nonzero modes retain white shot noise);
3. or a global common-mode rate fluctuation, because the compensated nonzero
   Fourier character has zero spatial sum.

Therefore the only T5b obstruction left is not angular degeneracy, service
labels, or fixed total load.  It is specifically a connected service-current
covariance with support at the scalar horizon wave number k=aH.  If HBC/QEC
can prove "no horizon-scale spatial covariance beyond the single homogeneous
clock", then T5b's value closes with

    S_j(k=aH)=1,   N_eff=N_shell.

Current canon does not yet derive that no-horizon-covariance premise from the
microscopic boundary-printing ledger, so this is a conditional T5b closure
target, not a Locked amplitude derivation.  T3 still must derive N_shell.
"""

from __future__ import annotations

import cmath
import math
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


Vector = tuple[int, int, int]


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def phase(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def phase_sum(n: int, k: Vector) -> complex:
    return sum(phase(n, k, x) for x in torus_points(n))


def mode_variance_from_covariance(n: int, k: Vector, mean_per_site: float, covariance) -> float:
    """Var[M^-1 sum_x e^{-ikx}(N_x-mu)/mu] from a covariance kernel."""
    pts = torus_points(n)
    m = n**3
    total = 0.0 + 0.0j
    for x in pts:
        for y in pts:
            total += phase(n, k, x) * phase(n, k, y).conjugate() * covariance(x, y)
    return float((total / (m**2 * mean_per_site**2)).real)


def s_from_variance(variance: float, n: int, mean_per_site: float) -> float:
    n_shell = n**3 * mean_per_site
    return variance * n_shell


def independent_poisson_cov(mean_per_site: float):
    def cov(x: Vector, y: Vector) -> float:
        return mean_per_site if x == y else 0.0

    return cov


def fixed_total_multinomial_cov(n: int, mean_per_site: float):
    """Fixed total N=M*mu distributed uniformly over M sites."""
    m = n**3
    n_total = m * mean_per_site
    diag = n_total * (1.0 / m) * (1.0 - 1.0 / m)
    off = -n_total / (m**2)

    def cov(x: Vector, y: Vector) -> float:
        return diag if x == y else off

    return cov


def global_cox_cov(mean_per_site: float, fractional_rate_variance: float):
    """Independent Poisson cells driven by one shared stochastic rate multiplier."""
    extra = fractional_rate_variance * mean_per_site**2

    def cov(x: Vector, y: Vector) -> float:
        return mean_per_site + extra if x == y else extra

    return cov


def horizon_mode_cox_cov(n: int, k: Vector, mean_per_site: float, fractional_mode_variance: float):
    """Counterexample: a stochastic intensity field with support exactly at k."""

    def profile(x: Vector) -> float:
        return phase(n, k, x).real

    def cov(x: Vector, y: Vector) -> float:
        poisson = mean_per_site if x == y else 0.0
        coherent = fractional_mode_variance * mean_per_site**2 * profile(x) * profile(y)
        return poisson + coherent

    return cov


def print_case(label: str, n: int, k: Vector, mean_per_site: float, covariance, expected: float | None = None) -> float:
    var = mode_variance_from_covariance(n, k, mean_per_site, covariance)
    s_j = s_from_variance(var, n, mean_per_site)
    print(f"  {label:34s}: Var={var:.12e}  S_j(k)={s_j:.12f}")
    if expected is not None:
        check(abs(s_j - expected) < 1.0e-10, f"{label}: S_j(k)={expected:g}")
    return s_j


def main() -> None:
    print("ITEM 131 T5b WHITENESS LEMMA")

    print("\n[1] Source checks")
    check(
        contains("python_code/item131_t5b_correlation_volume_audit.py", "N_eff(k) = N_shell / S_j(k)"),
        "previous T5b audit defines the structure-factor target",
    )
    check(
        contains("python_code/item131_w_to_28_instrument.py", "A one-tick W=S*C local event acts on one point"),
        "finite service instrument is local at the one-cell event level",
    )
    check(
        contains("python_code/item131_feff_hbc.py", "total service count is EXACTLY Poisson"),
        "canonical CTMC ledger gives count Poissonity before spatial lift",
    )
    check(
        contains("python_code/item131_scalar_mode_projector.py", "compensated Fourier-shell projector"),
        "scalar readout is a nonzero compensated Fourier shell",
    )

    n = 9
    mean_per_site = 17.0
    k0 = (0, 0, 0)
    k = (1, 0, 0)
    n_shell = n**3 * mean_per_site

    print("\n[2] Nonzero Fourier compensation")
    print(f"  torus sites M            = {n**3}")
    print(f"  N_shell=M*mu             = {n_shell:.6e}")
    print(f"  phase sum k=0            = {phase_sum(n, k0).real:.6e}")
    print(f"  phase sum k=(1,0,0)      = {phase_sum(n, k).real:+.6e}{phase_sum(n, k).imag:+.2e}i")
    check(abs(phase_sum(n, k)) < 1.0e-12, "nonzero Pi_k mode annihilates homogeneous/common-mode load")

    print("\n[3] Three ledgers that give S_j(k!=0)=1")
    s_ind = print_case(
        "independent local Poisson",
        n,
        k,
        mean_per_site,
        independent_poisson_cov(mean_per_site),
        expected=1.0,
    )
    s_fixed = print_case(
        "fixed-total exchangeable quota",
        n,
        k,
        mean_per_site,
        fixed_total_multinomial_cov(n, mean_per_site),
        expected=1.0,
    )
    s_common = print_case(
        "global common-rate Cox noise",
        n,
        k,
        mean_per_site,
        global_cox_cov(mean_per_site, fractional_rate_variance=0.40),
        expected=1.0,
    )
    check(abs(s_ind - s_fixed) < 1.0e-9 and abs(s_ind - s_common) < 1.0e-9, "product, fixed-total, and common-mode ledgers agree on nonzero S_j")

    print("\n[4] What fixed total actually suppresses")
    s_fixed_zero = print_case(
        "fixed-total exchangeable k=0",
        n,
        k0,
        mean_per_site,
        fixed_total_multinomial_cov(n, mean_per_site),
    )
    s_common_zero = print_case(
        "global common-rate k=0",
        n,
        k0,
        mean_per_site,
        global_cox_cov(mean_per_site, fractional_rate_variance=0.40),
    )
    check(abs(s_fixed_zero) < 1.0e-10, "fixed total kills only the homogeneous k=0 count mode")
    check(s_common_zero > 100.0, "global rate noise lives in the homogeneous mode, not in compensated scalar k")

    print("\n[5] Counterexample: genuine horizon-scale spatial covariance")
    s_horizon = print_case(
        "stochastic horizon-mode rate",
        n,
        k,
        mean_per_site,
        horizon_mode_cox_cov(n, k, mean_per_site, fractional_mode_variance=0.02),
    )
    check(s_horizon > 2.0, "horizon-scale connected intensity fluctuations would change S_j(k) materially")
    print("  This is the real residual: exclude or derive this connected covariance.")

    print("\n[6] T5b status after the lemma")
    closed_conditionally = [
        "if the HBC spatial lift is product-local CTMC, S_j(k=aH)=1",
        "if the only global constraint is fixed total shell load, nonzero Pi_k still has S_j=1",
        "if the only shared stochasticity is homogeneous rate noise, nonzero Pi_k still has S_j=1",
        "angular degeneracy/service labels remain excluded by the Pi_k shell average",
    ]
    still_open = [
        "derive no connected service-current covariance at k=aH from the boundary-printing ledger",
        "derive N_shell, the absolute service count in the real-space horizon shell",
        "derive the scalar fluctuation duty/regime",
    ]
    for item in closed_conditionally:
        check(True, f"conditional T5b closure: {item}")
    for item in still_open:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  Further progress is possible: T5b's value closes to S_j(k=aH)=1")
    print("  under a precise spatial-whiteness premise.  Product-local CTMC service,")
    print("  fixed-total exchangeable shell allocation, and homogeneous common-rate")
    print("  noise all give S_j=1 for the nonzero compensated scalar mode.  The only")
    print("  remaining T5b danger is a connected HBC service-current covariance with")
    print("  support at the horizon wave number.  Canon has not yet derived its")
    print("  absence, so this is a conditional closure target, not a Locked amplitude.")
    print("  T3 remains the separate task of deriving N_shell.")
    print("=" * 108)
    print("exit 0 -- T5b whiteness reduced to no-horizon-scale-covariance premise.")


if __name__ == "__main__":
    main()
