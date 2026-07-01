#!/usr/bin/env python3
r"""Continuum string-tension scaling audit for the record-grammar confinement rung.

This is not a new simulation.  It is the guardrail around the current
confinement result:

* the finite bond-bipyramid Wilson loops give the leading strong-coupling
  lattice-unit coefficient

      sigma_sc(beta) = -log(beta/18),

  which is a correct strong-coupling statement but not a physical string
  tension;
* a dimensionful ``sqrt(sigma)`` requires a weak-coupling scale-setting bridge:
  measure ``a sqrt(sigma)(beta)`` and show it scales consistently with the
  thermal scale from ``beta_c(N_t)``;
* the existing framework-bulk SU(3) Polyakov data are demonstration-grade
  ``N_t=2,4`` points, not such a bridge.

The script encodes the closure criterion and checks the numbers already in the
paper.  It also records a useful lead: the current bond-bipyramid SU(3)
transition pair is wildly inconsistent with two-loop scaling if its printed
``beta`` is read as the standard Wilson beta, but a simple beta/3 conversion
puts the *two-point ratio* close to the expected factor of two.  With only two
small-volume points, and at non-asymptotic couplings, that is a candidate
normalisation, not a derivation.

The next real computation is therefore precise: finite-size-scale the full
SU(3) bond-bipyramid bulk at several ``N_t`` values, measure large-loop Creutz
ratios at the same beta values, and test whether

      T_c/sqrt(sigma) = 1/(N_t a sqrt(sigma))

stabilises.  Until then, no dimensionful QCD string tension is claimed.
"""

from __future__ import annotations

import math


N = 3
B0 = 11.0 * N / (48.0 * math.pi**2)
B1 = 34.0 * N**2 / (3.0 * (16.0 * math.pi**2) ** 2)


def two_loop_a_lambda(beta: float) -> float:
    """Two-loop asymptotic lattice scaling factor, up to a scheme constant.

    Wilson convention: beta = 2N/g^2, N=3.  Only ratios are used.
    """

    g2 = 2.0 * N / beta
    return math.exp(-1.0 / (2.0 * B0 * g2)) * (B0 * g2) ** (-B1 / (2.0 * B0 * B0))


def two_loop_ratio(beta_nt2: float, beta_nt4: float, scale: float = 1.0) -> float:
    return two_loop_a_lambda(scale * beta_nt2) / two_loop_a_lambda(scale * beta_nt4)


def sigma_strong(beta: float) -> float:
    return -math.log(beta / 18.0)


def tc_over_sqrt_sigma_strong(beta_c: float, nt: int) -> float:
    return 1.0 / (nt * math.sqrt(sigma_strong(beta_c)))


def assert_close_fraction(name: str, value: float, target: float, frac: float) -> None:
    err = abs(value - target) / abs(target)
    print(f"  {name:<72s} value={value:.12g} target={target:.12g} frac_err={err:.3%}")
    if err > frac:
        raise AssertionError(name)


def assert_far_fraction(name: str, value: float, target: float, frac: float) -> None:
    err = abs(value - target) / abs(target)
    print(f"  {name:<72s} value={value:.12g} target={target:.12g} frac_err={err:.3%}")
    if err <= frac:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<72s} value={value}")
    if not value:
        raise AssertionError(name)


def solve_scale_for_ratio(beta_nt2: float, beta_nt4: float, target: float = 2.0) -> float:
    lo, hi = 0.01, 2.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if two_loop_ratio(beta_nt2, beta_nt4, mid) > target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main() -> None:
    print("Continuum string-tension scaling audit")
    print("=" * 96)
    print(f"  pure SU(3) two-loop coefficients: b0={B0:.12g}, b1={B1:.12g}")

    # Known hypercubic SU(3) deconfinement couplings used as a sanity check for
    # the two-loop ratio convention.  These are not framework inputs.
    hyper_nt2 = 5.09
    hyper_nt4 = 5.6925

    # Demonstration-grade full-SU(3) bond-bipyramid bulk values recorded by
    # record_grammar_tch_su3_bulk_polyakov_torus.py / the confinement paper.
    bulk_nt2 = 4.3
    bulk_nt4 = 6.5

    print("\n[1] Thermal scale-ratio sanity check")
    hyper_ratio = two_loop_ratio(hyper_nt2, hyper_nt4)
    bulk_raw_ratio = two_loop_ratio(bulk_nt2, bulk_nt4)
    bulk_third_ratio = two_loop_ratio(bulk_nt2, bulk_nt4, scale=1.0 / 3.0)
    best_scale = solve_scale_for_ratio(bulk_nt2, bulk_nt4)
    best_ratio = two_loop_ratio(bulk_nt2, bulk_nt4, scale=best_scale)

    print(f"  hypercubic SU(3): beta_c(2)={hyper_nt2:.4f}, beta_c(4)={hyper_nt4:.4f}")
    print(f"  bond-bipyramid SU(3 demo): beta_c(2)~{bulk_nt2:.3f}, beta_c(4)~{bulk_nt4:.3f}")
    assert_close_fraction("hypercubic pair gives expected a2/a4 ~= 2", hyper_ratio, 2.0, 0.05)
    assert_far_fraction("raw bond-bipyramid printed beta does not scale as standard beta", bulk_raw_ratio, 2.0, 0.50)
    assert_close_fraction("candidate beta_eff=beta/3 gives two-point factor ~= 2", bulk_third_ratio, 2.0, 0.05)
    assert_close_fraction("best two-point scale is close to 1/3", best_scale, 1.0 / 3.0, 0.03)
    assert_close_fraction("best-scale ratio lands by construction", best_ratio, 2.0, 1e-12)

    print("\n[2] Strong-coupling sigma is not a scale-setting sigma")
    t2 = tc_over_sqrt_sigma_strong(bulk_nt2, 2)
    t4 = tc_over_sqrt_sigma_strong(bulk_nt4, 4)
    spread = abs(t2 - t4) / (0.5 * (t2 + t4))
    print(f"  using sigma_sc=-log(beta/18): T_c/sqrt(sigma_sc) at Nt=2 = {t2:.6f}")
    print(f"  using sigma_sc=-log(beta/18): T_c/sqrt(sigma_sc) at Nt=4 = {t4:.6f}")
    print(f"  fractional two-point spread = {spread:.1%}")
    assert_true("strong-coupling sigma fails the scale-setting constancy test", spread > 0.20)

    print(
        f"""
VERDICT:
  FRONTIER SHARPENED, NOT CLOSED.

  The existing Wilson-loop coefficient sigma=-log(beta/18) is a valid leading
  strong-coupling lattice-unit result.  It is not the physical QCD string
  tension: when combined with the current beta_c(N_t) demonstration values it
  fails the elementary T_c/sqrt(sigma) constancy test.

  The current full-SU(3) bond-bipyramid transition pair also cannot be read as a
  standard Wilson-beta continuum scaling result: the raw two-loop scale ratio is
  {bulk_raw_ratio:.3f}, not 2.  There is one useful lead, however: beta_eff=beta/3
  gives ratio {bulk_third_ratio:.3f}, and the exact two-point fitted scale is
  {best_scale:.6f}.  This is a candidate action-normalisation bridge, not a
  closure, because two small-volume points and non-asymptotic beta values cannot
  establish weak-coupling scaling.

  Closure requires a real SU(3) bulk scale-setting run:
    (i) beta_c(N_t) for several N_t with finite-size scaling;
   (ii) a sqrt(sigma)(beta) from large Wilson/Creutz loops at the same betas;
  (iii) a stable T_c/sqrt(sigma) and a derived beta-normalisation convention.

  Until that exists, the confinement paper's string tension remains leading
  strong-coupling in lattice units, and the dimensionful physical sigma is a
  live hard frontier.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
