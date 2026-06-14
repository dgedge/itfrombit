#!/usr/bin/env python3
r"""ITEM 131 / ITEM 56 constructive gap-closure harness.

Question:
    Does the shared 1/28 spectral gap genuinely derive both

        n_s = 27/28
        w0  = -27/28

    without using the failed depolarizing/FDT route?

Answer encoded by the checks:

* The 28 count is real: the [8,4,4] extended Hamming code has 14
  weight-4 words, and 2 transverse spatial modes gives 28 channels.

* The natural graph operator on the 14 weight-4 words does not have gap 1/28.
  Joining two weight-4 words when their supports overlap in 2 positions
  gives the cocktail-party graph K_14 minus 7 complementary-pair edges.
  Its normalized random-walk/Laplacian relaxation gap is 1.  A gap 1/28
  can be produced by a lazy/rate-rescaled channel, but then 1/28 is the
  inserted clock rate, not an emergent code spectral gap.  This is the
  real leg-1 obstruction.

* There is, however, a natural non-depolarizing substrate-clock candidate:
  a serial absorbing QEC erasure clock on the 28 channels.  Each substrate
  tick services exactly one of the 28 active channels and clears it if it
  carries a syndrome.  The vacuum is absorbing/low-entropy, not maximally
  mixed.  The active-mode spectrum is 1 - m/28, so the first decay gap is
  exactly 1/28.  Companion audits derive uniformity by covariance, formalise
  the one-jump finite-bandwidth scheduler, and construct the exact W=S*C
  -> 2x14 service instrument by incidence refinement.

* The ordinary "3D Fourier transform of an exponential correlation"
  route does NOT derive a constant tilt 1 - n_s = Delta.  It gives a
  Lorentzian-squared spectrum whose dimensionless log-slope is

      3 - 4 (k xi)^2 / (1 + (k xi)^2),

  which is scale-dependent and equals -Delta only at one chosen pivot.

* The early n_s leg can be improved if the 28-channel operator is a
  transfer operator in logarithmic scale, not physical radius:

      Delta_R^2(k) = A_s (k/k_*)^{-Delta}.

  Then n_s - 1 = -Delta exactly.  This is an anomalous-dimension/RG
  statement.  It avoids the naive Fourier failure.  The companion audit
  item131_primordial_tilt_logscale.py sharpens the exact theorem form:
  use the continuous Markov generator Q=P-I, saturated horizon printing
  d ln k=d ln a, and density/power-level action.  A literal discrete
  multiplier, Hubble drift, or amplitude-level action shifts the coefficient.
  The companion item131_early_logscale_lift.py further shows that log-scale
  clocking follows from scale-covariant Markov horizon-shell printing:
  T_{lambda mu}=T_lambda T_mu forces tau=q ln lambda.

* The late w(a) leg is not implied by "a gap Delta" alone.  The continuity
  equation says

      d ln rho_DE / d ln a = -3(1+w).

  The item-131 law w(a)=-1+a Delta is equivalent to

      rho_DE(a) = rho0 exp[3 Delta (1-a)].

  A constant gap per ln a would instead give w=-1+Delta/3.  To get
  w(a)=-1+a Delta one needs a specific late response law:

      d ln rho_DE / d ln a = -3 a Delta,

  i.e. the ordinary FRW volume-continuity factor multiplied by a linearly
  activated R4/boundary fraction a.  The companion audit
  item131_late_activation.py sharpens this: a comoving d-dimensional active
  support gives f_d(a)=a^d and local CPL slope w_a=-d Delta, so item 131's
  w_a=-Delta uniquely selects d=1.  Companion audits now prove the finite R4
  support-dimension lemma, the in-instrument no-extra homogeneous ledger, and
  no independent generation-covariant R5 inside the current 8-bit physical
  register.  The residual is outside-sector/non-R4-coupling completeness, not
  re-deriving the FRW factor 3.

Therefore the best honest upgrade is:
    "28-channel sector verified; natural code graph gap is 1; 1/28 can be
     made natural by a serial absorbing QEC clock with an incidence-refined
     W-to-28 service instrument; the early log-scale map closes under the
     scale-covariant Markov HBC / saturated horizon / power-ledger reading;
     the late R4 support map closes at finite support level."
"""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np


DELTA = 1.0 / 28.0
N_CHANNELS = 28


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def codewords() -> list[int]:
    # A standard generator for the [8,4,4] extended Hamming code.
    generator = [0b10000111, 0b01001011, 0b00101101, 0b00011110]
    words = []
    for bits in itertools.product([0, 1], repeat=4):
        cw = 0
        for bit, row in zip(bits, generator):
            if bit:
                cw ^= row
        words.append(cw)
    return words


def verify_28_count() -> list[int]:
    words = codewords()
    weights = Counter(int.bit_count(word) for word in words)
    print(f"  weight distribution: {dict(sorted(weights.items()))}")
    check(weights[4] == 14, "extended Hamming [8,4,4] has 14 weight-4 words")
    check(2 * weights[4] == 28, "2 transverse modes x 14 logical words = 28 channels")
    return [word for word in words if int.bit_count(word) == 4]


def supports(words: list[int]) -> list[set[int]]:
    return [{i for i in range(8) if word & (1 << i)} for word in words]


def multiplicities(vals: np.ndarray, decimals: int = 8) -> dict[float, int]:
    rounded = np.round(vals, decimals)
    return dict(sorted(Counter(float(v) for v in rounded).items()))


def verify_natural_code_operator(weight4_words: list[int]) -> None:
    supp = supports(weight4_words)
    n = len(supp)
    adjacency = np.zeros((n, n))
    complement_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            overlap = len(supp[i] & supp[j])
            if overlap == 2:
                adjacency[i, j] = adjacency[j, i] = 1.0
            elif overlap == 0:
                complement_pairs += 1
            else:
                raise AssertionError(f"unexpected overlap {overlap} for pair {(i, j)}")

    degree = adjacency.sum(axis=1)
    check(n == 14, "natural weight-4 sector has 14 vertices")
    check(complement_pairs == 7, "the 14 words split into 7 complementary pairs")
    check(np.all(degree == 12), "overlap-2 graph is 12-regular: K14 minus a perfect matching")

    adjacency_eigs = np.linalg.eigvalsh(adjacency)
    check(multiplicities(adjacency_eigs) == {-2.0: 6, 0.0: 7, 12.0: 1}, "adjacency spectrum is {12, 0^7, (-2)^6}")

    transition = adjacency / 12.0
    laplacian = np.eye(n) - transition
    lap_eigs = np.linalg.eigvalsh(laplacian)
    check(multiplicities(lap_eigs) == {0.0: 1, 1.0: 7, 1.16666667: 6}, "normalized random-walk Laplacian spectrum is {0, 1^7, (7/6)^6}")
    check(abs(lap_eigs[1] - 1.0) < 1e-12, "natural code relaxation gap is 1, not 1/28")

    epsilon = DELTA
    lazy_transition = (1.0 - epsilon) * np.eye(n) + epsilon * transition
    lazy_gap = np.sort(np.linalg.eigvalsh(np.eye(n) - lazy_transition))[1]
    check(abs(lazy_gap - DELTA) < 1e-12, "a lazy/rate-rescaled natural walk can have gap 1/28")
    check(True, "but the 1/28 is exactly the inserted laziness/rate, not derived by the graph")

    # A non-depolarizing reset chain shows that "gap 1/28 => depolarizing" is too strong.
    # Row-stochastic convention: each state stays put with 1-epsilon and resets to state 0 with epsilon.
    reset = (1.0 - epsilon) * np.eye(28)
    reset[:, 0] += epsilon
    reset_eigs = np.linalg.eigvals(reset)
    nontrivial = sorted([abs(v) for v in reset_eigs if abs(v - 1.0) > 1e-10], reverse=True)
    reset_gap = 1.0 - nontrivial[0]
    stationary = np.zeros(28)
    stationary[0] = 1.0
    check(abs(reset_gap - DELTA) < 1e-12, "non-depolarizing pure-vacuum reset channel can also have gap 1/28")
    check(np.allclose(stationary @ reset, stationary), "reset channel fixed point is pure/vacuum-like, not maximally mixed")
    check(True, "therefore the obstruction is rate insertion, not a theorem that 1/28 forces depolarization")


def verify_serial_absorbing_qec_clock() -> None:
    """Absorbing erasure clock reduced by the number m of active syndromes.

    State m means m of the 28 channels still carry an active syndrome.  Each
    tick chooses one of the 28 channels uniformly.  If it is active, m->m-1;
    otherwise m stays fixed.  The transition matrix is triangular with
    eigenvalues 1-m/28.  The vacuum m=0 is absorbing.
    """
    n = N_CHANNELS
    transition = np.zeros((n + 1, n + 1))
    for m in range(n + 1):
        transition[m, m] = 1.0 - m / n
        if m > 0:
            transition[m, m - 1] = m / n
    check(np.allclose(transition.sum(axis=1), 1.0), "serial erasure clock is stochastic")
    vacuum = np.zeros(n + 1)
    vacuum[0] = 1.0
    check(np.allclose(vacuum @ transition, vacuum), "vacuum/no-syndrome state is absorbing")
    eigs = np.linalg.eigvals(transition)
    expected = np.array([1.0 - m / n for m in range(n + 1)])
    check(np.allclose(np.sort(np.real(eigs)), np.sort(expected)), "serial clock spectrum is {1-m/28}")
    gap = 1.0 - sorted([abs(v) for v in eigs if abs(v - 1.0) > 1e-10], reverse=True)[0]
    check(abs(gap - DELTA) < 1e-12, "serial absorbing QEC clock has natural first gap 1/28")

    # Expected active syndrome count decays by the same first gap.
    active_expectation = np.arange(n + 1, dtype=float)
    next_expectation = transition @ active_expectation
    ratios = next_expectation[1:] / active_expectation[1:]
    check(np.allclose(ratios, 1.0 - DELTA), "expected active-syndrome count decays by factor 1-1/28 per tick")


def exponential_correlation_dimensionless_slope(kxi: np.ndarray) -> np.ndarray:
    # C(r)=exp(-r/xi) in 3D gives P(k) proportional to (1+k^2 xi^2)^-2.
    # Dimensionless spectrum k^3 P(k) has the log-slope below.
    q = kxi * kxi
    return 3.0 - 4.0 * q / (1.0 + q)


def log_transfer_slope(k_over_kstar: np.ndarray, delta: float) -> np.ndarray:
    power = k_over_kstar ** (-delta)
    return np.gradient(np.log(power), np.log(k_over_kstar))


def w_from_rho(a: np.ndarray, rho: np.ndarray) -> np.ndarray:
    dlnrho_dln_a = np.gradient(np.log(rho), np.log(a))
    return -1.0 - dlnrho_dln_a / 3.0


def rho_for_activation_exponent(a: np.ndarray, delta: float, exponent: float) -> np.ndarray:
    """Density law implied by 1+w(a)=delta*a^exponent, normalized rho(1)=1."""
    if abs(exponent) < 1e-12:
        return a ** (-3.0 * delta)
    return np.exp(-(3.0 * delta / exponent) * (a**exponent - 1.0))


def main() -> None:
    print("ITEM 131 / 56 SHARED-GAP CLOSURE HARNESS")
    print(f"Delta = 1/28 = {DELTA:.8f}")

    print("\n[1] The 28-channel count is real")
    weight4_words = verify_28_count()

    print("\n[2] Shared-operator leg: natural code gap is not 1/28")
    verify_natural_code_operator(weight4_words)

    print("\n[3] Shared-operator leg: serial absorbing QEC clock can make 1/28 natural")
    verify_serial_absorbing_qec_clock()

    print("\n[4] Early n_s leg: ordinary 3D Fourier exponential route fails")
    kxi = np.logspace(-2, 2, 2000)
    slope = exponential_correlation_dimensionless_slope(kxi)
    target = -DELTA
    pivot_index = int(np.argmin(np.abs(slope - target)))
    pivot = kxi[pivot_index]
    check(slope.max() > 2.9 and slope.min() < -0.99, "exponential-correlation slope runs from blue to red; it is not constant")
    check(abs(slope[pivot_index] - target) < 5e-3, f"slope equals -Delta only at a tuned pivot k xi ~= {pivot:.3f}")
    check(np.std(slope[(kxi > 0.3) & (kxi < 3.0)]) > 0.8, "no broad scale window has constant n_s-1=-Delta")

    print("\n[5] Early n_s leg: log-scale transfer operator would close the tilt")
    kk = np.logspace(-3, 3, 5000)
    log_slope = log_transfer_slope(kk, DELTA)
    check(np.max(np.abs(log_slope[10:-10] + DELTA)) < 1e-6, "Delta_R^2(k)=(k/k*)^-Delta gives n_s-1=-Delta exactly")
    ns = 1.0 - DELTA
    check(abs(ns - 27.0 / 28.0) < 1e-15, "n_s = 1 - Delta = 27/28")

    print("\n[6] Late w(a) leg: continuity equation requirements")
    a = np.logspace(-3, 0, 5000)
    rho_item131 = np.exp(3.0 * DELTA * (1.0 - a))
    w_item131 = w_from_rho(a, rho_item131)
    idx_today = -20
    check(abs(w_item131[idx_today] - (-1.0 + a[idx_today] * DELTA)) < 1e-6, "rho=rho0 exp[3 Delta (1-a)] gives w(a)=-1+a Delta")
    w0 = -1.0 + DELTA
    wa = -DELTA
    check(abs(w0 + 27.0 / 28.0) < 1e-15 and abs(wa + 1.0 / 28.0) < 1e-15, "CPL values are w0=-27/28, wa=-1/28")

    rho_constant_gap = a ** (-DELTA)
    w_constant_gap = w_from_rho(a, rho_constant_gap)
    check(abs(np.median(w_constant_gap) - (-1.0 + DELTA / 3.0)) < 1e-8, "a bare constant gap per ln a gives w=-1+Delta/3, not item 131")
    check(abs((-1.0 + DELTA / 3.0) - w0) > 0.02, "the FRW continuity factor is load-bearing; a bare gap alone is insufficient")

    rho_three_channel_constant = a ** (-3.0 * DELTA)
    w_three_channel = w_from_rho(a, rho_three_channel_constant)
    check(abs(np.median(w_three_channel) - w0) < 1e-8, "three spatial channels can give w0=-1+Delta")
    check(np.std(w_three_channel) < 1e-10, "but without linear activation it gives wa=0, not w(a)=-1+aDelta")

    print("\n[7] Late w(a) leg: activation exponent selects the parameterisation")
    for exponent, label in [(-1.0, "area/volume"), (0.0, "constant"), (1.0, "1D R4 string length"), (2.0, "area"), (3.0, "volume")]:
        rho = rho_for_activation_exponent(a, DELTA, exponent)
        w = w_from_rho(a, rho)
        today = w[-20]
        wa_eff = (w[-20] - w[-200]) / (a[-20] - a[-200])
        print(f"  exponent {exponent:+.0f} ({label:18s}): w(near today)={today:.6f}, local dw/da={wa_eff:.6f}")
    rho_linear = rho_for_activation_exponent(a, DELTA, 1.0)
    w_linear = w_from_rho(a, rho_linear)
    check(abs(w_linear[-20] - (-1.0 + a[-20] * DELTA)) < 1e-6, "1D activation f(a)=a gives item-131 w(a)")
    rho_area = rho_for_activation_exponent(a, DELTA, 2.0)
    w_area = w_from_rho(a, rho_area)
    check(abs(w_area[-20] - (-1.0 + (a[-20] ** 2) * DELTA)) < 1e-6, "2D area activation would give w(a)=-1+a^2 Delta, not item 131")

    print("\n[8] Duality check")
    check(abs(w0 + ns) < 1e-15, "w0=-n_s follows once both formula forms are assumed")
    for d in [1 / 28, 1 / 29, 1 / 56, 0.05]:
        check(abs((-1.0 + d) + (1.0 - d)) < 1e-15, f"w0=-n_s is automatic for Delta={d:.5f}")

    print("\n" + "=" * 96)
    print("CONDITIONAL IMPROVEMENT:")
    print("  * Do not use depolarizing/FDT.  It is the wrong physical limit.")
    print("  * Leg 1 is sharper after pushback: the natural weight-4 code graph has")
    print("    gap 1, not 1/28.  Exact 1/28 currently enters as a clock/rate factor.")
    print("  * Progress: a serial absorbing QEC erasure clock on the 28 channels is")
    print("    non-depolarizing, vacuum-preserving, and has first gap 1/28.  Companion")
    print("    audits derive covariance-uniform weights, the finite-bandwidth scheduler,")
    print("    and the W-to-28 boundary instrument by incidence refinement.")
    print("  * Early leg: companion audit item131_primordial_tilt_logscale.py shows")
    print("    the exact theorem form.  The serial clock must be read as the")
    print("    continuous generator Q=P-I acting on the scalar power ledger in")
    print("    saturated horizon scale d ln k=d ln a.  Then n_s=27/28 exactly.")
    print("    Companion item131_early_logscale_lift.py upgrades log scale to a")
    print("    scale-covariant HBC consequence, with radial q=1 normalization.")
    print("  * Late leg: companion audits reduce the response to support dimension")
    print("    and prove the finite R4 boundary-QEC support is 1D.  With the")
    print("    homogeneous comoving line-ledger lift, f(a)=a and w(a)=-1+a/28.")
    print("  * Remaining theorem burden is now the physical cosmology premise:")
    print("    saturated self-similar radial HBC for the early branch, plus")
    print("    outside-sector/non-R4-coupling completeness. Companion audit")
    print("    item131_no_r5_instrument_completeness.py closes independent")
    print("    generation-covariant R5 inside the current 8-bit instrument.")
    print("=" * 96)


if __name__ == "__main__":
    main()
