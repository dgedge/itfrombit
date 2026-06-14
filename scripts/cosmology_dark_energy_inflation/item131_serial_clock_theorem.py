#!/usr/bin/env python3
r"""ITEM 131 serial-clock theorem candidate.

Goal:
    Prove the operator-theory part of the 28-channel absorbing QEC clock:

        if one substrate tick services exactly one of N=28 orthogonal
        boundary channels, uniformly, and clears it when active,

    then the resulting non-depolarizing absorbing clock has first gap 1/28.

What this proves:
    The 1/28 gap does not require depolarization and is not an arbitrary
    lazy rate once "one-channel-per-tick" serial service is supplied.

What it does NOT prove by itself:
    That the canonical substrate walk W=S*C supplies the 28 logical/transverse
    service labels.  Companion scripts derive uniformity by covariance and
    formalise the one-jump finite-bandwidth scheduler; the explicit
    W=S*C -> 2x14 boundary service instrument is constructed separately in
    item131_w_to_28_instrument.py.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np


N = 28


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def reduced_serial_matrix(n: int) -> np.ndarray:
    """Transition matrix on m=number of active syndromes.

    Each tick chooses one channel uniformly.  If the channel is active,
    m decreases by one; otherwise m is unchanged.
    """
    transition = np.zeros((n + 1, n + 1))
    for m in range(n + 1):
        transition[m, m] = 1.0 - m / n
        if m:
            transition[m, m - 1] = m / n
    return transition


def full_serial_matrix(n: int) -> np.ndarray:
    """Full 2^n subset transition matrix for small n checks."""
    states = list(range(1 << n))
    index = {state: i for i, state in enumerate(states)}
    transition = np.zeros((len(states), len(states)))
    for state in states:
        row = index[state]
        for channel in range(n):
            cleared = state & ~(1 << channel)
            transition[row, index[cleared]] += 1.0 / n
    return transition


def full_parallel_matrix(n: int) -> np.ndarray:
    """Parallel service: every active channel is cleared in one tick."""
    states = list(range(1 << n))
    transition = np.zeros((len(states), len(states)))
    transition[:, 0] = 1.0
    return transition


def full_independent_matrix(n: int, q: float) -> np.ndarray:
    """Each active channel is independently cleared with probability q."""
    states = list(range(1 << n))
    index = {state: i for i, state in enumerate(states)}
    transition = np.zeros((len(states), len(states)))
    for state in states:
        active = [i for i in range(n) if state & (1 << i)]
        inactive_mask = state & ~sum(1 << i for i in active)
        row = index[state]
        for keep_bits in itertools.product([0, 1], repeat=len(active)):
            prob = 1.0
            next_state = inactive_mask
            for bit, channel in zip(keep_bits, active):
                if bit:
                    prob *= 1.0 - q
                    next_state |= 1 << channel
                else:
                    prob *= q
            transition[row, index[next_state]] += prob
    return transition


def spectral_gap_row_stochastic(transition: np.ndarray) -> float:
    eigs = np.linalg.eigvals(transition)
    nontrivial = sorted([abs(v) for v in eigs if abs(v - 1.0) > 1e-10], reverse=True)
    if not nontrivial:
        return 1.0
    return 1.0 - nontrivial[0]


def multiplicities(vals: np.ndarray, decimals: int = 8) -> dict[float, int]:
    return dict(sorted(Counter(float(v) for v in np.round(vals, decimals)).items()))


def main() -> None:
    print("ITEM 131 SERIAL ABSORBING CLOCK THEOREM")

    print("\n[1] Reduced N=28 theorem")
    reduced = reduced_serial_matrix(N)
    eigs = np.linalg.eigvals(reduced)
    expected = np.array([1.0 - m / N for m in range(N + 1)])
    check(np.allclose(reduced.sum(axis=1), 1.0), "serial reduced clock is stochastic")
    vacuum = np.zeros(N + 1)
    vacuum[0] = 1.0
    check(np.allclose(vacuum @ reduced, vacuum), "vacuum/no-syndrome state is absorbing")
    check(np.allclose(np.sort(np.real(eigs)), np.sort(expected)), "reduced spectrum is exactly {1-m/28}")
    check(abs(spectral_gap_row_stochastic(reduced) - 1.0 / N) < 1e-12, "first gap is 1/28")

    print("\n[2] Full subset-space check for small n")
    for n in [3, 4, 5, 6]:
        full = full_serial_matrix(n)
        gap = spectral_gap_row_stochastic(full)
        eigs = np.linalg.eigvals(full)
        expected_mult = {round(1.0 - r / n, 8): int(math.comb(n, r)) for r in range(n + 1)}
        check(abs(gap - 1.0 / n) < 1e-12, f"full subset serial clock n={n} has gap 1/{n}")
        check(multiplicities(np.real(eigs)) == expected_mult, f"full subset spectrum n={n} has binomial multiplicities")

    print("\n[3] Schedule alternatives")
    serial_gap = spectral_gap_row_stochastic(full_serial_matrix(6))
    parallel_gap = spectral_gap_row_stochastic(full_parallel_matrix(6))
    independent_gap = spectral_gap_row_stochastic(full_independent_matrix(6, 1.0 / 6.0))
    check(abs(serial_gap - 1.0 / 6.0) < 1e-12, "exclusive serial service derives 1/n from channel count")
    check(abs(parallel_gap - 1.0) < 1e-12, "parallel service gives gap 1, not 1/n")
    check(abs(independent_gap - 1.0 / 6.0) < 1e-12, "independent per-channel service can match 1/n only by setting q=1/n")

    print("\n[4] Proof status")
    print("  Proven finite theorem:")
    print("    uniform exclusive service over N orthogonal channels => absorbing non-depolarizing gap 1/N.")
    print("  For N=28:")
    print("    gap = 1/28 and vacuum is absorbing.")
    print("  Companion bridge:")
    print("    item131_w_to_28_instrument.py realises the 2x14 logical/transverse")
    print("    service alphabet as an incidence refinement of the 8 single-bit Kraus events.")
    print("  Existing support:")
    print("    discrete tick, single-bit Kraus jumps, bipartite P^2 clock, sequential k-index")
    print("    microsteps.")
    print("  Companion cosmology lifts:")
    print("    item131_primordial_tilt_logscale.py gives the exact log-scale")
    print("    generator theorem form for n_s; item131_late_activation.py and")
    print("    item131_r4_support_dimension.py handle the late w(a) branch.")
    print("\nexit 0 -- serial-clock theorem proved; W-to-28 instrument lift handled separately.")


if __name__ == "__main__":
    main()
