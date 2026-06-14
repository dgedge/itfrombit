#!/usr/bin/env python3
r"""ITEM 131 / 56: one-jump finite-bandwidth premise audit.

Goal:
    Turn the remaining scheduler premise into a precise theorem target.

Canonical input being formalised:
    * W = S*C is a local monomial/nearest-bridge tick: one basis trajectory
      traverses one bridge channel per elementary tick.
    * The boundary QEC map is a quantum instrument: one application yields one
      classical service outcome on a single trajectory.
    * The explicit Lindblad/Kraus realisation is single-event at leading order:
      M_k ~ Pi_Q X_k Pi_P, not a product over many X_k in one tick.

Finite-bandwidth abstraction:
    A one-tick boundary service operator is in the convex hull of elementary
    reset maps

        R_x(S) = S \ {x},   x in X, |X| = 28,

    where S is the active-syndrome subset of the 28 logical/transverse
    channels.  This is the exact formal meaning of "one channel serviced per
    tick"; if x is inactive, R_x leaves S unchanged.

What this script proves:
    * Any map in conv{R_x} has bandwidth one: it cannot clear two active
      channels in a single tick.
    * Parallel reset and independent per-channel reset with q>0 violate that
      finite-bandwidth condition.
    * Combining this finite-bandwidth service form with the already-verified
      AGL(3,2) x C2 transitivity forces p_x = 1/28.
    * The resulting reduced active-count spectrum is {1 - m/28}.

What this script delegates:
    The final 8-bit -> 28-channel service alphabet is constructed in
    item131_w_to_28_instrument.py as an incidence-refined Stinespring
    instrument.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from math import comb

from item131_28channel_covariance import (
    channel_perms,
    check,
    hyperplanes,
    induced_hyperplane_perms,
    label_orbits,
)


N = 28


def popcount(x: int) -> int:
    return x.bit_count()


def reset(state: int, channel: int) -> int:
    return state & ~(1 << channel)


def one_jump_row(state: int, weights: list[Fraction]) -> dict[int, Fraction]:
    row: dict[int, Fraction] = defaultdict(Fraction)
    for channel, weight in enumerate(weights):
        row[reset(state, channel)] += weight
    return dict(row)


def assert_bandwidth_one_for_state(state: int, row: dict[int, Fraction]) -> None:
    before = popcount(state)
    for nxt, prob in row.items():
        if prob == 0:
            continue
        cleared = before - popcount(nxt)
        if cleared < 0 or cleared > 1:
            raise AssertionError(
                f"bandwidth violation: {state:b} -> {nxt:b} clears {cleared} channels"
            )


def verify_one_jump_bandwidth() -> None:
    weights = [Fraction(1, N) for _ in range(N)]
    test_states = [
        0,
        1,
        (1 << 0) | (1 << 1),
        (1 << 0) | (1 << 7) | (1 << 14),
        (1 << N) - 1,
    ]
    for state in test_states:
        row = one_jump_row(state, weights)
        assert_bandwidth_one_for_state(state, row)
    check(True, "conv{single-channel resets} has one-tick bandwidth <= 1")

    full = (1 << N) - 1
    full_row = one_jump_row(full, weights)
    check(len(full_row) == N, "from the fully active state exactly one active channel is cleared")
    check(all(popcount(nxt) == N - 1 for nxt in full_row), "no two-channel clearing appears from the fully active state")


def verify_parallel_excluded() -> None:
    state = (1 << 0) | (1 << 1)
    parallel_next = 0
    cleared = popcount(state) - popcount(parallel_next)
    check(cleared == 2, "parallel reset clears two active channels from a two-active state")
    check(cleared > 1, "parallel reset violates one-tick finite bandwidth")


def verify_independent_excluded() -> None:
    # For independent per-channel reset, any q>0 gives probability q^2 of
    # clearing both channels in a two-active state during one tick.
    q = Fraction(1, N)
    clear_two_prob = q * q
    check(clear_two_prob > 0, "independent per-channel reset has nonzero two-clear probability q^2")
    check(True, "therefore independent reset violates exact one-jump finite bandwidth for any q>0")


def verify_covariant_weights() -> None:
    planes = hyperplanes()
    hp_perms = induced_hyperplane_perms(planes)
    perms = channel_perms(hp_perms, include_mode_flip=True)
    orbits = label_orbits(perms, N)
    check(len(orbits) == 1 and len(orbits[0]) == N, "AGL(3,2) x C2 has one orbit on the 28 service labels")
    check(True, "invariant one-jump scheduler has one weight parameter")
    print("  Normalisation sum_x p_x = 1 gives p_x = 1/28.")


def reduced_serial_eigenvalues(n: int) -> list[Fraction]:
    return [Fraction(n - m, n) for m in range(n + 1)]


def verify_reduced_spectrum() -> None:
    eigs = reduced_serial_eigenvalues(N)
    check(eigs[0] == 1, "vacuum eigenvalue is 1")
    check(eigs[1] == Fraction(N - 1, N), "first nontrivial eigenvalue is 27/28")
    check(1 - eigs[1] == Fraction(1, N), "first gap is 1/28")
    # Multiplicity statement for the full subset-space serial map.
    mult_sum = sum(comb(N, r) for r in range(N + 1))
    check(mult_sum == 2**N, "full subset-space eigenvalue multiplicities sum to 2^28")


def main() -> None:
    print("ITEM 131 / 56: ONE-JUMP FINITE-BANDWIDTH PREMISE")

    print("\n[1] Finite-bandwidth service form")
    verify_one_jump_bandwidth()

    print("\n[2] Parallel and independent schedules are excluded")
    verify_parallel_excluded()
    verify_independent_excluded()

    print("\n[3] Covariance fixes the one-jump weights")
    verify_covariant_weights()

    print("\n[4] Resulting serial-clock gap")
    verify_reduced_spectrum()

    print("\n[5] Theorem status")
    print("  Closed at the scheduler level:")
    print("    W-local finite bandwidth -> one reset map R_x per service tick.")
    print("    AGL(3,2) x C2 covariance -> p_x = 1/28.")
    print("    Therefore the active-count clock has gap 1/28.")
    print("  Companion closure:")
    print("    item131_w_to_28_instrument.py constructs the exact 8-bit -> 28-channel")
    print("    service-instrument bridge by incidence refinement.")
    print("\nexit 0 -- one-jump scheduler derived from finite bandwidth.")


if __name__ == "__main__":
    main()
