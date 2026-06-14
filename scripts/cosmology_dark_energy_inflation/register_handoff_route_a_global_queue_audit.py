#!/usr/bin/env python3
r"""REGISTER-HANDOFF ROUTE A — can canon derive the global address queue?

Context:
    Clause B showed that a literal one-edge-per-tick sweep is not tick-equivalent
    to the 9-state queue.  Route A asks whether §5.2 nevertheless licenses the
    *global address queue* used by register_handoff_form_selection.py.

What is already canon-supported:
    * §5.2 says the 12 Q3 strain checks are measured simultaneously, so the
      edge-serial sweep is not the literal §5.2 readout.
    * record_content_from_syndrome.py proves that a single-site syndrome
      increment is an 8-ary vertex address.
    * the finite-bandwidth scheduler supplies one correction/service per tick.

The load-bearing gap:
    The tick-13 queue is stronger than simultaneous snapshot readout.  It acts
    like an event-address FIFO: all generated fault addresses are retained as
    individual records, and one active address is serviced per tick.  A full
    12-bit syndrome snapshot plus one finite service does not generally know
    individual addresses at finite occupancy.

This script tests three exact 256-state readings:
    1. global address FIFO / one-jump queue (control; should match tick 13);
    2. simultaneous full-syndrome snapshot + one decoded-bit toggle;
    3. an optimistic snapshot service that only clears decoded bits overlapping
       the true active set.

Verdict:
    Route A is partially supported but not closed.  §5.2 rules out edge-serial
    measurement and supplies the address alphabet, but it does not by itself
    derive the event-address FIFO needed for the exact transient CC hit.
"""

from __future__ import annotations

import math
from itertools import combinations

from register_handoff_form_selection import (
    BASE_GAMMA,
    rho_ratio_from_q1,
    solve_target,
    transient_readouts,
)


EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def syndrome(state: int) -> tuple[int, ...]:
    return tuple(((state >> a) & 1) ^ ((state >> b) & 1) for a, b in EDGES)


def build_min_weight_decoder() -> dict[tuple[int, ...], list[int]]:
    by_syndrome: dict[tuple[int, ...], list[int]] = {}
    for state in range(256):
        by_syndrome.setdefault(syndrome(state), []).append(state)
    decoder: dict[tuple[int, ...], list[int]] = {}
    for syn, states in by_syndrome.items():
        min_weight = min(state.bit_count() for state in states)
        decoder[syn] = [state for state in states if state.bit_count() == min_weight]
    return decoder


DECODER = build_min_weight_decoder()


def failure_weight(count: int) -> float:
    if count <= 3:
        return 0.0
    if count == 4:
        return 0.5
    return 1.0


def q1(dist: dict[int, float]) -> float:
    return sum(prob * failure_weight(state.bit_count()) for state, prob in dist.items())


def mean_active_fraction(dist: dict[int, float]) -> float:
    return sum(prob * state.bit_count() for state, prob in dist.items()) / 8.0


def generation_step(dist: dict[int, float]) -> dict[int, float]:
    out: dict[int, float] = {}
    for state, prob in dist.items():
        clean = [bit for bit in range(8) if not (state >> bit) & 1]
        for mask in range(1 << len(clean)):
            next_state = state
            weight = prob
            for j, bit in enumerate(clean):
                if (mask >> j) & 1:
                    weight *= BASE_GAMMA
                    next_state |= 1 << bit
                else:
                    weight *= 1.0 - BASE_GAMMA
            out[next_state] = out.get(next_state, 0.0) + weight
    return out


def service_global_fifo(dist: dict[int, float]) -> dict[int, float]:
    out: dict[int, float] = {}
    for state, prob in dist.items():
        active = [bit for bit in range(8) if (state >> bit) & 1]
        if not active:
            out[state] = out.get(state, 0.0) + prob
            continue
        share = prob / len(active)
        for bit in active:
            next_state = state & ~(1 << bit)
            out[next_state] = out.get(next_state, 0.0) + share
    return out


def service_snapshot_toggle(dist: dict[int, float]) -> dict[int, float]:
    """Use the full syndrome snapshot, pick a min-weight decoded bit, and toggle it."""

    out: dict[int, float] = {}
    for state, prob in dist.items():
        decoded_sets = DECODER[syndrome(state)]
        for decoded in decoded_sets:
            decoded_bits = [bit for bit in range(8) if (decoded >> bit) & 1]
            if not decoded_bits:
                out[state] = out.get(state, 0.0) + prob / len(decoded_sets)
                continue
            share = prob / (len(decoded_sets) * len(decoded_bits))
            for bit in decoded_bits:
                next_state = state ^ (1 << bit)
                out[next_state] = out.get(next_state, 0.0) + share
    return out


def service_snapshot_clear_overlap(dist: dict[int, float]) -> dict[int, float]:
    """Optimistic snapshot service: clear decoded bits only if truly active."""

    out: dict[int, float] = {}
    for state, prob in dist.items():
        decoded_sets = DECODER[syndrome(state)]
        for decoded in decoded_sets:
            decoded_bits = [bit for bit in range(8) if (decoded >> bit) & 1]
            if not decoded_bits:
                out[state] = out.get(state, 0.0) + prob / len(decoded_sets)
                continue
            share = prob / (len(decoded_sets) * len(decoded_bits))
            for bit in decoded_bits:
                if (state >> bit) & 1:
                    next_state = state & ~(1 << bit)
                else:
                    next_state = state
                out[next_state] = out.get(next_state, 0.0) + share
    return out


def run(service, ticks: int) -> dict[int, float]:
    dist: dict[int, float] = {0: 1.0}
    for _tick in range(ticks):
        dist = generation_step(dist)
        dist = service(dist)
    return dist


def report(label: str, dist: dict[int, float], q1_target: float) -> tuple[float, float, float]:
    value = q1(dist)
    rho = rho_ratio_from_q1(value, q1_target)
    mean = mean_active_fraction(dist)
    print(f"  {label:<34} q1={value:.12e}  rho/rho_obs={rho:.3e}  <n>/8={mean:.6f}")
    return value, rho, mean


def birth_multiplicity_pressure() -> None:
    print("\n[3] Why snapshot readout is weaker than an address FIFO")
    # At tick 1 from a clean cell, multiple births are already non-negligible
    # under the binomial generation law used by the successful queue.
    probs = {}
    for k in range(9):
        probs[k] = math.comb(8, k) * BASE_GAMMA**k * (1.0 - BASE_GAMMA) ** (8 - k)
    p_multi = sum(probs[k] for k in range(2, 9))
    p_four_plus = sum(probs[k] for k in range(4, 9))
    print(f"  clean-cell P(two or more births in one tick) = {p_multi:.6f}")
    print(f"  clean-cell P(four or more births in one tick) = {p_four_plus:.6e}")
    check(p_multi > 0.15, "multi-birth ticks are not a negligible sub-percent effect")

    # Show the full syndrome increment can decode up to weight 3 but not the
    # whole batch distribution.
    for weight in range(4):
        seen = set()
        for bits in combinations(range(8), weight):
            state = sum(1 << bit for bit in bits)
            syn = syndrome(state)
            seen.add(syn)
        check(len(seen) == math.comb(8, weight), f"weight-{weight} increments are syndrome-distinct")
    weight4_classes = 0
    for bits in combinations(range(8), 4):
        state = sum(1 << bit for bit in bits)
        comp = state ^ 0xFF
        if syndrome(state) == syndrome(comp):
            weight4_classes += 1
    print(f"  weight-4 complement-collision checks = {weight4_classes}/70")
    check(weight4_classes == 70, "all weight-4 increment sets collide with their complements")


def main() -> None:
    print("REGISTER-HANDOFF ROUTE A: GLOBAL ADDRESS-QUEUE AUDIT")
    _p_target, q1_target = solve_target()

    print("\n[1] Canon support that Route A can use")
    print("  §5.2 support: the 12 strain checks are simultaneous QND measurements,")
    print("  not a literal one-edge-per-tick measurement sweep.")
    print("  Record-content support: a single-site syndrome increment is a vertex")
    print("  address; the service variable is 'where', not raw strain count.")
    print("  Finite-bandwidth support: one correction/service action per tick.")
    check(len(EDGES) == 12, "Q3 has the 12 strain edges used by §5.2")

    print("\n[2] Exact tick-13 comparison")
    global_dist = run(service_global_fifo, 13)
    snapshot_dist = run(service_snapshot_toggle, 13)
    optimistic_dist = run(service_snapshot_clear_overlap, 13)

    q_global, rho_global, _ = report("global address FIFO", global_dist, q1_target)
    q_snapshot, rho_snapshot, _ = report("snapshot min-weight toggle", snapshot_dist, q1_target)
    q_optimistic, rho_optimistic, _ = report("snapshot clear-overlap oracle", optimistic_dist, q1_target)

    _pre13, post13 = transient_readouts(BASE_GAMMA, 1, 13)
    check(abs(q_global - post13) < 1e-15, "global FIFO is exactly the 9-state tick-13 queue")
    check(abs(rho_global - 1.0006397433046454) < 1e-10, "global FIFO lands on the transient CC hit")
    check(rho_snapshot > 1e40, "plain simultaneous snapshot decoding does not reproduce the queue")
    check(rho_optimistic > 1e35, "even optimistic snapshot service misses by many orders")

    birth_multiplicity_pressure()

    print("\n[4] Route A verdict")
    print(
        f"""
  Route A gets partway but does not close.

  Positive: canon's simultaneous QND readout and the derived address-increment
  theorem support a global address layer over the failed literal edge sweep.
  That is why the edge-sweep failure should not automatically kill the queue
  route.

  Negative: the exact tick-13 queue is not derived by a simultaneous syndrome
  snapshot.  Snapshot+one-service gives q1={q_snapshot:.6e} (or
  q1={q_optimistic:.6e} under an optimistic overlap-clearing oracle), many orders
  above the target q1={q1_target:.6e}.  The successful model is specifically an
  event-address FIFO: generated fault addresses are retained individually and
  one true active address is serviced per tick.

  The remaining theorem is therefore narrower and cleaner:
      §5.2 must derive an event-address FIFO/demultiplexer, not merely
      simultaneous 12-bit syndrome snapshots.

  If that FIFO is derived, the transient branch is back in play.  If not, the
  exact tick-13 landing remains demoted and the live non-horizon CC route is the
  stationary one-jump result plus the 0.22% gamma-correction branch, with the
  count-slaving candidate still separate class-2.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
