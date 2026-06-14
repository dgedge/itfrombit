#!/usr/bin/env python3
r"""REGISTER-HANDOFF CLAUSE B — explicit edge-sweep microdynamics.

The prior scheduler-closure audit narrowed the exact transient CC route to:

    one covariant 12-edge strain sweep + one stroboscopic recovery/handoff latch.

Clause B is the missing equivalence test:

    Does an explicit one-edge-per-tick Q3 strain scheduler reproduce the 9-state
    one-jump queue's tick-13 post-service residual?

This script answers with exact finite Markov chains over the 8-bit fault subset.

Models compared:
    * global one-jump queue:
        the previous 9-state model, lifted to 256 subsets as a control.
    * fixed round-robin incident-edge oracle:
        each tick services one specified Q3 edge; if either endpoint is faulty,
        one faulty endpoint is cleared.  This is optimistic because it grants
        endpoint localization from the active edge.
    * parity-safe edge service:
        each tick clears only an odd-parity active edge (exactly one faulty
        endpoint); even-parity double faults on the edge are not localized.
    * adaptive best-edge oracle:
        a stronger diagnostic than any fixed covariant sweep in spirit: at each
        tick, choose the edge minimizing the next-step q1.  It is not the physical
        scheduler, but it checks whether the edge-local restriction is already too
        slow even with favourable edge choices.

Verdict:
    The explicit edge-local schedulers do not reproduce the tick-13 one-jump
    residual.  The exact transient CC hit is therefore not established by
    "12-edge sweep + latch" alone; the framework must either derive the global
    address-queue abstraction, derive a new edge-sweep correction factor, or use
    the stationary/gamma-correction branch.
"""

from __future__ import annotations

import math
from itertools import product

from register_handoff_form_selection import (
    BASE_GAMMA,
    rho_ratio_from_q1,
    solve_target,
    transient_readouts,
)


VERTICES = list(product((0, 1), repeat=3))
VERTEX_INDEX = {v: i for i, v in enumerate(VERTICES)}


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def q3_edges() -> list[tuple[int, int, int]]:
    edges: list[tuple[int, int, int]] = []
    for i, vertex in enumerate(VERTICES):
        for axis in range(3):
            other = list(vertex)
            other[axis] ^= 1
            j = VERTEX_INDEX[tuple(other)]
            if i < j:
                edges.append((i, j, axis))
    return edges


EDGES = q3_edges()
ROUND_ROBIN = tuple(range(len(EDGES)))


def failure_weight(count: int) -> float:
    if count <= 3:
        return 0.0
    if count == 4:
        return 0.5
    return 1.0


def subset_q1(dist: dict[int, float]) -> float:
    return sum(prob * failure_weight(state.bit_count()) for state, prob in dist.items())


def mean_active_fraction(dist: dict[int, float]) -> float:
    return sum(prob * state.bit_count() for state, prob in dist.items()) / 8.0


def generation_step(dist: dict[int, float], gamma: float = BASE_GAMMA) -> dict[int, float]:
    out: dict[int, float] = {}
    for state, prob in dist.items():
        clean = [bit for bit in range(8) if not (state >> bit) & 1]
        for mask in range(1 << len(clean)):
            next_state = state
            weight = prob
            for j, bit in enumerate(clean):
                if (mask >> j) & 1:
                    weight *= gamma
                    next_state |= 1 << bit
                else:
                    weight *= 1.0 - gamma
            out[next_state] = out.get(next_state, 0.0) + weight
    return out


def service_global_queue(dist: dict[int, float]) -> dict[int, float]:
    """Clear one active bit anywhere, using uniform tie-breaking on subset states."""

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


def service_incident_edge_oracle(dist: dict[int, float], edge_index: int) -> dict[int, float]:
    """Optimistic edge service: clear one faulty endpoint if the current edge sees one."""

    a, b, _axis = EDGES[edge_index]
    out: dict[int, float] = {}
    for state, prob in dist.items():
        active = [bit for bit in (a, b) if (state >> bit) & 1]
        if not active:
            out[state] = out.get(state, 0.0) + prob
            continue
        share = prob / len(active)
        for bit in active:
            next_state = state & ~(1 << bit)
            out[next_state] = out.get(next_state, 0.0) + share
    return out


def service_parity_safe_edge(dist: dict[int, float], edge_index: int) -> dict[int, float]:
    """Clear only when the current edge has odd parity, leaving double faults intact."""

    a, b, _axis = EDGES[edge_index]
    out: dict[int, float] = {}
    for state, prob in dist.items():
        fa = (state >> a) & 1
        fb = (state >> b) & 1
        if fa ^ fb:
            bit = a if fa else b
            next_state = state & ~(1 << bit)
            out[next_state] = out.get(next_state, 0.0) + prob
        else:
            out[state] = out.get(state, 0.0) + prob
    return out


def run_fixed_policy(policy: str, ticks: int) -> dict[int, float]:
    dist: dict[int, float] = {0: 1.0}
    for tick in range(1, ticks + 1):
        dist = generation_step(dist)
        if policy == "global":
            dist = service_global_queue(dist)
        elif policy == "incident-oracle":
            dist = service_incident_edge_oracle(dist, ROUND_ROBIN[(tick - 1) % 12])
        elif policy == "parity-safe":
            dist = service_parity_safe_edge(dist, ROUND_ROBIN[(tick - 1) % 12])
        else:
            raise ValueError(policy)
    return dist


def run_adaptive_best_edge(ticks: int) -> tuple[dict[int, float], list[int]]:
    dist: dict[int, float] = {0: 1.0}
    chosen: list[int] = []
    for _tick in range(1, ticks + 1):
        generated = generation_step(dist)
        candidates = [
            (subset_q1(service_incident_edge_oracle(generated, edge_index)), edge_index)
            for edge_index in range(len(EDGES))
        ]
        _q, edge_index = min(candidates, key=lambda item: item[0])
        chosen.append(edge_index)
        dist = service_incident_edge_oracle(generated, edge_index)
    return dist, chosen


def report_policy(label: str, dist: dict[int, float], q1_target: float) -> tuple[float, float, float]:
    q1 = subset_q1(dist)
    rho = rho_ratio_from_q1(q1, q1_target)
    mean = mean_active_fraction(dist)
    print(f"  {label:<28} q1={q1:.12e}  rho/rho_obs={rho:.3e}  <n>/8={mean:.6f}")
    return q1, rho, mean


def main() -> None:
    print("REGISTER-HANDOFF CLAUSE B: EXPLICIT EDGE-SWEEP MICRODYNAMICS")
    _p_target, q1_target = solve_target()

    print("\n[1] Control: 256-state global queue reproduces the 9-state tick-13 result")
    global_dist = run_fixed_policy("global", 13)
    q_global, rho_global, _mean_global = report_policy("global one-jump", global_dist, q1_target)
    _pre13, post13 = transient_readouts(BASE_GAMMA, 1, 13)
    q13 = post13
    check(abs(q_global - q13) < 1e-15, "256-state global queue equals the 9-state post-service model")
    check(abs(rho_global - 1.0006397433046454) < 1e-10, "global queue reproduces the exact-grade tick-13 hit")

    print("\n[2] Fixed 12-edge round-robin microdynamics at the same tick")
    incident_dist = run_fixed_policy("incident-oracle", 13)
    parity_dist = run_fixed_policy("parity-safe", 13)
    q_incident, rho_incident, _mean_incident = report_policy(
        "round-robin incident oracle", incident_dist, q1_target
    )
    q_parity, rho_parity, _mean_parity = report_policy("round-robin parity-safe", parity_dist, q1_target)
    check(q_incident > 0.03, "even optimistic fixed edge service leaves far too much residual")
    check(rho_incident > 1e30, "optimistic fixed edge service misses the CC transient by many orders")
    check(q_parity > q_incident, "parity-only localization is slower than the optimistic oracle")

    print("\n[3] Stronger diagnostic: adaptive best-edge oracle")
    adaptive_dist, chosen = run_adaptive_best_edge(13)
    q_adapt, rho_adapt, _mean_adapt = report_policy("adaptive best-edge oracle", adaptive_dist, q1_target)
    print("  chosen edge sequence:", " ".join(str(edge) for edge in chosen))
    check(q_adapt > 0.01, "even adaptive edge-local service remains far above the target q1")
    check(rho_adapt > 1e20, "adaptive edge-local service still misses the transient by many orders")

    print("\n[4] Interpretation")
    print(
        f"""
  The global queue is a true one-fault-per-tick address queue: if any bit is
  active, one active bit is cleared.  That abstraction is exactly what produces
  q1={q_global:.6e} and rho={rho_global:.6f} at tick 13.

  A literal edge sweep is not tick-equivalent.  It only services faults incident
  to the current edge; finite-occupancy localization further slows it.  At tick
  13 the optimistic fixed edge sweep gives q1={q_incident:.6e}, while the
  parity-safe reading gives q1={q_parity:.6e}.  Even the nonphysical adaptive
  best-edge oracle gives q1={q_adapt:.6e}.

  Clause B therefore FAILS for the direct edge-sweep convention.  The exact
  transient cosmological-constant branch now needs a stronger scheduler theorem:
  either derive that §5.2's handoff uses a global address queue rather than a
  literal edge sweep, or compute a new edge-sweep/latch correction and show it
  lands independently.  Otherwise the finite-handoff route falls back to the
  stationary one-jump result (1.616 rho_obs) plus the 0.22% gamma-correction
  branch, while the count-slaving candidate remains separate class-2.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
