#!/usr/bin/env python3
r"""REGISTER-HANDOFF DEMUX DECISION AUDIT.

Question:
    Which route is now honest?

      A. promote the event-address FIFO/demultiplexer from §5.2/register
         mechanics and keep the tick-13 transient cosmological-constant hit live;
      B. keep tick 13 demoted and use the stationary one-jump result
         (1.616 rho_obs) plus the 0.22% gamma-correction theorem target.

Decision logic:
    The successful finite-handoff chain is not merely "one reset per tick".
    It is active-address aware: when faults are present, one true active address
    is serviced.  That has to be distinguished from two weaker readings that
    canon might otherwise suggest:

      * blind covariant service: choose one address/channel from the alphabet;
        if it is inactive, nothing is repaired;
      * syndrome-snapshot service: decode the current 12-bit strain word, with
        no event-address memory.

    This script compares those readings against the exact depth-6 q1 target.
    It also separates two logically distinct conclusions:

      * active-address demux is selected as the viable QEC recovery reading;
      * the clean-start tick-13 dwell is still not derived, because the
        edge-sweep/latch candidate has been killed and no replacement dwell
        theorem exists.

Verdict:
    Use the stationary active-demux one-jump branch as the live non-horizon
    cosmological-constant route.  Tick 13 remains demoted unless an independent
    dwell theorem is later derived.
"""

from __future__ import annotations

import math

from register_handoff_form_selection import (
    BASE_GAMMA,
    q1_from_dist,
    rho_ratio_from_q1,
    row_step,
    solve_gamma_for_stationary_queue,
    solve_target,
    stationary,
    transient_readouts,
)
from register_handoff_route_a_global_queue_audit import (
    q1 as q1_256,
    run as run_256,
    service_global_fifo,
    service_snapshot_clear_overlap,
    service_snapshot_toggle,
)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def transition_matrix(gamma: float, policy: str) -> list[list[float]]:
    """Transition on post-service active-fault count k=0..8.

    Generation is the same binomial-on-clean-bits law used by the finite
    handoff chain.  The service policy differs:

      active_demux:
        if pre-service count h>0, one true active address is cleared.
      blind_point:
        choose one of the 8 point addresses uniformly; it repairs only if active.
      blind_hyperplane:
        optimistic 28-channel coarse service.  A channel has four point
        preimages; if at least one preimage is active, count one repair.  This is
        intentionally generous to the blind-channel branch.
    """

    out = [[0.0 for _ in range(9)] for _ in range(9)]
    for k in range(9):
        clean = 8 - k
        for births in range(clean + 1):
            prob = (
                math.comb(clean, births)
                * gamma**births
                * (1.0 - gamma) ** (clean - births)
            )
            pre = k + births
            if policy == "active_demux":
                out[k][max(0, pre - 1)] += prob
            elif policy == "blind_point":
                if pre == 0:
                    out[k][0] += prob
                else:
                    p_hit = pre / 8.0
                    out[k][pre - 1] += prob * p_hit
                    out[k][pre] += prob * (1.0 - p_hit)
            elif policy == "blind_hyperplane":
                if pre == 0:
                    out[k][0] += prob
                else:
                    miss = math.comb(8 - pre, 4) / math.comb(8, 4) if 8 - pre >= 4 else 0.0
                    p_hit = 1.0 - miss
                    out[k][pre - 1] += prob * p_hit
                    out[k][pre] += prob * (1.0 - p_hit)
            else:
                raise ValueError(policy)
    return out


def read_stationary(policy: str, q1_target: float) -> tuple[float, float, float]:
    matrix = transition_matrix(BASE_GAMMA, policy)
    dist = stationary(matrix)
    q = q1_from_dist(dist)
    mean = sum(k * p for k, p in enumerate(dist)) / 8.0
    return q, rho_ratio_from_q1(q, q1_target), mean


def read_tick(policy: str, tick: int, q1_target: float) -> tuple[float, float, float]:
    matrix = transition_matrix(BASE_GAMMA, policy)
    dist = [0.0 for _ in range(9)]
    dist[0] = 1.0
    for _ in range(tick):
        dist = row_step(dist, matrix)
    q = q1_from_dist(dist)
    mean = sum(k * p for k, p in enumerate(dist)) / 8.0
    return q, rho_ratio_from_q1(q, q1_target), mean


def report_count_policy(label: str, policy: str, q1_target: float) -> None:
    q_stat, rho_stat, mean_stat = read_stationary(policy, q1_target)
    q_13, rho_13, mean_13 = read_tick(policy, 13, q1_target)
    print(
        f"  {label:<29} stationary rho={rho_stat:.3e} q1={q_stat:.6e} <n>/8={mean_stat:.4f}"
    )
    print(
        f"  {'':<29} tick-13    rho={rho_13:.3e} q1={q_13:.6e} <n>/8={mean_13:.4f}"
    )


def report_snapshot_policy(label: str, service, q1_target: float) -> None:
    dist = run_256(service, 13)
    q = q1_256(dist)
    rho = rho_ratio_from_q1(q, q1_target)
    mean = sum(prob * state.bit_count() for state, prob in dist.items()) / 8.0
    print(f"  {label:<29} tick-13    rho={rho:.3e} q1={q:.6e} <n>/8={mean:.4f}")


def main() -> None:
    print("REGISTER-HANDOFF DEMUX DECISION AUDIT")
    _p_target, q1_target = solve_target()

    print("\n[1] Premise separation")
    print("  closed: §5.2 gives simultaneous QND strain readout.")
    print("  closed: record-content audit gives a single-increment vertex address.")
    print("  closed: finite bandwidth allows at most one repair/reset per tick.")
    print("  not closed by those alone: active-address demultiplexing at finite occupancy.")
    print("  This audit asks which operational reading survives the q1 target.")

    print("\n[2] Scheduler readings against the exact q1 target")
    report_count_policy("active-address demux", "active_demux", q1_target)
    report_count_policy("blind 8-address reset", "blind_point", q1_target)
    report_count_policy("blind 28-channel optimistic", "blind_hyperplane", q1_target)
    report_snapshot_policy("snapshot min-weight", service_snapshot_toggle, q1_target)
    report_snapshot_policy("snapshot overlap oracle", service_snapshot_clear_overlap, q1_target)
    report_snapshot_policy("256-state active FIFO", service_global_fifo, q1_target)

    active_stat = read_stationary("active_demux", q1_target)[1]
    active_tick = read_tick("active_demux", 13, q1_target)[1]
    blind8 = read_stationary("blind_point", q1_target)[1]
    blind28 = read_stationary("blind_hyperplane", q1_target)[1]
    snapshot = rho_ratio_from_q1(q1_256(run_256(service_snapshot_toggle, 13)), q1_target)
    optimistic = rho_ratio_from_q1(q1_256(run_256(service_snapshot_clear_overlap, 13)), q1_target)

    check(1.4 < active_stat < 1.8, "active demux stationary branch is the factor-two live branch")
    check(abs(active_tick - 1.0006397433046454) < 1e-10, "active demux reproduces the tick-13 transient hit")
    check(blind8 > 1e60, "blind point-address reset is catastrophically too slow")
    check(blind28 > 1e14, "even optimistic blind 28-channel service is far too slow")
    check(snapshot > 1e40 and optimistic > 1e35, "syndrome snapshots do not derive the queue")

    print("\n[3] Stationary correction size")
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    print(f"  raw §5.2 gamma                         = {BASE_GAMMA:.12f}")
    print(f"  stationary exact gamma                 = {gamma_star:.12f}")
    print(f"  gamma_star / raw gamma                 = {gamma_star / BASE_GAMMA:.12f}")
    print(f"  barrier shift Delta[-ln gamma]         = {-math.log(gamma_star) + math.log(BASE_GAMMA):.12f}")
    check(abs(gamma_star / BASE_GAMMA - 0.997778) < 2e-6, "stationary closure needs the registered 0.22% correction")

    print("\n[4] Decision")
    print(
        """
  The active-address demultiplexer is operationally selected: every weaker
  reading that canon might otherwise license either fails catastrophically
  (blind reset) or by many orders (snapshot service).  The selection is not a
  new numerical fit; it is the QEC-recovery reading of the existing address
  ledger: repair acts on the decoder's tracked address, not on a blind channel
  draw or on a lossy current syndrome snapshot.

  But this does NOT resurrect the exact tick-13 branch.  The edge-sweep/latch
  dwell candidate is dead, and this audit derives no independent reason why
  the first cosmological handoff should be read exactly at clean-start tick 13.

  Decision:
      keep tick 13 demoted;
      use the stationary active-demux one-jump branch as the live route;
      make the remaining theorem the 0.22% gamma/barrier correction.

  A later independent dwell theorem could reopen the transient branch, but it
  is no longer the default route.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
