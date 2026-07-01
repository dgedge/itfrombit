#!/usr/bin/env python3
"""Substrate traffic-clock theorem attempt for the dressed-alpha bridge.

Question:
    Can the current substrate service algebra derive the missing bridge monitor
    clock gamma_mon without fitting alpha?

The available local theorem is the item-131 one-jump service clock:

    one covariant service event per substrate tick over N equivalent channels
    => per-channel service rate 1/N.

For a two-link Peierls readout one can also test a workload aggregation

    two active bridge links over an N-slot service sweep => gamma = 2/N.

This script checks only alphabets already present in the local canon: bridge
links/endpoints, cubic directions, Q3 vertices, Q3 strain edges, the
12-edge-plus-latch handoff candidate, and the 28 logical/transverse service
alphabet.  It reports near hits but does not promote them unless the same canon
also licenses the bridge-readout identification.  The bridge-local latch audit
now refutes the 12+latch identification under the current Peierls algebra.
"""
from __future__ import annotations

from dataclasses import dataclass

import dressed_alpha_bridge_web_open_system as bw


T_M = 1.0 / 3.0


@dataclass(frozen=True)
class ClockCandidate:
    name: str
    gamma_mon: float
    status: str
    note: str


def compute_target():
    h2, _pairs, idx, _bas = bw.build_pair_system()
    s_eta, _omega_max = bw.photon_form_factor()
    c_j, gap_unit = bw.response_coefficient(h2, idx)
    gamma_esc, mean_omega, _weight = bw.gamma_escape(
        h2,
        idx,
        s_eta,
        t_m=T_M,
        eta=bw.ETA_PIN,
        g_portal=bw.G_PORTAL,
    )
    target_gap = c_j * gamma_esc / bw.DELTA_TARGET
    gamma_target = gap_unit * T_M * T_M / target_gap
    return c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target


def delta_ratio(gamma_mon: float, gamma_target: float) -> float:
    """In this reduced response law delta is linear in gamma_mon."""
    return gamma_mon / gamma_target


def candidates() -> list[ClockCandidate]:
    return [
        ClockCandidate(
            "native service tick",
            1.0,
            "canon, wrong object",
            "one global substrate service event per tick; already too fast for this bridge clock",
        ),
        ClockCandidate(
            "Grover matter tick",
            T_M,
            "canon, wrong scale",
            "matter hopping scale t_m/t_web=1/3, not a monitor scheduler",
        ),
        ClockCandidate(
            "two bridge links",
            1.0 / 2.0,
            "underdefined",
            "serial service over only the two Peierls links",
        ),
        ClockCandidate(
            "four bridge endpoints",
            1.0 / 4.0,
            "underdefined",
            "serial service over the four endpoint incidences of the two links",
        ),
        ClockCandidate(
            "six cubic directions",
            1.0 / 6.0,
            "portable one-jump",
            "one direction channel out of six; no proof this is the bridge monitor alphabet",
        ),
        ClockCandidate(
            "seven incident hyperplanes",
            1.0 / 7.0,
            "wrong alphabet",
            "point-hyperplane incidence degree from the 8->28 instrument, not a bridge readout",
        ),
        ClockCandidate(
            "eight Q3 point labels",
            1.0 / 8.0,
            "canon alphabet",
            "single-bit Kraus alphabet before incidence refinement",
        ),
        ClockCandidate(
            "two links / 12 strain edges",
            2.0 / 12.0,
            "tempting workload",
            "two Peierls links inside the covariant 12-edge Q3 strain sweep",
        ),
        ClockCandidate(
            "two links / (12+latch)",
            2.0 / 13.0,
            "near hit, refuted",
            "bridge-local latch audit finds no one-bit Peierls latch",
        ),
        ClockCandidate(
            "twelve strain edges",
            1.0 / 12.0,
            "canon alphabet",
            "full Q3 strain-ledger sweep, single edge service",
        ),
        ClockCandidate(
            "twenty-eight service channels",
            1.0 / 28.0,
            "canon theorem",
            "item-131 logical/transverse service alphabet with one-jump gap 1/28",
        ),
    ]


def main() -> None:
    c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target = compute_target()
    target_single_n = 1.0 / gamma_target
    target_two_link_sweep = 2.0 / gamma_target

    print("DRESSED-ALPHA SUBSTRATE TRAFFIC-CLOCK AUDIT")
    print("=" * 96)
    print("Pinned bridge/web target:")
    print(f"  t_m/t_web                    = {T_M:.6f}")
    print(f"  Gamma_esc                    = {gamma_esc:.6e} Lambda")
    print(f"  <omega>_emit                 = {mean_omega:.6f} Lambda")
    print(f"  c_J                          = {c_j:.6f}")
    print(f"  unit bridge gap              = {gap_unit:.6f} Lambda")
    print(f"  required response gap        = {target_gap:.6f} Lambda")
    print(f"  required gamma_mon           = {gamma_target:.6f} Lambda")
    print(f"  as single-channel 1/N        = N_eff {target_single_n:.6f}")
    print(f"  as two-link workload 2/N     = N_eff {target_two_link_sweep:.6f}")
    print()
    print("Canon traffic-clock candidates:")
    print(f"  {'route':<30} {'gamma':>11} {'delta/target':>13} {'status':<22} note")

    exact_hits = []
    near_hits = []
    for cand in candidates():
        ratio = delta_ratio(cand.gamma_mon, gamma_target)
        if abs(ratio - 1.0) < 0.005 and "canon" in cand.status:
            exact_hits.append(cand.name)
        if abs(ratio - 1.0) < 0.05:
            near_hits.append(cand.name)
        print(f"  {cand.name:<30} {cand.gamma_mon:>11.6f} {ratio:>13.6f} {cand.status:<22} {cand.note}")

    print()
    print("Interpretation:")
    print("  The one-jump theorem can derive 1/N once the bridge monitor alphabet is")
    print("  fixed, but the required effective N is not a native integer: 6.416 for a")
    print("  single-channel reading, or 12.833 for a two-link workload reading.")
    print("  The closest portable integer clocks are N=6 and N=7; both miss by about")
    print("  seven to eight percent.  The closest workload candidate is two links over")
    print("  a 12-edge-plus-latch cycle, 2/13, but that latch is a cosmological handoff")
    print("  candidate and the bridge-local audit finds no analogous one-bit Peierls")
    print("  latch.  The eta branch is a static Wilson/Peierls phase, not a service")
    print("  slot.")
    print()
    print("VERDICT")
    print("  No substrate traffic-clock theorem currently derives gamma_mon.  The")
    print("  12+latch bridge-reading target has now been tested and refuted under the")
    print("  existing bridge algebra: the known latch belongs to the Q3/horizon")
    print("  coboundary kernel, not to the Peierls current readout.  Adopting 2/13,")
    print("  1/6, or 1/7 would therefore be another assignment choice.")

    assert not exact_hits
    assert "two links / (12+latch)" in near_hits
    assert abs((2.0 / 13.0) / gamma_target - 1.0) > 0.01
    assert 6.0 < target_single_n < 7.0
    assert 12.0 < target_two_link_sweep < 13.0
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
