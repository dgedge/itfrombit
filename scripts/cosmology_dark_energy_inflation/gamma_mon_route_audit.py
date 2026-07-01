#!/usr/bin/env python3
"""Audit candidate derivations for the dressed-alpha monitoring clock gamma_mon.

The bridge/web calculation has reduced the dressed-alpha magnitude to one scale:
the monitored bridge relaxation clock gamma_mon when t_m != t_web.  This script
checks the obvious routes and keeps the anti-numerology discipline explicit:

  * canon-derived clocks are tested first;
  * near-hit clocks are reported, but not adopted unless they are independently
    derived by the service algebra;
  * the forbidden route, choosing gamma_mon from observed alpha, is shown only as
    the target.

Uses the Peierls-current bridge/web calculation from
dressed_alpha_bridge_web_open_system.py.  No N1=31 loop-integral substitution.
"""
from __future__ import annotations

import math

import dressed_alpha_bridge_web_open_system as bw


T_M = 1.0 / 3.0
TARGET = bw.DELTA_TARGET


def compute_quantities():
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
    target_gap = c_j * gamma_esc / TARGET
    gamma_target = gap_unit * T_M * T_M / target_gap
    return c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target


def delta_for(gamma_mon: float, c_j: float, gap_unit: float, gamma_esc: float) -> float:
    gap = gap_unit * T_M * T_M / gamma_mon
    return c_j * gamma_esc / gap


def main() -> None:
    c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target = compute_quantities()
    s_record = math.log(8.0 * 137.0)

    candidates = [
        (
            "global service clock",
            1.0,
            "canon-derived",
            "one elementary monitor interrogation per lattice tick",
        ),
        (
            "matter-sector clock",
            T_M,
            "plausible but not enough",
            "monitor at the Grover matter-hop tick",
        ),
        (
            "dephased-rate preserving",
            T_M * T_M,
            "scaling convention",
            "keeps 2|H|^2/gamma at the unit-clock value",
        ),
        (
            "R14 firing clock",
            bw.ALPHA0,
            "wrong object",
            "rate of the distinguished firing projector, not the readout/dephasing traffic",
        ),
        (
            "half unit bridge gap",
            0.5 * gap_unit,
            "near hit, refuted",
            "not a Peierls-current half-linewidth; see gamma_mon_halfwidth_refutation.py",
        ),
        (
            "two links / 12+latch",
            2.0 / 13.0,
            "near hit, refuted",
            "no one-bit Peierls latch in current bridge algebra; see dressed_alpha_bridge_latch_audit.py",
        ),
        (
            "record entropy inverse",
            1.0 / s_record,
            "near hit, underived",
            "uses s1=ln(8*137) as a clock rather than an entropy per record",
        ),
        (
            "loop scale 1/(2pi)",
            1.0 / (2.0 * math.pi),
            "near hit, not service algebra",
            "generic loop scale; including it shows the density hazard",
        ),
        (
            "FORBIDDEN alpha-fit target",
            gamma_target,
            "fit if adopted",
            "computed backwards from observed delta",
        ),
    ]

    print("GAMMA_MON ROUTE AUDIT")
    print("=" * 96)
    print("Pinned bridge/web calculation:")
    print(f"  t_m/t_web      = {T_M:.6f}")
    print(f"  Gamma_esc      = {gamma_esc:.6e} Lambda")
    print(f"  <omega>_emit   = {mean_omega:.6f} Lambda")
    print(f"  c_J            = {c_j:.6f}")
    print(f"  unit gap       = {gap_unit:.6f} Lambda")
    print(f"  target gap     = {target_gap:.6f} Lambda")
    print(f"  target gamma   = {gamma_target:.6f} Lambda (not adoptable by itself)")
    print()
    print(f"{'route':<26} {'gamma':>11} {'delta':>11} {'/target':>9} {'status':<24} note")
    near_hits = []
    canon_hits = []
    for name, gamma_mon, status, note in candidates:
        delta = delta_for(gamma_mon, c_j, gap_unit, gamma_esc)
        ratio = delta / TARGET
        near_status = "underived" in status or "not service" in status or "refuted" in status
        if abs(ratio - 1.0) < 0.10 and near_status:
            near_hits.append(name)
        if abs(ratio - 1.0) < 0.10 and status == "canon-derived":
            canon_hits.append(name)
        print(f"{name:<26} {gamma_mon:>11.6f} {delta:>11.4e} {ratio:>9.3f} {status:<24} {note}")

    print("\nInterpretation:")
    print("  Canon-derived readout clocks do not land on the required shift.")
    print("  Several simple non-canonical clocks land near it, including half the unit")
    print("  bridge gap, two links over a 12+latch sweep, and 1/(2pi).  The first")
    print("  two have now been checked directly and fail as derivations.  That is")
    print("  exactly the assignment-density warning:")
    print("  a near numeric clock is not a derivation unless the service algebra selects")
    print("  it before looking at alpha.")
    print()
    print("Half-linewidth status:")
    print("  gamma_mon = gap_unit/2 has now been checked as a Peierls-current")
    print("  critical-damping/half-linewidth rule and fails that test.  It remains")
    print("  a graph-internal near-hit only.  The web-only bridge-plus-photon pole")
    print("  calculation also fails to supply the O(1) relaxation clock.")
    print("  The 12+latch traffic-clock route has also been checked: the real one-bit")
    print("  latch is the Q3/horizon coboundary kernel, while the Peierls bridge readout")
    print("  has no analogous one-bit local latch.  The live target is therefore a new")
    print("  substrate derivation of the bridge traffic clock, or a combined")
    print("  monitor-plus-web CPTP calculation.")

    assert not canon_hits
    assert len(near_hits) >= 2
    assert abs((0.5 * gap_unit) / gamma_target - 1.0) < 0.05
    print("\nALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
