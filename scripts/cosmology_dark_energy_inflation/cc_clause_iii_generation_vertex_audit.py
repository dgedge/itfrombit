#!/usr/bin/env python3
r"""CC CLAUSE III: can the 0.222% stationary gamma correction be derived?

Live route before this audit:
  stationary active-demux one-jump branch = 1.616 rho_obs
  exact closure needs gamma_star/gamma_Sec5.2 = 0.997778

Clause iii asks whether the small correction is an output of existing framework
machinery, especially the dressed-alpha / gauge-web vacuum-polarization sector.

This script separates three logically different objects:
  (1) the required generation-vertex barrier Delta(-ln gamma);
  (2) alpha0/pi-like radiative candidates, which are numerically close;
  (3) the actual canon dressed-alpha machinery, which K5 reclassified as a
      count-plus-ansatz rather than a loop integral.

exit 0 means the numerical target and the source-status verdict are reproduced.
It does NOT mean the cosmological constant is derived.
"""

from __future__ import annotations

import math

from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    queue_readouts,
    rho_ratio_from_q1,
    solve_gamma_for_stationary_queue,
    solve_target,
)


CANON_DRESSED_ALPHA_INV = 137.035999077
GENUINE_TCH_LOOP_PI = 0.06
COUNT_LOOP_SLOT = 31.0 / (2.0 * math.pi)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def rho_for_delta(delta: float, q1_target: float) -> float:
    """Stationary one-jump rho ratio after adding Delta to -ln(gamma)."""
    gamma = BASE_GAMMA * math.exp(-delta)
    _, q_post, _ = queue_readouts(gamma, 1)
    return rho_ratio_from_q1(q_post, q1_target)


def main() -> None:
    print("CC CLAUSE-III GENERATION-VERTEX AUDIT")

    _, q1_target = solve_target()
    _, q1_stat, mean_active = queue_readouts(BASE_GAMMA, 1)
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    delta_req = -math.log(gamma_star / BASE_GAMMA)
    c_req = delta_req / ALPHA0
    rho_raw = rho_ratio_from_q1(q1_stat, q1_target)

    print("\n[1] Target fixed by the stationary active-demux branch")
    print(f"  raw gamma_Sec5.2                       = {BASE_GAMMA:.12f}")
    print(f"  stationary post-service q1             = {q1_stat:.12e}")
    print(f"  mean active bit fraction               = {mean_active:.12f}")
    print(f"  raw stationary rho/rho_obs             = {rho_raw:.6f}")
    print(f"  gamma_star/gamma_Sec5.2                = {gamma_star / BASE_GAMMA:.12f}")
    print(f"  required Delta(-ln gamma)              = {delta_req:.12f}")
    print(f"  required coefficient C=Delta/alpha0    = {c_req:.12f}")
    check(1.55 < rho_raw < 1.70, "baseline stationary branch is the 1.616 rho_obs route")
    check(abs(gamma_star / BASE_GAMMA - 0.997778) < 5e-6, "closure requires the registered 0.222% gamma reduction")
    check(abs(delta_req - 0.0022245) < 5e-6, "Delta target matches the 256-state correction scan")

    print("\n[2] Candidate coefficients at the generation vertex")
    candidates = [
        ("alpha0/pi", ALPHA0 / math.pi, "near one-loop barrier form"),
        ("alpha0/3", ALPHA0 / 3.0, "nearest simple trivalent divisor"),
        ("q1_stat", q1_stat, "self-consistency feedback"),
        (
            "dressed-alpha fractional shift",
            (ALPHA0 - 1.0 / CANON_DRESSED_ALPHA_INV) / ALPHA0,
            "actual 137 -> 137.036 fractional alpha shift",
        ),
    ]
    print(f"  {'candidate':<32s} {'Delta':>13s} {'C=Delta/a0':>13s} {'rho/rho_obs':>13s}  status")
    for name, delta, note in candidates:
        rho = rho_for_delta(delta, q1_target)
        coeff = delta / ALPHA0
        miss = coeff / c_req - 1.0
        status = (
            "2% hit"
            if abs(math.log(rho)) < 0.02
            else "5% hit"
            if abs(math.log(rho)) < 0.05
            else "miss"
        )
        print(f"  {name:<32s} {delta:13.9f} {coeff:13.6f} {rho:13.6f}  {status}; {note}; C miss {miss:+.2%}")

    alpha_pi_rho = rho_for_delta(ALPHA0 / math.pi, q1_target)
    alpha_third_rho = rho_for_delta(ALPHA0 / 3.0, q1_target)
    q1_rho = rho_for_delta(q1_stat, q1_target)
    dressed_frac = (ALPHA0 - 1.0 / CANON_DRESSED_ALPHA_INV) / ALPHA0
    dressed_frac_rho = rho_for_delta(dressed_frac, q1_target)
    check(0.975 < alpha_pi_rho < 0.985, "alpha0/pi is the closest structured near-hit, but just outside a 2%-in-rho gate")
    check(0.94 < alpha_third_rho < 0.97, "alpha0/3 remains a second 5%-class near-hit")
    check(q1_rho < 0.96, "q1 self-consistency undershoots too much to close clause iii")
    check(dressed_frac < delta_req / 5.0, "the literal dressed-alpha fractional shift is far too small as the barrier")
    check(dressed_frac_rho > 1.45, "literal dressed-alpha shift leaves the CC branch high")

    print("\n[3] Source test for the dressed-alpha route")
    print(f"  canon dressed-alpha inverse             = {CANON_DRESSED_ALPHA_INV:.9f}")
    print(f"  fractional alpha shift                  = {dressed_frac:.9f}")
    print(f"  required barrier / fractional shift     = {delta_req / dressed_frac:.3f}")
    print(f"  count-loop slot 31/(2pi)                = {COUNT_LOOP_SLOT:.6f}")
    print(f"  genuine TCH loop Pi (K5 audit scale)    ~ {GENUINE_TCH_LOOP_PI:.3f}")
    print(f"  count slot / genuine loop scale         ~ {COUNT_LOOP_SLOT / GENUINE_TCH_LOOP_PI:.1f}")
    check(COUNT_LOOP_SLOT / GENUINE_TCH_LOOP_PI > 50.0, "K5 mismatch: count slot is not the genuine loop integral")

    print(
        """
VERDICT
  Clause iii is attacked but not closed.

  The alpha0/pi branch is the best structured candidate numerically:
      Delta = alpha0/pi gives rho = 0.979 rho_obs and
      C = 1/pi = 0.31831 versus C_req = 0.30475.
  But the existing dressed-alpha machinery does not derive this generation-vertex
  coefficient. Canon K5 says the 137.036 result is a count-plus-self-consistency
  construction, not an evaluated vacuum-polarization integral; the genuine loop
  is operator/mass dependent and sits in a different numerical slot. Reusing it
  here would be exactly the thermodynamic-cluster failure mode: a plausible
  coefficient without the observable map.

  The q1 self-consistency check also fails as a closure: it gives rho = 0.944
  rho_obs and has no derived feedback map.

  The remaining theorem target is now exact:
      derive a generation-vertex operator/current whose loop or local
      normalization shifts -ln(gamma) by
          Delta_vertex = 0.0022245 = alpha0 * 0.30475,
      not merely a dressed-alpha shift in alpha^{-1}.
  Until that operator and observable map are supplied, the live CC route remains
  stationary active-demux + open 0.222% barrier correction, not a Locked
  cosmological constant.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
