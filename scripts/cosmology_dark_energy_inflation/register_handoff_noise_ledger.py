#!/usr/bin/env python3
r"""REGISTER HANDOFF NOISE LEDGER — cosmological-constant p_c derivation attempt.

The live non-horizon rho_Lambda route is no longer a free O(1) coefficient:

    rho_Lambda = alpha_0 Lambda^4 * (21 q1)^32 / 21,

where q1 is the exact strain-decoder cell-failure probability at the boot handoff.
This script asks whether p_c can be derived from the local register thermodynamics
rather than inverted from rho_Lambda.

Outcome discipline:
  * The target p_c and its sensitivity are computed from the exact cell law.
  * Every canonical local ledger tested here is pre-declared and self-checking.
  * A near-hit attempt prefactor is recorded only as a theorem target, not adopted.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np


ALPHA0 = 1.0 / 137.0
LAMBDA_QCD_GEV = 0.332
M_P_GEV = 1.220890e19
HBAR_GEV_S = 6.582120e-25
H0_KM_S_MPC = 67.36
MPC_KM = 3.085678e19
OMEGA_L = 0.6847

PHI = (math.sqrt(5.0) - 1.0) / 2.0
KAPPA = 1.0 / (2.0 * PHI)  # ANCHOR §5.2 Group-II exponent.
BASE_DELTA = 3.0 * KAPPA
BASE_ODDS = math.exp(-BASE_DELTA)

VERTICES = range(8)
EDGES = [
    (u, v)
    for u in VERTICES
    for v in VERTICES
    if u < v and bin(u ^ v).count("1") == 1
]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def observed_rho_lambda() -> float:
    h = H0_KM_S_MPC / MPC_KM * HBAR_GEV_S
    return 3.0 * OMEGA_L * h * h * M_P_GEV * M_P_GEV / (8.0 * math.pi)


def cell_failure_weight(k: int) -> float:
    if k <= 3:
        return 0.0
    if k == 4:
        return 0.5
    return 1.0


def q1_iid(p: float) -> float:
    return sum(
        math.comb(8, k) * p**k * (1.0 - p) ** (8 - k) * cell_failure_weight(k)
        for k in range(9)
    )


def q_top_from_p(p: float) -> float:
    q1 = q1_iid(p)
    return (21.0 * q1) ** 32 / 21.0


def rho_ratio_from_p(p: float) -> float:
    return q_top_from_p(p) * ALPHA0 * LAMBDA_QCD_GEV**4 / observed_rho_lambda()


def solve_target_p() -> tuple[float, float, float]:
    q_target = observed_rho_lambda() / (ALPHA0 * LAMBDA_QCD_GEV**4)
    q1_target = math.exp((math.log(21.0) + math.log(q_target)) / 32.0 - math.log(21.0))
    lo, hi = 0.01, 0.30
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if q1_iid(mid) < q1_target:
            lo = mid
        else:
            hi = mid
    p_target = (lo + hi) / 2.0
    eps = 1.0e-7
    sensitivity = (
        32.0
        * p_target
        * (q1_iid(p_target + eps) - q1_iid(p_target - eps))
        / (2.0 * eps)
        / q1_iid(p_target)
    )
    return p_target, q1_target, sensitivity


def strain_frustration(mask: int) -> int:
    return sum(((mask >> u) & 1) ^ ((mask >> v) & 1) for u, v in EDGES)


def hamming_weight(mask: int) -> int:
    return int(mask).bit_count()


def stationary_ctmc(mode: str) -> tuple[float, float]:
    """Stationary distribution for exact 256-state local bit-fault ledgers.

    Modes:
      detailed_balance:
        reversible heat-bath flips with strain energy F/(2 phi). Complement
        symmetry is exact and forces p=1/2.
      isolated_service:
        independent bit queue. A clean bit activates with the isolated single-fault
        Boltzmann odds exp(-3/(2 phi)); an active fault is serviced at unit rate.
      positive_strain_birth:
        generation of a clean-bit fault is penalised only by positive incremental
        strain, exp[-kappa max(Delta F, 0)], while active faults are serviced at
        unit rate. This is the most generous local service-biased strain ledger.
    """

    q = np.zeros((256, 256), dtype=float)
    for state in range(256):
        f0 = strain_frustration(state)
        for bit in VERTICES:
            target = state ^ (1 << bit)
            df = strain_frustration(target) - f0
            occupied = (state >> bit) & 1
            if mode == "detailed_balance":
                rate = 1.0 / (1.0 + math.exp(KAPPA * df))
            elif mode == "isolated_service":
                rate = BASE_ODDS if not occupied else 1.0
            elif mode == "positive_strain_birth":
                rate = math.exp(-KAPPA * max(df, 0)) if not occupied else 1.0
            else:
                raise ValueError(mode)
            q[state, target] = rate
        q[state, state] = -q[state].sum()

    a = q.T.copy()
    b = np.zeros(256)
    a[-1, :] = 1.0
    b[-1] = 1.0
    pi = np.linalg.solve(a, b)
    p = sum(float(pi[state]) * hamming_weight(state) / 8.0 for state in range(256))
    q1 = sum(float(pi[state]) * cell_failure_weight(hamming_weight(state)) for state in range(256))
    return p, q1


def prefactor_scan(p_target: float) -> list[tuple[float, Fraction, float, float]]:
    hits: list[tuple[float, Fraction, float, float]] = []
    for den in range(1, 33):
        for num in range(1, 4 * den + 1):
            factor = Fraction(num, den)
            odds = BASE_ODDS * float(factor)
            p = odds / (1.0 + odds)
            ratio = rho_ratio_from_p(p)
            if 0.5 < ratio < 2.0:
                hits.append((abs(math.log(ratio)), factor, p, ratio))
    return sorted(set(hits))


def main() -> None:
    print("REGISTER HANDOFF NOISE LEDGER")

    print("\n[1] Target from the exact cell law")
    p_target, q1_target, sensitivity = solve_target_p()
    odds_target = p_target / (1.0 - p_target)
    delta_target = math.log((1.0 - p_target) / p_target)
    print(f"  q1_target                         = {q1_target:.12e}")
    print(f"  p_c target                         = {p_target:.12f}")
    print(f"  odds p/(1-p)                       = {odds_target:.12f}")
    print(f"  DeltaF/T = ln((1-p)/p)             = {delta_target:.12f}")
    print(f"  exact-law sensitivity dlnrho/dlnp  = {sensitivity:.3f}")
    print(f"  rho within x2 requires p to        = +/-{100.0 * math.log(2.0) / sensitivity:.3f}%")
    check(abs(p_target - 0.0972) < 6e-4, "target agrees with the registered p_c ~= 0.0972")
    check(118.0 < sensitivity < 122.0, "exact-law sensitivity is about 120, not an O(1) tolerance")

    print("\n[2] Canonical local ledgers")
    rows = []
    for mode in ("detailed_balance", "isolated_service", "positive_strain_birth"):
        p, q1 = stationary_ctmc(mode)
        rho_ratio = (21.0 * q1) ** 32 / 21.0 * ALPHA0 * LAMBDA_QCD_GEV**4 / observed_rho_lambda()
        rows.append((mode, p, q1, rho_ratio))
        print(f"  {mode:24s} p={p:.12f} q1={q1:.12e} rho/rho_obs={rho_ratio:.3e}")
    check(abs(rows[0][1] - 0.5) < 1e-12, "reversible strain thermodynamics is complement-symmetric and fails")
    check(0.080 < rows[1][1] < 0.082, "isolated service queue reproduces the normalized single-fault Boltzmann value")
    check(rows[2][1] > 0.20, "positive-strain local birth rule overshoots badly")

    raw_boltzmann_p = BASE_ODDS
    print("\n[3] Isolated single-fault Boltzmann readings")
    print(f"  ANCHOR §5.2 single fault: DeltaF=3, kappa=1/(2phi) -> exp(-3/(2phi)) = {BASE_ODDS:.12f}")
    print(f"  raw event-probability reading       p={raw_boltzmann_p:.12f} rho/rho_obs={rho_ratio_from_p(raw_boltzmann_p):.3e}")
    queue_p = BASE_ODDS / (1.0 + BASE_ODDS)
    print(f"  normalized occupancy/odds reading   p={queue_p:.12f} rho/rho_obs={rho_ratio_from_p(queue_p):.3e}")
    check(rho_ratio_from_p(raw_boltzmann_p) < 1e-4, "raw Boltzmann reading undershoots rho by at least four orders")
    check(rho_ratio_from_p(queue_p) < 1e-8, "normalized queue reading undershoots even more strongly")

    print("\n[4] Attempt-rate prefactor target")
    factor_target = odds_target / BASE_ODDS
    print(f"  required attempt factor A_req       = {factor_target:.12f}")
    print(f"  ln A_req                            = {math.log(factor_target):.12f}")
    candidates = [
        ("sqrt(3/2)   (walk-active/exiting Hessian target)", math.sqrt(3.0 / 2.0)),
        ("39/32       (nearer simple rational, no source)", 39.0 / 32.0),
        ("28/23       (simple rational, no source)", 28.0 / 23.0),
        ("11/9        (simple rational, no source)", 11.0 / 9.0),
        ("17/14       (simple rational, no source)", 17.0 / 14.0),
    ]
    for name, factor in candidates:
        odds = BASE_ODDS * factor
        p = odds / (1.0 + odds)
        print(f"  {name:46s} A={factor:.12f} p={p:.12f} rho/rho_obs={rho_ratio_from_p(p):.3f}")
    hits = prefactor_scan(p_target)
    print(f"  rationals p/q, q<=32 giving rho within x2: {len(hits)}")
    for _, factor, p, ratio in hits[:8]:
        print(f"    {str(factor):>5s}: p={p:.12f}, rho/rho_obs={ratio:.3f}")
    check(len(hits) >= 5, "nearby rational prefactors are dense enough to block coefficient adoption")

    sqrt_factor = math.sqrt(3.0 / 2.0)
    sqrt_odds = BASE_ODDS * sqrt_factor
    sqrt_p = sqrt_odds / (1.0 + sqrt_odds)
    check(0.5 < rho_ratio_from_p(sqrt_p) < 2.0, "sqrt(3/2) is a real target-grade near hit")
    check(any(factor == Fraction(39, 32) for _, factor, _, _ in hits), "a source-free rational beats or rivals the near hit")

    print(
        """
VERDICT
  The register-handoff ledger has not derived the cosmological constant.

  What is now closed:
    * the required local target is p_c = %.6f, DeltaF/T = %.6f;
    * exact-law sensitivity is d ln rho / d ln p = %.1f, so this is a
      half-percent target, not an adjustable coefficient;
    * reversible Q3 strain thermodynamics, isolated §5.2 Boltzmann, and the
      simplest service-biased 256-state ledgers all fail, with reasons.

  What survives as the sharp theorem target:
    * derive an attempt-rate prefactor A = %.6f multiplying the §5.2 isolated
      single-fault odds exp[-3/(2phi)].
    * sqrt(3/2) is structurally tempting because 3 walk-active legs and 2 exiting
      vacuum legs are already canonical, and it lands within a factor two in rho.
      It is NOT adopted: 39/32 and other small rationals occupy the same band.

  Promotion criterion:
    Build the microscopic service-current/Hessian theorem that fixes A before
    looking at rho_Lambda. If it gives A = A_req to ~0.6%%, the CC route closes;
    if it gives the uncorrected Boltzmann value or any generic local ledger above,
    the non-horizon cosmological-constant route fails by computation.
exit 0"""
        % (p_target, delta_target, sensitivity, factor_target)
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
