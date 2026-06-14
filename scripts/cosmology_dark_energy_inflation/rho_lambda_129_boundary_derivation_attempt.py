#!/usr/bin/env python3
"""Attempt to derive the 129/128 finite-depth correction from the boot ledger.

Question
--------
Does the actual boot Kraus/current ledger force

    r_eff = (2/9) * (129/128)

for the ell=6 residual-fault rho_Lambda route?

Answer from this audit
----------------------
No.  The available ledger gives an exact local pair probability r=2/9 and,
when tensorized over 2^6 malignant pair slots, no finite-depth correction to
the mean.  An open-boundary chain of adjacent pair windows does produce a real
finite-size current-covariance effect, but it changes the Fano factor, not the
mean current.  The 129/128 factor appears only if one adds an extra boundary
leg/half-interval by hand to the 128 leg slots.  That is a clean target for a
future boot-boundary theorem, not a derivation from the current ledger.
"""

from __future__ import annotations

import math
from fractions import Fraction


ELL = 6
N_PAIR = 2**ELL
N_LEG = 2 * N_PAIR
R_PAIR = Fraction(2, 9)
P_TH = Fraction(11, 100)
ALPHA0 = 1.0 / 137.0
LAMBDA_QCD_GEV = 0.332

M_P_GEV = 1.220890e19
HBAR_GEV_S = 6.582120e-25
H0_KM_S_MPC = 67.36
MPC_KM = 3.085678e19
OMEGA_L = 0.6847


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def observed_rho_lambda() -> float:
    h = H0_KM_S_MPC / MPC_KM * HBAR_GEV_S
    return 3.0 * OMEGA_L * h * h * M_P_GEV * M_P_GEV / (8.0 * math.pi)


def rho_from(r: float) -> float:
    return float(P_TH) * r**N_PAIR * ALPHA0 * LAMBDA_QCD_GEV**4


def main() -> None:
    print("RHO_LAMBDA 129/128 BOUNDARY-DERIVATION ATTEMPT")

    print("\n[1] Local pair Kraus ledger")
    walk_active = ("C0", "C1", "I3")
    exits = {"C0", "C1"}
    outcomes = []
    for singlet in (0, 1):
        for left in walk_active:
            for right in walk_active:
                event = singlet == 1 and left in exits and right in exits
                outcomes.append((singlet, left, right, event))
    total = len(outcomes)
    event_count = sum(1 for *_, event in outcomes if event)
    r = Fraction(event_count, total)
    print(f"  Kraus labels: singlet bit x 3 x 3 = {total}")
    print(f"  successful pair labels: {event_count}")
    print(f"  r = {event_count}/{total} = {r}")
    check(total == 18, "pair ledger has 18 exclusive labels")
    check(r == R_PAIR, "local boot pair probability is exactly 2/9")

    print("\n[2] Tensorized finite-depth boot ledger")
    tensor_correction = Fraction(1, 1)
    print(f"  ell = {ELL}; malignant pair slots N = 2^ell = {N_PAIR}; leg slots = {N_LEG}")
    print(f"  tensor trace ratio = r^N = (2/9)^{N_PAIR}; correction factor = {tensor_correction}")
    check(tensor_correction == 1, "independent tensor Kraus ledger gives no finite-depth mean correction")

    rho_obs = observed_rho_lambda()
    rho_base = rho_from(float(R_PAIR))
    print(f"  rho_base/rho_obs = {rho_base / rho_obs:.9f}")
    check(0.60 < rho_base / rho_obs < 0.61, "unmodified tensor ledger reproduces the known 0.607 rho_obs")

    print("\n[3] Open-chain current ledger: boundary covariance, not mean shift")
    # Pair windows share leg variables X_i on an open chain:
    # E_i = S_i X_i X_{i+1}; P(S_i)=1/2, P(X_i=exit)=2/3.
    p_x = Fraction(2, 3)
    p_s = Fraction(1, 2)
    mean_pair = p_s * p_x * p_x
    adjacent_second = p_s * p_s * p_x * p_x * p_x
    cov_adjacent = adjacent_second - mean_pair * mean_pair
    var_pair = mean_pair * (1 - mean_pair)
    total_mean = N_PAIR * mean_pair
    total_var = N_PAIR * var_pair + 2 * (N_PAIR - 1) * cov_adjacent
    fano = total_var / total_mean
    print(f"  open chain leg slots = N+1 = {N_PAIR + 1}")
    print(f"  mean per adjacent pair = {mean_pair}")
    print(f"  adjacent covariance    = {cov_adjacent}")
    print(f"  total mean             = {total_mean}")
    print(f"  total Fano             = {fano} = {float(fano):.9f}")
    check(mean_pair == R_PAIR, "open-boundary sharing leaves the mean pair probability at 2/9")
    check(cov_adjacent == Fraction(2, 81), "shared boundary leg creates a real positive adjacent covariance")
    check(fano == Fraction(287, 288), "for ell=6 the open-chain count has Fano factor 287/288")
    print("  The actual boundary effect found here is T4/covariance data, not a")
    print("  T3 mean-current correction. It cannot supply 129/128 for rho_Lambda.")

    print("\n[4] Where 129/128 enters")
    target = Fraction(129, 128)
    unsupported_leg_factor = Fraction(N_LEG + 1, N_LEG)
    unsupported_half_interval = Fraction(N_PAIR * 2 + 1, N_PAIR * 2)
    print(f"  (2N+1)/(2N) with N=64 = {unsupported_leg_factor}")
    print(f"  half-interval form     = {unsupported_half_interval}")
    check(unsupported_leg_factor == target, "129/128 is exactly one extra boundary leg on top of 128 leg slots")
    rho_129 = rho_from(float(R_PAIR * target))
    print(f"  rho with r*(129/128)  = {rho_129 / rho_obs:.9f} rho_obs")
    check(abs(rho_129 / rho_obs - 1.0) < 0.002, "the unsupported boundary leg would close the number")
    print("  But the local Kraus completeness relation supplies 2N leg labels, not")
    print("  2N+1, and the open-chain ledger supplies N+1 shared leg variables while")
    print("  preserving the same mean per pair. No current term in the constructed")
    print("  ledger licenses the extra '+1' leg.")

    print("\n[5] Verdict under T1-T9")
    print("  T1 event algebra: local pair Kraus labels built; finite-depth tensor")
    print("     and open-chain variants are explicit.")
    print("  T3 mean current: remains r=2/9. The 129/128 mean shift is NOT derived.")
    print("  T4 covariance: partially improved; open-chain boundary covariance gives")
    print("     Fano=287/288 at ell=6.")
    print("  T8 alternatives: unsupported extra boundary leg is identified as the")
    print("     only tested route that yields 129/128.")

    print("\n" + "=" * 100)
    print("VERDICT")
    print("  The attempted derivation fails in the strict sense.")
    print("  The actual boot Kraus/current ledger forces the local pair probability")
    print("  r=2/9 and tensorizes over 64 pair slots without a mean correction.")
    print("  Open finite-depth boundaries create a real covariance correction")
    print("  (Fano=287/288), but no 129/128 shift in the mean current.")
    print("  The 129/128 factor is equivalent to adding one unlicensed boundary leg")
    print("  to the 128 leg slots. It remains a target theorem: prove that extra")
    print("  boundary leg/half-interval exists in the microscopic boot process, or")
    print("  reject the exact rho_Lambda closure.")
    print("=" * 100)
    print("exit 0 -- 129/128 not derived; covariance target sharpened.")


if __name__ == "__main__":
    main()
