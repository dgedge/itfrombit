#!/usr/bin/env python3
"""Attack the remaining factor-1.65 in the boot residual-fault rho_Lambda route.

The existing item123_boot_residual_fault_theorem.py result is

    q_6 = 0.110 * (2/9)^64,
    rho_Lambda = q_6 * alpha_0 * Lambda_QCD^4 = 0.607 rho_obs.

This script asks whether the missing factor can be closed by a code-native
refinement rather than by an after-the-fact coefficient.  It does not promote a
new derivation.  It identifies the sharpest target: a finite-depth correction

    r_eff = (2/9) * (1 + 2^-7) = (2/9) * (129/128)

would close the gap to 0.15%, and the denominator 2^7 is adjacent to the ell=6
boot depth.  The correction is not derived here; it is a pre-registered target
for the missing boot-level Kraus/current ledger.
"""

from __future__ import annotations

import math
from fractions import Fraction


ALPHA0 = 1.0 / 137.0
LAMBDA_QCD_GEV = 0.332
P_TH = 0.110
R0 = 2.0 / 9.0
ELL = 6

M_P_GEV = 1.220890e19
HBAR_GEV_S = 6.582120e-25
H0_KM_S_MPC = 67.36
MPC_KM = 3.085678e19
OMEGA_L = 0.6847


def observed_rho_lambda() -> float:
    h = H0_KM_S_MPC / MPC_KM * HBAR_GEV_S
    return 3.0 * OMEGA_L * h * h * M_P_GEV * M_P_GEV / (8.0 * math.pi)


def rho_from(p_th: float = P_TH, r: float = R0, ell: int = ELL) -> float:
    return p_th * r ** (2**ell) * ALPHA0 * LAMBDA_QCD_GEV**4


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> None:
    print("RHO_LAMBDA BOOT GAP ATTACK")

    rho_obs = observed_rho_lambda()
    q_obs = rho_obs / (ALPHA0 * LAMBDA_QCD_GEV**4)
    rho_base = rho_from()
    gap = rho_obs / rho_base
    r_req = (q_obs / P_TH) ** (1.0 / (2**ELL))
    pth_req = q_obs / (R0 ** (2**ELL))

    print("\n[1] Baseline and exact target")
    print(f"  q_obs                         = {q_obs:.12e}")
    print(f"  rho_obs                       = {rho_obs:.12e} GeV^4")
    print(f"  rho_base = 0.110*(2/9)^64*aL4 = {rho_base:.12e} GeV^4")
    print(f"  rho_base/rho_obs              = {rho_base / rho_obs:.9f}")
    print(f"  missing factor rho_obs/base   = {gap:.9f}")
    print(f"  exact r at p_th=0.110         = {r_req:.12f}")
    print(f"  exact r/(2/9)-1               = {(r_req / R0 - 1.0) * 100:.6f}%")
    print(f"  exact p_th at r=2/9           = {pth_req:.12f}")
    print(f"  exact p_th/0.110              = {pth_req / P_TH:.9f}")
    check(0.60 < rho_base / rho_obs < 0.61, "baseline reproduces the recorded 0.607 rho_obs result")
    check(abs(gap - 1.6479565808) < 1e-6, "the factor-1.65 gap is reproduced")

    print("\n[2] Threshold routes")
    threshold_candidates = [
        ("canonical CSS/toric threshold", P_TH, "imported threshold; current canon"),
        ("1/C(4,2)", 1.0 / math.comb(4, 2), "code-native malignant-pair count, but not enough"),
        ("2/11", 2.0 / 11.0, "excellent numerically, no derivation located"),
        ("exact needed", pth_req, "target value, not a derivation"),
    ]
    for name, pth, note in threshold_candidates:
        rel = rho_from(p_th=pth) / rho_obs
        print(f"  {name:31s} p_th={pth:.12f}  rho/rho_obs={rel:.9f}  {note}")
    check(0.90 < rho_from(p_th=1.0 / 6.0) / rho_obs < 0.94, "the clean 1/C(4,2) threshold improves the result but leaves an O(8%) gap")
    check(abs(rho_from(p_th=2.0 / 11.0) / rho_obs - 1.0) < 0.005, "2/11 would close the gap numerically")
    print("  Verdict: p_th=2/11 is too good to ignore but too ungrounded to adopt.")

    print("\n[3] Finite-depth correction to r")
    exact_ratio = r_req / R0
    fd_rows = []
    for m in range(4, 11):
        for sign, label in ((+1, "1+2^-m"), (-1, "1-2^-m"), (+2, "1/(1-2^-m)")):
            if sign == +1:
                factor = 1.0 + 2.0 ** (-m)
            elif sign == -1:
                factor = 1.0 - 2.0 ** (-m)
            else:
                factor = 1.0 / (1.0 - 2.0 ** (-m))
            rel = rho_from(r=R0 * factor) / rho_obs
            fd_rows.append((abs(math.log(rel)), m, label, factor, rel))
    fd_rows.sort()
    for _, m, label, factor, rel in fd_rows[:8]:
        print(f"  {label:10s} m={m:2d}  r/(2/9)={factor:.12f}  rho/rho_obs={rel:.9f}")
    factor_129 = Fraction(129, 128)
    rho_129 = rho_from(r=R0 * float(factor_129))
    print(f"  exact r/(2/9)                 = {exact_ratio:.12f}")
    print(f"  129/128                       = {float(factor_129):.12f}")
    print(f"  relative r error              = {(float(factor_129) / exact_ratio - 1.0) * 100:.6f}%")
    print(f"  rho using (2/9)*(129/128)     = {rho_129:.12e} GeV^4")
    print(f"  rho_129/rho_obs               = {rho_129 / rho_obs:.9f}")
    check(abs(rho_129 / rho_obs - 1.0) < 0.002, "the 129/128 finite-depth correction closes rho_Lambda to 0.2%")
    check(abs(float(factor_129) / exact_ratio - 1.0) < 3e-5, "129/128 is within 3e-5 fractionally of the exact r correction")
    check(fd_rows[0][1] == 7, "the best finite-depth correction uses m=7 = ell+1")
    print("  Verdict: this is the sharpest constructive target.  It must be derived")
    print("  from the boot Kraus/current ledger as a boundary or finite-depth effect;")
    print("  otherwise it is just the target encoded as 1+2^-7.")

    print("\n[4] O(1) coefficient fitting null")
    simple_factors: list[tuple[float, str]] = [
        (5.0 / 3.0, "5/3"),
        (8.0 / 5.0, "8/5"),
        ((1.0 + math.sqrt(5.0)) / 2.0, "phi"),
        (math.pi / 2.0, "pi/2"),
        (math.sqrt(8.0 / 3.0), "sqrt(8/3)"),
    ]
    for value, name in simple_factors:
        print(f"  {name:9s} = {value:.9f}; gap/factor - 1 = {gap / value - 1.0:+.3%}")

    rational_hits = []
    for den in range(1, 17):
        for num in range(1, 4 * den + 1):
            value = num / den
            if abs(value / gap - 1.0) < 0.03:
                rational_hits.append(Fraction(num, den))
    rational_hits = sorted(set(rational_hits), key=float)
    print(f"  simple rationals p/q<=16 within 3% of the gap: {len(rational_hits)}")
    print("  examples:", ", ".join(str(x) for x in rational_hits[:12]))
    check(len(rational_hits) >= 8, "the standalone factor-1.65 alphabet is dense enough to reject coefficient shopping")

    print("\n[5] 2/9 status")
    print("  The 2/9 used here is not the retired item-86 universal banner.")
    print("  It is the local boot pair probability (1/2)*(2/3)^2:")
    print("    coin-singlet pair factor = 1/2")
    print("    walk-active exiting fraction = 2/3")
    print("    pair probability = 2/9")
    check(Fraction(1, 2) * Fraction(2, 3) ** 2 == Fraction(2, 9), "the boot 2/9 is an exact local ledger identity")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  No exact rho_Lambda derivation is closed here.")
    print("  The baseline non-horizon boot theorem remains real but conditional:")
    print("      q_6 = 0.110*(2/9)^64, rho = 0.607 rho_obs.")
    print("  Closing the coefficient by an arbitrary O(1) multiplier is rejected;")
    print("  the nearby simple-factor alphabet is too dense.")
    print("  Two concrete derivation targets survive:")
    print("    1. derive the logical-level threshold.  The code-native 1/C(4,2)")
    print("       improves the result to 0.919 rho_obs; 2/11 would close it but")
    print("       currently has no event-algebra derivation.")
    print("    2. derive a finite-depth correction to the pair ratio:")
    print("          r_eff = (2/9)*(129/128) = (2/9)*(1+2^-7),")
    print("       which gives 0.9985 rho_obs.  Since 7=ell+1, this is the")
    print("       highest-value target for the boot Kraus/current ledger.")
    print("  Required standard: derive the 129/128 or 2/11 from the microscopic")
    print("  boot ledger before comparing to rho_obs; otherwise the factor-1.65")
    print("  remains a coefficient fit.")
    print("=" * 104)
    print("exit 0 -- boot gap attacked; finite-depth 129/128 target identified, not adopted.")


if __name__ == "__main__":
    main()
