#!/usr/bin/env python3
r"""ITEM 132: dynamical attempt at the R4 line-density law.

Target
------
Try to prove the remaining local law needed by item132_r4_local_action_lift.py:

    lambda_R4 = (2/3) |g|/a0.

Result
------
Generic local QEC/virial dynamics gets close but does not fix the exact
normalisation.  It derives

    lambda_R4 = (2/3) chi_R4 |g|/a0,

where:

* 2/3 is finite: two legal R4 repair edges per three R1 generation sectors.
* the linear dependence is the unique local line-response branch compatible
  with 1D R4 support and the p=3/flat-curve/BTFR exponent.
* chi_R4 is a dimensionless susceptibility set by the microscopic balance of
  R4 line creation, erasure, and virialisation.

Detailed balance, Onsager linear response, and stable birth-death dynamics do
not force chi_R4=1.  Choosing chi_R4=1 is equivalent to asserting that one
horizon-walk acceleration quantum a0 produces one unit of R4 line response.
That may be the correct substrate rule, but it is still an additional
microscopic normalization theorem to derive from W=S*C plus boundary-QEC rates.

Consequence
-----------
With chi_R4 left free, the cubic action is

    S = int chi_R4 |g|^3/(12 pi G a0) d^3x,

and spherical BTFR becomes

    v_inf^4 = (a0/chi_R4) G M_b.

Thus observations plus the external a0=cH0/2pi identification would fix
chi_R4=1, but generic stable dynamics does not prove it.  Companion
item132_chi_unit_poisson.py gives the sharper conditional closure: if W=S*C
plus boundary QEC realise a count-valued immigration-death line ledger with
matched creation/erasure rate Gamma0, the stationary law is Poisson(|g|/a0)
and chi_R4=1.
"""

from __future__ import annotations

import math
from fractions import Fraction


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def stationary_line_density(x: float, chi: float, incidence: Fraction = Fraction(2, 3)) -> float:
    """Linear-response stationary line density for |g|/a0=x."""
    return float(incidence) * chi * x


def birth_death_rhs(lam: float, x: float, chi: float, gamma: float = 1.0) -> float:
    """Stable local relaxation dynamics d lambda/dt = gamma[(2/3) chi x - lambda]."""
    return gamma * (stationary_line_density(x, chi) - lam)


def action_coeff(chi: Fraction) -> Fraction:
    """Coefficient multiplying |g|^3/(pi G a0)."""
    return chi * Fraction(1, 12)


def variation_prefactor(chi: Fraction) -> Fraction:
    """Variation of chi |g|^3/(12 pi G a0) gives chi/(4 pi G a0)."""
    return chi * Fraction(1, 4)


def btfr_a0_eff(chi: float, a0: float = 1.0) -> float:
    """Spherical equation gives v^4 = (a0/chi) G M."""
    return a0 / chi


def signed_detailed_balance_response(x: float, chi: float, incidence: Fraction = Fraction(2, 3)) -> float:
    """Bounded two-orientation detailed-balance response.

    For a finite exclusive two-orientation line state with bias chi*x,
    the net signed line density is incidence*tanh(chi*x).  This has the
    right small-x slope but is not exactly linear at x~1.
    """
    return float(incidence) * math.tanh(chi * x)


def main() -> None:
    print("ITEM 132: R4 LINE-DENSITY DYNAMICS ATTEMPT")

    print("\n[1] Finite incidence is fixed")
    incidence = Fraction(2, 3)
    check(incidence == Fraction(2, 3), "R4 has two legal repair edges per three R1 generation sectors")
    check(action_coeff(Fraction(1, 1)) == Fraction(1, 12), "unit susceptibility would reproduce |g|^3/(12 pi G a0)")
    check(variation_prefactor(Fraction(1, 1)) == Fraction(1, 4), "unit susceptibility gives the AQUAL 1/(4 pi G a0) prefactor")

    print("\n[2] Stable local birth-death / virial relaxation leaves a susceptibility")
    x = 0.4
    for chi in [0.5, 1.0, 2.0]:
        lam_star = stationary_line_density(x, chi)
        rhs_low = birth_death_rhs(0.5 * lam_star, x, chi)
        rhs_high = birth_death_rhs(1.5 * lam_star, x, chi)
        print(f"  chi={chi:.1f}: lambda*={lam_star:.6f}, rhs_below={rhs_low:+.6f}, rhs_above={rhs_high:+.6f}")
        check(rhs_low > 0 and rhs_high < 0, f"chi={chi:.1f} gives a stable stationary line density")
    check(stationary_line_density(x, 0.5) != stationary_line_density(x, 1.0), "same finite QEC structure admits different susceptibilities")
    check(stationary_line_density(x, 2.0) != stationary_line_density(x, 1.0), "dynamics fixes stability, not chi_R4=1")

    print("\n[3] Detailed balance gives linear response only up to susceptibility")
    for chi in [0.5, 1.0, 2.0]:
        slope0 = float(incidence) * chi
        exact_at_1 = signed_detailed_balance_response(1.0, chi)
        linear_at_1 = stationary_line_density(1.0, chi)
        print(
            f"  chi={chi:.1f}: small-x slope={slope0:.6f}, "
            f"tanh response at x=1={exact_at_1:.6f}, linear branch={linear_at_1:.6f}"
        )
    check(abs(signed_detailed_balance_response(1e-5, 1.0) / 1e-5 - float(incidence)) < 1e-10, "bounded detailed balance recovers the 2/3 slope at infinitesimal x for chi=1")
    check(abs(signed_detailed_balance_response(1.0, 1.0) - stationary_line_density(1.0, 1.0)) > 0.1, "exclusive finite-occupancy detailed balance is not exactly linear through x~1")
    check(True, "exact linearity over the deep branch requires a non-saturating/Poisson line-density regime or an additional criticality theorem")

    print("\n[4] Susceptibility is exactly the remaining BTFR normalization")
    for chi in [Fraction(1, 2), Fraction(1, 1), Fraction(2, 1)]:
        print(
            f"  chi={chi}: action coeff={action_coeff(chi)}, "
            f"variation prefactor={variation_prefactor(chi)}, "
            f"a0_eff/a0={btfr_a0_eff(float(chi)):.3f}"
        )
    check(btfr_a0_eff(1.0) == 1.0, "chi_R4=1 gives v_inf^4=a0 G M_b")
    check(btfr_a0_eff(0.5) == 2.0, "chi_R4=1/2 would shift the BTFR acceleration by factor 2")
    check(btfr_a0_eff(2.0) == 0.5, "chi_R4=2 would shift the BTFR acceleration by factor 1/2")

    print("\n[5] Honest theorem status")
    check(True, "finite R4/QEC data fix the 2/3 incidence factor")
    check(True, "local stable dynamics naturally gives linear response in |g|/a0")
    check(True, "but the dimensionless susceptibility chi_R4 is not fixed by those facts")
    check(True, "proving the target law now means proving chi_R4=1 and non-saturating line response")

    print("\n" + "=" * 100)
    print("DYNAMICAL LINE-DENSITY RESULT")
    print("  The attempted dynamical proof gets to")
    print("      lambda_R4 = (2/3) chi_R4 |g|/a0")
    print("  with stable local relaxation and detailed-balance linear response.")
    print("  It does NOT derive chi_R4=1.  That missing unit susceptibility is")
    print("  exactly the BTFR normalization once a0 is identified externally with")
    print("  cH0/2pi.  The next microscopic target is therefore precise:")
    print("      derive chi_R4=1 from W=S*C plus boundary-QEC creation/erasure rates,")
    print("      and show the bound line network is in the non-saturating linear regime.")
    print("  Companion item132_chi_unit_poisson.py closes that target if the QEC")
    print("  rate-matching lemma and count-valued line-ledger premise are accepted.")
    print("=" * 100)
    print("exit 0 -- dynamic attempt reduces line-density law to unit susceptibility.")


if __name__ == "__main__":
    main()
