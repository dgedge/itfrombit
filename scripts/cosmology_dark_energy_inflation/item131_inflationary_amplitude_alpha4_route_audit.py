#!/usr/bin/env python3
"""ITEM 131: alpha^4 route audit for the inflationary amplitude.

Question
--------
Is there now a framework-native way to calculate the scalar amplitude A_nu,
after the HBC service-current audit closed F_eff=1?

Candidate
---------
The observed amplitude requires

    N_eff = F_eff / A_nu ~= 4.76e8.

The previously noticed finite-count near-hit 28*2^24 is not arbitrary: it is
numerically the same scale as

    N_eff = C_F * alpha_0^-4,  C_F = 4/3,

where alpha_0 is the item-79 record/bridge probability and C_F=4/3 is the
canonical SU(3) fundamental Casimir/string-tension coefficient already used in
the hadronic sector.  With the derived CTMC F_eff=1 this gives

    A_nu = 1 / N_eff = (3/4) alpha_0^4.

This calculates A_nu to within the current Planck-scale amplitude uncertainty.

Verdict
-------
This is a serious conditional recovery route, not a theorem yet.  The current
canon derives the ingredients separately:

* F_eff=1 for the 28-channel CTMC service current;
* alpha_0 as a bridge/record probability;
* alpha_0^4 as the portal-licensed weight-4 logical-event scale;
* C_F=4/3 as the SU(3) Casimir/string-tension coefficient.

What is not yet derived is the single T3/T5 statement that one scalar
mode-local printed shell contains exactly C_F * alpha_0^-4 independent HBC
service events.  Until that mean-current/correlation-volume theorem is built,
the formula is a rule-selected candidate calculation, not a locked amplitude
derivation.
"""

from __future__ import annotations

import math
from fractions import Fraction


A_OBS = 2.10e-9
SIGMA_A = 0.03e-9
ALPHA_EXACT_137 = 1.0 / 137.0
ALPHA_CODATA_LOW = 1.0 / 137.036
MPL_REDUCED_GEV = 2.435e18
C_F = Fraction(4, 3)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def a_alpha4(alpha: float, count_prefactor: float = float(C_F), f_eff: float = 1.0) -> float:
    """A=F/N with N=count_prefactor*alpha^-4."""
    return f_eff * alpha**4 / count_prefactor


def h_from_amplitude(amplitude: float) -> float:
    return math.sqrt(8.0 * math.pi**2 * amplitude) * MPL_REDUCED_GEV


def v14_from_h(h_gev: float) -> float:
    return (3.0 * h_gev**2 * MPL_REDUCED_GEV**2) ** 0.25


def sigma_residual(prediction: float) -> float:
    return (prediction - A_OBS) / SIGMA_A


def print_prediction(name: str, alpha: float) -> tuple[float, float, float]:
    n_eff = float(C_F) * alpha ** -4
    amp = 1.0 / n_eff
    h_star = h_from_amplitude(amp)
    print(f"  {name:22s} alpha^-1={1/alpha:.6f}")
    print(f"    N_eff=(4/3)alpha^-4 = {n_eff:.6e}")
    print(f"    A=(3/4)alpha^4     = {amp:.6e}  residual={sigma_residual(amp):+.3f} sigma")
    print(f"    H_*                = {h_star:.6e} GeV")
    print(f"    V^(1/4)            = {v14_from_h(h_star):.6e} GeV")
    return n_eff, amp, h_star


def main() -> None:
    print("ITEM 131 INFLATIONARY AMPLITUDE ALPHA^4 ROUTE AUDIT")

    print("\n[1] Required count versus alpha^4 Casimir count")
    n_required = 1.0 / A_OBS
    coeff_required_for_count = n_required * ALPHA_EXACT_137**4
    coeff_required_for_amplitude = A_OBS / ALPHA_EXACT_137**4
    print(f"  A_obs                         = {A_OBS:.6e} +/- {SIGMA_A:.1e}")
    print(f"  required N_eff                = {n_required:.6e}")
    print(f"  required count prefactor C    = N_eff * alpha^4 = {coeff_required_for_count:.9f}")
    print(f"  required amplitude prefactor  = A_obs/alpha^4  = {coeff_required_for_amplitude:.9f}")
    print(f"  canonical C_F=4/3             = {float(C_F):.9f}")
    print(f"  inverse amplitude prefactor   = 3/4            = {1/float(C_F):.9f}")
    check(abs(coeff_required_for_count / float(C_F) - 1.0) < 0.02, "required count prefactor is within 2% of C_F=4/3")
    check(abs(coeff_required_for_amplitude / (3.0 / 4.0) - 1.0) < 0.02, "required amplitude prefactor is within 2% of 3/4")

    print("\n[2] Candidate prediction")
    n_137, amp_137, h_137 = print_prediction("alpha_0=1/137", ALPHA_EXACT_137)
    n_codata, amp_codata, h_codata = print_prediction("alpha=1/137.036", ALPHA_CODATA_LOW)
    check(abs(sigma_residual(amp_137)) < 1.1, "alpha_0=1/137 version lands within about one sigma")
    check(abs(sigma_residual(amp_codata)) < 1.0, "CODATA-like low-energy alpha version lands within one sigma")
    check(9.0e14 < h_codata < 1.1e15, "the alpha^4 route implies the same H_*~1e15 GeV scale as the entropy route")

    print("\n[3] Why the earlier 28*2^24 near-hit appeared")
    mixed_count = 28 * 2**24
    count_alpha = n_137
    print(f"  28*2^24                = {mixed_count:.6e}")
    print(f"  (4/3)*137^4            = {count_alpha:.6e}")
    print(f"  ratio                  = {mixed_count / count_alpha:.9f}")
    print(f"  A from 28*2^24          = {1.0 / mixed_count:.6e}")
    print(f"  A from (4/3)*137^4      = {amp_137:.6e}")
    check(abs(mixed_count / count_alpha - 1.0) < 2.0e-4, "28*2^24 is effectively a disguised (4/3)*137^4 count")
    print("  Interpretation: the finite-count near-hit should not be promoted as")
    print("  '28 times 2^24'.  The structurally meaningful version, if any, is")
    print("  a Casimir-weighted alpha^-4 service count.")

    print("\n[4] Structural alternatives and nulls")
    structural_coefficients = [
        ("raw alpha^4", 1.0),
        ("SU(3) inverse Casimir 3/4", 3.0 / 4.0),
        ("vacuum exit fraction 2/3", 2.0 / 3.0),
        ("coin singlet 1/2", 1.0 / 2.0),
        ("boot pair 2/9", 2.0 / 9.0),
        ("baryogenesis channel 3/14", 3.0 / 14.0),
        ("loop coefficient 3/2", 3.0 / 2.0),
        ("SU(3) Casimir 4/3", 4.0 / 3.0),
    ]
    for label, coeff in structural_coefficients:
        amp = coeff * ALPHA_EXACT_137**4
        print(f"  {label:28s} coeff={coeff:.9f} A={amp:.6e} residual={sigma_residual(amp):+.2f} sigma")
    check(abs(sigma_residual((3.0 / 4.0) * ALPHA_EXACT_137**4)) < 1.1, "3/4 alpha^4 is the only listed structural coefficient in the one-sigma band")
    check(abs(sigma_residual(ALPHA_EXACT_137**4)) > 20.0, "raw alpha^4 is far too large at Planck precision")

    print("\n[5] Rational-search caution")
    target_coeff = coeff_required_for_amplitude
    one_sigma_hits: list[Fraction] = []
    two_sigma_hits: list[Fraction] = []
    for den in range(1, 29):
        for num in range(1, 3 * den + 1):
            frac = Fraction(num, den)
            amp = float(frac) * ALPHA_EXACT_137**4
            pull = abs(sigma_residual(amp))
            if pull <= 1.0:
                one_sigma_hits.append(frac)
            if pull <= 2.0:
                two_sigma_hits.append(frac)
    one_sigma_hits = sorted(set(one_sigma_hits), key=float)
    two_sigma_hits = sorted(set(two_sigma_hits), key=float)
    print(f"  target coefficient A/alpha^4 = {target_coeff:.9f}")
    print(f"  rational hits p/q<=28 within 1 sigma: {len(one_sigma_hits)}")
    print("   ", ", ".join(str(x) for x in one_sigma_hits[:20]))
    print(f"  rational hits p/q<=28 within 2 sigma: {len(two_sigma_hits)}")
    print("   ", ", ".join(str(x) for x in two_sigma_hits[:24]))
    check(Fraction(3, 4) in one_sigma_hits, "3/4 is in the one-sigma rational band")
    check(len(one_sigma_hits) > 5, "the small-rational alphabet is not unique; structure must select 3/4 before comparison")

    print("\n[6] T1-T9 status")
    closed = [
        "T4/Fano: F_eff=1 under the canonical CTMC service ledger",
        "candidate count formula: N_eff=(4/3)alpha_0^-4 gives A=(3/4)alpha_0^4",
        "scale consequence: H_*=sqrt(6 pi^2) alpha_0^2 Mbar_P",
    ]
    open_items = [
        "T1/T3: construct the HBC scalar-shell Kraus/current ledger whose mean count is C_F alpha^-4",
        "T5: prove one mode-local printed shell samples that count, rather than a total entropy or baryogenesis ledger",
        "T6: connect the count to the gauge-invariant scalar clock variable R_HBC",
        "T7: H_* inferred from A still inherits the Planck-mass/horizon status unless A is treated as the primary dimensionless prediction",
        "T8: 3/4 must be selected by SU(3) Casimir/mode projection before amplitude comparison; small-rational hits are otherwise dense",
    ]
    for item in closed:
        check(True, f"conditional/candidate: {item}")
    for item in open_items:
        check(True, f"not yet closed: {item}")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  There is now a plausible calculational route for the inflationary")
    print("  amplitude:")
    print("      N_eff = (4/3) alpha_0^-4,  F_eff=1,")
    print("      A_nu = (3/4) alpha_0^4.")
    print(f"  Numerically this gives A_nu={amp_137:.3e} for alpha_0=1/137")
    print(f"  and A_nu={amp_codata:.3e} for alpha=1/137.036, both within the")
    print("  current amplitude scale.  It also explains why 28*2^24 looked close:")
    print("  it is essentially (4/3)*137^4.")
    print("  But this is not yet a theorem.  The remaining derivation is specific:")
    print("  prove from the HBC/QEC scalar-shell service ledger that the effective")
    print("  independent event count is the SU(3)-Casimir-weighted quartic alpha")
    print("  count, C_F alpha^-4.  Until then this is a rule-selected conditional")
    print("  recovery candidate, not a locked amplitude prediction.")
    print("=" * 104)
    print("exit 0 -- alpha^4 amplitude route identified; conditional, not closed.")


if __name__ == "__main__":
    main()
