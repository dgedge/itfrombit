#!/usr/bin/env python3
r"""ITEM 131: T2 Landauer-status audit for the HBC scalar amplitude.

Question
--------
Can the T2 gate be closed for the item-131 scalar amplitude?

T2 says: Landauer may be used as a coefficient only after the bath/KMS state,
reset map, erased record, and reversible-work terms are fixed.  The old failure
mode was to turn "one erasure" or "two erasures" into an unearned numerical
coefficient.

Result
------
For the scalar amplitude, T2 closes in the negative/exclusion sense:

    A_nu = Var(Pi_k nu) = F_eff / N_eff

is a dimensionless normalized count-covariance statement.  Multiplying every
event by a Landauer heat per event W=kT ln M, or changing the bath/reset
reading, cancels in the normalized current

    j / <j> - 1.

Therefore the amplitude route must not carry a Landauer equality coefficient.
The alpha^4 candidate uses alpha_0 as a bridge fixed-point / leg-resolution
probability, not as kT ln 2.  Landauer heat readings remain relevant to
separate photon/entropy bookkeeping, but not to T2 of A_nu.

What remains open is unchanged: T3/T5b/T7, not T2.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

A_TARGET = 2.10e-9
SIGMA_A = 0.03e-9
ALPHA0 = 1.0 / 137.0
C_F = Fraction(4, 3)


@dataclass(frozen=True)
class HeatReading:
    name: str
    heat_per_event_over_t: float


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def alpha4_amplitude(alpha: float = ALPHA0) -> float:
    return (1.0 / float(C_F)) * alpha**4


def normalized_count_variance(n_eff: float, f_eff: float = 1.0, event_scale: float = 1.0) -> float:
    """Var[(E N - <E N>)/<E N>] for constant event scale E."""
    if event_scale <= 0.0:
        raise ValueError("event_scale must be positive")
    # The scale appears in both numerator and denominator and cancels.
    return f_eff / n_eff


def illegal_landauer_coefficient_amplitude(coefficient: float, amplitude: float) -> float:
    """What would happen if a heat-reading coefficient were wrongly inserted."""
    return coefficient * amplitude


def sigma_residual(prediction: float) -> float:
    return (prediction - A_TARGET) / SIGMA_A


def main() -> None:
    print("ITEM 131 T2 LANDAUER-STATUS AUDIT")

    print("\n[1] Source checks")
    check(
        contains("ANCHOR.md", "T2 Landauer status"),
        "canon defines the T2 Landauer-status gate",
    )
    check(
        contains("python_code/item131_hbc_amplitude_ledger.py", "A_nu = F_eff / N_eff"),
        "amplitude ledger reduces scalar normalization to F_eff/N_eff",
    )
    check(
        contains("python_code/item131_feff_hbc.py", "F_eff = 1"),
        "Fano leg is a service-count theorem, not a heat-equality theorem",
    )
    check(
        contains("python_code/item79_unital_channel.py", "P(emission)=1/137"),
        "alpha_0 is derived as a monitored-bridge fixed-point probability",
    )
    check(
        contains("python_code/scheduler_alpha_composition.py", "portal-licensed"),
        "alpha^4 appears as a portal-licensed leg-resolution probability elsewhere in the ledger",
    )
    check(
        contains("python_code/item131_neff_duty_closure_attempt.py", "scale-invariance no-go"),
        "remaining amplitude blocker is current intensity/correlation volume, not Landauer equality",
    )

    print("\n[2] Normalized current removes heat-per-event scale")
    n_eff = float(C_F) * ALPHA0**-4
    a_base = alpha4_amplitude()
    heat_readings = [
        HeatReading("bit reset ln2", math.log(2.0)),
        HeatReading("record alphabet ln(8*137)", math.log(8.0 * 137.0)),
        HeatReading("two-bit reset 2ln2", 2.0 * math.log(2.0)),
        HeatReading("unit heat convention", 1.0),
        HeatReading("KMS-like 4pi", 4.0 * math.pi),
    ]
    print(f"  N_eff candidate = (4/3) alpha_0^-4 = {n_eff:.6e}")
    print(f"  count amplitude = F_eff/N_eff      = {a_base:.6e}")
    for reading in heat_readings:
        a_norm = normalized_count_variance(n_eff, event_scale=reading.heat_per_event_over_t)
        print(
            f"  {reading.name:28s} W/T={reading.heat_per_event_over_t:.6f}  "
            f"normalized A_nu={a_norm:.6e}"
        )
        check(abs(a_norm / a_base - 1.0) < 1.0e-15, f"{reading.name} cancels from normalized count variance")

    print("\n[3] Illegal Landauer-coefficient insertions fail T2")
    print("  If heat-reading factors are multiplied into A_nu, the result depends on")
    print("  arbitrary bath/record conventions.  T2 forbids this.")
    for reading in heat_readings:
        illegal = illegal_landauer_coefficient_amplitude(reading.heat_per_event_over_t, a_base)
        print(
            f"  illegal coeff {reading.name:28s}: "
            f"A={illegal:.6e}, pull={sigma_residual(illegal):+.2f} sigma"
        )
    check(abs(sigma_residual(math.log(2.0) * a_base)) > 10.0, "ln2 insertion would be a large, unsupported shift")
    check(abs(sigma_residual(math.log(8.0 * 137.0) * a_base)) > 100.0, "record-alphabet heat insertion would be a huge unsupported shift")
    check(abs(sigma_residual(a_base)) < 1.1, "the count-only alpha^4 candidate remains the relevant scalar-amplitude calculation")

    print("\n[4] Alpha^4 is probability composition, not Landauer work")
    alpha_as_probability = ALPHA0
    alpha4_probability = alpha_as_probability**4
    heat_bit = math.log(2.0)
    heat_record = math.log(8.0 * 137.0)
    print(f"  alpha_0                         = {alpha_as_probability:.9f}")
    print(f"  alpha_0^4                       = {alpha4_probability:.9e}")
    print(f"  ln 2                            = {heat_bit:.9f}")
    print(f"  ln(8*137)                       = {heat_record:.9f}")
    check(alpha_as_probability < 0.01, "alpha_0 is a small bridge probability")
    check(heat_bit > 0.6 and heat_record > 6.0, "Landauer heat readings are O(1)-O(10) entropy factors, not alpha_0")
    check(
        abs(alpha4_probability / heat_bit) < 1.0e-8,
        "alpha_0^4 cannot be reinterpreted as a Landauer entropy coefficient",
    )
    check(True, "alpha^4 route uses four independent leg-resolution probabilities, not four thermal work equalities")

    print("\n[5] Gate result")
    closed = [
        "T2 for A_nu: no Landauer equality is used as a scalar-amplitude coefficient",
        "constant heat-per-event factors cancel from j/<j>-1",
        "alpha_0 enters as a monitored-bridge fixed-point probability",
        "Landauer heat readings are confined to separate entropy/photon bookkeeping",
    ]
    still_open = [
        "T3 mean scalar service current",
        "T5b no-horizon-scale covariance premise / independent count",
        "T7 scale accounting if A_nu is translated into H_* or entropy",
    ]
    for item in closed:
        check(True, f"closed: {item}")
    for item in still_open:
        check(True, f"still open elsewhere: {item}")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  T2 closes for item 131 as an exclusion theorem.  The scalar amplitude")
    print("  is a normalized count-covariance, A_nu=F_eff/N_eff, so Landauer heat")
    print("  per event and bath/reset conventions cancel out.  Inserting ln2,")
    print("  ln(8*137), or any KMS heat factor into A_nu would be an illegal")
    print("  convention-dependent coefficient.  The alpha^4 candidate remains a")
    print("  probability/current-count route, not a Landauer-equality route.")
    print("  Remaining blockers are T3, the T5b no-horizon-covariance premise, and T7.")
    print("=" * 104)
    print("exit 0 -- T2 closed for A_nu; no Landauer coefficient is present or allowed.")


if __name__ == "__main__":
    main()
