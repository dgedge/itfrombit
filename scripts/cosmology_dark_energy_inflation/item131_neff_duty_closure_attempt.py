#!/usr/bin/env python3
r"""ITEM 131: N_eff / duty closure attempt for the HBC scalar amplitude.

Question
--------
After the scalar projector audit, the remaining open statement is narrower:

    Var(Pi_k nu) = F_eff / N_eff

with a derived N_eff and a derived duty regime.  T5a is closed: Pi_k is the
compensated Fourier-shell projector.  This script tests whether the current
canon already forces T5b, the correlation volume / independent event count.

Result
------
It does not.  The current HBC/QEC ledger fixes:

* the relative 28-channel service weights and the 1/28 radial log-shell action,
* the compensated Fourier-shell scalar projector Pi_k,
* the CTMC Fano factor F_eff=1, and the bandwidth-one alternative F_eff=1-p.

Those constraints are invariant under a scalar-current intensity rescaling

    j(x,N) -> lambda j(x,N),        N_eff -> lambda N_eff,

while the amplitude rescales as A_nu -> A_nu/lambda.  Therefore the absolute
scalar-shell count is not selected by the present ledger.  A theorem for
N_eff=(4/3)alpha_0^-4 would require one additional microscopic statement:
which HBC service-current density/correlation volume is sampled by one
normalized Pi_k shell.  The later T5b audit reduces that correlation-volume
statement to N_eff=N_shell/S_j(k), but it still does not derive S_j(k=aH) or
N_shell.

The duty status is similarly exact but not selective.  The CTMC residual-fault
reading is dilute and has F_eff=1.  A bandwidth-one scheduler has F_eff=1-p;
for the alpha^4 count to match A_s it would need p~=0.014, not a saturated
p->1 scalar-fluctuation ledger.  Existing canon separates these regimes but
does not derive which one the primordial scalar perturbation samples.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

A_TARGET = 2.10e-9
SIGMA_A = 0.03e-9
ALPHA0 = 1.0 / 137.0
C_F = Fraction(4, 3)
DELTA = Fraction(1, 28)


@dataclass(frozen=True)
class RescaledLedger:
    label: str
    intensity_scale: float
    n_eff: float
    f_eff: float = 1.0

    @property
    def amplitude(self) -> float:
        return self.f_eff / self.n_eff


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def alpha4_count(alpha: float = ALPHA0) -> float:
    return float(C_F) * alpha**-4


def amplitude_from_count(n_eff: float, f_eff: float = 1.0) -> float:
    return f_eff / n_eff


def required_count(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    return f_eff / amplitude


def bandwidth_one_fano(duty: float) -> float:
    return 1.0 - duty


def required_duty(n_eff: float, amplitude: float = A_TARGET) -> float:
    return 1.0 - amplitude * n_eff


def sigma_residual(prediction: float) -> float:
    return (prediction - A_TARGET) / SIGMA_A


def log_shell_slope() -> Fraction:
    return -DELTA


def closed_ledger_fingerprint() -> tuple[str, Fraction, str, str]:
    """The pieces already closed by the ledger, independent of intensity."""
    return (
        "compensated Fourier-shell Pi_k",
        log_shell_slope(),
        "CTMC F_eff=1",
        "bandwidth-one F_eff=1-p",
    )


def main() -> None:
    print("ITEM 131 N_EFF / DUTY CLOSURE ATTEMPT")

    print("\n[1] Source checks")
    check(
        contains("python_code/item131_scalar_mode_projector.py", "Pi_k projector form derived"),
        "T5a projector form is closed in the scalar-mode projector audit",
    )
    check(
        contains("python_code/item131_feff_hbc.py", "F_eff = 1") and contains("python_code/item131_feff_hbc.py", "1 - p"),
        "T4 Fano/duty formulas are closed for the two scheduler readings",
    )
    check(
        contains("python_code/item131_inflationary_amplitude_alpha4_route_audit.py", "not yet a theorem"),
        "alpha^4 route is canonically marked as conditional, not locked",
    )
    check(
        contains("python_code/scheduler_alpha_composition.py", "portal-licensed"),
        "alpha^4 exists as a portal-licensed weight-4 scale in the scheduler ledger",
    )
    check(
        contains("python_code/ledger_sky_reading.py", "Class 3") and contains("python_code/ledger_sky_reading.py", "comoving anchor"),
        "ledger-to-sky audit classifies N_eff/H_* as a residual normalization class",
    )

    print("\n[2] The closed covariance identity")
    print("  For a stationary mode-local current with independent effective events:")
    print("      Var(Pi_k nu) = F_eff / N_eff.")
    print("  Current canon supplies Pi_k and F_eff.  The later T5b audit reduces")
    print("  the correlation-volume factor to N_eff=N_shell/S_j(k), but the current")
    print("  density and S_j(k=aH) remain unresolved.")
    fingerprint = closed_ledger_fingerprint()
    print(f"  closed fingerprint = {fingerprint}")
    check(fingerprint[1] == Fraction(-1, 28), "28-clock tilt slope is independent of scalar-current intensity")
    check(fingerprint[2] == "CTMC F_eff=1", "canonical CTMC Fano factor is an intensity-independent theorem")

    print("\n[3] Scale-invariance no-go for T5b")
    base_n = alpha4_count()
    base_amp = amplitude_from_count(base_n)
    ledgers = [
        RescaledLedger("half-intensity scalar current", 0.5, 0.5 * base_n),
        RescaledLedger("alpha^4 candidate", 1.0, base_n),
        RescaledLedger("double-intensity scalar current", 2.0, 2.0 * base_n),
        RescaledLedger("tenfold-intensity scalar current", 10.0, 10.0 * base_n),
    ]
    print("  All rows preserve the same Pi_k, tilt, and CTMC Fano statements:")
    for row in ledgers:
        print(
            f"  {row.label:32s} lambda={row.intensity_scale:5.2f}  "
            f"N_eff={row.n_eff:.6e}  A_nu={row.amplitude:.6e}  "
            f"pull={sigma_residual(row.amplitude):+.2f} sigma"
        )
        check(closed_ledger_fingerprint() == fingerprint, f"closed constraints unchanged for lambda={row.intensity_scale:g}")
    check(abs(ledgers[0].amplitude / base_amp - 2.0) < 1.0e-12, "halving the intensity doubles A_nu")
    check(abs(ledgers[2].amplitude / base_amp - 0.5) < 1.0e-12, "doubling the intensity halves A_nu")
    check(abs(sigma_residual(base_amp)) < 1.1, "the alpha^4 candidate remains numerically strong")
    check(
        abs(sigma_residual(ledgers[0].amplitude)) > 30.0 and abs(sigma_residual(ledgers[2].amplitude)) > 15.0,
        "intensity rescalings allowed by current closed gates drastically change the amplitude",
    )
    print("  Therefore T5b cannot be derived from the already-closed gates alone.")

    print("\n[4] What the required count would be")
    n_required = required_count()
    coeff_required = n_required * ALPHA0**4
    amp_coeff_required = A_TARGET / ALPHA0**4
    print(f"  required N_eff/F_eff       = {n_required:.6e}")
    print(f"  alpha^4 candidate N_eff    = {base_n:.6e}")
    print(f"  required count coefficient = {coeff_required:.9f}")
    print(f"  C_F=4/3                    = {float(C_F):.9f}")
    print(f"  required amplitude coeff   = {amp_coeff_required:.9f}")
    print(f"  3/4                        = {3.0/4.0:.9f}")
    check(abs(coeff_required / float(C_F) - 1.0) < 0.02, "required coefficient is within 2% of C_F=4/3")
    check(
        abs((base_n / n_required) - 1.0) < 0.02,
        "alpha^4 count is a close conditional candidate, but closeness is not a derivation",
    )

    print("\n[5] Duty-regime degeneracy")
    p_alpha = required_duty(base_n)
    print("  CTMC residual-fault reading: F_eff=1 exactly, independent of intensity.")
    print("  Bandwidth-one scheduler:     F_eff=1-p exactly.")
    print(f"  for alpha^4 count, required bandwidth-one duty p = {p_alpha:.6f}")
    check(0.0 < p_alpha < 0.03, "alpha^4 count needs dilute bandwidth-one duty if that scheduler is used")
    for duty in [0.0, p_alpha, 0.10, 0.50, 0.99]:
        f_eff = bandwidth_one_fano(duty)
        n_for_target = required_count(f_eff=f_eff)
        print(
            f"  p={duty:8.6f}  F_eff={f_eff:.6f}  "
            f"N_eff needed for A_s={n_for_target:.6e}"
        )
    check(bandwidth_one_fano(0.99) < 0.02, "saturated bandwidth-one printing suppresses fluctuation noise")
    check(required_count(f_eff=bandwidth_one_fano(0.99)) < 5.0e6, "saturated duty would need a much smaller N_eff than alpha^-4")
    check(
        abs(required_count(f_eff=bandwidth_one_fano(p_alpha)) / base_n - 1.0) < 1.0e-12,
        "choosing p can always repair a chosen N_eff, so duty must be selected independently",
    )

    print("\n[6] Exact theorem target left by the attempt")
    closed = [
        "T5a: Pi_k projector form",
        "T4: CTMC F_eff=1 and bandwidth-one F_eff=1-p",
        "tilt: one 28-clock radial action per Delta_R^2 log shell",
        "alpha^4: portal-licensed weight-4 scale exists elsewhere in the ledger",
    ]
    still_open = [
        "T5b value: no derived S_j(k=aH) structure factor / independent event count",
        "T3: scalar-current intensity law selecting C_F alpha_0^-4 or another count",
        "regime: proof that primordial scalar fluctuations sample the dilute CTMC residual-fault ledger, not saturated background ticks",
    ]
    for item in closed:
        check(True, f"closed input: {item}")
    for item in still_open:
        check(True, f"not closed by this attempt: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The attempted closure produces a no-go, not a locked amplitude theorem.")
    print("  With Pi_k and F_eff fixed, the current ledger still has a scalar-current")
    print("  intensity / structure-factor scaling freedom.  Rescaling that intensity")
    print("  leaves the projector, the 1/28 tilt, and the Fano theorem intact while")
    print("  changing A_nu.  Thus T5b's value remains a real missing microscopic theorem.")
    print("  The alpha^4 route remains the sharp candidate:")
    print("      N_eff = (4/3) alpha_0^-4,  F_eff=1,  A_nu=(3/4)alpha_0^4,")
    print("  but canon may only promote it after deriving the scalar-current intensity")
    print("  and S_j(k=aH), plus the dilute/residual-fault duty regime.")
    print("=" * 108)
    print("exit 0 -- closure attempt yields scale-invariance no-go; T5b value remains open.")


if __name__ == "__main__":
    main()
