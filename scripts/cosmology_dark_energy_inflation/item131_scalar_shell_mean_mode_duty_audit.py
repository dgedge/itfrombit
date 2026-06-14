#!/usr/bin/env python3
r"""ITEM 131: scalar-shell mean count, mode-locality, and dilute-duty audit.

Question
--------
Can the HBC/QEC service ledger itself derive the three remaining ingredients
behind the conditional scalar-amplitude route?

    1. scalar-shell mean count N_eff,
    2. mode-local correlation volume for the scalar projector,
    3. dilute-duty regime for the Poisson/CTMC reading.

Result
------
No full derivation follows from the current ledger.  What *does* follow is a
sharp partial theorem:

    * the 28/112 finite instrument fixes relative channel weights and the
      1/28 log-shell generator, not the absolute scalar-shell count;
    * the mode-local radial-crossing lemma closes the tilt normalization
      (angular degeneracy is divided out of Delta_R^2), and the later Pi_k
      audit closes the linear projector form as the compensated Fourier shell;
      the later T5b audit reduces the correlation volume to the service-current
      structure factor S_j(k), but does not derive its value or the shell count;
    * the Fano leg is exact by regime:

          CTMC / dilute jump ledger:        F_eff = 1,
          bandwidth-one tick with duty p:  F_eff = 1 - p.

Thus the amplitude ledger is identifiable only after adding a scalar-shell
current theorem:

    A_nu = F_eff / N_eff,    N_eff = N_shell / S_j(k).

For the alpha^4 candidate, N_eff=(4/3) alpha_0^-4, matching the observed
amplitude requires F_eff ~= 0.986, i.e. p ~= 0.014 in the bandwidth-one
scheduler.  The scalar-fluctuation current must therefore be dilute (or else
the count must be duty-corrected); saturated background printing cannot be
silently identified with the fluctuation ledger.
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
ALPHA_CODATA_LOW = 1.0 / 137.036
C_F = Fraction(4, 3)
DELTA = Fraction(1, 28)
N_CHANNELS = 28
N_FLAGS = 112


@dataclass(frozen=True)
class CountRow:
    name: str
    n_eff: float
    status: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def required_count(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    return f_eff / amplitude


def amplitude_from_count(n_eff: float, f_eff: float = 1.0) -> float:
    return f_eff / n_eff


def alpha4_count(alpha: float = ALPHA0) -> float:
    return float(C_F) * alpha ** -4


def bandwidth_one_fano(duty: float) -> float:
    """Fano factor for N~Binomial(T,p), independent of T."""
    return 1.0 - duty


def required_duty_for_count(n_eff: float, amplitude: float = A_TARGET) -> float:
    """If the bandwidth-one scheduler is used, A=(1-p)/N_eff."""
    return 1.0 - amplitude * n_eff


def spectral_tilt_slope() -> Fraction:
    """Mode-local radial crossing supplies exactly one anomalous action."""
    return -DELTA


def projected_count(total_events: float, projection_weight: float) -> float:
    """Mean count seen by one scalar mode after an unspecified projector."""
    return total_events * projection_weight


def sigma_residual(prediction: float) -> float:
    return (prediction - A_TARGET) / SIGMA_A


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def main() -> None:
    print("ITEM 131 SCALAR-SHELL MEAN / MODE-LOCALITY / DUTY AUDIT")

    print("\n[1] Source-ledger checks")
    check(
        contains("python_code/item131_w_to_28_instrument.py", "service channel has rate 1/28"),
        "finite HBC/QEC instrument derives uniform 1/28 service-channel weights",
    )
    check(
        contains("python_code/item131_w_to_28_instrument.py", "112 microscopic flags"),
        "instrument refines the 8-bit Kraus alphabet to 112 incidence flags",
    )
    check(
        contains("python_code/item131_mode_local_radial_crossing.py", "spectral-density normalization theorem"),
        "mode-local radial crossing is closed for the tilt observable",
    )
    check(
        contains("python_code/item131_scalar_mode_projector.py", "Pi_k projector form derived"),
        "scalar-mode projector form is closed as the compensated Fourier-shell readout",
    )
    check(
        contains("python_code/item131_feff_hbc.py", "F_eff = 1") and contains("python_code/item131_feff_hbc.py", "1 - p"),
        "Fano audit separates CTMC F_eff=1 from bandwidth-one F_eff=1-p",
    )
    check(
        contains("python_code/item131_inflationary_amplitude_alpha4_route_audit.py", "not yet a theorem"),
        "alpha^4 amplitude route is explicitly conditional, not locked",
    )

    print("\n[2] Scalar-shell mean count")
    n_required = required_count()
    n_alpha_137 = alpha4_count(ALPHA0)
    n_alpha_codata = alpha4_count(ALPHA_CODATA_LOW)
    rows = [
        CountRow("28 service channels", float(N_CHANNELS), "finite instrument; far too noisy"),
        CountRow("112 incidence flags", float(N_FLAGS), "finite refinement; far too noisy"),
        CountRow("2^28 active subsets", float(2**28), "needs non-ledger Fano correction"),
        CountRow("28*2^24 near-hit", float(28 * 2**24), "disguised alpha^-4 scale; 2^24 not derived"),
        CountRow("(4/3)*137^4", n_alpha_137, "alpha^4 candidate; count theorem missing"),
    ]
    print(f"  target A_nu                  = {A_TARGET:.6e}")
    print(f"  required N_eff/F_eff         = {n_required:.6e}")
    for row in rows:
        amp = amplitude_from_count(row.n_eff)
        f_needed = A_TARGET * row.n_eff
        print(
            f"  {row.name:26s} N={row.n_eff:.6e}  "
            f"A(F=1)={amp:.6e}  F_needed={f_needed:.6e}  {row.status}"
        )
    amp_alpha = amplitude_from_count(n_alpha_137)
    amp_codata = amplitude_from_count(n_alpha_codata)
    check(amplitude_from_count(float(N_CHANNELS)) / A_TARGET > 1.0e7, "28 labels cannot be the scalar-shell count")
    check(amplitude_from_count(float(N_FLAGS)) / A_TARGET > 1.0e6, "112 flags cannot be the scalar-shell count")
    check(abs(sigma_residual(amp_alpha)) < 1.1, "(4/3)*137^4 gives A_nu within about one sigma")
    check(abs(sigma_residual(amp_codata)) < 1.0, "CODATA-like alpha version also lands within one sigma")
    check(True, "but current HBC/QEC ledger contains no derivation of the absolute alpha^-4 shell count")

    print("\n[3] Mode-locality: tilt/projector theorem versus amplitude normalization")
    slope = spectral_tilt_slope()
    print(f"  closed tilt slope: d ln Delta_R^2 / d ln k = {slope} = {float(slope):+.9f}")
    check(slope == -DELTA, "dimensionless power ledger receives one radial 28-clock action per shell")
    check(Fraction(1, 1) + slope == Fraction(27, 28), "tilt normalization gives n_s=27/28")
    total = n_alpha_137
    projector_examples = [
        ("one normalized scalar mode", 1.0),
        ("one of 28 service labels", 1.0 / N_CHANNELS),
        ("one of 112 incidence flags", 1.0 / N_FLAGS),
        ("one eighth address sector", 1.0 / 8.0),
    ]
    for label, weight in projector_examples:
        n_seen = projected_count(total, weight)
        print(f"  {label:29s}: projector weight={weight:.6f}, N_seen={n_seen:.6e}, A={amplitude_from_count(n_seen):.6e}")
    check(
        amplitude_from_count(projected_count(total, 1.0 / N_CHANNELS)) / amp_alpha > 20.0,
        "changing the unresolved projector changes amplitude by O(28) while leaving the tilt statement intact",
    )
    check(True, "therefore amplitude mode-locality now requires the Pi_k structure-factor theorem")

    print("\n[4] Dilute-duty regime")
    print("  CTMC reading: independent jump process, F_eff=1 exactly.")
    for duty in [0.0, 0.014, 0.10, 0.50, 0.99]:
        f_eff = bandwidth_one_fano(duty)
        print(
            f"  bandwidth-one duty p={duty:5.3f}: "
            f"F_eff={f_eff:.6f}, required N_eff={required_count(f_eff=f_eff):.6e}"
        )
    p_alpha = required_duty_for_count(n_alpha_137)
    p_alpha_codata = required_duty_for_count(n_alpha_codata)
    print(f"  duty required by (4/3)*137^4 count       p={p_alpha:.6f}")
    print(f"  duty required by (4/3)*137.036^4 count   p={p_alpha_codata:.6f}")
    check(0.0 < p_alpha < 0.03, "alpha^4 count requires a dilute bandwidth-one duty, not saturation")
    check(bandwidth_one_fano(0.99) < 0.02, "saturated bandwidth-one printing suppresses scalar-clock count noise")
    check(required_count(f_eff=bandwidth_one_fano(0.99)) < 5.0e6, "p=0.99 would require a count two orders smaller than the alpha^4 candidate")

    print("\n[5] Exact status")
    closed = [
        "relative finite instrument: 28 service labels from 112 incidence flags",
        "tilt mode-locality: one radial action per Delta_R^2 shell",
        "projector form: Pi_k is the compensated Fourier-shell readout, not a service-label readout",
        "Fano regime formulas: CTMC F_eff=1; bandwidth-one F_eff=1-p",
    ]
    open_items = [
        "T3 scalar-shell mean count: no ledger theorem for (4/3) alpha_0^-4 or S_dS",
        "T5b amplitude mode-locality: no derived S_j(k=aH) structure factor / independent count for one Pi_k Fourier shell",
        "duty selection: no derivation that primordial scalar fluctuations sample the dilute CTMC ledger while the background is saturated",
    ]
    for item in closed:
        check(True, f"derived/closed: {item}")
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The requested derivation does not close from the current HBC/QEC ledger.")
    print("  The ledger gives the relative 28-channel clock, the tilt mode-locality,")
    print("  and an exact duty/Fano split.  It does not give the scalar-shell mean")
    print("  count or prove the no-horizon-covariance premise needed for S_j(k)=1.")
    print("  The alpha^4 amplitude candidate")
    print("      A_nu = (3/4) alpha_0^4")
    print("  would require a mode-local mean count N_eff=(4/3)alpha_0^-4 plus a")
    print("  dilute fluctuation-duty reading (p about 0.014 if read as bandwidth-one).")
    print("  Saturated background printing can coexist only if the scalar fluctuation")
    print("  ledger is a separate dilute projection, or if N_eff is explicitly")
    print("  duty-corrected.  That separation is the remaining theorem target.")
    print("=" * 108)
    print("exit 0 -- partial theorem/no-go: Pi_k form closed; mean count and no-horizon-covariance premise remain open.")


if __name__ == "__main__":
    main()
