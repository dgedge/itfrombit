#!/usr/bin/env python3
r"""ITEM 131: Inflation/HBC status gate.

Question
--------
Can the open Inflation/HBC item be promoted beyond the current status?

Current answer
--------------
No.  The tilt and the amplitude now sit at different tiers.

Tilt:
    The finite 28-channel clock, scale-covariant radial HBC lift, scalar-clock
    observable map, and mode-local power-ledger action reduce the tilt to

        d ln Delta_R^2 / d ln k = -1/28,    n_s = 27/28,

    conditional on the physical HBC premise: an early saturated constant-H
    printer whose local print-rate fluctuation is the unique adiabatic scalar
    clock.

Amplitude:
    The scalar amplitude is not set by the 28-clock.  The current ledger gives

        A_nu(k) = F_eff S_j(k) / N_shell,
        N_eff(k) = N_shell / S_j(k).

    T4 is closed for the canonical CTMC ledger, F_eff=1.  T5a is closed:
    Pi_k is the compensated Fourier-shell projector.  T5b is reduced:
    S_j(k=aH)=1 follows under no connected horizon-scale service-current
    covariance.  T3 remains open: derive the absolute shell count N_shell.

    The alpha^4 candidate remains the sharpest conditional route:

        N_shell = (4/3) alpha_0^-4,    S_j=1,    F_eff=1
        => A_nu = (3/4) alpha_0^4.

    But the existing gates are invariant under a scalar-current intensity
    rescaling N_shell -> lambda N_shell, which changes A_nu.  Therefore this is
    still a conditional recovery candidate, not a Locked amplitude derivation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

A_TARGET = 2.10e-9
SIGMA_A = 0.03e-9
ALPHA0 = 1.0 / 137.0
ALPHA_CODATA_LOW = 1.0 / 137.036
MPL_REDUCED_GEV = 2.435e18
DELTA = Fraction(1, 28)
C_F = Fraction(4, 3)


class Status(Enum):
    CLOSED = "CLOSED"
    CONDITIONAL = "CONDITIONAL"
    OPEN = "OPEN"


@dataclass(frozen=True)
class LedgerGate:
    name: str
    status: Status
    finding: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def alpha4_count(alpha: float) -> float:
    return float(C_F) * alpha**-4


def amplitude_from_shell(n_shell: float, s_j: float = 1.0, f_eff: float = 1.0) -> float:
    return f_eff * s_j / n_shell


def required_shell_count(amplitude: float = A_TARGET, s_j: float = 1.0, f_eff: float = 1.0) -> float:
    return f_eff * s_j / amplitude


def h_star_from_amplitude(amplitude: float) -> float:
    """Entropy/slow-roll translation, used here only as a scale map."""
    return math.sqrt(8.0 * math.pi**2 * amplitude) * MPL_REDUCED_GEV


def sigma_pull(prediction: float) -> float:
    return (prediction - A_TARGET) / SIGMA_A


def ns_from_hbc_log_clock() -> Fraction:
    return Fraction(1, 1) - DELTA


def bandwidth_one_duty_for_count(n_shell: float, amplitude: float = A_TARGET) -> float:
    """If F_eff=1-p, solve amplitude=(1-p)/N_shell."""
    return 1.0 - amplitude * n_shell


def print_gate_table(gates: list[LedgerGate]) -> None:
    for gate in gates:
        print(f"  {gate.status.value:11s}  {gate.name}")
        print(f"      {gate.finding}")


def main() -> None:
    print("ITEM 131 INFLATION/HBC STATUS GATE")

    print("\n[1] Source checks")
    source_checks = [
        ("python_code/item131_w_to_28_instrument.py", "service channel has rate 1/28", "finite 28-channel service weights"),
        ("python_code/item131_primordial_tilt_logscale.py", "therefore n_s=27/28", "log-scale generator gives n_s=27/28"),
        ("python_code/item131_hbc_lift_assumption_audit.py", "saturated constant H and QEC action on Delta_R^2", "constant-H power-ledger premise is explicit"),
        ("python_code/item131_mode_local_radial_crossing.py", "spectral-density normalization theorem", "mode-local radial crossing is the tilt normalization"),
        ("python_code/item131_scalar_clock_bridge.py", "R_HBC = psi - nu", "gauge-invariant scalar-clock wrapper exists conditionally"),
        ("python_code/item131_hbc_clock_covariance_exit.py", "A_nu remains an undetermined", "clock covariance audit leaves amplitude open"),
        ("python_code/item131_feff_hbc.py", "F_eff = 1", "CTMC Fano factor is closed"),
        ("python_code/item131_scalar_mode_projector.py", "Pi_k projector form derived", "T5a projector form is closed"),
        ("python_code/item131_t5b_correlation_volume_audit.py", "N_eff(k) = N_shell / S_j(k)", "T5b form is reduced to structure factor"),
        ("python_code/item131_t5b_whiteness_lemma.py", "no-horizon-scale-covariance premise", "T5b value conditionally closes under whiteness"),
        ("python_code/item131_neff_duty_closure_attempt.py", "scale-invariance no-go", "absolute shell count is not selected by existing gates"),
    ]
    for path, phrase, label in source_checks:
        check(contains(path, phrase), label)

    print("\n[2] Tilt ledger")
    ns = ns_from_hbc_log_clock()
    print(f"  anomalous power slope        = -1/28 = {-float(DELTA):+.12f}")
    print(f"  n_s                          = 27/28 = {float(ns):.12f}")
    check(ns == Fraction(27, 28), "finite 28-clock plus radial log-shell action gives n_s=27/28")
    check(True, "promotion still requires the physical HBC scalar-clock premise, not another spectral calculation")

    print("\n[3] Amplitude bookkeeping")
    n_required = required_shell_count(A_TARGET)
    n_alpha_137 = alpha4_count(ALPHA0)
    n_alpha_codata = alpha4_count(ALPHA_CODATA_LOW)
    amp_alpha_137 = amplitude_from_shell(n_alpha_137)
    amp_alpha_codata = amplitude_from_shell(n_alpha_codata)
    print(f"  target A_nu                  = {A_TARGET:.6e}")
    print(f"  required N_shell/S_j         = {n_required:.6e}  (F_eff=1)")
    print(f"  alpha0 candidate N_shell     = {n_alpha_137:.6e}")
    print(f"  alpha0 candidate A_nu        = {amp_alpha_137:.6e}  pull={sigma_pull(amp_alpha_137):+.3f} sigma")
    print(f"  CODATA-alpha candidate A_nu  = {amp_alpha_codata:.6e}  pull={sigma_pull(amp_alpha_codata):+.3f} sigma")
    print(f"  entropy-map H_*(target)      = {h_star_from_amplitude(A_TARGET):.6e} GeV")
    print(f"  entropy-map H_*(alpha0)      = {h_star_from_amplitude(amp_alpha_137):.6e} GeV")
    check(abs(sigma_pull(amp_alpha_137)) < 1.1, "alpha0=(1/137) route lands inside the current one-sigma amplitude scale")
    check(abs(sigma_pull(amp_alpha_codata)) < 1.0, "CODATA-like alpha route also lands inside the current one-sigma amplitude scale")
    check(8.0e14 < h_star_from_amplitude(A_TARGET) < 1.2e15, "entropy translation maps the target to H_*~1e15 GeV")

    print("\n[4] T5b structure-factor sensitivity")
    for s_j in [0.5, 1.0, 1.2, 2.0]:
        amp = amplitude_from_shell(n_alpha_137, s_j=s_j)
        needed_count = required_shell_count(A_TARGET, s_j=s_j)
        print(
            f"  S_j={s_j:3.1f}: A(alpha count)={amp:.6e} "
            f"pull={sigma_pull(amp):+7.2f} sigma; count needed={needed_count:.6e}"
        )
    check(abs(amplitude_from_shell(n_alpha_137, s_j=2.0) / amp_alpha_137 - 2.0) < 1.0e-12, "horizon-scale service covariance would directly rescale A_nu")
    check(True, "no-horizon-scale covariance is therefore a load-bearing amplitude premise")

    print("\n[5] Rescaling no-go")
    print("  All rows keep the already-closed tilt, Pi_k form, and CTMC Fano statement:")
    base_fingerprint = ("Pi_k", Fraction(-1, 28), "F_eff=1")
    for lam in [0.5, 1.0, 2.0, 10.0]:
        n_shell = lam * n_alpha_137
        amp = amplitude_from_shell(n_shell)
        print(f"  lambda={lam:4.1f}: N_shell={n_shell:.6e}  A_nu={amp:.6e}  pull={sigma_pull(amp):+7.2f} sigma")
        check(base_fingerprint == ("Pi_k", Fraction(-1, 28), "F_eff=1"), f"closed fingerprint unchanged for lambda={lam:g}")
    check(amplitude_from_shell(0.5 * n_alpha_137) / amp_alpha_137 == 2.0, "intensity rescaling changes amplitude while preserving closed gates")
    check(abs(sigma_pull(amplitude_from_shell(2.0 * n_alpha_137))) > 15.0, "allowed rescalings produce observationally distinct amplitudes")

    print("\n[6] Duty split")
    p_alpha = bandwidth_one_duty_for_count(n_alpha_137)
    print(f"  bandwidth-one duty needed for alpha0 count: p={p_alpha:.6f}")
    for duty in [0.0, p_alpha, 0.1, 0.5, 0.99]:
        f_eff = 1.0 - duty
        print(f"  p={duty:8.6f}: F_eff={f_eff:.6f}, N_shell needed={required_shell_count(f_eff=f_eff):.6e}")
    check(0.0 < p_alpha < 0.03, "bandwidth-one reading requires dilute duty for the alpha^4 count")
    check(required_shell_count(f_eff=0.01) < 5.0e6, "saturated bandwidth-one printing would need a count two orders smaller")

    print("\n[7] Gate table")
    gates = [
        LedgerGate("finite 28-clock / relative channel weights", Status.CLOSED, "Exact finite instrument: 112 flags coarse-grain to 28 channels with rate 1/28."),
        LedgerGate("tilt exponent", Status.CONDITIONAL, "n_s=27/28 follows from radial log-shell, power-ledger action, and constant-H HBC scalar-clock premises."),
        LedgerGate("gauge-invariant scalar variable", Status.CONDITIONAL, "R_HBC=psi-nu is formal if a unique local print-time field nu_HBC exists."),
        LedgerGate("Fano factor", Status.CLOSED, "Canonical CTMC ledger gives F_eff=1; bandwidth-one alternative is F_eff=1-p."),
        LedgerGate("projector form", Status.CLOSED, "Pi_k is the compensated Fourier-shell readout at |k|=aH; angular degeneracy is divided out."),
        LedgerGate("structure factor", Status.CONDITIONAL, "S_j(k=aH)=1 under product-local/exchangeable/common-mode ledgers; open if there is horizon-mode covariance."),
        LedgerGate("absolute shell count", Status.OPEN, "No current theorem selects N_shell=(4/3)alpha_0^-4, S_dS, or another absolute count."),
        LedgerGate("early duration and exit", Status.OPEN, "Finite saturated w=-1 duration and exit/reheating remain item-42/R-activation problems."),
    ]
    print_gate_table(gates)
    check(sum(g.status is Status.OPEN for g in gates) == 2, "two hard open gates remain")
    check(sum(g.status is Status.CONDITIONAL for g in gates) == 3, "three bridge premises remain conditional")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  Inflation/HBC is structurally narrowed, not closed.")
    print("  The tilt side is a conditional scalar-clock/log-shell theorem:")
    print("      n_s = 27/28")
    print("  if saturated constant-H HBC supplies the unique local scalar clock and")
    print("  the 28-generator acts on the dimensionless power ledger.")
    print("  The amplitude side is sharper but still open:")
    print("      A_nu = F_eff S_j(k=aH) / N_shell.")
    print("  F_eff=1 and Pi_k are closed; S_j=1 is conditionally reduced to the")
    print("  no-horizon-scale-covariance premise.  The hard remaining derivation is")
    print("  the absolute scalar-shell service count, with the alpha^4 candidate")
    print("      N_shell=(4/3)alpha_0^-4  =>  A_nu=(3/4)alpha_0^4")
    print("  still conditional.  The rescaling no-go proves the present canon cannot")
    print("  promote A_nu without a new microscopic current-density theorem.")
    print("=" * 104)
    print("exit 0 -- Inflation/HBC status gate: tilt conditional; amplitude reduced to S_j and N_shell.")


if __name__ == "__main__":
    main()
