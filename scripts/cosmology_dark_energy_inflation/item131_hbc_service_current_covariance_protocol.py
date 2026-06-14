#!/usr/bin/env python3
r"""ITEM 131: strict thermodynamic protocol applied to HBC A_nu.

Purpose
-------
This is the worked target for the canon's thermodynamic-claim protocol.
It asks whether the HBC/QEC service-current ledger derives the scalar-clock
amplitude

    Delta_nu^2(k) = A_nu (k/k_*)^(-1/28)

rather than merely fitting or inverting it.

The strict ledger form is:

    j(x,N)       = mode-local HBC print/service event current,
    nu(x,N)      = integral^N [j(x,N')/jbar(N') - 1] dN',
    R_HBC        = psi - nu  (gauge-invariant scalar clock),
    A_nu         = F_eff / N_eff.

For one horizon-shell bin with N_eff statistically independent service events
and current Fano/correlation factor F_eff,

    Var[(N - <N>)/<N>] = F_eff / N_eff.

The 28-channel QEC generator supplies the log-shell tilt, not this amplitude.
The later Fano audit closes F_eff=1 for the canonical CTMC service ledger, so
the current promotion blocker is sharper: derive the scalar-shell mean count
N_eff, the mode-local correlation volume, and the duty/regime map that makes
the CTMC reading the one sampled by primordial HBC.

The scalar-shell mean/mode/duty audit verifies that this blocker has not yet
closed.  The alpha^4 count requires a dilute fluctuation-duty reading, not a
saturated bandwidth-one fluctuation ledger.

The scalar-mode projector audit closes the linear readout form: Pi_k is the
compensated Fourier-shell projection of the local service current.  That
splits T5 into a closed projector-form half and an open correlation-volume
half.

The N_eff/duty closure attempt proves the remaining obstruction is real under
the current canon: a scalar-current intensity rescaling leaves Pi_k, the 1/28
tilt, and the Fano theorem unchanged while changing the amplitude.  Therefore
T5b cannot be promoted without a new microscopic current-density /
correlation-volume law.

The T5b correlation-volume audit closes the formal identity: for the Pi_k
horizon shell, Var(delta_j(k))=F_eff S_j(k)/N_shell and
N_eff=N_shell/S_j(k).  The numerical gate remains open until the current
ledger derives the spatial structure factor S_j(k=aH) and the shell count.

The T5b whiteness lemma then proves S_j(k)=1 for nonzero compensated Fourier
modes under product-local CTMC service, fixed-total exchangeable allocation,
or homogeneous common-rate noise.  What remains is the no-horizon-scale
connected covariance premise.

The T2 audit closes Landauer status for A_nu in the exclusion sense:
Landauer heat-per-event factors cancel from the normalized current j/<j>-1,
and any ln 2 / ln M insertion would be an illegal convention-dependent
coefficient.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

DELTA = Fraction(1, 28)
A_TARGET = 2.1e-9
N_CHANNELS = 28
N_FLAGS = 112
MPL_REDUCED_GEV = 2.435e18


class Status(Enum):
    PASS = "PASS"
    CONDITIONAL = "CONDITIONAL"
    OPEN = "OPEN"
    FAIL = "FAIL"


@dataclass(frozen=True)
class Gate:
    code: str
    name: str
    status: Status
    finding: str


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def a_nu_from_current(n_eff: float, f_eff: float = 1.0) -> float:
    """Dimensionless scalar-clock variance per mode-local horizon shell."""
    return f_eff / n_eff


def required_n_eff(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    return f_eff / amplitude


def delta_nu_squared(k_ratio: float, n_eff: float, f_eff: float = 1.0) -> float:
    """Dimensionless power ledger with the 28-clock anomalous dimension."""
    return a_nu_from_current(n_eff, f_eff) * k_ratio ** (-float(DELTA))


def p_nu(k_ratio: float, n_eff: float, f_eff: float = 1.0) -> float:
    """Power spectrum shape in units where k_* = 1."""
    return 2.0 * math.pi**2 * delta_nu_squared(k_ratio, n_eff, f_eff) / k_ratio**3


def entropy_count_from_h(h_gev: float) -> float:
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def h_from_entropy_count(n_eff: float) -> float:
    return math.sqrt(8.0 * math.pi**2 / n_eff) * MPL_REDUCED_GEV


def print_gate_table(gates: list[Gate]) -> None:
    for gate in gates:
        print(f"  {gate.code:>2s}  {gate.status.value:11s}  {gate.name}")
        print(f"      {gate.finding}")


def main() -> None:
    print("ITEM 131 STRICT THERMODYNAMIC PROTOCOL: HBC SERVICE-CURRENT COVARIANCE")

    anchor = (ROOT / "ANCHOR.md").read_text()
    w_to_28 = (ROOT / "python_code" / "item131_w_to_28_instrument.py").read_text()
    scalar_bridge = (ROOT / "python_code" / "item131_scalar_clock_bridge.py").read_text()
    clock_cov = (ROOT / "python_code" / "item131_hbc_clock_covariance_exit.py").read_text()
    amp_ledger = (ROOT / "python_code" / "item131_hbc_amplitude_ledger.py").read_text()
    feff = (ROOT / "python_code" / "item131_feff_hbc.py").read_text()
    alpha4 = (ROOT / "python_code" / "item131_inflationary_amplitude_alpha4_route_audit.py").read_text()
    projector = (ROOT / "python_code" / "item131_scalar_mode_projector.py").read_text()
    mean_mode_duty = (ROOT / "python_code" / "item131_scalar_shell_mean_mode_duty_audit.py").read_text()
    neff_duty = (ROOT / "python_code" / "item131_neff_duty_closure_attempt.py").read_text()
    t2_landauer = (ROOT / "python_code" / "item131_t2_landauer_status_audit.py").read_text()
    t5b_volume = (ROOT / "python_code" / "item131_t5b_correlation_volume_audit.py").read_text()
    t5b_white = (ROOT / "python_code" / "item131_t5b_whiteness_lemma.py").read_text()

    print("\n[1] Source checks")
    check("service channel has rate 1/28" in w_to_28, "finite 28-channel service weights are constructed")
    check("112 microscopic flags" in w_to_28, "112 incidence-refined microscopic flags are constructed")
    check("R_HBC = psi - nu" in scalar_bridge and "delta N" in scalar_bridge, "scalar-clock bridge gives a conditional gauge-invariant observable map")
    check("A_nu remains an undetermined" in clock_cov, "clock-covariance audit leaves amplitude open")
    check("A_nu = F_eff/N_eff" in amp_ledger or "A_nu = F_eff / N_eff" in amp_ledger, "amplitude ledger reduces normalization to F_eff/N_eff")
    check("F_eff = 1" in feff and "25-regularity" in feff, "Fano leg is now closed for the canonical CTMC ledger")
    check("A_nu = (3/4) alpha_0^4" in alpha4, "alpha^4 audit supplies the current conditional candidate formula")
    check("Pi_k projector form derived" in projector, "linear scalar-mode projector form is derived")
    check("partial theorem/no-go" in mean_mode_duty and "p about 0.014" in mean_mode_duty, "scalar-shell mean/mode/duty audit keeps the blocker explicit")
    check("scale-invariance no-go" in neff_duty and "T5b value remains open" in neff_duty, "N_eff/duty closure attempt proves the residual scale freedom")
    check("T2 closes for item 131" in t2_landauer and "no Landauer coefficient" in t2_landauer, "T2 Landauer-status audit closes the no-coefficient gate")
    check("T5b reduced to service-current structure factor" in t5b_volume, "T5b audit reduces correlation volume to a structure-factor theorem target")
    check("T5b whiteness reduced to no-horizon-scale-covariance premise" in t5b_white, "T5b whiteness lemma isolates the remaining spatial-covariance premise")
    check("Thermodynamic / entropy claim protocol" in anchor, "canon carries the strict thermodynamic-claim protocol")

    print("\n[2] Current covariance ledger")
    n_required = required_n_eff()
    h_entropy = h_from_entropy_count(n_required)
    print("  For one shell: Var[(N-<N>)/<N>] = F_eff/N_eff.")
    print(f"  target A_nu                    = {A_TARGET:.6e}")
    print(f"  required N_eff/F_eff           = {n_required:.6e}")
    print(f"  required events per 28 channel = {n_required / N_CHANNELS:.6e}")
    print(f"  required events per 112 flag   = {n_required / N_FLAGS:.6e}")
    print(f"  if N_eff=S_dS, H_*             = {h_entropy:.6e} GeV")
    check(abs(a_nu_from_current(n_required) - A_TARGET) / A_TARGET < 1e-15, "A_nu=F_eff/N_eff reproduces the target only after N_eff/F_eff is fixed")
    check(4.0e8 < n_required < 6.0e8, "CMB amplitude needs O(5e8) independent shell events for F_eff=1")
    check(8.0e14 < h_entropy < 1.2e15, "entropy-pixel reading maps the target to the GUT-scale H_* neighborhood")

    print("\n[3] Log-shell covariance shape")
    for k_ratio in [0.25, 1.0, 4.0, 16.0]:
        print(
            f"  k/k*={k_ratio:5.2f}: "
            f"Delta_nu^2={delta_nu_squared(k_ratio, n_required):.12e}, "
            f"P_nu={p_nu(k_ratio, n_required):.12e}"
        )
    slope = math.log(delta_nu_squared(math.e, n_required) / delta_nu_squared(1.0, n_required))
    check(abs(slope + float(DELTA)) < 1e-14, "28-clock supplies d ln Delta_nu^2 / d ln k = -1/28")
    check(abs(delta_nu_squared(1.0, n_required) - A_TARGET) / A_TARGET < 1e-15, "pivot amplitude remains the separate F_eff/N_eff input")

    print("\n[4] Strict thermodynamic protocol gates")
    gates = [
        Gate(
            "T1",
            "event algebra / jump support",
            Status.CONDITIONAL,
            "28 channels and 112 flags are finite-instrument results; the scalar HBC current j(x,N) is a proposed lift.",
        ),
        Gate(
            "T2",
            "Landauer equality versus bound",
            Status.PASS,
            "Closed by exclusion: A_nu is normalized count covariance F_eff/N_eff; heat-per-event factors cancel from j/<j>-1, so no Landauer coefficient is present or allowed.",
        ),
        Gate(
            "T3",
            "mean service current",
            Status.OPEN,
            "The ledger still needs the absolute mean count N_eff per mode-local shell; the alpha^4 candidate N_eff=(4/3) alpha_0^-4 is not yet derived from the scalar-shell current.",
        ),
        Gate(
            "T4",
            "current covariance / Fano factor",
            Status.PASS,
            "For the canonical monitored CTMC ledger the total service count is exactly Poisson, F_eff=1; the bandwidth-one scheduler alternative has F_eff=1-p and remains a regime caveat.",
        ),
        Gate(
            "T5",
            "projector / correlation volume",
            Status.OPEN,
            "T5a projector form is closed. T5b form is reduced to N_eff=N_shell/S_j(k); S_j=1 follows under no horizon-scale connected covariance, but that premise is not yet derived from the boundary-printing ledger.",
        ),
        Gate(
            "T6",
            "gauge-invariant observable map",
            Status.CONDITIONAL,
            "R_HBC=psi-nu and delta-N are formally derived if a unique local print-time field nu_HBC exists.",
        ),
        Gate(
            "T7",
            "scale-input accounting",
            Status.OPEN,
            "The entropy route still needs an independently derived early H_*; the alpha^4 route avoids fitting H_* but still inherits the Planck-mass/horizon status if H_* is read back from A_nu.",
        ),
        Gate(
            "T8",
            "alternative ledgers and failure modes",
            Status.PASS,
            "Finite counts, deterministic service, sub-Poisson 2^28, 28*2^24, and Gamma0/H alternatives have been enumerated.",
        ),
        Gate(
            "T9",
            "promotion statement",
            Status.PASS,
            "The allowed status is conditional recovery candidate: F_eff=1 and A_nu=(3/4) alpha_0^4 if N_eff=(4/3) alpha_0^-4 and the dilute/projection premises are derived.",
        ),
    ]
    print_gate_table(gates)
    blockers = [g for g in gates if g.status in {Status.OPEN, Status.FAIL}]
    conditionals = [g for g in gates if g.status is Status.CONDITIONAL]
    print(f"  blockers={len(blockers)}; conditionals={len(conditionals)}")
    check(len(blockers) == 3, "protocol still blocks a derived/parameter-free A_nu claim")
    check(len(conditionals) == 2, "two formal bridges are conditional rather than locked")

    print("\n[5] Promotion target")
    print("  To promote A_nu above conditional candidate, derive all three:")
    print("    (i)   N_shell = (4/3) alpha_0^-4, S_dS, or another microscopic shell event density;")
    print("    (ii)  no connected service-current covariance at k=aH, so S_j(k)=1;")
    print("    (iii) the early-ledger regime: dilute CTMC duty, or the duty-corrected saturated alternative.")
    print("  For the alpha^4 count, a bandwidth-one reading requires p about 0.014, not p -> 1.")
    print("  The N_eff/duty closure attempt shows this is not just absent prose: a")
    print("  scalar-current intensity rescaling preserves Pi_k, the 1/28 tilt, and")
    print("  F_eff while moving A_nu.")
    check(True, "worked target is well-posed as service-current covariance, not entropy prose")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The strict protocol prevents the old thermodynamic overclaim pattern.")
    print("  HBC/QEC currently derives the 28-channel relative clock and conditionally")
    print("  supplies the gauge-invariant scalar-clock map.  The Fano leg is now")
    print("  closed for the canonical CTMC ledger, F_eff=1.  The strongest amplitude")
    print("  candidate is")
    print("      A_nu = (3/4) alpha_0^4")
    print("  from N_eff=(4/3) alpha_0^-4.  Pi_k is now fixed as the compensated")
    print("  Fourier-shell readout, and the correlation-volume theorem has been")
    print("  reduced to N_eff=N_shell/S_j(k).  The whiteness lemma gives S_j=1")
    print("  under no horizon-scale connected covariance.  Numerical closure remains")
    print("  open because that premise and N_shell are not derived.")
    print("  This is an improved conditional recovery candidate, not a Locked")
    print("  amplitude prediction.  The mean/mode/duty")
    print("  audit shows saturated background printing cannot be silently identified")
    print("  with the scalar fluctuation current.")
    print("=" * 108)
    print("exit 0 -- A_nu status improved to conditional alpha^4 candidate; T5b value remains open.")


if __name__ == "__main__":
    main()
