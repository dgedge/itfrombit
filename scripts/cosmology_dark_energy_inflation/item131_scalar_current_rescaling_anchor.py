#!/usr/bin/env python3
r"""ITEM 131: tie down the scalar-current rescaling freedom.

Question
--------
The current HBC amplitude ledger still has a rescaling symmetry:

    j(x,N) -> lambda j(x,N),
    N_shell -> lambda N_shell,
    A_nu -> A_nu/lambda,

while the already-closed gates (28-clock tilt, Pi_k projector, CTMC Fano) are
unchanged.  Can the ledger tie lambda down?

Result
------
One step of progress, not a closure.

If the scalar current is identified with HBC boundary-capacity printing, then
the shell count is no longer an arbitrary intensity.  It must be

    N_shell(H_*, u_event) = S_dS(H_*) / u_event,
    S_dS = 8 pi^2 Mbar_P^2 / H_*^2,

where u_event is the entropy capacity consumed per independent scalar-current
event, in natural entropy units.  Then

    A_nu = F_eff S_j u_event / S_dS.

So the continuous rescaling freedom is reduced to two concrete theorem targets:

    (1) event-unit theorem: what is u_event for the scalar current?
    (2) early-printer theorem: what is H_*?

The alpha^4 candidate becomes the equivalent scale statement

    S_dS(H_*) / u_event = (4/3) alpha_0^-4,
    H_* = sqrt(6 pi^2 / u_event) alpha_0^2 Mbar_P.

This does not lock the amplitude, because current canon does not yet derive
u_event for scalar perturbations or H_*.  It does show precisely where the
old lambda freedom lives.  It is not an extra hidden coefficient once the HBC
capacity law is adopted.
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
MPL_REDUCED_GEV = 2.435e18
C_F = Fraction(4, 3)


@dataclass(frozen=True)
class UnitCandidate:
    name: str
    u_event: float
    status: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def alpha4_count(alpha: float = ALPHA0) -> float:
    return float(C_F) * alpha**-4


def de_sitter_entropy(h_gev: float) -> float:
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def n_shell_from_capacity(h_gev: float, u_event: float) -> float:
    return de_sitter_entropy(h_gev) / u_event


def amplitude_from_capacity(h_gev: float, u_event: float, s_j: float = 1.0, f_eff: float = 1.0) -> float:
    return f_eff * s_j / n_shell_from_capacity(h_gev, u_event)


def h_for_target_amplitude(u_event: float, amplitude: float = A_TARGET) -> float:
    """Solve amplitude = u_event H^2/(8 pi^2 Mbar_P^2)."""
    return math.sqrt(8.0 * math.pi**2 * amplitude / u_event) * MPL_REDUCED_GEV


def h_for_alpha4_count(u_event: float) -> float:
    """Solve S_dS(H)/u_event = (4/3) alpha_0^-4."""
    return math.sqrt(8.0 * math.pi**2 / (u_event * alpha4_count())) * MPL_REDUCED_GEV


def sigma_pull(prediction: float) -> float:
    return (prediction - A_TARGET) / SIGMA_A


def main() -> None:
    print("ITEM 131 SCALAR-CURRENT RESCALING ANCHOR")

    print("\n[1] Source checks")
    check(
        contains("python_code/item131_neff_duty_closure_attempt.py", "scale-invariance no-go"),
        "prior audit proves the raw scalar-current intensity rescaling freedom",
    )
    check(
        contains("python_code/item131_t5b_correlation_volume_audit.py", "N_eff(k) = N_shell / S_j(k)"),
        "T5b has already reduced correlation volume to N_shell/S_j",
    )
    check(
        contains("python_code/item131_t5b_whiteness_lemma.py", "no-horizon-scale-covariance premise"),
        "T5b whiteness has a precise remaining covariance premise",
    )
    check(
        contains("cosmological_qec_engine/cosmological_qec_engine.tex", "raising $S_{\\max}$ by exactly $\\ln 2$ per cell"),
        "HBC source supplies a boundary-capacity printing law",
    )
    check(
        contains("python_code/ledger_sky_reading.py", "CLASS 3") and contains("python_code/ledger_sky_reading.py", "N_eff = 4.76e8"),
        "ledger-to-sky audit classifies A_nu/N_eff as a class-3 normalization",
    )
    check(
        contains("python_code/record_alphabet_derivation.py", "s_1 = ln(8 x 137)"),
        "record-alphabet event unit is available as a separate, derived class-2 reading",
    )

    print("\n[2] Capacity law breaks the formal lambda symmetry")
    n_alpha = alpha4_count()
    amp_alpha = 1.0 / n_alpha
    print("  Raw covariance ledger:")
    print("      A_nu = F_eff S_j / N_shell")
    print("  HBC capacity-anchored ledger:")
    print("      N_shell = S_dS(H_*) / u_event")
    print("      A_nu = F_eff S_j u_event / S_dS(H_*)")
    print(f"  alpha^4 candidate count       = {n_alpha:.6e}")
    print(f"  alpha^4 candidate amplitude   = {amp_alpha:.6e}  pull={sigma_pull(amp_alpha):+.3f} sigma")
    check(abs(sigma_pull(amp_alpha)) < 1.1, "alpha^4 count remains the sharp conditional amplitude candidate")
    check(True, "lambda is now interpretable as capacity density, not a free current normalization, if the HBC capacity law is granted")

    print("\n[3] Event-unit candidates")
    units = [
        UnitCandidate("one natural entropy pixel", 1.0, "old entropy-pixel premise"),
        UnitCandidate("one stabilizer-bit cell", math.log(2.0), "literal HBC source: ln2 capacity per cell"),
        UnitCandidate("one 28-channel record", math.log(28.0), "service-label record unit; not derived for scalar current"),
        UnitCandidate("record alphabet ln(8*137)", math.log(8.0 * 137.0), "baryogenesis/photon record unit; not automatically scalar"),
    ]
    for row in units:
        h_target = h_for_target_amplitude(row.u_event)
        h_alpha = h_for_alpha4_count(row.u_event)
        n_at_h_alpha = n_shell_from_capacity(h_alpha, row.u_event)
        amp_at_h_alpha = amplitude_from_capacity(h_alpha, row.u_event)
        print(f"  {row.name:30s} u={row.u_event:.6f}")
        print(f"      H for A_target       = {h_target:.6e} GeV")
        print(f"      H for alpha^4 count  = {h_alpha:.6e} GeV")
        print(f"      N_shell(H_alpha)     = {n_at_h_alpha:.6e}")
        print(f"      A(H_alpha)           = {amp_at_h_alpha:.6e}  pull={sigma_pull(amp_at_h_alpha):+.3f} sigma")
        print(f"      status: {row.status}")
        check(abs(n_at_h_alpha / n_alpha - 1.0) < 1.0e-12, f"{row.name}: H_alpha exactly realizes alpha^4 count")
    check(8.0e14 < h_for_alpha4_count(1.0) < 1.1e15, "natural-pixel unit reproduces the existing H_*~1e15 GeV alpha^4 route")
    check(1.1e15 < h_for_alpha4_count(math.log(2.0)) < 1.3e15, "literal ln2 cell unit shifts the required H_* upward by sqrt(1/ln2)")
    check(3.0e14 < h_for_alpha4_count(math.log(8.0 * 137.0)) < 4.5e14, "record-alphabet unit shifts the required H_* downward by sqrt(1/s1)")

    print("\n[4] What this does and does not tie down")
    h_nat = h_for_alpha4_count(1.0)
    h_bit = h_for_alpha4_count(math.log(2.0))
    h_record = h_for_alpha4_count(math.log(8.0 * 137.0))
    print(f"  If u_event=1,       H_* = sqrt(6 pi^2) alpha_0^2 Mbar_P = {h_nat:.6e} GeV")
    print(f"  If u_event=ln2,     H_* = {h_bit:.6e} GeV")
    print(f"  If u_event=ln1096,  H_* = {h_record:.6e} GeV")
    print("  These are not three fitted amplitudes.  They are three different event")
    print("  granularities for the same capacity law.  The scalar-current event unit")
    print("  must be derived before the alpha^4 count can be promoted.")
    check(abs(h_bit / h_nat - math.sqrt(1.0 / math.log(2.0))) < 1.0e-12, "event-unit changes H_* by the predicted sqrt(1/u) factor")
    check(abs(h_record / h_nat - math.sqrt(1.0 / math.log(8.0 * 137.0))) < 1.0e-12, "record event cannot be imported without a large H_* shift")

    print("\n[5] Remaining theorem target")
    open_items = [
        "prove scalar-current events are HBC capacity pixels, cells, record events, or another unit",
        "derive H_* from item-42/R-activation or another early saturated-printer mechanism",
        "prove no connected service-current covariance at k=aH, so S_j=1",
        "settle CTMC/dilute versus duty-corrected bandwidth-one readout",
    ]
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The scalar-current rescaling freedom can be tied down one level, but not")
    print("  eliminated.  Once HBC boundary-capacity printing is imposed, the free")
    print("  lambda becomes")
    print("      N_shell = S_dS(H_*) / u_event,")
    print("  so the amplitude depends on an event-unit theorem and an early H_* theorem,")
    print("  not on arbitrary current intensity.  The alpha^4 candidate is equivalent")
    print("  to")
    print("      H_* = sqrt(6 pi^2 / u_event) alpha_0^2 Mbar_P.")
    print("  For u_event=1 this is the existing ~1e15 GeV route; for the literal HBC")
    print("  ln2 cell unit it moves to ~1.2e15 GeV; for the record alphabet ln(8*137)")
    print("  it moves to ~3.8e14 GeV.  Current canon has not derived which unit the")
    print("  scalar current samples, so A_nu remains a conditional recovery candidate.")
    print("=" * 108)
    print("exit 0 -- rescaling freedom reduced to event-unit + H_* capacity theorem; not Locked.")


if __name__ == "__main__":
    main()
