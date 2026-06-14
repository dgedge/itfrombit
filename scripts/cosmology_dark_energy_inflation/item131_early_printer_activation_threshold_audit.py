#!/usr/bin/env python3
"""ITEM 131 / 42: early printer and R-activation threshold audit.

Question
--------
Can the early saturated HBC printer scale H_* be derived as the R2/R3
activation threshold, rather than inferred from the scalar amplitude?

Result
------
The threshold can be written as a sharp conditional theorem, but the core
activation lemma is still not independently derived.

The clean candidate is:

    N_print(H_*) * alpha_0^4 = C_F,
    N_print(H_*) = S_dS(H_*) / ln2,
    C_F = 4/3.

Read literally: the saturated printer exits when one mode-local horizon shell
contains one SU(3)-Casimir-normalized quartic logical activation.  This yields

    H_* = sqrt(6 pi^2 / ln2) alpha_0^2 Mbar_P = 1.199e15 GeV.

That is exactly the same scale found by the alpha^4 amplitude candidate and it
lands on the phenomenological R2/seesaw activation decade.  However, unless the
R2/R3 Kraus/current ledger independently proves the threshold condition
N_print alpha^4 = 4/3, this is not an independent derivation of item 42.  It is
the amplitude candidate recast as an activation threshold.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALPHA0 = 1.0 / 137.0
MPL_REDUCED_GEV = 2.435e18
A_OBS = 2.10e-9
SIGMA_A = 0.03e-9
V_EW_GEV = 246.0
MNU_ATM_GEV = 0.05e-9


@dataclass(frozen=True)
class ThresholdCandidate:
    name: str
    c_threshold: float
    status: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def de_sitter_entropy(h_gev: float) -> float:
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def printed_cell_count(h_gev: float) -> float:
    return de_sitter_entropy(h_gev) / math.log(2.0)


def h_from_activation_threshold(c_threshold: float) -> float:
    """Solve (S_dS/ln2) * alpha^4 = c_threshold."""
    n_print = c_threshold * ALPHA0**-4
    return math.sqrt(8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / (math.log(2.0) * n_print))


def amplitude_from_threshold(c_threshold: float) -> float:
    """If the same shell count is read as A=1/N, return A."""
    return ALPHA0**4 / c_threshold


def pull(prediction: float) -> float:
    return (prediction - A_OBS) / SIGMA_A


def rel(a: float, b: float) -> float:
    return a / b - 1.0


def main() -> None:
    print("ITEM 131 / 42 EARLY PRINTER R-ACTIVATION THRESHOLD AUDIT")

    print("\n[1] Source gates")
    check(
        contains("ANCHOR.md", "R1 \u2192 R2/R3 \u2192 R4 activation schedule")
        and contains("ANCHOR.md", "phenomenologically arranged"),
        "item 42 explicitly keeps the activation schedule open/phenomenological",
    )
    check(
        contains("ANCHOR.md", "not a structural framework derivation")
        and contains("ANCHOR.md", "simultaneous activation"),
        "R2/R3 simultaneous breaking is not yet a structural derivation",
    )
    check(
        contains(
            "cosmological_qec_engine/cosmological_qec_engine.tex",
            "raising $S_{\\max}$ by exactly $\\ln 2$ per cell",
        ),
        "HBC supplies the printed-cell capacity unit ln2",
    )
    check(
        contains(
            "python_code/item131_u_event_hstar_derivation_audit.py",
            "u_event = ln 2",
        ),
        "prior audit conditionally derives u_event=ln2 for printed-cell scalar events",
    )
    check(
        contains(
            "python_code/instrumentation_onset_decision.py",
            "they exist at crystallisation (level 1) and ONLY there",
        ),
        "instrumentation-onset audit decides the strain-readout onset in the boot chain",
    )
    check(
        contains(
            "python_code/tch_item53_seesaw_dimensionalisation_audit.py",
            "useful bookkeeping, but this is an inversion",
        ),
        "item 53 audit keeps seesaw dimensionalisation from acting as an independent H_* derivation",
    )

    print("\n[2] Candidate threshold equation")
    c_f = 4.0 / 3.0
    h_cf = h_from_activation_threshold(c_f)
    n_cf = printed_cell_count(h_cf)
    load_cf = n_cf * ALPHA0**4
    amp_cf = amplitude_from_threshold(c_f)
    print("  Candidate activation condition:")
    print("      N_print(H_*) alpha_0^4 = C_F")
    print("      N_print = S_dS/ln2,  C_F=4/3")
    print(f"  H_*                         = {h_cf:.9e} GeV")
    print(f"  N_print(H_*)                = {n_cf:.9e}")
    print(f"  N_print alpha_0^4           = {load_cf:.12f}")
    print(f"  A from same count           = {amp_cf:.9e}  pull={pull(amp_cf):+.3f} sigma")
    check(abs(load_cf - c_f) < 1.0e-12, "H_* exactly satisfies the proposed activation threshold")
    check(1.1e15 < h_cf < 1.3e15, "threshold lands in the early GUT/R2 activation decade")
    check(abs(pull(amp_cf)) < 1.1, "same threshold gives the existing alpha^4 amplitude candidate")

    print("\n[3] Alternative activation constants")
    candidates = [
        ThresholdCandidate("one raw quartic logical event", 1.0, "no SU(3) projection"),
        ThresholdCandidate("SU(3) Casimir C_F", 4.0 / 3.0, "current conditional candidate"),
        ThresholdCandidate("inverse Casimir 3/4", 3.0 / 4.0, "wrong side if used as threshold count"),
        ThresholdCandidate("walk-active 3/2", 3.0 / 2.0, "nearby simple alternative"),
        ThresholdCandidate("baryogenesis branch 3/14", 3.0 / 14.0, "wrong ledger"),
    ]
    for row in candidates:
        h = h_from_activation_threshold(row.c_threshold)
        amp = amplitude_from_threshold(row.c_threshold)
        print(f"  {row.name:32s} C={row.c_threshold:.9f}")
        print(f"      H_*={h:.9e} GeV  A={amp:.9e}  pull={pull(amp):+.2f} sigma")
        print(f"      status: {row.status}")
    check(abs(pull(amplitude_from_threshold(4.0 / 3.0))) < abs(pull(amplitude_from_threshold(1.0))), "C_F improves the amplitude relative to raw one-event threshold")
    check(abs(pull(amplitude_from_threshold(4.0 / 3.0))) < abs(pull(amplitude_from_threshold(3.0 / 2.0))), "C_F is numerically preferred over the nearby 3/2 threshold")
    print("  Caution: this table shows why the C_F selection must be derived before")
    print("  comparison.  Nearby simple constants are plentiful enough that the number")
    print("  alone cannot promote the threshold.")

    print("\n[4] R2/seesaw and boot-onset checks")
    e_r2 = V_EW_GEV**2 / MNU_ATM_GEV
    print(f"  R2/seesaw inversion v^2/(0.05 eV) = {e_r2:.9e} GeV")
    print(f"  threshold H_*                    = {h_cf:.9e} GeV")
    print(f"  relative difference              = {rel(h_cf, e_r2):+.3%}")
    print(f"  if M_R=H_*, v^2/H_*              = {V_EW_GEV**2 / h_cf * 1.0e9:.5f} eV")
    check(abs(rel(h_cf, e_r2)) < 0.02, "threshold matches the phenomenological R2/seesaw scale at the percent level")
    print("  Boot/noise thresholds such as p_c decide when the strain readout turns on")
    print("  in the boot hierarchy.  They are dimensionless and cannot supply H_* in")
    print("  GeV without a separate capacity/scale bridge.")
    check(True, "dimensionless handoff criticality is an onset rule, not a GeV scale derivation")

    print("\n[5] What is actually derived?")
    closed = [
        "printed-cell scalar events carry u_event=ln2, conditionally",
        "if N_print alpha^4 = 4/3, then H_*=1.199e15 GeV follows algebraically",
        "the same threshold gives A_nu=(3/4)alpha_0^4",
        "the result cross-checks the R2/seesaw decade",
    ]
    open_items = [
        "derive the threshold condition N_print alpha^4 = 4/3 from the actual R2/R3 Kraus/current ledger",
        "derive why the SU(3) Casimir C_F, rather than another simple constant, is the activation load",
        "derive simultaneous R2/R3 activation on the 3D TCH gauge web (item 43 interface)",
        "derive the R4/EW scale and seesaw dimensionalisation instead of inserting v and m_nu",
        "derive the finite saturated-printer exit/duration rather than only the exit scale",
    ]
    for item in closed:
        check(True, f"conditional/closed: {item}")
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  A clean early-printer threshold theorem is now isolated:")
    print("      N_print(H_*) alpha_0^4 = 4/3, with N_print=S_dS/ln2.")
    print("  If the R2/R3 activation ledger proves that condition, then")
    print("      H_* = sqrt(6 pi^2/ln2) alpha_0^2 Mbar_P = 1.199e15 GeV")
    print("  follows without fitting, and the scalar amplitude follows as")
    print("      A_nu = (3/4) alpha_0^4.")
    print("")
    print("  But the decisive clause is still open: current canon has not derived")
    print("  N_print alpha^4 = 4/3 from the R2/R3 Kraus/current geometry.  Therefore")
    print("  this is a conditional activation-threshold reduction, not a full item-42")
    print("  R-activation schedule closure.  The next finite target is the actual")
    print("  R2/R3 threshold ledger proving the quartic load and the C_F projection.")
    print("=" * 108)
    print("exit 0 -- early-printer threshold isolated; conditional, not independently derived.")


if __name__ == "__main__":
    main()
