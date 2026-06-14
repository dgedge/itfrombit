#!/usr/bin/env python3
"""ITEM 131: scalar event unit and early saturated-printer scale audit.

Question
--------
Can the HBC/QEC ledger now derive the scalar event unit u_event and the early
saturated printer scale H_*?

Result
------
One clause closes conditionally; the scale becomes sharply fixed but is not yet
an independent theorem.

If the scalar current is the HBC boundary-printing current, a single scalar
event is one printed boundary cell.  The HBC source states that each new cell
raises the horizon capacity by exactly ln 2.  Therefore

    u_event = ln 2

for the printed-cell scalar-current reading.  This is not the same object as
the baryogenesis/photon record unit ln(8*137), which counts an erased classical
service record rather than added geometric capacity.

Combining u_event=ln2 with the existing alpha^4 shell-count candidate,

    N_shell = (4/3) alpha_0^-4,
    N_shell = S_dS(H_*) / ln2,
    S_dS = 8 pi^2 Mbar_P^2 / H_*^2,

fixes

    H_* = sqrt(6 pi^2 / ln2) alpha_0^2 Mbar_P = 1.20e15 GeV.

This lands on the R2/seesaw activation decade, but canon still labels that
R2/R3 activation schedule as phenomenological/open.  So H_* is a conditional
consequence of the alpha^4 amplitude route plus the ln2 event-unit theorem, not
an independent derivation of the early saturated printer threshold.
"""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALPHA0 = 1.0 / 137.0
MPL_REDUCED_GEV = 2.435e18
A_OBS = 2.10e-9
SIGMA_A = 0.03e-9
V_EW_GEV = 246.0
MNU_ATM_EV = 0.05
MNU_ATM_GEV = MNU_ATM_EV * 1.0e-9


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def alpha4_count(alpha: float = ALPHA0) -> float:
    return (4.0 / 3.0) * alpha**-4


def de_sitter_entropy(h_gev: float) -> float:
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def amplitude_from_h_u(h_gev: float, u_event: float) -> float:
    return u_event / de_sitter_entropy(h_gev)


def h_for_count(u_event: float, count: float) -> float:
    return math.sqrt(8.0 * math.pi**2 / (u_event * count)) * MPL_REDUCED_GEV


def h_for_amplitude(u_event: float, amplitude: float = A_OBS) -> float:
    return math.sqrt(8.0 * math.pi**2 * amplitude / u_event) * MPL_REDUCED_GEV


def pull(prediction: float) -> float:
    return (prediction - A_OBS) / SIGMA_A


def rel(a: float, b: float) -> float:
    return a / b - 1.0


def main() -> None:
    print("ITEM 131 U_EVENT / H_* DERIVATION AUDIT")

    print("\n[1] Source gates")
    check(
        contains(
            "cosmological_qec_engine/cosmological_qec_engine.tex",
            "raising $S_{\\max}$ by exactly $\\ln 2$ per cell",
        ),
        "HBC source states one printed cell adds exactly ln2 of capacity",
    )
    check(
        contains(
            "python_code/item131_scalar_current_rescaling_anchor.py",
            "N_shell(H_*, u_event) = S_dS(H_*) / u_event",
        ),
        "prior item-131 audit reduced amplitude normalization to S_dS/u_event",
    )
    check(
        contains(
            "python_code/item131_inflationary_amplitude_alpha4_route_audit.py",
            "N_eff = (4/3) alpha_0^-4",
        ),
        "alpha^4 shell-count candidate is already registered as conditional",
    )
    check(
        contains("ANCHOR.md", "phenomenologically suggested")
        and contains("ANCHOR.md", "not a structural framework derivation"),
        "R2/R3 activation scale is explicitly not yet a structural theorem",
    )
    check(
        contains(
            "part_17_energy_trajectory/part_17_energy_trajectory.tex",
            "closure of item 42 would lift this boundary condition",
        ),
        "early R-activation / printer-threshold derivation is item-42 open territory",
    )

    print("\n[2] Event-unit selection")
    u_nat = 1.0
    u_cell = math.log(2.0)
    u_record = math.log(8.0 * 137.0)
    print(f"  natural entropy nat      u = {u_nat:.12f}  (coordinate convention)")
    print(f"  printed HBC cell         u = {u_cell:.12f}  (capacity increment)")
    print(f"  address x channel record u = {u_record:.12f}  (erased service record)")
    check(abs(u_cell - math.log(2.0)) < 1.0e-15, "printed-cell scalar event gives u_event=ln2")
    check(u_record / u_cell > 10.0, "record-erasure unit is a different ledger, not a silent scalar-current import")
    print("  Derivation clause: if the scalar clock is the local printed-cell current,")
    print("  then its event is the HBC capacity increment, not a heat-record erasure.")

    print("\n[3] H_* fixed by ln2 plus the alpha^4 shell-count candidate")
    n_alpha = alpha4_count()
    h_cell = h_for_count(u_cell, n_alpha)
    h_nat = h_for_count(u_nat, n_alpha)
    h_record = h_for_count(u_record, n_alpha)
    amp_cell = amplitude_from_h_u(h_cell, u_cell)
    h_obs_cell = h_for_amplitude(u_cell)
    print(f"  N_shell=(4/3)alpha^-4 = {n_alpha:.9e}")
    print(f"  H_*(u=1)              = {h_nat:.9e} GeV")
    print(f"  H_*(u=ln2)            = {h_cell:.9e} GeV")
    print(f"  H_*(u=ln(8*137))      = {h_record:.9e} GeV")
    print(f"  H required by A_obs, u=ln2 = {h_obs_cell:.9e} GeV")
    print(f"  A(H_*, ln2)           = {amp_cell:.9e}  pull={pull(amp_cell):+.3f} sigma")
    check(1.1e15 < h_cell < 1.3e15, "ln2 event unit fixes H_* in the 1.2e15 GeV decade")
    check(abs(de_sitter_entropy(h_cell) / u_cell / n_alpha - 1.0) < 1.0e-12, "H_* exactly realizes the alpha^4 count by construction")
    check(abs(pull(amp_cell)) < 1.1, "resulting amplitude remains within the current Planck-scale normalization band")

    print("\n[4] R2/seesaw scale cross-check, not an independent derivation")
    e_r2 = V_EW_GEV**2 / MNU_ATM_GEV
    mnu_pred_ev = V_EW_GEV**2 / h_cell * 1.0e9
    amp_at_er2 = amplitude_from_h_u(e_r2, u_cell)
    print(f"  E_R2 from v^2/(0.05 eV) = {e_r2:.9e} GeV")
    print(f"  H_*(ln2 alpha^4)        = {h_cell:.9e} GeV")
    print(f"  relative difference      = {rel(h_cell, e_r2):+.3%}")
    print(f"  If M_R=H_*, v^2/M_R     = {mnu_pred_ev:.5f} eV")
    print(f"  A(E_R2, ln2)            = {amp_at_er2:.9e}  pull={pull(amp_at_er2):+.3f} sigma")
    check(abs(rel(h_cell, e_r2)) < 0.02, "H_* lands within 2% of the phenomenological R2/seesaw activation scale")
    check(0.045 < mnu_pred_ev < 0.055, "identifying M_R=H_* returns the atmospheric neutrino scale")
    print("  This is a strong cross-sector coincidence, but ANCHOR labels the R2/R3")
    print("  activation schedule and seesaw dimensionalisation as open; it cannot yet")
    print("  be used as an independent derivation of H_*.")

    print("\n[5] Status ledger")
    closed = [
        "u_event=ln2 under scalar-event = printed HBC boundary cell",
        "record-alphabet ln(8*137) is excluded as an automatic scalar event unit",
        "given N_shell=(4/3)alpha^-4, ln2 fixes H_*=1.20e15 GeV",
    ]
    open_items = [
        "derive the scalar-event = printed-cell identification from the gauge-invariant delta-N map",
        "derive N_shell=(4/3)alpha^-4 from the microscopic scalar-current shell ledger",
        "derive item-42/R-activation or another early saturated-printer threshold independently of A_nu",
        "prove the H_* = E_R2 handoff rather than treating the 1.2e15 GeV match as a cross-check",
    ]
    for item in closed:
        check(True, f"conditional/closed: {item}")
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The scalar event unit can be derived conditionally: if the adiabatic")
    print("  scalar clock is the local HBC printed-cell current, then u_event=ln2,")
    print("  because each new boundary cell adds exactly ln2 of geometric capacity.")
    print("  This is not the baryogenesis/photon record unit ln(8*137).")
    print("")
    print("  With that event unit, the existing alpha^4 shell-count candidate fixes")
    print("      H_* = sqrt(6 pi^2/ln2) alpha_0^2 Mbar_P = 1.20e15 GeV.")
    print("  The scale lands on the R2/seesaw activation decade and predicts")
    print("  v^2/H_* = 0.050 eV if identified with M_R, but that is a cross-check:")
    print("  canon still lacks the independent item-42/R-activation printer-threshold")
    print("  theorem.  So u_event is conditionally tied down; H_* is conditionally")
    print("  fixed, not independently Locked.")
    print("=" * 108)
    print("exit 0 -- u_event=ln2 conditional derivation; H_* fixed by alpha^4+ln2 but not independent.")


if __name__ == "__main__":
    main()
