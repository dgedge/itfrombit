#!/usr/bin/env python3
"""ITEM 131: residual consolidation around the saturation premise.

The Item-131 cluster has accumulated several narrow residuals:

  * R4 d=1 support for the late w(a) branch;
  * finite-to-cosmological lift premises for n_s=27/28;
  * scalar-amplitude shell count / stop-rule premises.

This audit does not add a new derivation.  It checks the owner scripts and
records the current simplification: after the colour-budget-capacity result,
the tilt and amplitude defences share one named premise,

    saturated printer dynamics: lambda_shell = N_shell alpha0^4 = C_F.

Under that premise, the amplitude count and the constant-H/log-shell tilt are
not independent knobs.  Locking item 131 now means deriving that saturation
identity from the printer dynamics itself.
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1 / 137.0
C_F = 4 / 3


def run_script(name: str) -> str:
    proc = subprocess.run(
        [sys.executable, str(ROOT / "python_code" / name)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f"{name} failed with exit {proc.returncode}")
    return proc.stdout


def require(text: str, needle: str, label: str) -> None:
    ok = needle in text
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    if not ok:
        raise AssertionError(f"{label}: missing {needle!r}")


def main() -> None:
    print("ITEM 131 SATURATION RESIDUALS AUDIT")
    print("Rerunning owner scripts for the late, tilt, amplitude, and C_F legs.")

    outputs = {
        "r4_support": run_script("item131_r4_support_dimension.py"),
        "late_activation": run_script("item131_late_activation.py"),
        "early_logscale": run_script("item131_early_logscale_lift.py"),
        "radial_crossing": run_script("item131_mode_local_radial_crossing.py"),
        "t5b_whiteness": run_script("item131_t5b_whiteness_lemma.py"),
        "saturation": run_script("item131_saturation_closure.py"),
        "saturated_dynamics": run_script("item131_saturated_printer_dynamics_attempt.py"),
        "cf_stop": run_script("item131_cf_stop_rule_closure_attempt.py"),
    }

    print("\n[1] Owner-script gates")
    require(
        outputs["r4_support"],
        "exit 0 -- R4 one-dimensional support lemma proved at finite support level.",
        "R4 d=1 support is finite-closed",
    )
    require(
        outputs["late_activation"],
        "exit 0 -- late activation law closed at finite R4 support level.",
        "late f(a)=a / w(a) branch closes under homogeneous line-ledger lift",
    )
    require(
        outputs["early_logscale"],
        "exit 0 -- early log-scale lift reduced to scale-covariant HBC.",
        "finite-to-log lift is reduced to scale-covariant saturated HBC",
    )
    require(
        outputs["radial_crossing"],
        "exit 0 -- q=1 mode-local radial-crossing lemma derived as a spectral-density normalization theorem.",
        "mode-local q=1 radial crossing is formal-closed",
    )
    require(
        outputs["t5b_whiteness"],
        "no-horizon-scale-covariance premise",
        "T5b residue is no connected horizon-scale service-current covariance",
    )
    require(
        outputs["saturation"],
        "one shared, content-fixed",
        "amplitude and tilt now share one content-fixed saturation premise",
    )
    require(
        outputs["saturated_dynamics"],
        "saturated-printer dynamics reduced to one maximal-throughput/backlog-latch theorem; not Locked.",
        "printer dynamics attempt leaves the saturation equality conditional",
    )
    require(
        outputs["cf_stop"],
        "exit 0 -- C_F load structurally selected; stop-rule marginality remains open.",
        "C_F coefficient selected; printer stop-rule remains the open dynamical law",
    )

    n_shell = C_F / ALPHA0**4
    a_nu = 1 / n_shell

    print("\n[2] Shared saturation identity")
    print("  Define lambda_shell = N_shell alpha0^4.")
    print("  Sustainability of colour-coherent printing gives lambda_shell <= C_F.")
    print("  Saturated capacity printing gives lambda_shell >= C_F.")
    print("  Therefore, if the saturation-printer theorem is true:")
    print(f"      lambda_shell = C_F = {C_F:.12f}")
    print(f"      N_shell = C_F alpha0^-4 = {n_shell:.6e}")
    print(f"      A_nu = 1/N_shell = (3/4) alpha0^4 = {a_nu:.6e}")
    assert math.isclose(a_nu, 0.75 * ALPHA0**4, rel_tol=1e-15, abs_tol=0.0)

    print("\n[3] Residual collapse table")
    rows = [
        (
            "R4 d=1 support",
            "FINITE-CLOSED",
            "R4 repair support is a 1-chain; cosmology still assumes homogeneous R4 line-ledger lift.",
        ),
        (
            "late w(a)",
            "CONDITIONAL-LIFT",
            "f(a)=a and w(a)=-1+a/28 follow under saturated homogeneous line support and no outside positive channel.",
        ),
        (
            "tilt n_s",
            "CONDITIONAL-LIFT",
            "finite 28-clock + q=1 radial log shell give 27/28 when the saturated printer supplies the scalar clock.",
        ),
        (
            "amplitude A_nu",
            "CONDITIONAL-LIFT",
            "F_eff=1 and S_j=1 leave N_shell; saturation-capacity gives N_shell=C_F alpha0^-4.",
        ),
        (
            "finite-to-cosmological lift",
            "SINGLE-RESIDUAL",
            "derive saturated printer dynamics: capacity, backlog feedback, constant-H phase, and exit.",
        ),
    ]
    for name, status, note in rows:
        print(f"  {name:28s} {status:18s} {note}")

    print("\n[4] The Lock target")
    print("  Prove a saturated-printer dynamics theorem from HBC/QEC mechanics:")
    print("    1. the boundary printer is a colour-load regulator with capacity C_F;")
    print("    2. coherent printing is admitted only for lambda_shell <= C_F;")
    print("    3. the early saturated phase is backlog-free maximal printing, so lambda_shell >= C_F;")
    print("    4. hence lambda_shell = C_F while saturated, fixing H*, n_s, and A_nu together;")
    print("    5. exit occurs when the printer leaves this capacity-saturated branch.")
    print("  The present dynamics attempt derives the capacity inequality and the")
    print("  conditional maximum-throughput corollary, but not this latch law.  The")
    print("  cleanup is that this remains the one named residual under both tilt and")
    print("  amplitude.")

    print("\nVERDICT")
    print("  Item 131 is tidier but not Locked.  R4 d=1 is finite-closed; tilt and")
    print("  amplitude both reduce to the same saturation premise; deriving saturation")
    print("  from printer dynamics would Lock both at once.")
    print("exit 0 -- item-131 residuals consolidated to the saturated-printer theorem target.")


if __name__ == "__main__":
    main()
