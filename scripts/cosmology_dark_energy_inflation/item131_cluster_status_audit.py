#!/usr/bin/env python3
"""Consolidated status audit for the item-131 HBC/R4 cluster.

This is intentionally a meta-audit: it reruns the small, targeted scripts that
own each leg and records the boundary between finite-instrument theorems,
conditional cosmology lifts, and genuinely open stop-rule/threshold lemmas.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_script(name: str) -> str:
    path = ROOT / "python_code" / name
    proc = subprocess.run(
        [sys.executable, str(path)],
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
    if needle not in text:
        raise AssertionError(f"{label}: missing {needle!r}")
    print(f"  [PASS] {label}")


def main() -> None:
    print("ITEM 131 CLUSTER STATUS AUDIT")
    print("Rerunning owner scripts for the three live subclaims.")

    scripts = {
        "r4_support": "item131_r4_support_dimension.py",
        "late_activation": "item131_late_activation.py",
        "early_logscale": "item131_early_logscale_lift.py",
        "radial_crossing": "item131_mode_local_radial_crossing.py",
        "scalar_clock": "item131_scalar_clock_bridge.py",
        "physical_coupling": "item131_hbc_physical_coupling_audit.py",
        "clock_covariance_exit": "item131_hbc_clock_covariance_exit.py",
        "inflation_gate": "item131_inflation_hbc_status_gate.py",
        "saturation": "item131_saturation_closure.py",
        "saturated_dynamics": "item131_saturated_printer_dynamics_attempt.py",
        "cf_stop_rule": "item131_cf_stop_rule_closure_attempt.py",
    }

    outputs = {key: run_script(script) for key, script in scripts.items()}

    print("\n[1] Late R4 / w(a) branch")
    require(
        outputs["r4_support"],
        "exit 0 -- R4 one-dimensional support lemma proved at finite support level.",
        "R4 one-dimensional support is finite-proved",
    )
    require(
        outputs["late_activation"],
        "exit 0 -- late activation law closed at finite R4 support level.",
        "late f(a)=a branch follows under homogeneous line-ledger lift",
    )
    print("  status: CLOSED at finite support level; cosmology lift assumes the homogeneous R4 line ledger and no outside-sector positive channel.")

    print("\n[2] Primordial tilt / n_s branch")
    require(
        outputs["early_logscale"],
        "exit 0 -- early log-scale lift reduced to scale-covariant HBC.",
        "scale-covariant Markov horizon ratios force the log clock",
    )
    require(
        outputs["radial_crossing"],
        "exit 0 -- q=1 mode-local radial-crossing lemma derived as a spectral-density normalization theorem.",
        "q=1 radial mode-local normalization is formal-closed",
    )
    require(
        outputs["scalar_clock"],
        "exit 0 -- scalar-clock wrapper closed conditionally; HBC clock covariance and early-stage duration remain open.",
        "gauge-invariant R_HBC/delta-N wrapper is conditional",
    )
    require(
        outputs["physical_coupling"],
        "saturated-H remains a premise",
        "physical HBC scalar-clock coupling remains a premise",
    )
    print("  status: CONDITIONAL THEOREM: finite 28-clock + radial log-shell + power-ledger action give n_s=27/28 if saturated HBC supplies the unique scalar clock.")

    print("\n[3] Scalar amplitude / A_nu branch")
    require(outputs["inflation_gate"], "F_eff=1", "CTMC Fano leg remains closed")
    require(outputs["inflation_gate"], "absolute shell count", "absolute shell count remains a hard gate")
    require(
        outputs["saturation"],
        "one shared, content-fixed",
        "saturation-capacity identification collapses amplitude and tilt onto one premise",
    )
    require(
        outputs["saturated_dynamics"],
        "not Locked",
        "saturated-printer dynamics attempt keeps the equality conditional",
    )
    require(
        outputs["cf_stop_rule"],
        "exit 0 -- C_F load structurally selected; stop-rule marginality remains open.",
        "C_F coefficient is selected but the marginal stop rule is open",
    )
    print("  status: CONDITIONAL CANDIDATE: A_nu=(3/4)alpha0^4 follows if saturated printer dynamics enforces lambda_shell=N_shell alpha0^4=C_F.")

    print("\n[4] Remaining theorem targets")
    rows = [
        (
            "T_amp_stop",
            "SHARED-OPEN",
            "derive saturated-printer dynamics: lambda_shell = E[N_topological commits per printed shell] = C_F",
        ),
        (
            "T_scalar_clock",
            "CONDITIONAL",
            "derive the local HBC print-time field and its covariance from the event ledger",
        ),
        (
            "T_constant_H_exit",
            "OPEN",
            "derive finite saturated w=-1 duration and exit/reheating from the same printer dynamics",
        ),
        (
            "T_late_R4",
            "FINITE-CLOSED",
            "R4 support dimension d=1; residual is outside-sector/non-R4-coupling completeness",
        ),
    ]
    for name, status, target in rows:
        print(f"  {name:18s} {status:14s} {target}")

    print("\nVERDICT")
    print("  The item-131 cluster is not uniformly open.")
    print("  R4 one-dimensional support is finite-closed; n_s=27/28 is a conditional")
    print("  HBC scalar-clock/log-shell theorem; A_nu is now conditional on the same")
    print("  saturation-capacity premise.  The remaining Lock target is one theorem:")
    print("  derive maximal-throughput/backlog-latch printer dynamics from HBC/QEC,")
    print("  including constant-H duration, no-horizon covariance, and exit.")
    print("exit 0 -- item-131 cluster status consolidated.")


if __name__ == "__main__":
    main()
