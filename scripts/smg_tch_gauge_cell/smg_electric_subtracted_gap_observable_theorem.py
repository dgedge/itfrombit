#!/usr/bin/env python3
r"""Observable theorem: the SMG/TCH continuum gate must use Delta_mir.

The continuum question is whether a charged mirror pole survives after the
pure gauge sector is sent to the Wilson-axis continuum.  A charged sector in a
gauge theory is not just "mirror matter": Gauss law forces a static flux/string
background.  The raw transfer-matrix gap therefore decomposes as

    Delta_raw(beta,L,cut)
      = E_full(charged) - E_full(vacuum)
      = E_string_min(beta,L) + Delta_mir(beta,L,cut).

The pure gauge term E_string_min is present even when the mirror block carries
no mass witness.  Therefore a theorem about mirror decoupling cannot use
Delta_raw.  It must use the quotient observable

    Delta_mir = Delta_raw - E_string_min.

This script checks the proof logic against the existing finite rows.  The
"smoke" row is the decisive counterexample: its raw charged gap is positive
and exactly C3/beta, while the electric-subtracted mirror offset is zero.  Any
continuum theorem that cited the raw gap would incorrectly count a pure
electric string as a mirror mass.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUN_DIR = ROOT / "smg_dmrg_runs"
C3 = 4.0 / 3.0


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def load_json(name: str):
    return json.loads((RUN_DIR / name).read_text(encoding="utf-8"))


def electric_ref(row: dict) -> float:
    return float(row.get("electric_ref", C3 / float(row["beta"])))


def offset(row: dict) -> float:
    return float(row["charged_gap"]) - electric_ref(row)


def raw_gap_can_lie(rows: list[dict]) -> None:
    print("[2] finite counterexample: raw charged gap can be pure string energy")
    max_abs_offset = 0.0
    for row in rows:
        beta = float(row["beta"])
        raw = float(row["charged_gap"])
        string = electric_ref(row)
        off = raw - string
        max_abs_offset = max(max_abs_offset, abs(off))
        print(
            f"  beta={beta:6.3f}: Delta_raw={raw:.12f}, "
            f"E_string=C3/beta={string:.12f}, Delta_mir={off:+.3e}"
        )
        check(raw > 0.0, "raw charged gap is positive")
        check(abs(off) < 1.0e-9, "electric-subtracted mirror offset is zero")
    check(max_abs_offset < 1.0e-9, "smoke row is pure electric string, not mirror mass")


def positive_witness_rows(rows: list[dict], label: str, floor: float) -> None:
    print(f"[3] positive witness rows: {label}")
    offsets = []
    for row in rows:
        beta = float(row["beta"])
        raw = float(row["charged_gap"])
        string = electric_ref(row)
        off = raw - string
        offsets.append(off)
        print(
            f"  beta={beta:6.3f}: Delta_raw={raw:.6f}, "
            f"E_string={string:.6f}, Delta_mir={off:.6f}"
        )
    print(f"  minimum electric-subtracted mirror offset = {min(offsets):.6f}")
    check(min(offsets) > floor, f"{label} has positive mirror offset above floor")


def theorem_text() -> str:
    return r"""
Derivation.

1. Let Q be the external charged/Gauss sector used to test a charged mirror
   excitation.  The transfer-matrix energy in that sector is

       E_full(Q) = E_gauge_string(Q) + E_residual(Q).

   The first term is forced by Gauss law: a non-singlet endpoint must be
   attached to gauge flux or to a Wilson-line dressing.  It is already present
   in the pure gauge theory.

2. Mirror decoupling is a statement about the residual mirror pole.  Hence the
   relevant mass gap is the charged energy above the pure gauge static
   background:

       Delta_mir = [E_full(Q)-E_full(0)] - [E_gauge(Q)-E_gauge(0)].

   On the finite electric rows, [E_gauge(Q)-E_gauge(0)] is the minimal electric
   string energy E_string_min=C3/beta.

3. Any raw-gap criterion is not invariant under changes of the pure gauge
   string energy with the mirror sector held fixed.  In the limiting
   countermodel E_residual(Q)=0 but E_gauge_string(Q)>0, Delta_raw is positive
   although there is no mirror mass.  Therefore Delta_raw is not a valid
   continuum mirror-pole observable.

4. The electric-subtracted gap is the quotient by this pure gauge background.
   A uniform positive lower bound on Delta_mir excludes a charged mirror pole;
   a positive raw Delta_raw alone does not.
"""


def main() -> int:
    print("SMG/TCH ELECTRIC-SUBTRACTED MIRROR-GAP OBSERVABLE THEOREM")
    print("=" * 86)
    print("[1] theorem")
    print(theorem_text().strip())

    smoke = load_json("smg_region2_deep_charged_gap_smoke_2x1.json")
    raw_gap_can_lie(smoke["rows"])

    cell = load_json("smg_region2_beta_sweep_cell336.json")
    positive_witness_rows(cell["rows"], "336-state magnetic cell", floor=3.0)

    two = load_json("smg_region2_2plaq_magnetic_gap.json")[0]
    positive_witness_rows(two["rows"], "2-plaquette magnetic cell", floor=2.0)

    print("[4] decision")
    print(
        "  The continuum theorem needs the electric-subtracted Delta_mir because "
        "Delta_raw is contaminated by the obligatory pure-gauge static string. "
        "The smoke row proves necessity: raw>0 but Delta_mir=0. The magnetic "
        "witness rows are meaningful only after the same subtraction."
    )
    print("ALL ASSERTIONS PASSED -- electric-subtracted mirror gap is the continuum observable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
