#!/usr/bin/env python3
r"""SMG/TCH weak-coupling frontier audit.

The finite TCH gauge-cell results are strong, but the continuum beta extension
is still a frontier.  This script is an observable-hygiene gate: it separates
the raw finite electric string gap from the electric-subtracted mirror/SMG gap,
then states the theorem target that would promote the finite rows to a
weak-coupling continuum result.

Inputs are already-generated lightweight JSON tables:

  * smg_region2_beta_sweep_cell336.json
      336-state non-reduced magnetic cell, sector-resolved.

  * smg_region2_deep_charged_gap_smoke_2x1.json
      43k-state smoke row.  This row is intentionally classified as an
      electric-string sanity check: its gap is C3/beta to numerical precision,
      so it is not the mirror-mass observable.

The 106,460-state extended rows are recorded directly from the committed
scripts tch_2plaq_extended_deep.py and tch_2plaq_extended_hop.py because those
scripts print the table rather than writing JSON.

Exit 0 means:
  1. the best current vacuum local-polymer certificate still stops at beta=0.661;
  2. the smoke row is rejected as a mirror-gap witness;
  3. all finite TCH rows that do include the mirror/SMG offset have positive
     electric-subtracted gaps;
  4. the remaining continuum theorem is named precisely, with the
     electric-subtracted lower-bound problem open from beta=1 upward.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RUN_DIR = ROOT / "smg_dmrg_runs"

C3 = 4.0 / 3.0
PHI = (1.0 + math.sqrt(5.0)) / 2.0
BETA_CERT = 0.661155716805
Z_FP = 0.06630807677


def z_activity(beta: float) -> float:
    return beta * beta * PHI / (8.0 * C3)


def load_json(name: str) -> dict:
    path = RUN_DIR / name
    return json.loads(path.read_text())


def electric_subtracted_rows(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        beta = float(row["beta"])
        raw = float(row["charged_gap"])
        electric = float(row.get("electric_ref", C3 / beta))
        out.append(
            {
                "beta": beta,
                "raw": raw,
                "electric": electric,
                "mirror_offset": raw - electric,
            }
        )
    return out


def main() -> int:
    print("SMG/TCH WEAK-COUPLING FRONTIER AUDIT")
    print()

    print("[0] certified-domain endpoint")
    beta_one_ratio = z_activity(1.0) / Z_FP
    print(f"    FP/independent-neighbourhood beta_cert = {BETA_CERT:.12f}")
    print(f"    z(beta=1)/z_cert = {beta_one_ratio:.3f}")
    assert 0.661 < BETA_CERT < 0.662
    assert 2.28 < beta_one_ratio < 2.30
    print("    beta=1 remains outside the cheap local-polymer certificate.")
    print()

    print("[1] 336-state gauge-cell region-II sweep: subtract the electric string")
    cell = load_json("smg_region2_beta_sweep_cell336.json")
    cell_offsets = electric_subtracted_rows(cell["rows"])
    min_cell = min(cell_offsets, key=lambda r: r["mirror_offset"])
    for row in cell_offsets:
        print(
            f"    beta={row['beta']:5.3f}: raw={row['raw']:8.4f}, "
            f"C3/beta={row['electric']:8.4f}, offset={row['mirror_offset']:8.4f}"
        )
    print(
        f"    finite-cell electric-subtracted floor = {min_cell['mirror_offset']:.4f} "
        f"at beta={min_cell['beta']:.3f}"
    )
    assert min_cell["mirror_offset"] > 3.9
    print("    This is the finite object that actually bears on mirror survival.")
    print()

    print("[2] smoke 2x1 row: classify, do not over-cite")
    smoke = load_json("smg_region2_deep_charged_gap_smoke_2x1.json")
    smoke_offsets = electric_subtracted_rows(smoke["rows"])
    max_smoke_offset = max(abs(row["mirror_offset"]) for row in smoke_offsets)
    for row in smoke_offsets:
        print(
            f"    beta={row['beta']:5.3f}: raw={row['raw']:8.4f}, "
            f"C3/beta={row['electric']:8.4f}, offset={row['mirror_offset']:+.3e}"
        )
    assert max_smoke_offset < 1e-3
    print(
        "    Verdict: this row is an electric finite-volume/string sanity check, "
        "not a mirror-mass witness."
    )
    print()

    print("[3] 106,460-state extended rows: finite-cell stability with the right offset")
    # From tch_2plaq_extended_deep.py and tch_2plaq_extended_hop.py.
    ext_rows = [
        {"label": "ext magnetic t=0", "beta": 0.5, "gap": 6.671182, "t": 0.0},
        {"label": "ext magnetic t=0", "beta": 1.0, "gap": 4.929975, "t": 0.0},
        {"label": "ext hop t=1", "beta": 0.5, "gap": 6.017000, "t": 1.0},
        {"label": "ext hop t=1", "beta": 1.0, "gap": 4.320000, "t": 1.0},
    ]
    min_ext = None
    for row in ext_rows:
        offset = row["gap"] - C3 / row["beta"]
        row["offset"] = offset
        min_ext = row if min_ext is None or offset < min_ext["offset"] else min_ext
        print(
            f"    {row['label']:18s} beta={row['beta']:3.1f}: "
            f"gap={row['gap']:8.4f}, C3/beta={C3/row['beta']:8.4f}, "
            f"offset={offset:8.4f}"
        )
    assert min_ext is not None and min_ext["offset"] > 2.9
    print(
        f"    finite extended-row floor = {min_ext['offset']:.4f} "
        f"({min_ext['label']}, beta={min_ext['beta']:.1f})."
    )
    print()

    print("[4] theorem target")
    print(
        "    Define Delta_mirror(L,cut,beta,t) = "
        "E0(gauge-invariant charged sector) - E0(vacuum) - E_string_min(beta)."
    )
    print(
        "    Finite evidence: Delta_mirror is positive in the 336-state cell "
        "and in the 106,460 extended two-plaquette rows."
    )
    print(
        "    Local-defect lower bound certifies the electric-subtracted offset on "
        "beta_cert <= beta < 1."
    )
    print(
        "    Remaining proof: show inf_{L,cut,beta>=1} Delta_mirror > 0 "
        "along the pure fundamental Wilson SU(3) axis, with no induced adjoint/mirror "
        "couplings and no bulk transition."
    )
    print(
        "    This is the beta-extension frontier.  It is not closed by another "
        "finite row unless that row participates in a cutoff/volume lower-bound "
        "argument or a genuine RG/continuum theorem."
    )

    print()
    print("[VERDICT]")
    print(
        "  The finite-domain TCH/Peter-Weyl/CG construction is strong and internally "
        "consistent.  The weak-coupling/continuum extension is now reduced to a "
        "specific electric-subtracted mirror-gap theorem plus the ordinary pure "
        "Wilson-axis no-bulk-transition/analyticity input.  Do not cite raw "
        "C3/beta finite electric gaps as SMG survival."
    )
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
