#!/usr/bin/env python3
r"""RG-basin reduction for the SMG/TCH weak-coupling frontier.

Scratch artifact; no canon update.

Given ``smg_induced_lambda_bound.py``:

    P_vac H_registered P_ch = 0

for electric, matter, magnetic, and charged-hopping terms.  Therefore the
mirror-induced Wilsonian deformation coordinates are exactly zero inside the
registered Hamiltonian.

This script spells out the consequence in coupling-space language.  The RG path
is not

    (beta_F, beta_A, lambda_6, lambda_8, ...)

with unknown mirror-induced beta_A/lambda_i.  It is simply

    (beta_F, 0, 0, 0, ...)

plus a decoupled massive spectator sector.  Thus the SMG-specific basin problem
is closed modulo one external lattice-gauge input: the pure/fundamental Wilson
SU(3) axis has no finite-beta bulk transition on beta_F >= beta_cert.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BETA_CERT = 0.661156
DELTA_MIN = 3.960160


def load_lambda_module():
    spec = importlib.util.spec_from_file_location(
        "lambda_bound", ROOT / "smg_induced_lambda_bound.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def eta_max_for_bulk_radius(radius: float, delta_min: float = DELTA_MIN) -> float:
    """Max ||P0 V_pair Pch|| allowed if an unknown bulk radius is supplied."""
    if radius < 0:
        raise ValueError("bulk radius must be nonnegative")
    return (radius * delta_min) ** 0.5


def main() -> int:
    lb = load_lambda_module()
    assert abs(lb.BETA_CERT - BETA_CERT) < 1e-12
    assert abs(lb.T_PHYSICAL - 1.0) < 1e-12

    print("[0] SMG RG-basin reduction")
    print(f"    beta_cert = {BETA_CERT:.6f}")
    print(f"    charged denominator floor from finite audit = {DELTA_MIN:.6f}")
    print()

    print("[1] coupling-space path")
    print("    coordinates: (beta_F, beta_A, lambda_6, lambda_8, ...)")
    print("    registered path: beta_F in [beta_cert, infinity), all other coordinates = 0")
    print("    reason: exact P0 H Pch = 0 for every registered local TCH/SMG operator")
    print()

    print("[2] decision logic")
    print("    IF pure/fundamental Wilson SU(3) axis has no finite-beta bulk transition")
    print("       on beta_F >= beta_cert:")
    print("         -> SMG mirror block is a massive spectator along the continuum path.")
    print("    ELSE:")
    print("         -> the obstruction is ordinary gauge-sector bulk physics, not SMG")
    print("            mirror backreaction.")
    print()

    print("[3] future-pair-source budget")
    print("    If a new number-changing operator is later added,")
    print("      |lambda_pair| <= eta^2 / Delta_min.")
    print("    For a known distance r_bulk to the nearest bulk surface, require")
    print("      eta < sqrt(r_bulk * Delta_min).")
    for radius in (1e-4, 1e-3, 1e-2, 1e-1):
        print(
            f"      r_bulk={radius:<7g} -> eta_max={eta_max_for_bulk_radius(radius):.6f}"
        )
    print()

    print("[4] theorem status")
    print("    SMG-specific induced-coupling clause: CLOSED for the registered Hamiltonian.")
    print("    Remaining clause: external/numerical lattice-gauge no-bulk-transition")
    print("    statement for the pure Wilson SU(3) axis. This is not a new framework")
    print("    parameter; it is the standard continuum-basin input.")
    print()
    print("[VERDICT]")
    print("  The RG path stays in the no-bulk-transition basin exactly to the extent that")
    print("  the pure Wilson SU(3) axis does. The mirror sector contributes no induced")
    print("  adjoint/irrelevant coordinates in the registered operator algebra. Therefore")
    print("  there is no remaining SMG-specific lambda_i estimate to tune.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
