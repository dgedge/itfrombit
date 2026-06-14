#!/usr/bin/env python3
r"""Bound mirror-induced local gauge couplings for the SMG/TCH RG route.

This is a scratch/audit artifact, not a canon update.

The continuum-decoupling theorem needs a bound on induced local gauge couplings

    S_eff[U] = S_Wilson[U; beta_R] + sum_i lambda_i O_i[U]

after the massive mirror sector is integrated out.  In the finite TCH operator
that has actually been built, this question has a sharp algebraic answer:

  * the charged/mirror count N_ch is an exact block label;
  * electric, matter, magnetic plaquette, and charged hopping all commute with
    N_ch;
  * therefore P_vac H P_ch = 0 exactly, and Schrieffer-Wolff induced pure-gauge
    lambdas from mirror virtual excursions vanish exactly inside the registered
    Hamiltonian.

This does not prove the full continuum no-bulk-transition theorem.  It says the
mirror block is a spectator unless an additional number-changing/pair-creation
operator is introduced.  For such an extra perturbation V_pair with
||P_vac V_pair P_ch|| <= eta, the conservative second-order bound is

    |lambda_ind| <= eta^2 / Delta_ch,

where Delta_ch is the charged-sector gap.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parent
BETA_CERT = 0.661156
T_PHYSICAL = 1.0


def load_magnetic_module():
    spec = importlib.util.spec_from_file_location(
        "mag", ROOT / "tch_nonreduced_magnetic_plaquette.py"
    )
    module = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        assert spec.loader is not None
        spec.loader.exec_module(module)
    return module


def sector_indices(basis):
    counts = np.array([sum(1 for f in fs if f != "1") for (_rs, fs, _ms) in basis])
    return counts, {q: np.where(counts == q)[0] for q in sorted(set(counts))}


def commutator_norm(counts, op):
    d = counts.astype(float)
    return float(np.linalg.norm((d[:, None] - d[None, :]) * op))


def offblock_norm(op, src, dst):
    if len(src) == 0 or len(dst) == 0:
        return 0.0
    return float(np.linalg.norm(op[np.ix_(src, dst)], 2))


def block_ground(op, idx):
    return float(np.linalg.eigvalsh(op[np.ix_(idx, idx)])[0])


def hamiltonian(module, beta, hopping):
    wsym = module.W + module.W.conj().T
    hsym = module.HOP + module.HOP.conj().T
    return module.ELE / beta + module.MAT - 0.5 * beta * wsym + hopping * hsym


def charged_gap(module, sectors, beta, hopping):
    h = hamiltonian(module, beta, hopping)
    e0 = block_ground(h, sectors[0])
    gaps = []
    for q, idx in sectors.items():
        if q == 0:
            continue
        gaps.append((q, block_ground(h, idx) - e0))
    return min(gaps, key=lambda row: row[1])


def minimized_charged_gap(module, sectors, hopping, beta_bounds=(BETA_CERT, 50.0)):
    lo, hi = map(math.log, beta_bounds)

    def objective(log_beta):
        _q, gap = charged_gap(module, sectors, math.exp(log_beta), hopping)
        return gap

    grid = np.linspace(lo, hi, 121)
    vals = np.array([objective(x) for x in grid])
    i = int(np.argmin(vals))
    left = grid[max(0, i - 2)]
    right = grid[min(len(grid) - 1, i + 2)]
    res = minimize_scalar(objective, bounds=(left, right), method="bounded")
    beta = math.exp(float(res.x))
    q, gap = charged_gap(module, sectors, beta, hopping)
    return gap, beta, q


def main() -> int:
    mag = load_magnetic_module()
    counts, sectors = sector_indices(mag.basis)
    q0 = sectors[0]
    qcharged = np.concatenate([idx for q, idx in sectors.items() if q > 0])

    ops = {
        "electric": mag.ELE,
        "matter": mag.MAT,
        "magnetic W+Wdag": mag.W + mag.W.conj().T,
        "hopping H+Hdag": mag.HOP + mag.HOP.conj().T,
    }

    print("[0] SMG induced-lambda bound")
    print(f"    beta_cert={BETA_CERT:.6f}, physical hopping t={T_PHYSICAL:.1f}")
    print(f"    sectors: " + ", ".join(f"q={q}: {len(idx)}" for q, idx in sectors.items()))
    print()

    print("[1] exact block-coupling audit")
    max_off = 0.0
    for name, op in ops.items():
        comm = commutator_norm(counts, op)
        off = offblock_norm(op, q0, qcharged)
        max_off = max(max_off, off)
        print(f"    {name:<18s} ||[N_ch, O]||={comm:.3e}   ||P0 O Pch||={off:.3e}")
        assert comm < 1e-10
        assert off < 1e-10
    print("    result: registered local operators induce no vacuum-to-mirror excursions.")
    print()

    print("[2] charged-sector denominators")
    q_cert, gap_cert = charged_gap(mag, sectors, BETA_CERT, T_PHYSICAL)
    gap_min, beta_min, q_min = minimized_charged_gap(mag, sectors, T_PHYSICAL)
    print(f"    gap at beta_cert: Delta_ch={gap_cert:.6f} in q={q_cert}")
    print(
        "    min over beta in [beta_cert, 50]: "
        f"Delta_ch={gap_min:.6f} at beta={beta_min:.6f}, q={q_min}"
    )
    assert gap_cert > 0.0
    assert gap_min > 0.0
    print()

    print("[3] induced lambda bound")
    print("    registered Hamiltonian:")
    print(f"      lambda_induced <= ||P0 H Pch||^2 / Delta_ch = {max_off**2 / gap_min:.3e}")
    print("      algebraic reading: lambda_induced = 0 exactly within this operator set.")
    print("    if a future pair-creation perturbation has ||P0 V_pair Pch|| <= eta:")
    for eta in (1e-3, 1e-2, 1e-1, 1.0):
        print(f"      eta={eta:<7g} -> |lambda_pair| <= {eta * eta / gap_min:.6e}")
    print()

    print("[4] basin reading")
    print("    Since lambda_induced=0 for the registered number-conserving mirror block,")
    print("    the RG path is the pure Wilson/fundamental TCH gauge path plus a decoupled")
    print("    massive spectator sector. The finite-cell data do not have to push an")
    print("    adjoint or irrelevant coupling through a bulk-transition surface; there is")
    print("    no induced coupling to push. The remaining non-finite input is therefore")
    print("    the standard gauge-sector statement: the fundamental Wilson SU(3) path from")
    print("    beta_cert to beta=infinity has no finite-beta bulk transition.")
    print()
    print("[VERDICT]")
    print("  Within the registered TCH/SMG Hamiltonian, mirror-induced local gauge")
    print("  lambdas are exactly zero by charged-number block diagonalization. A nonzero")
    print("  lambda requires adding a number-changing pair source not present in the")
    print("  current operator set; if added with strength eta, the conservative bound is")
    print("  eta^2/Delta_ch with Delta_ch >= the value above. Thus the RG/basin problem")
    print("  reduces to the ordinary pure-gauge no-bulk-transition input, not a mirror")
    print("  backreaction problem.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
