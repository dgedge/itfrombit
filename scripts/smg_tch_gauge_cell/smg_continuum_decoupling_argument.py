#!/usr/bin/env python3
r"""SMG continuum decoupling / no-phase-transition argument audit.

This is deliberately NOT a new ANCHOR development.  It tests whether the
remaining beta > 0.661 frontier can be advanced by a continuum-decoupling
argument rather than another strong-coupling polymer bound.

The target shape is:

  1. finite algebra: the mirror/charged-count sectors are superselected;
  2. finite evidence: no sector crossing appears in the non-reduced cell through
     the weak-coupling side for physical/moderate hopping;
  3. continuum premise: if the TCH SU(3) gauge path has no bulk transition and
     the gauge-invariant charged sector stays confining/positive, the mirror
     gap adiabatically continues past the polymer certificate.

The script proves (1), measures (2), and refuses to promote (3).  A no-transition
argument closes only conditional on that external RG/continuum premise.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parent


def load_magnetic_module():
    """Load the validated 336-state non-reduced magnetic cell quietly."""
    spec = importlib.util.spec_from_file_location(
        "mag", ROOT / "tch_nonreduced_magnetic_plaquette.py"
    )
    module = importlib.util.module_from_spec(spec)
    with contextlib.redirect_stdout(io.StringIO()):
        spec.loader.exec_module(module)
    return module


def sector_indices(basis):
    counts = np.array([sum(1 for f in Fs if f != "1") for (_Rs, Fs, _ms) in basis])
    return counts, {q: np.where(counts == q)[0] for q in sorted(set(counts))}


def commutator_norm(diag_values, op):
    d = diag_values.astype(float)
    return float(np.linalg.norm((d[:, None] - d[None, :]) * op))


def sector_ground(H, idx):
    if len(idx) == 0:
        return None
    return float(np.linalg.eigvalsh(H[np.ix_(idx, idx)])[0])


def scan_sector_gaps(module, counts, sectors, betas, hoppings):
    Wsym = module.W + module.W.conj().T
    Hsym = module.HOP + module.HOP.conj().T
    rows = []
    for t in hoppings:
        best = (float("inf"), None, None, None)
        for beta in betas:
            H = module.ELE / beta + module.MAT + t * Hsym - (beta / 2.0) * Wsym
            e_vac = sector_ground(H, sectors[0])
            charged = []
            for q in sorted(q for q in sectors if q > 0):
                e_q = sector_ground(H, sectors[q])
                if e_q is not None:
                    charged.append((q, e_q - e_vac))
            q_min, gap_min = min(charged, key=lambda x: x[1])
            if gap_min < best[0]:
                best = (gap_min, beta, q_min, charged)
        rows.append((t, *best))
    return rows


def block_operators(module, sectors):
    Wsym = module.W + module.W.conj().T
    Hsym = module.HOP + module.HOP.conj().T
    blocks = {}
    for q, idx in sectors.items():
        blocks[q] = {
            "A": module.ELE[np.ix_(idx, idx)],
            "B": module.MAT[np.ix_(idx, idx)],
            "C": -0.5 * Wsym[np.ix_(idx, idx)],
            "D": Hsym[np.ix_(idx, idx)],
        }
    return blocks


def ground_from_blocks(blocks, q, beta, hopping):
    blk = blocks[q]
    H = blk["A"] / beta + blk["B"] + beta * blk["C"] + hopping * blk["D"]
    return float(np.linalg.eigvalsh(H)[0])


def charged_gap_continuous(blocks, beta, hopping):
    e_vac = ground_from_blocks(blocks, 0, beta, hopping)
    charged = []
    for q in sorted(q for q in blocks if q > 0):
        charged.append((q, ground_from_blocks(blocks, q, beta, hopping) - e_vac))
    return min(charged, key=lambda x: x[1])


def minimized_gap(blocks, hopping, beta_bounds=(0.05, 50.0)):
    lo, hi = map(math.log, beta_bounds)

    def objective(log_beta):
        _q, gap = charged_gap_continuous(blocks, math.exp(log_beta), hopping)
        return gap

    # Use a coarse pass to avoid a local minimum of the lowest sector envelope.
    grid = np.linspace(lo, hi, 121)
    vals = np.array([objective(x) for x in grid])
    best_i = int(np.argmin(vals))
    left = grid[max(0, best_i - 2)]
    right = grid[min(len(grid) - 1, best_i + 2)]
    res = minimize_scalar(objective, bounds=(left, right), method="bounded", options={"xatol": 1e-10})
    beta = math.exp(float(res.x))
    q, gap = charged_gap_continuous(blocks, beta, hopping)
    return gap, beta, q


def main() -> int:
    mag = load_magnetic_module()
    counts, sectors = sector_indices(mag.basis)
    n_op = np.diag(counts)

    print("[0] SMG continuum/no-transition audit")
    print("    scope: 336-state non-reduced square-cell Hamiltonian; proof clauses only")
    print("    no ANCHOR update; no smg_dmrg execution")
    print()

    print("[1] finite algebra: charged-count superselection")
    ops = {
        "electric": mag.ELE,
        "matter": mag.MAT,
        "magnetic W+Wdag": mag.W + mag.W.conj().T,
        "hopping H+Hdag": mag.HOP + mag.HOP.conj().T,
    }
    for name, op in ops.items():
        norm = commutator_norm(counts, op)
        print(f"    ||[N_charged, {name}]|| = {norm:.3e}")
        assert norm < 1e-10
    print("    result: vacuum and charged/mirror sectors are exact blocks of the finite Hamiltonian.")
    print()

    print("[2] finite non-reduced scan across and beyond the certified domain")
    betas = np.concatenate(
        [
            np.linspace(0.25, 1.25, 41),
            np.linspace(1.30, 3.00, 18),
            np.linspace(3.25, 10.0, 28),
        ]
    )
    rows = scan_sector_gaps(mag, counts, sectors, betas, hoppings=(0.2, 1.0, 4.0))
    print("    t       min charged-sector gap     at beta    sector    reading")
    for t, gap, beta, q, charged in rows:
        if gap > 0:
            reading = "no finite crossing in scan"
        else:
            reading = "sector crossing: stress-only warning"
        print(f"    {t:<5.1f}   {gap:>+12.6f}              {beta:<7.3f}   q={q:<1d}      {reading}")
    assert rows[0][1] > 0.0
    assert rows[1][1] > 0.0
    assert rows[2][1] < 0.0
    print("    interpretation: physical/moderate hopping (t<=1) is adiabatically smooth in")
    print("    this cell through beta=10; the unphysical t=4 stress crossing proves that")
    print("    an unconditional all-hopping theorem would be false.")
    print()

    print("[3] continuous finite-cell envelope")
    blocks = block_operators(mag, sectors)
    dense_rows = []
    print("    t       min charged-sector gap     beta*      sector")
    for t in (0.2, 1.0, 2.0, 3.0, 4.0):
        gap, beta, q = minimized_gap(blocks, t)
        dense_rows.append((t, gap, beta, q))
        print(f"    {t:<5.1f}   {gap:>+12.6f}              {beta:<8.4f}  q={q:<1d}")

    def threshold_objective(t):
        gap, _beta, _q = minimized_gap(blocks, t)
        return gap

    gap_t1 = threshold_objective(1.0)
    gap_t4 = threshold_objective(4.0)
    assert gap_t1 > 0.0
    assert gap_t4 < 0.0
    t_crit = float(brentq(threshold_objective, 1.0, 4.0, xtol=1e-9, rtol=1e-9))
    gap_crit, beta_crit, q_crit = minimized_gap(blocks, t_crit)
    print()
    print(
        "    finite-cell stress threshold: "
        f"t_c = {t_crit:.6f} at beta* = {beta_crit:.6f}, q={q_crit}, "
        f"gap={gap_crit:+.3e}"
    )
    assert dense_rows[0][1] > 0.0
    assert dense_rows[1][1] > 0.0
    assert dense_rows[-1][1] < 0.0
    assert 1.0 < t_crit < 4.0
    print("    interpretation: the finite cell has a measured stress boundary. The")
    print("    physical/moderate t<=1 region is separated from the first finite-cell")
    print("    sector crossing by an O(1) hopping margin, but no all-hopping theorem exists.")
    print()

    print("[4] continuum decoupling clauses")
    print("    C1 finite block structure: PROVED above.")
    print("    C2 finite-row support: positive for t<=1 in the non-reduced cell; existing")
    print("       106,460 rows independently keep beta=1 positive under the registered t-scan.")
    print("    C3 continuum/no-bulk-transition premise: NOT proved here. Required statement:")
    print("       the TCH SU(3) gauge path from beta≈0.661 to weak coupling has no bulk")
    print("       transition in the gauge-invariant sector, and the charged meson/string")
    print("       sector remains positive/confining.")
    print("    C4 EFT decoupling: conditional on C3, the massive CSS-SMG mirror block is")
    print("       a gapped gauge-invariant sector; integrating it out can renormalise local")
    print("       gauge terms but cannot close the mirror gap without a sector crossing.")
    print()
    print("[VERDICT]")
    print("  A continuum decoupling/no-phase-transition route is coherent and sharper than")
    print("  further local polymer tinkering, but it is conditional.  The finite algebra")
    print("  supplies exact sector decoupling, finite rows show no crossing for the")
    print(f"  physical/moderate hopping range, and the stress crossing is at t_c={t_crit:.6f};")
    print("  however, the load-bearing continuum premise is external/RG-level.")
    print("  This advances the frontier to a named theorem:")
    print("    prove no gauge-invariant charged-sector transition from beta=0.661 to the")
    print("    continuum TCH SU(3) limit, or derive the corresponding RG flow.")
    print("  Until that theorem exists, this is not a closure.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
