#!/usr/bin/env python3
r"""Audit the black-hole 10/27 transfer premise.

Question checked
----------------
Does V_cell/V_Sch itself force the 27-slot Moore service alphabet, or is the
10/27 Hawking attempt-rate coefficient conditional on transferring item-120's
Landauer-Moore alphabet to horizon severing?

Result
------
V_cell/V_Sch does not force the Moore alphabet.  It fixes:

  * the Q3 coboundary syndrome record;
  * the one global-complement/vacuum latch;
  * orthogonal Schwarzschild shell-cell addresses in V_Sch.

It does not contain a spatial service-neighbour label.  Any orthogonal spatial
alphabet A can be tensored onto the same isometry without affecting
injectivity.  Therefore the 26-slot Moore shell is an extra local-service
principle, not a consequence of V_cell alone.

The transfer theorem is still strong: item 120 already derives the monitored
Landauer-Moore alphabet for an isotropic 3D local erasure event.  If horizon
severing is the same local Landauer service class, then the alphabet is the
26-slot Moore shell plus the V_cell latch, and the Schwarzschild normal selects
the outward 3x3 face plus the latch:

    (9 + 1) / (26 + 1) = 10/27.

This script is deliberately a boundary audit: it proves what is and is not
carried by V_cell/V_Sch.
"""

from __future__ import annotations

from collections import Counter
import itertools
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANCHOR = (ROOT / "ANCHOR.md").read_text(encoding="utf-8")

ALPHA0 = 1.0 / 137.0
PHI = (math.sqrt(5.0) - 1.0) / 2.0
EPS_F = 1.0 / (2.0 * PHI)
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
ALL = (1 << 8) - 1
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def syndrome_word(s: int) -> int:
    out = 0
    for n, (i, j) in enumerate(EDGES):
        out |= (bit(s, i) ^ bit(s, j)) << n
    return out


def complement_representatives() -> list[int]:
    return [s for s in range(256) if s < (s ^ ALL)]


def v_cell_rows() -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = []
    for rep in complement_representatives():
        for gamma in (0, 1):
            s = rep ^ (ALL if gamma else 0)
            rows.append((syndrome_word(s), gamma))
    return rows


def source_ladder(beta_eff: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta_eff * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return dict(sorted((f, v / z) for f, v in raw.items()))


def mean_f(beta_eff: float) -> float:
    return sum(f * p for f, p in source_ladder(beta_eff).items())


def required_gamma_for_stefan(beta_eff: float = 1.0) -> float:
    rho_lambda = beta_eff / (2.0 * math.pi * EPS_F)
    target_c = 1.0 / (15360.0 * math.pi)
    denom = (27.0 * math.pi / 32.0) * EPS_F * mean_f(beta_eff) * rho_lambda**4
    return target_c / denom


def model_over_target(prefactor: float) -> float:
    return prefactor * ALPHA0 / required_gamma_for_stefan()


def moore_shell() -> list[tuple[int, int, int]]:
    return [v for v in itertools.product((-1, 0, 1), repeat=3) if v != (0, 0, 0)]


def main() -> None:
    print("BLACK-HOLE 10/27 TRANSFER-PREMISE AUDIT")
    print("=" * 86)

    print("[1] What V_cell actually carries")
    rows = v_cell_rows()
    counts = Counter(rows)
    check(len(rows) == 256, "V_cell has 128 complement classes times 2 latch values")
    check(len(counts) == 256, "syndrome+latch rows are injective")
    check(all(isinstance(row[0], int) and row[1] in (0, 1) for row in rows), "rows are syndrome plus latch only")
    print("    V_cell row type: (12-edge syndrome word, latch bit)")
    print("    No row contains a spatial neighbour / Moore offset label.")

    print("\n[2] Underdetermination: arbitrary spatial alphabets tensor without changing V")
    # If a spatial service label is an orthogonal output factor, V_cell stays an
    # isometry for any alphabet size.  Hence V_cell cannot select 26 or 27.
    for n_alpha in (1, 6, 8, 12, 26, 27, 45):
        lifted = [(row, a) for row in rows for a in range(n_alpha)]
        ok = len(lifted) == len(set(lifted))
        print(f"    alphabet size {n_alpha:2d}: lifted rows={len(lifted):5d}, injective={ok}")
        check(ok, f"V_cell remains isometric with orthogonal alphabet size {n_alpha}")
    print("    Therefore the Moore count is not a theorem of V_cell/V_Sch alone.")

    print("\n[3] What the Landauer-Moore transfer adds")
    shell = moore_shell()
    outward = [v for v in shell if v[2] == 1]
    inward = [v for v in shell if v[2] == -1]
    tangent = [v for v in shell if v[2] == 0]
    check("uniform $1/26$" in ANCHOR and "monitored-channel theorem" in ANCHOR, "item 120 derives the monitored uniform Moore-shell measure")
    check(len(shell) == 26, "3D Moore shell has 26 spatial service slots")
    check(len(outward) == 9 and len(inward) == 9 and len(tangent) == 8, "a local horizon normal splits Moore shell as 9 outward / 9 inward / 8 tangent")
    pref = (len(outward) + 1) / (len(shell) + 1)
    print(f"    transferred alphabet: 26 Moore slots + 1 latch = {len(shell) + 1}")
    print(f"    emitted subset:       9 outward slots + 1 latch = {len(outward) + 1}")
    print(f"    prefactor:            {pref:.12f} = 10/27")
    print(f"    P/P_SB:               {model_over_target(pref):.9f}")
    check(abs(pref - 10.0 / 27.0) < 1e-15, "Landauer-Moore transfer gives exactly 10/27")
    check(abs(model_over_target(pref) - 1.0) < 4e-3, "10/27 alpha lands within 0.4 percent of Stefan-Hawking coefficient")

    print("\n[4] Boundary verdict")
    print(
        """
    Confirmed:
      V_cell/V_Sch alone fixes the latch and the shell-cell direct sum, but not
      the spatial service alphabet.  It is compatible with many orthogonal
      service alphabets, so it cannot by itself prove 26+1 slots.

    Conditional closure:
      If horizon severing is an isotropic local Landauer service event of the
      same class as item 120, the item-120 Moore-shell theorem supplies the
      missing alphabet.  Then the horizon normal plus partner asymmetry selects
      outward face plus latch, and the 10/27 scheduler follows.

    If that cross-sector transfer is rejected, the 10/27 coefficient remains a
    strong stencil candidate rather than a locked horizon theorem.
ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
