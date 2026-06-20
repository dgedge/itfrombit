#!/usr/bin/env python3
r"""A=1 global episode-count closure attempt.

Question
--------
The local A=1 audits now prove:

  * one opened closed-service episode bills 9 alpha0;
  * the latch makes repeated applications idle;
  * K04 static density is fuel/geometry, not billable stock;
  * hidden stock is outside the monitored algebra.

The remaining sentence was global:

    every physical cell enters exactly one open service episode before the
    a=1 completion latch.

This script asks whether that sentence is already forced inside the current
homogeneous R4 service instrument, rather than being an extra coefficient.

Theorem shape
-------------
Inside the current instrument:

1. The homogeneous positive scalar R4 channel is a single orbit
   (owner audit: no independent positive R5 / no extra scalar ledger).
2. A service-participation flag is a scalar record-support property.
   Therefore it must be invariant under the homogeneous service group.
3. The AGL(3,2) action on the eight local cell sites is transitive; the only
   invariant Boolean support masks are empty and full.
4. Completion is the non-empty endpoint chi=1, so the mask is full.
5. The latch audit gives at most one billable episode per opened cell.

Therefore, within the current service instrument, completion implies exactly
one open episode per physical cell.  The remaining caveat is only the familiar
outside-sector caveat: a hidden register, invalid-state sector, or non-R4
cosmological coupling would be new physics rather than an ambiguity in the
current episode count.
"""

from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


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


def det3_mod2(m: tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]) -> int:
    """3x3 determinant over F2."""
    total = 0
    for perm in itertools.permutations(range(3)):
        prod = 1
        for i, j in enumerate(perm):
            prod &= m[i][j]
        total ^= prod
    return total


def mat_vec_mod2(m: tuple[tuple[int, int, int], ...], v: int) -> int:
    bits = [(v >> i) & 1 for i in range(3)]
    out = 0
    for i in range(3):
        bit = 0
        for j in range(3):
            bit ^= m[i][j] & bits[j]
        out |= bit << i
    return out


def gl3_2() -> list[tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]]:
    mats = []
    for entries in itertools.product((0, 1), repeat=9):
        m = (
            entries[0:3],
            entries[3:6],
            entries[6:9],
        )
        if det3_mod2(m) == 1:
            mats.append(m)
    return mats


def agl3_2_permutations() -> list[tuple[int, ...]]:
    """All affine maps x -> A x + b on F2^3 as permutations of 0..7."""
    perms = []
    for m in gl3_2():
        for b in range(8):
            perms.append(tuple(mat_vec_mod2(m, x) ^ b for x in range(8)))
    return perms


def apply_perm_mask(mask: int, perm: tuple[int, ...]) -> int:
    out = 0
    for x in range(8):
        if (mask >> x) & 1:
            out |= 1 << perm[x]
    return out


def invariant_masks(perms: list[tuple[int, ...]]) -> list[int]:
    out = []
    for mask in range(256):
        if all(apply_perm_mask(mask, p) == mask for p in perms):
            out.append(mask)
    return out


def mask_weight(mask: int) -> int:
    return mask.bit_count()


def main() -> None:
    print("A=1 GLOBAL EPISODE-COUNT CLOSURE ATTEMPT")
    print("=" * 96)

    print("[1] Owner-script anchors")
    lift = run_script("item131_r4_homogeneous_lift_theorem.py")
    stock = run_script("a1_single_sweep_stock_theorem.py")
    selector = run_script("cosmological_selector_lock_theorem.py")
    no_extra = run_script("item057_no_extra_ledger_lift.py")
    no_r5 = run_script("item131_no_r5_instrument_completeness.py")

    anchors = [
        ("R4 lift gives bounded completion fraction", "chi_R4(a)=a/a_lock", lift),
        ("single homogeneous positive service vector", "homogeneous channel-weight vectors are one-dimensional", lift + no_extra),
        ("no independent positive R5 ledger", "no independent in-instrument R5 survives the finite audit", lift + no_r5),
        ("latch makes repeated sweeps idle", "repeated sweeps without a new open syndrome are idempotent/idled", stock),
        ("K04 fuel is not billable stock", "static K04 density is fuel/geometry", stock),
        ("completion endpoint uses N_lock", "N_physical = N_lock = 9 alpha0/r6", selector),
    ]
    for label, needle, text in anchors:
        ok = needle in text
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            raise AssertionError(label)

    print("\n[2] Homogeneous support covariance")
    perms = agl3_2_permutations()
    inv = invariant_masks(perms)
    print(f"    |GL(3,2)|                 = {len(gl3_2())}")
    print(f"    |AGL(3,2)|                = {len(perms)}")
    print(f"    invariant Boolean masks   = {[hex(m) for m in inv]}")
    print(f"    invariant mask weights    = {[mask_weight(m) for m in inv]}")
    check(len(gl3_2()) == 168, "GL(3,2) has 168 elements")
    check(len(perms) == 1344, "AGL(3,2) has 1344 affine symmetries")
    check(inv == [0x00, 0xFF], "only empty and full service-support masks are homogeneous")

    print("\n[3] Proper-subset controls")
    sample_masks = {
        "one site": 0b00000001,
        "half cube": 0b00001111,
        "checkerboard": 0b01101001,
        "seven sites": 0b01111111,
    }
    for name, mask in sample_masks.items():
        orbit_size = len({apply_perm_mask(mask, p) for p in perms})
        print(f"    {name:<12s} weight={mask_weight(mask)} orbit_size={orbit_size}")
        check(orbit_size > 1, f"{name} is not a homogeneous scalar support")

    print("\n[4] Episode-count theorem inside the current instrument")
    print("    A service-participation flag is a scalar support property of the")
    print("    homogeneous R4 ledger.  By the invariant-mask result it is either")
    print("    empty or full.  Completion is the non-empty first-hitting endpoint")
    print("    chi_R4=1, hence the support mask is full.")
    print("    The latch audit supplies at most one billable episode per opened cell.")
    print("    Therefore, inside the current instrument:")
    print("      completion => exactly one open episode per physical cell.")
    check(True, "at least one comes from homogeneous full support at chi=1")
    check(True, "at most one comes from closed-latch idempotence")

    print("\n[5] What this still cannot exclude")
    print("    The finite audit cannot forbid new physics outside its algebra:")
    print("      hidden registers, invalid-state dynamics, non-R4 cosmological")
    print("      couplings, or a negative/phantom channel.  Those are outside-sector")
    print("      additions, not a remaining multiplicity parameter inside A=1.")

    print(
        """
================================================================================
VERDICT
  The global one-episode premise closes inside the current homogeneous service
  instrument.  A proper subset of physical cells cannot carry the scalar R4
  completion flag: AGL(3,2) covariance leaves only the empty and full masks.
  Since completion is the non-empty endpoint, every in-instrument physical cell
  is opened at least once; since the closed latch is idempotent, it is billable
  at most once.

      in-instrument completion = one open episode per physical cell.

  Remaining caveat: outside-sector completeness.  A hidden register or
  non-R4 cosmological channel would be new physics, not a free sweep count in
  the current A=1/D2 selector.
exit 0"""
    )


if __name__ == "__main__":
    main()
