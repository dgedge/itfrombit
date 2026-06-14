#!/usr/bin/env python3
r"""Item 115: can the existing Q3/Final-Boss 3D Bloch core pin the matter loop?

The withdrawn item-115 Chadha-Nielsen mechanism was the wrong diagram: the
anisotropic T_1u/E_g object is neutral, so the only viable correction to the
gauge-boson velocity is the charged-matter vacuum-polarization bubble

    Pi_mn(q) = Tr[V_m G_F V_n G_F].

Before such a bubble is meaningful, the candidate matter Green function G_F must
come from a charge-conserving, three-dimensional, physical matter Hamiltonian.
This script tests the only currently explicit 3D Q3/Bloch core in canon
(`finalboss_3dbloch.py`) against that precondition.

Verdict target:
  * If the core preserves the 48 physical codewords, commutes with the item-116
    electric charge, and retains rank-3 spatial support, it can seed the loop.
  * If not, the Q3 matter loop is still unpinned; computing a number from this
    kernel would be a gauge-noninvariant proxy, not vacuum polarization.

This is a pre-loop audit, not the loop itself.
"""

from __future__ import annotations

import itertools as it
from fractions import Fraction

import numpy as np


BIT = {
    0: "G0",
    1: "G1",
    2: "LQ",
    3: "C0",
    4: "C1",
    5: "I3",
    6: "chi",
    7: "W",
}
REF_BITS = {3, 4, 5}  # finalboss_3dbloch.py reference hop: C0,C1,I3.


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def octant_perm(axis_perm: tuple[int, int, int], signs: tuple[int, int, int]) -> tuple[int, ...]:
    sigma = [0] * 8
    for v in range(8):
        c = [((v >> a) & 1) * 2 - 1 for a in range(3)]
        cp = [signs[j] * c[axis_perm[j]] for j in range(3)]
        sigma[v] = sum(((cp[j] + 1) // 2) << j for j in range(3))
    return tuple(sigma)


def det(axis_perm: tuple[int, int, int], signs: tuple[int, int, int]) -> int:
    p = np.zeros((3, 3))
    for j in range(3):
        p[j, axis_perm[j]] = signs[j]
    return round(np.linalg.det(p))


def proper_octahedral_rotations() -> list[tuple[int, ...]]:
    rots = []
    for axis_perm in it.permutations(range(3)):
        for signs in it.product((1, -1), repeat=3):
            if det(axis_perm, signs) == 1:
                rots.append(octant_perm(axis_perm, signs))
    check(len(rots) == 24, "proper octahedral subgroup has 24 rotations")
    return rots


def bits_to_mask(bits: set[int]) -> int:
    return sum(1 << b for b in bits)


def build_finalboss_hops() -> tuple[dict[int, int], dict[int, np.ndarray]]:
    rots = proper_octahedral_rotations()
    flips: dict[int, int] = {}
    dirs: dict[int, np.ndarray] = {}
    for f in range(8):
        sigma = next(
            sg
            for sg in rots
            if sg[0] == f and 0 not in {sg[b] for b in REF_BITS}
        )
        flips[f] = bits_to_mask({sigma[b] for b in REF_BITS})
        dirs[f] = np.array(
            [
                2 * ((f >> 2) & 1) - 1,
                2 * ((f >> 1) & 1) - 1,
                2 * (f & 1) - 1,
            ],
            dtype=float,
        )
    return flips, dirs


def is_codeword(n: int) -> bool:
    b = lambda i: (n >> i) & 1
    if b(0) and b(1):
        return False
    if b(7) != b(6):
        return False
    colour = (b(3), b(4))
    if b(2) == 0 and colour != (0, 0):
        return False
    if b(2) == 1 and colour == (0, 0):
        return False
    return True


def charge(n: int) -> Fraction:
    """Item-116/2.8 electric charge on physical matter codewords."""
    b = lambda i: (n >> i) & 1
    z_f = 1 if b(5) == 0 else -1
    sum_z_c = -3 if (b(3), b(4)) == (0, 0) else -1
    return Fraction(1, 2) * z_f + Fraction(1, 3) * sum_z_c + Fraction(1, 2)


def p_block(mask: int, physical: list[int], index: dict[int, int], keep_charge: bool) -> np.ndarray:
    out = np.zeros((len(physical), len(physical)))
    for col, n in enumerate(physical):
        m = n ^ mask
        if m not in index:
            continue
        if keep_charge and charge(m) != charge(n):
            continue
        out[index[m], col] = 1.0
    return out


def h_pp(
    k: np.ndarray,
    flips: dict[int, int],
    dirs: dict[int, np.ndarray],
    physical: list[int],
    index: dict[int, int],
    keep_charge: bool,
) -> np.ndarray:
    h = np.zeros((len(physical), len(physical)))
    for f in range(8):
        h += p_block(flips[f], physical, index, keep_charge) * np.cos(float(np.dot(k, dirs[f])))
    return h


def main() -> None:
    flips, dirs = build_finalboss_hops()
    physical = [n for n in range(256) if is_codeword(n)]
    invalid = [n for n in range(256) if not is_codeword(n)]
    index = {n: i for i, n in enumerate(physical)}
    q_diag = np.array([float(charge(n)) for n in physical])
    q_op = np.diag(q_diag)

    print("\n[1] Physical matter space and item-116 charge")
    check(len(physical) == 48, "R1/R2/R3 physical matter space has 48 codewords")
    check(len(invalid) == 208, "invalid complement has 208 codewords")
    charges = sorted(set(charge(n) for n in physical))
    charge_counts = {str(q): sum(charge(n) == q for n in physical) for q in charges}
    print(f"  physical charges: {charge_counts}")
    check(charge_counts == {"-1": 6, "-1/3": 18, "0": 6, "2/3": 18}, "charge multiplicities match item 116")

    print("\n[2] Candidate 3D Bloch core versus physical closure and U(1) charge")
    rows = []
    total_pp = total_pq = total_charge_ok = total_charge_bad = 0
    for f in range(8):
        mask = flips[f]
        pp = 0
        pq = 0
        charge_ok = 0
        charge_bad = 0
        for n in physical:
            m = n ^ mask
            if m in index:
                pp += 1
                if charge(m) == charge(n):
                    charge_ok += 1
                else:
                    charge_bad += 1
            else:
                pq += 1
        total_pp += pp
        total_pq += pq
        total_charge_ok += charge_ok
        total_charge_bad += charge_bad
        flip_names = ",".join(BIT[b] for b in range(8) if (mask >> b) & 1)
        rows.append((f, tuple(int(x) for x in dirs[f]), flip_names, pp, pq, charge_ok, charge_bad))

    print("  f   dir              flips          P->P  P->Q  charge-ok  charge-bad")
    for f, direction, flip_names, pp, pq, charge_ok, charge_bad in rows:
        print(f"  {f:<2}  {str(direction):<15} {flip_names:<13} {pp:>4} {pq:>5} {charge_ok:>10} {charge_bad:>11}")
    print(f"  totals: P->P={total_pp}, P->Q={total_pq}, charge-ok={total_charge_ok}, charge-bad={total_charge_bad}")
    check(total_pq > 0, "the Final-Boss core is not closed on the 48 physical matter states")
    check(total_charge_bad > 0, "the physical P-block contains electric-charge-changing hops")
    check(total_charge_ok == 32, "only 32 of 384 one-hop physical attempts are charge-preserving")

    sample_k = np.array([0.37, -0.21, 0.44])
    h_phys = h_pp(sample_k, flips, dirs, physical, index, keep_charge=False)
    comm = h_phys @ q_op - q_op @ h_phys
    comm_frob = float(np.linalg.norm(comm))
    comm_max = float(np.max(np.abs(comm)))
    print("\n[3] Charge-commutator gate")
    print(f"  sample k = {sample_k.tolist()}")
    print(f"  ||[P H(k) P, Q]||_F = {comm_frob:.6f}")
    print(f"  max |[P H(k) P, Q]_ij| = {comm_max:.6f}")
    check(comm_frob > 1.0 and comm_max > 0.1, "candidate physical Hamiltonian violates U(1) charge conservation")

    print("\n[4] Charge-preserving projection: does it leave a 3D matter kernel?")
    charge_projected_basis = []
    active = []
    for f in range(8):
        block = p_block(flips[f], physical, index, keep_charge=True)
        if np.count_nonzero(block):
            active.append((f, tuple(int(x) for x in dirs[f]), int(np.count_nonzero(block))))
        charge_projected_basis.append(block.reshape(-1))
    operator_rank = int(np.linalg.matrix_rank(np.stack(charge_projected_basis)))
    direction_rank = int(np.linalg.matrix_rank(np.array([dirs[f] for f, _, _ in active]))) if active else 0
    print(f"  active charge-preserving directions: {active}")
    print(f"  independent charge-preserving B_f operators: {operator_rank}")
    print(f"  span rank of active body-diagonal directions: {direction_rank}")
    check(operator_rank == 1, "charge-preserving projection collapses to one independent hop operator")
    check(direction_rank == 1, "charge-preserving projection is one-dimensional, not a 3D fermion kernel")

    h_q = h_pp(sample_k, flips, dirs, physical, index, keep_charge=True)
    comm_q = float(np.linalg.norm(h_q @ q_op - q_op @ h_q))
    check(comm_q < 1e-12, "charge-projected Hamiltonian commutes with Q, but only after losing 3D support")

    print(
        r"""
================================================================================
VERDICT (item 115 / Q3 matter vacuum-polarization precondition)

The existing 256-register Final-Boss Bloch core cannot be used as the charged
matter Green function G_F for item 115's vacuum-polarization loop:

  * It is not closed on the 48 physical matter codewords; 256 of 384 one-hop
    physical attempts leave the physical space.
  * Its physical P-block does not commute with the item-116 electric charge; 96
    physical P->P hops change Q, so a bubble built from it would violate the
    U(1) Ward/gauge precondition before any Pi_mn(q) integral is attempted.
  * Hard-projecting to charge-preserving hops leaves only f=4, a single
    body-diagonal direction and one independent operator. That is a legal
    charge-conserving toy, but not a 3D Q3 matter-fermion kernel and not a
    vacuum-polarization mechanism capable of smoothing the T_1u/E_g anisotropy.

So item 115 is still half-recovered only: K7's tree-level sqrt(2/3) anisotropy
is real, but the Q3 matter-fermion loop is NOT pinned by the current 3D Bloch
core. The next theorem target is now exact: construct H_matter(k) on the
physical 48 (or 3 x 16) matter states such that P H P = H, [H,Q]=0, the Peierls
current satisfies the lattice Ward identity, and the charge-preserving spatial
support has rank 3. Only then is Pi_mn(q) a meaningful item-115 computation.
================================================================================
"""
    )
    print("exit 0 -- pre-loop audit passed: current Q3 core fails the charged-matter loop gates.")


if __name__ == "__main__":
    main()
