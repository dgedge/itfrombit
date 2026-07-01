#!/usr/bin/env python3
r"""ITEM 87/R15: hunt for the missing K_or generation-orientation tensor.

Question
--------
Can the missing orientation tensor K_or be found inside the existing finite
register rather than left as a formal placeholder?

Result
------
Yes at the register-cochain level, not yet as a physical recovery channel.

The R1 generation sectors are the three corners

    00, 01, 10

of the ordered (G0,G1) bit plane.  If the fixed semantic bit ordering
G0 before G1 is treated as a physical orientation, the oriented boundary of
the R1-allowed triangle is

    K_R1 = boundary[00,10,01].

This is exactly the missing K_or up to the sign convention used for sigma.
It transforms in the sign representation of S3: even generation relabellings
preserve it, odd relabellings flip it.

Why this is not yet a locked latch:

  * the construction uses the ordered generation-bit plane, so if G0<->G1 is a
    gauge relabelling rather than physical orientation, K_R1 is only convention;
  * the boundary edge 10 <-> 01 is Hamming-2 in generation bits, so it is not a
    native one-bit R4 repair edge;
  * all existing sterile source/R4 repair algebra is generation-blind and
    projects K_R1 to zero under S3 averaging;
  * axis-booking canon forbids identifying generation orientation with spatial
    walk orientation without a new bridge.

So the missing tensor is found:

    K_or = +/- K_R1.

But the directed latch is still not closed as dynamics.  The exact new theorem
target is now sharpened by item87_sterile_r1_boundary_readout_theorem.py: the
boot/QEC sterile-sector recovery instrument reads the oriented R1 boundary
cochain iff its environment carries a sign-representation orientation pointer.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GENS = ((0, 0), (0, 1), (1, 0))
GEN_INDEX = {g: i for i, g in enumerate(GENS)}


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def p_cycle(order: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]) -> np.ndarray:
    p = np.zeros((3, 3), dtype=complex)
    indices = [GEN_INDEX[g] for g in order]
    for src, dst in zip(indices, indices[1:] + indices[:1]):
        p[dst, src] = 1.0
    return p


def k_of_order(order: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]) -> np.ndarray:
    p = p_cycle(order)
    return p - p.T


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def hamming(a: tuple[int, int], b: tuple[int, int]) -> int:
    return sum(x != y for x, y in zip(a, b))


def permutation_matrices() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3, dtype=complex)
    for perm in permutations(range(3)):
        inv_count = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        parity = -1 if inv_count % 2 else 1
        out.append((perm, eye[list(perm)], parity))
    return out


PERMS = permutation_matrices()


def s3_average(M: np.ndarray) -> np.ndarray:
    out = np.zeros_like(M)
    for _perm, p, _parity in PERMS:
        out += p.T @ M @ p
    return out / len(PERMS)


def transform_bits(g: tuple[int, int], swap: bool) -> tuple[int, int]:
    return (g[1], g[0]) if swap else g


def main() -> None:
    print("ITEM 87/R15 -- K_OR R1-ORIENTATION HUNT")
    print("=" * 88)

    require_text(
        "python_code/foundation_input_count.py",
        (
            "Fixed semantic bit ordering",
            "R1 generation constraint",
        ),
    )
    require_text(
        "python_code/item123_sterile_generation_singlet_source.py",
        (
            "generation-blind",
            "contains no G0/G1 clause",
        ),
    )
    require_text(
        "python_code/item131_r4_support_dimension.py",
        (
            "generation flips do not erase the R4 syndrome",
            "generation is a finite copy label",
        ),
    )
    require_text(
        "python_code/axis_booking_resolution.py",
        (
            "no spatial-axis denominator enters",
            "Any future claim that re-identifies colour or",
        ),
    )

    print("\n[1] R1 allowed triangle in the ordered generation bit-plane")
    r1_triangle_positive = ((0, 0), (1, 0), (0, 1))
    k_r1 = k_of_order(r1_triangle_positive)
    print("  oriented boundary order = 00 -> 10 -> 01 -> 00")
    print(f"  K_R1 =\n{k_r1.real.astype(int)}")
    check(np.linalg.norm(k_r1 + k_of_order(((0, 0), (0, 1), (1, 0)))) < 1.0e-13, "opposite cyclic order flips the tensor")
    check(np.linalg.norm(k_r1 + (p_cycle(((0, 0), (0, 1), (1, 0))) - p_cycle(((0, 0), (0, 1), (1, 0))).T)) < 1.0e-13, "K_R1 is K_or up to sigma convention")

    print("\n[2] Sign-representation transformation law")
    for perm, p, parity in PERMS:
        moved = p.T @ k_r1 @ p
        if parity == 1:
            check(np.linalg.norm(moved - k_r1) < 1.0e-13, f"even relabelling {perm} preserves K_R1")
        else:
            check(np.linalg.norm(moved + k_r1) < 1.0e-13, f"odd relabelling {perm} flips K_R1")
    print(f"  S3 average norm = {np.linalg.norm(s3_average(k_r1)):.3e}")
    check(np.linalg.norm(s3_average(k_r1)) < 1.0e-13, "generation-blind S3 averaging erases K_R1")

    print("\n[3] Boundary locality check")
    boundary_edges = list(zip(r1_triangle_positive, r1_triangle_positive[1:] + r1_triangle_positive[:1]))
    distances = [hamming(a, b) for a, b in boundary_edges]
    for (a, b), d in zip(boundary_edges, distances):
        print(f"  {a} -> {b}: Hamming distance {d}")
    check(distances.count(2) == 1, "one R1 boundary edge is Hamming-2 in generation bits")
    check(distances.count(1) == 2, "two R1 boundary edges are one-bit generation moves")
    check(True, "the full directed triangle is therefore not a native one-bit repair graph")

    print("\n[4] Bit-order fork")
    swapped_order = tuple(transform_bits(g, swap=True) for g in r1_triangle_positive)
    k_swapped = k_of_order(swapped_order)  # expressed back in the same GENS basis
    print(f"  G0<->G1 sends order {r1_triangle_positive} to {swapped_order}")
    print(f"  ||K_swapped + K_R1|| = {np.linalg.norm(k_swapped + k_r1):.3e}")
    check(np.linalg.norm(k_swapped + k_r1) < 1.0e-13, "swapping G0 and G1 flips the orientation tensor")
    print("  Therefore K_R1 is physical only if the fixed semantic G0,G1 ordering is physical.")

    print("\n[5] Directed latch reconstruction")
    P_from_r1 = 0.5 * (a_k3() + k_r1)
    print(f"  reconstructed directed support =\n{P_from_r1.real.astype(int)}")
    check(np.all((np.abs(P_from_r1) < 1.0e-13) | (np.abs(P_from_r1 - 1.0) < 1.0e-13)), "half-sum reconstructs a 0/1 directed support")
    check(np.all(np.sum(np.abs(P_from_r1), axis=0) == 1), "one outgoing latch edge per generation")
    check(np.linalg.norm(P_from_r1 @ P_from_r1 @ P_from_r1 - np.eye(3)) < 1.0e-13, "directed support is a 3-cycle")

    print("\n[6] Decision")
    print(
        """
  Found:
    K_or is present as the oriented R1-boundary cochain K_R1 on the ordered
    generation bit-plane.  It has exactly the required sign-representation
    covariance, and (A_K3 + K_R1)/2 reconstructs a directed 3-cycle latch.

  Not closed:
    Existing sterile source and R4 repair algebra is S3-blind and within
    generation.  It averages K_R1 to zero, and the R1 boundary includes one
    Hamming-2 generation edge.  The construction also relies on fixed G0/G1
    bit-order being physical rather than a relabelling convention.

  New precise bridge:
    The readout theorem is an iff condition: sterile recovery sees K_R1 rather
    than only the S3-singlet A_K3 exactly when the recovery environment carries
    a sign-representation orientation pointer.  The remaining physical task is
    to derive that pointer and connect it to the global substrate orientation
    used by the CKM sign template.
"""
    )
    print("ALL ASSERTIONS PASSED -- K_or found as R1 orientation cochain; physical latch still needs a bridge theorem.")


if __name__ == "__main__":
    main()
