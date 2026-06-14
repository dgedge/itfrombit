#!/usr/bin/env python3
"""
Can current TCH canon prove the neutrino phase-lift/principal-branch rule?

Target lemma:
    d_eff/N = 3/9  ==>  arg(B_nu) lies in the principal Voronoi cell
                         around 1/3, excluding remote IH branches.

Result:
    Not from the current documented operator.

Reason:
    The actual Part 18 Feshbach operator has H_QP H_PQ = 2 I_3 in generation
    space. Its generation-circulant off-diagonal coefficient B_nu is exactly
    zero, so arg(B_nu) is undefined. It supplies degeneracy/rank structure,
    not a phase-lift.

    The Part 05 Koide mass formula does contain a nonzero B_nu, but only
    after delta_nu has already been inserted into the scalar ansatz. Extracting
    arg(B_nu) from that matrix simply returns the input phase.

Run:
    python3 verify_neutrino_phase_lift_obstruction.py
"""

from __future__ import annotations

import cmath
import math


OMEGA = complex(math.cos(2.0 * math.pi / 3.0), math.sin(2.0 * math.pi / 3.0))


def matmul(a: list[list[complex]], b: list[list[complex]]) -> list[list[complex]]:
    rows, inner, cols = len(a), len(b), len(b[0])
    return [
        [sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)]
        for i in range(rows)
    ]


def transpose(a: list[list[complex]]) -> list[list[complex]]:
    return [list(row) for row in zip(*a)]


def make_part18_h_pq(g: float = 1.0) -> list[list[complex]]:
    """Six active rows, three nu_R columns; two parallel entries per generation."""

    h = [[0j for _ in range(3)] for _ in range(6)]
    for gen in range(3):
        h[2 * gen][gen] = complex(g)
        h[2 * gen + 1][gen] = complex(g)
    return h


def circulant_coefficients_from_first_row(row: list[complex]) -> tuple[complex, complex, complex]:
    """
    For a 3x3 circulant with first row [c0,c1,c2], eigenvalues are
      lambda_n = c0 + c1 omega^n + c2 omega^(2n).
    For a real/Hermitian nearest-neighbour circulant, c2 = conj(c1), and the
    Koide-like off-diagonal coefficient is B = c1.
    """

    return row[0], row[1], row[2]


def inverse_dft_to_first_row(eigs: list[complex]) -> list[complex]:
    """First row of the circulant whose eigenvalues are eigs[n]."""

    row = []
    for k in range(3):
        row.append(sum(eigs[n] * OMEGA ** (-n * k) for n in range(3)) / 3.0)
    return row


def koide_amplitudes(R: float, delta: float) -> list[float]:
    return [1.0 + R * math.cos(delta + 2.0 * math.pi * n / 3.0) for n in range(3)]


def angle(z: complex) -> float | None:
    if abs(z) < 1e-14:
        return None
    a = math.atan2(z.imag, z.real)
    return a % (2.0 * math.pi)


def format_complex(z: complex) -> str:
    if abs(z.imag) < 1e-12:
        return "%.12g" % z.real
    return "%.12g%+.12gj" % (z.real, z.imag)


def print_matrix(label: str, matrix: list[list[complex]]) -> None:
    print(label)
    for row in matrix:
        print("   " + " ".join("%22s" % format_complex(x) for x in row))


def main() -> None:
    print("=" * 78)
    print("Neutrino phase-lift obstruction verifier")
    print("=" * 78)

    h_pq = make_part18_h_pq()
    h_qp_h_pq = matmul(transpose(h_pq), h_pq)
    print_matrix("A. Actual Part 18 generation-space H_QP H_PQ:", h_qp_h_pq)
    first_row = h_qp_h_pq[0]
    c0, c1, c2 = circulant_coefficients_from_first_row(first_row)
    print("   circulant coefficients from first row:")
    print("     A  = %s" % format_complex(c0))
    print("     B  = %s" % format_complex(c1))
    print("     B* = %s" % format_complex(c2))
    print("     arg(B) = %s" % ("undefined" if angle(c1) is None else "%.12f" % angle(c1)))
    assert abs(c1) < 1e-14 and abs(c2) < 1e-14
    print("   conclusion: actual Part 18 Feshbach block has no nonzero B_nu.")
    print()

    delta = 1.0 / 3.0
    eigs = [complex(x) for x in koide_amplitudes(1.0, delta)]
    row = inverse_dft_to_first_row(eigs)
    a, b, bstar = circulant_coefficients_from_first_row(row)
    b_arg = angle(b)
    # Depending on DFT convention this is delta or -delta modulo 2pi. Fold to
    # the smaller absolute representative to compare with the input.
    folded = min(b_arg or 0.0, 2.0 * math.pi - (b_arg or 0.0))
    print_matrix("B. Part 05 Koide ansatz circulant reconstructed from inserted delta=1/3:", [
        row,
        [row[2], row[0], row[1]],
        [row[1], row[2], row[0]],
    ])
    print("   extracted |B|  = %.12f" % abs(b))
    print("   extracted arg = %.12f folded to %.12f" % (b_arg or 0.0, folded))
    print("   input delta   = %.12f" % delta)
    print("   conclusion: nonzero B_nu appears only after the Koide phase is inserted.")
    print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  Current documented operators do not prove the phase-lift rule.")
    print("  Part 18 gives B_nu = 0 in generation space, so there is no arg(B_nu)")
    print("  to bound near 3/9. Part 05 gives a nonzero B_nu only by assuming")
    print("  delta_nu in the Koide ansatz, so extracting arg(B_nu) is circular.")
    print()
    print("  To prove the NH theorem, the framework needs a new explicit")
    print("  generation-circulant neutrino operator with nonzero B_nu derived from")
    print("  the substrate, then a bound |arg(B_nu)-1/3| < 1/18.")


if __name__ == "__main__":
    main()
