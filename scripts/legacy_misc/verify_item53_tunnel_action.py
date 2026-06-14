#!/usr/bin/env python3
"""
Verify/refute the proposed item-53 tunnelling-action lemma:

    S_tunnel = (# nonzero H_PQ channels) x d_code = 6 x 4 = 24.

Result:
    Refuted for the actual Part 18 Feshbach operator.

Reason:
    The six H_PQ entries are parallel couplings into three nu_R columns,
    not six serial barriers. Feshbach elimination uses

        Sigma(E) = H_PQ (E I - H_QQ)^-1 H_QP,

    so the relevant channel algebra is H_QP H_PQ = 2 I_3.

    If each H_PQ entry has attenuation g = exp(-d_code), then each
    nonzero self-energy eigenvalue scales as 2 g^2, with action

        S_eigen = -log(2 exp(-2 d_code)) = 2 d_code - log(2),

    not 6 d_code.

Run:
    python3 verify_item53_tunnel_action.py
"""

from __future__ import annotations

import math


D_CODE = 4
N_Q = 3
CHANNELS_PER_Q = 2
N_CHANNELS = N_Q * CHANNELS_PER_Q


def make_h_pq(g: float = 1.0) -> list[list[float]]:
    """6x3 Part-18 H_PQ block: two active lepton rows per nu_R column."""
    matrix = [[0.0 for _ in range(N_Q)] for _ in range(N_CHANNELS)]
    for q in range(N_Q):
        matrix[2 * q][q] = g
        matrix[2 * q + 1][q] = g
    return matrix


def transpose(matrix: list[list[float]]) -> list[list[float]]:
    return [list(row) for row in zip(*matrix)]


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    rows, inner, cols = len(a), len(b), len(b[0])
    return [
        [sum(a[i][k] * b[k][j] for k in range(inner)) for j in range(cols)]
        for i in range(rows)
    ]


def diagonal(matrix: list[list[float]]) -> list[float]:
    return [matrix[i][i] for i in range(min(len(matrix), len(matrix[0])))]


def nonzero_eigs_hqphpq(g: float) -> list[float]:
    # Because columns are mutually orthogonal and each has two entries g.
    return [2.0 * g * g] * N_Q


def pretty_matrix(matrix: list[list[float]]) -> str:
    return "\n".join("   " + " ".join("%10.6g" % x for x in row) for row in matrix)


def main() -> None:
    h_pq = make_h_pq(1.0)
    h_qp_h_pq = matmul(transpose(h_pq), h_pq)
    expected = [[2.0 if i == j else 0.0 for j in range(N_Q)] for i in range(N_Q)]
    assert h_qp_h_pq == expected

    g = math.exp(-D_CODE)
    eig = nonzero_eigs_hqphpq(g)[0]
    s_eigen = -math.log(eig)
    s_proposed = N_CHANNELS * D_CODE
    s_determinant = -math.log(math.prod(nonzero_eigs_hqphpq(g)))
    s_trace = -math.log(sum(nonzero_eigs_hqphpq(g)))

    assert abs(s_eigen - (2 * D_CODE - math.log(2))) < 1e-12
    assert abs(s_proposed - s_eigen) > 1.0

    print("=" * 78)
    print("Item-53 tunnel-action verifier")
    print("=" * 78)
    print("H_PQ has 6 nonzero entries, arranged as two parallel couplings per nu_R:")
    print(pretty_matrix(h_pq))
    print()
    print("H_QP H_PQ =")
    print(pretty_matrix(h_qp_h_pq))
    print()
    print("This equals 2 I_3, so the physical nonzero Feshbach eigenvalues are")
    print("parallel-channel sums, not products of six serial barriers.")
    print()
    print("With per-entry attenuation g = exp(-d_code), d_code = %d:" % D_CODE)
    print("  eigenvalue scale       = 2 exp(-2 d_code)")
    print("  eigenvalue action      = 2 d_code - log(2) = %.6f" % s_eigen)
    print("  trace action           = %.6f" % s_trace)
    print("  determinant action     = %.6f" % s_determinant)
    print("  proposed 6 d_code      = %.6f" % s_proposed)
    print()
    print("VERDICT: refuted as stated.")
    print("The factor 6 is a channel count/rank fact. It does not enter the")
    print("single-neutrino mass eigenvalue exponent as 6*d_code.")


if __name__ == "__main__":
    main()
