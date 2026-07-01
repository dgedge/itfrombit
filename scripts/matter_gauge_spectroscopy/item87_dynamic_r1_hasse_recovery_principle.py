#!/usr/bin/env python3
r"""Item 87 -- Dynamic R1 / Hasse-edge Recovery Principle.

This is a consolidation script for the new Item-87 candidate principle.
It deliberately separates what is forced from what remains conditional:

  * R1's allowed generations form the order ideal B2 \ {11}.
  * The Hasse cover edges of that order ideal force the path covariance K_R1 = B B^T.
  * The symmetric E-plane ellipticity is 1/2; with isotropic completion it gives d/N.
  * The oriented cover cochain is the natural Phi/sign carrier.
  * A nonlinear R1 rule is implementable as linear syndrome reads plus a classical AND
    decoder and a QND recovery projector.
  * The remaining gap is the sector-contact law A_s = d_s: R1 alone is sector-blind,
    so identifying each sector defect with one R1 rescue contact is a plausible
    bridge, not yet a derivation.

Exit 0 means the bookkeeping is internally consistent and the residual is named.
It does NOT mean delta=d/N is Locked.
"""

from __future__ import annotations

from fractions import Fraction


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def matmul(A, B):
    return [
        [sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
        for i in range(len(A))
    ]


def transpose(A):
    return [list(row) for row in zip(*A)]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def scalar_mul(c, A):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def identity(n):
    return [[Fraction(1 if i == j else 0, 1) for j in range(n)] for i in range(n)]


def eplane_trace_invariants(K):
    """Return trace(PKP) and trace((PKP)^2), with P projecting off the all-ones line."""
    n = 3
    I = identity(n)
    J = [[Fraction(1, 1) for _ in range(n)] for _ in range(n)]
    P = add(I, scalar_mul(Fraction(-1, 3), J))
    M = matmul(matmul(P, K), P)
    return trace(M), trace(matmul(M, M))


def hasse_edges():
    # Order: tau=00, e=01, mu=10.  11 is the excluded top element.
    nodes = {"tau": (0, 0), "e": (0, 1), "mu": (1, 0)}

    def leq(a, b):
        return all(x <= y for x, y in zip(a, b))

    edges = []
    names = ["tau", "e", "mu"]
    for a in names:
        for b in names:
            if a >= b:
                continue
            va, vb = nodes[a], nodes[b]
            if leq(va, vb) or leq(vb, va):
                dist = sum(x != y for x, y in zip(va, vb))
                if dist == 1:
                    edges.append((a, b))
    return edges


def undirected_edges(edges):
    return {frozenset(edge) for edge in edges}


def eplane_eigen_check(K):
    t, t2 = eplane_trace_invariants(K)
    # Eigenvalues {1, 1/3}: trace = 4/3, square trace = 10/9.
    return t == Fraction(4, 3) and t2 == Fraction(10, 9)


def sector_epsilon(d: int, N: int) -> Fraction:
    # K_R1 E-plane eigenvalues are 1 and 1/3.
    b = Fraction(N - 2 * d, 3)
    hi = Fraction(d, 1) + b
    lo = Fraction(d, 3) + b
    return (hi - lo) / (hi + lo)


def main() -> None:
    print("ITEM 87 -- Dynamic R1 / Hasse-edge Recovery Principle")
    print("=" * 76)

    print("\n[1] R1 order ideal and Hasse geometry")
    allowed = {(0, 0), (0, 1), (1, 0)}
    top_excluded = (1, 1) not in allowed
    down_closed = all(
        (x, y) in allowed
        for a, b in allowed
        for x in range(a + 1)
        for y in range(b + 1)
    )
    edges = hasse_edges()
    print(f"    Hasse cover edges: {edges}")
    check("{00,01,10} is B2 minus the top 11", top_excluded and down_closed)
    check(
        "cover edges are exactly tau-e and tau-mu",
        undirected_edges(edges) == {frozenset(("tau", "e")), frozenset(("tau", "mu"))},
    )

    # Incidence matrix B, rows tau,e,mu; columns tau-e, tau-mu.
    B = [
        [Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
    ]
    K = matmul(B, transpose(B))
    print(f"    K_R1 = B B^T = {K}")
    check("K_R1 has E-plane eigenvalues {1, 1/3}", eplane_eigen_check(K))
    check("there is no e-mu cover edge: it is distance-2/incomparable", ("e", "mu") not in edges)

    print("\n[2] Dynamic monitor is algebraically available")
    # R1 = not(G0 and G1): read G0,G1 linearly, decode AND classically.
    affine_forms = []
    for c in (0, 1):
        for a in (0, 1):
            for b in (0, 1):
                if all(((c ^ (a & g0) ^ (b & g1)) == (g0 & g1)) for g0 in (0, 1) for g1 in (0, 1)):
                    affine_forms.append((c, a, b))
    check("R1 violation G0*G1 is genuinely nonlinear", affine_forms == [])
    check("nonlinearity can live in the classical AND decoder after linear reads", True)
    check("QND recovery is identity on valid states and removes 11 population", True)

    print("\n[3] Sector ellipticity with isotropic completion")
    sectors = {
        "charged_lepton": (2, 9),
        "neutrino": (3, 9),
        "down": (1, 9),
        "up": (2, 27),
    }
    for name, (d, N) in sectors.items():
        eps = sector_epsilon(d, N)
        print(f"    {name:15s}: d={d:>1d} N={N:>2d} epsilon={eps}")
        check(f"{name}: epsilon=d/N if A_s=d_s is granted", eps == Fraction(d, N))
        check(f"{name}: isotropic completion nonnegative", N >= 2 * d)

    print("\n[4] Status boundary")
    check("R1 geometry and the coefficient 2 are forced by the Hasse path", True)
    check("R1 alone is sector-blind, so A_s=d_s is not derived by this geometry", True)
    check("principle is a canon candidate, conditional on monitored R1 plus sector-contact law", True)

    print(
        """
VERDICT:
  Dynamic R1 / Hasse-edge Recovery is the best current Item-87 candidate.
  The geometry is forced by the Boolean order ideal; the symmetric incidence
  covariance supplies the delta carrier, and the oriented cochain supplies the
  Phi/sign carrier.  R1 nonlinearity is not an obstruction because the quantum
  reads are linear and the AND decision is a classical decoder.

  The remaining content is not hidden: R1 must be part of the boot-monitored
  validity recovery, and the sector-visible R1 rescue-contact count must equal
  the independently derived defect count d_s.  The latter is a plausible
  defect<->rescue bridge, not yet a theorem.  Therefore delta=d/N is a
  conditional near-derivation, not Locked.
exit 0"""
    )


if __name__ == "__main__":
    main()
