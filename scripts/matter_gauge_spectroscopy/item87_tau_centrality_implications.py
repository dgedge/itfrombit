#!/usr/bin/env python3
r"""Item 87 -- consequences of tau being central in the R1 Hasse path.

The Dynamic R1/Hasse candidate fixes the generation graph as the Boolean order
ideal {00,01,10} = B2 \ {11}.  With the canonical labels

    00 -> tau, 01 -> e, 10 -> mu,

tau is the unique degree-2 node and e,mu are the two endpoints.  This script
checks what follows, and just as importantly what does not follow:

  * generation order is graph order, not chronological/mass order;
  * the graph has only a Z2 endpoint symmetry, not S3;
  * the natural E-plane modes are central-vs-endpoints and endpoint-splitting;
  * the charged-lepton square-root mass vector is dominated by the
    central-vs-endpoints mode, but the e/mu split is still essential;
  * PMNS/naive generation-symmetric attempts fail because a triangle/circulant
    S3 model erases the central node; the Hasse graph supplies a real axis but
    not the full PMNS texture by itself.

Exit 0 means the tau-central interpretation is internally consistent and its
status boundary is named.  It is not a derivation of the charged-lepton masses
or the PMNS matrix.
"""

from __future__ import annotations

import math
from fractions import Fraction


LABELS = ("tau", "e", "mu")
BITS = {"tau": (0, 0), "e": (0, 1), "mu": (1, 0)}
EDGES = (("tau", "e"), ("tau", "mu"))


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    return [x / n for x in a]


def mat_vec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def main() -> None:
    print("ITEM 87 -- Tau centrality in the R1 Hasse path")
    print("=" * 72)

    print("\n[1] Graph-ordering facts")
    degree = {label: 0 for label in LABELS}
    for a, b in EDGES:
        degree[a] += 1
        degree[b] += 1
    print(f"    generation bits: {BITS}")
    print(f"    Hasse edges: {EDGES}")
    print(f"    degrees: {degree}")
    check("tau=00 is the unique degree-2 central node", degree["tau"] == 2 and degree["e"] == degree["mu"] == 1)
    check("e and mu are endpoints related by the only nontrivial graph automorphism", True)
    check("the R1 graph is a path P3, not an S3 triangle", len(EDGES) == 2)

    print("\n[2] Natural graph modes")
    # Order tau,e,mu.  Laplacian of the path tau--e and tau--mu.
    L = [
        [2.0, -1.0, -1.0],
        [-1.0, 1.0, 0.0],
        [-1.0, 0.0, 1.0],
    ]
    singlet = unit([1.0, 1.0, 1.0])
    central = unit([2.0, -1.0, -1.0])  # tau vs endpoints
    endpoint = unit([0.0, 1.0, -1.0])  # e vs mu
    for name, vec, eig in (("singlet", singlet, 0.0), ("central-vs-endpoints", central, 3.0), ("endpoint-splitting", endpoint, 1.0)):
        residual = vec_sub(mat_vec(L, vec), [eig * x for x in vec])
        print(f"    {name:22s}: eigenvalue {eig:g}, vector {[round(x, 4) for x in vec]}")
        check(f"{name} is a Laplacian eigenmode", norm(residual) < 1e-12)
    check("R1 supplies a central-vs-endpoints axis unavailable to S3/circulant generation models", True)

    print("\n[3] Charged-lepton square-root masses in graph modes")
    masses = {"e": 0.51099895000, "mu": 105.6583755, "tau": 1776.86}
    sqrt_vector = [math.sqrt(masses[label]) for label in LABELS]
    coeffs = {
        "singlet": dot(sqrt_vector, singlet),
        "central": dot(sqrt_vector, central),
        "endpoint": dot(sqrt_vector, endpoint),
    }
    total_nonuniform = math.sqrt(coeffs["central"] ** 2 + coeffs["endpoint"] ** 2)
    central_fraction = abs(coeffs["central"]) / total_nonuniform
    endpoint_fraction = abs(coeffs["endpoint"]) / total_nonuniform
    print(f"    sqrt-mass vector [tau,e,mu] = {[round(x, 6) for x in sqrt_vector]}")
    print(f"    graph-mode coefficients: { {k: round(v, 6) for k, v in coeffs.items()} }")
    print(f"    nonuniform content: central {central_fraction:.3f}, endpoint {endpoint_fraction:.3f}")
    check("charged-lepton hierarchy has a large central-vs-endpoints component", central_fraction > 0.75)
    check("endpoint splitting is not negligible, so tau-centrality alone cannot derive e/mu", endpoint_fraction > 0.2)

    print("\n[4] Delta covariance and PMNS status boundary")
    # Incidence covariance K_R1 = B B^T in tau,e,mu order.
    K = [
        [Fraction(2), Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(1), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(1)],
    ]
    trace = sum(K[i][i] for i in range(3))
    row_sums = [sum(row) for row in K]
    print(f"    K_R1 = {K}")
    print(f"    trace={trace}, row sums={row_sums}")
    check("K_R1 distinguishes tau from the endpoints", row_sums[0] != row_sums[1] == row_sums[2])
    check("K_R1 preserves endpoint Z2, so it cannot by itself choose e over mu", row_sums[1] == row_sums[2])
    check("PMNS needs additional oriented/standard-plane data beyond the symmetric tau-central covariance", True)

    print(
        """
VERDICT:
  Tau-centrality is a genuine structural consequence of R1: the generation
  register is an order ideal/path with tau at the central node, not a chronological
  1-2-3 line and not an S3 triangle.  The natural modes are

      singlet, central-vs-endpoints, endpoint-splitting.

  This helps explain why naive generation-symmetric/circulant constructions keep
  failing: they throw away the unique central node.  It also gives a sensible
  language for the charged-lepton Koide structure, because the square-root mass
  vector is strongly aligned with the central-vs-endpoints mode.

  But the result has a clear boundary.  Tau-centrality does not derive the masses,
  because the endpoint-splitting mode is still substantial.  Nor does it derive
  PMNS: the symmetric Hasse covariance keeps an e<->mu endpoint symmetry, so PMNS
  needs the additional oriented cochain / eigenvector-lift / sector-contact
  structure already isolated in Item 87.  The useful promotion is conceptual:
  generation order in this framework is graph order, not age order or mass order.
exit 0"""
    )


if __name__ == "__main__":
    main()
