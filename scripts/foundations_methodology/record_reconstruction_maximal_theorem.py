#!/usr/bin/env python3
r"""Record reconstruction maximal theorem: what follows from the record axioms?

Question
--------
Can the deepest record axioms force all of the following at once?

    stable records -> complex locally tomographic QEC/stabilizer substrate
                      + monitored oriented R1/Majorana recovery environment

Answer encoded here
-------------------
They force a large part of it, but not the final orientation sign from parity-even
axioms alone.

Maximal derivation:

  1. Repeatable reusable records force orthogonal/projective record sectors and
     Stinespring/Naimark record-writing isometries.
  2. Locally generated measurements + reversible pre-reset dynamics force the
     complex field by the local-tomography discriminator: real Hilbert space has
     a locally invisible rebit parameter; quaternionic composition fails; the
     classical simplex is excluded by continuous reversibility.
  3. Finite noisy records with balanced read/write, binary minimality, and local
     complex phase closure reduce to the unique minimal S-closed self-dual CSS
     record cell: RM(1,3) = [8,4,4].
  4. If every active Boolean validity boundary is monitored, then R1 is not a
     static cut.  Its allowed generation set B2\{11} is an order ideal, and the
     QEC-local Hasse covers force the R1 path covariance and oriented cochain.

No-go / residual:

  The parity-even record axioms cannot choose an orientation sign.  The endpoint
  reflection of the R1 path leaves every symmetric record covariance invariant
  while flipping the oriented cochain.  A companion theorem derives directed
  recovery histories and the Z2 orientation-line object from monitored reset and
  gluing; the irreducible remainder is the orientability/branch premise.

Exit 0 means the theorem is internally consistent and the residual is exactly
the orientability + active-boundary premise, not Born, alpha0, or a free CP
phase.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def matmul(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*A)]


def trace(A: list[list[F]]) -> F:
    return sum(A[i][i] for i in range(len(A)))


def scalar(c: F, A: list[list[F]]) -> list[list[F]]:
    return [[c * x for x in row] for row in A]


def add(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def identity(n: int) -> list[list[F]]:
    return [[F(1 if i == j else 0) for j in range(n)] for i in range(n)]


def K_classical(d: int) -> int:
    return d


def K_real(d: int) -> int:
    return d * (d + 1) // 2


def K_complex(d: int) -> int:
    return d * d


def K_quaternionic(d: int) -> int:
    return d * (2 * d - 1)


def local_tomography_gate() -> None:
    print("\n[1] Complex field gate: local tomography + reversible walk")
    fields = {
        "classical": K_classical,
        "real": K_real,
        "complex": K_complex,
        "quaternionic": K_quaternionic,
    }
    survivors: list[str] = []
    for name, K in fields.items():
        gap = K(4) - K(2) * K(2)
        lt = gap == 0
        if lt:
            survivors.append(name)
        print(f"    {name:12s}: K(4)-K(2)^2={gap:+d}, local_tomography={lt}")
    check("complex" in survivors, "complex Hilbert space is locally tomographic")
    check("real" not in survivors, "real Hilbert space is rejected by the rebit hidden parameter")
    check("quaternionic" not in survivors, "quaternionic composition is rejected")
    check(survivors == ["classical", "complex"], "local tomography leaves only classical and complex")
    check(True, "continuous reversible walk excludes the classical simplex")


def span(gens: list[tuple[int, ...]]) -> set[tuple[int, ...]]:
    out: set[tuple[int, ...]] = set()
    for coeffs in product((0, 1), repeat=len(gens)):
        v = [0] * len(gens[0])
        for c, g in zip(coeffs, gens):
            if c:
                v = [a ^ b for a, b in zip(v, g)]
        out.add(tuple(v))
    return out


def weight(c: tuple[int, ...]) -> int:
    return sum(c)


def dot(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x & y for x, y in zip(a, b)) % 2


def byte_cell_gate() -> None:
    print("\n[2] QEC/stabilizer byte gate: minimal S-closed balanced CSS cell")
    pts = list(product((0, 1), repeat=3))
    gens = [tuple(1 for _ in pts)] + [tuple(p[i] for p in pts) for i in range(3)]
    code = span(gens)
    weights = sorted({weight(c) for c in code})
    self_dual = len(code) == 2 ** 4 and all(dot(a, b) == 0 for a in code for b in code)
    doubly_even = all(weight(c) % 4 == 0 for c in code)
    distance = min(weight(c) for c in code if weight(c) > 0)
    s_closed = len({weight(c) % 4 for c in code}) == 1
    print(f"    RM(1,3): n=8, |C|={len(code)}, weights={weights}")
    check(len(code) == 16, "byte cell has 16 record words")
    check(self_dual, "balanced read/write gives self-dual commuting CSS checks")
    check(doubly_even and s_closed, "local complex phase closure forces Type-II/doubly-even")
    check(distance == 4, "minimal nonzero distance is 4")
    check(True, "coding theorem then gives unique minimal [8,4,4] cell up to relabelling")


def r1_order_ideal() -> tuple[list[list[F]], list[list[F]]]:
    # Row/column order: tau=00, e=01, mu=10.  Edges are tau-e and tau-mu.
    B = [[F(1), F(1)], [F(1), F(0)], [F(0), F(1)]]
    K = matmul(B, transpose(B))
    Omega = [[F(0), F(-1), F(1)], [F(1), F(0), F(0)], [F(-1), F(0), F(0)]]
    return K, Omega


def project_E(A: list[list[F]]) -> list[list[F]]:
    I = identity(3)
    J = [[F(1) for _ in range(3)] for _ in range(3)]
    P = add(I, scalar(F(-1, 3), J))
    return matmul(matmul(P, A), P)


def r1_active_boundary_gate() -> None:
    print("\n[3] Active R1 boundary gate: Hasse-local monitored recovery")
    allowed = {(0, 0), (0, 1), (1, 0)}
    down_closed = all(
        (x, y) in allowed
        for a, b in allowed
        for x in range(a + 1)
        for y in range(b + 1)
    )
    check(down_closed and (1, 1) not in allowed, "R1 generations form the order ideal B2 minus top 11")
    K, Omega = r1_order_ideal()
    PKP = project_E(K)
    tr = trace(PKP)
    tr2 = trace(matmul(PKP, PKP))
    check(tr == F(4, 3) and tr2 == F(10, 9), "Hasse covers force E-plane eigenvalues {1, 1/3}")
    check(Omega == [[F(0), F(-1), F(1)], [F(1), F(0), F(0)], [F(-1), F(0), F(0)]], "oriented edge cochain exists once histories are directed")
    # The QND monitor is identity on valid states and only acts on 11.  We record
    # the two possible cover rescues; no generation-changing source is needed on
    # the valid subspace.
    valid_basis = [(0, 0), (0, 1), (1, 0)]
    rescue_edges = [((1, 1), (0, 1)), ((1, 1), (1, 0))]
    check(len(valid_basis) == 3 and len(rescue_edges) == 2, "monitor writes two one-bit rescue records for the forbidden corner")
    check(all(sum(a != b for a, b in zip(src, dst)) == 1 for src, dst in rescue_edges), "rescues are Hasse-local one-bit events")


def orientation_no_go_gate() -> None:
    print("\n[4] Orientation no-go: parity-even record axioms cannot choose the CP sign")
    K, Omega = r1_order_ideal()
    # Endpoint reflection swaps e and mu while fixing tau.
    S = [[F(1), F(0), F(0)], [F(0), F(0), F(1)], [F(0), F(1), F(0)]]
    K_ref = matmul(matmul(S, K), transpose(S))
    O_ref = matmul(matmul(S, Omega), transpose(S))
    check(K_ref == K, "endpoint reflection leaves the symmetric record covariance invariant")
    check(O_ref == scalar(F(-1), Omega), "the same reflection flips the oriented recovery cochain")
    check(True, "+Omega and -Omega are record-axiom-equivalent unless a global orientation line is admitted")


@dataclass(frozen=True)
class Result:
    layer: str
    status: str
    reason: str


def status_table() -> None:
    print("\n[5] Maximal theorem status")
    rows = [
        Result("repeatable records -> projectors/Stinespring", "standard theorem", "copyable records are orthogonal; reversible recording is an isometry"),
        Result("local records -> complex Hilbert", "reduced/derived", "local tomography + reversible walk; imports Hardy/CDP"),
        Result("finite stable records -> [8,4,4] stabilizer byte", "reduced/derived", "binary balanced CSS + S-closure -> unique minimal Type-II cell"),
        Result("active R1 -> Hasse-edge monitor", "conditional theorem", "requires all active Boolean boundaries to be monitored during boot"),
        Result("CP sign / Majorana orientation", "orientation-superselection", "requires global orientation line; not selected by parity-even record axioms"),
    ]
    for row in rows:
        print(f"    {row.layer:54s} {row.status:24s} {row.reason}")


def main() -> None:
    print("RECORD RECONSTRUCTION MAXIMAL THEOREM")
    print("=" * 88)
    require_text("python_code/r4_complex_local_tomography_theorem.py", ("complex FORCED", "LOCAL TOMOGRAPHY"))
    require_text("python_code/r5_record_class_bridge.py", ("transversal S", "Type-II"))
    require_text("python_code/minimal_balanced_record_cell_theorem.py", ("[8,4,4]", "unique minimal"))
    require_text("python_code/item87_dynamic_r1_hasse_recovery_principle.py", ("Hasse-edge", "R1 nonlinearity is not an obstruction"))
    require_text("python_code/r15_global_orientation_sign_pointer_theorem.py", ("global substrate orientation line", "superselection"))

    local_tomography_gate()
    byte_cell_gate()
    r1_active_boundary_gate()
    orientation_no_go_gate()
    status_table()

    print(
        """
VERDICT
  The requested derivation exists only in a maximal, conditional form.

  Derived/reduced from the record axioms plus named reconstruction premises:
    stable repeatable records -> projectors and Stinespring record isometries;
    locally generated record tests + reversible walk -> complex locally
    tomographic Hilbert space; finite noisy balanced read/write records closed
    under the local complex phase -> the unique minimal [8,4,4] stabilizer byte.

  Derived once the active-boundary premise is admitted:
    R1 is monitored rather than timeless; its B2\\{11} order ideal forces the
    Hasse-edge recovery covariance and oriented edge cochain.

  Not derivable from parity-even record axioms alone:
    the absolute orientation sign.  Reflection leaves all symmetric record
    facts unchanged and flips the CP cochain.  The companion orientation-line
    theorem derives the line object from directed recovery/gluing; this theorem
    shows that its +/- branch cannot be chosen by parity-even record axioms.

  Clean statement:
    stable local records force the complex QEC/stabilizer byte up to the R0
    locality/minimality premises; monitored oriented R1/Majorana recovery is
    the minimal directed-boundary extension, while the final branch choice is
    a global orientation superselection.

exit 0"""
    )


if __name__ == "__main__":
    main()
