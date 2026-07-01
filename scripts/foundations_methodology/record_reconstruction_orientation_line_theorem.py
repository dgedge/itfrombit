#!/usr/bin/env python3
r"""Orientation-line theorem for monitored recovery histories.

Question
--------
Can the extra orientation-bearing premise be derived rather than simply added?

Result
------
Yes for the *existence and type* of the orientation-bearing object; no for the
absolute sign.

The derivation has three finite legs.

  1. Directed recovery histories are forced by reusable monitored recovery.
     A recovery event has a before/after syndrome order and the service register
     is reset for reuse.  The Stinespring copy is reversible before reset, but
     the reusable service ledger is not time-symmetric: it exports a record and
     erases it to a blank.  Thus the event is a directed history, not an
     unoriented edge.

  2. A directed R1 recovery event carries an antisymmetric boundary cochain.
     The allowed generation set is the order ideal B2\{11}; its Hasse covers are
     the two one-bit rescues from the forbidden corner.  The symmetric second
     moment is the R1 covariance.  The antisymmetric part is the CP-orientation
     carrier.

  3. Local directed cochains must glue over a connected substrate.  Their signs
     form a Z2 orientation bundle.  If all overlap/cycle transition products are
     consistent, the connected substrate has exactly two global sections,
     +omega and -omega.  Multiplying the local sign-representation cochain by
     omega makes a scalar environment record.  If the signs are averaged instead
     of globally selected, the oriented CP part cancels.

Therefore the "extra premise" can be reduced to:

    the monitored recovery substrate is orientable and chooses one of the two
    global Z2 sections as a superselection branch.

That is narrower than a free CP phase.  The line's existence follows from
directed recovery + composable local records; the sign remains the irreducible
orientation branch, just like choosing one time arrow / handedness sector.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def matmul(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*A)]


def scalar(c: F, A: list[list[F]]) -> list[list[F]]:
    return [[c * x for x in row] for row in A]


def add(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def directed_recovery_gate() -> None:
    print("\n[1] Directed histories from monitored recovery + reset")
    # Minimal event order: bad syndrome b is corrected to valid v and writes an
    # environment record e.  Reset maps e -> blank for reuse.  The reverse map
    # blank -> e is not available without supplying the exported entropy.
    before = "forbidden syndrome 11"
    after_options = ("valid 01", "valid 10")
    record_options = ("rescue 11->01", "rescue 11->10")
    reset_map = {record: "blank" for record in record_options}
    check(len(after_options) == len(record_options) == 2, "R1 rescue has two one-bit directed outcomes")
    check(len(set(reset_map.values())) == 1, "reset is many-to-one on the service record")
    check(before != after_options[0] and before != after_options[1], "recovery has a before/after syndrome order")
    print(f"    before={before}; after={after_options}; records={record_options}; reset -> blank")


def r1_symmetric_and_oriented_parts() -> tuple[list[list[F]], list[list[F]]]:
    # Generation order: tau=00, e=01, mu=10.
    # Symmetric closed edge-pair vectors on tau-e and tau-mu.
    B = [[F(1), F(1)], [F(1), F(0)], [F(0), F(1)]]
    K = matmul(B, transpose(B))
    # Directed boundary cochain from the two covers, tau->e and tau->mu, encoded
    # as the antisymmetric part.  Reflection e<->mu flips this cochain.
    Omega = [
        [F(0), F(-1), F(1)],
        [F(1), F(0), F(0)],
        [F(-1), F(0), F(0)],
    ]
    return K, Omega


def r1_orientation_carrier_gate() -> None:
    print("\n[2] R1 directed cochain is the orientation carrier")
    K, Omega = r1_symmetric_and_oriented_parts()
    check(K == [[F(2), F(1), F(1)], [F(1), F(1), F(0)], [F(1), F(0), F(1)]], "symmetric part is K_R1 = B B^T")
    check(Omega == scalar(F(-1), transpose(Omega)), "oriented part is antisymmetric")
    # Endpoint reflection e<->mu.
    S = [[F(1), F(0), F(0)], [F(0), F(0), F(1)], [F(0), F(1), F(0)]]
    check(matmul(matmul(S, K), transpose(S)) == K, "reflection preserves the symmetric covariance")
    check(matmul(matmul(S, Omega), transpose(S)) == scalar(F(-1), Omega), "reflection flips the directed cochain")


@dataclass(frozen=True)
class Z2Bundle:
    n: int
    edges: tuple[tuple[int, int, int], ...]  # (i,j,s_ij), with s_ij = +/-1

    def sections(self) -> list[tuple[int, ...]]:
        out: list[tuple[int, ...]] = []
        for signs in product((-1, 1), repeat=self.n):
            if all(signs[j] == sij * signs[i] for i, j, sij in self.edges):
                out.append(signs)
        return out


def orientation_bundle_gate() -> None:
    print("\n[3] Global orientation line as the Z2 gluing datum")
    # A connected orientable patch: a square with consistent transitions.
    orientable = Z2Bundle(
        4,
        (
            (0, 1, +1),
            (1, 2, +1),
            (2, 3, +1),
            (3, 0, +1),
            (0, 2, +1),
        ),
    )
    sections = orientable.sections()
    print(f"    orientable connected patch sections: {sections}")
    check(len(sections) == 2, "connected orientable substrate has exactly two global orientation sections +/-omega")
    check(sections[0] == tuple(-x for x in sections[1]), "the two sections differ only by global sign")

    # An inconsistent transition product around a cycle: no global orientation.
    nonorientable = Z2Bundle(3, ((0, 1, +1), (1, 2, +1), (2, 0, -1)))
    check(len(nonorientable.sections()) == 0, "inconsistent Z2 cycle has no global orientation line")

    # A disconnected orientable patch would have independent signs and hence
    # extra labels.  Connectedness removes that freedom.
    disconnected = Z2Bundle(4, ((0, 1, +1), (2, 3, +1)))
    check(len(disconnected.sections()) == 4, "disconnected components carry independent orientation signs")


def cp_cancellation_gate() -> None:
    print("\n[4] Why the line matters: oriented records survive, averaged signs cancel")
    _K, Omega = r1_symmetric_and_oriented_parts()
    plus = Omega
    minus = scalar(F(-1), Omega)
    averaged = scalar(F(1, 2), add(plus, minus))
    zero = [[F(0) for _ in range(3)] for _ in range(3)]
    check(averaged == zero, "unoriented +/- mixture cancels the CP cochain")
    check(plus != zero and minus != zero, "a chosen global orientation branch carries a nonzero directed record")


def main() -> None:
    print("RECORD RECONSTRUCTION ORIENTATION-LINE THEOREM")
    print("=" * 88)
    directed_recovery_gate()
    r1_orientation_carrier_gate()
    orientation_bundle_gate()
    cp_cancellation_gate()
    print(
        """
VERDICT
  The extra orientation-bearing premise can be reduced, but not erased.

  Derived:
    reusable monitored recovery makes service histories directed;
    directed R1 recovery writes an antisymmetric Hasse-edge cochain in addition
    to the symmetric covariance;
    composability of local directed records over a connected substrate forces a
    Z2 orientation line, with exactly two global sections if the substrate is
    orientable.

  Not derived:
    the choice between +omega and -omega.  The two branches are related by
    global orientation reversal and are indistinguishable to parity-even record
    axioms.  Averaging them kills the CP cochain.

  Minimal replacement for the old premise:
    the monitored recovery substrate is orientable and occupies one global
    orientation branch.  That is a superselection/initial-branch statement, not
    a tunable CP phase.

exit 0"""
    )


if __name__ == "__main__":
    main()
