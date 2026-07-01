#!/usr/bin/env python3
r"""ITEM 87/R15: when can sterile recovery read the oriented R1 boundary?

Question
--------
Can the framework prove that the sterile recovery instrument reads the
oriented R1-boundary cochain K_or, rather than only the S3-singlet Majorana
support A_K3?

Result
------
Yes, but only as a conditional representation theorem.

Let K_R1 be the oriented boundary cochain on the three R1 generation corners.
Under generation relabelling it transforms as the sign representation:

    pi . K_R1 = sgn(pi) K_R1.

Therefore a generation-blind Stinespring environment, whose pointer transforms
trivially under S3, cannot read K_R1.  The ordinary S3 average kills it and the
only surviving off-diagonal operator is the singlet A_K3.

The minimal bridge is also exact.  If the sterile recovery instrument carries
a pseudoscalar orientation pointer, i.e. a one-dimensional environment record
o with

    pi . o = sgn(pi) o,

then o K_R1 is an S3 scalar and the oriented boundary is readable relative to
that pointer.  Equivalently, the sign-twisted group average recovers K_R1.

Thus the theorem is:

    Sterile recovery reads the oriented R1 boundary iff the recovery
    environment contains a sign-representation orientation pointer, or an
    equivalent fixed physical orientation of the ordered (G0,G1) boundary.

Current-canon consequence:

    The existing sterile source algebra explicitly has no G0/G1 clause and no
    generation-labelled environment record.  It supplies a generation-singlet
    source port.  So the present canon refutes an unconditional K_or readout;
    it identifies the exact new service primitive needed to promote the CP
    latch from conditional to derived.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
GENS = ((0, 0), (0, 1), (1, 0))
GEN_INDEX = {g: i for i, g in enumerate(GENS)}
TOL = 1.0e-12


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def permutation_matrices() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3, dtype=complex)
    for perm in permutations(range(3)):
        inv_count = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        parity = -1 if inv_count % 2 else 1
        out.append((perm, eye[list(perm)], parity))
    return out


PERMS = permutation_matrices()


def p_cycle(order: tuple[tuple[int, int], tuple[int, int], tuple[int, int]]) -> np.ndarray:
    p = np.zeros((3, 3), dtype=complex)
    indices = [GEN_INDEX[g] for g in order]
    for src, dst in zip(indices, indices[1:] + indices[:1]):
        p[dst, src] = 1.0
    return p


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def k_r1() -> np.ndarray:
    """Positive R1 boundary convention: 00 -> 10 -> 01 -> 00."""

    order = ((0, 0), (1, 0), (0, 1))
    p = p_cycle(order)
    return p - p.T


def ordinary_average(M: np.ndarray) -> np.ndarray:
    out = np.zeros_like(M)
    for _perm, p, _parity in PERMS:
        out += p.T @ M @ p
    return out / len(PERMS)


def sign_twisted_average(M: np.ndarray) -> np.ndarray:
    out = np.zeros_like(M)
    for _perm, p, parity in PERMS:
        out += parity * (p.T @ M @ p)
    return out / len(PERMS)


def intertwiner_basis(
    domain_reps: list[np.ndarray],
    target_reps: list[np.ndarray],
) -> list[np.ndarray]:
    """Basis of T satisfying U_target(g) T = T U_domain(g)."""

    if len(domain_reps) != len(target_reps):
        raise ValueError("domain and target representations must have same group size")
    target_dim = target_reps[0].shape[0]
    domain_dim = domain_reps[0].shape[0]
    rows: list[np.ndarray] = []

    def unknown_index(row: int, col: int) -> int:
        return row * domain_dim + col

    for ud, ut in zip(domain_reps, target_reps):
        for i in range(target_dim):
            for j in range(domain_dim):
                row = np.zeros(target_dim * domain_dim, dtype=float)
                for r in range(target_dim):
                    row[unknown_index(r, j)] += float(np.real_if_close(ut[i, r]))
                for c in range(domain_dim):
                    row[unknown_index(i, c)] -= float(np.real_if_close(ud[c, j]))
                rows.append(row)

    constraints = np.vstack(rows)
    _u, svals, vh = np.linalg.svd(constraints)
    rank = int(np.sum(svals > TOL))
    null = vh[rank:]
    return [vec.reshape(target_dim, domain_dim) for vec in null]


def max_abs_column(basis: list[np.ndarray], col: int) -> float:
    if not basis:
        return 0.0
    return max(float(np.max(np.abs(b[:, col]))) for b in basis)


def main() -> None:
    print("ITEM 87/R15 -- STERILE R1-BOUNDARY READOUT THEOREM")
    print("=" * 92)

    require_text(
        "python_code/item123_sterile_generation_singlet_source.py",
        (
            "contains no G0/G1 clause",
            "generation-labeled environment record",
            "the source is one coherent port",
        ),
    )
    require_text(
        "python_code/item87_kor_r1_orientation_hunt.py",
        (
            "K_or is present as the oriented R1-boundary cochain",
            "It averages K_R1 to zero",
        ),
    )
    require_text(
        "python_code/axis_booking_resolution.py",
        (
            "Any future claim that re-identifies colour or",
            "generation with the TCH translation axes must be promoted as a new bridge",
        ),
    )

    A = a_k3()
    K = k_r1()
    P_plus = 0.5 * (A + K)

    print("\n[1] Transformation law")
    for perm, p, parity in PERMS:
        moved_a = p.T @ A @ p
        moved_k = p.T @ K @ p
        check(np.linalg.norm(moved_a - A) < TOL, f"A_K3 is invariant under {perm}")
        if parity == 1:
            check(np.linalg.norm(moved_k - K) < TOL, f"even relabelling {perm} preserves K_R1")
        else:
            check(np.linalg.norm(moved_k + K) < TOL, f"odd relabelling {perm} flips K_R1")

    print("\n[2] S3-blind averaging")
    print(f"  ||<A_K3>_S3 - A_K3||       = {np.linalg.norm(ordinary_average(A) - A):.3e}")
    print(f"  ||<K_R1>_S3||              = {np.linalg.norm(ordinary_average(K)):.3e}")
    print(f"  ||<P_+>_S3 - A_K3/2||      = {np.linalg.norm(ordinary_average(P_plus) - 0.5 * A):.3e}")
    check(np.linalg.norm(ordinary_average(A) - A) < TOL, "singlet support survives scalar recovery")
    check(np.linalg.norm(ordinary_average(K)) < TOL, "oriented boundary is erased by scalar recovery")
    check(np.linalg.norm(ordinary_average(P_plus) - 0.5 * A) < TOL, "directed latch averages to unoriented support")

    print("\n[3] Intertwiner no-go and minimal bridge")
    # Domain basis is (A_K3, K_R1): trivial plus sign representation.
    domain_reps = [np.diag([1.0, float(parity)]) for _perm, _p, parity in PERMS]
    trivial_reps = [np.array([[1.0]]) for _perm, _p, _parity in PERMS]
    sign_reps = [np.array([[float(parity)]]) for _perm, _p, parity in PERMS]
    scalar_plus_sign_reps = [np.diag([1.0, float(parity)]) for _perm, _p, parity in PERMS]

    trivial_basis = intertwiner_basis(domain_reps, trivial_reps)
    sign_basis = intertwiner_basis(domain_reps, sign_reps)
    combined_basis = intertwiner_basis(domain_reps, scalar_plus_sign_reps)

    print(f"  dim Hom_S3(trivial+sign, trivial pointer)     = {len(trivial_basis)}")
    print(f"  max K_R1 column in scalar-pointer maps        = {max_abs_column(trivial_basis, 1):.3e}")
    print(f"  dim Hom_S3(trivial+sign, sign pointer)        = {len(sign_basis)}")
    print(f"  max A_K3 column in sign-pointer maps          = {max_abs_column(sign_basis, 0):.3e}")
    print(f"  dim Hom_S3(trivial+sign, trivial+sign pointer)= {len(combined_basis)}")
    check(len(trivial_basis) == 1, "scalar environment has one equivariant readout channel")
    check(max_abs_column(trivial_basis, 1) < TOL, "scalar environment cannot read the sign component K_R1")
    check(len(sign_basis) == 1, "sign environment has one equivariant readout channel")
    check(max_abs_column(sign_basis, 0) < TOL, "sign environment reads the orientation component, not A_K3")
    check(len(combined_basis) == 2, "trivial plus sign pointer can read both A_K3 and K_R1")
    for b in combined_basis:
        check(abs(b[0, 1]) < TOL and abs(b[1, 0]) < TOL, "combined pointer does not mix scalar and sign channels")

    print("\n[4] Sign-pointer reconstruction")
    print(f"  ||<K_R1>_sign - K_R1||      = {np.linalg.norm(sign_twisted_average(K) - K):.3e}")
    print(f"  ||<A_K3>_sign||             = {np.linalg.norm(sign_twisted_average(A)):.3e}")
    print(f"  ||<P_+>_sign - K_R1/2||     = {np.linalg.norm(sign_twisted_average(P_plus) - 0.5 * K):.3e}")
    check(np.linalg.norm(sign_twisted_average(K) - K) < TOL, "sign-twisted pointer recovers K_R1")
    check(np.linalg.norm(sign_twisted_average(A)) < TOL, "sign-twisted pointer rejects the scalar support")
    check(np.linalg.norm(sign_twisted_average(P_plus) - 0.5 * K) < TOL, "directed latch decomposes into scalar plus sign readouts")

    print("\n[5] Decision")
    print(
        """
  Theorem found:
    The sterile recovery instrument can read the oriented R1 boundary exactly
    when its Stinespring environment carries a sign-representation orientation
    pointer, or when the ordered (G0,G1) R1 boundary is otherwise fixed as a
    physical orientation.  Then o*K_R1 is an S3 scalar readout.

  Current-canon no-go:
    The documented sterile source algebra is generation-blind, has no G0/G1
    clause, and has one coherent Stinespring source port.  Such a scalar
    environment reads A_K3 but provably averages K_R1 to zero.

  New primitive, stated minimally:
    Add or derive a sterile-recovery orientation pointer transforming as the
    sign representation of S3.  That is the exact service primitive required
    before the directed CP latch and its absolute orientation can be called
    derived rather than conditional.
"""
    )
    print("ALL ASSERTIONS PASSED -- orientation readout iff sign-pointer bridge; current scalar source refutes unconditional K_or readout.")


if __name__ == "__main__":
    main()
