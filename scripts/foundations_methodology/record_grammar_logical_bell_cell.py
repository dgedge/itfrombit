#!/usr/bin/env python3
"""Logical Bell records inside the [8,4,4] byte cell.

This is the second executable record-grammar toy theorem.  The first script
(`record_grammar_bell_machine.py`) used abstract qubits.  Here the endpoint
records are embedded into the framework's actual byte cell:

    C = RM(1,3) = [8,4,4],

with a two-record slice |0_L> = |00000000>, |1_L> = |11111111>.  This is not
claimed to be the whole [[8,0,4]] syndrome-record cell as a free logical qubit:
the full CSS cell is a complete stabilizer record, not a quantum memory with
logical dimension two.  The purpose is narrower and explicit:

    Does the endpoint/relational/detector/reset story survive when the endpoints
    are byte-level records and single-site errors are checked by the byte's
    parity syndrome?

The answer is yes, at this pedagogical slice level:

  * the two endpoint byte records are locally blank in the logical Bell state;
  * the relational content is the logical stabilizer pair X_L X_L and Z_L Z_L;
  * detector readout in the logical Z and logical X bases behaves like the
    abstract Bell machine;
  * a single physical bit flip is detected by the 12 edge Z_i Z_j syndrome and
    identifies the flipped vertex;
  * reset entropy remains ln(2) nats for one fair logical detector record.

This is a bridge, not a new framework result: it shows how the "record grammar"
language can be represented in the existing [8,4,4] machinery without treating
entanglement as a spatial wire.

Implementation note: the Bell calculation is performed in the explicit
two-record slice {|0_L>, |1_L>} and the byte-level syndrome calculation is
performed on the actual eight physical bits.  We do not materialise the full
256 x 256 cell algebra or the 65,536-dimensional two-cell Hilbert space because
that would obscure the theorem and waste memory.
"""

from __future__ import annotations

import math

import numpy as np


N = 8
ZERO = 0
ONE = (1 << N) - 1

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def bits(state: int) -> tuple[int, ...]:
    return tuple((state >> i) & 1 for i in range(N))


def ket(index: int, dim: int = 2) -> np.ndarray:
    v = np.zeros(dim, dtype=complex)
    v[index] = 1.0
    return v


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def rm13_words() -> set[int]:
    # Same convention as minimal_balanced_record_cell_theorem.py:
    # rows are {1, b2, b1, b0}.
    gens = [
        0b11111111,
        0b11110000,
        0b11001100,
        0b10101010,
    ]
    out = {0}
    for g in gens:
        out |= {x ^ g for x in list(out)}
    return out


CODE = rm13_words()


def cube_edges() -> list[tuple[int, int]]:
    return [(u, v) for u in range(8) for v in range(u + 1, 8) if (u ^ v).bit_count() == 1]


EDGES = cube_edges()


def syndrome_word(state: int) -> tuple[int, ...]:
    b = bits(state)
    return tuple(b[u] ^ b[v] for u, v in EDGES)


ZERO_SYNDROME = syndrome_word(ZERO)
INCIDENT = {v: tuple(1 if v in e else 0 for e in EDGES) for v in range(8)}


def partial_trace_two_logical(rho: np.ndarray, keep: str) -> np.ndarray:
    """Trace a two-logical-record density matrix over one 2-state endpoint."""

    tensor = rho.reshape(2, 2, 2, 2)  # row_A,row_B,col_A,col_B
    if keep == "A":
        return np.einsum("abcb->ac", tensor)
    if keep == "B":
        return np.einsum("abad->bd", tensor)
    raise ValueError("keep must be A or B")


def projectors_z() -> list[tuple[int, np.ndarray]]:
    return [
        (+1, np.outer(ket(0), ket(0).conj())),
        (-1, np.outer(ket(1), ket(1).conj())),
    ]


def projectors_x() -> list[tuple[int, np.ndarray]]:
    plus = (ket(0) + ket(1)) / math.sqrt(2)
    minus = (ket(0) - ket(1)) / math.sqrt(2)
    return [(+1, np.outer(plus, plus.conj())), (-1, np.outer(minus, minus.conj()))]


def kron(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.kron(a, b)


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def joint_distribution(rho: np.ndarray, pa: list[tuple[int, np.ndarray]], pb: list[tuple[int, np.ndarray]]) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for va, p_a in pa:
        for vb, p_b in pb:
            out[(va, vb)] = float(np.real_if_close(np.trace(kron(p_a, p_b) @ rho)))
    return out


def corr(dist: dict[tuple[int, int], float]) -> float:
    return sum(a * b * p for (a, b), p in dist.items())


def shannon_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0)


def assert_close(name: str, value: float, target: float, tol: float = 1e-10) -> None:
    err = abs(value - target)
    print(f"  {name:<52s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-10) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<52s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def single_bit_flip(state: int, vertex: int) -> int:
    return state ^ (1 << vertex)


def main() -> None:
    print("Logical Bell pair in the [8,4,4] record cell")
    print("=" * 76)

    print("\n[1] Byte-cell checks")
    print(f"  |C| = {len(CODE)}; contains 0^8? {ZERO in CODE}; contains 1^8? {ONE in CODE}")
    print(f"  cube edge checks = {len(EDGES)}")
    assert len(CODE) == 16 and ZERO in CODE and ONE in CODE and len(EDGES) == 12

    # |Phi_L+> = (|0_L 0_L> + |1_L 1_L>)/sqrt(2), represented in the
    # two-record slice where basis index 0 means byte word 0^8 and basis index 1
    # means byte word 1^8.
    psi = np.zeros(4, dtype=complex)
    psi[0] = 1 / math.sqrt(2)
    psi[3] = 1 / math.sqrt(2)
    rho = density(psi)

    print("\n[2] Endpoint cells are locally blank inside the logical record slice")
    rho_a = partial_trace_two_logical(rho, "A")
    rho_b = partial_trace_two_logical(rho, "B")
    assert_matrix_close("Tr_B rho = I_L/2", rho_a, I2 / 2.0)
    assert_matrix_close("Tr_A rho = I_L/2", rho_b, I2 / 2.0)
    assert_close("<Z_L A>", expect(rho, kron(Z, I2)), 0.0)
    assert_close("<X_L A>", expect(rho, kron(X, I2)), 0.0)
    assert_close("<Z_L A Z_L B>", expect(rho, kron(Z, Z)), 1.0)
    assert_close("<X_L A X_L B>", expect(rho, kron(X, X)), 1.0)
    assert_close("<Y_L A Y_L B>", expect(rho, kron(Y, Y)), -1.0)

    print("\n[3] Logical detector readout")
    for label, pa, pb, target in [
        ("Z_L/Z_L", projectors_z(), projectors_z(), 1.0),
        ("X_L/X_L", projectors_x(), projectors_x(), 1.0),
        ("Z_L/X_L", projectors_z(), projectors_x(), 0.0),
    ]:
        dist = joint_distribution(rho, pa, pb)
        print(f"  {label:<7s} dist={{{', '.join(f'{k}:{v:.3f}' for k, v in sorted(dist.items()))}}}")
        assert_close(f"{label} correlation", corr(dist), target)

    print("\n[4] Single-site physical errors are detected and addressed by the 12 edge checks")
    assert syndrome_word(ZERO) == ZERO_SYNDROME
    assert syndrome_word(ONE) == ZERO_SYNDROME
    for base_label, base in [("0_L", ZERO), ("1_L", ONE)]:
        for vertex in range(8):
            bad = single_bit_flip(base, vertex)
            delta = tuple(a ^ b for a, b in zip(syndrome_word(base), syndrome_word(bad)))
            assert delta == INCIDENT[vertex]
            assert bad not in CODE
        print(f"  every one-bit flip of {base_label} leaves C and carries a unique incident-edge address")
    assert len(set(INCIDENT.values())) == 8
    print("  the error syndrome alphabet is exactly the eight vertex addresses")

    print("\n[5] Reset ledger at logical-detector level")
    zz = joint_distribution(rho, projectors_z(), projectors_z())
    p_a_plus = zz[(+1, +1)] + zz[(+1, -1)]
    p_a_minus = zz[(-1, +1)] + zz[(-1, -1)]
    h_a = shannon_nats([p_a_plus, p_a_minus])
    h_ab = shannon_nats(list(zz.values()))
    assert_close("H(A logical detector)", h_a, math.log(2.0))
    assert_close("H(A,B correlated logical detectors)", h_ab, math.log(2.0))
    print(f"  one fair logical detector reset costs at least ln2 = {math.log(2.0):.12g} nats")
    print("  one nat is natural-log information; one fair bit is ln2 nats")

    print(
        """
Verdict
-------
The abstract Bell-record theorem embeds cleanly into a two-record slice of the
[8,4,4] byte cell.  The endpoint facts are byte records, the relational content
is the logical stabilizer constraint, the detector hears a chosen compatible
logical basis, and the byte's 12 edge checks detect and address any one-site
physical error.

Boundary
--------
This is not a claim that the full [[8,0,4]] record cell stores a free qubit, and
it is not yet a Wilson-dressed or dynamical scattering example.  It is the first
framework-native bridge: generic Bell algebra -> byte-level endpoint records ->
syndrome-addressed QEC protection -> reset ledger.
"""
    )


if __name__ == "__main__":
    main()
