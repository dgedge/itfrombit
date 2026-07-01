#!/usr/bin/env python3
"""Byte-level Bell endpoint records with an explicit Wilson link.

This is the fourth small "record grammar" executable:

  1. record_grammar_bell_machine.py
       abstract Bell records;
  2. record_grammar_logical_bell_cell.py
       Bell records embedded in the [8,4,4] byte-cell slice;
  3. record_grammar_wilson_dressed_pair.py
       a two-endpoint Wilson-dressed readout;
  4. this script
       a byte-level Bell endpoint pair whose phase-sensitive readout requires a
       Wilson link.

The state used here is the one-excitation Bell record

    |Psi_L(delta)> = ( |1_L 0_L> + exp(i delta) |0_L 1_L> ) / sqrt(2),

where |0_L> is the byte word 00000000 and |1_L> is the byte word 11111111
inside the RM(1,3) = [8,4,4] record cell.

Why this state rather than |Phi+>?
----------------------------------
The |Phi+> state is good for showing pure relational stabilizers.  A Wilson link
is most transparent when there is a charged record that can sit at endpoint A
or endpoint B.  Then the bare exchange coherence

    |1_L 0_L><0_L 1_L| + h.c.

is not gauge-readable: local U(1) phases at A and B change it.  The detector
can hear the phase-sensitive exchange only as the dressed operator

    X_U = U_AB |1_L 0_L><0_L 1_L| + U_AB^* |0_L 1_L><1_L 0_L|.

That is the smallest exact example of "endpoint records plus Wilson dressing".

Boundary
--------
This is still a pedagogical finite theorem, not a full lattice gauge theory.  It
does not derive Maxwell dynamics, charges, or scattering.  It only verifies the
operator grammar:

    byte endpoint records are protected by [8,4,4] syndrome checks;
    bare endpoint coherence is gauge-dependent;
    the Wilson-dressed endpoint relation is gauge-invariant and detector-readable.
"""

from __future__ import annotations

import math

import numpy as np


# ---------------------------------------------------------------------------
# Byte-cell data: RM(1,3) = [8,4,4] record words and the 12 cube-edge checks.
# These are the framework-native pieces.  The logical two-state slice below is
# deliberately only {|0_L>, |1_L>} inside the 16-word byte alphabet.
# ---------------------------------------------------------------------------

N_BITS = 8
ZERO_WORD = 0
ONE_WORD = (1 << N_BITS) - 1


def rm13_words() -> set[int]:
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


def bits(state: int) -> tuple[int, ...]:
    return tuple((state >> i) & 1 for i in range(N_BITS))


def cube_edges() -> list[tuple[int, int]]:
    return [
        (u, v)
        for u in range(N_BITS)
        for v in range(u + 1, N_BITS)
        if (u ^ v).bit_count() == 1
    ]


CODE = rm13_words()
EDGES = cube_edges()


def syndrome_word(state: int) -> tuple[int, ...]:
    b = bits(state)
    return tuple(b[u] ^ b[v] for u, v in EDGES)


ZERO_SYNDROME = syndrome_word(ZERO_WORD)
INCIDENT = {v: tuple(1 if v in e else 0 for e in EDGES) for v in range(N_BITS)}


# ---------------------------------------------------------------------------
# Logical two-endpoint Hilbert space.
#
# Basis order:
#   0 = |0_L 0_L>
#   1 = |0_L 1_L>
#   2 = |1_L 0_L>
#   3 = |1_L 1_L>
#
# The local U(1) toy gauge says the charged logical record |1_L> at endpoint A
# picks up exp(i lambda_A), while |1_L> at endpoint B picks up exp(i lambda_B).
# ---------------------------------------------------------------------------

DIM = 4
I4 = np.eye(DIM, dtype=complex)
I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def ket(index: int) -> np.ndarray:
    v = np.zeros(DIM, dtype=complex)
    v[index] = 1.0
    return v


K00 = ket(0)
K01 = ket(1)
K10 = ket(2)
K11 = ket(3)


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def psi_one_excitation(delta: float) -> np.ndarray:
    """(|1_L0_L> + exp(i delta)|0_L1_L>)/sqrt(2)."""

    return (K10 + np.exp(1j * delta) * K01) / math.sqrt(2.0)


def gauge_unitary(lambda_a: float, lambda_b: float) -> np.ndarray:
    """Local U(1) phase on charged logical records at A and B."""

    phases = [
        1.0,
        np.exp(1j * lambda_b),
        np.exp(1j * lambda_a),
        np.exp(1j * (lambda_a + lambda_b)),
    ]
    return np.diag(phases).astype(complex)


def link(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def transform_link(u_ab: complex, lambda_a: float, lambda_b: float) -> complex:
    """Wilson link orientation AB: U_AB -> exp(i(lambda_A-lambda_B)) U_AB."""

    return complex(np.exp(1j * (lambda_a - lambda_b)) * u_ab)


# Bare exchange coherence between the two one-excitation records.
T_AB = np.outer(K10, K01.conj())  # |1_L0_L><0_L1_L|
T_BA = np.outer(K01, K10.conj())  # |0_L1_L><1_L0_L|
X_BARE = T_AB + T_BA


def dressed_exchange(u_ab: complex) -> np.ndarray:
    """Gauge-dressed exchange readout X_U."""

    return u_ab * T_AB + np.conj(u_ab) * T_BA


def expect(psi: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.vdot(psi, op @ psi)))


def partial_trace_two_logical(rho: np.ndarray, keep: str) -> np.ndarray:
    tensor = rho.reshape(2, 2, 2, 2)  # row_A,row_B,col_A,col_B
    if keep == "A":
        return np.einsum("abcb->ac", tensor)
    if keep == "B":
        return np.einsum("abad->bd", tensor)
    raise ValueError("keep must be A or B")


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<58s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<58s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def shannon_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0)


def group_average_open_phase(n: int = 64) -> complex:
    """Average exp(i(lambda_A-lambda_B)) over a finite U(1)xU(1) grid."""

    grid = [2.0 * math.pi * k / n for k in range(n)]
    return complex(sum(np.exp(1j * (la - lb)) for la in grid for lb in grid) / (n * n))


def main() -> None:
    print("Byte-level Bell endpoint records with a Wilson link")
    print("=" * 82)

    print("\n[1] Byte-cell substrate checks")
    print(f"  RM(1,3) code size = {len(CODE)}; contains 0^8 and 1^8 = {ZERO_WORD in CODE and ONE_WORD in CODE}")
    print(f"  edge-check count = {len(EDGES)}")
    assert len(CODE) == 16
    assert ZERO_WORD in CODE and ONE_WORD in CODE
    assert len(EDGES) == 12

    print("\n[2] One-physical-bit errors are still byte-syndrome-addressed")
    assert syndrome_word(ZERO_WORD) == ZERO_SYNDROME
    assert syndrome_word(ONE_WORD) == ZERO_SYNDROME
    for label, base in [("0_L", ZERO_WORD), ("1_L", ONE_WORD)]:
        for vertex in range(N_BITS):
            bad = base ^ (1 << vertex)
            delta_s = tuple(a ^ b for a, b in zip(syndrome_word(base), syndrome_word(bad)))
            assert delta_s == INCIDENT[vertex]
            assert bad not in CODE
        print(f"  every one-bit flip of {label} leaves C and gives a unique vertex-address syndrome")
    assert len(set(INCIDENT.values())) == 8

    print("\n[3] One-excitation Bell record: local endpoints blank, relation sharp")
    delta = 0.41
    theta = -0.76
    psi = psi_one_excitation(delta)
    rho = density(psi)
    rho_a = partial_trace_two_logical(rho, "A")
    rho_b = partial_trace_two_logical(rho, "B")
    assert_matrix_close("Tr_B rho = I_L/2", rho_a, I2 / 2.0)
    assert_matrix_close("Tr_A rho = I_L/2", rho_b, I2 / 2.0)
    z_a_z_b = np.kron(Z, Z)
    assert_close("<Z_A Z_B> one-excitation anticorrelation", expect(psi, z_a_z_b), -1.0)
    assert_close("bare exchange <X_AB>", expect(psi, X_BARE), math.cos(delta))

    print("\n[4] Bare exchange is not a gauge-readable record")
    lambda_a = 1.31
    lambda_b = -0.52
    g = gauge_unitary(lambda_a, lambda_b)
    psi_g = g @ psi
    bare = expect(psi, X_BARE)
    bare_g = expect(psi_g, X_BARE)
    assert_close("bare before gauge", bare, math.cos(delta))
    assert_close("bare after gauge", bare_g, math.cos(delta + lambda_b - lambda_a))
    if abs(bare - bare_g) < 1e-6:
        raise AssertionError("chosen gauge transformation accidentally preserved bare exchange")
    avg = group_average_open_phase()
    print(f"  Gauss average of the open exchange phase = {avg.real:+.3e}{avg.imag:+.3e}i")
    if abs(avg) > 1e-12:
        raise AssertionError("open exchange phase survived Gauss averaging")
    print("  -> the relative endpoint phase exists in the state, but a detector cannot hear it naked.")

    print("\n[5] Wilson-dressed exchange is gauge-invariant")
    u = link(theta)
    u_g = transform_link(u, lambda_a, lambda_b)
    x_u = dressed_exchange(u)
    x_ug = dressed_exchange(u_g)
    assert_matrix_close("X(U') = G X(U) G^dag", x_ug, g @ x_u @ g.conj().T)
    dressed = expect(psi, x_u)
    dressed_g = expect(psi_g, x_ug)
    assert_close("dressed exchange before gauge", dressed, math.cos(delta + theta))
    assert_close("dressed exchange after gauge", dressed_g, dressed)
    for la, lb in [(0.0, 0.0), (0.2, -0.9), (2.4, 0.1), (-1.7, 3.3)]:
        gg = gauge_unitary(la, lb)
        assert_close(
            f"dressed invariant at ({la:.1f},{lb:.1f})",
            expect(gg @ psi, dressed_exchange(transform_link(u, la, lb))),
            dressed,
        )
    print("  -> the hearable exchange record is endpoint relation plus Wilson link.")

    print("\n[6] Reset ledger remains a logical detector fact")
    # A local occupancy detector on endpoint A has probabilities 1/2,1/2.
    h_a = shannon_nats([0.5, 0.5])
    assert_close("H(A occupancy detector)", h_a, math.log(2.0))
    print("  one fair logical detector bit is ln2 nats; reset, not correlation, is the irreversible step")

    print(
        """
Verdict
-------
The byte-level Bell endpoint pair now carries an explicit Wilson link.

  The nouns are byte records: |0_L>=00000000 and |1_L>=11111111 inside the
  [8,4,4] cell.  The code's 12 edge checks still detect and address any one-bit
  physical error.  The verb "compare the endpoint phase" is not gauge-readable
  by itself: the bare exchange coherence changes under local endpoint phases and
  vanishes under Gauss averaging.  The grammatical, detector-readable sentence is

      endpoint exchange + Wilson holonomy.

This is the first combined record/QEC/gauge toy: endpoint records, relational
entanglement, Wilson dressing, and reset ledger in one finite operator model.
"""
    )


if __name__ == "__main__":
    main()
