#!/usr/bin/env python3
"""GHZ record machine: a three-endpoint relational constraint.

The Bell scripts show that two endpoints can carry a relational record that is
not reducible to two local facts.  The GHZ state shows a stronger point:

    some record content is genuinely three-endpoint content.

For

    |GHZ(delta)> = ( |000> + exp(i delta)|111> ) / sqrt(2),

every one-endpoint reduced state is maximally mixed and every two-endpoint
reduced state is the same classical correlated mixture,

    rho_AB = (|00><00| + |11><11|)/2,

independent of the phase delta.  The phase is not stored in any single endpoint
or in any pair of endpoints.  It is stored only in the three-way relation:

    <XXX> = cos(delta),        <YYY> = -sin(delta).

At delta = 0 this becomes the standard GHZ/Mermin stabilizer record:

    XXX = +1,
    XYY = YXY = YYX = -1,

which cannot be reproduced by assigning pre-existing local +/-1 values to X and
Y at the three endpoints.

Record-grammar reading
----------------------
This is the first example where "what can be heard" is not merely basis
dependent but arity dependent:

    one endpoint hears no fact;
    two endpoints hear only the classical parity/correlation shell;
    three endpoints can hear the phase-sensitive relational record.

The script uses only ordinary finite Hilbert-space algebra.  It is a clean toy
theorem, not a new dynamics.
"""

from __future__ import annotations

import math
from itertools import product

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_all(*ops: np.ndarray) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def ket(index: int, dim: int = 8) -> np.ndarray:
    v = np.zeros(dim, dtype=complex)
    v[index] = 1.0
    return v


def ghz(delta: float) -> np.ndarray:
    return (ket(0) + np.exp(1j * delta) * ket(7)) / math.sqrt(2.0)


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def index_from_bits(bits: tuple[int, ...]) -> int:
    out = 0
    for bit in bits:
        out = (out << 1) | bit
    return out


def bits_from_index(index: int, n: int) -> tuple[int, ...]:
    return tuple((index >> (n - 1 - i)) & 1 for i in range(n))


def partial_trace(rho: np.ndarray, keep: tuple[int, ...], n: int = 3) -> np.ndarray:
    """Small explicit partial trace for n qubits.

    `keep` lists endpoint indices retained in the output, in that order.
    The implementation uses loops on purpose: for n=3 it is clearer than a
    clever einsum, and the note is pedagogical.
    """

    keep = tuple(keep)
    trace = tuple(i for i in range(n) if i not in keep)
    out_dim = 2 ** len(keep)
    out = np.zeros((out_dim, out_dim), dtype=complex)

    for kept_row in product((0, 1), repeat=len(keep)):
        for kept_col in product((0, 1), repeat=len(keep)):
            total = 0.0 + 0.0j
            for traced_bits in product((0, 1), repeat=len(trace)):
                row = [0] * n
                col = [0] * n
                for pos, bit in zip(keep, kept_row):
                    row[pos] = bit
                for pos, bit in zip(keep, kept_col):
                    col[pos] = bit
                for pos, bit in zip(trace, traced_bits):
                    row[pos] = bit
                    col[pos] = bit
                total += rho[index_from_bits(tuple(row)), index_from_bits(tuple(col))]
            out[index_from_bits(kept_row), index_from_bits(kept_col)] = total
    return out


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def projectors(obs: np.ndarray) -> list[tuple[int, np.ndarray]]:
    return [(+1, (I2 + obs) / 2.0), (-1, (I2 - obs) / 2.0)]


def nonselective_measure(rho: np.ndarray, endpoint: int, obs: np.ndarray) -> np.ndarray:
    out = np.zeros_like(rho)
    for _, p in projectors(obs):
        ops = [I2, I2, I2]
        ops[endpoint] = p
        k = kron_all(*ops)
        out += k @ rho @ k
    return out


def zzz_distribution(rho: np.ndarray) -> dict[tuple[int, int, int], float]:
    out: dict[tuple[int, int, int], float] = {}
    for values in product((+1, -1), repeat=3):
        ops = []
        for value in values:
            ops.append((I2 + value * Z) / 2.0)
        out[values] = float(np.real_if_close(np.trace(kron_all(*ops) @ rho)))
    return out


def shannon_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0)


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


def main() -> None:
    print("GHZ record machine: three-endpoint relational constraint")
    print("=" * 82)

    delta = 0.73
    psi = ghz(delta)
    rho = density(psi)
    rho0 = density(ghz(0.0))
    rho_shifted = density(ghz(delta + 0.41))

    print("\n[1] One-endpoint records are blank")
    for endpoint in range(3):
        assert_matrix_close(f"rho_{endpoint} = I/2", partial_trace(rho, (endpoint,)), I2 / 2.0)
        assert_close(f"<Z_{endpoint}>", expect(rho, kron_all(*(Z if i == endpoint else I2 for i in range(3)))), 0.0)

    print("\n[2] Two-endpoint records hear only the classical parity shell")
    pair_shell = np.diag([0.5, 0.0, 0.0, 0.5]).astype(complex)
    for pair in [(0, 1), (0, 2), (1, 2)]:
        assert_matrix_close(f"rho_{pair} = (|00><00|+|11><11|)/2", partial_trace(rho, pair), pair_shell)
    for label, op in [
        ("<Z Z I>", kron_all(Z, Z, I2)),
        ("<Z I Z>", kron_all(Z, I2, Z)),
        ("<I Z Z>", kron_all(I2, Z, Z)),
    ]:
        assert_close(label, expect(rho, op), 1.0)

    print("\n[3] The GHZ phase is invisible to every proper subset")
    for endpoint in range(3):
        assert_matrix_close(
            f"one-endpoint phase-blind {endpoint}",
            partial_trace(rho, (endpoint,)),
            partial_trace(rho_shifted, (endpoint,)),
        )
    for pair in [(0, 1), (0, 2), (1, 2)]:
        assert_matrix_close(
            f"two-endpoint phase-blind {pair}",
            partial_trace(rho, pair),
            partial_trace(rho_shifted, pair),
        )
    print("  -> changing delta changes no one- or two-endpoint readable record.")

    print("\n[4] The phase is a genuine three-endpoint record")
    assert_close("<X X X>", expect(rho, kron_all(X, X, X)), math.cos(delta))
    assert_close("<Y Y Y>", expect(rho, kron_all(Y, Y, Y)), -math.sin(delta))
    assert_close("<X Y Y>", expect(rho, kron_all(X, Y, Y)), -math.cos(delta))
    assert_close("<Y X Y>", expect(rho, kron_all(Y, X, Y)), -math.cos(delta))
    assert_close("<Y Y X>", expect(rho, kron_all(Y, Y, X)), -math.cos(delta))

    print("\n[5] GHZ/Mermin contradiction at delta=0")
    values = {
        "XXX": expect(rho0, kron_all(X, X, X)),
        "XYY": expect(rho0, kron_all(X, Y, Y)),
        "YXY": expect(rho0, kron_all(Y, X, Y)),
        "YYX": expect(rho0, kron_all(Y, Y, X)),
    }
    assert_close("XXX", values["XXX"], +1.0)
    assert_close("XYY", values["XYY"], -1.0)
    assert_close("YXY", values["YXY"], -1.0)
    assert_close("YYX", values["YYX"], -1.0)
    product_of_three_negative_equations = values["XYY"] * values["YXY"] * values["YYX"]
    print(f"  product(XYY,YXY,YYX) = {product_of_three_negative_equations:+.0f}, but XXX = {values['XXX']:+.0f}")
    print("  -> predetermined local X/Y values would force these to agree; the GHZ record is genuinely nonclassical.")

    print("\n[6] Remote readout still cannot signal to a retained pair")
    rho_ab = partial_trace(rho, (0, 1))
    for label, obs in [("C measures Z", Z), ("C measures X", X), ("C measures Y", Y)]:
        after = nonselective_measure(rho, endpoint=2, obs=obs)
        assert_matrix_close(label, partial_trace(after, (0, 1)), rho_ab)

    print("\n[7] Reset ledger")
    dist = zzz_distribution(rho)
    nonzero = {k: v for k, v in dist.items() if v > 1e-12}
    print(f"  ZZZ distribution nonzero outcomes = {nonzero}")
    h_one = shannon_nats([0.5, 0.5])
    h_joint = shannon_nats(list(dist.values()))
    assert_close("H(one local Z detector)", h_one, math.log(2.0))
    assert_close("H(three correlated Z detectors)", h_joint, math.log(2.0))
    print(f"  resetting three independent local devices without using correlation can cost {3 * h_one:.12g} nats")
    print(f"  resetting the joint correlated ZZZ record costs {h_joint:.12g} nats")

    print(
        """
Verdict
-------
The GHZ record is an exact three-endpoint lesson:

  no endpoint contains a local fact;
  no pair contains the phase;
  the phase is carried by the three-way operator record;
  the Mermin stabilizers give a nonclassical all-or-nothing constraint;
  remote non-selective readout still cannot signal;
  reset cost belongs to detector records, not to the existence of correlation.

This is the arity version of the record grammar: what can be heard depends not
only on basis and gauge dressing, but on how many endpoints the detector record
is allowed to jointly address.
"""
    )


if __name__ == "__main__":
    main()
