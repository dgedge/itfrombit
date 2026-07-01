#!/usr/bin/env python3
"""Executable theorem: the smallest record-grammar machine is a Bell pair.

Purpose
-------
This is the first deliberately small "record grammar" example.  It cashes out
the language

    endpoint records + shared relational constraint + detector readout + reset

as an exact two-qubit calculation, without introducing any new physics.

Theorem checked by this script
------------------------------
For the Bell record |Phi+> = (|00> + |11>)/sqrt(2):

  1. The state is not a pair of local facts.  The local endpoint states are
     maximally mixed, while the joint stabilizers XX and ZZ are sharp records.

  2. What a detector can "hear" depends on the readout basis.  Z/Z and X/X
     readouts are perfectly correlated; Y/Y is perfectly anti-correlated; mixed
     Z/X readouts are uncorrelated.

  3. A remote non-selective readout cannot signal.  Measuring B in X, Y, Z, or a
     tilted basis does not change A's reduced record.

  4. The relational record is not a classical shared-bit model.  The same fixed
     state violates CHSH with 2*sqrt(2).

  5. The thermodynamic arrow is in reset, not in correlation.  A fair local
     detector record carries ln(2) nats of Shannon entropy; the two correlated
     Z-detectors have joint entropy ln(2), or 2 ln(2) if reset separately with
     no use of their correlation.

Interpretation
--------------
This is intentionally modest.  It does not derive the Born rule; it uses the
ordinary Hilbert/stabilizer calculus already reconstructed elsewhere in canon.
Its job is to demonstrate that the proposed "records as endpoint constraints"
language has an exact algebraic model:

    local endpoint = readable subsystem,
    relational support = stabilizer constraint,
    detector = projective readout channel,
    reset = irreversible service step.

That is the first executable bridge from metaphor to operator grammar.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*ops: np.ndarray) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def partial_trace_two_qubits(rho: np.ndarray, keep: str) -> np.ndarray:
    """Partial trace of a 2-qubit density matrix.

    keep = "A" returns Tr_B rho; keep = "B" returns Tr_A rho.
    """

    tensor = rho.reshape(2, 2, 2, 2)  # row_A,row_B,col_A,col_B
    if keep == "A":
        return np.einsum("abcb->ac", tensor)
    if keep == "B":
        return np.einsum("abad->bd", tensor)
    raise ValueError("keep must be 'A' or 'B'")


def projectors_for_observable(obs: np.ndarray) -> list[tuple[int, np.ndarray]]:
    """Projectors for a +/-1 Hermitian observable."""

    return [(+1, (I2 + obs) / 2.0), (-1, (I2 - obs) / 2.0)]


def joint_distribution(rho: np.ndarray, obs_a: np.ndarray, obs_b: np.ndarray) -> dict[tuple[int, int], float]:
    probs: dict[tuple[int, int], float] = {}
    for va, pa in projectors_for_observable(obs_a):
        for vb, pb in projectors_for_observable(obs_b):
            p = np.trace(kron(pa, pb) @ rho)
            probs[(va, vb)] = float(np.real_if_close(p))
    return probs


def corr_from_distribution(dist: dict[tuple[int, int], float]) -> float:
    return sum(va * vb * p for (va, vb), p in dist.items())


def marginal(dist: dict[tuple[int, int], float], side: str) -> dict[int, float]:
    idx = 0 if side == "A" else 1
    out = {+1: 0.0, -1: 0.0}
    for key, p in dist.items():
        out[key[idx]] += p
    return out


def shannon_nats(probs: list[float]) -> float:
    return -sum(p * math.log(p) for p in probs if p > 0)


def nonselective_measure_b(rho: np.ndarray, obs_b: np.ndarray) -> np.ndarray:
    out = np.zeros_like(rho)
    for _, pb in projectors_for_observable(obs_b):
        k = kron(I2, pb)
        out += k @ rho @ k
    return out


def unit_axis(theta: float) -> np.ndarray:
    """A tilted X/Z-plane Pauli observable cos(theta) Z + sin(theta) X."""

    return math.cos(theta) * Z + math.sin(theta) * X


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<48s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(f"{name}: {value} != {target}")


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<48s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Record-grammar Bell machine: executable theorem")
    print("=" * 72)

    phi_plus = np.zeros(4, dtype=complex)
    phi_plus[0] = 1 / math.sqrt(2)
    phi_plus[3] = 1 / math.sqrt(2)
    rho = density(phi_plus)

    print("\n[1] Endpoint records are locally blank; relational stabilizers are sharp")
    rho_a = partial_trace_two_qubits(rho, "A")
    rho_b = partial_trace_two_qubits(rho, "B")
    assert_matrix_close("Tr_B rho = I/2", rho_a, I2 / 2)
    assert_matrix_close("Tr_A rho = I/2", rho_b, I2 / 2)
    for label, op, target in [
        ("<X_A>", kron(X, I2), 0.0),
        ("<Z_A>", kron(Z, I2), 0.0),
        ("<X_A X_B>", kron(X, X), 1.0),
        ("<Z_A Z_B>", kron(Z, Z), 1.0),
        ("<Y_A Y_B>", kron(Y, Y), -1.0),
    ]:
        assert_close(label, expect(rho, op), target)

    print("\n[2] What can be heard depends on detector basis")
    readouts = [
        ("Z/Z", Z, Z, +1.0),
        ("X/X", X, X, +1.0),
        ("Y/Y", Y, Y, -1.0),
        ("Z/X", Z, X, 0.0),
    ]
    for label, oa, ob, corr_target in readouts:
        dist = joint_distribution(rho, oa, ob)
        corr = corr_from_distribution(dist)
        ma = marginal(dist, "A")
        mb = marginal(dist, "B")
        print(f"  {label:<4s} dist={{{', '.join(f'{k}:{v:.3f}' for k, v in sorted(dist.items()))}}}")
        assert_close(f"{label} correlation", corr, corr_target)
        assert_close(f"{label} A marginal +", ma[+1], 0.5)
        assert_close(f"{label} B marginal +", mb[+1], 0.5)

    print("\n[3] No signalling: remote non-selective readout cannot alter A")
    for label, obs in [
        ("B measures Z", Z),
        ("B measures X", X),
        ("B measures Y", Y),
        ("B measures tilted pi/7", unit_axis(math.pi / 7.0)),
    ]:
        rho_after = nonselective_measure_b(rho, obs)
        assert_matrix_close(label, partial_trace_two_qubits(rho_after, "A"), rho_a)

    print("\n[4] The relational record is not a classical shared-bit model: CHSH")
    # For |Phi+>, choose A0=Z, A1=X, B0=(Z+X)/sqrt(2), B1=(Z-X)/sqrt(2).
    a0 = Z
    a1 = X
    b0 = (Z + X) / math.sqrt(2)
    b1 = (Z - X) / math.sqrt(2)
    chsh = (
        expect(rho, kron(a0, b0))
        + expect(rho, kron(a0, b1))
        + expect(rho, kron(a1, b0))
        - expect(rho, kron(a1, b1))
    )
    assert_close("CHSH", chsh, 2.0 * math.sqrt(2.0), tol=1e-12)
    print("  classical local-record bound = 2; Bell record reaches 2*sqrt(2)")

    print("\n[5] Reset ledger: correlation is reversible; reset is the arrow")
    zz = joint_distribution(rho, Z, Z)
    h_a = shannon_nats([marginal(zz, "A")[+1], marginal(zz, "A")[-1]])
    h_b = shannon_nats([marginal(zz, "B")[+1], marginal(zz, "B")[-1]])
    h_ab = shannon_nats(list(zz.values()))
    mutual = h_a + h_b - h_ab
    assert_close("H(A)", h_a, math.log(2.0))
    assert_close("H(B)", h_b, math.log(2.0))
    assert_close("H(A,B) correlated pair", h_ab, math.log(2.0))
    assert_close("I(A:B)", mutual, math.log(2.0))
    print(f"  reset one local detector costs at least {h_a:.12g} nats")
    print(f"  reset both jointly with correlation costs {h_ab:.12g} nats")
    print(f"  reset both locally without using correlation costs {h_a + h_b:.12g} nats")

    print(
        """
Verdict
-------
The smallest record-grammar machine is exact:

  endpoint records are locally random;
  the physical content is the shared stabilizer constraint;
  detector basis selects which part of the relation is heard;
  remote readout cannot signal;
  the relational record violates the classical shared-bit bound;
  the irreversible service cost enters when finite detector records are reset.

This is not yet a new dynamics.  It is the first executable theorem showing
that "what can be said / heard / remembered" has a precise operator-algebra
meaning in the stabilizer record grammar.
"""
    )


if __name__ == "__main__":
    main()
