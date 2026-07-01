#!/usr/bin/env python3
"""Magic-state warning: stabilizer records are not the whole quantum state.

This script is a guardrail for the record-grammar programme.

The earlier examples use stabilizers heavily because stabilizers are the
natural language of protected records and QEC syndromes.  But stabilizer
bookkeeping alone is not universal quantum mechanics.  It misses non-Clifford
phase resources, usually called "magic" in quantum information.

Toy model
---------
Start from the stabilizer state

    |+> = (|0> + |1>)/sqrt(2).

Apply the T gate

    T = diag(1, exp(i pi/4)),

to obtain

    |T> = (|0> + exp(i pi/4)|1>)/sqrt(2).

This state has Bloch vector

    (<X>, <Y>, <Z>) = (1/sqrt(2), 1/sqrt(2), 0).

It is not a Pauli stabilizer state: no Pauli expectation is +/-1.  It is also
outside the convex hull of the six one-qubit stabilizer states, because the
single-qubit stabilizer polytope is the octahedron

    |x| + |y| + |z| <= 1,

whereas |T> has |x|+|y|+|z| = sqrt(2).

Record-grammar reading
----------------------
Stabilizers are the protected record skeleton.  Magic is the non-stabilizer
phase resource carried on that skeleton.  A record grammar that silently
projects everything to stabilizers would lose detector-readable information.

This is a one-qubit warning example, not a theory of magic-state distillation.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {"X": X, "Y": Y, "Z": Z}


def ket(index: int) -> np.ndarray:
    v = np.zeros(2, dtype=complex)
    v[index] = 1.0
    return v


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def projector(psi: np.ndarray) -> np.ndarray:
    return density(psi / np.linalg.norm(psi))


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def bloch(rho: np.ndarray) -> np.ndarray:
    return np.array([expect(rho, X), expect(rho, Y), expect(rho, Z)])


def fidelity_pure(rho: np.ndarray, psi: np.ndarray) -> float:
    return float(np.real_if_close(np.vdot(psi, rho @ psi)))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<70s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<70s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def stabilizer_states() -> dict[str, np.ndarray]:
    zero = ket(0)
    one = ket(1)
    plus = (zero + one) / math.sqrt(2.0)
    minus = (zero - one) / math.sqrt(2.0)
    plus_i = (zero + 1j * one) / math.sqrt(2.0)
    minus_i = (zero - 1j * one) / math.sqrt(2.0)
    return {
        "|0>": zero,
        "|1>": one,
        "|+>": plus,
        "|->": minus,
        "|+i>": plus_i,
        "|-i>": minus_i,
    }


def t_gate() -> np.ndarray:
    return np.diag([1.0, np.exp(1j * math.pi / 4.0)])


def main() -> None:
    print("Magic-state warning: stabilizer skeleton is not magic-complete")
    print("=" * 88)

    zero = ket(0)
    one = ket(1)
    plus = (zero + one) / math.sqrt(2.0)
    t = t_gate()
    psi_t = t @ plus
    rho_t = density(psi_t)
    p_t = projector(psi_t)

    print("\n[1] T is non-Clifford: it does not map Pauli X to a Pauli")
    tx = t @ X @ t.conj().T
    target_tx = (X + Y) / math.sqrt(2.0)
    assert_matrix_close("T X T^dag = (X+Y)/sqrt(2)", tx, target_tx)
    distances = {name: float(np.linalg.norm(tx - op)) for name, op in PAULIS.items()}
    for name, dist in distances.items():
        print(f"  distance to {name}: {dist:.12g}")
    if min(distances.values()) < 1e-6:
        raise AssertionError("T accidentally mapped X to a Pauli")
    print("  -> Clifford summaries are not closed under the T operation.")

    print("\n[2] The T state sits between Pauli stabilizer axes")
    r_t = bloch(rho_t)
    assert_close("<X>_T", r_t[0], 1.0 / math.sqrt(2.0))
    assert_close("<Y>_T", r_t[1], 1.0 / math.sqrt(2.0))
    assert_close("<Z>_T", r_t[2], 0.0)
    max_pauli = float(np.max(np.abs(r_t)))
    assert_close("max Pauli expectation", max_pauli, 1.0 / math.sqrt(2.0))
    if max_pauli >= 1.0 - 1e-12:
        raise AssertionError("T state was accidentally a Pauli stabilizer state")
    print("  -> no Pauli stabilizer has a definite +/-1 value on |T>.")

    print("\n[3] The T state is outside the stabilizer octahedron")
    l1 = float(np.sum(np.abs(r_t)))
    assert_close("|x|+|y|+|z|", l1, math.sqrt(2.0))
    violation = l1 - 1.0
    assert_close("stabilizer-polytope violation", violation, math.sqrt(2.0) - 1.0)
    print("  -> |T> is not even a probabilistic mixture of one-qubit stabilizer states.")

    print("\n[4] Stabilizer replacements lose detector-readable phase information")
    stabs = stabilizer_states()
    fidelities = {label: fidelity_pure(projector(psi), psi_t) for label, psi in stabs.items()}
    for label, fid in sorted(fidelities.items(), key=lambda item: -item[1]):
        print(f"  magic-detector probability for {label:<4s}: {fid:.12g}")
    best_label, best_fid = max(fidelities.items(), key=lambda item: item[1])
    expected_best = (1.0 + 1.0 / math.sqrt(2.0)) / 2.0
    assert_close("best pure stabilizer overlap", best_fid, expected_best)
    print(f"  nearest pure stabilizer replacement: {best_label}, still misses {1.0 - best_fid:.6f}")

    # The best face projection of the Bloch vector onto the stabilizer octahedron
    # along the same X/Y direction is the mixed stabilizer state
    #   1/2 |+><+| + 1/2 |+i><+i|,
    # with Bloch vector (1/2, 1/2, 0).  It keeps the direction but shrinks the
    # magic vector back to the stabilizer boundary.
    rho_face = 0.5 * projector(stabs["|+>"]) + 0.5 * projector(stabs["|+i>"])
    r_face = bloch(rho_face)
    assert_matrix_close("best face stabilizer mixture Bloch vector", r_face, np.array([0.5, 0.5, 0.0]))
    face_magic_prob = float(np.real_if_close(np.trace(p_t @ rho_face)))
    assert_close("magic detector on face-projected stabilizer mixture", face_magic_prob, expected_best)
    print("  -> even the natural stabilizer-mixture surrogate loses the same T-basis signal.")

    print("\n[5] A rotated detector hears exactly what the stabilizer skeleton loses")
    p_magic_on_magic = float(np.real_if_close(np.trace(p_t @ rho_t)))
    p_magic_on_plus = float(np.real_if_close(np.trace(p_t @ projector(stabs["|+>"]))))
    p_magic_on_plusi = float(np.real_if_close(np.trace(p_t @ projector(stabs["|+i>"]))))
    assert_close("Pr(T-detector | T-state)", p_magic_on_magic, 1.0)
    assert_close("Pr(T-detector | |+>)", p_magic_on_plus, expected_best)
    assert_close("Pr(T-detector | |+i>)", p_magic_on_plusi, expected_best)
    print("  -> the non-Clifford phase is not metaphysical: the right detector basis reads it.")

    print("\n[6] Stabilizer measurements alone underdetermine the magic resource")
    # X and Y Pauli readout probabilities are ordinary detector facts.  They
    # show partial information, but neither stabilizer axis captures the whole
    # phase.  The T detector is the aligned readout.
    px_plus = (1.0 + r_t[0]) / 2.0
    py_plus = (1.0 + r_t[1]) / 2.0
    assert_close("Pr(X=+1 on |T>)", px_plus, expected_best)
    assert_close("Pr(Y=+1 on |T>)", py_plus, expected_best)
    print("  -> Pauli readouts see shadows of the phase; the magic projector sees the record exactly.")

    print(
        """
Verdict
-------
The stabilizer/QEC skeleton is essential, but it is not magic-complete:

  T maps a Pauli to a non-Pauli axis, so it is non-Clifford;
  |T> has no definite Pauli stabilizer;
  |T> lies outside the one-qubit stabilizer octahedron;
  stabilizer replacements lose at least 14.6% in the aligned T detector;
  a rotated non-Clifford detector reads the missing phase exactly.

Record-grammar rule:

  keep stabilizers as the protected record skeleton,
  but track non-Clifford/magic phases as named resources.

Discarding magic is not a harmless compression.  It destroys detector-readable
quantum information.
"""
    )


if __name__ == "__main__":
    main()
