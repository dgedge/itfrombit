#!/usr/bin/env python3
r"""Finite record-grammar compression certificate.

This script is the first integrated demonstration rather than a single toy
identity.  It builds a small full Hilbert-space state, chooses an admissible
record algebra, projects the full state onto that algebra, and certifies what
is preserved and what is deliberately lost.

Hilbert space
-------------
Use four qubits:

    A,B : two endpoint qubits with Bell/stabilizer relational records;
    R   : one route/loop qubit with a closed T holonomy record;
    D   : one finite detector branch bit.

The record basis is the tensor product

    Bell_AB basis  x  T-loop_R basis  x  detector-Z_D basis.

The record projection is the conditional expectation

    E_rec(rho) = sum_i P_i rho P_i,

where P_i are the orthogonal one-dimensional projectors onto that joint record
basis.  Equivalently, it is dephasing in the record basis.

Certificate
-----------
The script verifies:

  1. the record projectors are orthogonal and complete;
  2. E_rec is trace preserving, positive on the test state, and idempotent;
  3. every observable in the record algebra is preserved exactly;
  4. off-record witnesses are changed, but the change is bounded by
     ||rho - E_rec(rho)||_1;
  5. the trace-norm bound is saturated by the sign witness, so the residual is
     a real, measured object, not swept under the rug.

Record-grammar reading
----------------------
This is the executable version of the compression claim:

    full Hilbert state
      -> protected/stable record algebra
      -> exact preservation of detector-readable records
      -> explicit residual certificate for what was dropped.

The point is not that all correlations are small.  The point is that every
omitted correlation must be either outside the admissible detector algebra or
bounded by a measured residual.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def ket(index: int, dim: int = 2) -> np.ndarray:
    v = np.zeros(dim, dtype=complex)
    v[index] = 1.0
    return v


def dagger(a: np.ndarray) -> np.ndarray:
    return a.conj().T


def density(psi: np.ndarray) -> np.ndarray:
    psi = psi / np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def kron(*ops: np.ndarray) -> np.ndarray:
    out = np.array([[1.0]], dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def trace_norm(a: np.ndarray) -> float:
    return float(np.sum(np.linalg.svd(a, compute_uv=False)))


def operator_norm(a: np.ndarray) -> float:
    return float(np.linalg.svd(a, compute_uv=False)[0])


def bell_basis() -> list[tuple[str, np.ndarray]]:
    z0 = ket(0)
    z1 = ket(1)
    return [
        ("Phi+", (np.kron(z0, z0) + np.kron(z1, z1)) / math.sqrt(2.0)),
        ("Phi-", (np.kron(z0, z0) - np.kron(z1, z1)) / math.sqrt(2.0)),
        ("Psi+", (np.kron(z0, z1) + np.kron(z1, z0)) / math.sqrt(2.0)),
        ("Psi-", (np.kron(z0, z1) - np.kron(z1, z0)) / math.sqrt(2.0)),
    ]


def t_loop_basis() -> list[tuple[str, np.ndarray]]:
    """Closed-loop T basis on the route qubit."""

    z0 = ket(0)
    z1 = ket(1)
    phi = math.pi / 4.0
    return [
        ("T+", (z0 + phase(phi) * z1) / math.sqrt(2.0)),
        ("T-", (z0 - phase(phi) * z1) / math.sqrt(2.0)),
    ]


def detector_basis() -> list[tuple[str, np.ndarray]]:
    return [("D0", ket(0)), ("D1", ket(1))]


def record_unitary() -> tuple[np.ndarray, list[str]]:
    """Columns are Bell_AB x T_R x Z_D record basis vectors."""

    cols: list[np.ndarray] = []
    labels: list[str] = []
    for bell_label, bell in bell_basis():
        for t_label, route in t_loop_basis():
            for d_label, det in detector_basis():
                cols.append(np.kron(np.kron(bell, route), det))
                labels.append(f"{bell_label}:{t_label}:{d_label}")
    return np.column_stack(cols), labels


def conditional_expectation(rho: np.ndarray, u_record: np.ndarray) -> np.ndarray:
    """Project rho onto the diagonal algebra in the record basis."""

    sigma = dagger(u_record) @ rho @ u_record
    return u_record @ np.diag(np.diag(sigma)) @ dagger(u_record)


def record_observable(weights: np.ndarray, u_record: np.ndarray) -> np.ndarray:
    """An arbitrary observable in the commutative record algebra."""

    return u_record @ np.diag(weights.astype(complex)) @ dagger(u_record)


def record_projectors(u_record: np.ndarray) -> list[np.ndarray]:
    return [density(u_record[:, i]) for i in range(u_record.shape[1])]


def test_state(u_record: np.ndarray) -> np.ndarray:
    """A positive full state with small but real off-record coherences.

    In the record basis, the state is mostly diagonal plus a coherent
    off-record perturbation.  The convex construction guarantees positivity.
    """

    dim = u_record.shape[0]
    indices = np.arange(dim, dtype=float)
    base = 1.0 + 0.35 * np.sin(1.7 * indices) + 0.22 * np.cos(0.9 * indices)
    probs = np.maximum(base, 0.05)
    probs = probs / float(np.sum(probs))

    eta = np.exp(1j * (0.37 * indices + 0.11 * indices * indices)) / math.sqrt(dim)
    epsilon = 0.09
    sigma = (1.0 - epsilon) * np.diag(probs) + epsilon * density(eta)
    return u_record @ sigma @ dagger(u_record)


def haar_pure_state(rng: np.random.Generator, dim: int) -> np.ndarray:
    """Draw one complex Haar-distributed pure state in C^dim."""

    z = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return z / np.linalg.norm(z)


def sign_witness(delta: np.ndarray) -> np.ndarray:
    """Hermitian unit-norm witness saturating |Tr(delta O)| <= ||delta||_1."""

    vals, vecs = np.linalg.eigh(delta)
    return vecs @ np.diag(np.sign(vals)) @ dagger(vecs)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def named_record_observables() -> dict[str, np.ndarray]:
    """Detector-admissible observables generated by the chosen record algebra."""

    t_plus = t_loop_basis()[0][1]
    p_t = density(t_plus)
    p_d1 = density(ket(1))
    return {
        "Bell stabilizer ZZ_AB": kron(Z, Z, I2, I2),
        "Bell stabilizer XX_AB": kron(X, X, I2, I2),
        "T-loop projector P_T+": kron(np.eye(4, dtype=complex), p_t, I2),
        "detector branch Z_D": kron(I2, I2, I2, Z),
        "joint P_T+ and D1": kron(np.eye(4, dtype=complex), p_t, p_d1),
    }


def off_record_witnesses() -> dict[str, np.ndarray]:
    """Incompatible or off-record witnesses deliberately not preserved."""

    return {
        "detector coherence X_D": kron(I2, I2, I2, X),
        "route computational Z_R": kron(I2, I2, Z, I2),
        "local endpoint Z_A": kron(Z, I2, I2, I2),
        "local endpoint X_A": kron(X, I2, I2, I2),
    }


def main() -> None:
    print("Record-grammar compression certificate")
    print("=" * 92)

    u_record, labels = record_unitary()
    dim = u_record.shape[0]
    print(f"  Hilbert dimension: {dim}")
    print(f"  Record sectors:    {len(labels)}")
    print(f"  First sectors:     {', '.join(labels[:4])}, ...")

    print("\n[1] Record projectors are orthogonal and complete")
    assert_close("record basis unitarity", float(np.linalg.norm(dagger(u_record) @ u_record - np.eye(dim))), 0.0)
    projectors = record_projectors(u_record)
    identity_from_records = sum(projectors)
    assert_close("sum_i P_i = I", float(np.linalg.norm(identity_from_records - np.eye(dim))), 0.0)
    max_overlap = 0.0
    for i, pi in enumerate(projectors):
        assert_close(f"P_{i} idempotent", float(np.linalg.norm(pi @ pi - pi)), 0.0)
        for j, pj in enumerate(projectors):
            if i != j:
                max_overlap = max(max_overlap, float(np.linalg.norm(pi @ pj)))
    assert_close("max off-diagonal projector overlap", max_overlap, 0.0)
    print("  -> the admissible records are genuine orthogonal projectors.")

    print("\n[2] Conditional expectation is a well-behaved compression map")
    rho = test_state(u_record)
    rho_rec = conditional_expectation(rho, u_record)
    rho_rec_2 = conditional_expectation(rho_rec, u_record)
    assert_close("Tr rho", float(np.real_if_close(np.trace(rho))), 1.0)
    assert_close("Tr E(rho)", float(np.real_if_close(np.trace(rho_rec))), 1.0)
    assert_close("Hermiticity of rho", float(np.linalg.norm(rho - dagger(rho))), 0.0)
    assert_close("idempotence E(E(rho))=E(rho)", float(np.linalg.norm(rho_rec_2 - rho_rec)), 0.0)
    assert_greater("min eigenvalue rho", float(np.min(np.linalg.eigvalsh(rho))), -1e-12)
    assert_greater("min eigenvalue E(rho)", float(np.min(np.linalg.eigvalsh(rho_rec))), -1e-12)
    print("  -> the map keeps a valid state and removes only off-record coherences.")

    print("\n[3] Every observable in the record algebra is preserved")
    for name, op in named_record_observables().items():
        before = expect(rho, op)
        after = expect(rho_rec, op)
        assert_close(name, after, before)

    rng = np.random.default_rng(20260629)
    for trial in range(5):
        weights = rng.uniform(-1.0, 1.0, size=dim)
        op = record_observable(weights, u_record)
        assert_close(f"random record observable {trial}", expect(rho_rec, op), expect(rho, op))
    print("  -> detector-readable record expectations survive compression exactly.")

    print("\n[4] Off-record witnesses are changed, but bounded")
    delta = rho - rho_rec
    residual_trace_norm = trace_norm(delta)
    residual_frobenius = float(np.linalg.norm(delta))
    print(f"  ||rho - E(rho)||_1 = {residual_trace_norm:.12g}")
    print(f"  ||rho - E(rho)||_F = {residual_frobenius:.12g}")
    assert_greater("nonzero off-record residual", residual_trace_norm, 0.05)
    for name, op in off_record_witnesses().items():
        diff = abs(expect(delta, op))
        norm = operator_norm(op)
        print(f"  {name:<34s} |delta expectation|={diff:.12g}  ||O||={norm:.3g}")
        assert_less(f"{name} trace-norm bound", diff, residual_trace_norm * norm + 1e-12)
    print("  -> the dropped terms are visible to incompatible witnesses, not to the record algebra.")

    print("\n[5] The residual bound is sharp")
    witness = sign_witness(delta)
    saturation = abs(expect(delta, witness))
    assert_close("||sign(delta)||_op", operator_norm(witness), 1.0, tol=5e-12)
    assert_close("sign-witness saturation", saturation, residual_trace_norm, tol=5e-12)
    print("  -> the residual norm is a certificate: it is the exact worst-case lost expectation.")

    print("\n[6] Haar pure-state residual distribution")
    residuals = []
    for _ in range(64):
        psi = haar_pure_state(rng, dim)
        sigma = density(psi)
        sigma_rec = conditional_expectation(sigma, u_record)
        residuals.append(trace_norm(sigma - sigma_rec))
    qs = np.quantile(residuals, [0.0, 0.1, 0.5, 0.9, 1.0])
    print(
        "  ||rho - E(rho)||_1 quantiles over 64 Haar pure states: "
        f"min={qs[0]:.6f}, p10={qs[1]:.6f}, median={qs[2]:.6f}, "
        f"p90={qs[3]:.6f}, max={qs[4]:.6f}"
    )
    assert_greater("Haar residual median nonzero", float(qs[2]), 0.5)
    assert_less("Haar residual max bounded by trace distance diameter", float(qs[4]), 2.0 + 1e-12)
    print("  -> the selected test state is tame; generic pure states carry much larger off-record tails.")

    print("\n[7] Compression accounting")
    full_real_parameters = dim * dim - 1
    record_real_parameters = dim - 1
    print(f"  full density matrix parameters: {full_real_parameters}")
    print(f"  record probability parameters: {record_real_parameters}")
    print(f"  parameter compression factor:  {full_real_parameters / record_real_parameters:.1f}x")
    assert_close("full parameter count", float(full_real_parameters), 255.0)
    assert_close("record parameter count", float(record_real_parameters), 15.0)
    print("  -> for this finite network, the record algebra is much smaller than the full state.")

    print(
        """
Verdict
-------
The compression certificate demonstrates the record-grammar claim in a finite
Hilbert space:

  full state rho
    -> conditional expectation E_rec(rho)
    -> exact preservation of Bell/stabilizer records, T-loop records, and
       detector-branch records
    -> explicit trace-norm certificate for every off-record coherence removed.

So the grammar is not "drop correlations and hope."  For this network, the
omitted part is a measured residual.  Any detector observable inside the record
algebra is preserved exactly; any incompatible observable is bounded by
||rho - E_rec(rho)||_1, and the sign witness shows that this bound is sharp.

Boundary
--------
This is a finite four-qubit certificate, not a theorem for arbitrary many-body
matter.  Its purpose is to show what a serious record-compression claim must
look like: specify the admissible projectors, apply the conditional expectation,
preserve all allowed readouts, and report the residual instead of hiding it.
"""
    )


if __name__ == "__main__":
    main()
