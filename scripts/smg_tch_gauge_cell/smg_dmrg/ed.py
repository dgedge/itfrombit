"""Small exact-diagonalisation guardrails for the SMG prototype."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import ModelSpec, mode_index
from .terms import FermionOp, FermionTerm


@dataclass(frozen=True)
class ExactResult:
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    gap: float
    hermiticity_residual: float
    mirror_bilinear_m2: float


def parity_before(state: int, site: int) -> int:
    mask = (1 << site) - 1
    return (state & mask).bit_count() & 1


def apply_single_op(state: int, op: FermionOp) -> tuple[int, complex] | None:
    occupied = (state >> op.site) & 1
    if op.kind == "number":
        return (state, 1.0) if occupied else None

    sign = -1.0 if parity_before(state, op.site) else 1.0
    if op.kind == "annihilate":
        if not occupied:
            return None
        return state & ~(1 << op.site), sign
    if op.kind == "create":
        if occupied:
            return None
        return state | (1 << op.site), sign
    raise ValueError(f"unknown op kind: {op.kind}")


def apply_ops(state: int, ops: tuple[FermionOp, ...]) -> tuple[int, complex] | None:
    coeff = 1.0 + 0j
    current = state
    for op in reversed(ops):
        applied = apply_single_op(current, op)
        if applied is None:
            return None
        current, factor = applied
        coeff *= factor
    return current, coeff


def build_dense_hamiltonian(total_modes: int, terms: list[FermionTerm]) -> np.ndarray:
    if total_modes > 14:
        raise ValueError("dense ED guardrail refuses >14 modes")
    dim = 1 << total_modes
    hamiltonian = np.zeros((dim, dim), dtype=complex)
    for col in range(dim):
        for term in terms:
            applied = apply_ops(col, term.ops)
            if applied is None:
                continue
            row, factor = applied
            hamiltonian[row, col] += term.coeff * factor
    return hamiltonian


def many_body_gap(eigenvalues: np.ndarray, tol: float = 1e-9) -> float:
    levels: list[float] = []
    for value in np.sort(np.real(eigenvalues)):
        if not levels or value - levels[-1] > tol:
            levels.append(float(value))
        if len(levels) == 2:
            return levels[1] - levels[0]
    return 0.0


def build_pair_annihilator(total_modes: int, pairs: list[tuple[int, int]]) -> np.ndarray:
    dim = 1 << total_modes
    out = np.zeros((dim, dim), dtype=complex)
    pair_terms = [
        FermionTerm(1.0, (FermionOp("annihilate", i), FermionOp("annihilate", j)))
        for i, j in pairs
    ]
    for col in range(dim):
        for term in pair_terms:
            applied = apply_ops(col, term.ops)
            if applied is None:
                continue
            row, factor = applied
            out[row, col] += factor
    return out


def mirror_bilinear_m2(spec: ModelSpec, ground_state: np.ndarray) -> float:
    total = 0.0 + 0j
    ops = []
    for cell in range(spec.cells):
        pairs = []
        for flavor_a in range(spec.flavors):
            for flavor_b in range(flavor_a + 1, spec.flavors):
                pairs.append(
                    (
                        mode_index(spec, cell, "mirror", flavor_a),
                        mode_index(spec, cell, "mirror", flavor_b),
                    )
                )
        ops.append(build_pair_annihilator(spec.total_modes, pairs))

    for ox in ops:
        for oy in ops:
            total += np.vdot(ground_state, ox.conj().T @ (oy @ ground_state))
    return float(np.real(total) / (spec.cells * spec.cells))


def solve_exact(spec: ModelSpec, terms: list[FermionTerm]) -> ExactResult:
    hamiltonian = build_dense_hamiltonian(spec.total_modes, terms)
    hermiticity = float(np.linalg.norm(hamiltonian - hamiltonian.conj().T))
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    ground_state = eigenvectors[:, 0]
    return ExactResult(
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        gap=many_body_gap(eigenvalues),
        hermiticity_residual=hermiticity,
        mirror_bilinear_m2=mirror_bilinear_m2(spec, ground_state),
    )

