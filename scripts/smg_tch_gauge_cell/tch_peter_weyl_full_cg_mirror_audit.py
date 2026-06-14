#!/usr/bin/env python3
"""Full-Clebsch-Gordan TCH plaquette audit plus mirror-block stress test.

This is the next artifact after ``tch_peter_weyl_extended_lock_audit.py``.
That script included the low SU(3) channels

    1, 3, 3bar, 6, 6bar, 8

but used a normalized fusion-graph adjacency as the Wilson plaquette proxy.
Here the Wilson entries are replaced by actual Peter-Weyl/Clebsch-Gordan
matrix elements for the retained channels.

What is exact here:
  * CG embeddings are solved as SU(3) Lie-algebra intertwiners:
        (3 or 3bar) tensor R_source -> R_target.
  * Link multiplication uses the normalized Peter-Weyl coefficient
        sqrt(dim(source)/dim(target)) C_tail C_head^*.
  * The plaquette Wilson matrix is obtained by contracting the four vertex
    singlet intertwiners, not by a fusion adjacency.
  * The old 1,3,3bar exact script is reproduced spectrally in the overlap
    truncation, which is the regression guard.

What remains bounded:
  * The mirror-block coupling below is a finite hybrid stress test. It couples
    the full-CG plaquette oscillator to the already-built charged mirror-Fock
    SMG open-chain hopping channel. This tests whether a modest noncommuting
    gauge-plaquette modulation immediately destroys the mirror gap. It is not a
    final 3+1D SU(3) charged-hopping construction.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PYTHON_CODE = ROOT / "python_code"
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_peter_weyl_full_cg_mirror_verdict.json"

MATTER_REPS = ("3", "3b", "3", "3b")
REPS = ("1", "3", "3b", "6", "6b", "8")
TOL = 1e-10


def load_module(name: str, relative_path: str):
    path = PYTHON_CODE / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


EXT = load_module("tch_peter_weyl_extended_lock_mod", "tch_peter_weyl_extended_lock_audit.py")
GENS = dict(EXT.GENS)
# Older versions of the extended audit kept the adjoint action as real
# antisymmetric matrices. Multiplying by i puts the 8 in the same
# Hermitian-generator convention as 3, 3bar, 6, and 6bar. Keep this conditional
# so the full-CG audit remains stable after the source artifact is corrected.
if max(np.linalg.norm(generator - generator.conj().T) for generator in GENS["8"]) > 1e-10:
    GENS["8"] = [1j * generator for generator in GENS["8"]]


def dim(rep: str) -> int:
    return GENS[rep][0].shape[0]


def dual(rep: str) -> str:
    return EXT.dual(rep)


def casimir(rep: str) -> float:
    return EXT.casimir(rep)


def fuse_targets(factor_rep: str, source_rep: str) -> list[str]:
    if factor_rep == "3":
        return EXT.fuse_by_fundamental(source_rep, dagger=False)
    if factor_rep == "3b":
        return EXT.fuse_by_fundamental(source_rep, dagger=True)
    raise ValueError(f"unsupported factor rep {factor_rep!r}")


def generator_sum(reps: tuple[str, ...]) -> list[np.ndarray]:
    dims = [dim(rep) for rep in reps]
    total_dim = int(np.prod(dims))
    out = []
    for gen_index in range(8):
        total = np.zeros((total_dim, total_dim), complex)
        for pos, rep in enumerate(reps):
            factors = [
                GENS[rep][gen_index] if idx == pos else np.eye(dim(rep_here), dtype=complex)
                for idx, rep_here in enumerate(reps)
            ]
            term = factors[0]
            for factor in factors[1:]:
                term = np.kron(term, factor)
            total += term
        out.append(total)
    return out


def null_vector(stacked: np.ndarray, expected_nullity: int = 1) -> np.ndarray:
    _, singular, vh = np.linalg.svd(stacked)
    nullity = int(np.sum(singular < 1e-9))
    assert nullity == expected_nullity, (nullity, singular[-8:])
    return vh.conj().T[:, -1]


@lru_cache(maxsize=None)
def invariant_tensor(reps: tuple[str, str, str]) -> np.ndarray:
    blocks = generator_sum(reps)
    stacked = np.vstack(blocks)
    vec = null_vector(stacked, expected_nullity=1)
    vec = vec / np.linalg.norm(vec)
    pivot = np.argmax(np.abs(vec))
    if abs(vec[pivot]) > 0:
        vec *= np.conj(vec[pivot]) / abs(vec[pivot])
    return vec.reshape(tuple(dim(rep) for rep in reps))


@lru_cache(maxsize=None)
def invariant_multiplicity(reps: tuple[str, str, str]) -> int:
    stacked = np.vstack(generator_sum(reps))
    singular = np.linalg.svd(stacked, compute_uv=False)
    return int(np.sum(singular < 1e-9))


@lru_cache(maxsize=None)
def cg_embedding(factor_rep: str, source_rep: str, target_rep: str) -> np.ndarray:
    """Return Q[(factor,source), target] embedding target into factor tensor source."""
    assert target_rep in fuse_targets(factor_rep, source_rep)
    factor_dim = dim(factor_rep)
    source_dim = dim(source_rep)
    target_dim = dim(target_rep)
    product_dim = factor_dim * source_dim

    product_gens = generator_sum((factor_rep, source_rep))
    target_gens = GENS[target_rep]
    equations = []
    for product_gen, target_gen in zip(product_gens, target_gens):
        equations.append(
            np.kron(np.eye(target_dim, dtype=complex), product_gen)
            - np.kron(target_gen.T, np.eye(product_dim, dtype=complex))
        )
    stacked = np.vstack(equations)
    vec = null_vector(stacked, expected_nullity=1)
    q = vec.reshape((product_dim, target_dim), order="F")

    gram = q.conj().T @ q
    evals, evecs = np.linalg.eigh(gram)
    assert float(np.min(evals)) > 1e-12
    q = q @ (evecs @ np.diag(1 / np.sqrt(evals)) @ evecs.conj().T)

    pivot = np.argmax(np.abs(q))
    if abs(q.flat[pivot]) > 0:
        q *= np.conj(q.flat[pivot]) / abs(q.flat[pivot])

    iso_error = float(np.linalg.norm(q.conj().T @ q - np.eye(target_dim)))
    intertwiner_error = max(
        float(np.linalg.norm(product_gen @ q - q @ target_gen))
        for product_gen, target_gen in zip(product_gens, target_gens)
    )
    assert iso_error < 1e-10
    assert intertwiner_error < 1e-10
    return q


def local_link_basis(rep: str) -> list[tuple[str, int, int]]:
    return [(rep, tail, head) for tail in range(dim(rep)) for head in range(dim(rep))]


@lru_cache(maxsize=None)
def u_components(i: int, j: int, source: tuple[str, int, int]) -> tuple[tuple[tuple[str, int, int], complex], ...]:
    source_rep, source_tail, source_head = source
    out = []
    for target_rep in fuse_targets("3", source_rep):
        q = cg_embedding("3", source_rep, target_rep)
        source_dim = dim(source_rep)
        target_dim = dim(target_rep)
        norm = np.sqrt(source_dim / target_dim)
        left_row = i * source_dim + source_tail
        right_row = j * source_dim + source_head
        for target_tail in range(target_dim):
            left_coeff = q[left_row, target_tail]
            if abs(left_coeff) < 1e-14:
                continue
            for target_head in range(target_dim):
                coeff = norm * left_coeff * np.conj(q[right_row, target_head])
                if abs(coeff) > 1e-14:
                    out.append(((target_rep, target_tail, target_head), coeff))
    return tuple(out)


@lru_cache(maxsize=None)
def udag_components(i: int, j: int, source: tuple[str, int, int]) -> tuple[tuple[tuple[str, int, int], complex], ...]:
    """Projected adjoint of U_{ji}, matching the earlier exact script convention."""
    out = []
    for candidate_rep in REPS:
        for candidate in local_link_basis(candidate_rep):
            for target, coeff in u_components(j, i, candidate):
                if target == source:
                    out.append((candidate, np.conj(coeff)))
    return tuple(out)


@lru_cache(maxsize=None)
def link_map(kind: str, i: int, j: int, target_rep: str, source_rep: str) -> np.ndarray:
    target_dim = dim(target_rep)
    source_dim = dim(source_rep)
    matrix = np.zeros((target_dim, target_dim, source_dim, source_dim), complex)
    for source_tail in range(source_dim):
        for source_head in range(source_dim):
            source = (source_rep, source_tail, source_head)
            components = u_components(i, j, source) if kind == "u" else udag_components(i, j, source)
            for target, coeff in components:
                rep, target_tail, target_head = target
                if rep == target_rep:
                    matrix[target_tail, target_head, source_tail, source_head] += coeff
    return matrix


@lru_cache(maxsize=None)
def identity_link_map(target_rep: str, source_rep: str) -> np.ndarray:
    target_dim = dim(target_rep)
    source_dim = dim(source_rep)
    matrix = np.zeros((target_dim, target_dim, source_dim, source_dim), complex)
    if target_rep == source_rep:
        for tail in range(target_dim):
            for head in range(target_dim):
                matrix[tail, head, tail, head] = 1.0
    return matrix


@lru_cache(maxsize=None)
def vertex_tensors(link_reps: tuple[str, str, str, str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    r0, r1, r2, r3 = link_reps
    return (
        invariant_tensor((MATTER_REPS[0], r0, r3)),
        invariant_tensor((MATTER_REPS[1], dual(r0), r1)),
        invariant_tensor((MATTER_REPS[2], dual(r1), dual(r2))),
        invariant_tensor((MATTER_REPS[3], dual(r3), r2)),
    )


def contract_network(
    left_reps: tuple[str, str, str, str],
    right_reps: tuple[str, str, str, str],
    maps: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> complex:
    l0, l1, l2, l3 = vertex_tensors(left_reps)
    r0, r1, r2, r3 = vertex_tensors(right_reps)
    m0, m1, m2, m3 = maps
    return np.einsum(
        "xAD,yBC,zEF,wHG,ABab,CEce,GFgf,DHdh,xad,ybc,zef,whg->",
        l0.conj(),
        l1.conj(),
        l2.conj(),
        l3.conj(),
        m0,
        m1,
        m2,
        m3,
        r0,
        r1,
        r2,
        r3,
        optimize=True,
    )


def vertex_triples(link_reps: tuple[str, str, str, str]) -> tuple[tuple[str, str, str], ...]:
    r0, r1, r2, r3 = link_reps
    return (
        (MATTER_REPS[0], r0, r3),
        (MATTER_REPS[1], dual(r0), r1),
        (MATTER_REPS[2], dual(r1), dual(r2)),
        (MATTER_REPS[3], dual(r3), r2),
    )


def allowed_assignments(reps: tuple[str, ...]) -> list[tuple[str, str, str, str]]:
    out = []
    for link_reps in itertools.product(reps, repeat=4):
        if all(invariant_multiplicity(triple) == 1 for triple in vertex_triples(link_reps)):
            out.append(tuple(link_reps))
    return out


def gram_matrix(assignments: list[tuple[str, str, str, str]]) -> np.ndarray:
    gram = np.zeros((len(assignments), len(assignments)), complex)
    for row, left in enumerate(assignments):
        for col, right in enumerate(assignments):
            maps = tuple(identity_link_map(left[idx], right[idx]) for idx in range(4))
            gram[row, col] = contract_network(left, right, maps)
    herm = float(np.linalg.norm(gram - gram.conj().T))
    assert herm < 1e-10
    return gram


def full_cg_wilson_matrix(assignments: list[tuple[str, str, str, str]]) -> np.ndarray:
    wilson = np.zeros((len(assignments), len(assignments)), complex)
    for row, left in enumerate(assignments):
        for col, right in enumerate(assignments):
            total = 0.0j
            for a, b, c, d in itertools.product(range(3), repeat=4):
                maps = (
                    link_map("u", a, b, left[0], right[0]),
                    link_map("u", b, c, left[1], right[1]),
                    link_map("udag", c, d, left[2], right[2]),
                    link_map("udag", d, a, left[3], right[3]),
                )
                if any(not np.any(np.abs(entry) > 1e-14) for entry in maps):
                    continue
                total += contract_network(left, right, maps) / 3
            wilson[row, col] = total
    return wilson


def orthonormalize(matrix: np.ndarray, gram: np.ndarray) -> np.ndarray:
    evals, evecs = np.linalg.eigh(gram)
    assert float(np.min(evals)) > 1e-10
    g_inv_sqrt = evecs @ np.diag(1 / np.sqrt(evals)) @ evecs.conj().T
    out = g_inv_sqrt @ matrix @ g_inv_sqrt
    herm = float(np.linalg.norm(out - out.conj().T))
    assert herm < 1e-9
    return out


def geometry_lock_matrix(assignments: list[tuple[str, str, str, str]], strength: float) -> np.ndarray:
    lock = np.zeros((len(assignments), len(assignments)), complex)
    route_a = ("1", "3", "1", "3b")
    route_b = ("3b", "1", "3", "1")
    if route_a in assignments and route_b in assignments:
        i = assignments.index(route_a)
        j = assignments.index(route_b)
        lock[i, j] = -strength
        lock[j, i] = -strength
    return lock


def electric_matrix(assignments: list[tuple[str, str, str, str]], beta: float, gram: np.ndarray) -> np.ndarray:
    electric = np.zeros_like(gram)
    for idx, reps in enumerate(assignments):
        electric[idx, idx] = (1 / beta) * sum(casimir(rep) for rep in reps) * gram[idx, idx]
    return electric


def hamiltonian(
    assignments: list[tuple[str, str, str, str]],
    gram: np.ndarray,
    wilson: np.ndarray,
    beta: float,
    lock_strength: float,
) -> np.ndarray:
    raw = electric_matrix(assignments, beta, gram) - (beta / 2) * (wilson + wilson.conj().T)
    h = orthonormalize(raw, gram)
    if lock_strength:
        h += geometry_lock_matrix(assignments, lock_strength)
    herm = float(np.linalg.norm(h - h.conj().T))
    assert herm < 1e-9
    return h


def gap(eigvals: np.ndarray) -> tuple[float, int]:
    levels = []
    degeneracies = []
    for value in eigvals:
        if not levels or abs(value - levels[-1]) > 1e-8:
            levels.append(float(value))
            degeneracies.append(1)
        else:
            degeneracies[-1] += 1
    return levels[1] - levels[0], degeneracies[0]


def regression_against_minimal_exact() -> dict[str, Any]:
    old = load_module("css_peter_weyl_su3_truncation_mod", "css_peter_weyl_su3_truncation.py")
    assignments = allowed_assignments(("1", "3", "3b"))
    gram = gram_matrix(assignments)
    wilson = full_cg_wilson_matrix(assignments)
    rows = []
    max_spectral_error = 0.0
    for beta in [0.25, 0.5, 1.0]:
        new_h = hamiltonian(assignments, gram, wilson, beta=beta, lock_strength=0.0)
        old_h, _, _, _ = old.hamiltonian(beta=beta, tch_lock=0.0)
        new_eigs = np.linalg.eigvalsh(new_h)
        old_eigs = np.linalg.eigvalsh(old_h)
        spectral_error = float(np.linalg.norm(new_eigs - old_eigs))
        max_spectral_error = max(max_spectral_error, spectral_error)
        rows.append(
            {
                "beta": beta,
                "new_eigvals": [float(x) for x in new_eigs],
                "old_eigvals": [float(x) for x in old_eigs],
                "spectral_error": spectral_error,
            }
        )
    assert max_spectral_error < 1e-8
    return {
        "status": "full_cg_reproduces_exact_minimal_truncation",
        "max_spectral_error": max_spectral_error,
        "rows": rows,
    }


def full_cg_gauge_audit() -> dict[str, Any]:
    assignments = allowed_assignments(REPS)
    gram = gram_matrix(assignments)
    wilson = full_cg_wilson_matrix(assignments)
    gram_error = float(np.linalg.norm(gram - np.eye(len(assignments))))
    wilson_norm = float(np.linalg.norm(wilson))
    wilson_antiherm = float(np.linalg.norm(wilson - wilson.conj().T))

    scans: dict[str, list[dict[str, Any]]] = {}
    summary: dict[str, float] = {}
    for label, lock_strength in [("wilson_only", 0.0), ("with_geometry_lock", 1.0)]:
        rows = []
        min_gap = float("inf")
        for beta in [0.25, 0.5, 0.75, 1.0]:
            h = hamiltonian(assignments, gram, wilson, beta=beta, lock_strength=lock_strength)
            eigvals = np.linalg.eigvalsh(h)
            row_gap, degeneracy = gap(eigvals)
            min_gap = min(min_gap, row_gap)
            rows.append(
                {
                    "beta": beta,
                    "gap": row_gap,
                    "ground_degeneracy": degeneracy,
                    "lowest": [float(x) for x in eigvals[:8]],
                }
            )
        scans[label] = rows
        summary[f"{label}_min_gap"] = min_gap

    assert gram_error < 1e-9
    assert wilson_norm > 1e-8
    assert summary["with_geometry_lock_min_gap"] > 1.0

    return {
        "status": "full_cg_extended_plaquette_built",
        "assignments": [list(row) for row in assignments],
        "dimension": len(assignments),
        "gram_error_vs_identity": gram_error,
        "wilson_frobenius_norm": wilson_norm,
        "wilson_nonhermitian_norm_expected": wilson_antiherm,
        "scans": scans,
        "gap_summary": summary,
    }


def normalized_hermitian_part(matrix: np.ndarray) -> np.ndarray:
    herm = (matrix + matrix.conj().T) / 2
    norm = float(np.linalg.norm(herm))
    if norm < 1e-12:
        return np.zeros_like(herm)
    return herm / norm


def coupled_mirror_stress_audit(gauge: dict[str, Any]) -> dict[str, Any]:
    chain = load_module("css_spatial_mirror_chain_mod", "css_spatial_mirror_fock_gauge_gap.py")

    assignments = [tuple(row) for row in gauge["assignments"]]
    gram = gram_matrix(assignments)
    wilson = full_cg_wilson_matrix(assignments)
    h_gauge = hamiltonian(assignments, gram, wilson, beta=0.5, lock_strength=1.0)
    x_gauge = normalized_hermitian_part(wilson)

    h_mirror, charge_counts, components = chain.build_chain_hamiltonian(
        length=2,
        beta=0.5,
        hopping=0.2,
        states_per_charge=4,
    )
    h_hopping = h_mirror - components["matter"] - components["electric"]
    x_mirror = normalized_hermitian_part(h_hopping)

    identity_g = np.eye(h_gauge.shape[0], dtype=complex)
    identity_m = np.eye(h_mirror.shape[0], dtype=complex)
    rows = []
    min_full_gap = float("inf")
    min_matter_gap = float("inf")
    for gamma in [0.0, 0.1, 0.2, 0.5]:
        total = (
            np.kron(h_gauge, identity_m)
            + np.kron(identity_g, h_mirror)
            + gamma * np.kron(x_gauge, x_mirror)
        )
        herm = float(np.linalg.norm(total - total.conj().T))
        assert herm < 1e-9
        eigvals, eigvecs = np.linalg.eigh(total)
        full_gap, degeneracy = gap(eigvals)
        matter_total = np.kron(identity_g, components["matter"])
        ground_matter = float(np.real(eigvecs[:, 0].conj() @ matter_total @ eigvecs[:, 0]))
        matter_gap = float("inf")
        matter_delta = 0.0
        for idx in range(1, len(eigvals)):
            matter_expect = float(np.real(eigvecs[:, idx].conj() @ matter_total @ eigvecs[:, idx]))
            if matter_expect - ground_matter > 0.5:
                matter_gap = float(eigvals[idx] - eigvals[0])
                matter_delta = matter_expect - ground_matter
                break
        min_full_gap = min(min_full_gap, full_gap)
        min_matter_gap = min(min_matter_gap, matter_gap)
        rows.append(
            {
                "gamma": gamma,
                "dimension": int(total.shape[0]),
                "full_gap": full_gap,
                "matter_dominated_gap": matter_gap,
                "matter_expectation_delta": matter_delta,
                "ground_degeneracy": degeneracy,
            }
        )

    assert min_full_gap > 0.1
    assert min_matter_gap > 1.0
    return {
        "status": "hybrid_full_cg_plaquette_mirror_block_gap_survives",
        "scope": (
            "finite hybrid stress test: full-CG plaquette oscillator tensor-coupled "
            "to the charged mirror-Fock SMG chain hopping channel; not a final "
            "3+1D SU(3) charged-hopping construction"
        ),
        "mirror_case": {
            "length": 2,
            "beta": 0.5,
            "hopping": 0.2,
            "states_per_charge": 4,
            "charge_counts": {str(k): int(v) for k, v in charge_counts.items()},
        },
        "minimum_full_gap": min_full_gap,
        "minimum_matter_gap": min_matter_gap,
        "rows": rows,
    }


def build_verdict() -> dict[str, Any]:
    regression = regression_against_minimal_exact()
    gauge = full_cg_gauge_audit()
    mirror = coupled_mirror_stress_audit(gauge)
    return {
        "verdict": "full_cg_gate_closed_hybrid_mirror_stress_passes",
        "scope": (
            "The normalized fusion-adjacency Wilson proxy is replaced by full numerical "
            "Clebsch-Gordan/Peter-Weyl matrix elements for 6, 6bar, and 8. The mirror "
            "coupling is still a finite hybrid stress test, so continuum/chiral closure "
            "is not claimed."
        ),
        "checks": {
            "minimal_regression": regression,
            "full_cg_gauge": gauge,
            "mirror_stress": mirror,
        },
        "remaining_gate": (
            "promote the hybrid hopping modulation to a true SU(3) charged mirror-Fock "
            "hopping operator with CG link endpoints on a multi-plaquette TCH cell, then "
            "perform finite-size/cutoff scaling"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    gauge = verdict["checks"]["full_cg_gauge"]
    mirror = verdict["checks"]["mirror_stress"]
    print("TCH full-CG Peter-Weyl + mirror-block audit")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    print(f"  minimal regression max spectral error: {verdict['checks']['minimal_regression']['max_spectral_error']:.3e}")
    print(f"  full-CG plaquette dimension: {gauge['dimension']}")
    print("  full-CG min gaps:")
    for key, value in gauge["gap_summary"].items():
        print(f"    {key}: {value:.6g}")
    print(f"  hybrid mirror stress min full gap: {mirror['minimum_full_gap']:.6g}")
    print(f"  hybrid mirror stress min matter gap: {mirror['minimum_matter_gap']:.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
