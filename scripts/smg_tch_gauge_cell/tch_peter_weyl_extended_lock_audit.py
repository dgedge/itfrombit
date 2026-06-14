#!/usr/bin/env python3
"""Extended Peter-Weyl audit for the TCH plaquette/intertwiner lock.

The previous exact Peter-Weyl script kept only link irreps

    1, 3, 3bar

and therefore projected out the two channels that matter for the next check:

    3 x 3 -> 6 + 3bar
    3bar x 3 -> 1 + 8

This script adds the missing low SU(3) irreps

    1, 3, 3bar, 6, 6bar, 8

and audits the TCH-lock proposal at the fusion/intertwiner level.

What is exact here:
  * SU(3) generators for 3, 3bar, 6, 6bar, and 8.
  * Singlet/intertwiner multiplicities at each vertex, computed as nullspaces
    of the total Lie-algebra action.
  * The electric Casimir spectrum for the retained link irreps.
  * The square-plaquette geometry selecting the two electric-minimal boundary
    routings: links (1,3) versus links (0,2). These are the two ways the
    bipartite colour flux can run around opposite sides of the square face; the
    truncated-cube X-plaquette/Bianchi term is precisely the local move that
    flips one routing into the other.

What is not exact:
  * The full Clebsch-Gordan Wilson matrix for 6/6bar/8 is not built. The Wilson
    operator below is a fusion-graph adjacency with conservative normalization.
    It is a convergence/gap-stability audit, not a replacement for a full SU(3)
    tensor-network link Hilbert space.

Verdict logic:
  * If the extended channels close the gap even with the geometry-selected
    O(1) lock, the lock route is dead.
  * If the O(1) lock keeps a robust gap while Wilson-only remains soft, the next
    required artifact is a full-CG implementation deriving the same lock
    coefficient from the truncated-cube gauge-cell X-plaquette.
"""

from __future__ import annotations

import argparse
import itertools
import json
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "python_code" / "smg_dmrg_runs" / "tch_peter_weyl_extended_lock_verdict.json"

REPS = ("1", "3", "3b", "6", "6b", "8")
REP_LABEL = {
    "1": (0, 0),
    "3": (1, 0),
    "3b": (0, 1),
    "6": (2, 0),
    "6b": (0, 2),
    "8": (1, 1),
}
MATTER_REPS = ("3", "3b", "3", "3b")


def gell_mann() -> list[np.ndarray]:
    mats = [
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]],
        [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]],
        [[1, 0, 0], [0, -1, 0], [0, 0, 0]],
        [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
        [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]],
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
        [[1 / np.sqrt(3), 0, 0], [0, 1 / np.sqrt(3), 0], [0, 0, -2 / np.sqrt(3)]],
    ]
    return [np.array(m, dtype=complex) / 2 for m in mats]


FUND_GENS = gell_mann()
ANTI_GENS = [-g.conj() for g in FUND_GENS]


def symmetric_basis() -> np.ndarray:
    """Columns embed the 6 symmetric tensors into 3 x 3."""
    cols = []
    for i in range(3):
        vec = np.zeros(9, complex)
        vec[3 * i + i] = 1.0
        cols.append(vec)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        vec = np.zeros(9, complex)
        vec[3 * i + j] = 1 / np.sqrt(2)
        vec[3 * j + i] = 1 / np.sqrt(2)
        cols.append(vec)
    return np.column_stack(cols)


SYM_BASIS = symmetric_basis()


def symmetric_generators(base_gens: list[np.ndarray]) -> list[np.ndarray]:
    out = []
    for generator in base_gens:
        full = np.kron(generator, np.eye(3)) + np.kron(np.eye(3), generator)
        out.append(SYM_BASIS.conj().T @ full @ SYM_BASIS)
    return out


def adjoint_generators() -> list[np.ndarray]:
    # Orthonormal traceless matrix basis B_a = lambda_a / sqrt(2), with
    # Tr(B_a^dag B_b)=delta_ab. The commutator action is real antisymmetric;
    # multiply by i so the adjoint uses the same Hermitian-generator convention
    # as 3, 3bar, 6, and 6bar.
    basis = [2 * g / np.sqrt(2) for g in FUND_GENS]
    out = []
    for generator in FUND_GENS:
        mat = np.zeros((8, 8), complex)
        for row, b_row in enumerate(basis):
            for col, b_col in enumerate(basis):
                comm = generator @ b_col - b_col @ generator
                mat[row, col] = -1j * np.trace(b_row.conj().T @ comm)
        out.append(1j * mat)
    return out


GENS = {
    "1": [np.zeros((1, 1), complex) for _ in range(8)],
    "3": FUND_GENS,
    "3b": ANTI_GENS,
    "6": symmetric_generators(FUND_GENS),
    "6b": symmetric_generators(ANTI_GENS),
    "8": adjoint_generators(),
}


def dim(rep: str) -> int:
    return GENS[rep][0].shape[0]


def dual(rep: str) -> str:
    return {"1": "1", "3": "3b", "3b": "3", "6": "6b", "6b": "6", "8": "8"}[rep]


def triality(rep: str) -> int:
    p, q = REP_LABEL[rep]
    return (p - q) % 3


def casimir(rep: str) -> float:
    p, q = REP_LABEL[rep]
    return (p * p + q * q + p * q + 3 * p + 3 * q) / 3


@lru_cache(maxsize=None)
def invariant_multiplicity(reps: tuple[str, str, str]) -> int:
    dims = [dim(rep) for rep in reps]
    total_dim = int(np.prod(dims))
    blocks = []
    for gen_index in range(8):
        total = np.zeros((total_dim, total_dim), complex)
        for pos, rep in enumerate(reps):
            factors = []
            for idx, rep_here in enumerate(reps):
                factors.append(GENS[rep][gen_index] if idx == pos else np.eye(dim(rep_here)))
            term = factors[0]
            for factor in factors[1:]:
                term = np.kron(term, factor)
            total += term
        blocks.append(total)
    stacked = np.vstack(blocks)
    singular = np.linalg.svd(stacked, compute_uv=False)
    return int(np.sum(singular < 1e-9))


def vertex_triples(link_reps: tuple[str, str, str, str]) -> list[tuple[str, str, str]]:
    r0, r1, r2, r3 = link_reps
    return [
        (MATTER_REPS[0], r0, r3),
        (MATTER_REPS[1], dual(r0), r1),
        (MATTER_REPS[2], dual(r1), dual(r2)),
        (MATTER_REPS[3], dual(r3), r2),
    ]


def allowed_assignments(reps: tuple[str, ...]) -> list[dict[str, Any]]:
    out = []
    for link_reps in itertools.product(reps, repeat=4):
        mults = [invariant_multiplicity(triple) for triple in vertex_triples(link_reps)]
        if all(mult > 0 for mult in mults):
            out.append(
                {
                    "links": link_reps,
                    "vertex_multiplicities": mults,
                    "multiplicity_product": int(np.prod(mults)),
                    "electric_energy": float(sum(casimir(rep) for rep in link_reps)),
                    "center_flux": int((triality(link_reps[0]) + triality(link_reps[1])
                                        - triality(link_reps[2]) - triality(link_reps[3])) % 3),
                }
            )
    return out


def fuse_by_fundamental(rep: str, dagger: bool = False) -> list[str]:
    # Exact SU(3) product with 3 or 3bar, restricted to REPS.
    if not dagger:
        raw = {
            "1": ["3"],
            "3": ["6", "3b"],
            "3b": ["1", "8"],
            "6": ["8"],
            "6b": ["3b"],
            "8": ["3", "6b"],
        }
    else:
        raw = {
            "1": ["3b"],
            "3": ["1", "8"],
            "3b": ["6b", "3"],
            "6": ["3"],
            "6b": ["8"],
            "8": ["3b", "6"],
        }
    return [target for target in raw[rep] if target in REPS]


def wilson_matrix(assignments: list[dict[str, Any]]) -> np.ndarray:
    index = {row["links"]: pos for pos, row in enumerate(assignments)}
    dim_basis = len(assignments)
    w = np.zeros((dim_basis, dim_basis), complex)
    for col, row in enumerate(assignments):
        links = row["links"]
        target_lists = [
            fuse_by_fundamental(links[0], dagger=False),
            fuse_by_fundamental(links[1], dagger=False),
            fuse_by_fundamental(links[2], dagger=True),
            fuse_by_fundamental(links[3], dagger=True),
        ]
        norm = np.sqrt(np.prod([max(len(targets), 1) for targets in target_lists]))
        for target_links in itertools.product(*target_lists):
            row_index = index.get(target_links)
            if row_index is None:
                continue
            # Conservative fusion-graph normalization; full CG coefficients are
            # intentionally left to the next artifact.
            w[row_index, col] += 1.0 / norm
    return w


def electric_matrix(assignments: list[dict[str, Any]]) -> np.ndarray:
    return np.diag([row["electric_energy"] for row in assignments]).astype(complex)


def minimal_routings(assignments: list[dict[str, Any]]) -> list[int]:
    min_energy = min(row["electric_energy"] for row in assignments)
    return [
        idx for idx, row in enumerate(assignments)
        if abs(row["electric_energy"] - min_energy) < 1e-12
    ]


def geometric_lock_matrix(assignments: list[dict[str, Any]], strength: float) -> np.ndarray:
    """Square-plaquette X lock between the two electric-minimal boundary routings."""
    lock = np.zeros((len(assignments), len(assignments)), complex)
    # These are the two opposite-side routings around the oriented square face:
    # path e1+e3 and path e0+e2. They are exactly the two low spin networks in
    # the previous explicit 1,3,3bar script.
    routing_a = ("1", "3", "1", "3b")
    routing_b = ("3b", "1", "3", "1")
    index = {row["links"]: pos for pos, row in enumerate(assignments)}
    if routing_a in index and routing_b in index:
        i = index[routing_a]
        j = index[routing_b]
        lock[i, j] = -strength
        lock[j, i] = -strength
    return lock


def unique_gap(eigvals: np.ndarray, tol: float = 1e-9) -> tuple[float, int]:
    levels = []
    degeneracies = []
    for value in eigvals:
        if not levels or abs(value - levels[-1]) > tol:
            levels.append(float(value))
            degeneracies.append(1)
        else:
            degeneracies[-1] += 1
    return levels[1] - levels[0], degeneracies[0]


def scan(assignments: list[dict[str, Any]], lock_strength: float) -> list[dict[str, Any]]:
    e = electric_matrix(assignments)
    w = wilson_matrix(assignments)
    lock = geometric_lock_matrix(assignments, lock_strength)
    rows = []
    for beta in [0.25, 0.5, 0.75, 1.0]:
        h = (1 / beta) * e - (beta / 2) * (w + w.conj().T) + lock
        herm = float(np.linalg.norm(h - h.conj().T))
        assert herm < 1e-10
        eigvals = np.linalg.eigvalsh(h)
        gap, degeneracy = unique_gap(eigvals)
        rows.append(
            {
                "beta": beta,
                "gap": gap,
                "ground_degeneracy": degeneracy,
                "ground_energy": float(eigvals[0]),
                "lowest": [float(x) for x in eigvals[:8]],
            }
        )
    return rows


def build_verdict() -> dict[str, Any]:
    minimal_reps = ("1", "3", "3b")
    extended_reps = REPS
    minimal_assignments = allowed_assignments(minimal_reps)
    extended_assignments = allowed_assignments(extended_reps)

    max_vertex_mult = max(
        max(row["vertex_multiplicities"]) for row in extended_assignments
    )
    min_routing_rows = [extended_assignments[i] for i in minimal_routings(extended_assignments)]

    minimal_wilson = scan(minimal_assignments, lock_strength=0.0)
    minimal_lock = scan(minimal_assignments, lock_strength=1.0)
    extended_wilson = scan(extended_assignments, lock_strength=0.0)
    extended_lock = scan(extended_assignments, lock_strength=1.0)

    minimal_wilson_min = min(row["gap"] for row in minimal_wilson)
    minimal_lock_min = min(row["gap"] for row in minimal_lock)
    extended_wilson_min = min(row["gap"] for row in extended_wilson)
    extended_lock_min = min(row["gap"] for row in extended_lock)

    assert len(minimal_assignments) == 3
    assert len(extended_assignments) > len(minimal_assignments)
    assert minimal_wilson_min < 0.2
    assert minimal_lock_min > 1.0
    assert extended_lock_min > 1.0

    return {
        "verdict": "extended_channels_do_not_kill_o1_geometry_lock",
        "scope": (
            "fusion/intertwiner-level audit; Wilson matrix uses normalized fusion adjacency, "
            "not full 6/8 Clebsch-Gordan coefficients"
        ),
        "basis": {
            "minimal_reps": minimal_reps,
            "extended_reps": extended_reps,
            "minimal_allowed_assignments": len(minimal_assignments),
            "extended_allowed_assignments": len(extended_assignments),
            "max_vertex_intertwiner_multiplicity": max_vertex_mult,
            "electric_minimal_routings": min_routing_rows,
        },
        "geometry_lock": {
            "strength": 1.0,
            "meaning": (
                "unit-normalized truncated-cube square-face X-plaquette move coupling the "
                "two electric-minimal opposite-side boundary routings ('1','3','1','3b') "
                "and ('3b','1','3','1')"
            ),
        },
        "scans": {
            "minimal_1_3_3b_wilson_only": minimal_wilson,
            "minimal_1_3_3b_with_geometry_lock": minimal_lock,
            "extended_1_3_3b_6_6b_8_wilson_only": extended_wilson,
            "extended_1_3_3b_6_6b_8_with_geometry_lock": extended_lock,
        },
        "gap_summary": {
            "minimal_wilson_only_min_gap": minimal_wilson_min,
            "minimal_with_lock_min_gap": minimal_lock_min,
            "extended_wilson_only_min_gap": extended_wilson_min,
            "extended_with_lock_min_gap": extended_lock_min,
        },
        "next_required_artifact": (
            "replace fusion-adjacency Wilson entries by full Clebsch-Gordan tensors for 6, "
            "6bar, and 8; then couple this extended plaquette basis to the charged mirror-Fock "
            "SMG block"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    verdict = build_verdict()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")

    print("Extended Peter-Weyl TCH-lock audit")
    print(f"  verdict: {verdict['verdict']}")
    print(f"  wrote: {args.out.relative_to(ROOT)}")
    print(f"  minimal allowed assignments: {verdict['basis']['minimal_allowed_assignments']}")
    print(f"  extended allowed assignments: {verdict['basis']['extended_allowed_assignments']}")
    print(f"  max vertex intertwiner multiplicity: {verdict['basis']['max_vertex_intertwiner_multiplicity']}")
    print("  min gaps:")
    for key, value in verdict["gap_summary"].items():
        print(f"    {key}: {value:.6g}")
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
