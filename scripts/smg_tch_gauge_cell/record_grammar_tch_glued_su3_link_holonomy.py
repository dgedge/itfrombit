#!/usr/bin/env python3
r"""Non-commuting SU(3) link holonomies on a glued TCH shell.

Purpose
-------
This is the next rung after ``record_grammar_tch_glued_surface_selector.py``.
That script proved a scalar shell surface-action statement: on a glued chain of
truncated-cube cells, a local strong-coupling surface sum selects the minimal
shell surface and needs only one collective inside/outside promotion.

This script puts actual non-commuting SU(3) link matrices on the same glued
shell.  It checks the object a detector can read:

    W(C) = (1/3) Re Tr P prod_{e in C} U_e ,

where the product is path ordered around the closed boundary cycle C.

What is proved here
-------------------
1. Random near-identity SU(3) link records on the glued shell are genuinely
   non-commuting.
2. Boundary Wilson traces are invariant under local vertex gauge rotations.
3. For a one-face patch, the face holonomy and boundary Wilson loop agree.
4. For multi-face patches, the tempting "multiply independent face scalars"
   shortcut fails.  Non-abelian face records require based/transported
   holonomies or a path-integral surface action; they are not independent
   commuting weights.

Boundary
--------
This is a link-holonomy/detector certificate, not a confinement proof.  It
keeps the previous minimal-shell geometry as the correct surface ledger, but it
shows that the non-abelian dynamics cannot be reduced to independent face
scalars.  The next rung is a true SU(3) Wilson action/path-integral average on
this shell or a small TCH bulk patch.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass
from itertools import combinations

import numpy as np

from record_grammar_tch_glued_surface_selector import (
    SURFACE_Q,
    GluedComplex,
    build_glued_complex,
    connected_face_patches,
    surface_row,
)


N_COLOR = 3
RNG_SEED = 20260629
LINK_EPS = 0.38


@dataclass(frozen=True)
class Cycle:
    vertices: tuple[int, ...]

    @property
    def perimeter(self) -> int:
        return len(self.vertices) - 1


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def gell_mann_basis() -> list[np.ndarray]:
    zero = 0.0 + 0.0j
    one = 1.0 + 0.0j
    i = 1.0j
    return [
        np.array([[zero, one, zero], [one, zero, zero], [zero, zero, zero]], dtype=complex),
        np.array([[zero, -i, zero], [i, zero, zero], [zero, zero, zero]], dtype=complex),
        np.array([[one, zero, zero], [zero, -one, zero], [zero, zero, zero]], dtype=complex),
        np.array([[zero, zero, one], [zero, zero, zero], [one, zero, zero]], dtype=complex),
        np.array([[zero, zero, -i], [zero, zero, zero], [i, zero, zero]], dtype=complex),
        np.array([[zero, zero, zero], [zero, zero, one], [zero, one, zero]], dtype=complex),
        np.array([[zero, zero, zero], [zero, zero, -i], [zero, i, zero]], dtype=complex),
        (1.0 / np.sqrt(3.0)) * np.array(
            [[one, zero, zero], [zero, one, zero], [zero, zero, -2.0 * one]],
            dtype=complex,
        ),
    ]


GELL_MANN = gell_mann_basis()
I3 = np.eye(N_COLOR, dtype=complex)


def su3_exp_from_algebra(coeffs: np.ndarray, eps: float) -> np.ndarray:
    hermitian = sum(float(c) * basis for c, basis in zip(coeffs, GELL_MANN, strict=True))
    hermitian -= np.trace(hermitian) * I3 / N_COLOR
    evals, evecs = np.linalg.eigh(hermitian)
    unitary = evecs @ np.diag(np.exp(1j * eps * evals)) @ evecs.conj().T
    det_phase = np.linalg.det(unitary) ** (1.0 / N_COLOR)
    unitary = unitary / det_phase
    return unitary


def random_su3_near_identity(rng: np.random.Generator, eps: float = LINK_EPS) -> np.ndarray:
    coeffs = rng.normal(size=len(GELL_MANN))
    return su3_exp_from_algebra(coeffs, eps)


def random_su3_haar(rng: np.random.Generator) -> np.ndarray:
    z = rng.normal(size=(N_COLOR, N_COLOR)) + 1j * rng.normal(size=(N_COLOR, N_COLOR))
    q, r = np.linalg.qr(z)
    phases = np.diag(r) / np.abs(np.diag(r))
    q = q @ np.diag(phases.conj())
    det_phase = np.linalg.det(q) ** (1.0 / N_COLOR)
    return q / det_phase


def check_su3_links(links: list[np.ndarray]) -> float:
    max_unitarity = 0.0
    max_det = 0.0
    for link in links:
        max_unitarity = max(max_unitarity, float(np.linalg.norm(link.conj().T @ link - I3)))
        max_det = max(max_det, abs(np.linalg.det(link) - 1.0))
    assert_less("max link unitarity residual", max_unitarity, 1e-12)
    assert_less("max |det(U)-1| residual", max_det, 1e-12)
    return max_unitarity


def edge_matrix(complex_: GluedComplex, links: list[np.ndarray], left: int, right: int) -> np.ndarray:
    edge = (left, right) if left < right else (right, left)
    link = links[complex_.edge_index[edge]]
    if edge == (left, right):
        return link
    return link.conj().T


def cycle_holonomy(complex_: GluedComplex, links: list[np.ndarray], cycle: Cycle) -> np.ndarray:
    out = I3.copy()
    for left, right in zip(cycle.vertices[:-1], cycle.vertices[1:], strict=True):
        out = out @ edge_matrix(complex_, links, left, right)
    return out


def normalized_trace(matrix: np.ndarray) -> float:
    return float(np.real(np.trace(matrix)) / N_COLOR)


def patch_boundary_mask(complex_: GluedComplex, patch: tuple[int, ...]) -> int:
    mask = 0
    for face_i in patch:
        mask ^= complex_.face_masks[face_i]
    return mask


def cycles_from_mask(complex_: GluedComplex, mask: int) -> list[Cycle]:
    adj: dict[int, list[int]] = defaultdict(list)
    for edge_i, (left, right) in enumerate(complex_.edges):
        if mask & (1 << edge_i):
            adj[left].append(right)
            adj[right].append(left)

    if not adj:
        return []
    bad_degrees = {v: len(neighbours) for v, neighbours in adj.items() if len(neighbours) != 2}
    if bad_degrees:
        raise AssertionError(f"boundary is not a disjoint simple-cycle set: {bad_degrees}")

    seen_edges: set[tuple[int, int]] = set()
    cycles: list[Cycle] = []
    for start in sorted(adj):
        for first_next in sorted(adj[start]):
            edge = tuple(sorted((start, first_next)))
            if edge in seen_edges:
                continue
            path = [start, first_next]
            seen_edges.add(edge)
            prev = start
            cur = first_next
            while cur != start:
                candidates = [nxt for nxt in sorted(adj[cur]) if nxt != prev]
                if not candidates:
                    raise AssertionError("open path while walking closed boundary")
                nxt = candidates[0]
                edge = tuple(sorted((cur, nxt)))
                if edge in seen_edges and nxt != start:
                    raise AssertionError("self-intersecting boundary cycle encountered")
                seen_edges.add(edge)
                path.append(nxt)
                prev, cur = cur, nxt
            cycles.append(Cycle(vertices=tuple(path)))
    return cycles


def face_cycle(complex_: GluedComplex, face_i: int) -> Cycle:
    vertices = complex_.faces[face_i].vertices
    return Cycle(vertices=vertices + (vertices[0],))


def gauge_transform_links(
    complex_: GluedComplex,
    links: list[np.ndarray],
    gauges: list[np.ndarray],
) -> list[np.ndarray]:
    out = []
    for edge_i, (left, right) in enumerate(complex_.edges):
        out.append(gauges[left] @ links[edge_i] @ gauges[right].conj().T)
    return out


def max_adjacent_commutator_norm(complex_: GluedComplex, links: list[np.ndarray]) -> float:
    values = []
    for face in complex_.faces:
        edge_indices = []
        vertices = face.vertices
        for left, right in zip(vertices, vertices[1:] + vertices[:1], strict=True):
            edge = (left, right) if left < right else (right, left)
            edge_indices.append(complex_.edge_index[edge])
        for a, b in zip(edge_indices, edge_indices[1:] + edge_indices[:1], strict=True):
            values.append(float(np.linalg.norm(links[a] @ links[b] - links[b] @ links[a])))
    return max(values)


def one_component_patches(complex_: GluedComplex, max_size: int) -> list[tuple[int, ...]]:
    patches = []
    for patch in connected_face_patches(complex_, max_size=max_size):
        if len(patch) == 0:
            continue
        cycles = cycles_from_mask(complex_, patch_boundary_mask(complex_, patch))
        if len(cycles) == 1:
            patches.append(patch)
    return patches


def boundary_trace_for_patch(
    complex_: GluedComplex,
    links: list[np.ndarray],
    patch: tuple[int, ...],
) -> tuple[float, np.ndarray]:
    cycles = cycles_from_mask(complex_, patch_boundary_mask(complex_, patch))
    if len(cycles) != 1:
        raise AssertionError("expected a one-component boundary")
    hol = cycle_holonomy(complex_, links, cycles[0])
    return normalized_trace(hol), hol


def face_product_holonomy(
    complex_: GluedComplex,
    links: list[np.ndarray],
    patch: tuple[int, ...],
) -> np.ndarray:
    out = I3.copy()
    # Deterministic order.  This is deliberately the naive face-product
    # shortcut.  For non-abelian links it is not a gauge-covariant surface
    # transport rule, and the script verifies that it fails for multi-face
    # patches.
    for face_i in sorted(patch):
        out = out @ cycle_holonomy(complex_, links, face_cycle(complex_, face_i))
    return out


def main() -> None:
    print("Glued TCH shell with non-commuting SU(3) link holonomies")
    print("=" * 104)
    print(f"  random seed={RNG_SEED}; link algebra scale={LINK_EPS}; surface q={SURFACE_Q}")

    rng = np.random.default_rng(RNG_SEED)
    complex_ = build_glued_complex(3)
    links = [random_su3_near_identity(rng) for _ in complex_.edges]

    print("\n[1] Link-field sanity and non-commutativity")
    check_su3_links(links)
    comm = max_adjacent_commutator_norm(complex_, links)
    assert_greater("max adjacent SU(3) link commutator norm", comm, 0.02)
    print(
        f"  shell: cells={complex_.cells}, V={complex_.vertices}, "
        f"E={len(complex_.edges)}, F={len(complex_.faces)}"
    )

    print("\n[2] Gauge-invariant boundary Wilson traces")
    patches = one_component_patches(complex_, max_size=3)
    test_patches = [patch for patch in patches if 1 <= len(patch) <= 3][:40]
    gauges = [random_su3_haar(rng) for _ in range(complex_.vertices)]
    transformed = gauge_transform_links(complex_, links, gauges)
    max_gauge_delta = 0.0
    max_surface_geometry_delta = 0.0
    for patch in test_patches:
        before, _ = boundary_trace_for_patch(complex_, links, patch)
        after, _ = boundary_trace_for_patch(complex_, transformed, patch)
        max_gauge_delta = max(max_gauge_delta, abs(before - after))

        # The scalar surface selector remains a geometry ledger: A_min is still
        # the minimal face count for the same boundary.  It is not a replacement
        # for the ordered SU(3) trace.
        row = surface_row(complex_, patch, SURFACE_Q)
        if row is None:
            raise AssertionError("missing surface row for nonempty patch")
        max_surface_geometry_delta = max(
            max_surface_geometry_delta,
            abs(row.min_area - min(len(patch), len(complex_.faces) - len(patch))),
        )
    assert_less("max Wilson trace change under local SU(3) gauge rotations", max_gauge_delta, 1e-11)
    assert_less("surface geometry still returns the same minimal shell area", max_surface_geometry_delta, 1e-12)

    print("\n[3] One-face identity: face holonomy equals boundary holonomy")
    one_face_patches = [patch for patch in patches if len(patch) == 1][:12]
    max_one_face_delta = 0.0
    for patch in one_face_patches:
        boundary_trace, _ = boundary_trace_for_patch(complex_, links, patch)
        face_trace = normalized_trace(face_product_holonomy(complex_, links, patch))
        max_one_face_delta = max(max_one_face_delta, abs(boundary_trace - face_trace))
    assert_less("max one-face boundary-vs-face trace delta", max_one_face_delta, 1e-12)

    print("\n[4] Multi-face scalar shortcut fails for non-commuting links")
    multi_patches = [patch for patch in patches if 2 <= len(patch) <= 3][:60]
    mismatches = []
    for patch in multi_patches:
        boundary_trace, _ = boundary_trace_for_patch(complex_, links, patch)
        naive_face_trace = normalized_trace(face_product_holonomy(complex_, links, patch))
        mismatches.append(abs(boundary_trace - naive_face_trace))

    mismatch_median = float(np.median(mismatches))
    mismatch_max = float(np.max(mismatches))
    mismatch_count = sum(delta > 1e-3 for delta in mismatches)
    print(
        f"  tested multi-face patches={len(multi_patches)}, "
        f"median |delta trace|={mismatch_median:.6f}, "
        f"max |delta trace|={mismatch_max:.6f}, "
        f"failures above 1e-3={mismatch_count}"
    )
    assert_greater("median multi-face scalar-shortcut mismatch", mismatch_median, 1e-3)
    assert_true("at least 90% of tested multi-face shortcuts fail visibly", mismatch_count >= int(0.9 * len(multi_patches)))

    print("\n[5] Ledger interpretation")
    print(
        f"  detector record: one ordered SU(3) trace per boundary cycle; "
        f"sampled one-component cycles={len(test_patches)}"
    )
    print(
        "  surface record: the glued-shell minimal-area selector remains the "
        "geometry ledger, but a non-abelian Wilson action/path integral is "
        "needed to turn it into an averaged area law."
    )

    print(
        """
VERDICT:
  PASS.  Non-commuting SU(3) link matrices can be placed directly on the glued
  TCH shell.  The detector-readable object is the path-ordered boundary Wilson
  trace, and it is locally gauge-invariant.  The one-face check closes exactly.

  The multi-face check deliberately fails for the naive independent-face-scalar
  shortcut.  That is the useful result: once SU(3) links are non-commuting, the
  previous scalar surface selector remains a geometry/action ledger, not a
  substitute for a based non-abelian surface transport or a Wilson-action
  path-integral average.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
