#!/usr/bin/env python3
r"""Closed-shell TCH colour singlets from Gauss law and triality.

Purpose
-------
This is the exact, non-perturbative companion to the glued-shell Wilson-action
average.

The Wilson-action script proves a leading strong-coupling area law on the glued
TCH shell.  This script proves the simpler but deeper selection rule that does
not depend on strong coupling:

    a closed oriented colour shell has no gauge-invariant record with net
    SU(3) triality.

In centre language, a fundamental quark has triality +1, an antiquark has
triality -1, and an epsilon baryon vertex combines three fundamentals to
triality 0.  A compact shell has no boundary where net colour flux can escape,
so Gauss's law requires total triality 0.

What is proved here
-------------------
1. The glued TCH shell incidence matrix has rank V-1 over GF(3).  Therefore a
   source vector is in the Gauss-law image iff its total triality is zero.
2. A lone fundamental source is forbidden: no Z3 flux assignment solves Gauss.
3. A q qbar meson has an explicit finite open-line record.
4. A q q q epsilon baryon has an explicit finite Y/triangle record.
5. The SU(3) invariant tensors agree with the triality result:
      3 has no singlet, 3 x 3bar has delta, 3 x 3 x 3 has epsilon.
6. The finite-shell cap numbers are purely geometric.  For N glued cells,
      F = 12N + 2, half-shell = 6N + 1,
   while the local epsilon baryon cap is one triangular TCH face.  The resulting
   declared meson:baryon unit-cap ratio is (6N+1):1 under this shell convention.

Boundary
--------
This proves colour-singlet confinement as a closed-cell record-selection rule:
nonzero triality cannot be a gauge-invariant record on the compact shell.
It does not compute the physical proton mass, the baryon spectrum, or the QCD
string tension.  The area-ratio line is a finite-shell geometry diagnostic, not
a claim about physical hadron masses.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_glued_surface_selector import build_glued_complex
from record_grammar_tch_glued_su3_link_holonomy import random_su3_haar


TRIALITY = 3
RNG_SEED = 20260630


@dataclass(frozen=True)
class FluxRecord:
    charges: np.ndarray
    flux: np.ndarray
    support_edges: int


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def assert_false(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if value:
        raise AssertionError(name)


def gf3_rank(matrix: np.ndarray) -> int:
    a = np.array(matrix, dtype=int) % TRIALITY
    rows, cols = a.shape
    rank = 0
    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if a[row, col] % TRIALITY:
                pivot = row
                break
        if pivot is None:
            continue
        if pivot != rank:
            a[[rank, pivot]] = a[[pivot, rank]]
        inv = 1 if a[rank, col] == 1 else 2
        a[rank] = (inv * a[rank]) % TRIALITY
        for row in range(rows):
            if row != rank and a[row, col] % TRIALITY:
                a[row] = (a[row] - a[row, col] * a[rank]) % TRIALITY
        rank += 1
        if rank == rows:
            break
    return rank


def incidence_matrix(complex_) -> np.ndarray:
    mat = np.zeros((complex_.vertices, len(complex_.edges)), dtype=int)
    for edge_i, (left, right) in enumerate(complex_.edges):
        # Canonical orientation left -> right.  Divergence convention:
        # outgoing = -1, incoming = +1, all modulo 3.
        mat[left, edge_i] = 2
        mat[right, edge_i] = 1
    return mat


def in_gauss_image(complex_, charges: np.ndarray) -> bool:
    b = incidence_matrix(complex_)
    target = (-charges) % TRIALITY
    return gf3_rank(b) == gf3_rank(np.column_stack([b, target]))


def divergence(complex_, flux: np.ndarray) -> np.ndarray:
    return (incidence_matrix(complex_) @ flux) % TRIALITY


def adjacency(complex_) -> list[list[tuple[int, int]]]:
    adj = [[] for _ in range(complex_.vertices)]
    for edge_i, (left, right) in enumerate(complex_.edges):
        adj[left].append((right, edge_i))
        adj[right].append((left, edge_i))
    return adj


def shortest_path(complex_, start: int, stop: int) -> list[int]:
    adj = adjacency(complex_)
    parent: dict[int, int | None] = {start: None}
    queue = deque([start])
    while queue:
        cur = queue.popleft()
        if cur == stop:
            break
        for nxt, _ in adj[cur]:
            if nxt not in parent:
                parent[nxt] = cur
                queue.append(nxt)
    if stop not in parent:
        raise AssertionError("glued shell graph is disconnected")
    out = [stop]
    while out[-1] != start:
        prev = parent[out[-1]]
        if prev is None:
            raise AssertionError("bad path parent")
        out.append(prev)
    return list(reversed(out))


def add_path_flux(complex_, flux: np.ndarray, path: list[int], amount: int = 1) -> None:
    for left, right in zip(path[:-1], path[1:], strict=True):
        edge = (left, right) if left < right else (right, left)
        edge_i = complex_.edge_index[edge]
        if edge == (left, right):
            flux[edge_i] = (flux[edge_i] + amount) % TRIALITY
        else:
            flux[edge_i] = (flux[edge_i] - amount) % TRIALITY


def meson_record(complex_, quark: int, antiquark: int) -> FluxRecord:
    charges = np.zeros(complex_.vertices, dtype=int)
    charges[quark] = 1
    charges[antiquark] = 2
    flux = np.zeros(len(complex_.edges), dtype=int)
    add_path_flux(complex_, flux, shortest_path(complex_, quark, antiquark), amount=1)
    return FluxRecord(charges=charges, flux=flux, support_edges=int(np.count_nonzero(flux)))


def baryon_record(complex_, quarks: tuple[int, int, int], epsilon_vertex: int) -> FluxRecord:
    charges = np.zeros(complex_.vertices, dtype=int)
    for vertex in quarks:
        charges[vertex] = (charges[vertex] + 1) % TRIALITY
    flux = np.zeros(len(complex_.edges), dtype=int)
    for vertex in quarks:
        if vertex != epsilon_vertex:
            add_path_flux(complex_, flux, shortest_path(complex_, vertex, epsilon_vertex), amount=1)
    return FluxRecord(charges=charges, flux=flux, support_edges=int(np.count_nonzero(flux)))


def center_projector_value(total_triality: int) -> complex:
    omega = np.exp(2j * np.pi / TRIALITY)
    return sum(omega ** (k * (total_triality % TRIALITY)) for k in range(TRIALITY)) / TRIALITY


def epsilon_tensor() -> np.ndarray:
    eps = np.zeros((3, 3, 3), dtype=complex)
    for perm in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        eps[perm] = 1.0
    for perm in ((0, 2, 1), (2, 1, 0), (1, 0, 2)):
        eps[perm] = -1.0
    return eps / np.sqrt(6.0)


def invariant_tensor_checks() -> tuple[float, float, float]:
    rng = np.random.default_rng(RNG_SEED)
    u = random_su3_haar(rng)

    delta = np.eye(3, dtype=complex) / np.sqrt(3.0)
    delta_transformed = u @ delta @ u.conj().T
    delta_err = float(np.linalg.norm(delta_transformed - delta))

    eps = epsilon_tensor()
    eps_transformed = np.einsum("ia,jb,kc,abc->ijk", u, u, u, eps)
    eps_err = float(np.linalg.norm(eps_transformed - eps))

    lone = np.array([1.0, 0.0, 0.0], dtype=complex)
    lone_err = float(np.linalg.norm(u @ lone - lone))
    return delta_err, eps_err, lone_err


def first_triangle_vertices(complex_) -> tuple[int, int, int]:
    for face in complex_.faces:
        if ":T" in face.label:
            return face.vertices
    raise AssertionError("no triangular TCH face found")


def main() -> None:
    print("Closed-shell TCH colour-singlet Gauss/triality certificate")
    print("=" * 104)

    print("\n[1] Closed-shell Gauss law over the SU(3) centre")
    for cells in range(1, 7):
        complex_ = build_glued_complex(cells)
        b = incidence_matrix(complex_)
        assert_equal(f"N={cells} GF(3) incidence rank", gf3_rank(b), complex_.vertices - 1)

        lone = np.zeros(complex_.vertices, dtype=int)
        lone[0] = 1
        assert_false(f"N={cells} lone fundamental is Gauss-solvable", in_gauss_image(complex_, lone))

        neutral = np.zeros(complex_.vertices, dtype=int)
        neutral[0] = 1
        neutral[-1] = 2
        assert_true(f"N={cells} q qbar source is Gauss-solvable", in_gauss_image(complex_, neutral))

        tri = np.zeros(complex_.vertices, dtype=int)
        for vertex in first_triangle_vertices(complex_):
            tri[vertex] = 1
        assert_true(f"N={cells} q q q source is Gauss-solvable", in_gauss_image(complex_, tri))

    print("\n[2] Explicit finite confined records")
    complex_ = build_glued_complex(4)
    meson = meson_record(complex_, quark=0, antiquark=complex_.vertices - 1)
    baryon_vertices = first_triangle_vertices(complex_)
    baryon = baryon_record(complex_, baryon_vertices, epsilon_vertex=baryon_vertices[0])

    assert_equal("meson total triality", int(np.sum(meson.charges) % TRIALITY), 0)
    assert_equal("meson Gauss residual", int(np.max((divergence(complex_, meson.flux) + meson.charges) % TRIALITY)), 0)
    assert_true("meson finite open-line record has nonzero flux support", meson.support_edges > 0)

    assert_equal("baryon total triality", int(np.sum(baryon.charges) % TRIALITY), 0)
    assert_equal("baryon Gauss residual", int(np.max((divergence(complex_, baryon.flux) + baryon.charges) % TRIALITY)), 0)
    assert_true("baryon epsilon/Y record has finite flux support", baryon.support_edges > 0)
    print(f"  meson path support edges: {meson.support_edges}")
    print(f"  baryon epsilon support edges: {baryon.support_edges}")

    print("\n[3] SU(3) singlet tensor gate")
    center_lone = center_projector_value(1)
    center_meson = center_projector_value(1 + 2)
    center_baryon = center_projector_value(1 + 1 + 1)
    assert_less("global centre projection of lone fundamental", abs(center_lone), 1e-12)
    assert_less("global centre projection of q qbar minus 1", abs(center_meson - 1.0), 1e-12)
    assert_less("global centre projection of q q q minus 1", abs(center_baryon - 1.0), 1e-12)

    delta_err, eps_err, lone_err = invariant_tensor_checks()
    assert_less("delta_ij meson singlet invariance", delta_err, 1e-12)
    assert_less("epsilon_ijk baryon singlet invariance", eps_err, 1e-12)
    assert_true("generic SU(3) rotation moves a lone fundamental vector", lone_err > 1e-2)

    print("\n[4] Pure shell cap geometry")
    print("  cells   faces=12N+2   half-shell=6N+1   baryon T-face cap   meson:baryon")
    for cells in (1, 2, 3, 4, 6, 8, 10):
        complex_n = build_glued_complex(cells)
        faces = len(complex_n.faces)
        half_shell = faces // 2
        baryon_cap = 1
        assert_equal(f"N={cells} face count", faces, 12 * cells + 2)
        assert_equal(f"N={cells} half-shell cap", half_shell, 6 * cells + 1)
        print(f"  {cells:>5d}   {faces:>10d}   {half_shell:>14d}   {baryon_cap:>16d}   {half_shell}:1")

    print(
        """
VERDICT:
  PASS.  On a closed glued TCH colour shell, the GF(3) incidence rank is V-1,
  so Gauss's law admits exactly the total-triality-zero source sectors.  A lone
  fundamental source has no gauge-invariant record.  The q qbar meson and the
  q q q epsilon baryon have explicit finite records, and the SU(3) delta and
  epsilon invariant tensors verify the same statement in the full group.

  This is exact closed-cell colour-singlet confinement: orientation-complete
  means triality-neutral.  It is independent of the strong-coupling expansion.
  The cap-ratio line is a shell-geometry diagnostic, not a hadron-mass claim.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
