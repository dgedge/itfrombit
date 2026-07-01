#!/usr/bin/env python3
r"""Bulk baryon Y-string versus meson string on the bond-bipyramid graph.

Purpose
-------
The closed-shell Gauss/triality script proves the kinematic colour rule:
only total-triality-zero source sectors have gauge-invariant records.  The
bond-bipyramid static-potential scripts then license a genuine finite bulk
geometry and a leading strong-coupling string tension.

This script combines those two facts for the first baryon-in-bulk check.

It uses the bond graph underlying the corrected bond-centred bipyramid bulk.
A static meson is a q qbar line.  A static baryon is a q q q epsilon junction:
three fundamental Wilson arms meet at one SU(3) epsilon vertex.  The geometry
question is then a three-terminal Steiner problem on the same graph.

For the tri-axial source family

    J  = (1, 1, 1),
    q1 = J + n xhat,
    q2 = J + n yhat,
    q3 = J + n zhat,

the exact graph result is

    L_M = d(q1, q2) = 2n,
    L_Y = min_J' [d(q1,J') + d(q2,J') + d(q3,J')] = 3n,
    L_Y / L_M = 3/2.

The leading strong-coupling Wilson cost uses the same string tension per unit
length as the meson:

    V_M = 4 sigma L_M,
    V_B = 4 sigma L_Y,
    sigma = -log(beta/18).

Boundary
--------
This is a finite graph/leading-strong-coupling result.  It proves that the
bulk record grammar supports a baryonic epsilon/Y operator with the same string
tension as the mesonic line, and gives an N-independent tri-axial ratio 3/2.
It does not compute the proton mass, spin/flavour wavefunction, Coulomb term,
Luscher term, weak-coupling scaling, or the continuum equilateral ratio.  On
this axis-aligned finite graph the Y length equals the half-perimeter length,
so this geometry does not by itself distinguish the continuum Y ansatz from a
Delta ansatz; the Y character comes from the epsilon-junction operator.
"""

from __future__ import annotations

import heapq
import math
from collections import deque
from dataclasses import dataclass

import numpy as np

from record_grammar_tch_bond_bipyramid_bulk import AXES, Bond, bonds_box, build_bond_complex
from record_grammar_tch_bond_bipyramid_static_potential import BETA, u_fundamental
from record_grammar_tch_colour_singlet_gauss import epsilon_tensor
from record_grammar_tch_glued_su3_link_holonomy import random_su3_haar


Vertex = tuple[int, int, int]
MAX_N = 4
BOX_DIMS = (MAX_N + 2, MAX_N + 2, MAX_N + 2)
JUNCTION = (1, 1, 1)
RNG_SEED = 20260630


@dataclass(frozen=True)
class StringRow:
    n: int
    meson_length: int
    baryon_y_length: int
    half_perimeter: int
    witness_junction: Vertex
    meson_potential: float
    baryon_potential: float


def assert_equal(name: str, value: int | Vertex, target: int | Vertex) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_true(name: str, value: bool) -> None:
    print(f"  {name:<78s} value={value}")
    if not value:
        raise AssertionError(name)


def add_axis(vertex: Vertex, axis: int, amount: int) -> Vertex:
    out = list(vertex)
    out[axis] += amount
    return tuple(out)  # type: ignore[return-value]


def source_triple(n: int) -> tuple[Vertex, Vertex, Vertex]:
    return tuple(add_axis(JUNCTION, axis, n) for axis in AXES)  # type: ignore[return-value]


def adjacency_from_bonds(bonds: list[Bond]) -> dict[Vertex, set[Vertex]]:
    adj: dict[Vertex, set[Vertex]] = {}
    for coord, axis in bonds:
        left = coord
        right = add_axis(coord, axis, 1)
        adj.setdefault(left, set()).add(right)
        adj.setdefault(right, set()).add(left)
    return adj


def shortest_distances(adj: dict[Vertex, set[Vertex]], start: Vertex) -> dict[Vertex, int]:
    dist = {start: 0}
    queue: deque[Vertex] = deque([start])
    while queue:
        cur = queue.popleft()
        for nxt in adj[cur]:
            if nxt not in dist:
                dist[nxt] = dist[cur] + 1
                queue.append(nxt)
    return dist


def steiner_three_terminal(adj: dict[Vertex, set[Vertex]], terminals: tuple[Vertex, Vertex, Vertex]) -> tuple[int, Vertex]:
    """Exact three-terminal Steiner length by subset dynamic programming."""

    vertices = tuple(sorted(adj))
    full = (1 << len(terminals)) - 1
    inf = 10**9
    dist = {mask: {vertex: inf for vertex in vertices} for mask in range(1, full + 1)}
    for i, terminal in enumerate(terminals):
        dist[1 << i][terminal] = 0

    for mask in range(1, full + 1):
        sub = (mask - 1) & mask
        while sub:
            other = mask ^ sub
            if other:
                for vertex in vertices:
                    candidate = dist[sub][vertex] + dist[other][vertex]
                    if candidate < dist[mask][vertex]:
                        dist[mask][vertex] = candidate
            sub = (sub - 1) & mask

        heap = [(value, vertex) for vertex, value in dist[mask].items() if value < inf]
        heapq.heapify(heap)
        while heap:
            value, vertex = heapq.heappop(heap)
            if value != dist[mask][vertex]:
                continue
            for nxt in adj[vertex]:
                candidate = value + 1
                if candidate < dist[mask][nxt]:
                    dist[mask][nxt] = candidate
                    heapq.heappush(heap, (candidate, nxt))

    witness = min(vertices, key=lambda vertex: dist[full][vertex])
    return dist[full][witness], witness


def epsilon_junction_invariance() -> float:
    rng = np.random.default_rng(RNG_SEED)
    g = random_su3_haar(rng)
    eps = epsilon_tensor()
    transformed = np.einsum("ia,jb,kc,abc->ijk", g, g, g, eps)
    return float(np.linalg.norm(transformed - eps))


def main() -> None:
    print("Bulk baryon Y-string versus meson string")
    print("=" * 100)
    sigma = -math.log(u_fundamental(BETA))
    unit_cost = 4.0 * sigma
    print(f"  bond graph box={BOX_DIMS}; beta={BETA:.3f}; sigma=-ln(beta/18)={sigma:.9f}")
    print(f"  Wilson time-step face cost per spatial graph edge = 4 sigma = {unit_cost:.9f}")

    bonds = bonds_box(*BOX_DIMS)
    complex_ = build_bond_complex("bond-baryon-y-box", bonds)
    adj = adjacency_from_bonds(bonds)
    assert_equal("bond graph vertex count matches box", len(adj), BOX_DIMS[0] * BOX_DIMS[1] * BOX_DIMS[2])
    assert_equal("bond-bipyramid cell count matches bond graph edges", len(complex_.cell_masks), len(bonds))
    assert_true("junction lies in bulk graph", JUNCTION in adj)

    print("\n[1] SU(3) epsilon junction")
    eps_error = epsilon_junction_invariance()
    assert_close("epsilon_ijk invariant under common SU(3) junction rotation", eps_error, 0.0, tol=1e-12)
    print("  interpretation: three Wilson arms can meet at one gauge-invariant epsilon junction.")

    print("\n[2] Exact graph strings")
    print("  n  L_M=d(q1,q2)  L_Y  half-perimeter  junction      V_M        V_B       ratio")
    rows: list[StringRow] = []
    for n in range(1, MAX_N + 1):
        q1, q2, q3 = source_triple(n)
        terminals = (q1, q2, q3)
        distances = {terminal: shortest_distances(adj, terminal) for terminal in terminals}
        d12 = distances[q1][q2]
        d23 = distances[q2][q3]
        d31 = distances[q3][q1]
        half_perimeter = (d12 + d23 + d31) // 2
        steiner_length, witness = steiner_three_terminal(adj, terminals)

        expected_meson = 2 * n
        expected_y = 3 * n
        assert_equal(f"n={n} q1-q2 meson graph length", d12, expected_meson)
        assert_equal(f"n={n} pairwise graph lengths equal", d23, expected_meson)
        assert_equal(f"n={n} pairwise graph lengths equal", d31, expected_meson)
        assert_equal(f"n={n} exact baryon Y/Steiner length", steiner_length, expected_y)
        assert_equal(f"n={n} exact Y witness junction", witness, JUNCTION)
        assert_equal(f"n={n} half-perimeter control", half_perimeter, expected_y)

        meson_potential = unit_cost * d12
        baryon_potential = unit_cost * steiner_length
        assert_close(f"n={n} V_B/V_M", baryon_potential / meson_potential, 1.5)
        rows.append(
            StringRow(
                n=n,
                meson_length=d12,
                baryon_y_length=steiner_length,
                half_perimeter=half_perimeter,
                witness_junction=witness,
                meson_potential=meson_potential,
                baryon_potential=baryon_potential,
            )
        )
        print(
            f"  {n:>1d}"
            f" {d12:>13d}"
            f" {steiner_length:>4d}"
            f" {half_perimeter:>15d}"
            f" {str(witness):>12s}"
            f" {meson_potential:>9.3f}"
            f" {baryon_potential:>9.3f}"
            f" {baryon_potential / meson_potential:>8.3f}"
        )

    print("\n[3] Leading strong-coupling form")
    for row in rows:
        assert_close(f"n={row.n} meson V/L", row.meson_potential / row.meson_length, unit_cost)
        assert_close(f"n={row.n} baryon V/L_Y", row.baryon_potential / row.baryon_y_length, unit_cost)

    print(
        """
VERDICT:
  PASS.  The licensed bond-bipyramid bulk supports the expected baryonic
  epsilon/Y record: three fundamental Wilson arms meet at a gauge-invariant
  SU(3) epsilon junction.  The exact three-terminal Steiner calculation on
  the same bond graph gives L_Y=3n for the tri-axial source family, while the
  matched meson length is L_M=2n.  Therefore V_B/V_M=3/2 for this finite
  graph family, independent of n.

  More structurally, the per-unit string tension is the same in the meson and
  baryon sectors: V_M=4 sigma L_M and V_B=4 sigma L_Y.  This is the finite
  record-grammar analogue of the lattice-QCD Y-ansatz form
  V_3Q = sigma L_min + ... .

  Boundary: this is not a proton mass calculation.  The tri-axial graph metric
  also makes the half-perimeter control equal to L_Y, so this source family
  does not distinguish the continuum Y ansatz from a Delta ansatz.  The
  continuum/equilateral, Coulomb, Luscher, spin/flavour, and weak-coupling
  corrections remain outside this finite leading-order gate.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
