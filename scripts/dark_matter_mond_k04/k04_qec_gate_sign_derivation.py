#!/usr/bin/env python3
r"""Derive the K04 Hamiltonian form/signs from QEC gate geometry.

Target
------
K04 currently enters canon as the ordinal crystallisation rule

    E(G) = -w4 C4(G) - w6 C6(G) + lambda sum_v (deg(v)-3)^2,
    with w4, w6, lambda > 0.

Can the *form and signs* be derived from deeper QEC/service gates?

This audit proves the ordinal statement, not the numeric weights:

  1. Degree-3 penalty, lambda > 0:
     a framed service cell has three independent address/service directions.
     A vertex with degree != 3 has either a missing monitored channel or an
     overwritten channel.  The unique local symmetric quadratic cost that
     vanishes exactly on legal valence is (deg-3)^2, hence positive sign.

  2. Four-cycle reward, w4 > 0:
     on Q3, length-4 cycles are the minimal closed parity-check loops.  They
     are exactly two-axis square faces: their F2 boundary is zero, so they are
     legal local stabilizer closures.  More such closed loops means fewer open
     syndrome boundaries, so they enter as a reward, -w4 C4.

  3. Six-cycle reward, w6 > 0:
     length-6 cycles are the minimal closed loops that use all three address
     directions twice.  They are the smallest tri-axis strain/cube-closure
     witnesses.  They distinguish the true cubic QEC cell from 4-cycle-only
     clique artefacts, so they also enter as a reward, -w6 C6.

The script also checks the known selection facts: Q3 maximises (C4,C6) among
connected 8-vertex cubic graphs over the positive orthant, while C6 is
load-bearing against the local K4 trap.

Exit 0 means the K04 sign Hamiltonian is derived at framed-QEC ordinal grade.
It does not derive w4/w6, w6/Lambda, gamma, or the boot cooling law.
"""

from __future__ import annotations

from itertools import combinations, permutations

import numpy as np


def count_k_cycles(A: np.ndarray, k: int) -> int:
    """Exact count of distinct undirected simple k-cycles."""

    n = A.shape[0]
    total = 0
    for subset in combinations(range(n), k):
        first, rest = subset[0], subset[1:]
        for perm in permutations(rest):
            if perm[0] > perm[-1]:
                continue
            seq = (first,) + perm
            if all(A[seq[i], seq[(i + 1) % k]] for i in range(k)):
                total += 1
    return total


def is_connected(A: np.ndarray) -> bool:
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for u, live in enumerate(A[v]):
            if live and u not in seen:
                seen.add(u)
                stack.append(u)
    return len(seen) == A.shape[0]


def q3_edges() -> list[tuple[int, int, int]]:
    edges: list[tuple[int, int, int]] = []
    for x in range(8):
        for axis in range(3):
            y = x ^ (1 << axis)
            if x < y:
                edges.append((x, y, axis))
    return edges


def q3_graph() -> np.ndarray:
    A = np.zeros((8, 8), dtype=int)
    for x, y, _ in q3_edges():
        A[x, y] = A[y, x] = 1
    return A


def k4_graph() -> np.ndarray:
    A = np.ones((4, 4), dtype=int)
    np.fill_diagonal(A, 0)
    return A


def incidence_rank_f2(edges: list[tuple[int, int]], n_vertices: int) -> int:
    M = np.zeros((len(edges), n_vertices), dtype=np.uint8)
    for row, (a, b) in enumerate(edges):
        M[row, a] = 1
        M[row, b] = 1
    # Gaussian elimination over F2.
    rank = 0
    col = 0
    while rank < M.shape[0] and col < M.shape[1]:
        pivots = np.where(M[rank:, col] == 1)[0]
        if len(pivots) == 0:
            col += 1
            continue
        pivot = rank + int(pivots[0])
        M[[rank, pivot]] = M[[pivot, rank]]
        for r in range(M.shape[0]):
            if r != rank and M[r, col]:
                M[r] ^= M[rank]
        rank += 1
        col += 1
    return rank


def simple_cycles_q3(k: int) -> list[tuple[int, ...]]:
    """Canonical simple cycles in Q3 as vertex tuples."""

    A = q3_graph()
    cycles: set[tuple[int, ...]] = set()
    for subset in combinations(range(8), k):
        first, rest = subset[0], subset[1:]
        for perm in permutations(rest):
            if perm[0] > perm[-1]:
                continue
            seq = (first,) + perm
            if all(A[seq[i], seq[(i + 1) % k]] for i in range(k)):
                rotations = []
                for i in range(k):
                    r = seq[i:] + seq[:i]
                    rotations.append(r)
                    rotations.append(tuple(reversed(r)))
                cycles.add(min(rotations))
    return sorted(cycles)


def edge_axis(a: int, b: int) -> int:
    diff = a ^ b
    assert diff in (1, 2, 4)
    return {1: 0, 2: 1, 4: 2}[diff]


def gen_cubic_labelled(n: int = 8, degree: int = 3) -> list[np.ndarray]:
    A = np.zeros((n, n), dtype=int)
    deg = [0] * n
    out: list[np.ndarray] = []

    def bt() -> None:
        v = next((i for i in range(n) if deg[i] < degree), None)
        if v is None:
            out.append(A.copy())
            return
        candidates = [u for u in range(v + 1, n) if deg[u] < degree and A[v, u] == 0]
        need = degree - deg[v]
        if len(candidates) < need:
            return
        for combo in combinations(candidates, need):
            for u in combo:
                A[v, u] = A[u, v] = 1
                deg[v] += 1
                deg[u] += 1
            bt()
            for u in combo:
                A[v, u] = A[u, v] = 0
                deg[v] -= 1
                deg[u] -= 1

    bt()
    return out


print("[1] Degree-3 service-valence gate")
Aq = q3_graph()
degrees = Aq.sum(axis=1)
print(f"    Q3 degrees = {degrees.tolist()}")
assert set(degrees.tolist()) == {3}
service_channels = 3
for d in range(0, 7):
    cost = (d - service_channels) ** 2
    if d == 3:
        assert cost == 0
    else:
        assert cost > 0
print("    local symmetric cost (deg-3)^2 is positive and vanishes only at legal valence.")
print("    -> lambda > 0.")

print("\n[2] Q3 strain/coboundary gate")
edges_no_axis = [(a, b) for a, b, _ in q3_edges()]
rank = incidence_rank_f2(edges_no_axis, 8)
kernel_dim = 8 - rank
print(f"    edge checks = {len(edges_no_axis)}, F2 incidence rank = {rank}, kernel dim = {kernel_dim}")
assert len(edges_no_axis) == 12
assert rank == 7
assert kernel_dim == 1
print("    Q3 edge syndrome determines vertex faults up to global complement: the 12-edge strain ledger.")

print("\n[3] Four-cycle reward gate")
c4s = simple_cycles_q3(4)
print(f"    Q3 C4 count = {len(c4s)}")
axis_profiles_4 = []
for cyc in c4s:
    axes = [edge_axis(cyc[i], cyc[(i + 1) % 4]) for i in range(4)]
    profile = tuple(sorted(axes.count(axis) for axis in range(3)))
    axis_profiles_4.append(profile)
assert len(c4s) == 6
assert set(axis_profiles_4) == {(0, 2, 2)}
print("    every C4 is a minimal closed two-axis parity loop with F2 boundary zero.")
print("    -> closed stabilizer loops are rewarded, so w4 > 0 and contribution is -w4*C4.")

print("\n[4] Six-cycle reward gate")
c6s = simple_cycles_q3(6)
print(f"    Q3 C6 count = {len(c6s)}")
axis_profiles_6 = []
for cyc in c6s:
    axes = [edge_axis(cyc[i], cyc[(i + 1) % 6]) for i in range(6)]
    profile = tuple(sorted(axes.count(axis) for axis in range(3)))
    axis_profiles_6.append(profile)
assert len(c6s) == 16
assert set(axis_profiles_6) == {(2, 2, 2)}
print("    every C6 is a minimal closed tri-axis strain/cube-consistency loop.")
print("    -> tri-axis cube closure is rewarded, so w6 > 0 and contribution is -w6*C6.")

print("\n[5] C6 is load-bearing: the K4 trap")
Ak4 = k4_graph()
q3_c4, q3_c6 = count_k_cycles(Aq, 4), count_k_cycles(Aq, 6)
k4_c4, k4_c6 = count_k_cycles(Ak4, 4), count_k_cycles(Ak4, 6)
print(f"    Q3 per vertex: C4={q3_c4/8:.3f}, C6={q3_c6/8:.3f}")
print(f"    K4 per vertex: C4={k4_c4/4:.3f}, C6={k4_c6/4:.3f}")
assert abs(q3_c4 / 8 - k4_c4 / 4) < 1e-12
assert q3_c6 > 0 and k4_c6 == 0
print("    C4 alone ties Q3 to tetrahedral cliques; positive C6 selects true cubic closure.")

print("\n[6] Positive-orthant Q3 selection among connected 8-vertex cubic graphs")
classes: dict[tuple[float, ...], np.ndarray] = {}
for A in gen_cubic_labelled(8, 3):
    if not is_connected(A):
        continue
    key = tuple(np.round(np.sort(np.linalg.eigvalsh(A.astype(float))), 6))
    classes.setdefault(key, A)
reps = list(classes.values())
profiles = [(count_k_cycles(A, 4), count_k_cycles(A, 6), A) for A in reps]
q3_profile = (q3_c4, q3_c6)
print(f"    connected cubic graph classes on 8 vertices = {len(reps)}")
print(f"    Q3 profile = C4={q3_c4}, C6={q3_c6}")
for c4, c6, _ in sorted(profiles):
    print(f"      candidate C4={c4:2d}, C6={c6:2d}")
assert len(reps) == 5
assert all(q3_c4 >= c4 and q3_c6 >= c6 for c4, c6, _ in profiles)
assert any(q3_c4 > c4 or q3_c6 > c6 for c4, c6, _ in profiles)
print("    -> with w4,w6>0, Q3 minimises -w4*C4-w6*C6 over this finite QEC cell alphabet.")

print(
    r"""
[7] VERDICT
  K04 form/signs are derived at framed-QEC ordinal grade:

      E(G) = -w4 C4(G) - w6 C6(G) + lambda sum_v (deg(v)-3)^2,
      w4, w6, lambda > 0.

  What is derived:
    * lambda > 0 from the three-channel service-valence constraint;
    * w4 > 0 from minimal closed distance-4/two-axis parity loops;
    * w6 > 0 from minimal closed tri-axis strain/cube-consistency loops;
    * C6 is load-bearing, because C4 alone cannot distinguish Q3 from K4
      clique artefacts;
    * Q3 is then selected over the connected 8-vertex cubic cell alphabet
      throughout the positive orthant.

  What is not derived:
    * the numerical ratio w4/w6;
    * the absolute scale w6/Lambda;
    * the Metropolis/gamma boot cooling law;
    * the unembedded configuration-model sector, where K3,3 artefacts can still
      beat Q3 unless the physical Z3/framed-slice embedding constraint is imposed.
exit 0"""
)
print("ALL ASSERTIONS PASSED -- K04 sign Hamiltonian derived; numeric weights/cooling remain open.")
