#!/usr/bin/env python3
r"""RT first-law -> linearized-gravity audit on the verified substrate min-cut model.

Motivation
----------
The current gravity sector is still largely horizon-input / Dirac-class.  The
holographic-bridge result is stronger: in the link-state limit the substrate
already carries Ryu-Takayanagi structure,

    S(A) = min-cut(A|Ac).

The Jacobson/Faulkner route to gravity has a precise shape:

    RT area variation      delta S_A = delta Area(gamma_A)/(4G)
    entanglement first law delta S_A = delta <K_A>
    local modular K_A      delta <K_A> = integral_A kernel * delta T_00

for all ball-shaped regions A.  In continuum AdS/CFT, the resulting all-ball
identity is equivalent to the linearized Einstein equation.

This script performs the first finite-substrate test:

1. Verify the RT first variation on the same min-cut lattice: small edge-capacity
   perturbations h_e give delta S_A = sum_{e in gamma_A} h_e.
2. Test whether boundary interval entropies see the full linearized edge metric.
   A single smooth RT branch sees only part of the edge space, but the full
   symmetric RT subgradient family has full column rank.
3. Build the finite first-law constraint operator.  Given a candidate modular
   energy kernel B_{A,i}, first law requires M h in range(B).  Equivalently
       E h = 0,  E = (left-nullspace of B)^T M.
   This is the graph analog of the linearized Einstein constraint.

Exit 0 means the RT/linearized-area side is strong enough and the first-law
constraint operator is explicit.  It also rejects a naive boundary-site modular
kernel as too small: it over-constrains the graph to the zero perturbation.  It
does NOT mean Einstein gravity is fully derived.  The actual QEC modular
Hamiltonian/stress-tensor map B, the 1/4G normalization, and the dS/HBC horizon
continuum lift remain the named open gates.
"""

from __future__ import annotations

import math
from collections import defaultdict, deque
from dataclasses import dataclass

import numpy as np


R = 4
INF = 1.0e9


@dataclass(frozen=True)
class Region:
    start: int
    length: int
    nodes: frozenset[int]


def build_disk(radius: int = R):
    sites = [
        (x, y)
        for x in range(-radius, radius + 1)
        for y in range(-radius, radius + 1)
        if x * x + y * y <= radius * radius
    ]
    idx = {site: i for i, site in enumerate(sites)}
    edges: list[tuple[int, int]] = []
    for x, y in sites:
        for dx, dy in ((1, 0), (0, 1)):
            other = (x + dx, y + dy)
            if other in idx:
                edges.append((idx[(x, y)], idx[other]))
    deg: defaultdict[int, int] = defaultdict(int)
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    boundary = sorted(
        (i for i in range(len(sites)) if deg[i] < 4),
        key=lambda i: math.atan2(sites[i][1], sites[i][0]),
    )
    return sites, edges, boundary


SITES, EDGES, BOUNDARY = build_disk()
N = len(SITES)
NB = len(BOUNDARY)
NE = len(EDGES)


def deterministic_caps(seed: int = 0) -> np.ndarray:
    """Unit capacities with a tiny deterministic tie-breaker.

    The tie-breaker selects one representative RT surface when the square grid
    has exact degeneracies.  It is not physical input; it only makes the
    derivative test well-posed.  The rank/tomography result also holds in the
    unit-capacity probe used in substrate_rt_wedge.py for this radius.
    """

    tags = np.array([
        math.sin(
            (seed + 1) * 0.31 * (a + 1)
            + (seed + 2) * 0.47 * (b + 1)
            + (seed + 3) * 0.19 * (k + 1)
        )
        for k, (a, b) in enumerate(EDGES)
    ])
    return 1.0 + 1.0e-3 * tags


def maxflow_mincut(region: frozenset[int], capacities: np.ndarray) -> tuple[float, list[int]]:
    source, sink = N, N + 1
    total = N + 2
    cap: list[defaultdict[int, float]] = [defaultdict(float) for _ in range(total)]
    for eid, (a, b) in enumerate(EDGES):
        c = float(capacities[eid])
        cap[a][b] += c
        cap[b][a] += c
    for bnode in BOUNDARY:
        if bnode in region:
            cap[source][bnode] = INF
        else:
            cap[bnode][sink] = INF

    flow = 0.0
    while True:
        parent = [-1] * total
        parent[source] = source
        q: deque[int] = deque([source])
        while q:
            u = q.popleft()
            for v, c in cap[u].items():
                if c > 1.0e-12 and parent[v] == -1:
                    parent[v] = u
                    q.append(v)
        if parent[sink] == -1:
            break
        bottleneck = math.inf
        v = sink
        while v != source:
            bottleneck = min(bottleneck, cap[parent[v]][v])
            v = parent[v]
        v = sink
        while v != source:
            cap[parent[v]][v] -= bottleneck
            cap[v][parent[v]] = cap[v].get(parent[v], 0.0) + bottleneck
            v = parent[v]
        flow += bottleneck

    reachable = {source}
    q = deque([source])
    while q:
        u = q.popleft()
        for v, c in cap[u].items():
            if c > 1.0e-12 and v not in reachable:
                reachable.add(v)
                q.append(v)

    cut_edges = [
        eid
        for eid, (a, b) in enumerate(EDGES)
        if (a in reachable) != (b in reachable)
    ]
    return flow, cut_edges


def contiguous_regions() -> list[Region]:
    regions: list[Region] = []
    for start in range(NB):
        for length in range(1, NB):
            nodes = frozenset(BOUNDARY[(start + j) % NB] for j in range(length))
            regions.append(Region(start, length, nodes))
    return regions


def mincut_matrix(regions: list[Region], capacities: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.zeros((len(regions), NE))
    entropies = np.zeros(len(regions))
    for r, region in enumerate(regions):
        value, cut = maxflow_mincut(region.nodes, capacities)
        entropies[r] = value
        matrix[r, cut] = 1.0
    return matrix, entropies


def mincut_matrix_family(regions: list[Region], seeds: int = 8) -> tuple[np.ndarray, list[int]]:
    rows: list[np.ndarray] = []
    branch_ranks: list[int] = []
    for seed in range(seeds):
        matrix, _ = mincut_matrix(regions, deterministic_caps(seed))
        rows.append(matrix)
        branch_ranks.append(int(np.linalg.matrix_rank(matrix, tol=1.0e-10)))
    return np.vstack(rows), branch_ranks


def toy_modular_kernel(regions: list[Region]) -> np.ndarray:
    """Discrete ball-modular-Hamiltonian placeholder.

    For a continuum CFT ball/interval, K_A is a local weighted integral of T_00
    with a positive kernel that vanishes at the boundary.  This parabolic
    boundary kernel is not adopted as canon; it is a stand-in that lets the
    first-law constraint operator be built and rank-tested.
    """

    kernel = np.zeros((len(regions), NB))
    for r, region in enumerate(regions):
        m = region.length
        for j in range(m):
            # Positive midpoint parabolic weight, zero at the interval ends in
            # the continuum limit.
            x = (j + 0.5) / m
            weight = x * (1.0 - x) * m
            kernel[r, (region.start + j) % NB] = weight
    return kernel


def orth_projector_to_column_space(matrix: np.ndarray, tol: float = 1.0e-10) -> tuple[np.ndarray, int, np.ndarray]:
    u, s, _ = np.linalg.svd(matrix, full_matrices=True)
    rank = int(np.sum(s > tol))
    projector = u[:, :rank] @ u[:, :rank].T
    left_null = u[:, rank:]
    return projector, rank, left_null


def main() -> None:
    caps = deterministic_caps(0)
    regions = contiguous_regions()
    M_branch, S0 = mincut_matrix(regions, caps)
    M_unit, _ = mincut_matrix(regions, np.ones(NE))
    M_family, branch_ranks = mincut_matrix_family(regions, seeds=8)
    rank_branch = int(np.linalg.matrix_rank(M_branch, tol=1.0e-10))
    rank_unit = int(np.linalg.matrix_rank(M_unit, tol=1.0e-10))
    rank_family = int(np.linalg.matrix_rank(M_family, tol=1.0e-10))
    covered = int(np.sum(np.sum(M_family, axis=0) > 0))

    print("RT FIRST-LAW -> LINEARIZED-GRAVITY AUDIT")
    print("=" * 96)
    print("[0] substrate RT graph")
    print(f"    disk radius R={R}: {N} sites, {NE} bulk edges, {NB} ordered boundary sites")
    print(f"    contiguous boundary intervals: {len(regions)}")
    print(f"    one generic RT branch rank: {rank_branch}/{NE}")
    print(f"    symmetric unit-capacity representative rank: {rank_unit}/{NE}")
    print(f"    RT subgradient-family ranks by branch: {branch_ranks}")
    print(f"    RT subgradient-family combined rank: {rank_family}/{NE}; covered edges: {covered}/{NE}")
    assert rank_branch < NE
    assert rank_unit == NE
    assert rank_family == NE
    assert covered == NE

    rng = np.random.default_rng(151)
    h = rng.normal(size=NE)
    h /= np.linalg.norm(h)
    eps = 1.0e-8
    M1, S1 = mincut_matrix(regions, caps + eps * h)
    linear = M_branch @ h
    fd = (S1 - S0) / eps
    fd_err = float(np.max(np.abs(fd - linear)))
    changed = int(np.sum(np.any(M1 != M_branch, axis=1)))

    print("\n[1] RT first variation")
    print("    For a fixed RT surface gamma_A, delta S_A is the sum of edge-capacity")
    print("    perturbations on gamma_A: delta S_A = sum_{e in gamma_A} h_e.")
    print(f"    finite-difference max error = {fd_err:.2e}; changed cuts = {changed}/{len(regions)}")
    assert changed == 0
    assert fd_err < 1.0e-5

    h_family = rng.normal(size=NE)
    h_family /= np.linalg.norm(h_family)
    delta_s_family = M_family @ h_family
    h_hat, *_ = np.linalg.lstsq(M_family, delta_s_family, rcond=None)
    rec_err = float(np.linalg.norm(h_hat - h_family) / np.linalg.norm(h_family))
    print("\n[2] linearized entanglement tomography")
    print("    A single smooth RT branch does not see the whole edge space.  The")
    print("    symmetric link-state point is degenerate, however, and its RT subgradient")
    print("    family does span the whole edge space.")
    print(f"    reconstruction relative error from family delta S_A = M_family h: {rec_err:.2e}")
    assert rec_err < 1.0e-12

    B = toy_modular_kernel(regions)
    PB, rank_B, left_null = orth_projector_to_column_space(B)
    constraint = left_null.T @ M_unit
    rank_constraint = int(np.linalg.matrix_rank(constraint, tol=1.0e-10))

    delta_s_unit = M_unit @ h
    generic_residual = float(np.linalg.norm((np.eye(len(regions)) - PB) @ delta_s_unit) / np.linalg.norm(delta_s_unit))
    null_dim = NE - rank_constraint

    print("\n[3] finite first-law constraint operator")
    print("    Given a modular-energy kernel B_{A,i}, first law requires M h in range(B).")
    print("    Equivalently E h = 0 with E=(left-nullspace B)^T M. This is the graph")
    print("    analog of the all-ball linearized Einstein constraint.")
    print(f"    modular-kernel rank: {rank_B}/{NB} boundary directions")
    print(f"    independent edge constraints from naive first law: {rank_constraint}/{NE}")
    print(f"    first-law-compatible edge directions under this naive B: {null_dim}")
    print(f"    generic RT perturbation first-law residual: {generic_residual:.3f}")
    assert rank_B == NB
    assert rank_constraint == NE
    assert null_dim == 0
    assert generic_residual > 0.05

    print("""
[4] VERDICT
  PASSED: The verified substrate RT structure has the right linearized area
    calculus.  On a fixed nondegenerate branch, delta S_A is exactly the
    min-cut area variation.  At the symmetric link-state point, the RT
    subgradient family spans all edge-capacity perturbations.  There is no
    finite-QEC information obstruction to using RT entanglement as a linearized
    metric probe.

  NEW FINITE OBJECT: For any proposed modular-energy kernel B, the first law
    becomes the explicit graph equation E h = 0.  This is the finite analog of
    the Faulkner/Jacobson "all balls imply linearized Einstein" step.

  NEGATIVE CONTROL: A naive boundary-site parabolic kernel is rejected: it
    over-constrains the graph and leaves no nonzero first-law-compatible edge
    perturbation.  Thus RT alone is not enough, and a hand-waved boundary
    modular Hamiltonian does not close gravity.

  NOT CLOSED: The framework must derive the real QEC modular Hamiltonian /
    stress tensor from the service ledger, including its bulk/gravitational
    dressing, fix the 1/(4G) normalization intrinsically, and lift the disk
    cross-section to the HBC/dS horizon continuum.  Gravity is advanced from
    "horizon-input only" to a sharp intrinsic-gravity theorem target, not yet
    to a derived Einstein equation.
exit 0""")
    print("ALL ASSERTIONS PASSED -- RT first variation, full-rank tomography, and first-law constraint operator verified.")


if __name__ == "__main__":
    main()
