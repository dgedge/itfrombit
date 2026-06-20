#!/usr/bin/env python3
"""K04 numerical-residual checkout.

This is a fast companion to the heavier K04 instruments.  It checks the four
residuals left after the framed-QEC sign-Hamiltonian derivation:

  1. the numerical ratio w4/w6;
  2. the absolute bridge w6/Lambda;
  3. the Metropolis/gamma boot-cooling interpretation;
  4. the unembedded configuration-model K3,3 artefact sector.

Exit 0 means the arithmetic and graph-theoretic exclusions used in the canon
summary are internally consistent.  The w4/w6 ratio is now conditionally
promoted by k04_w4_w6_orbit_derivation.py; the failed historical 208 route is
kept separate.
"""

from __future__ import annotations

import math
from itertools import combinations


ALPHA0 = 1.0 / 137.035999084
GAMMA_PAPER = 0.995
T_START_OVER_W6 = 6.0
T_END_OVER_W6 = 0.5
LAMBDA_QCD_GEV = 0.332
DEBRIS_EXCESS_PER_VERTEX_W6 = 2.17


def metropolis_accept(delta_e: float, temperature: float) -> float:
    return min(1.0, math.exp(-delta_e / temperature))


def assert_metropolis_scale_invariance() -> None:
    for delta_e, temperature, scale in [
        (1.0, 6.0, 11.0),
        (4.0, 2.8, 0.37),
        (208.0, 3.15, 29.0),
    ]:
        base = metropolis_accept(delta_e, temperature)
        scaled = metropolis_accept(scale * delta_e, scale * temperature)
        assert abs(base - scaled) < 1e-15, (delta_e, temperature, scale)


def q3_edges() -> set[tuple[int, int]]:
    return {
        (u, v)
        for u in range(8)
        for v in range(u + 1, 8)
        if (u ^ v).bit_count() == 1
    }


def k33_edges() -> set[tuple[int, int]]:
    return {(i, 3 + j) for i in range(3) for j in range(3)}


def paths(adj: list[list[int]], start: int, goal: int, length: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    stack = [(start, (start,))]
    while stack:
        vertex, path = stack.pop()
        if len(path) == length:
            if goal in adj[vertex]:
                out.append(path)
            continue
        for nxt in adj[vertex]:
            if nxt not in path and nxt != goal:
                stack.append((nxt, path + (nxt,)))
    return out


def canonical_cycle(cycle: tuple[int, ...]) -> tuple[int, ...]:
    best: tuple[int, ...] | None = None
    for offset in range(len(cycle)):
        rot = cycle[offset:] + cycle[:offset]
        for candidate in (rot, tuple(reversed(rot))):
            if best is None or candidate < best:
                best = candidate
    assert best is not None
    return best


def count_cycles(edges: set[tuple[int, int]], n: int, length: int) -> int:
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    cycles: set[tuple[int, ...]] = set()
    for a, b in edges:
        for path in paths(adj, b, a, length - 1):
            cycles.add(canonical_cycle((a,) + path))
    return len(cycles)


def assert_k33_artifact_theorem() -> tuple[float, int]:
    q3 = q3_edges()
    k33 = k33_edges()
    q3_c4, q3_c6 = count_cycles(q3, 8, 4), count_cycles(q3, 8, 6)
    k33_c4, k33_c6 = count_cycles(k33, 6, 4), count_cycles(k33, 6, 6)
    assert (q3_c4, q3_c6) == (6, 16)
    assert (k33_c4, k33_c6) == (9, 6)

    # e_Q3/v = -(0.75 r + 2), e_K33/v = -(1.5 r + 1).
    # K33 beats Q3 iff 1.5 r + 1 > 0.75 r + 2.
    threshold = 4.0 / 3.0
    assert abs((1.5 * threshold + 1.0) - (0.75 * threshold + 2.0)) < 1e-15
    assert (1.5 * 1.5 + 1.0) > (0.75 * 1.5 + 2.0)

    # Z^3 nearest-neighbour graph exclusions:
    # - parity colouring makes it bipartite, hence triangle-free and K4-free;
    # - any two distinct vertices share at most two common nearest neighbours,
    #   while a same-side pair in K3,3 shares all three.
    units = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    max_common = 0
    for dx in range(-2, 3):
        for dy in range(-2, 3):
            for dz in range(-2, 3):
                if (dx, dy, dz) == (0, 0, 0):
                    continue
                common = sum(
                    1
                    for ux, uy, uz in units
                    if abs(ux - dx) + abs(uy - dy) + abs(uz - dz) == 1
                )
                max_common = max(max_common, common)
    assert max_common == 2

    # Direct finite sanity check: no triangle in a small Z^3 patch.
    patch = [
        (x, y, z)
        for x in range(-1, 2)
        for y in range(-1, 2)
        for z in range(-1, 2)
    ]
    patch_index = {p: i for i, p in enumerate(patch)}
    patch_edges: set[tuple[int, int]] = set()
    for a, b in combinations(patch, 2):
        if sum(abs(x - y) for x, y in zip(a, b)) == 1:
            patch_edges.add(tuple(sorted((patch_index[a], patch_index[b]))))
    for a, b, c in combinations(range(len(patch)), 3):
        tri = {
            tuple(sorted((a, b))),
            tuple(sorted((a, c))),
            tuple(sorted((b, c))),
        }
        assert not tri <= patch_edges

    return threshold, max_common


def main() -> None:
    assert_metropolis_scale_invariance()
    threshold, max_common = assert_k33_artifact_theorem()

    w6_over_lambda_ramp_start = 1.0 / T_START_OVER_W6
    w6_gev = w6_over_lambda_ramp_start * LAMBDA_QCD_GEV
    debris_scale = DEBRIS_EXCESS_PER_VERTEX_W6 * w6_over_lambda_ramp_start

    r_paper = math.log(T_START_OVER_W6 / T_END_OVER_W6) / -math.log(GAMMA_PAPER)
    gamma_alpha_bit = math.exp(-ALPHA0 * math.log(2.0))
    r_alpha_bit = math.log(T_START_OVER_W6 / T_END_OVER_W6) / (ALPHA0 * math.log(2.0))

    print("[0] K04 numerical-residual checkout")
    print("    Metropolis scale invariance: PASS")
    print("    cycle counters: Q3=(C4=6,C6=16), K3,3=(C4=9,C6=6): PASS")
    print()

    print("[1] w4/w6 ratio")
    print("    status: CONDITIONALLY DERIVED as w4/w6=2.")
    print("    source: Q3 loop-orbit edge-ledger normal ordering, not the old")
    print("    historical 208 fingerprint.")
    print("    negative check: the 208 fingerprint does not recover a consistent")
    print("    integer (w4,w6,lambda) under the cycle-optimal premise.")
    print("    coefficient audit: python_code/k04_w4_w6_orbit_derivation.py")
    print("    heavy audit: python_code/k04_weight_recovery.py")
    print()

    print("[2] w6/Lambda bridge")
    print("    status: CONDITIONALLY SELECTED, not K04-derived.")
    print("    premise: the boot bath/ramp-start temperature is Lambda.")
    print(f"    T_start = 6 w6 -> w6/Lambda = {w6_over_lambda_ramp_start:.6f}")
    print(f"    w6 = {w6_gev:.6f} GeV for Lambda = {LAMBDA_QCD_GEV:.3f} GeV")
    print(f"    debris scale ~= 2.17 w6 = {debris_scale:.3f} Lambda")
    print("    transition anchoring remains a different physical premise.")
    print()

    print("[3] Metropolis/gamma boot cooling")
    print("    status: GAMMA=0.995 IS A CANON-PINNED PROTOCOL PARAMETER.")
    print(f"    R_paper = ln(12)/[-ln(0.995)] = {r_paper:.3f} sweeps")
    print(f"    gamma_alpha_bit = exp(-alpha0 ln2) = {gamma_alpha_bit:.9f}")
    print(f"    R_alpha_bit = {r_alpha_bit:.3f} sweeps")
    print("    verdict: the numerical closeness is a theorem target, not a closure;")
    print("    the registered entropy-spike gamma gate failed at percent grade.")
    print()

    print("[4] unembedded configuration-model K3,3 sector")
    print("    status: CLOSED NEGATIVELY AS AN ARTIFACT.")
    print(f"    K3,3 beats Q3 in the unembedded toy iff w4/w6 > {threshold:.6f}.")
    print("    This explains why the old toy sweep range 1.5..2.5 was contaminated.")
    print("    Physical embedding excludes the phase: Z3 is bipartite (no K4),")
    print(f"    and max common neighbours of two Z3 sites is {max_common}, while K3,3 needs 3.")
    print("    heavy audit: python_code/debris_dark_matter_audit.py")
    print()

    print("[5] canon verdict")
    print("    Derived/conditional: K04 ordinal signs; w4/w6=2 under edge-ledger")
    print("    normal ordering; K3,3 artifact exclusion by embedding.")
    print("    Conditional: w6/Lambda ramp-start bridge.")
    print("    Open: a first-principles boot-cooling law.")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
