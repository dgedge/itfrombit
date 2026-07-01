#!/usr/bin/env python3
"""Bridge-local latch audit for the dressed-alpha traffic-clock route.

Question:
    Can the Peierls-current bridge readout be promoted to a two-link workload
    over a 12-edge Q3 strain sweep plus one bridge-local latch?

The current canon has one clean one-bit latch theorem: the Q3 coboundary
delta: F2^8 -> F2^12 has kernel {0, ALL}.  The missing bit is the
global-complement/vacuum latch used in V_cell and horizon ledgers.

This audit asks whether the Peierls bridge readout has an analogous local
one-bit latch.  It does not test whether a new primitive could be invented; it
tests whether the existing bridge/current algebra already contains the latch
needed to license the 2/13 traffic-clock assignment.
"""
from __future__ import annotations

from itertools import product

import numpy as np

import bridge_web_lindblad_keldysh_poles as poles
import dressed_alpha_bridge_web_open_system as bw


T_M = 1.0 / 3.0


def gf2_rank(rows: list[int], n_cols: int) -> int:
    """Rank of a binary matrix represented by integer row masks."""
    pivots = list(rows)
    rank = 0
    for col in reversed(range(n_cols)):
        pivot_at = next((i for i in range(rank, len(pivots)) if (pivots[i] >> col) & 1), None)
        if pivot_at is None:
            continue
        pivots[rank], pivots[pivot_at] = pivots[pivot_at], pivots[rank]
        for i in range(len(pivots)):
            if i != rank and ((pivots[i] >> col) & 1):
                pivots[i] ^= pivots[rank]
        rank += 1
    return rank


def kernel_words(rows: list[int], n_cols: int) -> list[int]:
    """Brute-force kernel words for small GF(2) maps."""
    out = []
    for word in range(1 << n_cols):
        if all(((row & word).bit_count() % 2) == 0 for row in rows):
            out.append(word)
    return out


def q3_edges() -> list[tuple[int, int]]:
    vertices = list(product((0, 1), repeat=3))
    index = {v: i for i, v in enumerate(vertices)}
    edges: list[tuple[int, int]] = []
    for i, vertex in enumerate(vertices):
        for axis in range(3):
            other = list(vertex)
            other[axis] ^= 1
            j = index[tuple(other)]
            if i < j:
                edges.append((i, j))
    return edges


def q3_latch_reference() -> tuple[int, int, list[int], list[tuple[int, int]]]:
    edges = q3_edges()
    rows = [(1 << i) | (1 << j) for i, j in edges]
    rank = gf2_rank(rows, 8)
    ker = kernel_words(rows, 8)
    return rank, 8 - rank, ker, edges


def peierls_endpoint_parity() -> tuple[int, int, int, int]:
    """Rank/kernel of the two Peierls bridge-link parity rows.

    The four endpoint map is the most generous bridge-local version:
    columns are (0, 1, 8, 9), rows are links (0,8) and (1,9).
    The sixteen-vertex map embeds the same two rows in the full bridge.
    """
    endpoint_rows = [
        (1 << 0) | (1 << 2),  # endpoint labels 0 and 8
        (1 << 1) | (1 << 3),  # endpoint labels 1 and 9
    ]
    full_rows = [
        (1 << 0) | (1 << 8),
        (1 << 1) | (1 << 9),
    ]
    endpoint_rank = gf2_rank(endpoint_rows, 4)
    full_rank = gf2_rank(full_rows, 16)
    return endpoint_rank, 4 - endpoint_rank, full_rank, 16 - full_rank


def pair_current_diagnostics(
    h2: np.ndarray,
    idx: dict[tuple[int, int], int],
    s_eta,
) -> dict[str, float | int]:
    j_minus = bw.current_portal(idx, -1)
    j_plus = bw.current_portal(idx, +1)
    current_projector = np.outer(j_minus, j_minus)
    current_rank = int(np.linalg.matrix_rank(current_projector, tol=1e-12))
    rows, rates, _weights, web_dark_dim = poles.grouped_web_poles(h2, idx, s_eta)
    positive = rates[rates > poles.RATE_TOL]
    strongest = float(positive.max()) if len(positive) else 0.0
    return {
        "pair_dim": len(idx),
        "current_rank": current_rank,
        "current_dark_dim": len(idx) - current_rank,
        "j_minus_dot_j_plus": float(j_minus @ j_plus),
        "j_minus_h_j_minus": float(j_minus @ h2 @ j_minus),
        "j_plus_h_j_plus": float(j_plus @ h2 @ j_plus),
        "web_bright_eigenspaces": int(np.sum(rates > poles.RATE_TOL)),
        "web_dark_dim": web_dark_dim,
        "web_strongest_pole": strongest,
        "web_rows": len(rows),
    }


def traffic_target(h2: np.ndarray, idx: dict[tuple[int, int], int], s_eta):
    c_j, gap_unit = bw.response_coefficient(h2, idx)
    gamma_esc, mean_omega, _weight = bw.gamma_escape(
        h2,
        idx,
        s_eta,
        t_m=T_M,
        eta=bw.ETA_PIN,
        g_portal=bw.G_PORTAL,
    )
    target_gap = c_j * gamma_esc / bw.DELTA_TARGET
    gamma_target = gap_unit * T_M * T_M / target_gap
    return c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target


def main() -> None:
    q3_rank, q3_kernel_dim, q3_ker, q3_edge_list = q3_latch_reference()
    endpoint_rank, endpoint_kernel_dim, full_rank, full_kernel_dim = peierls_endpoint_parity()

    h2, _pairs, idx, _bas = bw.build_pair_system()
    s_eta, _omega_max = bw.photon_form_factor()
    current = pair_current_diagnostics(h2, idx, s_eta)
    c_j, gap_unit, gamma_esc, mean_omega, target_gap, gamma_target = traffic_target(h2, idx, s_eta)

    two_over_12 = 2.0 / 12.0
    two_over_13 = 2.0 / 13.0
    exact_two_link_slots = 2.0 / gamma_target
    fractional_latch_weight = exact_two_link_slots - 12.0

    print("DRESSED-ALPHA BRIDGE-LOCAL LATCH AUDIT")
    print("=" * 92)
    print("[1] Reference latch: Q3 coboundary delta")
    print(f"  Q3 strain edges              = {len(q3_edge_list)}")
    print(f"  rank(delta)                  = {q3_rank}")
    print(f"  kernel dimension             = {q3_kernel_dim}")
    print(f"  kernel words                 = {[format(k, '08b') for k in q3_ker]}")
    print("  This is the real one-bit latch: the global complement/intercept a0.")

    print("\n[2] Peierls bridge endpoint parity")
    print("  rows                         = link parities (0,8), (1,9)")
    print(f"  endpoint-local rank          = {endpoint_rank} of 4")
    print(f"  endpoint-local kernel dim    = {endpoint_kernel_dim}")
    print(f"  full 16-vertex rank          = {full_rank} of 16")
    print(f"  full 16-vertex kernel dim    = {full_kernel_dim}")
    print("  This is not a Q3-style single blind latch; even the endpoint-restricted")
    print("  bridge map leaves a two-bit kernel.")

    print("\n[3] Peierls current and integrated-web poles")
    print(f"  pair Hilbert dimension       = {current['pair_dim']}")
    print(f"  raw current-projector rank   = {current['current_rank']}")
    print(f"  raw current dark dimension   = {current['current_dark_dim']}")
    print(f"  <J_-|J_+>                    = {current['j_minus_dot_j_plus']:.3e}")
    print(f"  <J_-|H_pair|J_->             = {current['j_minus_h_j_minus']:.3e}")
    print(f"  <J_+|H_pair|J_+>             = {current['j_plus_h_j_plus']:.3e}")
    print(f"  web bright eigenspaces       = {current['web_bright_eigenspaces']}")
    print(f"  web dark bridge dimension    = {current['web_dark_dim']} of {current['pair_dim']}")
    print(f"  strongest web pole           = {current['web_strongest_pole']:.6e} Lambda")
    print("  The eta=-1 branch is a Wilson/Peierls phase choice already consumed by")
    print("  the current vertex and form factor.  It is not an extra service slot.")

    print("\n[4] Traffic-clock numerics")
    print(f"  Gamma_esc                    = {gamma_esc:.6e} Lambda")
    print(f"  <omega>_emit                 = {mean_omega:.6f} Lambda")
    print(f"  response c_J                 = {c_j:.6f}")
    print(f"  unit bridge gap              = {gap_unit:.6f} Lambda")
    print(f"  required response gap        = {target_gap:.6f} Lambda")
    print(f"  required gamma_mon           = {gamma_target:.6f} Lambda")
    print(f"  two links / 12               = {two_over_12:.6f} ({two_over_12 / gamma_target:.6f} x target)")
    print(f"  two links / 13               = {two_over_13:.6f} ({two_over_13 / gamma_target:.6f} x target)")
    print(f"  exact two-link slot count    = {exact_two_link_slots:.6f}")
    print(f"  implied latch weight over 12 = {fractional_latch_weight:.6f}")

    print("\nVERDICT")
    print("  Refuted as a current theorem.  The existing substrate has a genuine")
    print("  one-bit Q3/horizon latch, but the Peierls bridge readout does not inherit")
    print("  it.  The two-link current is low-rank, leaves large dark subspaces, and")
    print("  its eta branch is a static flux phase rather than a monitored latch.")
    print("  Therefore 2/13 remains a near-hit assignment, not a derived bridge-local")
    print("  traffic clock.  A new bridge-local service primitive could be postulated,")
    print("  but the present bridge-plus-web algebra does not prove it.")

    assert len(q3_edge_list) == 12
    assert q3_rank == 7
    assert q3_kernel_dim == 1
    assert q3_ker == [0, 0xFF]
    assert endpoint_rank == 2
    assert endpoint_kernel_dim == 2
    assert full_rank == 2
    assert full_kernel_dim == 14
    assert current["current_rank"] == 1
    assert current["current_dark_dim"] == current["pair_dim"] - 1
    assert current["web_dark_dim"] > 0
    assert abs(current["j_minus_dot_j_plus"]) < 1e-12
    assert abs(current["j_minus_h_j_minus"]) < 1e-12
    assert abs(current["j_plus_h_j_plus"]) < 1e-12
    assert 12.0 < exact_two_link_slots < 13.0
    assert abs((two_over_13 / gamma_target) - 1.0) > 0.01
    assert fractional_latch_weight < 1.0
    print("ALL ASSERTS PASSED")


if __name__ == "__main__":
    main()
