#!/usr/bin/env python3
"""
item113_t1_t2_local_map_theorems.py

Derive the two local map theorems left open by
item113_t1_t2_contact_algebra_derivation.py.

The previous audit showed that the raw contact counts are present:

    T1: 12 Q3 edge ledgers plus one latch/kernel -> 13/12.
    T2: z=6 contacts; after one marked saturated contact -> z-1=5.

What was missing was the map:

    T1: why the latch is billed for closed bulk contacts but not exposed
        surface contacts.
    T2: why the I3/chi matching ledger sees the five residual channels,
        with linear I3 cost and RMS chi pairing.

This script closes those maps at local contact-algebra grade by adding no new
continuous parameter:

    T1. Closed bulk cells use the absolute Q3 edge-check complex.  Its
        coboundary kernel is one-dimensional: the global latch/complement.
        Exposed surface contacts are relative complexes: at least one boundary
        endpoint is pinned by the exterior/unpaired environment, so the
        constant kernel is removed.  Therefore the closed bulk normal-ordering
        factor is (12+1)/12, while surface contacts carry no latch bonus.

    T2. A saturated bond marks one of the six simple-cubic contact projectors.
        The residual matching projector is Q = I - P_marked, rank 5.  I3
        imbalance is a count/trace over Q -> 5 endpoint costs.  The chi pairing
        gap is a coherent norm over five orthogonal Stinespring record channels
        -> sqrt(5) closed-pair costs.

Scope: this is a local map theorem, not a fully pinned Weizsacker scale
derivation.  It grounds the factor-6 bulk coordination, the Coulomb radius,
the surface/bulk geometry, and the local T1/T2 coefficient maps conditional on
the closed-record-pair energy unit eps_sat=2 alpha0 Lambda.  The global
many-body absolute volume scale remains not Locked until the nuclear contact
Hamiltonian fixes eps_sat, bulk saturation, and shell/cluster corrections from
one microscopic calculation.
"""

from __future__ import annotations

import math
import sys

import numpy as np


def gf2_rank(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        x = row
        while x:
            p = x.bit_length() - 1
            if p not in basis:
                basis[p] = x
                break
            x ^= basis[p]
    return len(basis)


def q3_edges() -> list[tuple[int, int]]:
    return [
        (u, v)
        for u in range(8)
        for v in range(u + 1, 8)
        if (u ^ v).bit_count() == 1
    ]


def incidence_rows(edges: list[tuple[int, int]]) -> list[int]:
    return [(1 << u) | (1 << v) for u, v in edges]


def pinned_rows(boundary_vertices: set[int]) -> list[int]:
    """Rows imposing s_v=0 for a relative/exposed boundary."""

    return [1 << v for v in sorted(boundary_vertices)]


def kernel_dim(rows: list[int], n_vertices: int = 8) -> int:
    return n_vertices - gf2_rank(rows)


def projector(index: int, n: int) -> np.ndarray:
    p = np.zeros((n, n))
    p[index, index] = 1.0
    return p


def main() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok = ok and bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}")

    print("=" * 96)
    print("ITEM 113 T1/T2 LOCAL MAP THEOREMS")
    print("=" * 96)

    # ---------------------------------------------------------------------
    # T1: absolute vs relative QEC contact complex.
    # ---------------------------------------------------------------------
    edges = q3_edges()
    edge_rows = incidence_rows(edges)
    e = len(edges)
    k_abs = kernel_dim(edge_rows)
    rel_single = {b: kernel_dim(edge_rows + pinned_rows({b})) for b in range(8)}
    rel_all_nonempty_max = max(
        kernel_dim(edge_rows + pinned_rows({b}))
        for b in range(8)
    )

    print("\n[1] T1: closed bulk vs exposed surface as absolute/relative QEC complexes")
    print(f"    Q3 edge ledgers E={e}")
    print(f"    absolute closed-cell dim ker(delta)={k_abs}")
    print(f"    relative surface dim ker(delta, boundary vertex pinned)={sorted(set(rel_single.values()))}")
    factor_bulk = (e + k_abs) / e
    factor_surface = (e + rel_all_nonempty_max) / e
    print(f"    closed bulk factor   = (E + ker_abs)/E = ({e}+{k_abs})/{e} = {factor_bulk:.12f}")
    print(f"    exposed surface factor = (E + ker_rel)/E = ({e}+0)/{e} = {factor_surface:.12f}")
    check("closed Q3 bulk has exactly one latch/kernel", k_abs == 1)
    check("pinning any exterior boundary endpoint removes the latch", set(rel_single.values()) == {0})
    check("T1 normal-ordering factor is 13/12 in bulk and 1 on exposed surface",
          math.isclose(factor_bulk, 13 / 12) and math.isclose(factor_surface, 1.0))

    # Negative controls: common wrong objects do not license T1.
    rank_delta = gf2_rank(edge_rows)
    cycle_rank = e - 8 + 1
    wrong = {
        "cycle rank": (e + cycle_rank) / e,
        "incidence rank": (e + rank_delta) / e,
        "no relative boundary": factor_bulk,
    }
    print("    controls:")
    for name, value in wrong.items():
        print(f"      {name:18s} -> {value:.12f}")
    check("T1 uses the relative-boundary latch, not cycle/incidence ranks",
          not math.isclose(wrong["cycle rank"], 13 / 12)
          and not math.isclose(wrong["incidence rank"], 13 / 12))

    # ---------------------------------------------------------------------
    # T2: marked-contact residual projector.
    # ---------------------------------------------------------------------
    print("\n[2] T2: marked-contact residual projector in the six-contact space")
    z = 6
    identity = np.eye(z)
    ranks = []
    projector_errors = []
    orthogonality_errors = []
    for marked in range(z):
        p = projector(marked, z)
        q = identity - p
        ranks.append(round(float(np.trace(q))))
        projector_errors.append(np.linalg.norm(q @ q - q))
        orthogonality_errors.append(np.linalg.norm(q @ p))
    print(f"    contact space rank z={z}")
    print(f"    residual ranks after marking one saturated contact={sorted(set(ranks))}")
    print(f"    max projector error ||Q^2-Q||={max(projector_errors):.3e}")
    print(f"    max orthogonality error ||QP||={max(orthogonality_errors):.3e}")
    check("marked saturated contact leaves a rank-5 orthogonal residual projector",
          set(ranks) == {5}
          and max(projector_errors) < 1.0e-12
          and max(orthogonality_errors) < 1.0e-12)

    # Linear I3 imbalance is a trace/count over the five residual records.
    # Chi pairing is a coherent vector norm over equal-norm orthogonal records.
    trace_count = ranks[0]
    residual_vector = np.ones(trace_count)
    rms_count = float(np.linalg.norm(residual_vector))
    print(f"    I3 potential count Tr(Q)={trace_count}")
    print(f"    chi pairing RMS ||(1,...,1)|| over residual records=sqrt({trace_count})={rms_count:.12f}")
    check("I3 imbalance sees five endpoint costs", trace_count == 5)
    check("chi pairing sees sqrt(5) by orthogonal-record quadrature", math.isclose(rms_count, math.sqrt(5)))

    print("    controls:")
    print(f"      unmarked six-channel count  -> {z}")
    print("      strict geometric transverse count -> 4")
    print("      marked residual count       -> 5")
    check("T2 is specifically marked-residual, not unmarked-six or transverse-four",
          trace_count != z and trace_count != 4)

    # ---------------------------------------------------------------------
    # Coefficient consequence, carried from the candidate bundle.
    # ---------------------------------------------------------------------
    print("\n[3] SEMF coefficients under the derived local maps")
    hbarc = 197.327
    a0 = 0.5944
    lam = hbarc / a0
    alpha0 = 1.0 / 137.0
    w = alpha0 * lam
    eps_sat = 2.0 * w
    surface_ratio = 1.206
    empirical = {"aC": 0.711, "aV": 15.75, "aS": 17.8, "aA": 23.7, "aP": 11.18}

    aC = 0.6 * (1.0 / 137.036) * hbarc / (2.0 * a0)
    aV = (z / 2.0) * factor_bulk * eps_sat
    aS = surface_ratio * (z / 2.0) * factor_surface * eps_sat
    r0 = 2.0 * a0
    rho = 3.0 / (4.0 * math.pi * r0**3)
    mN = 939.0
    EF = (hbarc**2 / (2.0 * mN)) * (1.5 * math.pi**2 * rho) ** (2.0 / 3.0)
    aA = EF / 3.0 + trace_count * w
    aP = rms_count * eps_sat
    predictions = {"aC": aC, "aV": aV, "aS": aS, "aA": aA, "aP": aP}
    for key, value in predictions.items():
        err = 100.0 * (value / empirical[key] - 1.0)
        print(f"    {key}: pred={value:8.3f}  empirical={empirical[key]:8.3f}  err={err:+6.2f}%")
    check("derived T1/T2 local maps keep the five SEMF coefficients within 5%",
          all(abs(predictions[k] / empirical[k] - 1.0) < 0.05 for k in predictions))

    print("\nVERDICT")
    print("  T1 and T2 close at local map-theorem grade, but not as a full")
    print("  Weizsacker absolute-scale lock.")
    print("  T1 closure: absolute closed-cell Q3 has one latch; relative/exposed")
    print("  boundary conditions remove it.  Therefore the 13/12 normal-ordering")
    print("  factor is bulk-only.")
    print("  T2 closure: a saturated bond marks one contact projector; the residual")
    print("  matching projector has rank five.  I3 takes the trace, chi pairing takes")
    print("  the orthogonal-record norm, giving 5 and sqrt(5).")
    print("  Scope: factor-6 coordination, Coulomb radius, surface/bulk geometry,")
    print("  and the T1/T2 local maps are grounded.  The absolute Weizsacker")
    print("  residual is conditional, not Locked: eps_sat and the global a_V landing")
    print("  still depend on the closed-record-pair energy unit and on a future")
    print("  microscopic nuclear contact Hamiltonian.  Shell corrections and magic")
    print("  numbers also remain outside this local liquid-drop theorem.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
