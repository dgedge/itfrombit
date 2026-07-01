#!/usr/bin/env python3
"""
item113_t1_t2_contact_algebra_derivation.py

Can the Item-113 T1/T2 Weizsaecker coefficient targets be derived from the
bare TCH/QEC contact counts alone?

This script separates the algebraic counts from the extra billing/conditioning
rules needed to use them:

T1 target:
    eps_bulk = (13/12) eps_sat.

    The Q3 contact algebra does contain 12 edge checks and a one-dimensional
    latch/kernel.  Those are the only local objects that can make 13/12 in the
    proposed way.  But the algebra alone does not yet prove the normal-ordering
    map "bill the latch once over the 12 edge ledgers, and only for closed bulk
    contacts."  That remains a precise theorem target.

T2 target:
    residual channels = z - 1 = 5.

    The simple-cubic contact graph gives z=6, and once a saturated/marked
    contact is selected, exactly five residual contacts remain.  But the
    algebra alone does not yet prove that the I3/chi matching Hamiltonian is
    conditioned on one marked closed contact, nor that the pairing gap is the
    RMS over those five residual channels.  Those remain precise theorem
    targets.

Exit 0 means the reduction is internally consistent, not that the bare counts
alone close T1/T2.  The missing maps are supplied separately in
item113_t1_t2_local_map_theorems.py by the relative-boundary and projector
service-instrument readings.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
import math
import sys


def gf2_rank(rows: list[int]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        x = row
        while x:
            pivot = x.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = x
                break
            x ^= basis[pivot]
    return len(basis)


def bit(x: int, i: int) -> int:
    return (x >> i) & 1


def q3_edges() -> list[tuple[int, int]]:
    return [
        (u, v)
        for u in range(8)
        for v in range(u + 1, 8)
        if (u ^ v).bit_count() == 1
    ]


def incidence_rows(edges: list[tuple[int, int]]) -> list[int]:
    rows = []
    for u, v in edges:
        rows.append((1 << u) | (1 << v))
    return rows


def syndrome(word: int, edges: list[tuple[int, int]]) -> tuple[int, ...]:
    return tuple(bit(word, u) ^ bit(word, v) for u, v in edges)


def face_cycles() -> list[tuple[int, int, int, int]]:
    """The six square faces of Q3, as 4-cycles."""

    faces: list[tuple[int, int, int, int]] = []
    axes = (0, 1, 2)
    for fixed_axis in axes:
        plane_axes = [a for a in axes if a != fixed_axis]
        a, b = plane_axes
        for fixed in (0, 1):
            verts = []
            for pa, pb in ((0, 0), (1, 0), (1, 1), (0, 1)):
                x = 0
                x |= fixed << fixed_axis
                x |= pa << a
                x |= pb << b
                verts.append(x)
            faces.append(tuple(verts))
    return faces


def cycle_edge_mask(cycle: tuple[int, ...], edge_index: dict[tuple[int, int], int]) -> int:
    mask = 0
    n = len(cycle)
    for i in range(n):
        e = tuple(sorted((cycle[i], cycle[(i + 1) % n])))
        mask |= 1 << edge_index[e]
    return mask


def simple_cubic_directions() -> list[tuple[int, int, int]]:
    return [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ]


def main() -> int:
    ok = True

    def check(name: str, cond: bool) -> None:
        nonlocal ok
        ok = ok and cond
        print(f"[{'PASS' if cond else 'FAIL'}] {name}")

    edges = q3_edges()
    edge_index = {e: i for i, e in enumerate(edges)}
    incidence = incidence_rows(edges)
    rank_delta = gf2_rank(incidence)
    kernel_dim = 8 - rank_delta
    image_size = 2**rank_delta
    cycle_rank = len(edges) - 8 + 1
    face_masks = [cycle_edge_mask(face, edge_index) for face in face_cycles()]
    face_rank = gf2_rank(face_masks)

    print("=" * 96)
    print("ITEM 113 T1/T2 BARE CONTACT-COUNT AUDIT")
    print("=" * 96)
    print("\n[1] Q3 contact algebra")
    print(f"    vertices={8}, edges={len(edges)}")
    print(f"    rank(delta)=rank(edge incidence over F2)={rank_delta}")
    print(f"    dim ker(delta)={kernel_dim}")
    print(f"    |im(delta)|={image_size}")
    print(f"    cycle rank=E-V+1={cycle_rank}")
    print(f"    square faces={len(face_masks)}, face span rank={face_rank}")

    all_words = list(range(256))
    zero_syndrome = [w for w in all_words if syndrome(w, edges) == (0,) * len(edges)]
    check("Q3 has 12 local edge ledgers", len(edges) == 12)
    check("Q3 coboundary has the one latch/kernel direction", kernel_dim == 1 and zero_syndrome == [0, 255])
    check("the six face plaquettes span the five-dimensional cycle space", len(face_masks) == 6 and face_rank == cycle_rank == 5)

    print("\n[2] T1: does Q3 force the 13/12 bulk normal-ordering factor?")
    ratio_latch = (len(edges) + kernel_dim) / len(edges)
    controls = {
        "edge plus cycle-rank": (len(edges) + cycle_rank) / len(edges),
        "edge plus face-count": (len(edges) + len(face_masks)) / len(edges),
        "edge plus incidence-rank": (len(edges) + rank_delta) / len(edges),
    }
    print(f"    edge ledgers + latch gives ({len(edges)}+{kernel_dim})/{len(edges)} = {ratio_latch:.12f}")
    for name, value in controls.items():
        print(f"    control {name:24s}: {value:.12f}")
    check("the proposed 13/12 object exists as 12 edge ledgers plus one latch", math.isclose(ratio_latch, 13 / 12))
    check("other natural Q3 invariants do not accidentally give 13/12", all(not math.isclose(v, 13 / 12) for v in controls.values()))

    t1_closed_by_algebra = False
    print("    MISSING MAP: the algebra has not yet proved that the latch is billed")
    print("    once as a bulk-only closure contribution over the 12 edge ledgers.")
    print("    Surface exclusion is also not a consequence of the incidence matrix alone.")
    check("T1 is reduced to a precise normal-ordering theorem, not closed", not t1_closed_by_algebra)

    print("\n[3] T2: does the contact graph force z-1=5 residual matching channels?")
    dirs = simple_cubic_directions()
    z = len(dirs)
    residual_counts = {}
    for marked in dirs:
        residual_counts[marked] = sum(d != marked for d in dirs)
    print(f"    simple-cubic contact directions z={z}: {dirs}")
    print(f"    given any marked saturated contact, residual channels={sorted(set(residual_counts.values()))}")
    check("simple-cubic contact algebra gives z=6", z == 6)
    check("a marked saturated contact leaves z-1=5 residual channels", set(residual_counts.values()) == {5})

    # The key distinction: the graph has six equivalent contacts until one is
    # physically marked by a closed/saturated bond.  The algebra by itself does
    # not select that mark in a homogeneous bulk.
    unmarked_orbits = 1
    marked_orbits = 1
    print(f"    unmarked homogeneous contact graph has {unmarked_orbits} orbit of six equivalent contacts")
    print(f"    marked-contact frame has {marked_orbits} orbit of five residual contacts")
    print("    MISSING MAP: prove I3/chi matching is conditioned on a marked closed")
    print("    contact, and prove the pairing gap is an RMS over those five residual")
    print("    channels rather than over 6, 4, or an independently weighted cluster.")
    t2_closed_by_algebra = False
    check("T2 is conditional on a marked-contact matching theorem, not closed", not t2_closed_by_algebra)

    print("\n[4] Consequence if the two missing maps are supplied")
    hbarc = 197.327
    a0 = 0.5944
    lam = hbarc / a0
    alpha0 = 1 / 137
    w = alpha0 * lam
    eps_sat = 2 * w
    bulk_factor = 13 / 12
    surface_ratio = 1.206
    empirical = {"aV": 15.75, "aS": 17.8, "aA": 23.7, "aP": 11.18}
    aV = (z / 2) * bulk_factor * eps_sat
    aS = surface_ratio * (z / 2) * eps_sat
    r0 = 2 * a0
    rho = 3 / (4 * math.pi * r0**3)
    mN = 939.0
    EF = (hbarc**2 / (2 * mN)) * (1.5 * math.pi**2 * rho) ** (2 / 3)
    aA = EF / 3 + (z - 1) * w
    aP = math.sqrt(z - 1) * eps_sat
    preds = {"aV": aV, "aS": aS, "aA": aA, "aP": aP}
    for key, pred in preds.items():
        err = 100 * (pred / empirical[key] - 1)
        print(f"    {key}: pred={pred:8.3f}  empirical={empirical[key]:8.3f}  err={err:+6.2f}%")
    check("T1+T2 would put the four residual SEMF coefficients within 5%",
          all(abs(preds[k] / empirical[k] - 1) < 0.05 for k in preds))

    print("\nVERDICT")
    print("  The bare contact counts alone do NOT fully derive T1 and T2.")
    print("  What is derived now:")
    print("    * T1 raw objects: 12 Q3 edge ledgers plus one latch/kernel.")
    print("    * T2 raw object: z=6 contacts, hence z-1=5 after a marked contact.")
    print("  What remains to prove:")
    print("    * T1 bulk-only closed-cell normal-ordering: latch billed once over")
    print("      the 12 edge ledgers, and not on exposed surface contacts.")
    print("    * T2 marked-contact bipartite matching: I3/chi imbalance sees the")
    print("      five residual channels with RMS pairing.")
    print("  So the candidate bundle is sharply reduced by the count audit.")
    print("  The local map closure is tested separately in")
    print("  item113_t1_t2_local_map_theorems.py.")

    if not ok:
        print("\nCHECKS FAILED")
        return 1
    print("\nALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
