#!/usr/bin/env python3
"""Derive the K04 weight ratio from Q3 loop-orbit normal ordering.

The sign-Hamiltonian audit derives only

    E = -w4*C4 - w6*C6 + lambda*sum_v(deg(v)-3)^2,  w4,w6,lambda>0.

This script tests the natural next coefficient principle:

  * closed QEC/service loops are billed from local edge records;
  * a single edge record must not be counted more than once inside the same
    symmetry orbit of loop witnesses;
  * therefore each loop in an orbit is weighted by the inverse of the number
    of loops from that orbit containing a fixed edge;
  * the published K04 Hamiltonian has one C6 coefficient, so the two Q3
    six-cycle orbits are coarse-grained to one effective C6 coefficient by
    their microcanonical loop-count average.

The result is

    w4/w6 = 2.

Scope: this is a conditional derivation of the numerical ratio under the
edge-ledger normal-ordering rule.  The script also prints the exact point at
which the rule enters, so the result does not masquerade as pure covariance.
"""

from __future__ import annotations

from itertools import combinations, permutations, product
from fractions import Fraction


Point = int
Edge = tuple[int, int]
Cycle = tuple[int, ...]
Perm = tuple[int, ...]


def canonical_cycle(cycle: Cycle) -> Cycle:
    candidates: list[Cycle] = []
    n = len(cycle)
    for i in range(n):
        rot = cycle[i:] + cycle[:i]
        candidates.append(rot)
        candidates.append(tuple(reversed(rot)))
    return min(candidates)


def q3_edges() -> set[Edge]:
    return {
        (u, v)
        for u in range(8)
        for v in range(u + 1, 8)
        if (u ^ v).bit_count() == 1
    }


def q3_adjacency() -> list[set[int]]:
    adj = [set() for _ in range(8)]
    for a, b in q3_edges():
        adj[a].add(b)
        adj[b].add(a)
    return adj


def simple_cycles(k: int) -> list[Cycle]:
    adj = q3_adjacency()
    out: set[Cycle] = set()
    for subset in combinations(range(8), k):
        first, rest = subset[0], subset[1:]
        for perm in permutations(rest):
            if perm[0] > perm[-1]:
                continue
            cyc = (first,) + perm
            if all(cyc[(i + 1) % k] in adj[cyc[i]] for i in range(k)):
                out.add(canonical_cycle(cyc))
    return sorted(out)


def automorphisms_q3() -> list[Perm]:
    """Aut(Q3) = translations F2^3 semidirect coordinate permutations."""

    autos: list[Perm] = []
    for shift in range(8):
        for axes in permutations((0, 1, 2)):
            image = []
            for p in range(8):
                q = 0
                for new_axis, old_axis in enumerate(axes):
                    if p & (1 << old_axis):
                        q |= 1 << new_axis
                image.append(q ^ shift)
            autos.append(tuple(image))
    assert len(set(autos)) == 48
    return sorted(set(autos))


def apply_perm(cycle: Cycle, perm: Perm) -> Cycle:
    return canonical_cycle(tuple(perm[p] for p in cycle))


def cycle_orbits(cycles: list[Cycle]) -> list[list[Cycle]]:
    autos = automorphisms_q3()
    remaining = set(cycles)
    orbits: list[list[Cycle]] = []
    while remaining:
        seed = next(iter(remaining))
        orbit = {apply_perm(seed, g) for g in autos}
        assert orbit <= set(cycles)
        orbits.append(sorted(orbit))
        remaining -= orbit
    return sorted(orbits, key=lambda xs: (len(xs), xs[0]))


def cycle_edges(cycle: Cycle) -> set[Edge]:
    return {
        tuple(sorted((cycle[i], cycle[(i + 1) % len(cycle)])))
        for i in range(len(cycle))
    }


def incidence_per_edge(orbit: list[Cycle]) -> int:
    edges = q3_edges()
    counts = []
    for edge in edges:
        counts.append(sum(edge in cycle_edges(cycle) for cycle in orbit))
    assert len(set(counts)) == 1, counts
    return counts[0]


def complement_distance(cycle: Cycle) -> int | None:
    complement = sorted(set(range(8)) - set(cycle))
    if len(complement) != 2:
        return None
    return (complement[0] ^ complement[1]).bit_count()


def orbit_summary(k: int) -> list[dict[str, object]]:
    rows = []
    for orbit in cycle_orbits(simple_cycles(k)):
        inc = incidence_per_edge(orbit)
        weight = Fraction(1, inc)
        distances = sorted({complement_distance(cycle) for cycle in orbit})
        rows.append(
            {
                "length": k,
                "size": len(orbit),
                "incidence": inc,
                "weight": weight,
                "complement_distances": distances,
            }
        )
    return rows


def main() -> None:
    print("[0] Q3 loop-orbit normal-ordering audit")
    c4 = simple_cycles(4)
    c6 = simple_cycles(6)
    assert len(c4) == 6
    assert len(c6) == 16
    print(f"    C4(Q3)={len(c4)}, C6(Q3)={len(c6)}")
    print(f"    |Aut(Q3)|={len(automorphisms_q3())}")
    print()

    print("[1] Symmetry orbits and edge overcount")
    rows4 = orbit_summary(4)
    rows6 = orbit_summary(6)
    for row in rows4 + rows6:
        dist = row["complement_distances"]
        dist_label = "-" if dist == [None] else ",".join(str(x) for x in dist)
        print(
            "    L={length}: orbit size={size:2d}, loops/edge={incidence}, "
            "normal-ordered weight={weight}, complement distance={dist}".format(
                length=row["length"],
                size=row["size"],
                incidence=row["incidence"],
                weight=row["weight"],
                dist=dist_label,
            )
        )

    assert len(rows4) == 1
    assert [row["size"] for row in rows6] == [4, 12]
    assert sorted(row["incidence"] for row in rows6) == [2, 6]
    print()

    print("[2] Effective K04 coefficients")
    w4 = rows4[0]["weight"]
    assert isinstance(w4, Fraction)
    # K04 uses one C6 count.  Coarse-grain the two C6 orbit weights by the
    # number of loop witnesses in each orbit.
    c6_total = sum(int(row["size"]) for row in rows6)
    w6 = sum(int(row["size"]) * row["weight"] for row in rows6) / c6_total
    assert isinstance(w6, Fraction)
    ratio = w4 / w6
    print(f"    w4 = {w4} from the single C4 orbit")
    print(f"    w6 = (sum_orbits |O|/incidence_O) / 16 = {w6}")
    print(f"    w4/w6 = {ratio}")
    assert ratio == 2
    print()

    print("[3] What would not be derived by covariance alone")
    print("    Pure covariance allows one C4 orbit and two C6 orbit coefficients.")
    print("    The coefficient 2 enters only after adding edge-ledger normal ordering")
    print("    and then collapsing the two six-cycle orbits to the single K04 C6 term.")
    print()

    print("[4] Verdict")
    print("    Conditional closure: w4/w6 = 2 follows from Q3 loop-orbit")
    print("    normal-ordering of the local service/current edge ledger.")
    print("    Remaining non-K04 bridges: w6/Lambda and physical boot cooling gamma.")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
