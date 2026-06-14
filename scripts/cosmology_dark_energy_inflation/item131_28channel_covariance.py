#!/usr/bin/env python3
r"""ITEM 131 / 56: 28-channel covariance and reset-instrument audit.

This script tests the proposed route for deriving the uniform serial
28-channel clock rather than postulating it.

Representation:
    The 14 weight-4 words of the [8,4,4] code are represented as affine
    hyperplanes in F_2^3:

        H(a,b) = {x in F_2^3 : a.x = b},  a != 0, b in {0,1}.

    The automorphism group is AGL(3,2), order 8*168 = 1344.  It acts
    transitively on the 14 hyperplanes.  Tensoring with the C2 exchange of
    the two transverse modes gives one orbit on 28 channels.

What closes:
    Under the one-jump reset-instrument constraint, a scheduler is just a
    non-negative weight vector p_x on the 28 reset outcomes.  Covariance
    under the transitive 28-channel group forces p_x to be constant, and
    normalisation then gives p_x = 1/28.

What does not close:
    The full covariant 28x28 operator commutant is not one-dimensional.
    It has six orbital classes, so covariance alone still permits many
    non-serial covariant operators.  The one-jump/finite-bandwidth
    instrument constraint is therefore load-bearing.  The companion script
    item131_one_jump_premise.py formalises that scheduler layer, and
    item131_w_to_28_instrument.py constructs the W=S*C -> 2x14
    service-instrument lift.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import product


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def parity(x: int) -> int:
    return x.bit_count() & 1


def dot(a: int, x: int) -> int:
    return parity(a & x)


def rank_f2(vectors: list[int]) -> int:
    """Rank of 3-bit vectors over F_2."""
    basis = [0, 0, 0]
    rank = 0
    for raw in vectors:
        v = raw
        while v:
            pivot = v.bit_length() - 1
            if basis[pivot]:
                v ^= basis[pivot]
            else:
                basis[pivot] = v
                rank += 1
                break
    return rank


def linear_maps_gl3() -> list[tuple[int, ...]]:
    """All GL(3,2) maps as point permutations on {0,...,7}."""
    maps: list[tuple[int, ...]] = []
    nonzero = list(range(1, 8))
    for c0, c1, c2 in product(nonzero, repeat=3):
        if rank_f2([c0, c1, c2]) != 3:
            continue
        image = []
        for x in range(8):
            y = 0
            if x & 1:
                y ^= c0
            if x & 2:
                y ^= c1
            if x & 4:
                y ^= c2
            image.append(y)
        maps.append(tuple(image))
    return maps


def affine_maps_agl32() -> list[tuple[int, ...]]:
    maps = []
    for lin in linear_maps_gl3():
        for t in range(8):
            maps.append(tuple(y ^ t for y in lin))
    return maps


@dataclass(frozen=True)
class Hyperplane:
    normal: int
    offset: int
    support: int


def hyperplanes() -> list[Hyperplane]:
    planes = []
    seen: set[int] = set()
    for normal in range(1, 8):
        for offset in [0, 1]:
            support = 0
            for x in range(8):
                if dot(normal, x) == offset:
                    support |= 1 << x
            planes.append(Hyperplane(normal, offset, support))
            seen.add(support)
    check(len(planes) == 14, "there are 14 affine hyperplanes H(a,b), a != 0")
    check(len(seen) == 14, "the 14 hyperplanes have distinct weight-4 supports")
    check(all(p.support.bit_count() == 4 for p in planes), "every hyperplane support has weight 4")
    return planes


def permute_support(point_perm: tuple[int, ...], support: int) -> int:
    out = 0
    for x in range(8):
        if support & (1 << x):
            out |= 1 << point_perm[x]
    return out


def induced_hyperplane_perms(planes: list[Hyperplane]) -> list[tuple[int, ...]]:
    support_to_index = {p.support: i for i, p in enumerate(planes)}
    perms = []
    for point_perm in affine_maps_agl32():
        image = []
        for plane in planes:
            image.append(support_to_index[permute_support(point_perm, plane.support)])
        perms.append(tuple(image))
    # AGL(3,2) acts faithfully on the 14 affine hyperplanes.
    return sorted(set(perms))


def channel_perms(hyperplane_perms: list[tuple[int, ...]], include_mode_flip: bool) -> list[tuple[int, ...]]:
    """Permutation action on channel labels (mode, hyperplane), encoded as 14*mode+h."""
    perms = []
    mode_maps = [(0, 1), (1, 0)] if include_mode_flip else [(0, 1)]
    for hp in hyperplane_perms:
        for mm in mode_maps:
            image = [0] * 28
            for mode in [0, 1]:
                for h in range(14):
                    image[14 * mode + h] = 14 * mm[mode] + hp[h]
            perms.append(tuple(image))
    return sorted(set(perms))


class DSU:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def groups(self) -> dict[int, list[int]]:
        out: dict[int, list[int]] = defaultdict(list)
        for i in range(len(self.parent)):
            out[self.find(i)].append(i)
        return out


def label_orbits(perms: list[tuple[int, ...]], n: int) -> list[list[int]]:
    dsu = DSU(n)
    for perm in perms:
        for i, pi in enumerate(perm):
            dsu.union(i, pi)
    return list(dsu.groups().values())


def ordered_pair_orbits(perms: list[tuple[int, ...]], n: int) -> list[list[tuple[int, int]]]:
    dsu = DSU(n * n)
    for perm in perms:
        for i in range(n):
            base = n * i
            pbase = n * perm[i]
            for j in range(n):
                dsu.union(base + j, pbase + perm[j])
    groups = []
    for members in dsu.groups().values():
        groups.append([(m // n, m % n) for m in members])
    return groups


def relation_14(i: int, j: int, planes: list[Hyperplane]) -> str:
    if i == j:
        return "same hyperplane"
    overlap = (planes[i].support & planes[j].support).bit_count()
    if overlap == 0:
        return "complement/parallel"
    if overlap == 2:
        return "intersecting"
    raise AssertionError(f"unexpected hyperplane overlap {overlap}")


def relation_28(i: int, j: int, planes: list[Hyperplane]) -> str:
    mi, hi = divmod(i, 14)
    mj, hj = divmod(j, 14)
    mode = "same mode" if mi == mj else "opposite mode"
    return f"{mode}; {relation_14(hi, hj, planes)}"


def orbit_relation_summary(orbits: list[list[tuple[int, int]]], relation_fn) -> Counter:
    summary: Counter[str] = Counter()
    for orbit in orbits:
        labels = {relation_fn(i, j) for i, j in orbit}
        check(len(labels) == 1, f"ordered-pair orbit has a single relation class: {labels}")
        summary[next(iter(labels))] += 1
    return summary


def invariant_weight_dimension(perms: list[tuple[int, ...]], n: int) -> int:
    """Dimension of p with p_i = p_{g i}; equals number of label orbits."""
    return len(label_orbits(perms, n))


def print_orbits(name: str, orbits: list[list[int]]) -> None:
    sizes = sorted(len(o) for o in orbits)
    print(f"  {name}: {len(orbits)} orbit(s), sizes {sizes}")


def main() -> None:
    print("ITEM 131 / 56: 28-CHANNEL COVARIANCE AUDIT")

    print("\n[1] Build [8,4,4] weight-4 hyperplane representation")
    planes = hyperplanes()
    gl = linear_maps_gl3()
    agl = affine_maps_agl32()
    hp_perms = induced_hyperplane_perms(planes)
    check(len(gl) == 168, "GL(3,2) has order 168")
    check(len(agl) == 1344, "AGL(3,2) has order 1344")
    check(len(hp_perms) == 1344, "AGL(3,2) acts faithfully on the 14 hyperplanes")

    print("\n[2] Transitivity")
    hp_orbits = label_orbits(hp_perms, 14)
    print_orbits("14 hyperplanes under AGL(3,2)", hp_orbits)
    check(len(hp_orbits) == 1 and len(hp_orbits[0]) == 14, "AGL(3,2) is transitive on the 14 weight-4 words")

    perms_no_flip = channel_perms(hp_perms, include_mode_flip=False)
    perms_with_flip = channel_perms(hp_perms, include_mode_flip=True)
    no_flip_orbits = label_orbits(perms_no_flip, 28)
    with_flip_orbits = label_orbits(perms_with_flip, 28)
    print_orbits("28 channels without transverse C2 exchange", no_flip_orbits)
    print_orbits("28 channels with transverse C2 exchange", with_flip_orbits)
    check(len(no_flip_orbits) == 2, "without mode exchange the two transverse sectors remain distinct")
    check(len(with_flip_orbits) == 1 and len(with_flip_orbits[0]) == 28, "AGL(3,2) x C2 is transitive on all 28 channels")

    print("\n[3] Full covariant operator commutant")
    pair14 = ordered_pair_orbits(hp_perms, 14)
    summary14 = orbit_relation_summary(pair14, lambda i, j: relation_14(i, j, planes))
    print(f"  14-channel ordered-pair orbit classes: {dict(summary14)}")
    check(len(pair14) == 3, "commutant dimension on the 14 hyperplanes is 3, not 1")

    pair28 = ordered_pair_orbits(perms_with_flip, 28)
    summary28 = orbit_relation_summary(pair28, lambda i, j: relation_28(i, j, planes))
    print(f"  28-channel ordered-pair orbit classes: {dict(summary28)}")
    check(len(pair28) == 6, "full 28x28 covariant operator commutant has dimension 6, not 1")

    print("\n[4] One-jump reset-instrument invariant weights")
    dim_no_flip = invariant_weight_dimension(perms_no_flip, 28)
    dim_with_flip = invariant_weight_dimension(perms_with_flip, 28)
    check(dim_no_flip == 2, "without transverse C2, one-jump weights have two independent mode rates")
    check(dim_with_flip == 1, "with transverse C2, one-jump weights have one invariant rate")
    print("  Normalisation of the one-jump instrument gives sum_x p_x = 1.")
    print("  Since the invariant weight space is one-dimensional, p_x = 1/28 for every channel.")

    print("\n[5] Theorem status")
    print("  Derived here:")
    print("    AGL(3,2) x C2 transitivity forces uniform one-jump reset weights.")
    print("    Under the serial one-jump constraint, the scheduler is p_x = 1/28.")
    print("  Not derived here:")
    print("    The finite-bandwidth and W-to-28 instrument layers; see companion scripts.")
    print("    The full covariant operator algebra still has six classes, so covariance alone")
    print("    does not exclude parallel or pair-relation-dependent channels.")
    print("\nexit 0 -- uniformity derived from covariance; finite-bandwidth handled separately.")


if __name__ == "__main__":
    main()
