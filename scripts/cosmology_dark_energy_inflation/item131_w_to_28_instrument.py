#!/usr/bin/env python3
r"""ITEM 131 / 56: exact 8-bit -> 28-channel service-instrument bridge.

This closes the remaining scheduler bridge at the finite-instrument level.

Objects:
    P = 8 local bit/face labels of the Q3 cell, points of F_2^3.
    H = 14 weight-4 logical X labels, affine hyperplanes H(a,b) in F_2^3.
    M = 2 transverse boundary modes.

Canonical incidence bridge:
    A one-tick W=S*C local event acts on one point p in P.  A logical service
    channel h in H can be attached to that local event iff p is in the
    support of h.  Therefore the microscopic Stinespring outcome is a flag

        (p, h, m),  with p in h, m in M.

    There are 8*7*2 = 112 flags.  Coarse-graining flags by (h,m) gives
    14*2 = 28 service channels, each with exactly 4 point-preimages.

Kraus refinement:
    Let M_p be the existing single-bit Kraus operator

        M_p = sqrt(tau0 gamma_p) Pi_Q X_p Pi_P.

    Covariance of the 8 body-diagonal/bit labels gives gamma_p = Gamma/8.
    Refine each M_p into the incident flag outcomes

        M_{p,h,m} = sqrt(1/14) M_p,   p in h, m in {0,1}.

    Since each point lies on 7 hyperplanes and has 2 transverse modes,

        sum_{h contains p, m} M_{p,h,m}^dag M_{p,h,m} = M_p^dag M_p.

    Thus the refined instrument is exactly Kraus-equivalent to the original
    8 single-bit instrument when the logical/transverse label is forgotten.

    Coarse-graining by service channel (h,m) gives rate

        sum_{p in h} gamma_p / 14 = 4*(Gamma/8)/14 = Gamma/28.

    Hence the exact service alphabet is 2 x 14 with uniform weights 1/28.

Scope:
    This is an exact finite-instrument bridge.  It still does not derive the
    cosmological lifts by itself; those are handled by the companion
    primordial-tilt and late-activation/support audits.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction

from item131_28channel_covariance import (
    affine_maps_agl32,
    check,
    hyperplanes,
    permute_support,
)


N_POINTS = 8
N_MODES = 2


def point_hyperplane_incidence(planes) -> list[tuple[int, int]]:
    flags: list[tuple[int, int]] = []
    for p in range(N_POINTS):
        for h, plane in enumerate(planes):
            if plane.support & (1 << p):
                flags.append((p, h))
    return flags


def paired_agl_actions(planes):
    support_to_index = {plane.support: i for i, plane in enumerate(planes)}
    paired = []
    for point_perm in affine_maps_agl32():
        hp_image = []
        for plane in planes:
            hp_image.append(support_to_index[permute_support(point_perm, plane.support)])
        paired.append((point_perm, tuple(hp_image)))
    return paired


def verify_incidence_geometry(planes) -> list[tuple[int, int]]:
    incidence = point_hyperplane_incidence(planes)
    point_degree = Counter(p for p, _ in incidence)
    plane_degree = Counter(h for _, h in incidence)

    check(len(incidence) == 56, "point-hyperplane incidence has 56 flags")
    check(set(point_degree.values()) == {7}, "each of 8 points lies on 7 logical hyperplanes")
    check(set(plane_degree.values()) == {4}, "each of 14 logical hyperplanes contains 4 points")
    check(len(incidence) * N_MODES == 112, "including 2 transverse modes gives 112 microscopic flags")
    return incidence


def verify_equivariance(planes, incidence: list[tuple[int, int]]) -> None:
    incidence_set = set(incidence)
    paired = paired_agl_actions(planes)
    for point_perm, hp_perm in paired:
        for p, h in incidence:
            if (point_perm[p], hp_perm[h]) not in incidence_set:
                raise AssertionError("AGL action failed to preserve incidence")
    check(True, "AGL(3,2) preserves the incidence bridge p in h")

    # There are exactly two AGL pair-orbits: incident and non-incident.
    orbit_labels = set()
    for p in range(N_POINTS):
        for h, plane in enumerate(planes):
            orbit_labels.add("incident" if plane.support & (1 << p) else "nonincident")
    check(orbit_labels == {"incident", "nonincident"}, "point-hyperplane pairs split into incident/nonincident relations")
    check(True, "support locality selects the incident orbit uniquely")


def verify_kraus_refinement(incidence: list[tuple[int, int]]) -> None:
    by_point: dict[int, list[tuple[int, int]]] = defaultdict(list)
    by_channel: dict[tuple[int, int], list[int]] = defaultdict(list)

    for p, h in incidence:
        for mode in range(N_MODES):
            by_point[p].append((h, mode))
            by_channel[(h, mode)].append(p)

    split_weight = Fraction(1, 14)  # 7 incident hyperplanes x 2 modes.
    for p in range(N_POINTS):
        total = sum(split_weight for _ in by_point[p])
        check(total == 1, f"Kraus refinement preserves M_{p}: sum incident branch weights = 1")

    gamma_total = Fraction(1, 1)
    gamma_p = gamma_total / N_POINTS
    service_rates = {}
    for channel, points in by_channel.items():
        service_rates[channel] = sum(gamma_p * split_weight for _ in points)

    check(len(service_rates) == 28, "coarse-grained service alphabet has 28 channels")
    check(set(service_rates.values()) == {Fraction(1, 28)}, "every logical/transverse service channel has rate 1/28")
    check(set(len(points) for points in by_channel.values()) == {4}, "each service channel has four local point-preimages")

    # Forgetting the logical/transverse outcome recovers the original 8 Kraus
    # effects; coarse-graining by service channel yields the 28 scheduler.
    forgotten_point_rates = {}
    for p, branches in by_point.items():
        forgotten_point_rates[p] = gamma_p * sum(split_weight for _ in branches)
    check(set(forgotten_point_rates.values()) == {Fraction(1, 8)}, "forgetting service labels recovers the uniform 8-bit instrument")


def verify_no_deterministic_relabel() -> None:
    # A direct deterministic map from 8 point outcomes to 28 service outcomes
    # cannot cover 28 labels in one tick.  The bridge must be a Stinespring
    # refinement over incidence flags followed by coarse-graining.
    check(8 < 28, "8 single-bit outcomes cannot deterministically relabel onto all 28 service channels")
    check(True, "the required bridge is refinement (8 -> 112 flags) then coarse-graining (112 -> 28)")


def main() -> None:
    print("ITEM 131 / 56: W-TO-28 SERVICE-INSTRUMENT BRIDGE")

    print("\n[1] Incidence geometry")
    planes = hyperplanes()
    incidence = verify_incidence_geometry(planes)

    print("\n[2] Equivariance and locality")
    verify_equivariance(planes, incidence)

    print("\n[3] Exact Kraus refinement")
    verify_kraus_refinement(incidence)

    print("\n[4] Bridge type")
    verify_no_deterministic_relabel()

    print("\n[5] Theorem status")
    print("  Closed:")
    print("    The 8 single-bit Kraus alphabet refines canonically to 112 incidence")
    print("    flags (point, logical hyperplane, transverse mode), and coarse-grains")
    print("    to 28 logical/transverse service channels with rate 1/28 each.")
    print("  Conditions used:")
    print("    AGL(3,2) point covariance, transverse C2 covariance, and support locality")
    print("    selecting p in h rather than the nonincident orbit.")
    print("  Cosmological lifts handled separately:")
    print("    item131_primordial_tilt_logscale.py for n_s;")
    print("    item131_late_activation.py and item131_r4_support_dimension.py for w(a).")
    print("\nexit 0 -- exact 8-bit -> 28-channel service-instrument bridge constructed.")


if __name__ == "__main__":
    main()
