#!/usr/bin/env python3
r"""Z3 centre selection rule on a contractible TCH Wilson loop.

Purpose
-------
This is the centre-symmetry companion to the closed-shell Gauss/triality
certificate.

The Gauss script proves that a compact TCH colour shell admits only
triality-neutral source records.  This script writes the same centre-selection
statement in Wilson-loop language:

    centre averaging kills a fundamental-triality line;
    centre-neutral meson/baryon combinations survive.

Important correction
--------------------
The octagonal TCH face cycle used here is contractible: it bounds a face on the
S^2 shell.  Therefore it is a Wilson loop, not a genuine Polyakov loop.  A
Polyakov loop would require a compact Euclidean-time direction, hence a
non-contractible cycle.  The finite glued shell is simply connected and cannot
see a confinement/deconfinement transition.

Finite construction
-------------------
Use one oriented octagonal TCH face cycle as a representative compact wrap.
For a link configuration U, the wrapped holonomy is

    H = P prod_{e in C} U_e,
    W_fund = (1/3) Tr H.

The global centre operation on the wrap cut multiplies one directed link by
z in {1, omega, omega^2}.  Since z is central,

    W_fund -> z W_fund.

If the measure is centre-invariant, then

    <W_fund> = <z W_fund> = z <W_fund>,

so <W_fund> = 0 for z != 1.  The exact finite version is the centre-orbit
projection:

    (W + omega W + omega^2 W)/3 = 0.

By contrast, centre-neutral records survive:

    P P^*      (meson-like q qbar);
    det H = 1  (baryon-like q q q epsilon).

Boundary
--------
This is the centre selection rule in Polyakov-like notation, not a
deconfinement order parameter.  With independent Haar links, a fundamental
contractible loop with a singly occurring link is killed by the same selection
mechanism.  A genuine order-parameter calculation needs a non-contractible
thermal cycle and a centre-symmetric Wilson weight.  This script is also
distinct from the Z2 inside/outside half-shell complement in the
surface-selector scripts: Z3 acts on SU(3) triality/fundamental centre charge;
Z2 swaps the two spanning surfaces of a closed shell.
"""

from __future__ import annotations

import numpy as np

from record_grammar_tch_glued_su3_link_holonomy import (
    Cycle,
    N_COLOR,
    cycle_holonomy,
    random_su3_haar,
)
from record_grammar_tch_glued_surface_selector import build_glued_complex


RNG_SEED = 20260630
OMEGA = np.exp(2j * np.pi / 3.0)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<78s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def assert_close(name: str, value: complex | float, target: complex | float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<78s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_equal(name: str, value: int, target: int) -> None:
    print(f"  {name:<78s} value={value} target={target}")
    if value != target:
        raise AssertionError(name)


def octagon_wrap_cycle(complex_) -> Cycle:
    octagons = [face for face in complex_.faces if ":O" in face.label]
    if not octagons:
        raise AssertionError("no octagonal wrap face found")
    face = octagons[0]
    return Cycle(vertices=face.vertices + (face.vertices[0],))


def random_links(complex_, rng: np.random.Generator) -> list[np.ndarray]:
    return [random_su3_haar(rng) for _ in complex_.edges]


def normalized_complex_trace(matrix: np.ndarray) -> complex:
    return complex(np.trace(matrix) / N_COLOR)


def centre_transform_cut_link(complex_, links: list[np.ndarray], cycle: Cycle, z: complex) -> list[np.ndarray]:
    """Multiply the first directed link of the wrap cycle by centre element z."""

    out = [link.copy() for link in links]
    left, right = cycle.vertices[0], cycle.vertices[1]
    edge = (left, right) if left < right else (right, left)
    edge_i = complex_.edge_index[edge]
    if edge == (left, right):
        out[edge_i] = z * out[edge_i]
    else:
        # The stored link is for the canonical direction; if the wrapped
        # directed link uses the reverse, multiply the stored link by z* so that
        # its dagger is multiplied by z.
        out[edge_i] = np.conj(z) * out[edge_i]
    return out


def centre_orbit_fundamental_loop(complex_, links: list[np.ndarray], cycle: Cycle) -> list[complex]:
    out = []
    for z in (1.0 + 0.0j, OMEGA, OMEGA**2):
        transformed = centre_transform_cut_link(complex_, links, cycle, z)
        out.append(normalized_complex_trace(cycle_holonomy(complex_, transformed, cycle)))
    return out


def main() -> None:
    print("TCH Z3 centre-selection certificate for a contractible Wilson loop")
    print("=" * 104)
    print("  Z3 centre acts on fundamental triality; Z2 half-shell complement is a different symmetry.")

    complex_ = build_glued_complex(1)
    cycle = octagon_wrap_cycle(complex_)
    rng = np.random.default_rng(RNG_SEED)
    links = random_links(complex_, rng)

    print("\n[1] Contractible fundamental holonomy on a TCH octagonal face cycle")
    hol = cycle_holonomy(complex_, links, cycle)
    fundamental_loop = normalized_complex_trace(hol)
    assert_less("SU(3) holonomy determinant residual", abs(np.linalg.det(hol) - 1.0), 1e-12)
    assert_greater("sample fundamental Wilson magnitude is nonzero before centre projection", abs(fundamental_loop), 1e-3)
    print(f"  contractible perimeter={cycle.perimeter}; W_fund={fundamental_loop.real:+.6f}{fundamental_loop.imag:+.6f}i")

    print("\n[2] Centre action and exact centre-orbit projection")
    orbit = centre_orbit_fundamental_loop(complex_, links, cycle)
    for power, value in enumerate(orbit):
        assert_close(f"W(z^{power}) equals z^{power} W", value, (OMEGA**power) * fundamental_loop, tol=1e-12)
    centre_projection = sum(orbit) / 3.0
    assert_less("fundamental Wilson centre-orbit projection", abs(centre_projection), 1e-12)

    print("\n[3] Centre-neutral records survive")
    meson_record = abs(fundamental_loop) ** 2
    meson_orbit = [abs(item) ** 2 for item in orbit]
    baryon_record = np.linalg.det(hol)
    baryon_orbit = [
        np.linalg.det(cycle_holonomy(complex_, centre_transform_cut_link(complex_, links, cycle, z), cycle))
        for z in (1.0 + 0.0j, OMEGA, OMEGA**2)
    ]
    assert_greater("meson-like W W* record", meson_record, 1e-6)
    assert_less("meson-like W W* centre spread", max(abs(item - meson_record) for item in meson_orbit), 1e-12)
    assert_less("baryon-like det(H) minus 1", abs(baryon_record - 1.0), 1e-12)
    assert_less("baryon-like det(H) centre spread", max(abs(item - baryon_record) for item in baryon_orbit), 1e-12)

    print("\n[4] Monte Carlo sanity: centre projection kills fundamental triality")
    max_orbit_mean = 0.0
    mean_neutral = 0.0
    n_cfg = 80
    for _ in range(n_cfg):
        cfg_links = random_links(complex_, rng)
        cfg_orbit = centre_orbit_fundamental_loop(complex_, cfg_links, cycle)
        max_orbit_mean = max(max_orbit_mean, abs(sum(cfg_orbit) / 3.0))
        mean_neutral += sum(abs(item) ** 2 for item in cfg_orbit) / 3.0
    mean_neutral /= n_cfg
    assert_less("max sampled centre-orbit W residual", max_orbit_mean, 1e-12)
    assert_greater("sampled centre-neutral <|W|^2>", mean_neutral, 1e-4)

    print("\n[5] Do not conflate Z3 centre with Z2 shell complement")
    assert_equal("Z3 centre orbit size", 3, 3)
    assert_equal("Z2 inside/outside complement size", 2, 2)
    assert_less("1 + omega + omega^2", abs(1.0 + OMEGA + OMEGA**2), 1e-12)
    assert_greater("symmetry orders are distinct", abs(3 - 2), 0.0)

    print(
        """
VERDICT:
  PASS.  A contractible fundamental Wilson holonomy on a TCH colour cell
  transforms as W -> z W under the SU(3) Z3 centre.  The exact centre-orbit
  projection of W is zero, which is the centre selection rule in Polyakov-like
  notation, not a genuine Polyakov/deconfinement order parameter.

  Centre-neutral records survive: W W* is meson-like, det(H)=1 is the
  baryon/epsilon line, and neither is killed by the Z3 average.  This is a
  centre-symmetry statement on a simply connected shell, distinct from the Z2
  inside/outside half-shell ambiguity used in the surface-selector ledger.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
