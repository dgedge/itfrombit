#!/usr/bin/env python3
"""Non-abelian path ordering: the first SU(2) Wilson-loop toy.

This script adds the first genuinely non-abelian rung to the record-grammar
ladder.

Abelian lesson
--------------
For U(1), link variables are phases and commute.  A closed loop is a phase:

    W = exp(i Phi).

Non-abelian lesson
------------------
For SU(2), link variables are matrices.  The order of links matters:

    U_x U_y != U_y U_x.

The smallest loop commutator is

    W = U_x U_y U_x^dag U_y^dag.

The raw matrix W is not itself a gauge-invariant number.  Under local gauge
rotations at the plaquette vertices, the loop based at vertex 0 transforms by
conjugation:

    W -> G_0 W G_0^dag.

Therefore trace and eigenvalues are the detector-readable loop records:

    Tr(W), spec(W).

Detector reading
----------------
For a two-route interferometer with unresolved internal SU(2) color and a
maximally mixed incoming color, the bright/dark probabilities are

    P_bright = 1/2 + Re Tr(W)/(2N),
    P_dark   = 1/2 - Re Tr(W)/(2N),       N = 2.

So the non-abelian detector hears the normalized Wilson-loop trace, not a raw
endpoint phase.

This is a finite SU(2) path-ordering toy.  It is not a Yang-Mills action and it
does not claim to derive non-abelian gauge dynamics.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def su2_exp(generator: np.ndarray, angle: float) -> np.ndarray:
    """exp(i angle sigma) for a Pauli generator sigma."""

    return math.cos(angle) * I2 + 1j * math.sin(angle) * generator


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<70s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<70s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<70s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def check_su2(name: str, u: np.ndarray) -> None:
    assert_matrix_close(f"{name}: U^dag U = I", u.conj().T @ u, I2)
    assert_close(f"{name}: det U", abs(np.linalg.det(u)), 1.0)


def loop_commutator(ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
    """Group commutator Ux Uy Ux^dag Uy^dag."""

    return ux @ uy @ ux.conj().T @ uy.conj().T


def plaquette_links(ux: np.ndarray, uy: np.ndarray) -> dict[tuple[int, int], np.ndarray]:
    """Square links 0->1->2->3->0 with lower route 0->3->2."""

    return {
        (0, 1): ux,
        (1, 2): uy,
        (2, 3): ux.conj().T,
        (3, 0): uy.conj().T,
    }


def gauge_transform_links(
    links: dict[tuple[int, int], np.ndarray],
    gauges: dict[int, np.ndarray],
) -> dict[tuple[int, int], np.ndarray]:
    """Local lattice-gauge transformation U_ij -> G_i U_ij G_j^dag."""

    return {(i, j): gauges[i] @ u @ gauges[j].conj().T for (i, j), u in links.items()}


def plaquette_loop(links: dict[tuple[int, int], np.ndarray]) -> np.ndarray:
    """Loop based at vertex 0: 0->1->2->3->0."""

    return links[(0, 1)] @ links[(1, 2)] @ links[(2, 3)] @ links[(3, 0)]


def upper_route(links: dict[tuple[int, int], np.ndarray]) -> np.ndarray:
    return links[(0, 1)] @ links[(1, 2)]


def lower_route(links: dict[tuple[int, int], np.ndarray]) -> np.ndarray:
    # The stored links include 3->0 and 2->3.  The lower 0->3->2 route uses
    # their Hermitian conjugates.
    return links[(3, 0)].conj().T @ links[(2, 3)].conj().T


def bright_dark_from_routes(upper: np.ndarray, lower: np.ndarray) -> tuple[float, float]:
    """Balanced two-route probabilities with unresolved maximally mixed color."""

    n = upper.shape[0]
    rho = np.eye(n, dtype=complex) / n
    bright_op = (upper + lower) / 2.0
    dark_op = (upper - lower) / 2.0
    p_bright = float(np.real_if_close(np.trace(bright_op @ rho @ bright_op.conj().T)))
    p_dark = float(np.real_if_close(np.trace(dark_op @ rho @ dark_op.conj().T)))
    return p_bright, p_dark


def sorted_eigenvalues(u: np.ndarray) -> np.ndarray:
    vals = np.linalg.eigvals(u)
    return np.array(sorted(vals, key=lambda z: (round(float(np.angle(z)), 12), round(float(abs(z)), 12))))


def main() -> None:
    print("Non-abelian path ordering: SU(2) Wilson-loop toy")
    print("=" * 88)

    theta = 0.67
    phi = -0.43
    ux = su2_exp(X, theta)
    uy = su2_exp(Y, phi)
    check_su2("Ux", ux)
    check_su2("Uy", uy)

    print("\n[1] Path order matters")
    xy = ux @ uy
    yx = uy @ ux
    comm_norm = float(np.linalg.norm(xy - yx))
    print(f"  ||Ux Uy - Uy Ux|| = {comm_norm:.12g}")
    if comm_norm < 1e-6:
        raise AssertionError("chosen SU(2) links accidentally commuted")
    print("  -> a non-abelian path is an ordered word, not just a sum of phases.")

    print("\n[2] The plaquette loop is a nontrivial SU(2) matrix")
    w = loop_commutator(ux, uy)
    check_su2("W", w)
    assert_matrix_close("W from routes = Uupper Ulower^dag", w, xy @ yx.conj().T)
    scalar_part = np.trace(w) / 2.0
    noncentral_norm = float(np.linalg.norm(w - scalar_part * I2))
    print(f"  ||W - Tr(W)I/2|| = {noncentral_norm:.12g}")
    if noncentral_norm < 1e-6:
        raise AssertionError("loop was accidentally central")
    print("  -> unlike U(1), the loop record is not merely a number before taking trace/eigenvalues.")

    print("\n[3] Local gauge transformations conjugate the loop, preserving trace/eigenvalues")
    links = plaquette_links(ux, uy)
    w_from_links = plaquette_loop(links)
    assert_matrix_close("plaquette links reproduce W", w_from_links, w)

    gauges = {
        0: su2_exp(Z, 0.22) @ su2_exp(X, -0.31),
        1: su2_exp(Y, 0.74),
        2: su2_exp(X, -0.58) @ su2_exp(Z, 0.19),
        3: su2_exp(Z, -0.41),
    }
    transformed = gauge_transform_links(links, gauges)
    w_g = plaquette_loop(transformed)
    assert_matrix_close("W' = G0 W G0^dag", w_g, gauges[0] @ w @ gauges[0].conj().T)
    assert_complex_close("Tr(W') = Tr(W)", np.trace(w_g), np.trace(w))
    eig = sorted_eigenvalues(w)
    eig_g = sorted_eigenvalues(w_g)
    assert_matrix_close("spec(W') = spec(W)", eig_g, eig, tol=1e-10)

    raw_change = float(np.linalg.norm(w_g - w))
    print(f"  ||W' - W|| = {raw_change:.12g}")
    if raw_change < 1e-6:
        raise AssertionError("chosen base gauge accidentally left raw W unchanged")
    print("  -> the raw based-loop matrix is gauge-covariant; trace/eigenvalues are gauge-invariant.")

    print("\n[4] Unresolved-color detector hears Re Tr(W), not a raw matrix element")
    upper = upper_route(links)
    lower = lower_route(links)
    p_bright, p_dark = bright_dark_from_routes(upper, lower)
    n = 2
    expected_bright = 0.5 + float(np.real(np.trace(w))) / (2.0 * n)
    expected_dark = 0.5 - float(np.real(np.trace(w))) / (2.0 * n)
    assert_close("P_bright", p_bright, expected_bright)
    assert_close("P_dark", p_dark, expected_dark)
    assert_close("P_bright + P_dark", p_bright + p_dark, 1.0)

    upper_g = upper_route(transformed)
    lower_g = lower_route(transformed)
    p_bright_g, p_dark_g = bright_dark_from_routes(upper_g, lower_g)
    assert_close("P_bright gauge invariant", p_bright_g, p_bright)
    assert_close("P_dark gauge invariant", p_dark_g, p_dark)
    print("  -> a color-blind detector reads the normalized Wilson-loop trace.")

    print("\n[5] The small-loop limit exposes the commutator/curvature term")
    eps_x = 1.0e-3
    eps_y = -1.7e-3
    ux_small = su2_exp(X, eps_x)
    uy_small = su2_exp(Y, eps_y)
    w_small = loop_commutator(ux_small, uy_small)
    # With exp(i eps_x X), exp(i eps_y Y):
    # [i eps_x X, i eps_y Y] = -eps_x eps_y [X,Y] = -2i eps_x eps_y Z.
    first_commutator = I2 - 2j * eps_x * eps_y * Z
    assert_matrix_close("W = I - 2i eps_x eps_y Z + higher order", w_small, first_commutator, tol=1e-8)
    print("  -> the first nontrivial loop content is the generator commutator.")

    print(
        """
Verdict
-------
The non-abelian rung changes the grammar:

  U(1):  closed loop = one phase.
  SU(2): closed loop = ordered matrix product; trace/eigenvalues are readable.

The raw path word is gauge-covariant, not gauge-invariant.  A detector that
does not resolve internal color hears Re Tr(W), while the small-loop expansion
hears the commutator.  This is the finite record-grammar version of why
non-abelian gauge fields carry path ordering and curvature.

Boundary
--------
This script is not a Yang-Mills action or a continuum gauge theory.  It only
shows the operator grammar that must exist before such dynamics can be posed.
"""
    )


if __name__ == "__main__":
    main()
