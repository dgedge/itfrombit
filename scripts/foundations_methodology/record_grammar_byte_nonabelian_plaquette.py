#!/usr/bin/env python3
r"""Byte-level non-abelian plaquette: SU(2) Wilson loop on finite records.

This script combines three earlier record-grammar rungs:

  * [8,4,4] byte endpoint records and local syndrome checks;
  * non-abelian SU(2) path ordering;
  * detector readout by the normalized Wilson-loop trace.

Toy geometry
------------
Use four byte endpoints on a square plaquette:

        A
      /   \
    S       D
      \   /
        B

The source byte S can reach detector byte D by an upper route S->A->D or a
lower route S->B->D.  Each plaquette vertex is a byte endpoint, so a one-site
physical bit fault still has a local [8,4,4] syndrome address.

The link records are now SU(2) matrices.  The open routes are ordered products,
not commuting phases.  The closed loop based at S is

    W = U_upper U_lower^dag,

and transforms by conjugation under local SU(2) gauge rotations at S.  Therefore
the gauge-readable data are Tr(W) and spec(W).  A color-blind detector with
maximally mixed incoming color reads

    P_bright = 1/2 + Re Tr(W)/(2N),     N = 2.

This is a finite QEC/gauge toy.  It is not an SU(2) Yang-Mills action or a
claim about the full framework's non-abelian sector.
"""

from __future__ import annotations

import math

import numpy as np


BYTE_BITS = 8
ZERO_L = 0
ONE_L = (1 << BYTE_BITS) - 1
NODES = ("S", "A", "B", "D")

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def rm13_words() -> set[int]:
    gens = [
        0b11111111,
        0b11110000,
        0b11001100,
        0b10101010,
    ]
    out = {0}
    for g in gens:
        out |= {x ^ g for x in list(out)}
    return out


CODE = rm13_words()


def bits(state: int) -> tuple[int, ...]:
    return tuple((state >> i) & 1 for i in range(BYTE_BITS))


def cube_edges() -> list[tuple[int, int]]:
    return [(u, v) for u in range(8) for v in range(u + 1, 8) if (u ^ v).bit_count() == 1]


EDGES = cube_edges()


def syndrome_word(state: int) -> tuple[int, ...]:
    b = bits(state)
    return tuple(b[u] ^ b[v] for u, v in EDGES)


ZERO_SYNDROME = syndrome_word(ZERO_L)
INCIDENT = {v: tuple(1 if v in edge else 0 for edge in EDGES) for v in range(8)}


def one_excitation_configuration(excited_node: str) -> dict[str, int]:
    if excited_node not in NODES:
        raise ValueError(f"unknown node {excited_node}")
    return {node: (ONE_L if node == excited_node else ZERO_L) for node in NODES}


def single_physical_fault(config: dict[str, int], node: str, byte_vertex: int) -> dict[str, int]:
    out = dict(config)
    out[node] ^= 1 << byte_vertex
    return out


def all_byte_syndromes(config: dict[str, int]) -> dict[str, tuple[int, ...]]:
    return {node: syndrome_word(word) for node, word in config.items()}


def su2_exp(generator: np.ndarray, angle: float) -> np.ndarray:
    return math.cos(angle) * I2 + 1j * math.sin(angle) * generator


def check_su2(name: str, u: np.ndarray) -> None:
    assert_matrix_close(f"{name}: U^dag U = I", u.conj().T @ u, I2)
    assert_close(f"{name}: |det U|", abs(np.linalg.det(u)), 1.0)


def byte_su2_links() -> dict[tuple[str, str], np.ndarray]:
    """A noncommuting SU(2) plaquette on byte endpoints."""

    ux = su2_exp(X, 0.61)
    uy = su2_exp(Y, -0.47)
    return {
        ("S", "A"): ux,
        ("A", "D"): uy,
        ("B", "D"): ux,
        ("S", "B"): uy,
    }


def dagger_reverse_links(links: dict[tuple[str, str], np.ndarray]) -> dict[tuple[str, str], np.ndarray]:
    out = dict(links)
    for (i, j), u in list(links.items()):
        out[(j, i)] = u.conj().T
    return out


def gauge_transform_links(
    links: dict[tuple[str, str], np.ndarray],
    gauges: dict[str, np.ndarray],
) -> dict[tuple[str, str], np.ndarray]:
    return {(i, j): gauges[i] @ u @ gauges[j].conj().T for (i, j), u in links.items()}


def upper_route(links: dict[tuple[str, str], np.ndarray]) -> np.ndarray:
    return links[("S", "A")] @ links[("A", "D")]


def lower_route(links: dict[tuple[str, str], np.ndarray]) -> np.ndarray:
    return links[("S", "B")] @ links[("B", "D")]


def loop_from_routes(links: dict[tuple[str, str], np.ndarray]) -> np.ndarray:
    return upper_route(links) @ lower_route(links).conj().T


def bright_dark_from_routes(up: np.ndarray, low: np.ndarray) -> tuple[float, float]:
    n = up.shape[0]
    rho = np.eye(n, dtype=complex) / n
    bright_op = (up + low) / 2.0
    dark_op = (up - low) / 2.0
    p_bright = float(np.real_if_close(np.trace(bright_op @ rho @ bright_op.conj().T)))
    p_dark = float(np.real_if_close(np.trace(dark_op @ rho @ dark_op.conj().T)))
    return p_bright, p_dark


def sorted_eigenvalues(u: np.ndarray) -> np.ndarray:
    vals = np.linalg.eigvals(u)
    return np.array(sorted(vals, key=lambda z: (round(float(np.angle(z)), 12), round(float(abs(z)), 12))))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<76s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Byte-level non-abelian plaquette: SU(2) loop on finite record cells")
    print("=" * 96)

    print("\n[1] The byte endpoint layer is still present")
    assert len(CODE) == 16
    assert len(EDGES) == 12
    assert ZERO_L in CODE and ONE_L in CODE
    for excited in NODES:
        config = one_excitation_configuration(excited)
        syn = all_byte_syndromes(config)
        assert all(s == ZERO_SYNDROME for s in syn.values())
    print("  all four one-excitation endpoint configurations are valid byte records")

    print("\n[2] Physical byte faults are separate from SU(2) link records")
    base = one_excitation_configuration("S")
    for node in NODES:
        for byte_vertex in range(8):
            bad = single_physical_fault(base, node, byte_vertex)
            delta = {
                n: tuple(a ^ b for a, b in zip(all_byte_syndromes(base)[n], all_byte_syndromes(bad)[n]))
                for n in NODES
            }
            for n in NODES:
                expected = INCIDENT[byte_vertex] if n == node else ZERO_SYNDROME
                if delta[n] != expected:
                    raise AssertionError(f"wrong byte syndrome for node={node}, vertex={byte_vertex}")
    print("  a physical bit fault identifies endpoint byte plus internal cube vertex")

    print("\n[3] SU(2) open routes are ordered products")
    links = byte_su2_links()
    links = dagger_reverse_links(links)
    for edge in [("S", "A"), ("A", "D"), ("S", "B"), ("B", "D")]:
        check_su2(f"U_{edge[0]}{edge[1]}", links[edge])
    up = upper_route(links)
    low = lower_route(links)
    order_gap = float(np.linalg.norm(up - low))
    print(f"  ||U_upper - U_lower|| = {order_gap:.12g}")
    if order_gap < 1e-6:
        raise AssertionError("chosen byte SU(2) routes accidentally matched")
    print("  -> the same endpoint pair can carry distinct ordered route records.")

    print("\n[4] The closed byte plaquette loop is gauge-covariant, not a scalar phase")
    w = loop_from_routes(links)
    check_su2("W", w)
    scalar_part = np.trace(w) / 2.0
    noncentral = float(np.linalg.norm(w - scalar_part * I2))
    print(f"  ||W - Tr(W)I/2|| = {noncentral:.12g}")
    if noncentral < 1e-6:
        raise AssertionError("byte SU(2) loop accidentally central")

    gauges = {
        "S": su2_exp(Z, 0.28) @ su2_exp(X, -0.19),
        "A": su2_exp(Y, -0.51),
        "D": su2_exp(X, 0.37) @ su2_exp(Z, 0.22),
        "B": su2_exp(Z, -0.44),
    }
    transformed = gauge_transform_links(links, gauges)
    w_g = loop_from_routes(transformed)
    assert_matrix_close("W' = G_S W G_S^dag", w_g, gauges["S"] @ w @ gauges["S"].conj().T)
    assert_complex_close("Tr(W') = Tr(W)", np.trace(w_g), np.trace(w))
    assert_matrix_close("spec(W') = spec(W)", sorted_eigenvalues(w_g), sorted_eigenvalues(w), tol=1e-10)
    if float(np.linalg.norm(w_g - w)) < 1e-6:
        raise AssertionError("chosen gauge accidentally left raw W unchanged")
    print("  -> trace/eigenvalues are readable; the based-loop matrix is only covariant.")

    print("\n[5] A color-blind byte detector reads normalized Re Tr(W)")
    p_bright, p_dark = bright_dark_from_routes(up, low)
    n = 2
    expected_bright = 0.5 + float(np.real(np.trace(w))) / (2.0 * n)
    expected_dark = 0.5 - float(np.real(np.trace(w))) / (2.0 * n)
    assert_close("P_bright", p_bright, expected_bright)
    assert_close("P_dark", p_dark, expected_dark)
    assert_close("P_bright + P_dark", p_bright + p_dark, 1.0)

    p_bright_g, p_dark_g = bright_dark_from_routes(upper_route(transformed), lower_route(transformed))
    assert_close("P_bright gauge invariant", p_bright_g, p_bright)
    assert_close("P_dark gauge invariant", p_dark_g, p_dark)
    print("  -> finite byte endpoints plus unresolved color give the same trace readout rule.")

    print("\n[6] Small-loop content is the same non-abelian commutator on byte records")
    eps_x = 9.0e-4
    eps_y = -1.4e-3
    small_links = dagger_reverse_links(
        {
            ("S", "A"): su2_exp(X, eps_x),
            ("A", "D"): su2_exp(Y, eps_y),
            ("B", "D"): su2_exp(X, eps_x),
            ("S", "B"): su2_exp(Y, eps_y),
        }
    )
    w_small = loop_from_routes(small_links)
    first_commutator = I2 - 2j * eps_x * eps_y * Z
    assert_matrix_close("W = I - 2i eps_x eps_y Z + higher order", w_small, first_commutator, tol=1e-8)
    print("  -> the non-abelian curvature signal lives on the Wilson loop, not inside byte faults.")

    print(
        """
Verdict
-------
The byte-level non-abelian toy combines the finite record cell and SU(2)
Wilson-loop grammar:

  each plaquette vertex is a [8,4,4] byte endpoint;
  one-site physical faults remain byte-syndrome events;
  SU(2) route records are ordered products;
  the closed loop is gauge-covariant, with trace/eigenvalues gauge-readable;
  a color-blind detector reads normalized Re Tr(W);
  the small-loop content is the generator commutator.

Boundary
--------
This is still a toy.  It does not derive Yang-Mills dynamics, color confinement,
or the Standard Model gauge group.  It only proves that the non-abelian Wilson
grammar can be placed on finite byte endpoint records without confusing link
curvature with local byte syndromes.
"""
    )


if __name__ == "__main__":
    main()
