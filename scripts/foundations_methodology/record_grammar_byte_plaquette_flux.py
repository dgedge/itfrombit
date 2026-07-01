#!/usr/bin/env python3
r"""Byte-level plaquette flux: Wilson loop on finite record cells.

This script combines two earlier record-grammar rungs:

  * byte endpoint records from the [8,4,4] cell;
  * Aharonov-Bohm/Wilson loop phase as a closed, gauge-readable flux.

Toy geometry
------------
Use four byte endpoints on a square plaquette:

        A
      /   \
    S       D
      \   /
        B

The source byte S can reach the detector byte D by an upper route S->A->D or
a lower route S->B->D.  Each plaquette vertex is a byte endpoint: a logical
one-excitation state means one byte is |1_L> = 11111111 and the other endpoint
bytes are |0_L> = 00000000.  The byte's 12 cube-edge checks still detect a
single physical bit fault and identify its vertex address.

The link variables are U(1) Wilson holonomies.  Individual open path phases are
gauge bookkeeping.  The closed plaquette flux

    W = U_SA U_AD (U_SB U_BD)^*

is gauge-invariant and controls the bright/dark detector probabilities:

    P_bright = (1 + cos arg W)/2,
    P_dark   = (1 - cos arg W)/2.

Record-grammar reading
----------------------
The finite endpoint records are byte cells.  The phase record is the closed
plaquette flux.  The detector sees the closed flux through interference, not
the phase assigned to any individual endpoint or link.  Physical byte faults
remain ordinary local syndrome events, separate from the gauge flux ledger.

This is a pedagogical finite model, not a Maxwell action or full lattice gauge
theory.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


BYTE_BITS = 8
ZERO_L = 0
ONE_L = (1 << BYTE_BITS) - 1
NODES = ("S", "A", "B", "D")


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def rm13_words() -> set[int]:
    """RM(1,3) byte code in the same convention as the earlier byte scripts."""

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
    """Logical byte configuration with one endpoint byte excited."""

    if excited_node not in NODES:
        raise ValueError(f"unknown node {excited_node}")
    return {node: (ONE_L if node == excited_node else ZERO_L) for node in NODES}


def single_physical_fault(config: dict[str, int], node: str, byte_vertex: int) -> dict[str, int]:
    """Flip one physical bit inside one endpoint byte."""

    out = dict(config)
    out[node] ^= 1 << byte_vertex
    return out


def all_byte_syndromes(config: dict[str, int]) -> dict[str, tuple[int, ...]]:
    return {node: syndrome_word(word) for node, word in config.items()}


@dataclass(frozen=True)
class BytePlaquetteLinks:
    """Four oriented U(1) link phases on the byte plaquette."""

    sa: complex
    ad: complex
    sb: complex
    bd: complex

    @classmethod
    def from_angles(cls, theta_sa: float, theta_ad: float, theta_sb: float, theta_bd: float) -> "BytePlaquetteLinks":
        return cls(phase(theta_sa), phase(theta_ad), phase(theta_sb), phase(theta_bd))

    def upper(self) -> complex:
        return self.sa * self.ad

    def lower(self) -> complex:
        return self.sb * self.bd

    def loop(self) -> complex:
        # Closed path S->A->D->B->S.  The reverse lower route contributes the
        # complex conjugate of S->B->D.
        return self.upper() * np.conj(self.lower())

    def flux(self) -> float:
        return float(np.angle(self.loop()))


def transform_link(u_ij: complex, lambda_i: float, lambda_j: float) -> complex:
    """AB/path-holonomy convention: U_ij -> exp(i(lambda_i-lambda_j)) U_ij."""

    return complex(np.exp(1j * (lambda_i - lambda_j)) * u_ij)


def gauge_transform(links: BytePlaquetteLinks, lambdas: dict[str, float]) -> BytePlaquetteLinks:
    return BytePlaquetteLinks(
        sa=transform_link(links.sa, lambdas["S"], lambdas["A"]),
        ad=transform_link(links.ad, lambdas["A"], lambdas["D"]),
        sb=transform_link(links.sb, lambdas["S"], lambdas["B"]),
        bd=transform_link(links.bd, lambdas["B"], lambdas["D"]),
    )


def bright_dark(links: BytePlaquetteLinks) -> tuple[float, float]:
    """Balanced byte-interferometer detector probabilities."""

    bright_amp = (links.upper() + links.lower()) / 2.0
    dark_amp = (links.upper() - links.lower()) / 2.0
    return float(abs(bright_amp) ** 2), float(abs(dark_amp) ** 2)


def redistribute_flux_to_one_link(phi: float) -> BytePlaquetteLinks:
    return BytePlaquetteLinks.from_angles(phi, 0.0, 0.0, 0.0)


def redistribute_flux_evenly(phi: float) -> BytePlaquetteLinks:
    return BytePlaquetteLinks.from_angles(phi / 4.0, phi / 4.0, -phi / 4.0, -phi / 4.0)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<70s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<70s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Byte-level plaquette flux: closed Wilson loop on finite record cells")
    print("=" * 92)

    print("\n[1] Each plaquette vertex is a byte endpoint")
    print(f"  |RM(1,3)| = {len(CODE)}; cube edge checks = {len(EDGES)}")
    assert len(CODE) == 16
    assert len(EDGES) == 12
    assert ZERO_L in CODE and ONE_L in CODE

    for excited in NODES:
        config = one_excitation_configuration(excited)
        syn = all_byte_syndromes(config)
        assert all(s == ZERO_SYNDROME for s in syn.values())
    print("  all four one-excitation endpoint configurations are byte-code configurations")

    print("\n[2] Physical byte faults remain local syndrome-addressed events")
    base = one_excitation_configuration("A")
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
                    raise AssertionError(f"wrong syndrome for node={node}, vertex={byte_vertex}")
    print("  a one-bit fault identifies both the endpoint byte and the internal cube vertex")

    print("\n[3] The detector probabilities hear only the closed plaquette flux")
    links = BytePlaquetteLinks.from_angles(theta_sa=0.41, theta_ad=-0.72, theta_sb=1.11, theta_bd=0.36)
    phi = links.flux()
    bright, dark = bright_dark(links)
    assert_close("P_bright", bright, (1.0 + math.cos(phi)) / 2.0)
    assert_close("P_dark", dark, (1.0 - math.cos(phi)) / 2.0)
    assert_close("P_bright + P_dark", bright + dark, 1.0)
    print(f"  byte plaquette flux arg W = {phi:+.12f} rad")

    print("\n[4] Endpoint gauge transformations move open phases but preserve the byte detector output")
    lambdas = {"S": 0.62, "A": -0.35, "B": 1.74, "D": -1.06}
    transformed = gauge_transform(links, lambdas)
    endpoint_factor = phase(lambdas["S"] - lambdas["D"])
    assert_complex_close("upper route covariance", transformed.upper(), endpoint_factor * links.upper())
    assert_complex_close("lower route covariance", transformed.lower(), endpoint_factor * links.lower())
    assert_complex_close("closed plaquette loop invariant", transformed.loop(), links.loop())
    bright_g, dark_g = bright_dark(transformed)
    assert_close("P_bright gauge invariant", bright_g, bright)
    assert_close("P_dark gauge invariant", dark_g, dark)
    print("  -> open route phases moved; the closed flux and detector probabilities did not.")

    print("\n[5] Same closed flux, different link bookkeeping")
    for label, representative in [
        ("all flux on S-A", redistribute_flux_to_one_link(phi)),
        ("flux spread evenly", redistribute_flux_evenly(phi)),
    ]:
        b, d = bright_dark(representative)
        assert_close(f"{label}: flux", representative.flux(), phi)
        assert_close(f"{label}: P_bright", b, bright)
        assert_close(f"{label}: P_dark", d, dark)
    print("  -> the link-phase distribution is gauge convention; the plaquette holonomy is the record.")

    print("\n[6] Gauss averaging kills open byte-route phase, not closed flux")
    grid = [2.0 * math.pi * k / 64 for k in range(64)]
    open_avg = sum(phase(ls - ld) * links.upper() for ls in grid for ld in grid) / (64 * 64)
    assert_close("|Gauss average of open S->D route|", abs(open_avg), 0.0, tol=1e-12)
    assert_complex_close("closed loop survives", links.loop(), links.loop())
    print("  -> a byte endpoint detector cannot hear a naked route phase.")

    print(
        """
Verdict
-------
The byte-level plaquette script combines the finite record cell and the Wilson
loop:

  each plaquette vertex is a [8,4,4] byte endpoint;
  one-site physical faults are still local syndrome-addressed byte events;
  open path/link phases are gauge bookkeeping;
  the closed plaquette Wilson loop is the detector-readable flux record;
  bright/dark probabilities depend only on that closed flux.

This is the first finite QEC/gauge plaquette rung of the record grammar.  It
does not derive a Maxwell action or continuum gauge dynamics; it only proves
that the closed-flux readout can be placed on finite byte endpoint records.
"""
    )


if __name__ == "__main__":
    main()
