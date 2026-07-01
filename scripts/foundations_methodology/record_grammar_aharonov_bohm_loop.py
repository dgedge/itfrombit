#!/usr/bin/env python3
r"""Aharonov-Bohm record machine: phase as closed flux, not endpoint property.

This script adds the next small gauge example to the record-grammar ladder.

Physical picture
----------------
Use the smallest interferometer graph:

        A
      /   \
    S       D
      \   /
        B

A charged record can travel from source S to detector D by an upper path
S->A->D or a lower path S->B->D.  Each oriented link carries a U(1) Wilson
holonomy U_ij = exp(i theta_ij).  Under a local gauge transformation,

    U_ij -> exp(i(lambda_i - lambda_j)) U_ij.

The individual open-path phases are not physical: each path amplitude transforms
by the same endpoint phase exp(i(lambda_S - lambda_D)).  The relative phase
between the two paths is the closed Wilson loop

    W = U_SA U_AD (U_SB U_BD)^*
      = exp(i Phi_AB),

which is gauge-invariant.  The detector probabilities are

    P_bright = (1 + cos Phi_AB)/2,
    P_dark   = (1 - cos Phi_AB)/2.

Record-grammar reading
----------------------
The open path phases are "unsayable" gauge bookkeeping.  The closed loop flux is
the readable record.  This is the clean Aharonov-Bohm lesson in the same
endpoint/holonomy language used by the Wilson-dressed Bell scripts.

This is a toy graph calculation, not a full Maxwell field calculation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


NODES = ("S", "A", "B", "D")


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


@dataclass(frozen=True)
class DiamondLinks:
    """Four oriented link phases for the diamond interferometer."""

    sa: complex
    ad: complex
    sb: complex
    bd: complex

    @classmethod
    def from_angles(cls, theta_sa: float, theta_ad: float, theta_sb: float, theta_bd: float) -> "DiamondLinks":
        return cls(phase(theta_sa), phase(theta_ad), phase(theta_sb), phase(theta_bd))

    def upper(self) -> complex:
        return self.sa * self.ad

    def lower(self) -> complex:
        return self.sb * self.bd

    def loop(self) -> complex:
        # Closed loop S->A->D->B->S.  The reverse lower path contributes the
        # complex conjugate because U_ji = U_ij^* for U(1).
        return self.upper() * np.conj(self.lower())

    def flux(self) -> float:
        return float(np.angle(self.loop()))


def transform_link(u_ij: complex, lambda_i: float, lambda_j: float) -> complex:
    return complex(np.exp(1j * (lambda_i - lambda_j)) * u_ij)


def gauge_transform(links: DiamondLinks, lambdas: dict[str, float]) -> DiamondLinks:
    return DiamondLinks(
        sa=transform_link(links.sa, lambdas["S"], lambdas["A"]),
        ad=transform_link(links.ad, lambdas["A"], lambdas["D"]),
        sb=transform_link(links.sb, lambdas["S"], lambdas["B"]),
        bd=transform_link(links.bd, lambdas["B"], lambdas["D"]),
    )


def bright_dark(links: DiamondLinks) -> tuple[float, float]:
    """Balanced recombiner probabilities.

    Amplitude to the bright port is (upper + lower)/2.  Since each path first
    receives a 1/sqrt(2) split and the recombiner contributes another 1/sqrt(2),
    this normalisation gives P_bright + P_dark = 1.
    """

    bright_amp = (links.upper() + links.lower()) / 2.0
    dark_amp = (links.upper() - links.lower()) / 2.0
    return float(abs(bright_amp) ** 2), float(abs(dark_amp) ** 2)


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<58s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<58s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def phase_difference(a: complex, b: complex) -> float:
    return float(np.angle(a * np.conj(b)))


def redistribute_flux_to_one_link(phi: float) -> DiamondLinks:
    """Gauge-equivalent representative with all loop phase on the upper SA link."""

    return DiamondLinks.from_angles(phi, 0.0, 0.0, 0.0)


def redistribute_flux_evenly(phi: float) -> DiamondLinks:
    """Gauge-equivalent representative with the same loop phase spread over links."""

    return DiamondLinks.from_angles(phi / 4.0, phi / 4.0, -phi / 4.0, -phi / 4.0)


def main() -> None:
    print("Aharonov-Bohm record machine: closed flux is hearable")
    print("=" * 78)

    # A deliberately uneven set of link phases.  The individual numbers have no
    # invariant meaning; only upper/lower relative phase, equivalently loop flux,
    # has a detector meaning.
    links = DiamondLinks.from_angles(theta_sa=0.31, theta_ad=-0.84, theta_sb=0.57, theta_bd=1.22)
    phi = links.flux()
    bright, dark = bright_dark(links)

    print("\n[1] Detector probabilities depend only on the closed loop flux")
    assert_close("P_bright", bright, (1.0 + math.cos(phi)) / 2.0)
    assert_close("P_dark", dark, (1.0 - math.cos(phi)) / 2.0)
    assert_close("P_bright + P_dark", bright + dark, 1.0)
    print(f"  loop flux Phi_AB = {phi:+.12f} rad")

    print("\n[2] Local gauge transformations move open phases but preserve flux and probabilities")
    lambdas = {"S": 0.73, "A": -1.10, "B": 2.40, "D": -0.22}
    transformed = gauge_transform(links, lambdas)
    bright_g, dark_g = bright_dark(transformed)

    # Both open paths pick up the same endpoint phase exp(i(lambda_S-lambda_D)).
    endpoint_factor = phase(lambdas["S"] - lambdas["D"])
    assert_complex_close("upper path covariance", transformed.upper(), endpoint_factor * links.upper())
    assert_complex_close("lower path covariance", transformed.lower(), endpoint_factor * links.lower())
    assert_complex_close("closed loop invariant", transformed.loop(), links.loop())
    assert_close("P_bright gauge invariant", bright_g, bright)
    assert_close("P_dark gauge invariant", dark_g, dark)
    print("  -> open path phases moved; the closed flux and detector probabilities did not.")

    print("\n[3] Same flux, different phase bookkeeping")
    one_link = redistribute_flux_to_one_link(phi)
    even = redistribute_flux_evenly(phi)
    for label, representative in [("all flux on one link", one_link), ("flux spread evenly", even)]:
        b, d = bright_dark(representative)
        assert_close(f"{label}: flux", representative.flux(), phi)
        assert_close(f"{label}: P_bright", b, bright)
        assert_close(f"{label}: P_dark", d, dark)
    print("  -> the distribution of link phases is gauge convention; the loop holonomy is the record.")

    print("\n[4] Open endpoint phase is not hearable after Gauss averaging")
    # Average an open S->D path phase over independent endpoint gauges.  A closed
    # loop has no endpoint and survives; an open path does not.
    grid = [2.0 * math.pi * k / 64 for k in range(64)]
    open_avg = sum(phase(ls - ld) * links.upper() for ls in grid for ld in grid) / (64 * 64)
    assert_close("|Gauss average of open path|", abs(open_avg), 0.0, tol=1e-12)
    assert_complex_close("Gauss average of closed loop", links.loop(), links.loop())
    print("  -> an open path phase is gauge bookkeeping, while the closed loop is gauge-readable.")

    print("\n[5] Relative phase equals the closed Wilson loop")
    rel = phase_difference(links.upper(), links.lower())
    assert_close("arg(upper lower*)", rel, phi)
    assert_complex_close("upper lower* = loop", links.upper() * np.conj(links.lower()), links.loop())

    print(
        """
Verdict
-------
The Aharonov-Bohm toy is the closed-loop version of the record grammar:

  individual endpoint/path phases are not physical records;
  local gauge transformations can move those phases around;
  the detector hears the gauge-invariant closed Wilson loop;
  the bright/dark interference probabilities depend only on the loop flux.

This is why phase is better described as a closed relational flux record than
as an endpoint property.  The next larger step would be to put this loop phase
on byte-level endpoints or to promote the graph toy to a small lattice plaquette.
"""
    )


if __name__ == "__main__":
    main()
