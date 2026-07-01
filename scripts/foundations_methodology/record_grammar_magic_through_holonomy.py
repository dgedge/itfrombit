#!/usr/bin/env python3
r"""Magic through a Wilson loop: non-Clifford phase as closed holonomy.

This script combines two guardrails:

  * Magic-state warning: a T phase is outside stabilizer/Clifford bookkeeping.
  * Wilson-loop grammar: a physical phase must be gauge-readable, not naked.

Toy model
---------
Use the same two-route diamond as the Aharonov-Bohm scripts:

        A
      /   \
    S       D
      \   /
        B

The route qubit is the two-dimensional space spanned by the upper and lower
routes.  The physical relative phase is the closed Wilson loop

    W = U_upper U_lower^* = exp(i Phi).

Set Phi = pi/4.  In a gauge-fixed route basis, the state is

    |psi_W> = (|upper> + exp(i Phi)|lower>)/sqrt(2),

which is the one-qubit T magic state when Phi = pi/4.

Record-grammar reading
----------------------
The non-Clifford phase is not stored as an endpoint property or open-link
phase.  A naked open phase is gauge bookkeeping and vanishes under Gauss
averaging.  The closed Wilson loop is gauge-invariant and can carry the magic
phase as detector-readable content.

This is a finite U(1) toy, not a theory of magic-state distillation or a
derivation of gauge dynamics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def ket(index: int) -> np.ndarray:
    v = np.zeros(2, dtype=complex)
    v[index] = 1.0
    return v


def density(psi: np.ndarray) -> np.ndarray:
    return np.outer(psi, psi.conj())


def route_state(phi: float) -> np.ndarray:
    """Gauge-fixed route qubit (|upper> + exp(i phi)|lower>)/sqrt(2)."""

    return (ket(0) + phase(phi) * ket(1)) / math.sqrt(2.0)


def projector(psi: np.ndarray) -> np.ndarray:
    psi = psi / np.linalg.norm(psi)
    return density(psi)


def expect(rho: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.trace(rho @ op)))


def bloch(rho: np.ndarray) -> np.ndarray:
    return np.array([expect(rho, X), expect(rho, Y), expect(rho, Z)])


@dataclass(frozen=True)
class DiamondLinks:
    """Four oriented U(1) link phases for the diamond interferometer."""

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


def route_detector_probability(phi_state: float, phi_detector: float) -> float:
    """Probability that |psi(phi_state)> clicks in projector |psi(phi_detector)>."""

    psi = route_state(phi_state)
    p_det = projector(route_state(phi_detector))
    return float(np.real_if_close(np.vdot(psi, p_det @ psi)))


def stabilizer_states() -> dict[str, float]:
    """Route-state phases for the equatorial Pauli stabilizer axes."""

    return {
        "|+>": 0.0,
        "|+i>": math.pi / 2.0,
        "|->": math.pi,
        "|-i>": -math.pi / 2.0,
    }


def group_average_open_link(links: DiamondLinks, n: int = 64) -> complex:
    """Average a naked S->A link over its endpoint gauges."""

    grid = [2.0 * math.pi * k / n for k in range(n)]
    total = 0.0j
    for ls in grid:
        for la in grid:
            total += transform_link(links.sa, ls, la)
    return complex(total / (n * n))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Magic through a Wilson loop: T phase as closed holonomy")
    print("=" * 92)

    phi_t = math.pi / 4.0
    links = DiamondLinks.from_angles(theta_sa=phi_t, theta_ad=0.0, theta_sb=0.0, theta_bd=0.0)

    print("\n[1] The closed loop carries the T phase")
    assert_close("arg W", links.flux(), phi_t)
    assert_complex_close("W = exp(i pi/4)", links.loop(), phase(phi_t))
    print("  -> the non-Clifford phase is attached to a closed Wilson loop.")

    print("\n[2] The gauge-fixed route state is a magic state")
    psi = route_state(links.flux())
    rho = density(psi)
    r = bloch(rho)
    assert_close("<X>", r[0], 1.0 / math.sqrt(2.0))
    assert_close("<Y>", r[1], 1.0 / math.sqrt(2.0))
    assert_close("<Z>", r[2], 0.0)
    l1 = float(np.sum(np.abs(r)))
    assert_close("|x|+|y|+|z|", l1, math.sqrt(2.0))
    assert_close("stabilizer-octahedron violation", l1 - 1.0, math.sqrt(2.0) - 1.0)
    print("  -> the Wilson-loop phase is magic: outside the one-qubit stabilizer polytope.")

    print("\n[3] A naked open-link phase is not the magic record")
    lambdas = {"S": 0.61, "A": -1.72, "B": 2.31, "D": -0.44}
    transformed = gauge_transform(links, lambdas)
    if abs(np.angle(transformed.sa) - np.angle(links.sa)) < 1e-6:
        raise AssertionError("chosen gauge accidentally left naked S-A link unchanged")
    assert_complex_close("closed loop invariant", transformed.loop(), links.loop())
    open_avg = group_average_open_link(links)
    assert_close("|Gauss average of naked S-A link|", abs(open_avg), 0.0, tol=1e-12)
    print("  -> the open-link phase moves under gauge and averages away; W does not.")

    print("\n[4] An aligned Wilson-loop detector reads the magic phase exactly")
    p_aligned = route_detector_probability(phi_t, phi_t)
    p_plus = route_detector_probability(phi_t, 0.0)
    p_plusi = route_detector_probability(phi_t, math.pi / 2.0)
    expected_stab = (1.0 + 1.0 / math.sqrt(2.0)) / 2.0
    assert_close("Pr(T-loop detector | T-loop state)", p_aligned, 1.0)
    assert_close("Pr(|+> stabilizer detector | T-loop state)", p_plus, expected_stab)
    assert_close("Pr(|+i> stabilizer detector | T-loop state)", p_plusi, expected_stab)
    print("  -> stabilizer-axis readouts hear shadows; the loop-aligned detector hears the full phase.")

    print("\n[5] Best stabilizer replacement still loses detector-readable signal")
    best_label = ""
    best_prob = -1.0
    for label, stab_phi in stabilizer_states().items():
        prob = route_detector_probability(stab_phi, phi_t)
        print(f"  T-loop detector probability for {label:<4s}: {prob:.12g}")
        if prob > best_prob:
            best_label = label
            best_prob = prob
    assert_close("best stabilizer replacement", best_prob, expected_stab)
    print(f"  nearest stabilizer route state is {best_label}, losing {1.0 - best_prob:.6f}")

    print("\n[6] Same loop magic after arbitrary link redistribution")
    redistributed = DiamondLinks.from_angles(phi_t / 4.0, phi_t / 4.0, -phi_t / 4.0, -phi_t / 4.0)
    assert_complex_close("redistributed W = original W", redistributed.loop(), links.loop())
    assert_close("redistributed T detector probability", route_detector_probability(redistributed.flux(), phi_t), 1.0)
    print("  -> where the link phases sit is gauge bookkeeping; the loop holonomy carries the resource.")

    print(
        """
Verdict
-------
The magic-through-holonomy toy keeps the two guardrails together:

  a T phase is non-Clifford magic and cannot be replaced by stabilizers;
  a naked open-link phase is gauge bookkeeping and cannot be the physical record;
  the closed Wilson loop can carry the T phase as gauge-readable content;
  a loop-aligned detector reads it exactly, while stabilizer readouts lose 14.6%.

Record-grammar rule:

  magic can be retained as a holonomy resource,
  but only when the phase is attached to gauge-readable loop data.

Boundary
--------
This is not a magic-state distillation protocol or a gauge-dynamics derivation.
It is the smallest finite example showing that non-Clifford content can live in
a closed Wilson record rather than in a naked endpoint phase.
"""
    )


if __name__ == "__main__":
    main()
