#!/usr/bin/env python3
"""Tiny Wilson-dressed scattering vertex.

This script connects the Wilson-link language to transition amplitudes.

Toy model
---------
Use a one-particle Hilbert space with two endpoint records |A>, |B>.  A local
U(1) gauge transformation acts as

    |A> -> exp(i lambda_A)|A>,
    |B> -> exp(i lambda_B)|B>.

A bare hopping vertex |B><A| is therefore not a gauge-covariant object by
itself.  The Wilson-dressed vertex is

    V_U = g ( U_BA |B><A| + U_BA^* |A><B| ),

where the oriented link U_BA transforms as

    U_BA -> exp(i(lambda_B - lambda_A)) U_BA.

This is the smallest scattering/transition example:

    amplitude(A -> B; t) = <B| exp(-i V_U t) |A>
                         = -i U_BA sin(g t).

The transition probability is sin^2(g t), while the amplitude carries the
Wilson phase.  With two possible weak vertices, amplitudes add before
probabilities, and the interference depends only on the closed relative Wilson
phase U_1 U_2^*.

Record-grammar reading
----------------------
A transition amplitude is not a detector fact yet.  It is the complex weight
assigned to an allowed endpoint rewrite.  Gauge grammar says that the rewrite
must be dressed by the link record.  A detector hears the final endpoint
projector, while interference hears only closed relative phase.

This is a toy graph calculation, not a full QFT scattering theory.
"""

from __future__ import annotations

import math

import numpy as np


A = np.array([1.0, 0.0], dtype=complex)
B = np.array([0.0, 1.0], dtype=complex)
T_BA = np.outer(B, A.conj())  # |B><A|, the A -> B rewrite.
T_AB = np.outer(A, B.conj())  # |A><B|, the reverse rewrite.


def phase(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def gauge_unitary(lambda_a: float, lambda_b: float) -> np.ndarray:
    return np.diag([np.exp(1j * lambda_a), np.exp(1j * lambda_b)])


def transform_link_ba(u_ba: complex, lambda_a: float, lambda_b: float) -> complex:
    """Transform the oriented link carried by |B><A|."""

    return complex(np.exp(1j * (lambda_b - lambda_a)) * u_ba)


def dressed_vertex(g: float, u_ba: complex) -> np.ndarray:
    """Hermitian two-endpoint Wilson-dressed hopping/scattering vertex."""

    return g * (u_ba * T_BA + np.conj(u_ba) * T_AB)


def expm_hermitian(h: np.ndarray, t: float) -> np.ndarray:
    """Small dense exp(-i H t) using the Hermitian eigensystem."""

    evals, evecs = np.linalg.eigh((h + h.conj().T) / 2.0)
    return evecs @ np.diag(np.exp(-1j * evals * t)) @ evecs.conj().T


def transition_amplitude(h: np.ndarray, t: float) -> complex:
    """<B| exp(-i H t) |A>."""

    return complex(np.vdot(B, expm_hermitian(h, t) @ A))


def born_two_route_amplitude(g1: float, u1: complex, g2: float, u2: complex) -> complex:
    """First-order weak transition amplitude through two alternative vertices."""

    return g1 * u1 + g2 * u2


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<66s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_complex_close(name: str, value: complex, target: complex, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<66s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<66s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def main() -> None:
    print("Tiny Wilson-dressed scattering vertex")
    print("=" * 88)

    g = 0.23
    t = 1.17
    theta = -0.64
    u = phase(theta)
    h = dressed_vertex(g, u)

    print("\n[1] The Wilson-dressed vertex gives the exact transition amplitude")
    amp = transition_amplitude(h, t)
    expected_amp = -1j * u * math.sin(g * t)
    assert_complex_close("<B|exp(-iVt)|A> = -i U sin(gt)", amp, expected_amp)
    assert_close("P(A -> B)", abs(amp) ** 2, math.sin(g * t) ** 2)
    print("  -> the amplitude carries the Wilson phase; the endpoint-click probability does not.")

    print("\n[2] The vertex is gauge-covariant only when the link transforms")
    lambda_a = 0.91
    lambda_b = -1.42
    gmat = gauge_unitary(lambda_a, lambda_b)
    u_g = transform_link_ba(u, lambda_a, lambda_b)
    h_g = dressed_vertex(g, u_g)
    assert_matrix_close("V(U') = G V(U) G^dag", h_g, gmat @ h @ gmat.conj().T)

    # Fixed-basis amplitudes transform covariantly.  If the external endpoint
    # states are transformed too, the physical matrix element is unchanged.
    amp_g_fixed_basis = transition_amplitude(h_g, t)
    assert_complex_close(
        "fixed-basis amplitude gains endpoint phase",
        amp_g_fixed_basis,
        phase(lambda_b - lambda_a) * amp,
    )
    amp_with_external_legs = complex(np.vdot(gmat @ B, expm_hermitian(h_g, t) @ (gmat @ A)))
    assert_complex_close("amplitude with transformed external legs is invariant", amp_with_external_legs, amp)
    print("  -> the Wilson link is the record that makes the endpoint rewrite gauge-consistent.")

    print("\n[3] A bare hopping vertex is not a gauge-covariant local rule")
    h_bare = dressed_vertex(g, 1.0 + 0j)
    bare_covariant_form = gmat @ h_bare @ gmat.conj().T
    mismatch = float(np.linalg.norm(bare_covariant_form - h_bare))
    print(f"  ||G V_bare G^dag - V_bare|| = {mismatch:.12g}")
    if mismatch < 1e-6:
        raise AssertionError("chosen gauge accidentally left the bare vertex unchanged")
    print("  -> a naked rewrite has no gauge-invariant meaning as a local service rule.")

    print("\n[4] Two weak vertices interfere through their closed relative Wilson phase")
    g1, g2 = 0.19, 0.31
    u1 = phase(0.25)
    u2 = phase(-1.13)
    two_route = born_two_route_amplitude(g1, u1, g2, u2)
    phi = float(np.angle(u1 * np.conj(u2)))
    expected_prob = g1 * g1 + g2 * g2 + 2.0 * g1 * g2 * math.cos(phi)
    assert_close("|g1 U1 + g2 U2|^2", abs(two_route) ** 2, expected_prob)

    u1_g = transform_link_ba(u1, lambda_a, lambda_b)
    u2_g = transform_link_ba(u2, lambda_a, lambda_b)
    two_route_g = born_two_route_amplitude(g1, u1_g, g2, u2_g)
    assert_complex_close("two-route amplitude covariant", two_route_g, phase(lambda_b - lambda_a) * two_route)
    assert_close("two-route probability gauge invariant", abs(two_route_g) ** 2, abs(two_route) ** 2)
    print(f"  relative closed phase arg(U1 U2*) = {phi:+.12f} rad")
    print("  -> interference hears the closed relative holonomy, not either path phase alone.")

    print("\n[5] Constructive and destructive limits are amplitude statements")
    constructive = born_two_route_amplitude(1.0, 1.0 + 0j, 1.0, 1.0 + 0j)
    destructive = born_two_route_amplitude(1.0, 1.0 + 0j, 1.0, -1.0 + 0j)
    assert_close("constructive |1+1|^2", abs(constructive) ** 2, 4.0)
    assert_close("destructive |1-1|^2", abs(destructive) ** 2, 0.0)
    print("  -> the complex transition weights add first; only then does a detector probability appear.")

    print(
        """
Verdict
-------
The tiny scattering vertex adds one new grammatical rule:

  an allowed endpoint rewrite carries a complex transition amplitude;
  a charged rewrite must be Wilson-dressed to be gauge-covariant;
  detector probabilities are squared amplitudes after the allowed rewrites add;
  interference depends only on closed relative Wilson phase.

This is still a toy.  It does not derive a Maxwell action, a continuum S-matrix,
or QED.  It only supplies the finite operator bridge from Wilson-link language
to transition amplitudes.
"""
    )


if __name__ == "__main__":
    main()
