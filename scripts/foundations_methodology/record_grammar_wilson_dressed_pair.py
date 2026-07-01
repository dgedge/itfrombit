#!/usr/bin/env python3
"""Wilson-dressed endpoint records: the smallest gauge version of the record grammar.

The previous two toy theorems established the Bell-record language:

    endpoint records + relational stabilizer + detector readout + reset.

This script adds the first gauge layer.  It demonstrates the elementary rule
that charged endpoint records are not readable by themselves.  A detector can
hear only a gauge-invariant combination:

    endpoint relation + Wilson link/holonomy.

Toy model
---------
Use a one-particle two-endpoint Hilbert space spanned by |A>, |B>.  The local
U(1) gauge unitary is

    G(lambda_A, lambda_B) = diag(exp(i lambda_A), exp(i lambda_B)).

The bare endpoint coherence T_AB = |A><B| is not gauge-invariant.  A link
holonomy U_AB = exp(i theta_AB) transforms as

    U_AB -> exp(i(lambda_A - lambda_B)) U_AB,

so the dressed hopping/readout operator

    X_U = U_AB |A><B| + U_AB^* |B><A|

is covariant as a whole and has gauge-invariant expectation values when the
state and link are transformed together.

This is not a full gauge theory.  It is the smallest executable theorem for the
phrase "detectors hear dressed records, not bare endpoint phases".
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
A = np.array([1.0, 0.0], dtype=complex)
B = np.array([0.0, 1.0], dtype=complex)
T_AB = np.outer(A, B.conj())  # |A><B|
T_BA = np.outer(B, A.conj())  # |B><A|
X_BARE = T_AB + T_BA


def state(delta: float) -> np.ndarray:
    """(|A> + exp(i delta)|B>)/sqrt(2)."""

    return (A + np.exp(1j * delta) * B) / math.sqrt(2.0)


def gauge_unitary(lambda_a: float, lambda_b: float) -> np.ndarray:
    return np.diag([np.exp(1j * lambda_a), np.exp(1j * lambda_b)])


def link(theta: float) -> complex:
    return complex(np.exp(1j * theta))


def transform_link(u_ab: complex, lambda_a: float, lambda_b: float) -> complex:
    return complex(np.exp(1j * (lambda_a - lambda_b)) * u_ab)


def dressed_x(u_ab: complex) -> np.ndarray:
    return u_ab * T_AB + np.conj(u_ab) * T_BA


def expect(psi: np.ndarray, op: np.ndarray) -> float:
    return float(np.real_if_close(np.vdot(psi, op @ psi)))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<55s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_matrix_close(name: str, value: np.ndarray, target: np.ndarray, tol: float = 1e-12) -> None:
    err = float(np.linalg.norm(value - target))
    print(f"  {name:<55s} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def group_average_bare_phase(n: int = 64) -> complex:
    phases = []
    grid = [2.0 * math.pi * k / n for k in range(n)]
    for la in grid:
        for lb in grid:
            phases.append(np.exp(1j * (la - lb)))
    return complex(sum(phases) / len(phases))


def triangle_loop(theta_ab: float, theta_bc: float, theta_ca: float) -> complex:
    return link(theta_ab) * link(theta_bc) * link(theta_ca)


def transform_triangle(
    theta_ab: float,
    theta_bc: float,
    theta_ca: float,
    lambda_a: float,
    lambda_b: float,
    lambda_c: float,
) -> tuple[complex, complex, complex]:
    """Transform oriented links AB, BC, CA."""

    u_ab = np.exp(1j * (lambda_a - lambda_b)) * link(theta_ab)
    u_bc = np.exp(1j * (lambda_b - lambda_c)) * link(theta_bc)
    u_ca = np.exp(1j * (lambda_c - lambda_a)) * link(theta_ca)
    return complex(u_ab), complex(u_bc), complex(u_ca)


def main() -> None:
    print("Wilson-dressed endpoint records: executable gauge-readout theorem")
    print("=" * 78)

    delta = 0.37
    theta = -0.91
    lambda_a = 1.23
    lambda_b = -0.44
    psi = state(delta)
    u = link(theta)
    g = gauge_unitary(lambda_a, lambda_b)
    psi_g = g @ psi
    u_g = transform_link(u, lambda_a, lambda_b)

    print("\n[1] Bare endpoint coherence is gauge-dependent")
    bare = expect(psi, X_BARE)
    bare_g = expect(psi_g, X_BARE)
    assert_close("bare <X> = cos(delta)", bare, math.cos(delta))
    assert_close("bare transformed <X>", bare_g, math.cos(delta + lambda_b - lambda_a))
    if abs(bare - bare_g) < 1e-6:
        raise AssertionError("chosen gauge transformation accidentally left bare readout unchanged")
    print("  -> a bare relative endpoint phase is not a physical readout.")

    print("\n[2] Wilson-dressed readout is gauge-invariant")
    x_u = dressed_x(u)
    x_ug = dressed_x(u_g)
    # Covariance of the operator family: X(U') = G X(U) G^dag.
    assert_matrix_close("X(U') = G X(U) G^dag", x_ug, g @ x_u @ g.conj().T)
    dressed = expect(psi, x_u)
    dressed_g = expect(psi_g, x_ug)
    assert_close("dressed <X_U>", dressed, math.cos(delta + theta))
    assert_close("transformed dressed <X_U'>", dressed_g, dressed)
    print("  -> detector-visible endpoint coherence is endpoint relation plus holonomy.")

    print("\n[3] Gauss averaging kills bare open operators")
    avg = group_average_bare_phase()
    print(f"  group average of bare open phase exp(i(lambda_A-lambda_B)) = {avg.real:+.3e}{avg.imag:+.3e}i")
    if abs(avg) > 1e-12:
        raise AssertionError("bare open phase survived gauge average")
    # The dressed expectation is unchanged for a sample of gauge transformations.
    for la, lb in [(0.1, 0.2), (2.2, -0.7), (-1.1, 3.0), (5.0, 0.33)]:
        gg = gauge_unitary(la, lb)
        psi2 = gg @ psi
        u2 = transform_link(u, la, lb)
        assert_close(f"dressed invariant at ({la:.2f},{lb:.2f})", expect(psi2, dressed_x(u2)), dressed)
    print("  -> Gauss projection retains dressed open lines with endpoints, not naked endpoint coherence.")

    print("\n[4] Closed Wilson loops are hearable flux")
    th_ab, th_bc, th_ca = 0.2, -0.6, 1.4
    loop = triangle_loop(th_ab, th_bc, th_ca)
    u_ab, u_bc, u_ca = transform_triangle(th_ab, th_bc, th_ca, 0.9, -1.4, 2.1)
    loop_g = u_ab * u_bc * u_ca
    assert_close("arg loop", float(np.angle(loop)), float(np.angle(loop_g)))
    assert_close("|loop|", abs(loop), 1.0)
    assert_close("|loop transformed|", abs(loop_g), 1.0)
    print(f"  Wilson loop phase = {np.angle(loop):+.6f} rad, invariant under local endpoint gauges")
    print("  -> open link phases are gauge bookkeeping; closed holonomy/flux is directly hearable.")

    print(
        """
Verdict
-------
The first gauge layer of the record grammar is exact:

  bare endpoint relations are gauge-dependent and vanish under Gauss averaging;
  endpoint + Wilson holonomy is gauge-covariant as an operator family and has
  gauge-invariant detector expectation values;
  closed Wilson loops are invariant flux records.

This is the gauge version of "what can be heard": a detector cannot hear a
naked endpoint phase.  It hears either a dressed endpoint relation or a closed
loop/flux.  The next step, if desired, is to combine this with the byte-level
Bell cell so the logical endpoint records carry an explicit Wilson dressing.
"""
    )


if __name__ == "__main__":
    main()
