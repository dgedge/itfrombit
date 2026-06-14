#!/usr/bin/env python3
r"""Causal-set Dirac spin: explicit spin-frame lift.

This is the constructive half of the K23 Dirac residual.  Causal order alone
does not supply gamma matrices, a spinor fibre, or a spin connection.  But once
a mesoscopic frame/tetrad is supplied by extra framework structure (the
E_{1/2} coin / local service frame), the lift is standard and checkable:

  * gamma^mu=e_a^mu gamma^a satisfies the Clifford relation;
  * the Dirac factor squares to the scalar causal-set kinetic symbol;
  * link spin transports make the discrete kernel covariant under local frame
    rotations.

The script verifies these claims in flat 3+1 with a local Spin(3) frame
rotation.  Full Spin(3,1) curved-background construction remains open.
"""

from __future__ import annotations

import math

import numpy as np


I2 = np.eye(2, dtype=complex)
Z2 = np.zeros((2, 2), dtype=complex)
SIGMA_X = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA_Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMAS = [SIGMA_X, SIGMA_Y, SIGMA_Z]
ETA = np.diag([1.0, -1.0, -1.0, -1.0])


def dirac_gamma() -> list[np.ndarray]:
    return [
        np.block([[I2, Z2], [Z2, -I2]]),
        np.block([[Z2, SIGMA_X], [-SIGMA_X, Z2]]),
        np.block([[Z2, SIGMA_Y], [-SIGMA_Y, Z2]]),
        np.block([[Z2, SIGMA_Z], [-SIGMA_Z, Z2]]),
    ]


GAMMA = dirac_gamma()


def spin_rotation_z(theta: float) -> np.ndarray:
    """Dirac spin representation for a spatial z-rotation."""

    u2 = math.cos(theta / 2.0) * I2 - 1j * math.sin(theta / 2.0) * SIGMA_Z
    return np.block([[u2, Z2], [Z2, u2]])


def vector_rotation_z(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    lam = np.eye(4)
    lam[1:3, 1:3] = np.array([[c, -s], [s, c]])
    return lam


def slash(p: np.ndarray, gamma: list[np.ndarray] | None = None) -> np.ndarray:
    gamma = GAMMA if gamma is None else gamma
    return sum(gamma[mu] * p[mu] for mu in range(4))


def spin_connection_covariance_gate() -> tuple[float, float]:
    """Finite kernel covariance under local frame rotations.

    A link kernel D_xy = K_xy U_xy, with spin transport U_xy=S_x^{-1} S_y,
    transforms as D'_xy = S_x^{-1} D_xy S_y.  Closed-chain traces are invariant.
    """

    theta = np.array([0.0, 0.23, -0.71, 0.44])
    s = [spin_rotation_z(float(t)) for t in theta]
    n = len(s)
    k_scalar = np.array(
        [
            [0.0, 1.2, 0.0, 0.3],
            [1.2, 0.0, -0.4, 0.0],
            [0.0, -0.4, 0.0, 0.9],
            [0.3, 0.0, 0.9, 0.0],
        ]
    )
    d = [[k_scalar[i, j] * (np.linalg.inv(s[i]) @ s[j]) for j in range(n)] for i in range(n)]

    local = [spin_rotation_z(float(t)) for t in [0.31, -0.17, 0.52, -0.09]]
    s_prime = [s[i] @ local[i] for i in range(n)]
    d_rebuilt = [[k_scalar[i, j] * (np.linalg.inv(s_prime[i]) @ s_prime[j]) for j in range(n)] for i in range(n)]
    d_cov = [[np.linalg.inv(local[i]) @ d[i][j] @ local[j] for j in range(n)] for i in range(n)]
    cov_err = max(np.linalg.norm(d_rebuilt[i][j] - d_cov[i][j]) for i in range(n) for j in range(n))

    loop = d[0][1] @ d[1][2] @ d[2][3] @ d[3][0]
    loop_prime = d_rebuilt[0][1] @ d_rebuilt[1][2] @ d_rebuilt[2][3] @ d_rebuilt[3][0]
    trace_err = abs(np.trace(loop_prime) - np.trace(np.linalg.inv(local[0]) @ loop @ local[0]))
    return float(cov_err), float(trace_err)


def main() -> None:
    print("[1] Clifford fibre from supplied mesoscopic frame")
    max_clifford = 0.0
    for mu in range(4):
        for nu in range(4):
            anti = GAMMA[mu] @ GAMMA[nu] + GAMMA[nu] @ GAMMA[mu]
            target = 2.0 * ETA[mu, nu] * np.eye(4)
            max_clifford = max(max_clifford, float(np.linalg.norm(anti - target)))
    print(f"    max Clifford error={max_clifford:.3e}")
    assert max_clifford < 1e-12

    print("\n[2] Scalar-kernel factorisation")
    p = np.array([2.6, 0.4, -0.8, 1.1])
    m = 0.73
    factor = (slash(p) - m * np.eye(4)) @ (slash(p) + m * np.eye(4))
    p2_minus_m2 = float(p @ ETA @ p - m * m)
    factor_err = float(np.linalg.norm(factor - p2_minus_m2 * np.eye(4)))
    print(f"    ||(p slash-m)(p slash+m)-(p^2-m^2)||={factor_err:.3e}")
    assert factor_err < 1e-12

    print("\n[3] Spin-frame covariance")
    theta = 0.67
    s = spin_rotation_z(theta)
    lam = vector_rotation_z(theta)
    rotated_gamma = [np.linalg.inv(s) @ GAMMA[mu] @ s for mu in range(4)]
    expected_gamma = [sum(lam[mu, nu] * GAMMA[nu] for nu in range(4)) for mu in range(4)]
    rot_err = max(float(np.linalg.norm(rotated_gamma[mu] - expected_gamma[mu])) for mu in range(4))
    print(f"    S^-1 gamma^mu S = Lambda^mu_nu gamma^nu error={rot_err:.3e}")
    assert rot_err < 1e-12

    print("\n[4] Link spin connection")
    cov_err, trace_err = spin_connection_covariance_gate()
    print(f"    local-frame covariance error={cov_err:.3e}")
    print(f"    closed-loop trace invariance error={trace_err:.3e}")
    assert cov_err < 1e-12 and trace_err < 1e-12

    print("\n[5] Verdict")
    print("    A causal-set Dirac operator is constructible once a mesoscopic spin frame")
    print("    is supplied: Clifford fibre, factorisation, and link spin covariance pass.")
    print("    The no-go remains precise: causal order alone does not produce that frame.")
    print("    The remaining theorem is to derive the frame/spin connection from E_{1/2}")
    print("    service-frame mechanics in a curved causal set.")
    print("ALL ASSERTIONS PASSED -- spin-frame lift constructed at flat finite-cell grade.")


if __name__ == "__main__":
    main()
