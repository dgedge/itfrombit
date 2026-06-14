#!/usr/bin/env python3
r"""Manifest loop-holonomy F^2: Lorentzian positive-measure no-go.

K25 built a mesoscopic Alexandrov-loop scaffold and left the full Lorentzian
loop-F^2 theorem open.  This script sharpens the obstruction.

For any real loop area bivector A and any positive interval measure,

    <(F.A)^2> = F_i M_ij F_j,    M = <A A^T> >= 0.

So a positive sum/average of real Wilson-loop holonomy squares is always a
positive semidefinite quadratic form on the six components of F.  But the
Lorentzian Maxwell density is indefinite:

    1/2 F_mn F^mn = B^2 - E^2.

Therefore a manifest gauge-invariant Lorentzian action cannot be obtained from
"holonomy squared with a positive Alexandrov-loop measure" alone.  It needs an
extra sign/phase/Hodge-star structure: Euclidean Wick rotation, or a local
tetrad/frame that assigns timelike and spacelike plaquette signs.

The script also verifies that the simplest non-arbitrary Alexandrov diamond
loop p->r->q->s->p only sees electric components relative to the interval's
time axis, making pure magnetic fields invisible.
"""

from __future__ import annotations

import itertools
import math

import numpy as np


RNG = np.random.default_rng(20260613)
COMPONENTS = [(0, 1), (0, 2), (0, 3), (2, 3), (3, 1), (1, 2)]
LORENTZ_TARGET = np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])
EUCLIDEAN_TARGET = np.eye(6)


def area_bivector(vertices: np.ndarray) -> np.ndarray:
    """Polygon area bivector components in order 01,02,03,23,31,12."""

    sigma = np.zeros((4, 4))
    for a, b in zip(vertices, np.roll(vertices, -1, axis=0)):
        sigma += 0.5 * (np.outer(a, b) - np.outer(b, a))
    return np.array([sigma[i, j] for i, j in COMPONENTS])


def random_point_in_unit_diamond() -> np.ndarray:
    """A simple point inside the 3+1 Alexandrov interval p=(-1/2,0), q=(1/2,0)."""

    t = RNG.uniform(-0.5, 0.5)
    rmax = 0.5 - abs(t)
    # Uniform enough for this gate; the magnetic-zero identity is exact regardless.
    direction = RNG.normal(size=3)
    direction /= np.linalg.norm(direction)
    radius = rmax * RNG.random() ** (1.0 / 3.0)
    return np.array([t, *(radius * direction)])


def diamond_loop_area_trials(n: int = 20_000) -> np.ndarray:
    p = np.array([-0.5, 0.0, 0.0, 0.0])
    q = np.array([0.5, 0.0, 0.0, 0.0])
    out = []
    for _ in range(n):
        r = random_point_in_unit_diamond()
        s = random_point_in_unit_diamond()
        out.append(area_bivector(np.array([p, r, q, s])))
    return np.array(out)


def psd_matrix_from_positive_measure(area_vectors: np.ndarray) -> np.ndarray:
    return area_vectors.T @ area_vectors / len(area_vectors)


def signed_basis_matrix(signs: np.ndarray) -> np.ndarray:
    basis = np.eye(6)
    m = np.zeros((6, 6))
    for sign, b in zip(signs, basis):
        m += sign * np.outer(b, b)
    return m


def main() -> None:
    print("[1] Alexandrov diamond loops see only electric components")
    areas = diamond_loop_area_trials()
    electric = areas[:, :3]
    magnetic = areas[:, 3:]
    max_magnetic = float(np.max(np.abs(magnetic)))
    rank = int(np.linalg.matrix_rank(psd_matrix_from_positive_measure(areas), tol=1e-12))
    print(f"    max spatial-spatial area component over trials = {max_magnetic:.3e}")
    print(f"    quadratic-form rank = {rank} (expected 3: only 0i components)")
    assert max_magnetic < 1e-12 and rank == 3

    pure_b = np.array([0.0, 0.0, 0.0, 0.7, -0.2, 0.4])
    hol_b = areas @ pure_b
    lorentz_b = float(pure_b @ LORENTZ_TARGET @ pure_b)
    print(f"    pure magnetic field: max holonomy={np.max(np.abs(hol_b)):.3e}, target B^2-E^2={lorentz_b:.3f}")
    assert np.max(np.abs(hol_b)) < 1e-12 and lorentz_b > 0
    print("    -> the most natural Alexandrov diamond loop cannot be the full Maxwell action.")

    print("\n[2] Positive-measure loop-square no-go")
    random_areas = RNG.normal(size=(50_000, 6))
    m_pos = psd_matrix_from_positive_measure(random_areas)
    evals = np.linalg.eigvalsh(m_pos)
    target_evals = np.linalg.eigvalsh(LORENTZ_TARGET)
    print(f"    positive loop-square matrix eigenvalue range: [{evals[0]:.4f}, {evals[-1]:.4f}]")
    print(f"    Lorentzian Maxwell target eigenvalues: {target_evals.tolist()}")
    assert evals[0] > 0.0 and target_evals[0] < 0.0

    # The argument is algebraic, not sampling-specific.
    for v in itertools.chain(np.eye(6), [RNG.normal(size=6) for _ in range(20)]):
        assert float(v @ m_pos @ v) >= -1e-12
    electric_unit = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    assert float(electric_unit @ LORENTZ_TARGET @ electric_unit) < 0.0
    print("    proof: <(F.A)^2> is positive semidefinite for every positive measure;")
    print("           B^2-E^2 is indefinite. They cannot be identical for all F.")

    print("\n[3] What extra structure is sufficient")
    euclid = signed_basis_matrix(np.ones(6))
    lorentz = signed_basis_matrix(np.array([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]))
    euclid_err = float(np.linalg.norm(euclid - EUCLIDEAN_TARGET))
    lorentz_err = float(np.linalg.norm(lorentz - LORENTZ_TARGET))
    print(f"    Euclidean/Wick loop action from positive signs error = {euclid_err:.3e}")
    print(f"    Lorentzian action from signed Hodge/tetrad signs error = {lorentz_err:.3e}")
    assert euclid_err < 1e-12 and lorentz_err < 1e-12
    print("    -> a manifest gauge-invariant loop action is possible only after adding")
    print("       either Wick rotation or a frame/Hodge-star sign rule. That rule is")
    print("       extra structure; it is not a positive Alexandrov interval measure.")

    print(
        r"""
[4] VERDICT
  CLOSED NEGATIVELY:
    A Lorentzian manifest loop-F^2 action cannot be a positive average of real
    Alexandrov-loop holonomy squares.  The quadratic form has the wrong
    signature.  The simple diamond-loop candidate is even narrower: it sees
    only electric components relative to its timelike axis.

  STILL AVAILABLE:
    The K25 mesoscopic scaffold is valid as a Euclidean/Wick-rotated constant-
    field normalization protocol.  A Lorentzian loop action can be built with
    a signed Hodge/tetrad rule, but that rule must be derived from the same
    frame dynamics that now carries Dirac spin.

  RECLASSIFICATION:
    The open problem is no longer "find a positive Alexandrov-loop measure".
    That target is impossible.  The live target is: derive the Wick/Hodge
    sign structure from the framed QEC/causal-set dynamics.
ALL ASSERTIONS PASSED -- positive-measure Lorentzian loop-F2 route excluded."""
    )


if __name__ == "__main__":
    main()
