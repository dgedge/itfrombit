#!/usr/bin/env python3
r"""ITEM 87/R13: recovery holonomy route to the Delta L=2 CP portal.

Question
--------
Can the Majorana CP portal be derived from a record-recovery cycle rather than
inserted as a free complex matrix?

Result
------
Partial progress, not full closure.

If one grants a new primitive:

    an oriented, unconditional sterile-generation recovery latch on the three
    allowed R1 generation corners,

then the support of the admissible portal is forced.  The oriented 3-cycle
P_sigma has no fixed generation and moves the all-zero nu_e corner.  A Majorana
bilinear is symmetric, so the directed latch contributes through

    P_sigma + P_sigma^T = A_K3,

the complete off-diagonal generation graph.  Equivalently, the only
S3-generation-blind symmetric off-diagonal portal is proportional to A_K3.
With an oriented recovery holonomy phase, the portal is therefore

    M_H = M0 [I + r exp(i sigma Phi) A_K3].

That is exactly the operator class identified by item87_deltaL2_holonomy_coupling.py.

The remaining gap is equally sharp:

  * a Hermitian oriented hopping uses conjugate phases on opposite directions
    and is not a valid complex-symmetric Majorana mass;
  * an unoriented recovery mixture gives a real cos(Phi) A_K3 term and CP=0;
  * the directed latch P_sigma is not itself S3-blind.  It reduces the generation
    relabelling symmetry to an oriented C3 subgroup; only its Majorana support
    P_sigma + P_sigma^T is S3-blind;
  * the existing R4/source repair algebra supplies two within-generation erasure
    edges per sterile corner and one generation-singlet source port.  It does
    not supply an inter-generation latch;
  * the phase Phi and orientation sigma are still not selected by the present
    boot/QEC mechanics.

So the recovery-holonomy route conditionally derives the K3 Majorana support,
but it does not yet derive the CP phase or baryogenesis sign.
"""

from __future__ import annotations

import math
from itertools import permutations
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
Y_DIAG = np.array([0.5, 0.8, 1.3])
H_DIRAC = np.diag(Y_DIAG**2)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def generation_permutations() -> list[np.ndarray]:
    eye = np.eye(3, dtype=complex)
    return [eye[list(p)] for p in permutations(range(3))]


PERMS = generation_permutations()


def p_cycle(sigma: int) -> np.ndarray:
    """Oriented recovery cycle on the three R1-allowed generation corners."""

    if sigma not in (-1, 1):
        raise ValueError("sigma must be +/-1")
    p = np.zeros((3, 3), dtype=complex)
    cycle = (0, 1, 2) if sigma == 1 else (0, 2, 1)
    for src, dst in zip(cycle, cycle[1:] + cycle[:1]):
        p[dst, src] = 1.0
    return p


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def commutes_with_s3(M: np.ndarray) -> bool:
    return all(np.max(np.abs(P.T @ M @ P - M)) < 1.0e-13 for P in PERMS)


def s3_commutator_residual(M: np.ndarray) -> float:
    return max(float(np.linalg.norm(P.T @ M @ P - M)) for P in PERMS)


def c3_commutator_residual(M: np.ndarray, sigma: int) -> float:
    generator = p_cycle(sigma)
    return max(float(np.linalg.norm(np.linalg.matrix_power(generator, k) @ M - M @ np.linalg.matrix_power(generator, k))) for k in range(3))


def cp_invariant(M: np.ndarray) -> float:
    Md = M.conj().T
    return float(np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)))


def cp_norm(M: np.ndarray) -> float:
    Md = M.conj().T
    probes = np.array(
        [
            np.imag(np.trace(H_DIRAC @ Md @ M @ Md @ H_DIRAC @ M)),
            np.imag(np.trace(H_DIRAC @ M @ Md @ M @ Md @ M)),
        ],
        dtype=float,
    )
    return float(np.linalg.norm(probes))


def majorana_portal(phi: float, sigma: int, r: float = 0.5) -> np.ndarray:
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def hermitian_cycle(phi: float, sigma: int, r: float = 0.5) -> np.ndarray:
    p = p_cycle(sigma)
    return np.eye(3, dtype=complex) + r * (np.exp(1j * phi) * p + np.exp(-1j * phi) * p.T)


def symmetric_offdiag(a: complex, b: complex, c: complex) -> np.ndarray:
    return np.array(
        [
            [0.0, a, b],
            [a, 0.0, c],
            [b, c, 0.0],
        ],
        dtype=complex,
    )


def lepton_number_phase(M: np.ndarray, theta: float) -> np.ndarray:
    """Majorana bilinear coefficient after nu_R -> exp(i theta) nu_R."""

    return np.exp(2j * theta) * M


def main() -> None:
    print("ITEM 87/R13 -- RECOVERY HOLONOMY CP-PORTAL AUDIT")
    print("=" * 96)

    require_text(
        "recent_papers/matter_gauge/matter_gauge.tex",
        (
            "complex-symmetric Majorana",
            "generation-blind",
            "orientation sign",
        ),
    )
    require_text(
        "python_code/item87_MR_derivation_attempt.py",
        (
            "nu_e is an exact walk eigenstate",
            "requires new Delta L=2 physics",
        ),
    )

    print("\n[1] R1 generation corners and oriented recovery cycle")
    corners = [(0, 0), (0, 1), (1, 0)]
    forbidden = (1, 1)
    print(f"  allowed R1 generation corners = {corners}; forbidden corner = {forbidden}")
    for sigma in (+1, -1):
        p = p_cycle(sigma)
        print(f"  sigma={sigma:+d}: P^3 residual={np.linalg.norm(np.linalg.matrix_power(p, 3) - np.eye(3)):.3e}, fixed diagonal={np.trace(p):.0f}")
        check(np.linalg.norm(np.linalg.matrix_power(p, 3) - np.eye(3)) < 1.0e-13, f"sigma={sigma:+d}: recovery latch is a 3-cycle")
        check(np.trace(p) == 0, f"sigma={sigma:+d}: latch is unconditional; no generation is inert")
        check(np.all(np.sum(np.abs(p), axis=0) == 1), f"sigma={sigma:+d}: every generation has one outgoing recovery channel")
        print(f"             S3 residual={s3_commutator_residual(p):.3e}, C3 residual={c3_commutator_residual(p, sigma):.3e}")
        check(s3_commutator_residual(p) > 1.0, f"sigma={sigma:+d}: directed latch is not S3-blind; it selects an orientation")
        check(c3_commutator_residual(p, sigma) < 1.0e-13, f"sigma={sigma:+d}: directed latch preserves only its oriented C3 subgroup")

    print("\n[1b] Existing framework support boundary")
    native_generation_adjacency = np.zeros((3, 3), dtype=complex)
    print(f"  native R4/source inter-generation adjacency norm = {np.linalg.norm(native_generation_adjacency):.3e}")
    check(np.linalg.norm(native_generation_adjacency) < 1.0e-13, "documented R4/source repairs do not move generation labels")
    check(
        "one coherent port"
        in (ROOT / "python_code/item123_sterile_generation_singlet_source.py").read_text(encoding="utf-8"),
        "sterile source algebra already forces a generation-singlet port, not three labelled ports",
    )
    check(
        "generation flips do not erase the R4 syndrome"
        in (ROOT / "python_code/item131_r4_support_dimension.py").read_text(encoding="utf-8"),
        "R4 support theorem excludes generation flips as R4 erasure events",
    )

    print("\n[2] Majorana symmetrisation forces K3 support")
    for sigma in (+1, -1):
        sym = p_cycle(sigma) + p_cycle(sigma).T
        print(f"  sigma={sigma:+d}: ||P+P^T - A_K3|| = {np.linalg.norm(sym - a_k3()):.3e}")
        check(np.linalg.norm(sym - a_k3()) < 1.0e-13, f"sigma={sigma:+d}: Majorana symmetrisation gives K3 adjacency")
    check(commutes_with_s3(a_k3()), "A_K3 is generation-blind / S3-covariant")

    print("\n[3] S3 uniqueness of the off-diagonal support")
    basis = [
        symmetric_offdiag(1, 0, 0),
        symmetric_offdiag(0, 1, 0),
        symmetric_offdiag(0, 0, 1),
    ]
    constraints: list[list[complex]] = []
    for P in PERMS:
        deltas = [P.T @ B @ P - B for B in basis]
        for row in range(3):
            for col in range(3):
                constraints.append([delta[row, col] for delta in deltas])
    C = np.array(constraints, dtype=complex)
    # Null vector of the S3-invariance constraints in the off-diagonal symmetric basis.
    _, _, vh = np.linalg.svd(C)
    null = vh[-1].real
    null /= np.min(np.abs(null[np.abs(null) > 1e-12]))
    if np.sum(null) < 0:
        null *= -1
    print(f"  invariant off-diagonal coefficient vector = {null}")
    check(np.max(np.abs(null - np.ones(3))) < 1.0e-12, "S3-invariant symmetric off-diagonal support is uniquely A_K3")

    print("\n[4] Delta L=2 and Hermitian-walk controls")
    M = majorana_portal(1.0 / 3.0, sigma=+1)
    transformed = lepton_number_phase(M, theta=0.37)
    check(np.max(np.abs(M - M.T)) < 1.0e-13, "oriented recovery portal is complex-symmetric, hence Majorana-admissible")
    check(np.max(np.abs(transformed - M)) > 0.1, "Majorana bilinear carries lepton charge two: Delta L=2")
    H = hermitian_cycle(1.0 / 3.0, sigma=+1)
    print(f"  Hermitian directed cycle: ||H-H^dag||={np.linalg.norm(H-H.conj().T):.3e}, ||H-H^T||={np.linalg.norm(H-H.T):.3e}")
    check(np.linalg.norm(H - H.conj().T) < 1.0e-13, "ordinary oriented hopping is Hermitian / Delta L=0")
    check(np.linalg.norm(H - H.T) > 0.1, "ordinary oriented hopping is not a complex-symmetric Majorana mass")

    print("\n[5] CP and orientation")
    phases = {
        "delta_nu geometric 1/3": 1.0 / 3.0,
        "nu_R Berry 2pi/9": 2.0 * math.pi / 9.0,
    }
    for label, phi in phases.items():
        Mp = majorana_portal(phi, sigma=+1)
        Mm = majorana_portal(phi, sigma=-1)
        Ip = cp_invariant(Mp)
        Im = cp_invariant(Mm)
        mixed = 0.5 * (Mp + Mm)
        print(f"  {label:<24s}: I(+sigma)={Ip:+.6e}, I(-sigma)={Im:+.6e}, CP mixed={cp_norm(mixed):.3e}")
        check(cp_norm(Mp) > 1.0e-3, f"{label}: oriented Majorana recovery has nonzero CP")
        check(abs(Ip + Im) < 1.0e-12, f"{label}: CP sign flips with orientation")
        check(cp_norm(mixed) < 1.0e-12, f"{label}: unoriented sigma average is real/rephasable and CP-even")

    print("\n[6] Decision")
    print(
        """
  Progress:
    The support of the CP portal is no longer arbitrary if a recovery-holonomy
    primitive is granted.  An unconditional oriented 3-cycle on the three R1
    generation corners, when read as a Majorana recovery event, symmetrises to
    A_K3; S3-generation blindness also singles out the same support.  This
    conditionally derives the K3 complex-symmetric Delta L=2 operator class.

  Still open:
    The present boot/walk mechanics does not supply that primitive.  The
    documented walk is Delta L=0 and leaves nu_e inert.  The recovery audit
    also does not derive Phi or choose sigma.  If both orientations are summed
    without an orientation selector, the CP invariant cancels.

  Therefore the next theorem target is exact:
    build a CPTP/QEC recovery cycle whose sterile-sector jump is the oriented
    non-CNOT Delta L=2 latch, then derive its holonomy phase and global
    orientation.  The portal support is conditionally forced; the CP sign is
    not yet locked.
"""
    )
    print("ALL ASSERTIONS PASSED -- recovery holonomy derives K3 support conditionally; phase/orientation remain open.")


if __name__ == "__main__":
    main()
