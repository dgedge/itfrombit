#!/usr/bin/env python3
r"""ITEM 87/R15: can the CP portal latch, phase, and orientation be closed?

Question
--------
The recovery-holonomy audit reduced the CP/baryogenesis-sign frontier to three
missing objects:

    1. an oriented unconditional sterile-generation latch P_sigma,
    2. its holonomy phase Phi,
    3. the absolute orientation sigma.

This script tests whether those three can be closed from the current canon
without adding a new primitive.

Result
------
No full closure from current canon.

What can be closed conditionally:

    If a new oriented sterile recovery latch is granted, its directed part is

        P_sigma = (A_K3 + sigma K_or)/2,

    where A_K3 is the already-forced S3-blind Majorana support and K_or is an
    antisymmetric generation-orientation tensor.  Thus the missing latch is
    exactly the missing sign-representation object K_or.

    If, further, the latch is required to carry a faithful one-dimensional
    character of its C3 recovery cycle, then the native character phase is

        Phi = 2*pi/3

    up to orientation conjugation.  This is a mathematically clean conditional
    phase closure, and it gives nonzero weak-basis CP with sign reversed by
    sigma.

What cannot be closed:

    Current R4/source repair algebra is S3-blind and within-generation.  It has
    no inter-generation adjacency and does not read the R1 sign-representation
    tensor K_or as a service channel.  The R1-boundary readout theorem sharpens
    this to an iff condition: K_or is readable only if the recovery environment
    carries a sign-representation orientation pointer.  The documented walk has
    a global spatial phase orientation, but axis-booking canon explicitly
    demotes generation=spatial-axis identification.  Therefore there is no
    current-canon bridge from global spatial handedness to the sterile
    generation latch orientation.

So the honest status is:

    * support A_K3: conditionally forced by the previous audit;
    * phase Phi: conditionally closable as 2*pi/3 only after adding the faithful
      C3-character premise for a new latch;
    * service-readable directed latch and absolute sigma: not closed by current
      canon.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
import math

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


def generation_permutations() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3, dtype=complex)
    for p in permutations(range(3)):
        inv_count = sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j])
        parity = -1 if inv_count % 2 else 1
        out.append((p, eye[list(p)], parity))
    return out


PERMS = generation_permutations()


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def p_cycle(sigma: int) -> np.ndarray:
    if sigma not in (-1, 1):
        raise ValueError("sigma must be +/-1")
    p = np.zeros((3, 3), dtype=complex)
    cycle = (0, 1, 2) if sigma == 1 else (0, 2, 1)
    for src, dst in zip(cycle, cycle[1:] + cycle[:1]):
        p[dst, src] = 1.0
    return p


def k_orientation() -> np.ndarray:
    return p_cycle(+1) - p_cycle(+1).T


def commutant_dimension() -> int:
    rows = []
    for _perm, p, _parity in PERMS:
        rows.append(np.kron(p.T, np.eye(3)) - np.kron(np.eye(3), p))
    mat = np.vstack(rows)
    svals = np.linalg.svd(mat, compute_uv=False)
    return int(np.sum(svals < 1.0e-12))


def projection_residual_to_s3_commutant(M: np.ndarray) -> float:
    """Project M onto span{I,J}, the S3 commutant on three generations."""

    basis = [np.eye(3, dtype=complex), np.ones((3, 3), dtype=complex)]
    B = np.stack([b.reshape(-1) for b in basis], axis=1)
    coeffs, *_ = np.linalg.lstsq(B, M.reshape(-1), rcond=None)
    proj = sum(c * b for c, b in zip(coeffs, basis))
    return float(np.linalg.norm(M - proj))


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


def main() -> None:
    print("ITEM 87/R15 -- DIRECTED LATCH / PHASE / ORIENTATION CLOSURE AUDIT")
    print("=" * 96)

    require_text(
        "python_code/item123_sterile_generation_singlet_source.py",
        (
            "the source is one coherent port",
            "generation-blind",
        ),
    )
    require_text(
        "python_code/item131_r4_support_dimension.py",
        (
            "generation flips do not erase the R4 syndrome",
            "generation is a finite copy label",
        ),
    )
    require_text(
        "python_code/axis_booking_resolution.py",
        (
            "Generation's \"3\" is the R1 count",
            "The old literal 'same three axes' reading is not a theorem",
        ),
    )
    require_text(
        "python_code/ckm_walk_signed_template.py",
        (
            "The sign is set by the GLOBAL walk-phase orientation",
            "not a from-nothing derivation",
        ),
    )

    print("\n[1] Current generation-blind algebra cannot make a directed latch")
    Pp = p_cycle(+1)
    A = a_k3()
    K = k_orientation()
    print(f"  dim commutant_S3(three-generation permutation rep) = {commutant_dimension()}")
    print(f"  residual of P_+ projected to S3 commutant span{{I,J}} = {projection_residual_to_s3_commutant(Pp):.6f}")
    print(f"  residual of A_K3 projected to S3 commutant span{{I,J}} = {projection_residual_to_s3_commutant(A):.6f}")
    check(commutant_dimension() == 2, "S3-blind operators are only the I/J two-parameter algebra")
    check(projection_residual_to_s3_commutant(Pp) > 0.9, "directed latch P_sigma is outside current S3-blind algebra")
    check(projection_residual_to_s3_commutant(A) < 1.0e-13, "Majorana support A_K3 is inside the S3-blind algebra")

    print("\n[2] The missing latch is exactly an orientation tensor")
    reconstructed = 0.5 * (A + K)
    print(f"  ||P_+ - (A_K3+K_or)/2|| = {np.linalg.norm(Pp - reconstructed):.3e}")
    check(np.linalg.norm(Pp - reconstructed) < 1.0e-13, "directed latch = symmetric support plus antisymmetric orientation")
    for perm, mat, parity in PERMS:
        moved = mat.T @ K @ mat
        if parity == 1:
            check(np.linalg.norm(moved - K) < 1.0e-13, f"even permutation {perm} preserves K_or")
        else:
            check(np.linalg.norm(moved + K) < 1.0e-13, f"odd permutation {perm} flips K_or")
    print("  K_or transforms in the sign representation.  It is the missing datum.")

    print("\n[3] Conditional C3-character phase closure")
    characters = {
        "trivial": 0.0,
        "faithful +": 2.0 * math.pi / 3.0,
        "faithful -": -2.0 * math.pi / 3.0,
    }
    for label, phi in characters.items():
        M = majorana_portal(abs(phi), sigma=+1 if phi >= 0 else -1)
        print(f"  {label:10s}: Phi={phi:+.9f}, CP norm={cp_norm(M):.6e}, I_CP={cp_invariant(M):+.6e}")
    check(cp_norm(majorana_portal(0.0, +1)) < 1.0e-12, "trivial C3 character is CP-even")
    check(cp_norm(majorana_portal(2.0 * math.pi / 3.0, +1)) > 1.0e-3, "faithful C3 character gives nonzero CP")
    check(
        abs(
            cp_invariant(majorana_portal(2.0 * math.pi / 3.0, +1))
            + cp_invariant(majorana_portal(2.0 * math.pi / 3.0, -1))
        )
        < 1.0e-12,
        "faithful conjugate characters reverse the CP sign",
    )

    print("\n[4] Absolute orientation remains unselected")
    mixed = 0.5 * (
        majorana_portal(2.0 * math.pi / 3.0, +1)
        + majorana_portal(2.0 * math.pi / 3.0, -1)
    )
    print(f"  unoriented faithful-character mixture CP norm = {cp_norm(mixed):.3e}")
    check(cp_norm(mixed) < 1.0e-12, "summing both orientations cancels CP")
    print("  The quark walk has a global spatial phase orientation, but canon currently")
    print("  does not identify generation labels with spatial axes.  Therefore that sign")
    print("  cannot yet select the sterile-generation latch orientation.")

    print("\n[5] Decision")
    print(
        """
  Not closed from current canon:
    The directed latch P_sigma needs K_or to be read as a service channel.
    The R1 boundary supplies such a sign-representation cochain, but existing
    R4/source/QEC repairs are generation-blind, within-generation, and have one
    singlet source port, so they average it away.  The readout theorem shows
    that a sign-representation environment pointer is the minimal missing
    primitive.

  Conditionally closable:
    If a new oriented sterile recovery latch is admitted and if its phase is
    required to be the faithful C3 character of that latch, then Phi=2*pi/3
    follows without fitting and gives the required CP-odd portal.

  Still not locked:
    The absolute orientation sigma remains a boundary/global-orientation datum
    unless a new theorem maps the substrate walk orientation to the generation
    orientation tensor K_or.  Without that selector, the two conjugate latches
    are equally allowed and their mixture is CP-even.

  Therefore the current framework cannot close all three targets.  It can
  reduce the missing theorem to one precise new bridge:

      derive the R1-boundary K_or cochain as a sterile-sector QEC recovery
      tensor and couple it to the global substrate orientation.
"""
    )
    print("ALL ASSERTIONS PASSED -- no full current-canon closure; C3 phase closes only under a new latch-character premise.")


if __name__ == "__main__":
    main()
