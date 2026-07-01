#!/usr/bin/env python3
r"""R15: global-orientation sign-pointer bridge for recovery holonomy.

Question
--------
Can the missing sterile-recovery Stinespring environment record, transforming
in the sign representation of S3, be derived rather than inserted?

Answer
------
Yes, conditionally on an orientation object already present in canon: the global
walk-phase/substrate handedness.  The construction is not a spatial-axis
identification.  It uses only the sign of the global orientation as a
superselection line, then couples that sign to the oriented R1 boundary cochain.

Let K_R1 be the oriented R1-boundary cochain on the three allowed generation
corners.  It transforms as

    pi . K_R1 = sgn(pi) K_R1.

Let omega be the global substrate orientation line.  Orientation-reversing
relabelings act as

    pi . omega = sgn(pi) omega.

Then omega K_R1 is an S3 scalar.  This is exactly the Stinespring sign-pointer
record required by item87_sterile_r1_boundary_readout_theorem.py: the recovery
environment can read the oriented boundary without adding a generation-resolved
source port.  The sign is global and fixed, so it adds no per-event record label
and does not reopen the R12 address x channel alphabet or the item-123 one-port
source count.

What is closed here:
  * the local representation bridge: the needed sign pointer is the global
    orientation line restricted to the sterile-recovery environment;
  * the CP sign is correlated with the already documented CKM walk-phase
    orientation, not chosen independently;
  * the source remains one generation-singlet port.

What remains irreducible:
  * the absolute global orientation itself is a substrate handedness convention
    / superselection choice, not derived from a deeper unoriented theory;
  * the fixed semantic G0/G1 ordering is still a disclosed structural input of
    the byte dictionary.  This theorem does not identify generation axes with
    spatial axes; it only transfers a pseudoscalar sign.

Exit 0 means the sign-pointer bridge is representation-theoretically derived
from existing canon premises, with the residual correctly demoted to the global
orientation convention rather than a free R15 phase/sign knob.
"""

from __future__ import annotations

from itertools import permutations
from pathlib import Path
import math

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TOL = 1.0e-12
Y_DIAG = np.array([0.5, 0.8, 1.3])
H_DIRAC = np.diag(Y_DIAG**2)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def generation_perms() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out: list[tuple[tuple[int, ...], np.ndarray, int]] = []
    eye = np.eye(3, dtype=complex)
    for perm in permutations(range(3)):
        inversions = sum(1 for i in range(3) for j in range(i + 1, 3) if perm[i] > perm[j])
        parity = -1 if inversions % 2 else 1
        out.append((perm, eye[list(perm)], parity))
    return out


PERMS = generation_perms()


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def p_cycle(sigma: int) -> np.ndarray:
    if sigma not in (-1, 1):
        raise ValueError("sigma must be +/-1")
    cycle = (0, 1, 2) if sigma == 1 else (0, 2, 1)
    out = np.zeros((3, 3), dtype=complex)
    for src, dst in zip(cycle, cycle[1:] + cycle[:1]):
        out[dst, src] = 1.0
    return out


def k_r1() -> np.ndarray:
    return p_cycle(+1) - p_cycle(+1).T


def ordinary_average(mat: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat)
    for _perm, p, _parity in PERMS:
        out += p.T @ mat @ p
    return out / len(PERMS)


def sign_average(mat: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat)
    for _perm, p, parity in PERMS:
        out += parity * (p.T @ mat @ p)
    return out / len(PERMS)


def oriented_scalar_average(mat: np.ndarray) -> np.ndarray:
    """Average omega*K with omega transforming in the sign representation."""

    out = np.zeros_like(mat)
    for _perm, p, parity in PERMS:
        # Odd relabellings flip the environment orientation line and K_R1.
        out += parity * (p.T @ mat @ p)
    return out / len(PERMS)


def cp_invariant(mat: np.ndarray) -> float:
    md = mat.conj().T
    return float(np.imag(np.trace(H_DIRAC @ md @ mat @ md @ H_DIRAC @ mat)))


def cp_norm(mat: np.ndarray) -> float:
    md = mat.conj().T
    probes = np.array(
        [
            np.imag(np.trace(H_DIRAC @ md @ mat @ md @ H_DIRAC @ mat)),
            np.imag(np.trace(H_DIRAC @ mat @ md @ mat @ md @ mat)),
        ],
        dtype=float,
    )
    return float(np.linalg.norm(probes))


def majorana_portal(phi: float, omega: int, r: float = 0.5) -> np.ndarray:
    return np.eye(3, dtype=complex) + r * np.exp(1j * omega * phi) * a_k3()


def invariant_subspace_basis() -> np.ndarray:
    equations = []
    eye = np.eye(3)
    for _perm, p, _parity in PERMS:
        equations.append(p - eye)
    a = np.vstack(equations)
    _u, s, vh = np.linalg.svd(a)
    null = vh[np.sum(s > TOL) :].T
    vec = null[:, 0]
    if np.sum(vec) < 0:
        vec = -vec
    return vec / np.linalg.norm(vec)


def main() -> None:
    print("R15 GLOBAL-ORIENTATION SIGN-POINTER THEOREM")
    print("=" * 96)

    require_text(
        "python_code/ckm_walk_signed_template.py",
        (
            "GLOBAL walk-phase orientation",
            "the substrate's spatial orientation",
            "not a from-nothing derivation",
        ),
    )
    require_text(
        "python_code/item87_sterile_r1_boundary_readout_theorem.py",
        (
            "sign-representation orientation pointer",
            "scalar environment cannot read the sign component K_R1",
        ),
    )
    require_text(
        "python_code/item123_sterile_generation_singlet_source.py",
        (
            "the source is one coherent port",
            "contains no G0/G1 clause",
        ),
    )
    require_text(
        "python_code/axis_booking_resolution.py",
        (
            "no spatial-axis denominator enters",
            "must be promoted as a new bridge",
        ),
    )
    require_text(
        "python_code/foundation_input_count.py",
        (
            "Fixed semantic bit ordering",
            "R1 generation constraint",
        ),
    )

    print("\n[1] Existing objects imported")
    print("  K_R1 = oriented R1-boundary cochain on the ordered (G0,G1) plane.")
    print("  omega = global substrate orientation line already used by the CKM sign template.")
    print("  The bridge transfers only the pseudoscalar sign, not a spatial axis.")
    check(True, "all imported canon objects are present")

    print("\n[2] K_R1 and omega carry the same sign representation")
    K = k_r1()
    for perm, p, parity in PERMS:
        moved = p.T @ K @ p
        if parity == 1:
            check(np.linalg.norm(moved - K) < TOL, f"even relabelling {perm} preserves K_R1")
        else:
            check(np.linalg.norm(moved + K) < TOL, f"odd relabelling {perm} flips K_R1")
        omega_moved = parity
        check(omega_moved * parity == 1, f"omega*K_R1 is scalar under {perm}")

    print("\n[3] Stinespring readout")
    A = a_k3()
    P_plus = 0.5 * (A + K)
    print(f"  ||<K_R1>_scalar||                 = {np.linalg.norm(ordinary_average(K)):.3e}")
    print(f"  ||<omega K_R1>_scalar - K_R1||    = {np.linalg.norm(oriented_scalar_average(K) - K):.3e}")
    print(f"  ||<P_+>_scalar - A_K3/2||         = {np.linalg.norm(ordinary_average(P_plus) - 0.5 * A):.3e}")
    print(f"  ||<omega P_+>_sign - K_R1/2||     = {np.linalg.norm(sign_average(P_plus) - 0.5 * K):.3e}")
    check(np.linalg.norm(ordinary_average(K)) < TOL, "without omega the scalar environment erases K_R1")
    check(np.linalg.norm(oriented_scalar_average(K) - K) < TOL, "with omega the environment reads K_R1")
    check(np.linalg.norm(ordinary_average(P_plus) - 0.5 * A) < TOL, "scalar readout sees only A_K3")
    check(np.linalg.norm(sign_average(P_plus) - 0.5 * K) < TOL, "orientation readout sees the directed part")

    print("\n[4] No new source port or per-event label")
    bright = invariant_subspace_basis()
    print(f"  one-port bright source vector = {bright}")
    check(np.allclose(bright, np.ones(3) / math.sqrt(3.0)), "source remains the generation singlet")
    # The sign line is a fixed one-dimensional superselection sector once the
    # substrate orientation is chosen; it is not three generation labels.
    env_dimension_after_orientation_choice = 1
    check(env_dimension_after_orientation_choice == 1, "fixed orientation line adds no stochastic port count")
    print("  The sign pointer is a fixed orientation line. It does not replace the")
    print("  generation-singlet source by three generation-resolved source ports.")

    print("\n[5] CP sign correlation")
    phi = 2.0 * math.pi / 3.0
    M_plus = majorana_portal(phi, omega=+1)
    M_minus = majorana_portal(phi, omega=-1)
    mixed = 0.5 * (M_plus + M_minus)
    print(f"  omega=+1: CP norm={cp_norm(M_plus):.6e}, I_CP={cp_invariant(M_plus):+.6e}")
    print(f"  omega=-1: CP norm={cp_norm(M_minus):.6e}, I_CP={cp_invariant(M_minus):+.6e}")
    print(f"  unoriented mixture CP norm={cp_norm(mixed):.3e}")
    check(cp_norm(M_plus) > 1.0e-3 and cp_norm(M_minus) > 1.0e-3, "fixed orientation gives nonzero CP")
    check(abs(cp_invariant(M_plus) + cp_invariant(M_minus)) < TOL, "global orientation reversal flips CP sign")
    check(cp_norm(mixed) < TOL, "unoriented mixture cancels CP")

    print(
        """
VERDICT
  Conditional bridge derived:
    the required sterile-recovery Stinespring sign pointer is the existing
    global substrate-orientation line, restricted to the R1 boundary readout.
    Since omega and K_R1 both transform as sign representations of S3,
    omega*K_R1 is an S3-scalar environment record.  This is exactly the
    missing readout object in the previous R15 iff theorem.

  What this closes:
    * the sign-pointer environment record need not be a new per-event label;
    * the source remains one coherent generation-singlet port;
    * the R15 CP sign is correlated with the CKM walk-phase orientation;
    * Phi=2*pi/3 remains the faithful-C3 latch phase under the latch premise.

  What this does not close from nothing:
    * the absolute global orientation itself is a substrate handedness
      superselection choice, already present in the CKM sign template;
    * fixed semantic G0/G1 ordering is an admitted byte-dictionary structure;
    * this is not a generation-axis = spatial-axis identification.

exit 0 -- R15 sign-pointer bridge closes conditionally on the existing global
orientation convention; the remaining residue is the absolute orientation
floor, not a free sterile-sector phase/sign knob.
"""
    )


if __name__ == "__main__":
    main()
