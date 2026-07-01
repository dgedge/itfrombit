#!/usr/bin/env python3
r"""R15 recovery-holonomy / CP-origin closure status.

Question
--------
Can the matter/antimatter CP sign be closed as a holonomy of record repair?

Answer
------
Not as a full current-canon derivation.  The problem has split into three
clauses:

  1. Operator support: closed conditionally.  The only S3-blind symmetric
     off-diagonal sterile-generation Majorana support is A_K3.
  2. Holonomy phase: closed conditionally.  If the missing latch is a faithful
     C3 recovery character, Phi = +/- 2*pi/3 is forced and gives nonzero CP.
  3. Physical orientation/readout: not closed by a scalar source.  The oriented R1 boundary
     cochain K_R1 exists, but a scalar/generation-blind Stinespring environment
     averages it to zero.  Sterile recovery reads it iff the environment
     carries a sign-representation pointer.  Source/R4 algebra does not
     contain that pointer; r15_global_orientation_sign_pointer_theorem.py
     supplies it conditionally from the global substrate orientation line.

So this script is the scalar-source iff/no-go boundary.  It is superseded, for
the oriented-substrate theory, by r15_global_orientation_sign_pointer_theorem.py:
the sign pointer is conditionally the global orientation line.  R15 still is not
a Locked from-nothing CP-origin theorem, because the absolute global orientation
is a substrate superselection convention.
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


def perms() -> list[tuple[tuple[int, ...], np.ndarray, int]]:
    out = []
    eye = np.eye(3)
    for p in permutations(range(3)):
        inv = sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j])
        out.append((p, eye[list(p)], -1 if inv % 2 else 1))
    return out


PERMS = perms()


def a_k3() -> np.ndarray:
    return np.ones((3, 3), dtype=complex) - np.eye(3, dtype=complex)


def p_cycle(sigma: int) -> np.ndarray:
    p = np.zeros((3, 3), dtype=complex)
    cycle = (0, 1, 2) if sigma == 1 else (0, 2, 1)
    for src, dst in zip(cycle, cycle[1:] + cycle[:1]):
        p[dst, src] = 1.0
    return p


def k_r1() -> np.ndarray:
    """Oriented R1-boundary cochain 00 -> 10 -> 01 -> 00."""

    return p_cycle(+1) - p_cycle(+1).T


def ordinary_average(M: np.ndarray) -> np.ndarray:
    out = np.zeros_like(M)
    for _perm, p, _parity in PERMS:
        out += p.T @ M @ p
    return out / len(PERMS)


def sign_average(M: np.ndarray) -> np.ndarray:
    out = np.zeros_like(M)
    for _perm, p, parity in PERMS:
        out += parity * (p.T @ M @ p)
    return out / len(PERMS)


def cp_invariant(M: np.ndarray) -> float:
    md = M.conj().T
    return float(np.imag(np.trace(H_DIRAC @ md @ M @ md @ H_DIRAC @ M)))


def cp_norm(M: np.ndarray) -> float:
    md = M.conj().T
    probes = np.array(
        [
            np.imag(np.trace(H_DIRAC @ md @ M @ md @ H_DIRAC @ M)),
            np.imag(np.trace(H_DIRAC @ M @ md @ M @ md @ M)),
        ],
        dtype=float,
    )
    return float(np.linalg.norm(probes))


def majorana_portal(phi: float, sigma: int, r: float = 0.5) -> np.ndarray:
    return np.eye(3, dtype=complex) + r * np.exp(1j * sigma * phi) * a_k3()


def s3_invariant_offdiag_null_vector() -> np.ndarray:
    basis = [
        np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex),
        np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex),
        np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex),
    ]
    rows = []
    for _perm, p, _parity in PERMS:
        for b in basis:
            pass
        deltas = [p.T @ b @ p - b for b in basis]
        for i in range(3):
            for j in range(3):
                rows.append([d[i, j] for d in deltas])
    mat = np.array(rows, dtype=complex)
    _u, _s, vh = np.linalg.svd(mat)
    vec = vh[-1].real
    vec /= np.min(np.abs(vec[np.abs(vec) > TOL]))
    if np.sum(vec) < 0:
        vec *= -1
    return vec


def main() -> None:
    require_text(
        "python_code/item87_lepton_cp_obstruction.py",
        (
            "Leptonic CP vanishes in every physical lepton sector",
            "existence gap",
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

    print("R15 RECOVERY-HOLONOMY / CP-ORIGIN CLOSURE STATUS")
    print("=" * 88)

    print("\n[1] Existing walk obstruction")
    print("  The documented lepton walk is real/Delta-L=0 and has no physical")
    print("  neutrino-sector CP invariant; this is the need for a new Majorana latch.")
    check(True, "existing item87 obstruction is imported as a proven no-go")

    print("\n[2] Operator support")
    vec = s3_invariant_offdiag_null_vector()
    print(f"  S3-invariant symmetric off-diagonal coefficient vector = {vec}")
    check(np.max(np.abs(vec - np.ones(3))) < TOL, "S3-blind symmetric support is uniquely A_K3")
    for sigma in (+1, -1):
        check(
            np.linalg.norm(p_cycle(sigma) + p_cycle(sigma).T - a_k3()) < TOL,
            f"sigma={sigma:+d}: Majorana symmetrisation of directed latch gives A_K3",
        )

    print("\n[3] Conditional phase closure")
    for sigma in (+1, -1):
        phi = 2.0 * math.pi / 3.0
        M = majorana_portal(phi, sigma=sigma)
        print(f"  sigma={sigma:+d}, Phi=2pi/3: CP norm={cp_norm(M):.6e}, I_CP={cp_invariant(M):+.6e}")
        check(cp_norm(M) > 1.0e-3, f"sigma={sigma:+d}: faithful C3 character gives nonzero CP")
    check(
        abs(cp_invariant(majorana_portal(2.0 * math.pi / 3.0, +1)) + cp_invariant(majorana_portal(2.0 * math.pi / 3.0, -1))) < TOL,
        "conjugate faithful characters reverse the CP sign",
    )

    print("\n[4] Orientation/readout obstruction")
    A = a_k3()
    K = k_r1()
    P = 0.5 * (A + K)
    print(f"  ||<K_R1>_S3||                  = {np.linalg.norm(ordinary_average(K)):.3e}")
    print(f"  ||<P_+>_S3 - A_K3/2||          = {np.linalg.norm(ordinary_average(P) - 0.5 * A):.3e}")
    print(f"  ||<K_R1>_sign - K_R1||         = {np.linalg.norm(sign_average(K) - K):.3e}")
    print(f"  ||<P_+>_sign - K_R1/2||        = {np.linalg.norm(sign_average(P) - 0.5 * K):.3e}")
    check(np.linalg.norm(ordinary_average(K)) < TOL, "scalar/generation-blind recovery erases K_R1")
    check(np.linalg.norm(sign_average(K) - K) < TOL, "sign pointer reads K_R1")
    check(np.linalg.norm(ordinary_average(P) - 0.5 * A) < TOL, "without sign pointer the latch becomes unoriented")

    print("\n[5] Absolute sign")
    mixed = 0.5 * (majorana_portal(2.0 * math.pi / 3.0, +1) + majorana_portal(2.0 * math.pi / 3.0, -1))
    print(f"  CP norm of unoriented faithful-character mixture = {cp_norm(mixed):.3e}")
    check(cp_norm(mixed) < TOL, "summing both orientations cancels CP")

    print(
        """
VERDICT
  R15 is not fully closable from a scalar/generation-blind source alone.

  Closed/conditional:
    * The only admissible S3-blind Majorana support is A_K3.
    * A faithful C3 recovery character conditionally fixes Phi=2*pi/3 and gives
      a nonzero CP invariant.

  Blocked:
    * Current sterile source/R4 repair algebra is scalar and generation-blind;
      it has one coherent source port and no G0/G1 orientation clause.
    * Such a scalar environment provably averages K_R1 to zero.
    * Without a sign-representation pointer, the two orientations are summed
      and the CP invariant cancels.

  Superseding oriented-substrate bridge:
    r15_global_orientation_sign_pointer_theorem.py constructs that sign pointer
    conditionally as the global substrate-orientation line.  The remaining
    floor is the absolute global orientation convention, not a free
    sterile-sector phase/sign knob.
exit 0 -- scalar-source no-go retained; sign-pointer bridge closes conditionally via global orientation.
"""
    )


if __name__ == "__main__":
    main()
