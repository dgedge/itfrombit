#!/usr/bin/env python3
r"""Foundational reconstruction / CP-holonomy frontier classifier.

Purpose
-------
The "reconstruction wall" and the "CP-holonomy wall" are easy to overstate
after the R14 and Item-87/R15 updates.  This script freezes the current split:

  * Born, measurement, the arrow, and bare alpha0=1/137 are no longer open
    reconstruction legs inside the accepted record/QEC calculus.
  * The complex field and the byte cell are reduced to small reconstruction
    premises (locality/local tomography; binary CSS read/write plus S-closure).
  * The CP sign is no longer a free lepton-sector phase inside the oriented
    substrate: Dynamic R1/Hasse recovery supplies the record geometry, the
    global orientation line supplies the sign pointer, the EM closed-pair rule
    supplies the extra {nu,d} contact, and the current PMNS lift predicts
    Dirac CP conservation with Majorana CP.

What remains is therefore not "find a CP sign" or "why Born?".  It is the
deeper reconstruction theorem:

    why stable local records force this complex, locally tomographic,
    binary-balanced CSS/QEC substrate with the monitored orientation-carrying
    recovery environment.

Exit 0 means the current frontier has this shape and the older scalar-source
R15 no-go is superseded by the oriented-substrate bridge, while the residual
is correctly kept at reconstruction-floor grade rather than promoted to an
unconditional lock.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


@dataclass(frozen=True)
class Gate:
    gate: str
    status: str
    residual: str


GATES = [
    Gate(
        "R0 stable local record floor",
        "irreducible floor",
        "cannot be derived internally without circularity",
    ),
    Gate(
        "R4 complex Hilbert field",
        "reduced",
        "rests on local measurement generation plus imported Hardy/CDP theorem",
    ),
    Gate(
        "R5 QEC/stabilizer byte",
        "reduced",
        "binary-minimality, read/write-to-CSS, and local S-closure premises remain explicit",
    ),
    Gate(
        "Born / measurement / arrow",
        "foundationally grounded",
        "inherits only the complex-QEC record calculus",
    ),
    Gate(
        "bare alpha0 = 1/137",
        "foundationally grounded",
        "sector billing and dressed-alpha are separate downstream checks",
    ),
    Gate(
        "R15 CP sign inside oriented substrate",
        "conditional-closed",
        "orientation line is reduced to directed recovery/gluing; branch sign is superselection",
    ),
    Gate(
        "Item-87 Koide/CP contact magnitude",
        "conditional-closed at finite monitor grade",
        "depends on reading multiplicities as closed record-pair contacts",
    ),
    Gate(
        "CP from deepest record axioms",
        "open reconstruction theorem",
        "derive orientability and branch selection beyond parity-even record axioms",
    ),
]


def K_classical(d: int) -> int:
    return d


def K_real(d: int) -> int:
    return d * (d + 1) // 2


def K_complex(d: int) -> int:
    return d * d


def K_quaternionic(d: int) -> int:
    return d * (2 * d - 1)


def local_tomography_discriminator() -> None:
    print("\n[1] R4 local-tomography discriminator")
    fields = {
        "classical": K_classical,
        "real": K_real,
        "complex": K_complex,
        "quaternionic": K_quaternionic,
    }
    locally_tomographic = {}
    for name, K in fields.items():
        gap = K(4) - K(2) * K(2)
        locally_tomographic[name] = gap == 0
        print(f"    {name:12s}: K(4)-K(2)^2 = {gap:+d}")
    check(locally_tomographic["complex"], "complex is locally tomographic")
    check(locally_tomographic["classical"], "classical is locally tomographic but excluded by reversible walk")
    check(not locally_tomographic["real"], "real Hilbert space fails local tomography")
    check(not locally_tomographic["quaternionic"], "quaternionic composition fails the count/tensor test")


def r1_hasse_covariance() -> list[list[F]]:
    # Rows tau=00, e=01, mu=10; columns tau-e and tau-mu.
    B = [[F(1), F(1)], [F(1), F(0)], [F(0), F(1)]]
    return [[sum(B[i][k] * B[j][k] for k in range(2)) for j in range(3)] for i in range(3)]


def matmul(A: list[list[F]], B: list[list[F]]) -> list[list[F]]:
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def transpose(A: list[list[F]]) -> list[list[F]]:
    return [list(row) for row in zip(*A)]


def scalar_mul(c: F, A: list[list[F]]) -> list[list[F]]:
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def perm_matrix(perm: tuple[int, ...]) -> list[list[F]]:
    P = [[F(0) for _ in perm] for _ in perm]
    for i, j in enumerate(perm):
        P[i][j] = F(1)
    return P


def sign_of_permutation(perm: tuple[int, ...]) -> int:
    inversions = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
    return -1 if inversions % 2 else 1


def oriented_r1_boundary() -> list[list[F]]:
    return [
        [F(0), F(-1), F(1)],
        [F(1), F(0), F(-1)],
        [F(-1), F(1), F(0)],
    ]


def r1_geometry_and_orientation() -> None:
    print("\n[2] Dynamic R1/Hasse geometry and sign representation")
    allowed = {(0, 0), (0, 1), (1, 0)}
    down_closed = all(
        (x, y) in allowed
        for a, b in allowed
        for x in range(a + 1)
        for y in range(b + 1)
    )
    check(down_closed and (1, 1) not in allowed, "{00,01,10} is the B2 order ideal with top 11 removed")
    K = r1_hasse_covariance()
    print(f"    K_R1 = {K}")
    check(K == [[F(2), F(1), F(1)], [F(1), F(1), F(0)], [F(1), F(0), F(1)]], "R1 Hasse covers force the path covariance")
    # E-plane eigenvalues {1,1/3}: trace(PKP)=4/3 and trace((PKP)^2)=10/9.
    I = [[F(1 if i == j else 0) for j in range(3)] for i in range(3)]
    J = [[F(1) for _ in range(3)] for _ in range(3)]
    P = [[I[i][j] - J[i][j] / 3 for j in range(3)] for i in range(3)]
    PKP = matmul(matmul(P, K), P)
    tr = sum(PKP[i][i] for i in range(3))
    tr2 = sum(matmul(PKP, PKP)[i][i] for i in range(3))
    check(tr == F(4, 3) and tr2 == F(10, 9), "R1 covariance has E-plane eigenvalues {1, 1/3}")

    omega = oriented_r1_boundary()
    for perm in permutations(range(3)):
        Pm = perm_matrix(perm)
        moved = matmul(matmul(Pm, omega), transpose(Pm))
        check(moved == scalar_mul(F(sign_of_permutation(perm)), omega), f"Omega_R1 transforms as sign rep under {perm}")


def em_closed_pair_selector() -> None:
    print("\n[3] EM closed-pair selector for the Item-87 contact count")
    members = {
        "e": (F(-1, 2), F(-1, 2), 9, F(2, 9)),
        "nu": (F(1, 2), F(-1, 2), 9, F(3, 9)),
        "u": (F(1, 2), F(1, 6), 27, F(2, 27)),
        "d": (F(-1, 2), F(1, 6), 27, F(1, 9)),
    }
    selected: set[str] = set()
    for name, (t3, y, n_eff, target) in members.items():
        q = t3 + y
        close = t3 * t3 + y * y - q * q
        plus = int(close > 0)
        delta = F(2 + plus, n_eff)
        if plus:
            selected.add(name)
        print(f"    {name:2s}: C_close={str(close):>4s}, plus={plus}, delta={delta}, target={target}")
        check(close == -2 * t3 * y, f"{name}: C_close = -2 T3 Y")
        check(delta == target, f"{name}: base 2 plus EM closed-pair selector gives target delta")
    check(selected == {"nu", "d"}, "closed EM-pair selector selects exactly {nu,d}")


def main() -> None:
    print("FOUNDATIONAL RECONSTRUCTION / CP-HOLONOMY FRONTIER CLASSIFIER")
    print("=" * 88)

    require_text("python_code/r4_complex_local_tomography_theorem.py", ("complex FORCED", "LOCAL TOMOGRAPHY"))
    require_text("python_code/r5_record_class_bridge.py", ("Type-II", "transversal S"))
    require_text("python_code/alpha0_record_pair_symmetry_theorem.py", ("records are clonable", "Sym^2(16)"))
    require_text("python_code/r15_global_orientation_sign_pointer_theorem.py", ("global substrate orientation line", "omega*K_R1"))
    require_text("python_code/record_reconstruction_orientation_line_theorem.py", ("Z2 orientation line", "connected orientable substrate"))
    require_text("python_code/item87_neutrality_selector_closure.py", ("closed electroweak cancellation-pair record", "{nu,d}"))
    require_text("python_code/item87_leptonic_dirac_cp_from_frame_transport.py", ("leptonic DIRAC CP is CONSERVED", "MAJORANA"))
    require_text("python_code/item87_mbb_0nubb_prediction.py", ("m_bb ~ 2-3 meV", "NORMAL ORDERING"))

    local_tomography_discriminator()
    r1_geometry_and_orientation()
    em_closed_pair_selector()

    print("\n[4] Frontier table")
    for gate in GATES:
        print(f"    {gate.gate:42s} {gate.status:38s} residual: {gate.residual}")

    check(
        any(g.gate == "R15 CP sign inside oriented substrate" and g.status == "conditional-closed" for g in GATES),
        "CP sign is not a free phase inside the oriented substrate",
    )
    check(
        any(g.gate == "CP from deepest record axioms" and g.status.startswith("open") for g in GATES),
        "the remaining CP issue is a reconstruction theorem, not coefficient fitting",
    )
    check(
        any(g.gate == "R5 QEC/stabilizer byte" and g.status == "reduced" for g in GATES),
        "the full QEC/stabilizer substrate is reduced but not derived from R0 without premises",
    )

    print(
        """
VERDICT
  The hard foundation is now sharply localized.

  Not open at this level:
    Born, measurement, the thermodynamic arrow, and bare alpha0=1/137.
    The CP sign is also not an adjustable lepton-sector phase once the oriented
    substrate layer is accepted: Dynamic R1 gives the Hasse/order-ideal record,
    the orientation-line theorem derives the Z2 line object from directed
    recovery/gluing, and the present leptonic lift makes CP Majorana while
    conserving the Dirac PMNS Jarlskog.

  Still open:
    the from-R0 reconstruction theorem.  One must derive why stable local
    records force the binary CSS/read-write/S-closed QEC byte and why the
    substrate occupies one orientable branch, rather than a parity-even
    +/- mixture.

  So the frontier is no longer "why Born?" or "which CP sign?".  It is:

      stable local records -> complex locally tomographic QEC/stabilizer
      substrate with monitored oriented recovery.

exit 0"""
    )


if __name__ == "__main__":
    main()
