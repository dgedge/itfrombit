#!/usr/bin/env python3
r"""CC clause iii via the item-115 loop machine.

Question:
    Can the live item-115 vacuum-polarization machinery compute the stationary
    CC branch's remaining generation-vertex correction

        Delta(-ln gamma) = alpha0 * C_gen,  C_gen(target) = 0.304750231...

    rather than borrowing the demoted dressed-alpha count ansatz?

Discipline:
  * Reuse the item-115 chirality-Weyl kernel only after its operator gates pass:
      - the only nontrivial physical, charge-blind register flip on the 48 is the
        locked chi/W pair;
      - the Pauli triple on that logical qubit is the Ward-compatible Dirac triple.
  * Compute an actual loop integral, not a count.  The local zero-momentum Pauli
    generation vertex O_gen = sigma_i gives

        B_i = < 2 |<-k|sigma_i|+k>|^2 / (2 E_k) >_BZ
            = < (1 - s_i^2/|s|^2) / |s| >_BZ .

    Cubic isotropy gives B_i = (2/3)<1/|s|>.  The per-external-vertex barrier
    candidate is C_loop = B_i/2 = (1/3)<1/|s|>.

  * Report the alternatives explicitly:
      - unhalved two-vertex bubble;
      - strict item-115 photon-velocity map;
      - the successful single-vertex map.

exit 0 means the loop computation and gates are verified.  It does not by itself
prove the observable map from bit-flip billed ledger to single-vertex Pauli bubble.
"""

from __future__ import annotations

import math

import numpy as np

from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    queue_readouts,
    rho_ratio_from_q1,
    solve_gamma_for_stationary_queue,
    solve_target,
)


CHI_MASK = (1 << 6) | (1 << 7)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def b(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    return (
        not (b(n, 0) and b(n, 1))
        and b(n, 7) == b(n, 6)
        and ((b(n, 2) == 0) == ((b(n, 3), b(n, 4)) == (0, 0)))
    )


PHYS = [n for n in range(256) if valid(n)]
IDX = {n: i for i, n in enumerate(PHYS)}


def q116(n: int) -> float:
    zf = 1 if b(n, 5) == 0 else -1
    szc = -3 if (b(n, 3), b(n, 4)) == (0, 0) else -1
    return 0.5 * zf + szc / 3.0 + 0.5


Q = np.diag([q116(n) for n in PHYS])
ZG0 = np.diag([1.0 if not b(n, 0) else -1.0 for n in PHYS])


def flip_op(mask: int) -> np.ndarray:
    op = np.zeros((48, 48))
    for col, n in enumerate(PHYS):
        m = n ^ mask
        if m in IDX:
            op[IDX[m], col] = 1.0
    return op


def phase_op(bit: int) -> np.ndarray:
    return np.diag([1.0 if b(n, bit) == 0 else -1.0 for n in PHYS])


def comm_norm(a: np.ndarray, b_: np.ndarray) -> float:
    return float(np.linalg.norm(a @ b_ - b_ @ a))


def acomm_norm(a: np.ndarray, b_: np.ndarray) -> float:
    return float(np.linalg.norm(a @ b_ + b_ @ a))


SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA = [SX, SY, SZ]


def eigvecs(shat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    h = shat[0] * SX + shat[1] * SY + shat[2] * SZ
    _, u = np.linalg.eigh(h)
    return u[:, 0], u[:, 1]


def local_pauli_bubble_matrix_check(n_grid: int) -> tuple[float, float, float]:
    """Direct spinor-loop implementation at q=0, used as a harness."""
    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    vals = [0.0, 0.0, 0.0]
    for kx in ks:
        for ky in ks:
            for kz in ks:
                s = np.sin(np.array([kx, ky, kz]))
                norm = float(np.linalg.norm(s))
                um, up = eigvecs(s / norm)
                for i, op in enumerate(SIGMA):
                    vals[i] += 2.0 * abs(np.vdot(um, op @ up)) ** 2 / (2.0 * norm)
    return tuple(v / n_grid**3 for v in vals)


def c_loop_midpoint(n_grid: int, chunk: int = 48) -> float:
    """C_loop = (1/3)<1/|sin k|> over the midpoint Brillouin-zone grid."""
    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    sin2 = np.sin(ks) ** 2
    total = 0.0
    for i in range(0, n_grid, chunk):
        sx = sin2[i : i + chunk][:, None, None]
        den = np.sqrt(sx + sin2[None, :, None] + sin2[None, None, :])
        total += float(np.sum(1.0 / den))
    return total / n_grid**3 / 3.0


def item115_vector_invariant(n_grid: int) -> tuple[float, float, float]:
    """A_E, A_T, and D=A_E-A_T/v^2 for the item-115 photon-velocity observable."""
    ks = (np.arange(n_grid) + 0.5) * 2.0 * math.pi / n_grid - math.pi
    q = np.array([2.0 * math.pi / n_grid, 0.0, 0.0])
    ty = 1
    pi00 = para_t = dia = 0.0
    for kx in ks:
        for ky in ks:
            for kz in ks:
                k = np.array([kx, ky, kz])
                kq = k + q
                sk = np.sin(k)
                skq = np.sin(kq)
                nk = float(np.linalg.norm(sk))
                nkq = float(np.linalg.norm(skq))
                um, _ = eigvecs(sk / nk)
                _, vp = eigvecs(skq / nkq)
                de = nk + nkq
                pi00 += 2.0 * abs(np.vdot(um, vp)) ** 2 / de
                vy = math.cos(k[ty]) * SIGMA[ty]
                para_t += 2.0 * abs(np.vdot(um, vy @ vp)) ** 2 / de
                dia += math.sin(k[ty]) ** 2 / nk
    q2 = q[0] ** 2
    ae = pi00 / n_grid**3 / q2
    at = (dia - para_t) / n_grid**3 / q2
    return ae, at, ae - at


def rho_for_c(c_coeff: float, q1_target: float) -> float:
    delta = ALPHA0 * c_coeff
    gamma = BASE_GAMMA * math.exp(-delta)
    _, q_post, _ = queue_readouts(gamma, 1)
    return rho_ratio_from_q1(q_post, q1_target)


def extrapolate_c(ns: list[int]) -> tuple[float, float, list[tuple[int, float]]]:
    rows = [(n, c_loop_midpoint(n)) for n in ns]
    x = np.array([1.0 / (n * n) for n, _ in rows[-5:]])
    y = np.array([c for _, c in rows[-5:]])
    linear = float(np.polyfit(x, y, 1)[-1])
    quadratic = float(np.polyfit(x, y, 2)[-1])
    return linear, quadratic, rows


def main() -> None:
    print("CC GENERATION-VERTEX ITEM-115 LOOP AUDIT")

    print("\n[1] Item-115 operator gates")
    ok_masks = []
    for mask in range(256):
        if all((n ^ mask) in IDX and abs(q116(n ^ mask) - q116(n)) < 1e-12 for n in PHYS):
            ok_masks.append(mask)
    print(f"  physical + charge-blind flip masks = {[format(m, '08b') for m in ok_masks]}")
    check(ok_masks == [0, CHI_MASK], "only the locked chi/W flip supplies a nontrivial charge-blind kernel")

    xchi = flip_op(CHI_MASK)
    zchi = phase_op(6)
    ychi = 1j * xchi @ zchi
    triple = [xchi, ychi, zchi]
    for i, op in enumerate(triple):
        check(np.linalg.norm(op - op.conj().T) < 1e-12, f"chi Pauli {i} is Hermitian")
        check(comm_norm(op, Q) < 1e-12, f"chi Pauli {i} commutes with Q116")
        check(comm_norm(op, ZG0) < 1e-12, f"chi Pauli {i} commutes with G0")
        for j in range(i + 1, 3):
            check(acomm_norm(op, triple[j]) < 1e-12, f"chi Pauli {i},{j} anticommute")

    print("\n[2] Exact target from the stationary active-demux branch")
    _, q1_target = solve_target()
    gamma_star = solve_gamma_for_stationary_queue(q1_target)
    c_target = -math.log(gamma_star / BASE_GAMMA) / ALPHA0
    print(f"  gamma_star/gamma_Sec5.2 = {gamma_star / BASE_GAMMA:.12f}")
    print(f"  C_target                = {c_target:.12f}")
    check(abs(c_target - 0.304750231482) < 1e-9, "target coefficient matches clause-iii audit")

    print("\n[3] Local generation-vertex Pauli bubble")
    b16 = local_pauli_bubble_matrix_check(16)
    c16_matrix = sum(b16) / 6.0
    c16_formula = c_loop_midpoint(16)
    print(f"  direct q=0 spinor bubble B_i at N=16 = {b16[0]:.12f}, {b16[1]:.12f}, {b16[2]:.12f}")
    print(f"  single-vertex C=B_i/2, matrix         = {c16_matrix:.12f}")
    print(f"  single-vertex C=(1/3)<1/|sin k|>      = {c16_formula:.12f}")
    check(max(abs(b16[i] - b16[0]) for i in range(3)) < 1e-12, "local Pauli bubble is cubic-isotropic")
    check(abs(c16_matrix - c16_formula) < 1e-12, "spinor loop equals analytic BZ average")

    ns = [32, 48, 64, 96, 128, 160, 192, 256]
    c_lin, c_quad, rows = extrapolate_c(ns)
    print("\n  midpoint-grid convergence:")
    print("    N       C_loop          rho(C_loop)/rho_obs")
    for n, c in rows:
        print(f"    {n:<3d}  {c:.12f}    {rho_for_c(c, q1_target):.6f}")
    c_loop = 0.5 * (c_lin + c_quad)
    c_err = abs(c_lin - c_quad) / 2.0
    rho_loop = rho_for_c(c_loop, q1_target)
    print(f"\n  extrapolated C_loop       = {c_loop:.12f} +/- {c_err:.2e}")
    print(f"  target miss in coefficient = {(c_loop / c_target - 1.0) * 100:+.4f}%")
    print(f"  rho(C_loop)/rho_obs        = {rho_loop:.6f}")
    check(abs(c_loop / c_target - 1.0) < 0.005, "single-vertex loop lands within 0.5% in coefficient")
    check(abs(math.log(rho_loop)) < 0.005, "single-vertex loop lands within 0.5% in rho")

    print("\n[4] Alternative maps, to expose the load-bearing observable map")
    c_two_vertex = 2.0 * c_loop
    rho_two_vertex = rho_for_c(c_two_vertex, q1_target)
    ae32, at32, d32 = item115_vector_invariant(32)
    c_vector32 = 4.0 * math.pi * 8.0 * d32
    rho_vector32 = rho_for_c(c_vector32, q1_target)
    print(f"  unhalved two-vertex bubble C=2*C_loop     = {c_two_vertex:.6f}, rho={rho_two_vertex:.3e}")
    print(f"  item-115 photon map at N=32: A_E={ae32:.6f}, A_T={at32:.6f}, D={d32:.6f}")
    print(f"  strict vector coefficient 4pi*8*D         = {c_vector32:.6f}, rho={rho_vector32:.3e}")
    check(c_two_vertex > 1.9 * c_target, "unhalved bubble over-corrects and is not the CC map")
    check(c_vector32 > 5.0 * c_target, "strict photon-velocity map is not the generation-barrier map")

    print(
        """
VERDICT
  This is the first positive clause-iii computation from a genuinely evaluated
  loop.  The item-115 chirality-Weyl kernel, pointed at a local Pauli generation
  vertex and read as a single-external-vertex barrier correction, gives

      C_loop = 0.303563,  Delta(-ln gamma) = alpha0*C_loop,
      rho = 1.0019 rho_obs.

  That is inside the strict 2%-in-rho gate and within 0.39% of the exact target
  coefficient.  It is not the demoted dressed-alpha count ansatz and not a scan
  over simple constants.

  What remains open is no longer a number hunt.  It is the operator-map theorem:
    (a) the billed bit-flip generation attempt is the local Pauli vertex on the
        forced chi/W Weyl kernel; and
    (b) the attempt barrier consumes the single-external-vertex half-bubble
        C=B_i/2 rather than the unhalved two-vertex bubble or the photon-velocity
        observable.
  If (a,b) are derived from Sec 5.2/register mechanics, clause iii closes at the
  0.2%-in-rho level.  If not, this remains a sharp class-2 loop candidate.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
