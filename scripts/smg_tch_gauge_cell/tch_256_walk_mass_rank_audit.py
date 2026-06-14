#!/usr/bin/env python3
"""
256-state TCH walk/mass-operator color-rank audit.

Question:
  Can C_light = rank-one and C_middle = C_heavy = I_3 be derived from the
  documented 256-state TCH Yukawa/walk mass operator?

Documented operator implemented here:
  Part 01 defines the 8-bit register

      (G0, G1, LQ, C0, C1, I3, chi, W)

  and the 256-state walk

      U = sum_{k=0}^7 A_k CNOT^(k),
      A_0 = sqrt(1-delta),
      A_k = sqrt(delta/7) exp(i k pi/4), k=1..7,

  where CNOT^(k) acts on

      ctrl=(2-k) mod 8, target=(5-k) mod 8.

  The paper names M1=U^dagger U and M2=(U^dagger U)^2. The CKM paper then
  uses the 9-state left-handed quark projection with a Boltzmann color-singlet
  vector w_c proportional to exp(-delta HW_c).

What this script checks:
  1. Reproduce the published color-projected H_up/H_down matrices. This guards
     that the implemented 256-state operator is the one used in the corpus.
  2. Inspect the unprojected 9x9 quark-sector operator as 3 generations x
     3 colors.
  3. Test whether the diagonal color block for the Koide light node is rank-one
     while middle/heavy are color identity blocks.

Result:
  Negative for the documented 256-state operator. All diagonal generation color
  blocks are full-rank and not proportional to I_3. The only rank-one color
  object in the documented CKM construction is the external Boltzmann color
  singlet projection, and it is applied to every generation, not only the
  first-generation mass node.
"""

from __future__ import annotations

import cmath
import math


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
DELTA = 2 / 9
N_STATE = 256
GENS = [(0, 0), (0, 1), (1, 0)]  # Koide labels: heavy, light, middle.
COLORS = [(0, 1), (1, 0), (1, 1)]  # R, G, B in the framework's color-bit basis.
PAPER_ORDER = [0, 2, 1]  # Part 04 prints generations as 00,10,01.
KOIDE_LABELS = ["heavy n=0", "light n=1", "middle n=2"]

Matrix = list[list[complex]]


def hd(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def bits_to_index(bits: list[int] | tuple[int, ...]) -> int:
    out = 0
    for bit in bits:
        out = (out << 1) | bit
    return out


def index_to_bits(index: int) -> tuple[int, ...]:
    return tuple((index >> (7 - bit)) & 1 for bit in range(8))


def cnot_k_state(index: int, k: int) -> int:
    bits = list(index_to_bits(index))
    ctrl = (2 - k) % 8
    target = (5 - k) % 8
    if bits[ctrl] == 1:
        bits[target] ^= 1
    return bits_to_index(bits)


def zeros(n: int, m: int) -> Matrix:
    return [[0j for _ in range(m)] for _ in range(n)]


def dagger(matrix: Matrix) -> Matrix:
    return [[matrix[j][i].conjugate() for j in range(len(matrix))] for i in range(len(matrix[0]))]


def matmul(a: Matrix, b: Matrix) -> Matrix:
    n = len(a)
    p = len(b)
    m = len(b[0])
    out = zeros(n, m)
    for i in range(n):
        for k in range(p):
            aik = a[i][k]
            if abs(aik) < 1e-14:
                continue
            for j in range(m):
                out[i][j] += aik * b[k][j]
    return out


def build_walk() -> Matrix:
    amplitudes = [0j for _ in range(8)]
    amplitudes[0] = math.sqrt(1 - DELTA)
    for k in range(1, 8):
        amplitudes[k] = math.sqrt(DELTA / 7) * cmath.exp(1j * k * math.pi / 4)

    walk = zeros(N_STATE, N_STATE)
    for source in range(N_STATE):
        for k, amplitude in enumerate(amplitudes):
            walk[cnot_k_state(source, k)][source] += amplitude
    return walk


def quark_basis(isospin: int) -> list[int]:
    basis = []
    for gen in GENS:
        for color in COLORS:
            state = [0 for _ in range(8)]
            state[G0], state[G1] = gen
            state[LQ] = 1
            state[C0], state[C1] = color
            state[I3] = isospin
            state[CHI] = 0
            state[W] = 0
            basis.append(bits_to_index(state))
    return basis


def restrict(matrix: Matrix, basis: list[int]) -> Matrix:
    return [[matrix[row][col] for col in basis] for row in basis]


def color_boltzmann_vector() -> list[float]:
    # R=(0,1), G=(1,0), B=(1,1) have color-bit Hamming weights 1,1,2.
    raw = [math.exp(-DELTA), math.exp(-DELTA), math.exp(-2 * DELTA)]
    norm = math.sqrt(sum(value * value for value in raw))
    return [value / norm for value in raw]


def project_color_to_generation(h9: Matrix) -> Matrix:
    w = color_boltzmann_vector()
    out = zeros(3, 3)
    for gen_a in range(3):
        for gen_b in range(3):
            out[gen_a][gen_b] = sum(
                w[color_a]
                * h9[3 * gen_a + color_a][3 * gen_b + color_b]
                * w[color_b]
                for color_a in range(3)
                for color_b in range(3)
            )
    return out


def reorder_generation(matrix: Matrix, order: list[int]) -> Matrix:
    return [[matrix[i][j] for j in order] for i in order]


def matrix_rank(matrix: Matrix, tol: float = 1e-10) -> int:
    rows = [row[:] for row in matrix]
    n = len(rows)
    m = len(rows[0])
    rank = 0
    for col in range(m):
        if rank >= n:
            break
        pivot = max(range(rank, n), key=lambda row: abs(rows[row][col]))
        if abs(rows[pivot][col]) < tol:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(n):
            if row == rank or abs(rows[row][col]) < tol:
                continue
            factor = rows[row][col]
            rows[row] = [rows[row][j] - factor * rows[rank][j] for j in range(m)]
        rank += 1
    return rank


def color_block(h9: Matrix, gen_a: int, gen_b: int) -> Matrix:
    return [
        [h9[3 * gen_a + color_a][3 * gen_b + color_b] for color_b in range(3)]
        for color_a in range(3)
    ]


def identity_residual(block: Matrix) -> tuple[float, complex]:
    scalar = sum(block[i][i] for i in range(3)) / 3
    residual = max(
        abs(block[i][j] - (scalar if i == j else 0))
        for i in range(3)
        for j in range(3)
    )
    return residual, scalar


def fmt_complex(value: complex) -> str:
    if abs(value.imag) < 5e-4:
        return f"{value.real:.3f}"
    return f"{value.real:.3f}{value.imag:+.3f}i"


def print_matrix(title: str, matrix: Matrix) -> None:
    print(f"\n  {title}")
    for row in matrix:
        print("    " + "  ".join(f"{fmt_complex(value):>13s}" for value in row))


def assert_published_matrix_match(up: Matrix, down: Matrix) -> None:
    # Part 04 values after Boltzmann projection and paper generation ordering.
    expected_up = [
        [1.071, 0.051, 0.002],
        [0.051, 1.086, 0.000],
        [0.002, 0.000, 1.199],
    ]
    expected_down_abs = [
        [1.169, 0.034, 0.002],
        [0.034, 1.152, 0.000],
        [0.002, 0.000, 1.249],
    ]
    max_up = max(abs(up[i][j].real - expected_up[i][j]) for i in range(3) for j in range(3))
    max_down = max(abs(abs(down[i][j]) - expected_down_abs[i][j]) for i in range(3) for j in range(3))
    phase = math.atan2(down[0][1].imag, down[0][1].real) / math.pi
    print(f"\n  max residual vs Part 04 H_up real entries       = {max_up:.3e}")
    print(f"  max residual vs Part 04 |H_down| entries       = {max_down:.3e}")
    print(f"  down-sector phase arg(H_12)/pi                 = {phase:.3f}")
    assert max_up < 5e-3
    assert max_down < 5e-3
    assert 0.80 < abs(phase) < 0.90


def audit_color_blocks(name: str, h9: Matrix) -> bool:
    print(f"\n  {name}: diagonal generation color blocks in Koide order")
    all_match_desired = True
    for gen_index, label in enumerate(KOIDE_LABELS):
        block = color_block(h9, gen_index, gen_index)
        rank = matrix_rank(block)
        residual, scalar = identity_residual(block)
        desired = "rank-one" if gen_index == 1 else "identity"
        desired_ok = (rank == 1) if gen_index == 1 else residual < 1e-10
        all_match_desired = all_match_desired and desired_ok
        print(
            f"    {label:<10s} gen={GENS[gen_index]} desired={desired:<9s} "
            f"rank={rank}  tau=Tr/3={fmt_complex(scalar):>8s}  "
            f"identity_residual={residual:.3e}  "
            f"{'OK' if desired_ok else 'FAIL'}"
        )
    return all_match_desired


def main() -> None:
    hd("A. Build the Documented 256-State Walk Operator")
    walk = build_walk()
    m1 = matmul(dagger(walk), walk)
    m2 = matmul(m1, m1)
    print("  Built U = sum_k A_k CNOT^(k) on all 256 bit states.")
    print("  Built M1 = U^dagger U and M2 = (U^dagger U)^2.")

    up9 = restrict(m2, quark_basis(isospin=0))
    down9 = restrict(m2, quark_basis(isospin=1))
    up3 = reorder_generation(project_color_to_generation(up9), PAPER_ORDER)
    down3 = reorder_generation(project_color_to_generation(down9), PAPER_ORDER)

    hd("B. Guardrail: Reproduce the Published Color-Projected Matrices")
    print(f"  Boltzmann color vector w = {[round(x, 6) for x in color_boltzmann_vector()]}")
    print_matrix("H_up projected to generation space, Part 04 order", up3)
    print_matrix("H_down projected to generation space, Part 04 order", down3)
    assert_published_matrix_match(up3, down3)

    hd("C. Unprojected Color-Rank Test")
    up_ok = audit_color_blocks("up-left sector", up9)
    down_ok = audit_color_blocks("down-left sector", down9)
    if up_ok or down_ok:
        raise SystemExit("unexpected desired rank pattern appeared; inspect operator convention")

    hd("D. Verdict")
    print(
        "  The desired pattern is absent in the documented 256-state operator:\n"
        "    C_light is not rank-one; it is full rank.\n"
        "    C_middle and C_heavy are not proportional to I_3; they are color-anisotropic.\n\n"
        "  The rank-one color object used by Part 04 is the external Boltzmann\n"
        "  color-singlet vector |w><w|, applied uniformly to all generations when\n"
        "  projecting the 9-state quark space down to a 3-state generation matrix.\n"
        "  It is not a first-generation-only mass-node block.\n\n"
        "  Therefore the requested derivation cannot be made from this actual\n"
        "  256-state walk/mass operator. To get C_light=rank-one and\n"
        "  C_middle=C_heavy=I_3, a new color-resolved Yukawa operator must be\n"
        "  specified beyond the scalar Koide operator and the color-projected CKM\n"
        "  operator documented in the repo."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
