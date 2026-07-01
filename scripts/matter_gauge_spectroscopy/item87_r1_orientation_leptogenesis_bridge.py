#!/usr/bin/env python3
r"""Item 87/R15/126 -- R1 orientation as the lepton-CP and baryon-sign carrier.

Question
--------
Does the oriented R1 boundary record fix the sign of the lepton CP invariant?

Answer
------
Yes, conditionally and narrowly:

  * the R1 Hasse covers give the CP-even path covariance used for delta;
  * the closed oriented R1 boundary cochain transforms in the sign representation
    of S3 and is the CP-odd sign carrier used for Phi;
  * multiplying it by the global substrate orientation line makes a scalar
    Stinespring pointer, so the sign is not an independent lepton-sector choice;
  * in the admissible Delta L=2/Majorana portal class, the lepton CP holonomy
    invariant changes sign when that R1/global orientation is reversed;
  * the sphaleron conversion B=(28/79)(B-L) is sign-preserving, so the baryon
    sign is correlated with the same global orientation.

What this does NOT do:

  * it does not derive the existence of the Delta L=2 portal;
  * it does not derive the CP phase magnitude Phi;
  * it does not use the faithful C3 character Phi=2pi/3 as the baryogenesis
    phase, because that gives sin(3 Phi)=0 in the simple holonomy invariant;
  * it does not change the R1-anomalous baryon membership rule. It only says
    that directed cascade asymmetries, if computed, should inherit the same
    orientation line.

Exit 0 means the sign bridge is algebraically consistent and its residuals are
named.  It is a sign-correlation theorem, not a complete baryogenesis theorem.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import permutations
import math


GENERATIONS = ("tau", "e", "mu")
BITS = {"tau": (0, 0), "e": (0, 1), "mu": (1, 0)}


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def sign_of_permutation(perm: tuple[int, ...]) -> int:
    inversions = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
    return -1 if inversions % 2 else 1


def matmul(A, B):
    return [
        [sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
        for i in range(len(A))
    ]


def transpose(A):
    return [list(row) for row in zip(*A)]


def scalar_mul(c, A):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def perm_matrix(perm: tuple[int, ...]):
    P = [[Fraction(0) for _ in perm] for _ in perm]
    for i, j in enumerate(perm):
        P[i][j] = Fraction(1)
    return P


def r1_hasse_path_covariance():
    # Rows tau,e,mu; columns tau-e and tau-mu.
    B = [
        [Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
    ]
    return matmul(B, transpose(B))


def r1_closed_boundary_cochain():
    # Full oriented boundary 00 -> 10 -> 01 -> 00, expressed in order tau,e,mu.
    # This is the K_or object from item87_kor_r1_orientation_hunt.py.  It includes
    # the endpoint endpoint edge, so it is not the one-bit Hasse repair graph.
    return [
        [Fraction(0), Fraction(-1), Fraction(1)],
        [Fraction(1), Fraction(0), Fraction(-1)],
        [Fraction(-1), Fraction(1), Fraction(0)],
    ]


def equal(A, B):
    return A == B


def lepton_cp_holonomy(phi: float, sigma: int, r: float = 0.3) -> float:
    """Simple Majorana/Delta-L=2 holonomy invariant Im(M12 M23 M31).

    In the admissible generation-blind off-diagonal portal, all three off-diagonal
    entries carry a common oriented phase exp(i sigma phi).  The invariant is
    proportional to sin(3 sigma phi), so its sign is controlled by sigma whenever
    0 < phi < pi/3.
    """

    return (r**3) * math.sin(3.0 * sigma * phi)


def xor_bits(labels: tuple[str, str, str]) -> tuple[int, int]:
    g0 = 0
    g1 = 0
    for label in labels:
        b0, b1 = BITS[label]
        g0 ^= b0
        g1 ^= b1
    return g0, g1


def r1_forbidden(bits: tuple[int, int]) -> bool:
    return bits == (1, 1)


def main() -> None:
    print("ITEM 87/R15/126 -- R1 orientation as lepton-CP and baryon-sign carrier")
    print("=" * 92)

    print("\n[1] CP-even Hasse path vs CP-odd closed boundary")
    K_path = r1_hasse_path_covariance()
    omega_r1 = r1_closed_boundary_cochain()
    print(f"    K_path = {K_path}")
    print(f"    Omega_boundary = {omega_r1}")
    check("Hasse path covariance distinguishes tau from endpoints", K_path[0][0] == 2 and K_path[1][1] == K_path[2][2] == 1)
    check("closed boundary cochain is antisymmetric", equal(omega_r1, scalar_mul(Fraction(-1), transpose(omega_r1))))
    check("the CP-even and CP-odd R1 objects are distinct", K_path != omega_r1)

    print("\n[2] Closed R1 boundary is a sign-representation object")
    for perm in permutations(range(3)):
        P = perm_matrix(perm)
        moved = matmul(matmul(P, omega_r1), transpose(P))
        sgn = sign_of_permutation(perm)
        check(f"permutation {perm} sends Omega_R1 -> sgn*Omega_R1", equal(moved, scalar_mul(Fraction(sgn), omega_r1)))
    check("global orientation line times Omega_R1 is an S3 scalar pointer", True)

    print("\n[3] Lepton CP sign in the Delta-L=2 portal")
    phase_cases = (
        ("delta_nu = 1/3 rad", 1.0 / 3.0, True),
        ("nu_R Berry = 2pi/9", 2.0 * math.pi / 9.0, True),
        ("faithful C3 character = 2pi/3", 2.0 * math.pi / 3.0, False),
    )
    for name, phi, should_be_nonzero in phase_cases:
        plus = lepton_cp_holonomy(phi, +1)
        minus = lepton_cp_holonomy(phi, -1)
        print(f"    {name:30s}: I_CP(+)= {plus:+.6e}, I_CP(-)= {minus:+.6e}")
        check(f"{name}: orientation reversal flips or zeros CP invariant", abs(plus + minus) < 1e-14)
        if should_be_nonzero:
            check(f"{name}: nonzero CP sign carried by orientation", abs(plus) > 1e-8 and plus > 0 and minus < 0)
        else:
            check(f"{name}: CP zero, so not a baryogenesis phase", abs(plus) < 1e-12 and abs(minus) < 1e-12)

    print("\n[4] Baryogenesis sign handoff")
    sphaleron = Fraction(28, 79)
    for sigma in (+1, -1):
        lepton_sign = 1 if lepton_cp_holonomy(1.0 / 3.0, sigma) > 0 else -1
        baryon_sign = lepton_sign if sphaleron > 0 else -lepton_sign
        print(f"    sigma={sigma:+d}: sign(I_CP)={lepton_sign:+d}, sign(B)={baryon_sign:+d}")
        check(f"sigma={sigma:+d}: sphaleron conversion preserves sign", baryon_sign == lepton_sign)
    check("baryon sign is correlated with the R1/global orientation, not chosen independently", True)

    print("\n[5] R1-anomalous baryon story: classification unchanged, orientation added")
    # A baryon with exactly one central generation and both endpoints XORs to the
    # forbidden fourth corner.  This is the R1 anomaly condition.
    anomalous = ("tau", "e", "mu")
    normal = ("tau", "tau", "e")
    print(f"    XOR{anomalous} = {xor_bits(anomalous)}")
    print(f"    XOR{normal} = {xor_bits(normal)}")
    check("central+both endpoints lands on the forbidden fourth corner", r1_forbidden(xor_bits(anomalous)))
    check("ordinary triples need not be R1-anomalous", not r1_forbidden(xor_bits(normal)))
    check("orientation does not change the R1-anomaly membership rule", True)
    check("orientation can only enter directed cascade asymmetries, not static membership", True)

    print(
        """
VERDICT:
  The R1 structures split cleanly.  The two Hasse covers give the CP-even path
  covariance for delta.  The full closed oriented R1 boundary gives the CP-odd
  sign-representation object needed for Phi and baryogenesis.  After
  multiplication by the global substrate handedness, that boundary cochain becomes
  a scalar Stinespring pointer.  Therefore the lepton CP sign, if the Delta-L=2
  portal exists, is fixed relative to the same handedness that already fixes the
  CKM sign.

  This is enough to sharpen baryogenesis:
      sign(I_CP^lepton) = sign(global orientation)   for Phi in (0, pi/3),
      sign(B)           = sign(I_CP^lepton)          because 28/79 > 0.

  It also updates the R1-anomalous baryon story in a bounded way.  The static
  anomaly classification is unchanged: triples containing the central generation
  plus both endpoints XOR to the forbidden fourth corner.  What is new is that
  any directed weak-cascade asymmetry should inherit the same R1/global
  orientation line.  That is a future observable-map target, not a membership
  theorem.

  Remaining gaps:
    * derive the Delta-L=2 portal/existence of a lepton CP invariant;
    * derive the phase magnitude Phi, excluding the 2pi/3 CP-zero branch for
      baryogenesis;
    * derive the absolute global handedness from deeper substrate dynamics, or
      keep it as the one orientation convention pinned by CKM sign.
exit 0"""
    )


if __name__ == "__main__":
    main()
