#!/usr/bin/env python3
r"""Item 87 -- explicit nonlinear R1 monitor and fourth-generation boundary.

Question
--------
Can the R1 "no fourth generation" rule be represented as a monitored QEC
boundary rather than only as a static missing label?

Construction
------------
Generation register basis:

    tau = 00, e = 01, mu = 10, F = 11.

R1 forbids F.  The nonlinear check is implemented by two linear bit reads
(G0,G1) followed by a classical AND decoder.  The monitored recovery channel is

    P_R1 = I - |F><F|,
    K_e  = |e><F| / sqrt(2),
    K_mu = |mu><F| / sqrt(2).

This is CPTP and identity on the valid three-dimensional subspace.  The
environment records the repair choice.  If the environment records the closed
Hasse edge in the allowed order ideal, the two records are

    tau-e and tau-mu,

and their unsigned incidence covariance is K_R1 = B B^T.  The oriented edge
record is the antisymmetric R1 cochain.

Status boundary
---------------
The monitor closes gates 1--3.  Gate 4 is stricter: R1 alone is sector-blind and
therefore cannot derive d_s=(2,3,1,2).  Those counts reproduce the Koide
ellipticities only if the independently-derived sector defect counts are
identified with R1 rescue-contact multiplicities.  This script keeps that
conditional reproduction separate from strict derivation.

Exit 0 means the explicit channel is valid and the residual is exactly named.
"""

from __future__ import annotations

from fractions import Fraction


VALID = ("tau", "e", "mu")
ALL = ("tau", "e", "mu", "F")
INDEX = {name: i for i, name in enumerate(ALL)}
SECTORS = {
    "charged_lepton": (2, 9),
    "neutrino": (3, 9),
    "down": (1, 9),
    "up": (2, 27),
}


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def zero(n: int, m: int):
    return [[Fraction(0) for _ in range(m)] for _ in range(n)]


def identity(n: int):
    return [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(row) for row in zip(*A)]


def matmul(A, B):
    return [
        [sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))]
        for i in range(len(A))
    ]


def add(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def scalar_mul(c, A):
    return [[c * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def trace(A):
    return sum(A[i][i] for i in range(len(A)))


def projector_valid():
    P = zero(4, 4)
    for name in VALID:
        P[INDEX[name]][INDEX[name]] = Fraction(1)
    return P


def kraus_endpoint(endpoint: str):
    # Store K^dag K with the 1/2 weight; amplitudes are never needed explicitly.
    K = zero(4, 4)
    K[INDEX[endpoint]][INDEX["F"]] = Fraction(1, 1)  # symbolic |endpoint><F|
    return K


def dagger_times(K):
    return matmul(transpose(K), K)


def restrict_valid(A):
    return [[A[i][j] for j in range(3)] for i in range(3)]


def channel_on_basis_pair(a: str, b: str):
    """Return E(|a><b|) for the explicit channel, as a 4x4 matrix."""
    rho = zero(4, 4)
    rho[INDEX[a]][INDEX[b]] = Fraction(1)
    P = projector_valid()
    Ke = kraus_endpoint("e")
    Km = kraus_endpoint("mu")
    # Endpoint Kraus amplitudes are 1/sqrt(2), so K rho K^T carries 1/2.
    return add(
        matmul(matmul(P, rho), P),
        add(
            scalar_mul(Fraction(1, 2), matmul(matmul(Ke, rho), transpose(Ke))),
            scalar_mul(Fraction(1, 2), matmul(matmul(Km, rho), transpose(Km))),
        ),
    )


def eplane_projector():
    I = identity(3)
    J = [[Fraction(1) for _ in range(3)] for _ in range(3)]
    return add(I, scalar_mul(Fraction(-1, 3), J))


def eplane_invariants(K):
    P = eplane_projector()
    M = matmul(matmul(P, K), P)
    return trace(M), trace(matmul(M, M))


def eplane_equiv(A, B):
    P = eplane_projector()
    return matmul(matmul(P, A), P) == matmul(matmul(P, B), P)


def ellipticity_from_scaled_K(d: int, N: int):
    # K_R1 has E-plane eigenvalues {1,1/3}; isotropic completion has coefficient b.
    b = Fraction(N - 2 * d, 3)
    hi = Fraction(d) + b
    lo = Fraction(d, 3) + b
    return (hi - lo) / (hi + lo)


def main() -> None:
    print("ITEM 87 -- explicit nonlinear R1 monitor / fourth-generation boundary")
    print("=" * 82)

    print("\n[1] Physical channel: monitored fourth corner, valid subspace untouched")
    P = projector_valid()
    Ke = kraus_endpoint("e")
    Km = kraus_endpoint("mu")
    completeness = add(P, add(scalar_mul(Fraction(1, 2), dagger_times(Ke)), scalar_mul(Fraction(1, 2), dagger_times(Km))))
    check("P_R1 + endpoint recoveries form a CPTP channel", completeness == identity(4))
    for a in VALID:
        for b in VALID:
            out = channel_on_basis_pair(a, b)
            expected = zero(4, 4)
            expected[INDEX[a]][INDEX[b]] = Fraction(1)
            check(f"valid operator |{a}><{b}| is unchanged", out == expected)
    out_F = channel_on_basis_pair("F", "F")
    check("forbidden population is removed in one tick", out_F[INDEX["F"]][INDEX["F"]] == 0)
    check("forbidden state is rescued to endpoint mixture e/mu", out_F[INDEX["e"]][INDEX["e"]] == Fraction(1, 2) and out_F[INDEX["mu"]][INDEX["mu"]] == Fraction(1, 2))

    print("\n[2] Environment covariance: closed Hasse-edge record")
    # Rows in valid order tau,e,mu.  Columns are closed cover records tau-e and tau-mu.
    B_unsigned = [
        [Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
    ]
    K_r1 = matmul(B_unsigned, transpose(B_unsigned))
    print(f"    K_R1 = B B^T = {K_r1}")
    check("closed Hasse-edge covariance equals B B^T", K_r1 == [[2, 1, 1], [1, 1, 0], [1, 0, 1]])
    check("K_R1 E-plane eigenvalues are {1,1/3}", eplane_invariants(K_r1) == (Fraction(4, 3), Fraction(10, 9)))

    # Literal forbidden-corner rescue endpoint record is not identical, but it has the same E-plane anisotropy.
    C_rescue = [
        [Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(1)],
    ]
    check("literal forbidden-corner endpoint record is not literally K_R1", C_rescue != K_r1)
    check("literal rescue record is E-plane-equivalent to the closed Hasse record up to isotropic/singlet terms", eplane_invariants(C_rescue) == eplane_invariants(K_r1))

    print("\n[3] Oriented environment record: R1 cochain")
    # Oriented cover-edge adjacency on tau--e and tau--mu.  Endpoint swap flips this cochain.
    omega = [
        [Fraction(0), Fraction(1), Fraction(-1)],
        [Fraction(-1), Fraction(0), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0)],
    ]
    omega_swap = [
        [omega[0][0], omega[0][2], omega[0][1]],
        [omega[2][0], omega[2][2], omega[2][1]],
        [omega[1][0], omega[1][2], omega[1][1]],
    ]
    print(f"    Omega_R1 = {omega}")
    check("oriented cochain is antisymmetric", omega == scalar_mul(Fraction(-1), transpose(omega)))
    check("endpoint swap sends Omega_R1 -> -Omega_R1", omega_swap == scalar_mul(Fraction(-1), omega))
    check("oriented record is invisible to the symmetric covariance", all((add(omega, transpose(omega))[i][j] == 0) for i in range(3) for j in range(3)))

    print("\n[4] Sector contact counts: reproduction vs derivation")
    strict_r1_contact_count = 2  # two Hasse cover records for every sector.
    print(f"    strict R1-only contact count = {strict_r1_contact_count} for every sector")
    strict_derives = all(d == strict_r1_contact_count for d, _N in SECTORS.values())
    check("strict R1-only monitor is sector-blind and cannot derive all d_s", not strict_derives)
    for name, (d, N) in SECTORS.items():
        eps = ellipticity_from_scaled_K(d, N)
        b = Fraction(N - 2 * d, 3)
        print(f"    {name:15s}: imported contacts A_s=d_s={d}, isotropic B_s={b}, epsilon={eps}")
        check(f"{name}: conditional reproduction gives epsilon=d/N", eps == Fraction(d, N))
        check(f"{name}: isotropic completion nonnegative", b >= 0)

    print(
        """
VERDICT:
  The nonlinear R1 monitor can be built explicitly.  It treats the fourth corner
  11 as a monitored forbidden state, not a missing label: the channel is CPTP,
  identity on the valid generation subspace, and rescues 11 to the two one-bit
  endpoints.  This supports the stronger reading of "no fourth generation" as an
  actively maintained QEC boundary.

  Gates 1--3 close at finite-channel grade:
    1. the physical channel is identity/generation-blind on the valid subspace;
    2. closed Hasse-edge environment records give K_R1 = B B^T;
    3. the oriented edge record gives the antisymmetric R1 cochain.

  Gate 4 does NOT strictly close.  R1 alone has the same two cover records in
  every sector, so it cannot derive d_s=(2,3,1,2).  The four Koide ellipticities
  are reproduced only after importing the independently-derived sector defect
  counts as R1 contact multiplicities, A_s=d_s, with isotropic completion.

  Net result:
    * fourth-generation exclusion is upgraded from static absence to monitored
      boundary candidate;
    * virtual fourth-corner algebraic completions are allowed as off-boundary
      recovery excursions, but not as asymptotic particles;
    * Koide/CP unification remains conditional on the defect<->R1-contact theorem.
exit 0"""
    )


if __name__ == "__main__":
    main()
