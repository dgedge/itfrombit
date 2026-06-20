#!/usr/bin/env python3
r"""ITEM 123: absolute nu_R / total dark density from a boot-QEC source law.

The zero-mode reservoir scripts close the pressureless form and the 4:1
R4-zero/nu_R split, but they still use either the existing nu_R anchor
omega_nuR h^2=0.024 or the observed total omega_c h^2=0.120 to set the absolute
normalisation.  This script tests the only compact boot/QEC source law visible in
the finite register:

    n_nuR / n_gamma = alpha_0 / |Q|,
    |Q| = 2^8 - 48 = 208.

Interpretation
--------------
The boot service bath has a 208-label complement Q outside the R1-R3 matter
code.  A sterile source release is one non-unitary alpha_0-billed event addressed
through that complement.  This is NOT the old, failed "205 = 208 - 3" gravity
claim: the sterile nu_R corners are codewords in P, not Q, so we do not subtract
them from Q.  The proposed source law uses Q as the service-complement address
alphabet, with one alpha_0 suppression per release.

If this source map is correct, the absolute density follows without inserting
omega_nuR or omega_c:

    m_nuR = alpha_0^2 Lambda_QCD,
    n_nuR = (alpha_0/208) n_gamma(T_0),
    omega_nuR h^2 = m_nuR n_nuR / rho_crit,h^2,
    omega_zero = 4 omega_nuR       (directed R4 service-edge incidence),
    omega_dark = 5 omega_nuR.

Result
------
The candidate lands at omega_nuR h^2=0.02418 and omega_dark h^2=0.12089, a
0.7% high prediction relative to the standard 0.024/0.120 values.  The result is
strong enough to promote as a candidate absolute-normalisation theorem target.
It is not Locked: the load-bearing clause is the microscopic source map.  The
finite denominator is structurally selected by the finite code, and companion
scripts now derive monitored-uniform Q addressing plus a generation-blind
one-port singlet release from the finite source algebra.  What remains is the
microscopic boot assertion that the sterile release is one alpha-billed,
non-unitary event addressed through that Q service complement.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math


ALPHA0 = 1.0 / 137.035999
LAMBDA_QCD_EV = 0.332e9
T_CMB = 2.7255
ZETA3 = 1.2020569031595943
KB_EV_PER_K = 8.617333262e-5
HBARC_EV_CM = 1.973269804e-5
RHO_CRIT_H2_EV_CM3 = 1.05371e4

OMEGA_B_H2 = 0.0224
OMEGA_NUR_REFERENCE = 0.024
OMEGA_C_REFERENCE = 0.120
N_EFF = 3.044

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    if c[LQ] == 0:
        return (c[C0], c[C1]) == (0, 0)
    return (c[C0], c[C1]) != (0, 0)


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def species(c: tuple[int, ...]) -> str:
    if c[LQ] == 0 and (c[C0], c[C1]) == (0, 0):
        if c[I3] == 0 and c[CHI] == 0:
            return "nu_L"
        if c[I3] == 1 and c[CHI] == 0:
            return "e_L"
        if c[I3] == 0 and c[CHI] == 1:
            return "nu_R"
        if c[I3] == 1 and c[CHI] == 1:
            return "e_R"
    return "other"


def photon_density_cm3(t_cmb: float = T_CMB) -> float:
    return 2.0 * ZETA3 / math.pi**2 * (KB_EV_PER_K * t_cmb / HBARC_EV_CM) ** 3


def omega_from_ratio(n_over_ngamma: float, mass_ev: float) -> float:
    return n_over_ngamma * photon_density_cm3() * mass_ev / RHO_CRIT_H2_EV_CM3


def omega_r_h2() -> float:
    omega_gamma = 2.469e-5 * (T_CMB / 2.7255) ** 4
    neutrino_factor = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * N_EFF
    return omega_gamma * neutrino_factor


def z_eq(omega_m_h2: float) -> float:
    return omega_m_h2 / omega_r_h2() - 1.0


@dataclass(frozen=True)
class Clause:
    name: str
    status: str
    reason: str


def main() -> None:
    print("ITEM 123: ABSOLUTE nu_R / TOTAL DARK DENSITY FROM BOOT-QEC SOURCE LAW")
    print("=" * 100)

    print("\n[1] Finite code alphabet")
    allwords = [tuple(c) for c in product((0, 1), repeat=8)]
    p48 = [c for c in allwords if valid_r123(c)]
    q208 = [c for c in allwords if not valid_r123(c)]
    nur_sources = [c for c in p48 if not r4(c)]
    print(f"  full register labels       = {len(allwords)}")
    print(f"  R1-R3 matter codewords P   = {len(p48)}")
    print(f"  service-complement labels Q= {len(q208)}")
    print(f"  R4 sterile source corners  = {len(nur_sources)}")
    check(len(allwords) == 256, "8-bit register has 256 labels")
    check(len(p48) == 48, "R1-R3 matter code has 48 labels")
    check(len(q208) == 208, "complement Q has 208 labels")
    check(len(nur_sources) == 3 and all(species(c) == "nu_R" for c in nur_sources), "R4 excludes exactly three nu_R source corners in P")
    check(all(c not in q208 for c in nur_sources), "nu_R source corners are not Q labels; this is not the old 208-3 claim")

    print("\n[2] Candidate boot source law")
    n_ratio = ALPHA0 / len(q208)
    m_nur = ALPHA0**2 * LAMBDA_QCD_EV
    n_gamma = photon_density_cm3()
    n_nur = n_ratio * n_gamma
    print("  source law: n_nuR/n_gamma = alpha0 / |Q|")
    print(f"  alpha0                 = {ALPHA0:.9e}")
    print(f"  |Q|                    = {len(q208)}")
    print(f"  n_nuR/n_gamma          = {n_ratio:.9e}")
    print(f"  n_gamma(T0={T_CMB} K)  = {n_gamma:.3f} cm^-3")
    print(f"  n_nuR                  = {n_nur:.6f} cm^-3")
    print(f"  m_nuR=alpha0^2 Lambda  = {m_nur/1e3:.4f} keV")
    check(400.0 < n_gamma < 420.0, "photon density computed from T_CMB")
    check(17.5 < m_nur / 1e3 < 17.9, "sterile mass is the existing alpha^2 Lambda value")

    print("\n[3] Absolute density prediction")
    omega_nur = omega_from_ratio(n_ratio, m_nur)
    omega_zero = 4.0 * omega_nur
    omega_dark = 5.0 * omega_nur
    omega_m = OMEGA_B_H2 + omega_dark
    zeq = z_eq(omega_m)
    print(f"  omega_nuR h^2           = {omega_nur:.6f}")
    print(f"  omega_zero h^2 = 4 nuR  = {omega_zero:.6f}")
    print(f"  omega_dark h^2 = 5 nuR  = {omega_dark:.6f}")
    print(f"  omega_m h^2             = {omega_m:.6f}")
    print(f"  z_eq                    = {zeq:.1f}")
    print(f"  vs reference nuR 0.024  = {(omega_nur/OMEGA_NUR_REFERENCE-1)*100:+.2f}%")
    print(f"  vs reference CDM 0.120  = {(omega_dark/OMEGA_C_REFERENCE-1)*100:+.2f}%")
    check(abs(omega_nur / OMEGA_NUR_REFERENCE - 1.0) < 0.01, "alpha/208 predicts the sterile share within 1%")
    check(abs(omega_dark / OMEGA_C_REFERENCE - 1.0) < 0.01, "4:1 reservoir incidence then predicts total dark density within 1%")
    check(3400.0 < zeq < 3450.0, "predicted matter-radiation equality is LCDM-like")

    print("\n[4] Controls and competitor denominators")
    denominators = [45, 48, 137, 205, 208, 210, 256]
    for den in denominators:
        omega = omega_from_ratio(ALPHA0 / den, m_nur)
        print(f"  alpha/{den:<3d}: omega_nuR={omega:.5f}, omega_dark(5x)={5*omega:.5f}")
    required_den = ALPHA0 / (OMEGA_NUR_REFERENCE * RHO_CRIT_H2_EV_CM3 / (m_nur * n_gamma))
    simple_hits = [
        den
        for den in range(1, 513)
        if abs(omega_from_ratio(ALPHA0 / den, m_nur) / OMEGA_NUR_REFERENCE - 1.0) < 0.01
    ]
    print(f"  observational inverse denominator for omega_nuR=0.024: {required_den:.3f}")
    print(f"  integer denominators <=512 within 1% of the sterile target: {simple_hits}")
    check(208 in simple_hits, "the code-complement denominator 208 is in the 1% band")
    check(210 in simple_hits, "nearby integers also land: the number alone is not evidence")
    check(205 not in simple_hits, "the old 205-like value is close but outside the 1% band")
    print("  T8 reading: the numerical hit is not unique.  The claim stands or falls on")
    print("  the structural source map alpha0/|Q|, not on closeness to 0.024.")

    print("\n[5] Clause ledger")
    clauses = [
        Clause(
            "finite denominator",
            "DERIVED",
            "|Q|=256-48=208 from R1-R3 code complement",
        ),
        Clause(
            "alpha billing",
            "CANON-GROUNDED",
            "one non-unitary source release carries one alpha0 erasure weight",
        ),
        Clause(
            "source map",
            "COMPANION-ADVANCED",
            "uniform Q addressing plus generation-blind one-port singlet release; still conditional on alpha-billed boot release",
        ),
        Clause(
            "sterile mass",
            "CANON-INPUT",
            "m_nuR=alpha0^2 Lambda_QCD from item 118",
        ),
        Clause(
            "4:1 zero-mode split",
            "COMPANION-CLOSED",
            "directed R4 service-edge incidence gives omega_zero=4 omega_nuR",
        ),
        Clause(
            "absolute dark density",
            "CONDITIONAL-PREDICTED",
            "omega_dark h^2=0.12089 follows without inserting observed omega_c",
        ),
    ]
    for clause in clauses:
        print(f"  {clause.name:24s} {clause.status:22s} {clause.reason}")

    print("\n[6] Verdict")
    print("  Candidate closure: if the boot source law n_nuR/n_gamma=alpha0/208 is")
    print("  derived from the service-complement algebra, the absolute sterile density")
    print("  and the total dark density follow: omega_nuR h^2=0.02418 and")
    print("  omega_dark h^2=0.12089.  This uses the CMB photon density as the bath")
    print("  normalization, but does not use the existing 0.024 anchor or observed")
    print("  omega_c h^2 as inputs.")
    print("  Not Locked: the microscopic boot-rate theorem is still load-bearing.")
    print("  Companion scripts now derive uniform Q addressing and a generation-blind")
    print("  one-port singlet release.  The remaining assertion is one alpha0-billed")
    print("  non-unitary sterile release through that Q complement.")
    print("  The old 205=208-3 gravity claim stays dead; here 208 is the")
    print("  service-complement address alphabet, not a set containing the sterile states.")
    print("exit 0 -- alpha/208 is a strong conditional boot-QEC density theorem target.")


if __name__ == "__main__":
    main()
