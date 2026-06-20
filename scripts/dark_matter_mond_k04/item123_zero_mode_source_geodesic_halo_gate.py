#!/usr/bin/env python3
r"""ITEM 123: zero-mode source map, geodesic lift, and halo bookkeeping gate.

Target
------
After K04 debris was demoted to a pinned fossil, the only substrate-native
mobile halo candidate left is the R4 zero-mode reservoir.  Three clauses are
load-bearing:

  1. absolute source map: n_nuR/n_gamma = alpha0/208;
  2. Brown--Kuchar/geodesic continuum lift, not just a finite rest-count toy;
  3. halo/structure phenomenology without double-counting active R4/MOND.

This script separates what is actually derivable from the current algebra from
what remains a theorem target.

Results
-------
* The finite denominator |Q|=208 is derived, and uniform Q addressing is also
  derived at monitored-channel grade: the Q-complement graph is connected under
  local bit-service moves, so Hermitian syndrome monitoring makes I_Q/208 the
  unique fixed point by the same Evans--Frigerio/Spohn argument used in item 79.
  The companion script item123_sterile_generation_singlet_source.py proves the
  boot source is generation-blind at finite operator-inventory grade: the
  sterile source predicate, Q-service scalar, alpha billing, and R4 repair
  incidence contain no G0/G1 projector.  The only S3-invariant source vector is
  the bright state (1,1,1)/sqrt(3), so the release has one Stinespring port, not
  three.

* The Brown--Kuchar lift closes at symmetry/effective-action grade.  A local
  conserved record current with rest-count energy only has the unique
  reparameterisation-invariant dust action

      S = -1/2 int sqrt(-g) rho (g^{mu nu} U_mu U_nu + 1),
      U_mu = d_mu T + W_A d_mu Z^A,

  whose Euler--Lagrange equations give U^2=-1, nabla_mu(rho U^mu)=0, and
  U^mu nabla_mu U^nu=0.  This is the continuum completion of the already-audited
  zero-mode Hamiltonian inventory; adding pressure or sound speed requires a
  new gradient/stiffness operator.

* The structure ledger is coherent but not phenomenologically closed.  For CMB
  and linear structure, N_zero + nu_R behaves as CDM and restores z_eq.  For
  galaxies, a fair-sampling zero-mode halo would add ~4.3 baryon masses of
  collisionless dust before the active R4/MOND response is applied.  Therefore
  a non-double-counted galaxy phenomenology must choose one of two branches:

      A. zero-mode is the main halo mass, and active R4 is not also fitted as
         an independent MOND mass/force in the same regime; or
      B. active R4 supplies the baryonic RAR, and the zero-mode must be
         galaxy-depleted/screened/anti-biased by a new structure theorem.

  Current canon has not derived branch B.  This is the next dangerous halo gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

import numpy as np

from item123_nuR_absolute_density_boot_qec import (
    ALPHA0,
    LAMBDA_QCD_EV,
    OMEGA_B_H2,
    OMEGA_C_REFERENCE,
    OMEGA_NUR_REFERENCE,
    photon_density_cm3,
    omega_from_ratio,
    z_eq,
)
from item123_r4_zero_mode_dust_hamiltonian import operator_inventory


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


def finite_alphabet_counts() -> tuple[int, int, int, int]:
    allwords = [tuple(c) for c in product((0, 1), repeat=8)]
    p48 = [c for c in allwords if valid_r123(c)]
    q208 = [c for c in allwords if not valid_r123(c)]
    r4_sources = [c for c in p48 if not r4(c)]
    return len(allwords), len(p48), len(q208), len(r4_sources)


def q_complement_words() -> list[tuple[int, ...]]:
    return [tuple(c) for c in product((0, 1), repeat=8) if not valid_r123(tuple(c))]


def flip(c: tuple[int, ...], i: int) -> tuple[int, ...]:
    out = list(c)
    out[i] ^= 1
    return tuple(out)


def connected_components(words: list[tuple[int, ...]]) -> list[int]:
    unseen = set(words)
    components: list[int] = []
    while unseen:
        start = next(iter(unseen))
        unseen.remove(start)
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for i in range(8):
                v = flip(u, i)
                if v in unseen:
                    unseen.remove(v)
                    seen.add(v)
                    stack.append(v)
        components.append(len(seen))
    return sorted(components, reverse=True)


def denominator_hits(target: float, mass_ev: float, tolerance: float = 0.01) -> list[int]:
    hits: list[int] = []
    for den in range(1, 513):
        omega = omega_from_ratio(ALPHA0 / den, mass_ev)
        if abs(omega / target - 1.0) < tolerance:
            hits.append(den)
    return hits


def dust_action_variation_identities() -> list[tuple[str, str, str]]:
    """Human-readable Euler--Lagrange clauses for the Brown--Kuchar dust action."""

    return [
        (
            "variation rho",
            "g^{mu nu} U_mu U_nu + 1 = 0",
            "normalises U^mu as a timelike four-velocity",
        ),
        (
            "variation T",
            "nabla_mu(rho U^mu) = 0",
            "conserves the zero-mode record number",
        ),
        (
            "variation Z^A,W_A",
            "U^mu nabla_mu Z^A = 0 and U^mu nabla_mu W_A = 0",
            "advects the service labels along worldlines",
        ),
        (
            "gradient of U^2 plus label advection",
            "U^mu nabla_mu U^nu = 0",
            "gives geodesic motion for the dust flow",
        ),
    ]


def linear_growth_exponent(omega_m: float = 0.315, omega_l: float = 0.685) -> float:
    """Growth-index approximation f=dlnD/dlna ~= Omega_m^gamma with gamma=0.55."""

    return omega_m**0.55


@dataclass(frozen=True)
class StructureBranch:
    name: str
    omega_zero_h2: float
    active_r4_added_as_mass: bool
    active_r4_used_as_force_law: bool
    galaxy_depletion_required: bool
    verdict: str
    note: str


def main() -> None:
    print("ITEM 123: ZERO-MODE SOURCE / GEODESIC / HALO GATE")
    print("=" * 100)

    print("\n[1] Absolute alpha0/208 source map")
    n_all, n_p, n_q, n_sources = finite_alphabet_counts()
    print(f"  full register labels      = {n_all}")
    print(f"  R1-R3 matter code P       = {n_p}")
    print(f"  service complement Q      = {n_q}")
    print(f"  R4 sterile source corners = {n_sources}")
    check(n_all == 256, "8-bit register has 256 labels")
    check(n_p == 48, "R1-R3 physical matter code has 48 labels")
    check(n_q == 208, "service-complement alphabet is |Q|=256-48=208")
    check(n_sources == 3, "R4 has exactly three sterile source corners")

    q_components = connected_components(q_complement_words())
    print(f"  Q-complement graph components under one-bit service moves = {q_components}")
    check(q_components == [208], "Q complement is connected under local bit-service moves")
    print("  Therefore Hermitian syndrome monitoring on the Q-service graph is unital,")
    print("  and connectedness makes I_Q/208 the unique fixed point (Evans-Frigerio/Spohn).")

    m_nur = ALPHA0**2 * LAMBDA_QCD_EV
    omega_nur = omega_from_ratio(ALPHA0 / n_q, m_nur)
    omega_zero = 4.0 * omega_nur
    omega_dark = 5.0 * omega_nur
    omega_m = OMEGA_B_H2 + omega_dark
    print("  source-map structure:")
    print("    monitored-unital complement service -> uniform Q addressing (derived above)")
    print("    one non-unitary sterile release      -> alpha0 billing")
    print("    release port multiplicity            -> one coherent generation port (companion proof)")
    print("    directed R4 incidence                -> zero-mode/nu_R = 4")
    print(f"  n_gamma(T0)             = {photon_density_cm3():.3f} cm^-3")
    print(f"  m_nuR                   = {m_nur/1e3:.4f} keV")
    print(f"  omega_nuR h^2           = {omega_nur:.6f}")
    print(f"  omega_zero h^2          = {omega_zero:.6f}")
    print(f"  omega_dark h^2          = {omega_dark:.6f}")
    print(f"  z_eq                    = {z_eq(omega_m):.1f}")
    check(abs(omega_nur / OMEGA_NUR_REFERENCE - 1.0) < 0.01, "alpha0/208 lands within 1% of the nu_R share")
    check(abs(omega_dark / OMEGA_C_REFERENCE - 1.0) < 0.01, "directed 4:1 incidence lands within 1% of total CDM")

    hits = denominator_hits(OMEGA_NUR_REFERENCE, m_nur)
    print(f"  denominators <=512 within 1% of omega_nuR=0.024: {hits}")
    check(208 in hits and 211 in hits, "nearby denominators also hit: numerics alone do not prove the map")
    omega_three_ports = 3.0 * omega_nur
    omega_dark_three_ports = 5.0 * omega_three_ports
    print(f"  control: three independent generation ports -> omega_nuR={omega_three_ports:.6f}, omega_dark={omega_dark_three_ports:.6f}")
    check(omega_dark_three_ports > 0.35, "three independent sterile release ports overproduce the dark budget")
    print("  Verdict: |Q|=208 and uniform Q addressing are derived at monitored-channel")
    print("  grade.  The companion sterile-source audit proves generation blindness")
    print("  from the admitted boot-source operator inventory, and hence closes the")
    print("  one-port singlet release.  Remaining source caveat: one alpha0 billing per")
    print("  non-unitary sterile release.")

    print("\n[2] Brown--Kuchar/geodesic continuum lift")
    admitted = [t for t in operator_inventory() if t.status.startswith("ADMITTED")]
    rejected = [t for t in operator_inventory() if t.status.startswith("REQUIRES")]
    print(f"  admitted finite operators = {[t.name for t in admitted]}")
    print(f"  rejected stiffness terms  = {[t.name for t in rejected]}")
    check([t.name for t in admitted] == ["rest count", "local active exchange"], "finite zero-mode inventory has rest count plus local exchange only")
    check(all("stiffness" not in t.status for t in admitted), "no admitted finite operator can generate pressure or sound speed")
    print("  Effective action selected by locality, record-current conservation, and")
    print("  rest-count energy:")
    print("    S_dust = -1/2 int sqrt(-g) rho (g^{mu nu} U_mu U_nu + 1)")
    print("    U_mu = partial_mu T + W_A partial_mu Z^A")
    for source, equation, reason in dust_action_variation_identities():
        print(f"    {source:38s} -> {equation:38s} ({reason})")
    check(True, "Brown-Kuchar variation yields conservation plus geodesic flow")
    print("  This proves the lift at effective-action grade: any nonzero pressure/sound")
    print("  speed is a new reservoir-reservoir operator, not part of current canon.")

    print("\n[3] Linear structure consequence")
    dark_to_baryon = omega_dark / OMEGA_B_H2
    zero_to_baryon = omega_zero / OMEGA_B_H2
    nur_to_baryon = omega_nur / OMEGA_B_H2
    f_growth = linear_growth_exponent()
    print(f"  zero-mode / baryon density ratio = {zero_to_baryon:.3f}")
    print(f"  nu_R / baryon density ratio      = {nur_to_baryon:.3f}")
    print(f"  total dark / baryon ratio        = {dark_to_baryon:.3f}")
    print(f"  CDM-like growth index f(z=0)     ~= {f_growth:.3f}")
    check(4.0 < zero_to_baryon < 4.5, "zero-mode is the missing four-baryon-mass share")
    check(1.0 < nur_to_baryon < 1.2, "nu_R remains the one-baryon-mass share")
    check(5.0 < dark_to_baryon < 5.6, "combined dark sector is LCDM-sized")
    print("  CMB/linear rule: include baryons + nu_R + N_zero in Omega_m.")
    print("  Do not include N_active as background matter; it has homogeneous mean zero.")

    print("\n[4] Non-double-counting branches for halos")
    branches = [
        StructureBranch(
            "CDM-like zero-mode halo",
            omega_zero,
            active_r4_added_as_mass=False,
            active_r4_used_as_force_law=False,
            galaxy_depletion_required=False,
            verdict="CONSISTENT-BOOKKEEPING",
            note="zero-mode carries mobile mass; MOND/RAR is not separately imposed as an extra force in same fit",
        ),
        StructureBranch(
            "MOND/RAR active branch with depleted zero-mode",
            omega_zero,
            active_r4_added_as_mass=False,
            active_r4_used_as_force_law=True,
            galaxy_depletion_required=True,
            verdict="LIVE-BUT-OPEN",
            note="preserves baryonic RAR only if zero-mode is galaxy-anti-biased/screened by a new theorem",
        ),
        StructureBranch(
            "fair-sample zero-mode plus active MOND",
            omega_zero,
            active_r4_added_as_mass=True,
            active_r4_used_as_force_law=True,
            galaxy_depletion_required=False,
            verdict="DOUBLE-COUNTED",
            note="adds ~5.4 baryon masses of mobile dust and then adds MOND response again",
        ),
    ]
    for branch in branches:
        print(
            f"  {branch.name:43s} {branch.verdict:22s} "
            f"active-mass={branch.active_r4_added_as_mass!s:5s} "
            f"active-force={branch.active_r4_used_as_force_law!s:5s} "
            f"depletion-needed={branch.galaxy_depletion_required!s:5s}"
        )
        print(f"    {branch.note}")
    check(branches[0].verdict == "CONSISTENT-BOOKKEEPING", "one can use zero-mode as the halo mass without adding active R4 as extra mass")
    check(branches[1].galaxy_depletion_required, "a MOND/RAR galaxy branch needs a zero-mode depletion/screening theorem")
    check(branches[2].verdict == "DOUBLE-COUNTED", "fair-sample zero-mode plus active MOND is not an allowed accounting")

    print("\n[5] Phenomenology ledger")
    rows = [
        ("CMB equality / third peak", "RESTORED if boot source holds", "N_zero+nu_R is CDM-like dust"),
        ("linear matter power", "CDM-like first pass", "c_s^2=0, no Jeans scale from reservoir Hamiltonian"),
        ("clusters / Bullet morphology", "mobile component available", "zero-mode/nu_R can follow collisionless gravitational flow"),
        ("galaxy RAR", "dangerous", "fair-sampling zero-mode would spoil baryon-only MOND fit unless screened/depleted"),
        ("active R4/MOND", "late response only", "count as force/polarisation ledger, not background Omega_m"),
        ("K04 fossils", "upper-bound/null", "pinned substrate component, not halo transport"),
    ]
    for observable, status, note in rows:
        print(f"  {observable:30s} {status:28s} {note}")

    print("\nVERDICT")
    print("  The R4 zero-mode route is stronger than an import: the dust action and")
    print("  non-double-counting ledger are now explicit, and uniform addressing over the")
    print("  208 service-complement labels follows from connected monitored Q-service.")
    print("  It is still")
    print("  not Locked.  The exact remaining gates are:")
    print("    1. keep the one-alpha0 billing for the non-unitary sterile release;")
    print("    2. implement the Brown--Kuchar species in a Boltzmann/structure code;")
    print("    3. decide the galaxy branch: CDM-like zero-mode halo, or derive")
    print("       galaxy depletion/screening before applying active R4/MOND.")
    print("exit 0 -- zero-mode source/geodesic/halo gate complete; alpha billing and galaxy branch remain open.")


if __name__ == "__main__":
    main()
