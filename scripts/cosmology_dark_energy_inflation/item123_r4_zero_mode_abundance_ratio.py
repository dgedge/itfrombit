#!/usr/bin/env python3
r"""ITEM 123: R4 zero-mode abundance ratio from finite support incidence.

The conserved-reservoir lift makes the CMB dust component mathematically
possible but leaves the conserved total N_tot free.  This script tests the
only canon-native relative-abundance hook visible in the finite R4 support.

Finite R4 support:
  * three sterile-nu_R forbidden source corners, one per R1 generation sector;
  * six undirected legal repair edges, two per source corner;
  * a service reservoir is directed, because a record has a billed entry and
    a billed exit orientation.

Therefore:

    directed R4 service-edge records / sterile nu_R source records
      = (2 * 6) / 3 = 4.

If the conserved zero-mode reservoir counts directed R4 service-edge records
while the sterile relic counts the three source corners, the R4 reservoir
abundance is fixed relative to the nu_R anchor:

    omega_zero = 4 omega_nuR.

With the existing cold nu_R anchor omega_nuR h^2 = 0.024, this gives
omega_zero h^2 = 0.096 and restores z_eq.  Equivalently, if the total cold dark
budget is taken from the CMB as 0.120, the finite support ratio partitions it
as 1/5 sterile and 4/5 R4 zero-mode.

This advances the abundance clause from "free number" to "relative split
derived under the directed-service-reservoir reading."  It still does NOT
derive the absolute total dark density; that remains either the observed CMB
budget or a future boot/initial-condition theorem.
"""

from __future__ import annotations

from fractions import Fraction

from item131_r4_support_dimension import CHI, I3, W, flip, gen, r4, species, valid_active, valid_r123
from r4_eos_cmb_resolution import (
    OMEGA_B_H2,
    OMEGA_C_H2_NEEDED,
    OMEGA_NUR_H2,
    OMEGA_R4_H2,
)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def forbidden_and_edges():
    forbidden = []
    for g0, g1 in [(0, 0), (0, 1), (1, 0)]:
        c = [0] * 8
        # (G0, G1, LQ, C0, C1, I3, chi, W)
        c[0], c[1] = g0, g1
        c[2], c[3], c[4], c[5], c[6], c[7] = 0, 0, 0, 0, 1, 1
        q = tuple(c)
        assert valid_r123(q) and not r4(q)
        forbidden.append(q)

    edges = []
    for q in forbidden:
        for label, target in (("I3", flip(q, I3)), ("chi/W", flip(q, CHI, W))):
            assert valid_active(target)
            edges.append((q, target, label))
    return forbidden, edges


def z_eq(omega_m_h2: float) -> float:
    omega_gamma_h2 = 2.469e-5
    n_eff = 3.044
    neutrino_factor = 1.0 + (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * n_eff
    omega_r_h2 = omega_gamma_h2 * neutrino_factor
    return omega_m_h2 / omega_r_h2 - 1.0


def main() -> None:
    print("ITEM 123: R4 ZERO-MODE ABUNDANCE RATIO")
    print("=" * 90)

    forbidden, edges = forbidden_and_edges()
    n_sources = len(forbidden)
    n_undirected_edges = len(edges)
    n_directed_edges = 2 * n_undirected_edges
    ratio = Fraction(n_directed_edges, n_sources)

    print("\n[1] Finite R4 support count")
    print(f"  sterile nu_R source corners      = {n_sources}")
    print(f"  legal undirected R4 repair edges = {n_undirected_edges}")
    print(f"  directed service-edge records    = {n_directed_edges}")
    print(f"  directed edge/source ratio       = {ratio}")
    check(n_sources == 3, "R4 has one forbidden sterile source per R1 generation sector")
    check(n_undirected_edges == 6, "R4 has two legal repair edges per source")
    check(ratio == 4, "directed R4 service-edge records give a 4:1 ratio to sterile sources")

    print("\n[2] Species and generation distribution")
    by_generation = {}
    for q in forbidden:
        by_generation[gen(q)] = []
    for q, target, label in edges:
        by_generation[gen(q)].append((label, species(target)))
    for generation, rows in sorted(by_generation.items()):
        print(f"  gen={generation}: source={species(next(q for q in forbidden if gen(q) == generation))}")
        for label, target_species in rows:
            print(f"    directed pair {label:5s}: nu_R <-> {target_species}")
    check(all(len(rows) == 2 for rows in by_generation.values()), "each generation has exactly two repair channels")

    print("\n[3] Relative abundance consequence")
    omega_from_nur = float(ratio) * OMEGA_NUR_H2
    total_from_nur = OMEGA_NUR_H2 + omega_from_nur
    sterile_fraction = OMEGA_NUR_H2 / total_from_nur
    r4_fraction = omega_from_nur / total_from_nur
    print(f"  omega_nuR h^2             = {OMEGA_NUR_H2:.4f}")
    print(f"  omega_zero h^2 = 4 nuR    = {omega_from_nur:.4f}")
    print(f"  total dark h^2            = {total_from_nur:.4f}")
    print(f"  fractions                 = sterile {sterile_fraction:.3f}, R4 zero-mode {r4_fraction:.3f}")
    check(abs(omega_from_nur - OMEGA_R4_H2) < 1.0e-12, "4:1 ratio maps the existing nu_R anchor to omega_zero=0.096")
    check(abs(total_from_nur - OMEGA_C_H2_NEEDED) < 1.0e-12, "nu_R plus directed R4 zero-mode gives omega_c=0.120")
    check(abs(sterile_fraction - 0.2) < 1.0e-12 and abs(r4_fraction - 0.8) < 1.0e-12, "the 20/80 split follows from 1 source : 4 directed service records")

    print("\n[4] Equality check")
    zeq_without = z_eq(OMEGA_B_H2 + OMEGA_NUR_H2)
    zeq_with = z_eq(OMEGA_B_H2 + OMEGA_NUR_H2 + omega_from_nur)
    print(f"  z_eq without zero-mode = {zeq_without:.0f}")
    print(f"  z_eq with zero-mode    = {zeq_with:.0f}")
    check(1000 < zeq_without < 1250, "without zero-mode equality is at recombination")
    check(3300 < zeq_with < 3500, "with ratio-fixed zero-mode equality is LCDM-like")

    print("\n[5] Tier boundary")
    print("  Derived here: the relative 4:1 R4-zero / nu_R split under the")
    print("  directed-service-reservoir reading.")
    print("  Not derived here: the absolute total dark density.  The script uses either")
    print("  the existing nu_R cold anchor omega_nuR=0.024 or the observed CMB total")
    print("  omega_c=0.120.  A full prediction still needs a boot theorem fixing the")
    print("  sterile source density or total dark sector density.")
    print("exit 0 -- relative abundance ratio closes conditionally; absolute density remains open.")


if __name__ == "__main__":
    main()
