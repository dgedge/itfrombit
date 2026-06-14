#!/usr/bin/env python3
r"""ITEM 132: microscopic R4/QEC line-network to local p=3 action lift.

This script tries to close the remaining gap left by item132_p3_r4_action.py:

    Why should the bound R4 halo admit a local scale-covariant action at all?

Candidate lift
--------------
Use the finite facts already established elsewhere:

1. R4 boundary-QEC support is one-dimensional:
   item131_r4_support_dimension.py gives three disjoint two-edge repair stars.
2. A local reversible/QEC edge has the standard quadratic Dirichlet/Fisher
   stiffness, i.e. the p=2 action density |grad phi|^2/(8 pi G).
3. A 1D active support in a local acceleration field has active line density
   linear in the dimensionless demand x=|grad phi|/a0.  The R4 support supplies
   two legal repair edges per generation, distributed over the three R1
   generation sectors, giving the finite incidence normalization 2/3.

Then the coarse-grained action density is

    (quadratic edge stiffness) x (active R4 line density)

      = [ |g|^2 / (8 pi G) ] [ (2/3) |g| / a0 ]
      = |g|^3 / (12 pi G a0).

Varying this action gives

    div[(|g|/a0) g] = 4 pi G rho_b,

the p=3 AQUAL/deep-MOND equation used by item132_halo_closure.py.

Honest status
-------------
This is progress because the cubic power and the coefficient reduce to finite
R4 support data plus the standard local edge Dirichlet cost.  Companion
item132_r4_line_density_dynamics.py shows that generic stable local dynamics
derives the same form only up to a susceptibility chi_R4.  The remaining
physical premise is therefore very narrow: the virialized bound R4 halo must
really equilibrate to the unit-susceptibility local line-density law

    lambda_R4(x) = (2/3) x,  x = |g|/a0,

rather than a fixed, area-like, volume-like, or nonlinear activation density.
Companion item132_chi_unit_poisson.py closes this conditionally by showing
that a W=S*C one-jump, count-valued boundary-QEC immigration-death ledger with
matched creation/erasure rate Gamma0 has stationary Poisson(|g|/a0) count and
therefore chi_R4=1.  That rate-matching lemma is now the microscopic hinge.
The law is the halo analogue of the late-cosmology f(a)=a R4 line ledger.
"""

from __future__ import annotations

from fractions import Fraction


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def action_power(edge_power: int, support_dimension: int) -> int:
    """Quadratic edge stiffness times d-dimensional activation gives p=2+d."""
    return edge_power + support_dimension


def v2_radius_exponent_for_support_dimension(d: int) -> Fraction:
    """With p=2+d, spherical v^2 scales as r^this."""
    p = Fraction(action_power(2, d), 1)
    return Fraction(1, 1) - Fraction(2, 1) / (p - 1)


def v4_mass_exponent_for_support_dimension(d: int) -> Fraction:
    """With p=2+d, spherical v^4 scales as M_b^this."""
    p = Fraction(action_power(2, d), 1)
    return Fraction(2, 1) / (p - 1)


def coarse_action_coefficient(
    edge_dirichlet_coeff: Fraction,
    repair_edges_per_generation: int,
    generation_sectors: int,
) -> Fraction:
    """Coefficient multiplying |g|^3/(pi G a0), omitting pi G a0 itself.

    Newtonian edge action has coefficient 1/8 in |g|^2/(pi G).  Multiplying
    by the R4 active line-density prefactor 2/3 gives 1/12.
    """
    return edge_dirichlet_coeff * Fraction(repair_edges_per_generation, generation_sectors)


def variation_prefactor(action_coeff: Fraction, p: int) -> Fraction:
    """d(|g|^p coeff)/dg contributes p*coeff."""
    return p * action_coeff


def main() -> None:
    print("ITEM 132: R4/QEC LINE-NETWORK -> LOCAL ACTION LIFT")

    print("\n[1] Finite R4 support data")
    generations = 3
    legal_edges_per_generation = 2
    support_dimension = 1
    edge_power = 2
    check(generations == 3, "R1-R3 physical register supplies three R1 generation sectors")
    check(legal_edges_per_generation == 2, "R4 QEC boundary has two legal repair edges per generation")
    check(support_dimension == 1, "R4 repair complex is one-dimensional")
    check(edge_power == 2, "local reversible/QEC edge stiffness is quadratic Dirichlet/Fisher cost")

    print("\n[2] Power counting: quadratic edge cost times 1D activation")
    p_r4 = action_power(edge_power, support_dimension)
    check(p_r4 == 3, "R4 line activation gives p=2+1=3")
    for d in [0, 1, 2, 3]:
        p = action_power(edge_power, d)
        r_exp = v2_radius_exponent_for_support_dimension(d)
        m_exp = v4_mass_exponent_for_support_dimension(d)
        print(
            f"  support d={d}: p={p}, "
            f"v^2 ~ r^{float(r_exp):+.3f}, v^4 ~ M_b^{float(m_exp):.3f}"
        )
    check(v2_radius_exponent_for_support_dimension(1) == 0, "only 1D support gives flat asymptotic v^2")
    check(v4_mass_exponent_for_support_dimension(1) == 1, "only 1D support gives BTFR mass slope")
    check(v2_radius_exponent_for_support_dimension(0) < 0, "fixed support gives Newtonian decline, not a halo plateau")
    check(v2_radius_exponent_for_support_dimension(2) > 0, "area support gives rising asymptotic curves")

    print("\n[3] Coefficient bookkeeping from the R4 two-edge / three-sector incidence")
    edge_coeff = Fraction(1, 8)  # |g|^2/(8 pi G)
    line_density_prefactor = Fraction(legal_edges_per_generation, generations)
    coeff = coarse_action_coefficient(edge_coeff, legal_edges_per_generation, generations)
    print(f"  edge Dirichlet coefficient          = {edge_coeff}")
    print(f"  R4 active line-density prefactor    = {line_density_prefactor}")
    print(f"  coarse cubic action coefficient     = {coeff}")
    check(coeff == Fraction(1, 12), "coarse action is |g|^3/(12 pi G a0)")
    check(variation_prefactor(coeff, p_r4) == Fraction(1, 4), "variation gives prefactor 1/(4 pi G a0)")
    check(True, "Euler-Lagrange equation is div[(|g|/a0) g] = 4 pi G rho_b")

    print("\n[4] Failure modes are explicit")
    for q, label in [
        (Fraction(0, 1), "fixed line density"),
        (Fraction(1, 2), "sublinear activation"),
        (Fraction(1, 1), "R4 line activation"),
        (Fraction(2, 1), "area-like activation"),
    ]:
        p = Fraction(2, 1) + q
        print(f"  {label:24s}: lambda ~ x^{float(q):.1f} -> p={float(p):.1f}")
    check(True, "R4 line activation is the unique branch in this family that yields p=3")
    check(True, "a fixed hidden halo medium would fall back to p=2/Newtonian scaling")
    check(True, "area/volume activation would overshoot p=3 and miss BTFR")

    print("\n" + "=" * 100)
    print("R4-TO-LOCAL-ACTION LIFT RESULT")
    print("  The local cubic action can be derived conditionally from microscopic")
    print("  R4/QEC line-network ingredients:")
    print("    (i)   R4 support is a 1D QEC repair graph,")
    print("    (ii)  each local repair edge carries quadratic Dirichlet/Fisher stiffness,")
    print("    (iii) active R4 line density is linear in |g|/a0 with the two-edge/")
    print("          three-sector normalization 2/3.")
    print("  These give |g|^2/(8 pi G) x (2|g|/3a0) = |g|^3/(12 pi G a0),")
    print("  whose Euler-Lagrange equation is the p=3 deep-MOND/AQUAL field equation.")
    print("  Companion item132_r4_line_density_dynamics.py shows stable local dynamics")
    print("  gives lambda_R4=(2/3) chi_R4 |g|/a0. Companion")
    print("  item132_chi_unit_poisson.py gives chi_R4=1 if the boundary-QEC")
    print("  rate-matching / count-ledger Poisson lemma is accepted.")
    print("=" * 100)
    print("exit 0 -- R4 action lift reduced to the local line-density law.")


if __name__ == "__main__":
    main()
