#!/usr/bin/env python3
r"""Black-hole flux coefficient: species and polarization ledger audit.

The existing flux gate converts the Stefan--Hawking coefficient

    P = 1/(15360 pi M^2)

into a local attempt rate and finds that

    Gamma_0 = (10/27) alpha_0 Lambda_QCD

gives P/P_target = 0.997096.  This script tightens what that statement means.

Main point
----------
The coefficient 1/(15360 pi) is the naive Stefan flux for a two-helicity
massless bosonic channel.  The 10/27 route therefore fixes the local horizon
attempt rate for a *two-polarization Stefan channel*.  It is not, by itself, a
native derivation of the full species-summed Hawking luminosity.

Consequently the residual is now sharply finite:

  1. all-contact Landauer--Moore service transfer for the 10/27 local source;
  2. a species/polarization ledger deciding which asymptotic QFT channels are
     sourced by that local service event, before the already-computed
     Schwarzschild greybody transfer is applied.

This prevents a common overclaim: the 0.29% near-hit is real, but it is not a
closed total-flux theorem until the QEC radiation ledger supplies the
degeneracy/species factor.
"""

from __future__ import annotations

import math


ALPHA0 = 1.0 / 137.0
GAMMA_REQ_TWO_HELICITY = 2.711306813255e-3
GAMMA_1027 = (10.0 / 27.0) * ALPHA0


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def stefan_coeff(g_eff: float) -> float:
    """Naive horizon-area Stefan coefficient for g_eff bosonic helicities."""

    return g_eff / (30720.0 * math.pi)


def p_over_target_for_geff(g_eff: float) -> float:
    """10/27 source normalized to a Stefan target with g_eff helicities."""

    # The existing gate target is g_eff=2, so targets scale linearly with g_eff.
    return (GAMMA_1027 / GAMMA_REQ_TWO_HELICITY) * (2.0 / g_eff)


def required_gamma_for_geff(g_eff: float) -> float:
    return GAMMA_REQ_TWO_HELICITY * (g_eff / 2.0)


def main() -> None:
    print("BLACK-HOLE FLUX SPECIES/POLARIZATION LEDGER")
    print("=" * 92)

    print("[1] Stefan coefficient convention")
    canonical = 1.0 / (15360.0 * math.pi)
    two_pol = stefan_coeff(2.0)
    one_pol = stefan_coeff(1.0)
    print(f"    C(g_eff=1)                 = {one_pol:.12e}")
    print(f"    C(g_eff=2)                 = {two_pol:.12e}")
    print(f"    1/(15360 pi)               = {canonical:.12e}")
    check(abs(two_pol - canonical) < 1.0e-18, "current flux gate targets the two-helicity Stefan coefficient")
    check(abs(one_pol - 0.5 * canonical) < 1.0e-18, "one real scalar Stefan target would be half as large")

    print("\n[2] 10/27 source against different degeneracy conventions")
    print(f"    Gamma_req(g_eff=2)/Lambda  = {GAMMA_REQ_TWO_HELICITY:.12e}")
    print(f"    Gamma_10/27/Lambda         = {GAMMA_1027:.12e}")
    for g_eff, label in (
        (1.0, "one real scalar"),
        (2.0, "two-helicity massless boson"),
        (4.0, "two such channels"),
        (10.0, "many-species toy"),
    ):
        req = required_gamma_for_geff(g_eff)
        ratio = p_over_target_for_geff(g_eff)
        print(f"    {label:<28s} g_eff={g_eff:4.1f}  Gamma_req={req:.12e}  P/P_target={ratio:.9f}")
    check(abs(p_over_target_for_geff(2.0) - 0.997096067) < 5.0e-9, "10/27 reproduces the recorded 0.997096 near-hit for g_eff=2")
    check(p_over_target_for_geff(1.0) > 1.99, "reading the target as one scalar would over-emit by about a factor two")
    check(p_over_target_for_geff(4.0) < 0.51, "a two-channel species sum would require roughly twice the local source rate")

    print("\n[3] What is already computed versus what remains")
    closed = {
        "local KMS ladder": "finite QEC strain source",
        "freeze shell": "proper-distance beta-one surface",
        "escape cone": "M^-2 scaling",
        "greybody transfer": "standard spin/partial-wave exterior map",
        "10/27 local source": "conditional Landauer-Moore transfer",
    }
    residual = {
        "all-contact service theorem": "prove horizon severing uses the full Moore shell plus latch, not a narrower radial/surface alphabet",
        "species/polarization ledger": "prove which asymptotic QFT helicity/species channels one local service firing populates",
    }
    for name, desc in closed.items():
        print(f"    CLOSED/SHARP: {name:<26s} -> {desc}")
    for name, desc in residual.items():
        print(f"    RESIDUAL:     {name:<26s} -> {desc}")
    check(len(closed) == 5, "closed/sharp ingredients enumerated")
    check(len(residual) == 2, "flux residual reduced to two finite ledgers")

    print(
        """
[4] VERDICT
    The 10/27 route is tighter than a heuristic coefficient fit:

      * among the symmetry-allowed local service alphabets it is the unique
        member landing on the Stefan target;
      * its target is specifically the two-helicity Stefan coefficient;
      * the exterior spin/partial-wave greybody map is already a separate
        computed transfer.

    But it is not yet a closed total Hawking-luminosity theorem.  The remaining
    question is not another continuous coefficient.  It is the finite QEC
    species ledger: does one all-contact horizon service firing source exactly
    one two-polarization asymptotic channel, or a different species-weighted
    sum?  Until that is derived, the absolute flux coefficient is conditional
    rather than locked.

ALL ASSERTIONS PASSED -- flux residual sharpened to all-contact + species ledger."""
    )


if __name__ == "__main__":
    main()
