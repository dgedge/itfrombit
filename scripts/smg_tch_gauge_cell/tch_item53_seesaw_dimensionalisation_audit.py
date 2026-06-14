#!/usr/bin/env python3
"""
Audit ANCHOR §15 item 53: seesaw-product dimensionalisation.

What item 53 has in hand:
  * Part 18 gives an exact dimensionless identity
      |E1| |E2| = (1 + sqrt(2)) (sqrt(2) - 1) = 1.

What item 53 still needs:
  * A substrate rule that injects a physical mass scale, instead of assigning
    the Higgs VEV v = 246 GeV by hand and/or inverting the measured neutrino
    mass to define M_R.

This script separates:
  A. exact dimensionless Part 18 algebra;
  B. standard seesaw inversion (useful, but not a derivation);
  C. independent-scale checks using GUT/proton-scale inputs;
  D. the only tempting internal numerical foothold found here:
       mu_nu / mu_l ~= exp(-24),
     where 24 can be read as 6 Feshbach P-Q couplings times [8,4,4] distance 4.

The exp(-24) line is explicitly labelled a tempting ansatz, not a theorem.
The companion script tch_item53_tunnel_action_refutation.py refutes it as a
consequence of the actual Part 18 H_PQ matrix: the six entries are parallel
channels, not serial barriers.
"""

from __future__ import annotations

import math


# Physical constants / repo conventions.
V_HIGGS_GEV = 246.0
M_TAU_GEV = 1.77686
DM31_GEV2 = 2.511e-3 * 1e-18  # eV^2 -> GeV^2
M3_PART05_GEV = 50.116103e-3 * 1e-9  # meV -> eV -> GeV
ALPHA_INV = 137.035999084
ALPHA = 1.0 / ALPHA_INV
MPL_GEV = 1.2209e19

R_LEPTON = math.sqrt(2.0)
DELTA_LEPTON = 2.0 / 9.0
R_NU = 1.0
DELTA_NU = 1.0 / 3.0


def koide_factors(R: float, delta: float) -> list[float]:
    """m_n / mu = (1 + R cos(delta + 2 pi n/3))^2."""

    return [(1.0 + R * math.cos(delta + 2.0 * math.pi * n / 3.0)) ** 2 for n in range(3)]


def charged_lepton_mu_from_tau() -> float:
    """Use the charged-lepton n=0 heavy/tau node to fix mu_l."""

    f_tau = koide_factors(R_LEPTON, DELTA_LEPTON)[0]
    return M_TAU_GEV / f_tau


def neutrino_mu_from_dm31() -> float:
    """Use Part 05 local NH assignment (raw n=1,n=2,n=0) and Delta m31^2."""

    f = koide_factors(R_NU, DELTA_NU)
    ordered = [f[1], f[2], f[0]]
    denom = ordered[2] ** 2 - ordered[0] ** 2
    return math.sqrt(DM31_GEV2 / denom)


def neutrino_masses_from_mu(mu_gev: float) -> list[float]:
    """Part 05 local NH masses in GeV."""

    f = koide_factors(R_NU, DELTA_NU)
    return [mu_gev * f[i] for i in (1, 2, 0)]


def valid_r123(state: tuple[int, ...]) -> bool:
    """Part 18 bit order: (G0,G1,chi,I3,LQ,C0,C1,W)."""

    G0, G1, chi, _I3, LQ, C0, C1, W = state
    if G0 == 1 and G1 == 1:
        return False
    if W != chi:
        return False
    if (LQ == 0) != ((C0, C1) == (0, 0)):
        return False
    return True


def violates_r4_nu_r(state: tuple[int, ...]) -> bool:
    """R4-forbidden right-handed neutrino pseudocodewords."""

    _G0, _G1, chi, I3, LQ, C0, C1, _W = state
    return LQ == 0 and I3 == 0 and chi == 1 and (C0, C1) == (0, 0)


def hamming(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x != y for x, y in zip(a, b))


def code_partition() -> tuple[list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Return P and Q in the Part 18 R1+R2+R3-valid 48-space."""

    r123 = [
        tuple((n >> (7 - i)) & 1 for i in range(8))
        for n in range(256)
        if valid_r123(tuple((n >> (7 - i)) & 1 for i in range(8)))
    ]
    q = [s for s in r123 if violates_r4_nu_r(s)]
    p = [s for s in r123 if not violates_r4_nu_r(s)]
    return p, q


def fmt_gev_as_mev(x: float) -> str:
    return "%.6f meV" % (x * 1e12)


def main() -> None:
    print("=" * 78)
    print("TCH item-53 seesaw-dimensionalisation audit")
    print("=" * 78)

    e1 = 1.0 + math.sqrt(2.0)
    e2 = math.sqrt(2.0) - 1.0
    print("A. Part 18 dimensionless algebra")
    print("   |E1| = 1 + sqrt(2) = %.12f" % e1)
    print("   |E2| = sqrt(2) - 1 = %.12f" % e2)
    print("   |E1| |E2| = %.16f" % (e1 * e2))
    print("   status: exact, but dimensionless.")
    print()

    mu_l = charged_lepton_mu_from_tau()
    mu_nu = neutrino_mu_from_dm31()
    masses = neutrino_masses_from_mu(mu_nu)
    ratio = mu_nu / mu_l
    s_required = -math.log(ratio)
    print("B. Scales implied by Part 05")
    print("   mu_l  from charged-lepton Koide/tau = %.12e GeV" % mu_l)
    print("   mu_nu from Delta m31^2             = %.12e GeV" % mu_nu)
    print("   mu_nu / mu_l                       = %.12e" % ratio)
    print("   -ln(mu_nu / mu_l)                  = %.6f" % s_required)
    print("   Part 05 masses                     = [%s, %s, %s]" % tuple(fmt_gev_as_mev(m) for m in masses))
    print()

    mr_from_m3 = V_HIGGS_GEV ** 2 / masses[2]
    mr_from_mu = V_HIGGS_GEV ** 2 / mu_nu
    print("C. Standard seesaw dimensional injection")
    print("   M_R inferred from m3  = v^2 / m3     = %.6e GeV" % mr_from_m3)
    print("   M_R inferred from mu  = v^2 / mu_nu  = %.6e GeV" % mr_from_mu)
    print("   status: useful bookkeeping, but this is an inversion unless M_R or v")
    print("   is derived independently from the substrate.")
    print()

    print("D. Independent heavy-scale checks")
    for label, mr in [
        ("Part-19 rough E_R3 = 1e15 GeV", 1e15),
        ("ANCHOR proton-scale M_X = 1e16 GeV", 1e16),
        ("needed from Part-05 m3", mr_from_m3),
    ]:
        m3 = V_HIGGS_GEV ** 2 / mr
        mu = m3 / koide_factors(R_NU, DELTA_NU)[0]
        print("   %-34s -> m3 = %s, mu_nu = %s" % (label, fmt_gev_as_mev(m3), fmt_gev_as_mev(mu)))
    print("   status: 1e15 GeV is close at order-one level; 1e16 GeV is not.")
    print("   The Part-19 convergence is therefore suggestive, not a derivation.")
    print()

    print("E. Internal exponential-suppression candidates")
    candidates = [
        ("6 H_PQ entries x [8,4,4] distance 4", 6 * 4),
        ("3 nu_R pseudocodewords x 8-bit register", 3 * 8),
        ("d_eff=3 x 9-plaquette sites", 3 * 9),
        ("2 x 9-plaquette sites", 2 * 9),
        ("single [8,4,4] distance", 4),
    ]
    for label, action in candidates:
        mu_pred = mu_l * math.exp(-action)
        pred_masses = neutrino_masses_from_mu(mu_pred)
        rel_mu = (mu_pred / mu_nu) - 1.0
        print("   S = %-2d  %-42s mu_nu = %s  rel = %+7.2f%%  m3 = %s" % (
            action,
            label,
            fmt_gev_as_mev(mu_pred),
            100.0 * rel_mu,
            fmt_gev_as_mev(pred_masses[2]),
        ))
    print("   exact required action S_req = %.6f, so S=24 misses by %.6f in exponent." % (
        s_required,
        24.0 - s_required,
    ))
    print()

    p_space, q_space = code_partition()
    nearest = []
    one_bit_pairs = 0
    for q in q_space:
        distances = sorted((hamming(q, p), p) for p in p_space)
        nearest.append(distances[0][0])
        one_bit_pairs += sum(1 for p in p_space if hamming(q, p) == 1)
    print("F. Literal code-partition distance audit")
    print("   R1+R2+R3-valid states: P = %d, Q(nu_R) = %d" % (len(p_space), len(q_space)))
    print("   nearest raw Hamming distances Q -> P = %s" % nearest)
    print("   raw one-bit Q-P neighbour pairs = %d" % one_bit_pairs)
    print("   status: the literal P/Q boundary is not distance 4. The S=6x4")
    print("   candidate must use [8,4,4] stabilizer/code distance as an action")
    print("   penalty, not raw Hamming distance to the nu_R pseudocodeword.")
    print()

    print("G. Planck/alpha collision check")
    mr_alpha2 = 2.0 * (ALPHA ** 2) * MPL_GEV
    m3_alpha2 = V_HIGGS_GEV ** 2 / mr_alpha2
    print("   2 alpha^2 M_Pl = %.6e GeV -> m3 = %s" % (mr_alpha2, fmt_gev_as_mev(m3_alpha2)))
    print("   status: another simple framework-looking expression lands nearby.")
    print("   That weakens any claim that exp(-24) alone is evidentially decisive.")
    print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  Closed: Part 18 supplies an exact dimensionless seesaw product.")
    print("  Not closed: item 53's dimensionalisation is still applied by hand if it")
    print("  simply injects v = 246 GeV or defines M_R by inverting m_nu.")
    print()
    print("  Numerical foothold: the required scale suppression is very close to")
    print("  exp(-24), and 24 has framework-looking decompositions such as 6x4.")
    print("  This predicts mu_nu at the 10-20% level and gives the right scale.")
    print("  However, the companion refutation shows that 6x4 does not follow from")
    print("  the actual Part 18 H_PQ Feshbach operator.")
    print()
    print("  Obstruction: Part 18's explicit H_PQ entries are one-bit couplings, so")
    print("  the code-distance-4 action is not yet derived; the raw P/Q boundary is")
    print("  Hamming-distance 1. Proving a stabilizer-action penalty, rather than a")
    print("  raw path-length penalty, is the real item-53 target. Without that,")
    print("  exp(-24) is a promising ansatz.")


if __name__ == "__main__":
    main()
