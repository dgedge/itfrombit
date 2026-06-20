#!/usr/bin/env python3
r"""Black-hole freeze-surface and greybody transfer gate.

The finite horizon channel now gives a local source ladder

    P(F) proportional to g_Q(F) exp[-beta_eff F] .

This script attacks the next observable-map question: how a local QCD-scale
service event can be seen at infinity as a Hawking-scale quantum.

Near a Schwarzschild horizon, a source at proper distance rho from the horizon
has

    sqrt(f) = rho / (2 r_s) + O((rho/r_s)^3),
    T_H    = 1 / (4 pi r_s).

If the local service energy is

    E_loc(F) = eps_F F Lambda_QCD,

then the energy at infinity is E_inf = sqrt(f) E_loc and therefore

    E_inf / T_H = 2 pi rho eps_F Lambda_QCD F.

Consequences:

  * a universal proper-distance freeze shell gives an M-independent Hawking
    line temperature;
  * matching a dimensionless KMS beta_eff fixes only the O(1) shell position
    rho Lambda_QCD = beta_eff / (2 pi eps_F);
  * adding the Schwarzschild escape cone supplies the extra near-horizon f
    factor needed for the total luminosity to scale as M^-2;
  * exterior greybody factors and absolute attempt rate remain separate maps.

Exit 0 means the freeze-surface relation is self-consistent across black-hole
masses and the remaining open pieces are correctly localized.
"""

from __future__ import annotations

from collections import Counter
import math


PHI = (math.sqrt(5.0) - 1.0) / 2.0
EPS_F = 1.0 / (2.0 * PHI)
LAMBDA_QCD_GEV = 0.332
M_PLANCK_GEV = 1.220890e19
GEV_PER_KG = 5.60959e26
M_SUN_GEV = 1.98892e30 * GEV_PER_KG
GEV_INV_TO_FM = 0.1973269804
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def horizon_radius(m_gev: float) -> float:
    return 2.0 * m_gev / M_PLANCK_GEV**2


def hawking_temperature(m_gev: float) -> float:
    return 1.0 / (4.0 * math.pi * horizon_radius(m_gev))


def sqrt_f_near_horizon(rho_gev_inv: float, m_gev: float) -> float:
    return rho_gev_inv / (2.0 * horizon_radius(m_gev))


def coordinate_offset_fraction(rho_gev_inv: float, m_gev: float) -> float:
    """Return (r-r_s)/r_s from the near-horizon proper distance relation.

    For Schwarzschild, rho ~= 2 sqrt(r_s (r-r_s)), so
    (r-r_s)/r_s ~= (rho/(2 r_s))^2.
    """

    return sqrt_f_near_horizon(rho_gev_inv, m_gev) ** 2


def horizon_cells(m_gev: float) -> float:
    rs = horizon_radius(m_gev)
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2


def escape_fraction_outward(rho_gev_inv: float, m_gev: float) -> float:
    """Near-horizon photon escape-cone fraction in the outward hemisphere.

    For a Schwarzschild static observer, sin(alpha_c)=b_c sqrt(f)/r with
    b_c=3 sqrt(3) M = (3 sqrt(3)/2) r_s.  The outward-hemisphere fraction is
    1-cos(alpha_c), which is (27/8) f near the horizon.
    """

    rs = horizon_radius(m_gev)
    sf = sqrt_f_near_horizon(rho_gev_inv, m_gev)
    sin_alpha = min(1.0, (3.0 * math.sqrt(3.0) / 2.0) * sf)
    if sin_alpha < 1.0e-4:
        return 0.5 * sin_alpha * sin_alpha
    cos_alpha = math.sqrt(max(0.0, 1.0 - sin_alpha * sin_alpha))
    return 1.0 - cos_alpha


def beta_from_shell(rho_lambda: float) -> float:
    return 2.0 * math.pi * rho_lambda * EPS_F


def rho_lambda_for_beta(beta_eff: float) -> float:
    return beta_eff / (2.0 * math.pi * EPS_F)


def source_ladder(beta_eff: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta_eff * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return dict(sorted((f, v / z) for f, v in raw.items()))


def mean_f(src: dict[int, float]) -> float:
    return sum(f * p for f, p in src.items())


def coarse_greybody_transfer(src: dict[int, float], eta: dict[int, float]) -> dict[str, float]:
    """A normalized example transfer map with an exterior escape factor.

    The eta_F factors stand for the spin/partial-wave greybody barrier.  They
    are not derived here; the point is to keep them explicit rather than hiding
    them in the local KMS ladder.
    """

    escaped = {f: src[f] * eta.get(f, 0.0) for f in src}
    absorbed = sum(src[f] * (1.0 - eta.get(f, 0.0)) for f in src)
    return {
        "low(F<=4)": sum(v for f, v in escaped.items() if f <= 4),
        "mid(5<=F<=7)": sum(v for f, v in escaped.items() if 5 <= f <= 7),
        "high(F>=8)": sum(v for f, v in escaped.items() if f >= 8),
        "barrier_absorbed": absorbed,
    }


def summarize_shell(beta_eff: float) -> tuple[float, float]:
    rho_lambda = rho_lambda_for_beta(beta_eff)
    rho_gev_inv = rho_lambda / LAMBDA_QCD_GEV
    return rho_lambda, rho_gev_inv * GEV_INV_TO_FM


def main() -> None:
    print("BLACK-HOLE FREEZE-SURFACE / GREYBODY TRANSFER GATE")
    print("=" * 92)

    print("[1] Near-horizon freeze-surface identity")
    beta_target = 1.0
    rho_lambda, rho_fm = summarize_shell(beta_target)
    print(f"    eps_F = 1/(2 phi) = {EPS_F:.12f}")
    print(f"    beta_eff target = {beta_target:.6f}")
    print(f"    required rho* Lambda_QCD = {rho_lambda:.12f} = phi/pi")
    print(f"    required rho* = {rho_fm:.6f} fm")
    check(abs(rho_lambda - PHI / math.pi) < 1.0e-15, "beta=1 shell is rho Lambda = phi/pi")

    a0_lambda = 1.0
    beta_a0 = beta_from_shell(a0_lambda)
    print(f"    if rho = a0 = 1/Lambda_QCD, beta_eff = {beta_a0:.6f} = pi/phi")
    check(abs(beta_a0 - math.pi / PHI) < 1.0e-15, "a0 shell gives beta_eff=pi/phi")

    print("\n[2] Mass-independence check")
    rho_gev_inv = rho_lambda / LAMBDA_QCD_GEV
    for m_solar in (1.0e-18, 3.0, 30.0, 4.3e6):
        m = m_solar * M_SUN_GEV
        rs = horizon_radius(m)
        th = hawking_temperature(m)
        sf = sqrt_f_near_horizon(rho_gev_inv, m)
        eps_coord = coordinate_offset_fraction(rho_gev_inv, m)
        ratios = []
        for f in (3, 4, 6, 9, 12):
            e_loc = EPS_F * f * LAMBDA_QCD_GEV
            e_inf = sf * e_loc
            ratios.append(e_inf / th / f)
        spread = max(ratios) - min(ratios)
        print(
            f"    M={m_solar:g} Msun: r_s={rs:.3e} GeV^-1, "
            f"sqrt(f)={sf:.3e}, (r-r_s)/r_s={eps_coord:.3e}, "
            f"E_inf/(T_H F)={ratios[0]:.12f}"
        )
        check(spread < 1.0e-12, "same shell gives F-independent beta_eff")
        check(abs(ratios[0] - beta_target) < 1.0e-12, "same shell gives M-independent beta_eff")

    print("\n[3] Source ladder and explicit greybody placeholder")
    src = source_ladder(beta_target)
    print(f"    source P(F), beta=1: {src}")
    transparent = {f: 1.0 for f in src}
    lowpass = {f: max(0.0, min(1.0, (10.0 - f) / 7.0)) for f in src}
    obs_trans = coarse_greybody_transfer(src, transparent)
    obs_lowpass = coarse_greybody_transfer(src, lowpass)
    print(f"    transparent transfer: {obs_trans}")
    print(f"    sample lowpass barrier: {obs_lowpass}")
    check(obs_trans != obs_lowpass, "local ladder and exterior barrier are distinct observable maps")

    print("\n[4] Flux scaling: escape cone supplies the Hawking M^-2 factor")
    e_loc_mean = EPS_F * mean_f(src) * LAMBDA_QCD_GEV
    rows = []
    for m_solar in (3.0, 30.0, 4.3e6):
        m = m_solar * M_SUN_GEV
        sf = sqrt_f_near_horizon(rho_gev_inv, m)
        f_schw = sf * sf
        n_h = horizon_cells(m)
        no_cone = n_h * e_loc_mean * f_schw
        cone = no_cone * escape_fraction_outward(rho_gev_inv, m)
        rows.append((m_solar, no_cone, cone, cone * m * m))
        print(
            f"    M={m_solar:g} Msun: no-cone power unit={no_cone:.6e}, "
            f"with escape cone={cone:.6e}, cone*M^2={cone*m*m:.6e}"
        )
    no_cone_spread = max(r[1] for r in rows) / min(r[1] for r in rows) - 1.0
    cone_m2_spread = max(r[3] for r in rows) / min(r[3] for r in rows) - 1.0
    check(no_cone_spread < 1.0e-12, "redshift-only horizon shell gives mass-independent power")
    check(cone_m2_spread < 1.0e-12, "escape-cone dressed flux scales as M^-2")

    print("\n[5] What is now closed, and what is not")
    print(
        """
    CLOSED BY THIS SCRIPT:
      A fixed proper-distance freeze shell converts local Lambda_QCD-scale
      service jumps to Hawking-temperature energies at infinity:

          E_inf/T_H = (2 pi rho Lambda_QCD eps_F) F.

      The required shell for beta_eff=1 is rho Lambda_QCD = phi/pi.  Choosing
      rho=a0 instead gives beta_eff=pi/phi; that is an O(1) scheduler/shell
      convention, not an M-dependent tuning.

      The Schwarzschild escape cone supplies the missing near-horizon f factor:
      a redshift-only shell gives mass-independent power, while the escape-cone
      dressed source has the Hawking M^-2 luminosity scaling.

    STILL OPEN:
      The QEC geometry must select the actual O(1) proper shell position
      rho_*; the exterior spin/partial-wave greybody barrier must be computed;
      and an absolute attempt-rate Gamma_0 is needed for the flux coefficient.
      The local KMS ladder alone cannot supply those upper maps.

ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
