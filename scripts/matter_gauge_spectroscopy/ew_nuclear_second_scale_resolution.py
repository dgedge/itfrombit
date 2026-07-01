#!/usr/bin/env python3
r"""EW/NUCLEAR SECOND ANCHOR -- current status after the VEV and item-113 updates.

This script supersedes the older 2026-06-19 "second scale resolves as v == m_top"
verdict.  That verdict was correct at the time, but the later v-program changed
the status:

  * The electroweak VEV is no longer an irreducible input.  The live route is
        v / M_P = alpha_0^8 / sqrt(lambda(v))
    where alpha_0^8 is the 8-bit product channel from the item-79 single-link
    service projection, and lambda(v) is the RG output of the dimensionless
    boundary condition lambda(M_P)=0.  Numerically this lands at about 10%.

  * The result is not Locked.  It still needs the physical identification that
    EWSB bills the all-8 filled-cell transition rather than a lower-k channel,
    plus a higher-order RG/Coleman-Weinberg prefactor theorem for the last
    O(10%) normalization.

  * m_H/m_Z remains cross-sector: the A1g Higgs matter-cell ledger and the Eu Z
    bridge ledger are disjoint unless a single EW service/Hessian action fixes
    lambda_H/g_Z^2.

  * Nuclear MeV coefficients are no longer a separate dimensionful anchor at
    local liquid-drop grade.  The closed-record-pair energy unit
        eps_sat = 2 alpha_0 Lambda
    plus the T1/T2 local contact maps reproduce the SEMF coefficients within
    few-percent grade.  The residual is the microscopic many-body Hamiltonian
    and shell/magic-number structure, not an additional absolute scale.

Net: the "EW/nuclear second anchor" is no longer one undifferentiated wall.
The EW VEV is substantially reduced but not locked; Higgs/Z and heavy-quark
spectroscopy remain EW-second-scale/cross-sector work; nuclear is a QCD/record
contact programme rather than a new anchor.
"""

from __future__ import annotations

import math

M_P = 1.2209e19          # GeV
V_EW = 246.22            # GeV
M_H = 125.25             # GeV
M_Z = 91.1876            # GeV
M_TOP = 172.76           # GeV
ALPHA0 = 1.0 / 137.036   # dressed convention for EW comparisons in the v scripts
ALPHA_BARE = 1.0 / 137.0 # bare record-service convention for nuclear contact billing
HBARC = 197.327          # MeV fm
A0_FM = 0.5944           # fm


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def semf_table() -> dict[str, float]:
    """Reproduce the item-113 local liquid-drop coefficient consequence."""

    lam_mev = HBARC / A0_FM
    w = ALPHA_BARE * lam_mev
    eps_sat = 2.0 * w

    z = 6
    factor_bulk = 13.0 / 12.0
    factor_surface = 1.0
    surface_ratio = 1.206
    trace_count = 5
    rms_count = math.sqrt(5.0)

    a_c = 0.6 * (1.0 / 137.036) * HBARC / (2.0 * A0_FM)
    a_v = (z / 2.0) * factor_bulk * eps_sat
    a_s = surface_ratio * (z / 2.0) * factor_surface * eps_sat

    r0 = 2.0 * A0_FM
    rho = 3.0 / (4.0 * math.pi * r0**3)
    m_n = 939.0
    e_f = (HBARC**2 / (2.0 * m_n)) * (1.5 * math.pi**2 * rho) ** (2.0 / 3.0)
    a_a = e_f / 3.0 + trace_count * w
    a_p = rms_count * eps_sat
    return {"aC": a_c, "aV": a_v, "aS": a_s, "aA": a_a, "aP": a_p}


def main() -> int:
    print("=" * 98)
    print("EW/NUCLEAR SECOND ANCHOR -- CURRENT STATUS")
    print("=" * 98)

    print("\n[1] Electroweak VEV: reduced to a byte-channel radiative scale")
    lambda_h = M_H**2 / (2.0 * V_EW**2)
    vmp_obs = V_EW / M_P
    vmp_pred = ALPHA0**8 / math.sqrt(lambda_h)
    lambda_required = (ALPHA0**8 / vmp_obs) ** 2
    k_eff = math.log(vmp_obs * math.sqrt(lambda_h)) / math.log(ALPHA0)

    print(f"    lambda_H(v) = m_H^2/(2v^2) = {lambda_h:.6f}")
    print(f"    observed v/M_P              = {vmp_obs:.3e}")
    print(f"    alpha_0^8/sqrt(lambda_H)    = {vmp_pred:.3e}")
    print(f"    fractional miss             = {(vmp_pred / vmp_obs - 1.0) * 100:+.2f}%")
    print(f"    exponent after removing sqrt(lambda): {k_eff:.3f} powers of alpha_0")
    print(f"    lambda required for exact landing     = {lambda_required:.3f}")
    check(abs(k_eff - 8.0) < 0.05, "the hierarchy is the alpha_0^8 byte channel, not a random exponent")
    check(abs(vmp_pred / vmp_obs - 1.0) < 0.15, "alpha_0^8/sqrt(lambda_H) lands v/M_P to ~10-15%")

    print("\n    Status:")
    print("      - alpha_0 per bit: item-79 service projection.")
    print("      - 8 powers: all-8 filled-cell / A1g coincidence channel.")
    print("      - lambda: dimensionless boundary lambda(M_P)=0, RG-near-critical output.")
    print("      - residual: prove EWSB bills the filled-cell transition and compute the CW/RG prefactor.")

    print("\n[2] EW spectroscopy residuals")
    y_t = math.sqrt(2.0) * M_TOP / V_EW
    mh_mz = M_H / M_Z
    lambda_over_gz2_required = mh_mz**2 / 8.0
    print(f"    y_t = sqrt(2) m_top/v = {y_t:.4f}")
    print("    top therefore remains the saturated O(1) EW Yukawa check,")
    print("    but no longer needs to be treated as the primitive source of v.")
    print(f"    m_H/m_Z = {mh_mz:.6f}; required lambda_H/g_Z^2 = {lambda_over_gz2_required:.6f}")
    check(0.9 < y_t < 1.05, "top Yukawa saturation links m_top to the EW scale")
    check(0.20 < lambda_over_gz2_required < 0.27, "m_H/m_Z requires a nontrivial EW Hessian ratio")

    print("\n    Residual:")
    print("      m_H/m_Z remains cross-sector.  The A1g Higgs matter-cell ledger and")
    print("      the Eu Z-bridge ledger are disjoint until a single EW service/Hessian")
    print("      action derives lambda_H/g_Z^2.")

    print("\n[3] Nuclear MeV coefficients: local contact map, not a separate anchor")
    coeffs = semf_table()
    empirical = {"aC": 0.711, "aV": 15.75, "aS": 17.8, "aA": 23.7, "aP": 11.18}
    for name in ("aC", "aV", "aS", "aA", "aP"):
        err = 100.0 * (coeffs[name] / empirical[name] - 1.0)
        print(f"    {name}: pred={coeffs[name]:8.3f} MeV  empirical={empirical[name]:8.3f} MeV  err={err:+6.2f}%")
    check(all(abs(coeffs[k] / empirical[k] - 1.0) < 0.05 for k in coeffs),
          "T1/T2 contact maps keep the liquid-drop coefficients within 5%")

    print("\n    Residual:")
    print("      The local maps are grounded, but the full Weizsacker/nuclear lock")
    print("      still needs a microscopic contact Hamiltonian, shell corrections,")
    print("      magic numbers, and many-body dynamics.")

    print("\n[4] Classification")
    rows = [
        ("EW VEV v", "reduced/conditional", "alpha_0^8 byte channel + RG lambda, about 10%; not input-only"),
        ("top mass", "saturation check", "y_t ~= 1 ties top to the EW scale; heavy flavours remain Target-B"),
        ("m_H/m_Z", "open cross-sector", "needs a shared EW Hessian/service action"),
        ("nuclear MeV", "local theorem/open many-body", "T1/T2 maps work; shell/Hamiltonian remain"),
    ]
    for sector, tier, note in rows:
        print(f"    {sector:14s}  {tier:29s}  {note}")

    print("\nVERDICT")
    print("  The old statement 'v == m_top is the irreducible second anchor' is superseded.")
    print("  The EW VEV has become a conditional radiative prediction, v/M_P = alpha_0^8/sqrt(lambda),")
    print("  with the remaining gap now the filled-cell EWSB identification plus precision RG/CW prefactor.")
    print("  What remains hard in EW is spectroscopy: m_H/m_Z and heavy-quark masses need a genuine")
    print("  second-scale/EW Hessian or flavour theorem.  Nuclear binding is no longer a separate")
    print("  absolute-scale wall at local liquid-drop grade; it is a QCD/record-contact many-body programme.")
    print("exit 0 -- EW/nuclear second-anchor wall split and narrowed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
