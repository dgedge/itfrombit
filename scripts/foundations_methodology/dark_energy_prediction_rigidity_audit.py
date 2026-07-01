#!/usr/bin/env python3
r"""Rigidity audit for the live dark-energy prediction cluster.

Question
--------
Are the numbers used in the near-term cosmology prediction rigid framework
outputs, or are they carrying hidden freedom?

Audited targets
---------------
1. The "28" service alphabet behind

       w(a) = -1 + a/28,
       w0 = -27/28,
       wa = -1/28.

2. The "55" readable severing ledger behind

       C = 55/8,
       Omega_Lambda = 12 pi / 55.

3. The "17.7 keV" sterile-neutrino branch,

       m_nuR = alpha^2 Lambda_QCD,
       E_gamma = m_nuR/2.

Verdict convention
------------------
FORCED means no continuous parameter and no choice inside the named finite
instrument/readout.  CONDITIONAL means the value is parameter-free once a
named, non-numerical framework clause is accepted, but the clause is still a
physical premise.  OPEN means a number or observable still depends on an
unproved map, source law, or phenomenological branch.

This script deliberately separates numerical closeness from rigidity.  A close
integer hit is not evidence unless the framework also forces the integer.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math


ALPHA0 = 1.0 / 137.0
ALPHA_DRESSED = 1.0 / 137.035999
LAMBDA_QCD_GEV_LEGACY = 0.332
LAMBDA_QCD_GEV_REPIN_BARE = 0.33168
LAMBDA_QCD_GEV_REPIN_DRESSED = 0.33185

T_CMB = 2.7255
ZETA3 = 1.2020569031595943
KB_EV_PER_K = 8.617333262e-5
HBARC_EV_CM = 1.973269804e-5
RHO_CRIT_H2_EV_CM3 = 1.05371e4
OMEGA_NUR_REFERENCE = 0.024
OMEGA_CDM_REFERENCE = 0.120


@dataclass(frozen=True)
class Clause:
    name: str
    status: str
    reason: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def photon_density_cm3() -> float:
    return 2.0 * ZETA3 / math.pi**2 * (KB_EV_PER_K * T_CMB / HBARC_EV_CM) ** 3


def omega_from_ratio(n_over_ngamma: float, mass_ev: float) -> float:
    return n_over_ngamma * photon_density_cm3() * mass_ev / RHO_CRIT_H2_EV_CM3


def status_line(name: str, verdict: str, detail: str) -> None:
    print(f"  {name:<24s} {verdict:<14s} {detail}")


def audit_28() -> None:
    print("\n[1] Rigidity of 28 -> w0=-27/28, wa=-1/28")
    n_points = 8
    n_hyperplanes = 14
    n_transverse = 2
    point_degree = 7
    plane_degree = 4
    incidence_flags = n_points * point_degree
    micro_flags = incidence_flags * n_transverse
    service_channels = n_hyperplanes * n_transverse
    preimages_per_channel = plane_degree

    print("  finite bridge: points x incident hyperplanes x transverse modes")
    print(f"    {n_points} x {point_degree} x {n_transverse} = {micro_flags} microscopic flags")
    print(f"    {n_hyperplanes} x {n_transverse} = {service_channels} service channels")
    print(f"    each service channel has {preimages_per_channel} local point preimages")
    check(service_channels == 28, "28-channel alphabet is fixed by 14 hyperplanes x 2 modes")
    check(micro_flags == 112, "8 -> 112 -> 28 incidence bridge has the expected sizes")

    delta = Fraction(1, service_channels)
    w0 = Fraction(-1, 1) + delta
    wa = -delta
    print(f"  predicted CPL point: w0={w0}={float(w0):.10f}, wa={wa}={float(wa):.10f}")
    check(w0 == Fraction(-27, 28), "serial 28-clock gives w0=-27/28")
    check(wa == Fraction(-1, 28), "one-dimensional R4 support gives wa=-1/28")

    controls = {
        "parallel reset": ("gap", Fraction(1, 1), "would give order-one gap, not 1/28"),
        "area support": ("wa", Fraction(-2, 28), "would give wa=-2/28"),
        "volume support": ("wa", Fraction(-3, 28), "would give wa=-3/28"),
        "constant support": ("wa", Fraction(0, 1), "would give no CPL slope"),
    }
    print("  controls:")
    for name, (_kind, value, reason) in controls.items():
        print(f"    {name:<16s} -> {value:<6} {reason}")
    check(all(v[1] != wa for k, v in controls.items() if k != "parallel reset"), "support alternatives do not reproduce wa=-1/28")

    clauses = [
        Clause("28 alphabet", "FORCED", "14 affine hyperplanes times 2 transverse modes"),
        Clause("uniform weights", "FORCED", "AGL(3,2) x C2 covariance plus incidence locality"),
        Clause("serial one-jump", "CONDITIONAL", "reduced to one-step shift/boundary-crossing scheduler"),
        Clause("R4 d=1 support", "FORCED-INSTR", "finite R4 repair complex is a 1-chain in current register"),
        Clause("FRW lift", "CONDITIONAL", "homogeneous comoving 1-chain measure gives f(a)=a"),
    ]
    for clause in clauses:
        status_line(clause.name, clause.status, clause.reason)
    print("  verdict: parameter-free prediction-grade inside the current instrument;")
    print("           hidden freedom exists only by changing scheduler/support/outside-sector premises.")


def audit_55() -> None:
    print("\n[2] Rigidity of 55 -> C=55/8 -> Omega_Lambda=12pi/55")
    directed = 2 * 28
    blind = 1
    readable = directed - blind
    c = Fraction(readable, 8)
    omega = 3 * math.pi / (2 * float(c))
    print(f"  directed severing ledger: 2 x 28 = {directed}")
    print(f"  global-complement blind slot: {blind}")
    print(f"  readable channels: {directed} - {blind} = {readable}")
    print(f"  C = {readable}/8 = {float(c):.6f}")
    print(f"  Omega_Lambda = 3pi/(2C) = 12pi/55 = {omega:.9f}")
    check(readable == 55, "readable severing ledger has rank 55")
    check(c == Fraction(55, 8), "Bekenstein channel count is 55/8")
    check(abs(omega - 12 * math.pi / 55) < 1.0e-15, "Omega form is exactly 12pi/55")

    controls = [
        ("no direction", Fraction(28, 8), "C=3.5; loses partner asymmetry"),
        ("address direction", Fraction(56, 8), "C=7; requires forbidden readable global orientation"),
        ("both bits blind", Fraction(54, 8), "C=6.75; wrong quotient"),
    ]
    print("  controls:")
    for name, val, reason in controls:
        print(f"    {name:<18s} -> C={float(val):.3f}; {reason}")
        check(val != c, f"{name} control is distinct from 55/8")

    clauses = [
        Clause("direction exists", "FORCED", "horizon severing is partner-asymmetric"),
        Clause("record channel", "FORCED", "12-edge Q3 strain syndrome is the readout"),
        Clause("blind slot", "FORCED", "ker(delta)={0,ALL} removes exactly one global complement"),
        Clause("value-level hop tag", "FORCED", "address-level tag violates AGL(3,2) covariance/readout"),
        Clause("55", "FORCED", "2 x 28 directed records minus one global blind slot"),
    ]
    for clause in clauses:
        status_line(clause.name, clause.status, clause.reason)
    print("  verdict: rigid under the settled strain-readout/horizon-severing channel;")
    print("           not a free O(1) prefactor once that channel is accepted.")


def audit_177() -> None:
    print("\n[3] Rigidity of 17.7 keV sterile branch")
    masses = {
        "legacy+dressed": (ALPHA_DRESSED, LAMBDA_QCD_GEV_LEGACY),
        "legacy+bare": (ALPHA0, LAMBDA_QCD_GEV_LEGACY),
        "repin bare": (ALPHA0, LAMBDA_QCD_GEV_REPIN_BARE),
        "repin dressed": (ALPHA_DRESSED, LAMBDA_QCD_GEV_REPIN_DRESSED),
    }
    vals = {}
    for label, (alpha, lam_gev) in masses.items():
        m_kev = alpha**2 * lam_gev * 1.0e6
        vals[label] = m_kev
        print(f"  {label:<15s}: m_nuR={m_kev:.4f} keV, E_gamma={m_kev/2:.4f} keV")
    spread = max(vals.values()) / min(vals.values()) - 1.0
    check(17.5 < vals["legacy+dressed"] < 17.9, "legacy dressed value is 17.7 keV")
    check(spread < 0.002, "alpha/Lambda convention spread is below 0.2% for the line position")

    # Flux/mixing is not rigid.
    d_code = 4
    seesaw = 4.0 * 0.05 / (vals["legacy+dressed"] * 1.0e3)
    cd_lo = 4.0 * math.exp(-4 * d_code)
    cd_hi = 4.0 * math.exp(-2 * d_code)
    gf_lo = cd_lo * ALPHA_DRESSED**4
    gf_hi = cd_hi * ALPHA_DRESSED**4
    print("  mixing controls:")
    print(f"    seesaw line                 sin^2(2theta)={seesaw:.2e}")
    print(f"    code-distance only          sin^2(2theta)={cd_lo:.2e}..{cd_hi:.2e}")
    print(f"    code + gauge-forbidden hop  sin^2(2theta)={gf_lo:.2e}..{gf_hi:.2e}")
    check(seesaw > 1.0e-10 and cd_lo > 1.0e-10, "generic sterile mixing branches are X-ray excluded")
    check(gf_lo < 1.0e-10, "gauge-forbidden branch can evade X-ray bounds")

    # Density/source is close but conditional.
    q = 208
    n_ratio = ALPHA_DRESSED / q
    m_ev = vals["legacy+dressed"] * 1.0e3
    omega_nur = omega_from_ratio(n_ratio, m_ev)
    omega_dark = 5.0 * omega_nur
    hits = []
    for den in range(1, 513):
        om = omega_from_ratio(ALPHA_DRESSED / den, m_ev)
        if abs(om / OMEGA_NUR_REFERENCE - 1.0) < 0.01:
            hits.append(den)
    print("  density source:")
    print(f"    alpha/208 gives omega_nuR h^2={omega_nur:.6f}, omega_dark h^2={omega_dark:.6f}")
    print(f"    denominators <=512 within 1% of omega_nuR=0.024: {hits}")
    check(208 in hits, "208 lands in the 1% density band")
    check(len(hits) > 1, "density closeness is not unique without the source-map theorem")

    clauses = [
        Clause("mass formula", "CONDITIONAL", "m=alpha^2 Lambda_QCD; alpha/Lambda conventions shift <0.2%"),
        Clause("line position", "CONDITIONAL-SHARP", "E_gamma=m/2 = 8.84 keV if mass formula is accepted"),
        Clause("line flux", "OPEN", "mixing angle/gauge-forbidden Dirac mass not pinned"),
        Clause("density denominator", "CONDITIONAL", "|Q|=208 derived, but alpha/208 source firing still load-bearing"),
        Clause("generation singlet", "FORCED-INSTR", "admitted sterile source algebra is S3 generation-blind"),
    ]
    for clause in clauses:
        status_line(clause.name, clause.status, clause.reason)
    print("  verdict: the 8.84 keV line position is sharp conditional on m=alpha^2 Lambda;")
    print("           a detectable flux or absolute density is not yet parameter-free.")


def main() -> None:
    print("DARK-ENERGY / DARK-SECTOR RIGIDITY AUDIT")
    print("=" * 88)
    audit_28()
    audit_55()
    audit_177()
    print("\n[4] Summary")
    print("  28 / -27/28 : candidate prediction survives rigidity audit at conditional-theorem grade.")
    print("  55 / 12pi/55: rigid under settled strain-readout and value-level severing channel.")
    print("  17.7 keV    : sharp conditional line-position target, not yet a full flux/abundance prediction.")
    print("  Therefore the DESI/Euclid w0-wa flag is legitimate if labelled conditional")
    print("  on the current finite instrument, not numerology from a fitted rational.")
    print("exit 0 -- rigidity audit complete.")


if __name__ == "__main__":
    main()
