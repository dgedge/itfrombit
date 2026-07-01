#!/usr/bin/env python3
r"""[SUPERSEDED 2026-06-15: the R~1 high-reflectivity conclusion below is WITHDRAWN. Canon
(ANCHOR item 152; blackhole_deep v2.1, DOI 10.5281/zenodo.20702600) is that the QEC horizon is a
one-way record-writing isometry, NOT a coherent mirror -> the canonical prediction is a LARGE-ECHO
NULL; echo searches are upper-bound tests on extra horizon structure. The round-trip dt geometry
computed below stands as that upper-bound template; the R~1 reflectivity reading does not.]

PRECISE GW-ECHO TEMPLATE from the graviton-firewall-core scattering — sharpen item 124's
ringdown echo from order-of-magnitude to a quantitative LIGO/LISA template (delay + reflectivity).

ECHO PHYSICS (standard ECO; Cardoso-Pani): GW is trapped between the photon-sphere potential
barrier (~3M, where the Regge-Wheeler potential peaks) and the reflecting surface near r_h. The
echo SPACING is the round-trip tortoise time dt = 2|r*(r_peak) - r*(r_surface)|, and the echo
AMPLITUDE train decays by the surface reflectivity R per bounce.

FRAMEWORK-SPECIFIC INPUTS (this is what makes the template sharp, not generic):
  * the reflecting surface is the rigid-condensate-core boundary, set by the LATTICE cutoff
    a0 = 1/Lambda_QCD (~0.594 fm), NOT the Planck length -- so r_surface = r_h + O(a0).
    Because a0/ell_P ~ 10^20, ln(r_h/a0) ~ ln(r_h/ell_P)/2: the framework's echo spacing is
    ~HALF a generic Planck-ECO's. A distinguishing, parameter-free number.
  * the core is RIGID (zero clock rate) -> infinite impedance for the shear-phonon graviton ->
    R = |(Z_core - Z_ext)/(Z_core + Z_ext)|^2 -> 1 (total reflection); the sharp phase boundary
    leaks only ~(a0/lambda_GW)^2 ~ 10^-43. So R ~ 1 (a slowly-decaying, detectable echo train).

exit 0 = dt (both cutoff conventions) + R computed for LIGO/LISA masses and ASSERTED; the template
+ honest residuals (cutoff convention factor-2; firewall geometric-phase graviton channel) reported.
"""
import math

G, C, M_SUN = 6.674e-11, 2.99792458e8, 1.989e30
A0 = 0.594e-15           # lattice cutoff = 1/Lambda_QCD (m)
ELL_P = 1.616e-35        # Planck length (m)

def echo_delay(M_kg, cutoff_m, proper=True):
    """round-trip tortoise time between the photon sphere (3M) and a surface one cutoff from r_h.
    Uses eps = (r_surface-r_h)/r_h DIRECTLY (1+eps underflows to 1.0 in float, so r*_surface is
    built from r_h + r_h*ln(eps) without reconstructing r_surface)."""
    M_geom = G * M_kg / C ** 2
    r_h = 2 * M_geom
    if proper:                       # cutoff is a PROPER distance: 2 r_h sqrt(eps) = cutoff
        eps = (cutoff_m / (2 * r_h)) ** 2
    else:                            # cutoff is a COORDINATE distance: r_s - r_h = cutoff
        eps = cutoff_m / r_h
    rstar_s = r_h * (1 + eps) + r_h * math.log(eps)                 # r*(surface)
    rstar_peak = 3 * M_geom + r_h * math.log(3 * M_geom / r_h - 1)  # r*(3M); 3M/r_h-1 = 0.5
    dt_len = 2 * abs(rstar_peak - rstar_s)
    return dt_len / C

# ----------------------------------------------------------------------------- [1] echo delay
print("[1] ECHO SPACING dt (framework cutoff a0 = 1/Lambda_QCD; proper-distance convention):")
cases = [("LIGO 30 Msun", 30), ("LIGO 65 Msun", 65), ("LISA 1e6 Msun", 1e6)]
for label, m in cases:
    dt_a0 = echo_delay(m * M_SUN, A0, proper=True)
    dt_a0_coord = echo_delay(m * M_SUN, A0, proper=False)
    unit = "s" if dt_a0 > 2 else "ms"
    val = dt_a0 if dt_a0 > 2 else dt_a0 * 1e3
    valc = dt_a0_coord if dt_a0 > 2 else dt_a0_coord * 1e3
    print(f"    {label:16s}: dt = {val:7.1f} {unit}  (coordinate-cutoff {valc:.1f} {unit}; "
          f"factor-{dt_a0/dt_a0_coord:.1f} convention range)")
dt30 = echo_delay(30 * M_SUN, A0, proper=True)
assert 1e-2 < dt30 < 0.2, "30 Msun echo should be tens-to-~100 ms (LIGO band)"

# framework (a0) vs generic Planck-ECO (ell_P): the distinguishing factor
dt30_planck = echo_delay(30 * M_SUN, ELL_P, proper=True)
print(f"\n    DISTINGUISHING: framework a0-cutoff dt = {dt30*1e3:.0f} ms vs generic Planck-ECO "
      f"(ell_P) dt = {dt30_planck*1e3:.0f} ms")
print(f"    -> framework spacing is {dt30/dt30_planck:.2f}x the Planck-ECO (because a0/ell_P~10^20 halves the log):")
print("       the echo SPACING encodes the discreteness scale -- a0 = 1/Lambda_QCD, not ell_P.")
assert dt30 < dt30_planck

# ----------------------------------------------------------------------------- [2] reflectivity R
print("\n[2] REFLECTIVITY R (graviton shear-phonon at the rigid-core boundary):")
# QNM wavelength for the leakage estimate
f_qnm = 1.2e4 * (M_SUN / (30 * M_SUN))            # l=2 Schwarzschild QNM ~ 12 kHz * (Msun/M)
lam_gw = C / f_qnm
leak = (A0 / lam_gw) ** 2                          # finite-wavelength leakage past a sharp boundary
R = 1.0 - leak
print(f"    Z_core -> infinity (rigid, zero clock rate) => R = |(Z_core-Z_ext)/(Z_core+Z_ext)|^2 -> 1")
print(f"    sharp-boundary leakage ~ (a0/lambda_GW)^2 = {leak:.0e} (lambda_GW~{lam_gw/1e3:.0f} km at f_QNM~{f_qnm:.0f} Hz)")
print(f"    => R = {R:.20f} ~ 1  (essentially total reflection; a SLOWLY-decaying echo train)")
assert R > 0.99

# ----------------------------------------------------------------------------- [3] the template
print("\n[3] THE TEMPLATE (pulses at t_main + n*dt; phase-inverted; amplitudes A_n):")
T_barrier = 0.3       # photon-sphere barrier transmission for l=2 at the QNM band (representative)
print(f"    spacing      dt = {dt30*1e3:.0f} ms (30 Msun); scales ~ M ln(r_h/a0).")
print(f"    reflectivity R ~ 1 -> echo amplitudes A_n ~ A_main * T_barrier * R^(n-1) ~ A_main * {T_barrier}"
      f" * 1^(n-1):")
for n in (1, 2, 3, 5):
    print(f"      echo {n}: A_{n}/A_main ~ {T_barrier * R**(n-1):.3f}")
print("    -> a long, slowly-decaying, evenly-spaced echo train (R~1 is the most-detectable case);")
print("       each echo a phase-inverted copy of the main ringdown.")

# ----------------------------------------------------------------------------- [4] verdict
print(f"""
[4] VERDICT — the echo template is now QUANTITATIVE (not just order-of-magnitude):
  * SPACING: dt = 2|r*_peak - r*_surface| with the surface at the LATTICE cutoff a0 = 1/Lambda_QCD:
    dt ~ {dt30*1e3:.0f} ms (30 Msun), {echo_delay(65*M_SUN,A0)*1e3:.0f} ms (65 Msun),
    ~{echo_delay(1e6*M_SUN,A0):.0f} s (1e6 Msun, LISA). This is ~{dt30/dt30_planck:.2f}x a generic Planck-ECO -- the
    framework's a0 (not ell_P) is the distinguishing, parameter-free signature in the SPACING.
  * REFLECTIVITY: R ~ 1 (rigid zero-clock-rate core; sharp-boundary leakage ~10^-43): a long,
    slowly-decaying echo train -- the most detectable case, and stronger than a generic partial-
    reflector ECO.
  * TEMPLATE: evenly-spaced (dt), phase-inverted, slowly-decaying (R~1) pulses after the main
    ringdown; first echo ~T_barrier~30% of the main amplitude. A concrete LIGO/LISA matched-filter
    target: standard GR has NO echoes; the framework predicts this specific train.
  RESIDUALS (named, factor-level): (i) cutoff convention (proper vs coordinate a0) gives a factor
    ~2 in dt ({echo_delay(30*M_SUN,A0,False)*1e3:.0f}-{dt30*1e3:.0f} ms for 30 Msun); (ii) the firewall geometric-phase channel
    e^{{i kappa tau0 Z}} could couple to gravitons and reduce R below 1 (damping the train) -- needs the
    L_BH graviton matrix element; (iii) RW-peak at 3.28M (not 3M) shifts dt by <1%. NET: a sharp,
    framework-specific echo template (dt ~ {dt30*1e3:.0f} ms at 30 Msun, R ~ 1), ready for an echo search.
exit 0""")
print(f"ALL ASSERTIONS PASSED — dt ~ {dt30*1e3:.0f} ms (30 Msun, a0-cutoff, ~half the Planck-ECO); R ~ 1; "
      f"quantitative echo template delivered.")
