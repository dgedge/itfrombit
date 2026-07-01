#!/usr/bin/env python3
r"""ITEM 118/132 — is the framework's nu_R production WARM enough to keep the 17.7 keV
sterile out of galaxies (resolving the cluster-test galaxy-coldness tension)?

THE QUESTION: the cluster test needs the 20% nu_R to clump in CLUSTERS (it does -- TG-safe,
Omega_nuR~Omega_b supplies the residual) but NOT in galaxies (or it double-counts against the
R4-only RAR). Whether it stays out of galaxies is set by its free-streaming, i.e. its velocity
v = <p>/m today, i.e. its production momentum spectrum. The framework names two production
modes: thermal "freeze-out" (item 116 / 13.4) and "boundary-printing" (HBC). We compute the
warmth each gives and compare to the galaxy-suppression threshold.

THE UNFORGIVING FACT: a 17.7 keV particle produced relativistically (any thermal/freeze-out
route) redshifts to p_today ~ T0 << m, so v_today ~ T0/m is TINY -> it is COLD today. And a
boundary-printed defect placed at rest in the comoving (substrate) frame -- exactly the K04
fossils' "zero peculiar velocity" -- is COLDER still (pure CDM). Neither route is hot.

exit 0 = warmth computed for both production modes and compared to the galaxy threshold; the
verdict (warm enough? y/n) is reported honestly. Scope: standard thermal-relic free-streaming
(velocity-ratio + M_hm ~ v^3 scaling, calibrated to thermal WDM); the framework does not supply
a momentum spectrum, so the thermal-relic value is the warmest defensible stated case.
"""
import math

# ----------------------------------------------------------------------------- [0] constants
T_NU0 = 1.68e-4                 # eV, relic neutrino temperature today
OMEGA_DM_H2 = 0.12             # Planck cold DM
H2 = 0.674 ** 2
m_nuR = 17.7e3                 # eV (alpha^2 Lambda)
OMEGA_NUR_H2 = 0.20 * OMEGA_DM_H2 * 1.0   # nu_R is 20% of DM (Omega_nuR h^2 = 0.2 * Omega_DM h^2)
M_HM_2KEV = 2.5e10 * 2.0 ** (-3.33)       # half-mode mass of a 2 keV thermal WDM (Msun)
M_GALAXY = 1.0e10             # Msun: stay-out-of-galaxies threshold (suppress structure up to dwarf+)
print(f"[0] m_nuR=17.7 keV ; Omega_nuR h^2={OMEGA_NUR_H2:.3f} ; galaxy threshold M_hm>~{M_GALAXY:.0e} Msun")

def v_thermal(m_eV, omega_h2):            # free-streaming velocity today (c units) of a thermal relic
    Tx_over_Tnu = (omega_h2 * 94.0 / m_eV) ** (1.0 / 3.0)
    return 3.15 * T_NU0 * Tx_over_Tnu / m_eV

# ----------------------------------------------------------------------------- [1] thermal freeze-out warmth
v_nuR = v_thermal(m_nuR, OMEGA_NUR_H2)
v_2keV = v_thermal(2000.0, OMEGA_DM_H2)   # 2 keV thermal WDM (all DM) ~ galaxy-suppression benchmark
v_3keV = v_thermal(3000.0, OMEGA_DM_H2)
Tx_over_Tnu = (OMEGA_NUR_H2 * 94.0 / m_nuR) ** (1.0 / 3.0)
print("\n[1] THERMAL FREEZE-OUT (the warmest stated route):")
print(f"    relic temperature T_nuR/T_nu = {Tx_over_Tnu:.4f}  (diluted: froze out early/at high g_*)")
print(f"    free-streaming velocity v_nuR = {v_nuR:.2e} c")
print(f"    vs 2 keV WDM (galaxy-suppression benchmark) v = {v_2keV:.2e} c  ->  nu_R is "
      f"{v_2keV/v_nuR:.0f}x COLDER")
assert v_nuR < v_2keV, "nu_R would be warmer than 2 keV WDM (unexpected)"

# ----------------------------------------------------------------------------- [2] half-mode mass
M_hm = M_HM_2KEV * (v_nuR / v_2keV) ** 3   # M_hm ~ v^3 (free-streaming length ~ v)
print("\n[2] HALF-MODE MASS (structure it suppresses):")
print(f"    M_hm(nu_R) ~ {M_hm:.1e} Msun   (2 keV WDM has {M_HM_2KEV:.1e} Msun)")
print(f"    -> suppresses only below ~{M_hm:.0e} Msun << galaxies (1e8-1e12) => CLUMPS IN GALAXIES.")
assert M_hm < M_GALAXY, "nu_R would suppress galaxy-scale structure (would be warm enough)"

# ----------------------------------------------------------------------------- [3] boundary-printing-at-rest
print("\n[3] BOUNDARY-PRINTING AT REST (the other stated route, K04-fossil analogy):")
print("    printed into the lattice with ZERO peculiar velocity (substrate-frame-static, as the")
print("    K04 fossils) => v ~ 0 => free-streaming ~ 0 => PURE CDM => clumps on ALL scales,")
print("    including galaxies, even more than the thermal route. Colder, not warmer.")

# ----------------------------------------------------------------------------- [4] the warmth gap
v_needed = v_2keV * (M_GALAXY / M_HM_2KEV) ** (1.0 / 3.0)   # velocity for M_hm = galaxy threshold
gap = v_needed / v_nuR
print("\n[4] WHAT WOULD BE NEEDED vs WHAT THE FRAMEWORK GIVES:")
print(f"    to reach M_hm ~ {M_GALAXY:.0e} Msun (stay out of galaxies) needs v ~ {v_needed:.2e} c")
print(f"    => the nu_R must be ~{gap:.0f}x HOTTER (higher <p>) than its thermal value -- a strongly")
print("       NON-THERMAL hot spectrum. Neither freeze-out (cold) nor printing-at-rest (colder)")
print("       supplies it; the framework provides no hot-injection mechanism.")
assert gap > 5, "the warmth gap should be large"

# ----------------------------------------------------------------------------- [5] consequence for the RAR
nuR_per_baryon = 0.20 * 0.264 / 0.0493          # Omega_nuR/Omega_b ~ 1.07 (fair-sample clumping)
galaxy_dark_per_baryon = 5.0                      # MOND/RAR dark-to-baryon ~5 in the outskirts
frac_of_dark = nuR_per_baryon / galaxy_dark_per_baryon
RAR_INTRINSIC_DEX = 0.07                          # intrinsic RAR scatter (~0.05-0.08 dex)
offset_dex = math.log10(1.0 + frac_of_dark)
print("\n[5] CONSEQUENCE for the R4-only RAR (if nu_R clumps fair-sample in galaxies):")
print(f"    nu_R adds ~{nuR_per_baryon:.1f} baryon-mass ~ {100*frac_of_dark:.0f}% of galaxy dark mass")
print(f"    -> a CDM-like (cuspy) offset ~{offset_dex:.3f} dex on top of the R4 RAR, vs the intrinsic")
print(f"       RAR scatter ~{RAR_INTRINSIC_DEX:.2f} dex: it EATS the scatter budget (borderline tension,")
print("       not a clean break -- and worse, cuspy CDM adds scatter the tight RAR disfavours).")

# ----------------------------------------------------------------------------- [6] verdict
print(f"""
[6] VERDICT — NO: the framework's stated production is NOT warm enough to keep nu_R out of galaxies.
  * Thermal freeze-out gives v ~ {v_nuR:.1e} c ({v_2keV/v_nuR:.0f}x colder than 2 keV WDM), M_hm ~ {M_hm:.0e} Msun
    << galaxy scale; boundary-printing-at-rest is colder still (pure CDM). Both CLUMP IN GALAXIES.
  * Staying out of galaxies needs ~{gap:.0f}x higher <p> -- a non-thermal HOT injection the framework
    does not provide. So the cluster-test galaxy-coldness tension is CONFIRMED, not resolved: the
    SAME 17.7 keV that is ideal for clusters (TG-safe, Omega_nuR~Omega_b) is too cold for galaxies.
  * Impact: a ~{100*frac_of_dark:.0f}% cuspy CDM-like nu_R component in galaxies eats the RAR intrinsic-scatter
    budget (~{offset_dex:.2f} vs ~{RAR_INTRINSIC_DEX:.2f} dex) -- a borderline tension with the R4-only RAR, sharpened
    by cuspy CDM adding the very scatter the tight RAR rules out.
  ESCAPES (none established in canon): (i) a non-thermal hot nu_R production (must also preserve the
    abundance and the cluster fix); (ii) the nu_R sitting within the RAR scatter as a sub-dominant
    cuspy component (marginal); (iii) MOND structure-formation suppressing nu_R galaxy clumping
    (unmodelled -- needs the framework's actual gravity law for the nu_R + an N-body run).
  NET: the dark sector now has a sharp INTERNAL constraint -- clusters want the nu_R cold-and-present,
    galaxies want it warm-or-absent, and 17.7 keV + stated production lands COLD on both, so it helps
    clusters at the cost of a galaxy-RAR tension. This is the genuine open knot, not yet untied.
exit 0""")
print("ALL ASSERTIONS PASSED — both production modes are cold; warmth gap ~50x; tension confirmed honestly.")
