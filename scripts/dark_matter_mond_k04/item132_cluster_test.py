#!/usr/bin/env python3
r"""ITEM 132 — the CLUSTER test: does the 80/20 split (80% R4 deep-MOND + 20% sterile
nu_R at 17.7 keV) carry the galaxy-CLUSTER mass budget, where pure MOND famously fails?

THE PROBLEM (robust, decades old): MOND applied to galaxy clusters still UNDER-predicts
the dynamical mass. The boost at cluster-core accelerations (g ~ a0) is mild, so MOND
leaves a residual ~ a factor 2 -- i.e. clusters need roughly ONE BARYON-MASS of extra,
unseen, collisionless matter that MOND alone does not supply (Sanders 1999/2003; Angus,
Famaey & Buote 2008; Pointecouteau & Silk 2005). The standard fix is a ~few-keV-to-keV
sterile neutrino with Omega_nu ~ Omega_b ("nuHDM", Angus 2009, m~11 keV).

THE FRAMEWORK'S CLAIM: the 20% sterile nu_R is exactly that fix. This script tests it on
four fronts -- mass budget, Tremaine-Gunn phase space, free-streaming, and the CMB -- and
reports honestly which pass and which are tensions.

exit 0 = the four tests are computed and the solid facts ASSERTED; the verdict (which pass)
is reported. Scope: order-of-magnitude cluster physics + published calibrations, NOT a
cluster-by-cluster fit or an N-body/Boltzmann run (the framework has none -- those remain open).
"""
import math

# ----------------------------------------------------------------------------- [0] constants
OMEGA_DM = 0.264          # Planck 2018 cold/dark matter density
OMEGA_B = 0.0493          # Planck 2018 baryons
OMEGA_M = 0.315
H = 0.674
SPLIT_NU = 0.20           # framework: 20% of the dark sector is particulate nu_R
m_nu_keV = (1/137.035999)**2 * 0.332e6   # alpha^2 Lambda_QCD in keV (= 17.7 keV)
A0 = 1.1e-10              # m/s^2 (framework cH0/2pi)
print(f"[0] Omega_DM={OMEGA_DM}, Omega_b={OMEGA_B} (ratio {OMEGA_DM/OMEGA_B:.2f}); "
      f"m_nu_R={m_nu_keV:.1f} keV; split=80/20")

# ----------------------------------------------------------------------------- [1] MOND cluster residual
def nu_boost(y):                          # g_obs/g_bar = nu(g_bar/a0), framework mu=1-e^-x reading
    return 1.0 / (1.0 - math.exp(-math.sqrt(y)))
# cluster-core baryonic acceleration is ~a0 (M~1e15 Msun within ~1 Mpc): y ~ 1
y_cluster = 1.0
boost = nu_boost(y_cluster)
print("\n[1] MOND CLUSTER RESIDUAL (the problem):")
print(f"    at cluster scale g_bar ~ a0 (y={y_cluster}): MOND boost nu = {boost:.2f} (mild)")
print(f"    observed clusters are ~fair samples: M_dyn/M_b ~ Omega_m/Omega_b = {OMEGA_M/OMEGA_B:.1f}")
residual_per_Mb = 1.0     # literature: MOND leaves a residual ~ 1 baryon-mass (factor-2 total)
print(f"    -> MOND under-predicts; the residual is ~{residual_per_Mb:.0f} baryon-mass of unseen")
print("       collisionless matter per cluster (Sanders 2003; Angus, Famaey & Buote 2008).")

# ----------------------------------------------------------------------------- [2] the nu_R mass budget
def nuR_per_baryon(split):                # fair-sample cluster: M_nuR/M_b = Omega_nuR/Omega_b
    return split * OMEGA_DM / OMEGA_B
supply = nuR_per_baryon(SPLIT_NU)
print("\n[2] nu_R SUPPLY vs the residual (the split's real work):")
print(f"    Omega_nuR = {SPLIT_NU}*Omega_DM = {SPLIT_NU*OMEGA_DM:.3f}  (Omega_nuR/Omega_b = {supply:.2f})")
print(f"    => in a fair-sample cluster nu_R supplies {supply:.2f} baryon-mass vs the ~{residual_per_Mb:.0f} needed")
for s, lab in ((0.05, "95/5"), (0.20, "80/20"), (0.50, "50/50")):
    r = nuR_per_baryon(s)
    tag = "  <-- framework" if abs(s - SPLIT_NU) < 1e-9 else ""
    print(f"      split {lab:>5s}: supplies {r:.2f} M_b  ({'too little' if r<0.5 else 'too much' if r>2 else 'MATCHES residual'}){tag}")
assert 0.5 < supply < 2.0, "80/20 split does NOT supply ~1 baryon-mass"
print("    -> the 80/20 split is SELECTED by the cluster residual: Omega_nuR ~ Omega_b.")

# ----------------------------------------------------------------------------- [3] Tremaine-Gunn phase space
# calibrate on dwarf spheroidals (the strongest bound): m >~ 0.4 keV at rho~0.1 Msun/pc^3, sigma~10 km/s
# (Tremaine & Gunn 1979; Boyarsky, Ruchayskiy & Iakubovskyi 2009). Scale by m_min ~ (rho/sigma^3)^(1/4).
TG_DWARF_keV, RHO_D, SIG_D = 0.4, 0.1, 10.0          # Msun/pc^3, km/s
RHO_CL, SIG_CL = 1.0e-3, 1000.0                       # cluster core (Msun/pc^3, km/s)
def tg_bound(rho, sig):
    return TG_DWARF_keV * ((rho / RHO_D) * (SIG_D / sig) ** 3) ** 0.25
tg_cl, tg_dw = tg_bound(RHO_CL, SIG_CL), TG_DWARF_keV
print("\n[3] TREMAINE-GUNN PHASE-SPACE BOUND (can 17.7 keV pack in?):")
print(f"    cluster bound ~ {tg_cl*1e3:.0f} eV  ; dwarf bound ~ {tg_dw*1e3:.0f} eV  (17.7 keV state)")
print(f"    17.7 keV beats cluster TG by {m_nu_keV/tg_cl:.0f}x  and dwarf TG by {m_nu_keV/tg_dw:.0f}x")
assert m_nu_keV > tg_cl, "fails the cluster phase-space bound"
print("    -> EASILY packs into clusters (cluster fix is phase-space-allowed); also clears the")
print("       dwarf bound -> it CAN clump in galaxies too (flagged in [4]).")

# ----------------------------------------------------------------------------- [4] free-streaming
# thermal half-mode mass M_hm ~ 2.5e10 Msun (m/keV)^-3.33 (Schneider/Lovell-class calibration)
M_hm = 2.5e10 * m_nu_keV ** (-3.33)
M_DWARF_HALO = 1.0e9
print("\n[4] FREE-STREAMING (does it stay OUT of galaxies, as nuHDM needs?):")
print(f"    thermal half-mode mass M_hm ~ {M_hm:.1e} Msun  vs dwarf-halo ~ {M_DWARF_HALO:.0e} Msun")
cold_in_galaxies = M_hm < M_DWARF_HALO
print(f"    -> M_hm << dwarf scale: a 17.7 keV THERMAL relic is COLD, clusters down to ~{M_hm:.0e} Msun")
print("       => it would clump in GALAXIES too (CDM-like), adding ~Omega_nuR/Omega_b ~ 1 baryon-mass")
print(f"          (~{100*supply/(OMEGA_DM/OMEGA_B):.0f}% of galaxy dark mass) ON TOP of the R4 MOND boost --")
print("          a double-count vs the R4-only RAR fit, unless production is non-thermal/warmer")
print("          (the framework's freeze-out velocity distribution is unspecified). TENSION.")

# ----------------------------------------------------------------------------- [5] the CMB
omega_cmb_needed = 0.26      # collisionless density the 3rd acoustic peak requires (~Omega_c)
omega_nuR = SPLIT_NU * OMEGA_DM
print("\n[5] CMB THIRD PEAK (the hard one):")
print(f"    needs collisionless Omega ~ {omega_cmb_needed:.2f}; nu_R provides {omega_nuR:.3f} -> "
      f"{omega_cmb_needed/omega_nuR:.0f}x SHORT")
assert omega_nuR < omega_cmb_needed
print("    -> clusters are fixed by Omega_nuR ~ Omega_b, but the CMB needs ~5x more collisionless")
print("       mass than the nu_R supplies. The CMB completion (R4 as effective collisionless")
print("       source?) is the separate, harder open problem -- no Boltzmann run exists.")

# ----------------------------------------------------------------------------- [6] verdict
print(f"""
[6] VERDICT — the 80/20 split DOES its cluster work, but the same coldness opens two tensions:
  PASS (the real work): Omega_nuR ~ Omega_b ({supply:.2f} M_b) naturally supplies the ~1 baryon-mass
    MOND cluster residual -- the split is SELECTED by clusters (95/5 too little, 50/50 too much),
    and 17.7 keV is comfortably Tremaine-Gunn-safe in clusters. This is a genuine, parameter-light
    success and reproduces the nuHDM cluster fix (Angus 2009) with the framework's own mass+fraction.
  TENSION 1 (galaxies): 17.7 keV is COLDER than canonical nuHDM (M_hm ~ {M_hm:.0e} Msun << dwarf),
    so for thermal-class production it clumps in galaxies too, adding ~20% to galaxy dark mass and
    double-counting against the R4-only RAR. Escape needs a warmer/non-thermal freeze-out spectrum,
    which the framework does not specify. (The X-ray sister-result keeps the line itself safe.)
  TENSION 2 (CMB): the third peak needs ~5x more collisionless mass than the nu_R provides; fixing
    it falls to the R4 sector behaving as effective collisionless mass in the CMB -- unbuilt.
  TIER: cluster mass-budget + TG = order-of-magnitude PASS (the split's defining job); the
    galaxy-coldness and CMB consistency are real, quantified tensions whose resolution needs the
    nu_R production spectrum + an N-body/Boltzmann treatment (neither exists). Net: clusters are
    where the 80/20 split is strongest; galaxies+CMB are where it must still be shown to hold.
exit 0""")
print("ALL ASSERTIONS PASSED — budget match, cluster TG-safe, galaxy-coldness + CMB shortfall quantified.")
