#!/usr/bin/env python3
r"""ITEM 118 — the 17.7 keV sterile nu_R X-ray decay line: pin what CAN be pinned
(line POSITION, stability, and the mixing-angle BRACKET) and confront XRISM/NuSTAR.

THE PARAMETER-FREE PART (pinned):
  m_s = alpha^2 * Lambda_QCD = 17.7 keV (item 118; alpha=1/137.036, Lambda=0.332 GeV).
  Radiative decay nu_s -> nu + gamma is two-body => a MONOCHROMATIC line at
      E_gamma = m_s/2 = 8.85 keV  (no free parameter).
  Decay rate (Pal-Wolfenstein, the convention the X-ray bounds are quoted in;
  Boyarsky-Ruchayskiy-Shaposhnikov 2019):
      Gamma = 1.38e-32 s^-1 * (sin^2 2theta / 1e-7) * (m_s/keV)^5.

THE UNPINNED PART (and why): the LINE FLUX scales with the active-sterile mixing
  sin^2(2theta), which item 118 says is GAUGE-FORBIDDEN (the nu_R is the zero-syndrome
  decoupled vertex; the 'equatorial-hop' Dirac mass is forbidden). The framework does
  NOT give a single mixing value (the neutrino-mixing sector is open: item 87, and the
  seesaw dimensionalisation is applied by hand). So we BRACKET sin^2(2theta) from the
  framework's own suppression laws and compare to the X-ray bound, rather than invent
  one number. Crucially, the 17.7 keV nu_R is NOT the active-mass seesaw partner
  (that is the M_R~1e15 GeV state, ANCHOR Q6/item 53: m_nu = v^2/M_R), so its mixing is
  NOT fixed by the active-neutrino seesaw line.

THE BRACKET (each a framework-grounded suppression of the Dirac mass m_D; theta=m_D/m_s):
  (a) seesaw line (IF it WERE the active partner): sin^2 2theta = 4 m_nu/m_s ~ 1e-5;
  (b) code-distance tunnelling alone (item ~53/Part 5: amplitude e^{-d}..e^{-2d}, d=4
      the [8,4,4] distance): sin^2 2theta ~ 1e-7 .. 1e-3;
  (c) code-distance PLUS the item-118 gauge-forbidden alpha^2 hop (extra alpha^2 on the
      amplitude => alpha^4 on sin^2 2theta): sin^2 2theta ~ 1e-15 .. 1e-11.

exit 0 = machinery verified (line position + stability + seesaw exclusion ASSERTED);
the physics verdict (which bracket survives the X-ray bound) is reported honestly.
Scope: representative published bounds (NuSTAR/INTEGRAL ~1e-10; XRISM frontier ~1e-11..-12),
not a re-reduction of X-ray data; the bracket is model-dependent on the OPEN item-87 sector.
"""
import math

# ----------------------------------------------------------------------------- [0] constants
ALPHA = 1.0 / 137.035999
LAMBDA_QCD_eV = 0.332e9                     # canon Lambda_QCD = 0.332 GeV
KEV = 1.0e3
m_s_eV = ALPHA**2 * LAMBDA_QCD_eV          # item 118: m_s = alpha^2 Lambda
m_s_keV = m_s_eV / KEV
E_gamma_keV = m_s_keV / 2.0                 # two-body decay line
M_NU_eV = 0.05                              # heaviest active nu (sqrt(dm2_atm), NH)
D_CODE = 4                                  # [8,4,4] code distance
T_HUBBLE_s = 4.35e17                        # age of universe ~13.8 Gyr
# representative published X-ray bounds at E_gamma ~ 8.85 keV (sin^2 2theta):
XRAY_BOUND = 1.0e-10                        # NuSTAR/INTEGRAL-class, conservative current
XRAY_FRONTIER = 1.0e-12                     # deep / XRISM-era frontier (target)

print("[0] inputs: alpha=1/137.036, Lambda_QCD=0.332 GeV, m_nu=0.05 eV, d_code=4")

# ----------------------------------------------------------------------------- [1] line position (pinned)
print("\n[1] PARAMETER-FREE LINE POSITION:")
print(f"    m_s   = alpha^2 Lambda_QCD = {m_s_keV:.2f} keV")
print(f"    E_gamma = m_s/2            = {E_gamma_keV:.2f} keV   (monochromatic; no free parameter)")
assert abs(m_s_keV - 17.7) < 0.2, "m_s != 17.7 keV"
assert abs(E_gamma_keV - m_s_keV / 2) < 1e-9
print("    -> a sterile-nu line from this state can ONLY sit at 8.85 keV; a detection")
print("       elsewhere would falsify the mass identification (this part IS sharp).")

# ----------------------------------------------------------------------------- [2] decay rate + stability
def Gamma_s(sin2_2th):                      # s^-1 (Pal-Wolfenstein; BRS19 normalisation)
    return 1.38e-32 * (sin2_2th / 1e-7) * (m_s_keV ** 5)
def tau_over_Hubble(sin2_2th):
    return 1.0 / Gamma_s(sin2_2th) / T_HUBBLE_s
print("\n[2] LIFETIME (must be >> Hubble for a viable DM relic):")
for s2 in (1e-7, XRAY_BOUND, XRAY_FRONTIER):
    print(f"    sin^2 2theta={s2:.0e}: tau = {1/Gamma_s(s2):.2e} s = {tau_over_Hubble(s2):.1e} x Hubble")
assert tau_over_Hubble(XRAY_BOUND) > 1e6, "not stable even at the X-ray bound"
print("    -> stable by >=1e8 Hubble times at/below the X-ray bound: a viable relic regardless.")

# ----------------------------------------------------------------------------- [3] X-ray bound
print(f"\n[3] X-RAY BOUND at E_gamma=8.85 keV: sin^2 2theta < {XRAY_BOUND:.0e} (NuSTAR/INTEGRAL),")
print(f"    deep/XRISM frontier ~ {XRAY_FRONTIER:.0e}.")

# ----------------------------------------------------------------------------- [4] mixing-angle bracket
def s2(theta):  # sin^2(2theta) ~ 4 theta^2 for small theta
    return 4.0 * theta * theta
seesaw = 4.0 * M_NU_eV / m_s_eV                                   # (a)
amp_2d, amp_1d = math.exp(-2 * D_CODE), math.exp(-D_CODE)         # code-distance amplitudes
cd_lo, cd_hi = s2(amp_2d), s2(amp_1d)                             # (b)
cdgf_lo, cdgf_hi = cd_lo * ALPHA**4, cd_hi * ALPHA**4            # (c) + gauge-forbidden alpha^2 on amp
print("\n[4] MIXING-ANGLE BRACKET (framework-grounded suppressions; theta=m_D/m_s):")
print(f"    (a) seesaw line 4 m_nu/m_s ............ sin^2 2theta = {seesaw:.1e}"
      f"   -> {seesaw/XRAY_BOUND:.0e}x the bound : EXCLUDED")
print(f"    (b) code-distance e^-d..e^-2d (d=4) ... sin^2 2theta = {cd_lo:.1e} .. {cd_hi:.1e}"
      f"   -> EXCLUDED ({cd_lo/XRAY_BOUND:.0e}x .. {cd_hi/XRAY_BOUND:.0e}x bound)")
print(f"    (c) code-distance + gauge-forbidden a^2  sin^2 2theta = {cdgf_lo:.1e} .. {cdgf_hi:.1e}"
      f"   -> {'VIABLE' if cdgf_lo < XRAY_BOUND else 'EXCLUDED'} (straddles the XRISM frontier)")
# the decisive structural statement: (a),(b) are X-ray excluded; only (c) survives.
assert seesaw > XRAY_BOUND and cd_lo > XRAY_BOUND, "expected (a),(b) above the bound"
assert cdgf_lo < XRAY_BOUND, "expected the gauge-forbidden bracket to dip below the bound"
# what suppression the bound demands of the Dirac mass:
theta_max = math.sqrt(XRAY_BOUND / 4.0)
mD_max_eV = theta_max * m_s_eV
print(f"    => the bound DEMANDS m_D = theta*m_s < {mD_max_eV:.2f} eV (theta < {theta_max:.1e}),")
print(f"       i.e. a Dirac mass {LAMBDA_QCD_eV*ALPHA/ (mD_max_eV):.0e}x below even alpha*Lambda --")
print("       sub-eV, far below any EW/equatorial Dirac mass: exactly 'gauge-forbidden'.")

# ----------------------------------------------------------------------------- [5] verdict
print(f"""
[5] VERDICT (machinery asserted; flux verdict honest):
  * PINNED, PARAMETER-FREE: the line sits at E_gamma = {E_gamma_keV:.2f} keV (m_s=alpha^2 Lambda),
    and the relic is stable by >=1e8 Hubble times for any X-ray-allowed mixing.
  * NOT PINNED: the line FLUX. The framework gives a GAUGE-FORBIDDEN Dirac mass, not a
    number, and the neutrino-mixing sector is open (item 87). So sin^2 2theta is bracketed,
    not predicted.
  * THE DECISIVE STRUCTURAL RESULT: a GENERIC 17.7 keV sterile neutrino is X-ray-EXCLUDED --
    the seesaw value (1e-5) by ~5 OOM, even bare code-distance tunnelling (1e-7..1e-3) by
    3-7 OOM. The framework SURVIVES only because item 118's Dirac mass is gauge-forbidden:
    the bound forces m_D < {mD_max_eV:.2f} eV, and the gauge-forbidden alpha^2 hop supplies
    exactly that, landing sin^2 2theta ~ 1e-15..1e-11 -- at/below the XRISM frontier.
  * Two genuine escapes from the standard keV-sterile X-ray death, both structural (not
    tuned): (i) production is freeze-out/boundary-printing, NOT Dodelson-Widrow, so the
    abundance is decoupled from the mixing; (ii) the Dirac mass is gauge-forbidden, so the
    mixing is far below the seesaw line.
  HONEST STATE: the X-ray line is X-ray-SAFE-BY-STRUCTURE, not a positive prediction --
    its POSITION (8.85 keV) is sharp, its FLUX is unpredicted until item 87 pins the
    residual Dirac-mass suppression. The optimistic end of the bracket (~1e-11) is a real
    near-term XRISM target; the pessimistic end (~1e-15) is invisible. So: 'derive the
    mixing angle' is GATED on item 87 -- but the X-ray bound now sharply CONSTRAINS that
    open sector (it must suppress m_D below ~0.1 eV), which is a concrete forward target.
exit 0""")
print("ALL ASSERTIONS PASSED — line position + stability + seesaw/code-distance exclusion verified.")
