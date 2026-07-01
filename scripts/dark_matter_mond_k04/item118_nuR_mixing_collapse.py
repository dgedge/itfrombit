#!/usr/bin/env python3
r"""ITEM 118 — COLLAPSE the 17.7 keV nu_R X-ray bracket to a single mixing angle by
deriving the residual gauge-forbidden Dirac-mass suppression from the framework's own
R4 enforcement scale; turn the 8.84 keV line into a single (falsifiable) XRISM target.

THE MOVE (sidesteps the open item-53 dimensionalisation): the X-ray observable is
  sin^2 2theta with theta = m_D/m_M -- a RATIO, so the unpinned dimensionful base of the
  code amplitudes CANCELS. We need theta, not m_D and m_M separately.

THE FRAMEWORK-NATIVE SUPPRESSION (lines 349 / 742): the neutrino Dirac mass is "governed
  by the R4 Feshbach pole v = 246 GeV" (the EW Higgs VEV = the R4 enforcement scale). For
  the 17.7 keV nu_R that Dirac coupling is GAUGE-FORBIDDEN by R4 (item 118: zero gauge
  charge, equatorial hop forbidden). A coupling forbidden by a symmetry that is enforced
  (restored) at scale v is suppressed, for a state of mass m_M, by the standard approximate-
  symmetry ratio
        theta = m_D / m_M  ~  m_M / v_R4 ,     v_R4 = 246 GeV,
  i.e. the gauge-forbidden mixing is the sterile mass measured in units of the scale of the
  gauge that forbids it. This is the leading dimensional form; it uses ONLY the framework's
  own R4 scale, no new input.

  => theta ~ m_M/v_R4,  m_D ~ m_M^2/v_R4,  sin^2 2theta ~ 4 (m_M/v_R4)^2.

exit 0 = the collapse is computed and confronted with XRISM; the TIER is reported honestly
(this is a leading-dimensional ESTIMATE built on the framework's R4 scale, NOT a microscopic
theorem -- the rigorous matrix element is gated on item 53 (dimensionalisation) + item 87
(mixing operator); alternative gauge-forbiddenness powers are bracketed and X-ray-tested).
"""
import math

# ----------------------------------------------------------------------------- [0] constants
ALPHA = 1.0 / 137.035999
LAMBDA_QCD_eV = 0.332e9
V_R4_eV = 246.0e9                       # EW Higgs VEV = R4 enforcement / Feshbach scale (lines 349/742)
KEV = 1.0e3
m_M_eV = ALPHA**2 * LAMBDA_QCD_eV       # nu_R Majorana mass (item 118)
m_M_keV = m_M_eV / KEV
E_gamma_keV = m_M_keV / 2.0
T_HUBBLE_s = 4.35e17
XRAY_BOUND = 1.0e-10                    # NuSTAR/INTEGRAL-class current bound at 8.85 keV
XRISM_FRONTIER = 1.0e-12               # deep XRISM-era reach (target)

def Gamma_s(s2):                        # Pal-Wolfenstein (BRS19 normalisation)
    return 1.38e-32 * (s2 / 1e-7) * (m_M_keV ** 5)

print("[0] m_M = alpha^2 Lambda = {:.2f} keV ; v_R4 = 246 GeV (R4 enforcement scale) ; "
      "E_gamma = {:.2f} keV".format(m_M_keV, E_gamma_keV))

# ----------------------------------------------------------------------------- [1] the collapse
theta = m_M_eV / V_R4_eV                # gauge-forbidden mixing = m_M / v_R4
m_D_eV = theta * m_M_eV                 # = m_M^2 / v_R4
s2_pred = 4.0 * theta * theta          # sin^2(2theta)
print("\n[1] BRACKET COLLAPSED (theta = m_M/v_R4, the R4-gauge-forbidden ratio):")
print(f"    theta        = m_M/v_R4      = {theta:.3e}")
print(f"    m_D          = m_M^2/v_R4    = {m_D_eV*1e3:.3f} meV   (sub-eV: the 'gauge-forbidden' regime)")
print(f"    sin^2 2theta = 4 theta^2     = {s2_pred:.2e}   <-- SINGLE prediction")
assert abs(m_M_keV - 17.7) < 0.2
assert 1e-15 < s2_pred < 1e-13         # lands in the deep-but-viable band

# ----------------------------------------------------------------------------- [2] lifetime
tau = 1.0 / Gamma_s(s2_pred)
print(f"\n[2] LIFETIME at the predicted mixing: tau = {tau:.2e} s = {tau/T_HUBBLE_s:.1e} x Hubble  (super-stable relic)")
assert tau / T_HUBBLE_s > 1e6

# ----------------------------------------------------------------------------- [3] confront X-ray / XRISM
print("\n[3] CONFRONT THE X-RAY DATA (E_gamma = 8.84 keV):")
print(f"    current bound  sin^2 2theta < {XRAY_BOUND:.0e} : prediction is {XRAY_BOUND/s2_pred:.0f}x below  -> VIABLE")
print(f"    XRISM frontier sin^2 2theta ~ {XRISM_FRONTIER:.0e} : prediction is {XRISM_FRONTIER/s2_pred:.0f}x below  -> "
      f"{'BELOW REACH -> NON-DETECTION PREDICTED' if s2_pred < XRISM_FRONTIER else 'IN REACH'}")
assert s2_pred < XRAY_BOUND, "prediction must respect the current X-ray bound"
detectable = s2_pred >= XRISM_FRONTIER

# ----------------------------------------------------------------------------- [4] alternative powers (honest bracket + X-ray selection)
print("\n[4] ALTERNATIVE gauge-forbiddenness powers (the X-ray bound SELECTS among them):")
alts = {
    "alpha^2 (one forbidden hop)      ": 4 * (ALPHA**2) ** 2,
    "m_M/v_R4 (R4-scale; ADOPTED)     ": s2_pred,
    "alpha^4 (double-forbidden)       ": 4 * (ALPHA**4) ** 2,
}
for name, s2 in alts.items():
    sel = "EXCLUDED by X-ray" if s2 > XRAY_BOUND else ("XRISM-visible" if s2 >= XRISM_FRONTIER else "X-ray-safe, XRISM-invisible")
    print(f"    {name}: sin^2 2theta = {s2:.1e}  -> {sel}")
assert alts["alpha^2 (one forbidden hop)      "] > XRAY_BOUND   # the naive alpha^2 reading is killed
print("    -> the naive alpha^2 reading is X-ray-EXCLUDED; the R4-scale ratio m_M/v_R4 is the")
print("       largest reading that survives, and it lands X-ray-safe but XRISM-invisible.")

# ----------------------------------------------------------------------------- [5] verdict
print(f"""
[5] VERDICT — bracket collapsed to a SINGLE estimate; the X-ray line becomes a falsifiable test:
  * PREDICTION (single number): the 8.84 keV line has sin^2 2theta ~ {s2_pred:.1e}
    (theta ~ m_M/v_R4 = {theta:.1e}, m_D ~ {m_D_eV*1e3:.2f} meV), from the gauge-forbidden
    Dirac mass suppressed by the framework's own R4 enforcement scale v = 246 GeV.
  * FALSIFIABLE CONSEQUENCE: this is ~{XRISM_FRONTIER/s2_pred:.0f}x BELOW the XRISM frontier ->
    the framework predicts XRISM will NOT detect the 8.84 keV line. A detection at
    sin^2 2theta >~ 1e-12 would be ~2 OOM too bright for this reading and would FALSIFY it
    (pointing to a less-suppressed, non-R4-scale Dirac mass). Non-detection confirms it.
  * The X-ray bound also does real work upstream: it EXCLUDES the naive alpha^2
    gauge-forbiddenness reading (sin^2 2theta ~ 1e-8), selecting the deeper R4-scale ratio.
  TIER (honest): leading-DIMENSIONAL ESTIMATE on the framework's R4 scale, not a microscopic
    theorem. The rigorous matrix element is gated on item 53 (the 'applied-by-hand'
    dimensionalisation) + item 87 (the inter-generation mixing operator); and canon still
    carries two nu_R-Majorana mechanisms to reconcile (alpha^2-Lambda, item 118, vs the
    e^{{-2d}} tunnelling of item 53). So: the bracket is collapsed to a single FALSIFIABLE
    estimate (XRISM should see nothing); promoting it to a theorem needs item 53/87.
exit 0""")
print("ALL ASSERTIONS PASSED — collapse computed, alpha^2 reading excluded, XRISM non-detection predicted.")
