#!/usr/bin/env python3
r"""Audit of the cosmology-pivot proposal #2: "scrub the prefactor, keep rho_Lambda ~ Lambda^3 H_0 as a
bounding theorem; it places dark energy at the correct order of magnitude NATIVELY (~1e-47 GeV^4); the 120-OOM
problem VANISHES because the substrate has no physics above Lambda_QCD."
Tests both concrete claims by direct computation, decomposing the 120 OOM into its pieces. All quoted numbers
asserted. (External inputs are standard PDG/Planck constants, cited; not framework-internal.)
"""
import sys, math

# -- standard constants (PDG/CODATA/Planck; external, cited) --
M_P   = 1.2209e19        # GeV, Planck mass
LAM   = 0.332            # GeV, framework chiral scale (ANCHOR §1.4)
hbar  = 6.582119e-25     # GeV*s
H0    = 67.4 * 1e3 / (3.0857e22) * hbar   # H0=67.4 km/s/Mpc -> s^-1 -> GeV  (Planck 2018)
alpha = 1/137.036
# observed dark-energy density: rho_Lambda^{1/4} ~ 2.3e-3 eV (Omega_L~0.69) -> GeV^4
rho_obs = (2.3e-3 * 1e-9)**4 * 0.69**0   # (2.3e-3 eV in GeV)^4 ; the 2.3e-3 eV already folds Omega_L
rho_obs = (2.3e-3 * 1e-9)**4              # GeV^4  ~ 2.8e-47
print(f"inputs: M_P={M_P:.3e} GeV, Lambda={LAM} GeV, H0={H0:.3e} GeV, rho_obs~{rho_obs:.2e} GeV^4")
assert 2e-47 < rho_obs < 3.5e-47, "observed rho_Lambda ~ 2.5e-47 GeV^4"

# -- CLAIM A: 'rho ~ Lambda^3 H0 lands at the correct OOM (1e-47) NATIVELY (prefactor scrubbed)' --
bare_scale = LAM**3 * H0                  # the proposal's 'bounding theorem', prefactor = 1
ratio = bare_scale / rho_obs
print(f"\n[CLAIM A: scaling places DE at ~1e-47 natively]")
print(f"   Lambda^3 * H0 (prefactor 1) = {bare_scale:.2e} GeV^4   vs observed {rho_obs:.2e} GeV^4")
print(f"   -> overshoots by x{ratio:.0f} = {math.log10(ratio):.1f} orders of magnitude. So 'natively ~1e-47' is FALSE:")
print(f"      the bare scaling sits at ~1e-44, NOT 1e-47. The missing 3.3 OOM ARE the scrubbed prefactor.")
assert 1e3 < ratio < 1e4, "bare Lambda^3 H0 overshoots observed by ~3 OOM"
# the canon prefactor (3/16pi)*alpha closes exactly that gap (item123_lambda_audit):
pref = (3/(16*math.pi))*alpha
print(f"   canon prefactor (3/16pi)*alpha = {pref:.2e};  pref * Lambda^3 H0 = {pref*bare_scale:.2e}  (~ observed)")
print(f"   1/pref = {1/pref:.0f} ~ the x{ratio:.0f} gap -> the scrubbed prefactor IS the 3.3 OOM. Can't drop it AND hit 1e-47.")
assert abs(pref*bare_scale/rho_obs - 1) < 0.5, "the (3/16pi)alpha prefactor reproduces observed to <50%"

# -- CLAIM B: 'the 120-OOM problem VANISHES because no physics above Lambda' --
# decompose the full M_P^4 -> rho_obs gap into (cutoff) + (horizon input) + (prefactor)
gap_total  = math.log10(M_P**4 / rho_obs)                 # the famous ~120
gap_cutoff = math.log10(M_P**4 / LAM**4)                  # M_P^4 -> Lambda^4 : the 'no UV modes' part
gap_horizon= math.log10(LAM**4 / bare_scale)              # Lambda^4 -> Lambda^3 H0 = Lambda/H0 : HORIZON input
gap_pref   = math.log10(bare_scale / rho_obs)             # Lambda^3 H0 -> observed : the scrubbed prefactor
print(f"\n[CLAIM B: 120 OOM 'vanishes' via the Lambda cutoff]  decomposition of the {gap_total:.0f}-OOM gap:")
print(f"   M_P^4 -> Lambda^4  (UV-cutoff at Lambda, 'no physics above Lambda'): {gap_cutoff:6.1f} OOM  <- genuine substrate insight")
print(f"   Lambda^4 -> Lambda^3 H0  (the factor H0/Lambda = horizon input):     {gap_horizon:6.1f} OOM  <- DIRAC COINCIDENCE, inserted via H0, NOT solved")
print(f"   Lambda^3 H0 -> observed  (the scrubbed O(1)/alpha prefactor):         {gap_pref:6.1f} OOM  <- the piece the proposal deletes")
print(f"   {'':36}sum = {gap_cutoff+gap_horizon+gap_pref:.1f} OOM (= the {gap_total:.0f})")
assert abs((gap_cutoff+gap_horizon+gap_pref) - gap_total) < 0.2, "OOM decomposition closes"
assert 75 < gap_cutoff < 82, "the cutoff explains ~78 OOM"
assert 38 < gap_horizon < 45, "the horizon factor H0/Lambda supplies ~41 OOM"
print(f"   => the cutoff explains ~{gap_cutoff:.0f} of {gap_total:.0f} OOM (real). The remaining ~{gap_horizon+gap_pref:.0f} is NOT 'gone':")
print(f"      ~{gap_horizon:.0f} OOM is the horizon H0 put in BY HAND (rho ~ ...*H0) -- i.e. 'why does DE track the horizon',")
print(f"      which IS the coincidence problem; ~{gap_pref:.0f} OOM is the deleted prefactor. So '120 vanishes' overstates:")
print(f"      it splits into ~78 (cutoff, genuine) + ~41 (horizon input, coincidence relocated) + ~3 (scrubbed prefactor).")

print(f"""
=========================================================================================
VERDICT (cosmology-pivot proposal #2: rho_Lambda ~ Lambda^3 H0 'bounding theorem'):
  GENUINE KERNEL: cutting the UV cutoff from M_P to Lambda_QCD ('no continuous physics above the confinement
  scale') legitimately removes the WORST ~78 of the ~123 OOM -- a real qualitative insight, worth leading with.
  Scrubbing the exact O(1)/alpha prefactor as a §16.3 fit is also correct (and already canon: item123 audit,
  62 competitors).
  BUT the two headline claims are FALSE as stated:
   * 'places DE at ~1e-47 natively': NO -- the prefactor-free scaling Lambda^3 H0 = 5e-44, off by 3.3 OOM;
     1e-47 needs exactly the prefactor the proposal deletes. Cannot both scrub it AND claim the right OOM.
   * 'the 120-OOM problem vanishes': only ~78 OOM go away via the cutoff. The other ~41 OOM are supplied by
     the horizon H0 INSERTED in rho ~ Lambda^3 H0 -- the Dirac (Lambda_QCD <-> horizon) coincidence, i.e. the
     coincidence problem relocated, NOT solved. (Same horizon-consuming finding as item123 / DRIFT.)
  HONEST FORM: 'rho_Lambda ~ Lambda^3 H0 is a horizon-INPUT dimensional bound that gets the cutoff-suppression
  (~78 OOM) right from the substrate scale, leaving an order-1e3 prefactor (§16.3) and an ~41-OOM horizon
  factor that is assumed, not derived.' That is a real partial result -- but it is NOT 'the 120 OOM vanish'
  and NOT 'the correct OOM natively'.
=========================================================================================""")
print("exit 0 -- proposal #2: cutoff insight genuine (~78 OOM); but 'native 1e-47' false by 3 OOM (=scrubbed prefactor), and ~41 OOM is horizon input (coincidence relocated, not solved).")
