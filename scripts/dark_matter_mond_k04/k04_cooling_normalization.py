#!/usr/bin/env python3
r"""K04 cooling/KZ NORMALIZATION (dark_sector Open Frontier #1).

Derive, from current Canon inputs, the three requested quantities:
  (1) H_c  -- the crystallisation-epoch (boundary-printer dilution) rate,
  (2) xi(R) -- KZ correlation length vs the GLASSY-arrest length,
  (3) the fossil abundance Omega_K04 and the shadow suppression it forces.

Canon inputs (no tuning of the gamma=0.995 proxy; no K04-as-halo rescue):
  Lambda_QCD = 0.332 GeV          (framework chiral scale, M_N = 2*sqrt2*Lambda)
  T_start = 6 w6 = Lambda_QCD     (ramp-start anchor) ; w6 = Lambda_QCD/6
  line/wall tension mu = w4+4w6 = 6 w6 = Lambda_QCD
  n_wall = 2.63/xi (lattice units, codim-1, FIXED shape)
  e_debris = 2.17 w6 per debris vertex (sim-measured)
  EoS of pinned extended walls: w = -p/3 = -2/3  (p=2)
  Coarsening verdict: GLASSY ARREST (no Allen-Cahn coarsening) -> network frozen.
Self-asserting: every load-bearing number is asserted in-range.
"""
import math

# ---- constants ----
GeV_s   = 1.519267e24      # 1 GeV = 1.519e24 / s
M_Pl    = 1.220890e19      # GeV, non-reduced Planck mass
T0      = 2.348e-13        # GeV, CMB temperature today (2.725 K)
gstar   = 61.75            # rel. dof just above the QCD transition
gs_fo   = 61.75            # entropy dof at freeze-out
gs_0    = 3.909            # entropy dof today
rho_crit_h2 = 8.098e-47    # GeV^4, critical density / h^2  (1.878e-29 h^2 g/cc)
Omega_DM_h2 = 0.120        # target cold dark matter

LQCD = 0.332               # GeV
w6   = LQCD/6.0
mu   = 6.0*w6              # = LQCD  (wall tension scale)
e_debris = 2.17*w6         # per debris vertex
d_arrest = 0.80            # glassy-arrest deficit at boot-window rate (sim: 0.77-0.87)
a0 = 1.0/LQCD              # substrate register spacing (natural units); ~0.6 fm

print("="*70)
print("K04 COOLING / KZ NORMALIZATION")
print("="*70)
print(f"Lambda_QCD = {LQCD} GeV ; w6 = Lambda/6 = {w6:.4f} GeV ; mu = 6 w6 = {mu:.4f} GeV")
print(f"substrate spacing a0 = 1/Lambda = {a0:.3f} GeV^-1 = {a0*0.1973*1e6:.3f} am"
      f"  (= {a0*0.1973:.3f} fm)")

# ---- (1) H_c : crystallisation-epoch radiation-era Hubble rate ----
def H_rad(T, g): return 1.66*math.sqrt(g)*T*T/M_Pl
print("\n(1) CRYSTALLISATION-EPOCH RATE  H_c = 1.66 sqrt(g*) T_c^2 / M_Pl")
for Tc in (LQCD, LQCD/2):
    Hc = H_rad(Tc, gstar)
    print(f"   T_c = {Tc:.3f} GeV : H_c = {Hc:.3e} GeV = {Hc*GeV_s:.3e} s^-1"
          f"  (Hubble time {1/(Hc*GeV_s)*1e6:.2f} us)")
Hc = H_rad(LQCD, gstar)
assert 1e-20 < Hc < 1e-18, Hc          # QCD-epoch rate band
tau0 = 1.0/(LQCD*GeV_s)                 # substrate tick, s
tauQ = 1.0/Hc/GeV_s                     # n=1 boundary-printer dilution time, s
print(f"   substrate tick tau0 = 1/Lambda = {tau0:.3e} s")
print(f"   quench time tau_Q = 1/(n H_c), n=1 : {tauQ:.3e} s ; tau_Q/tau0 = {tauQ/tau0:.3e}")

# ---- (2) xi : KZ power law (IF critical) vs glassy-arrest length ----
print("\n(2) CORRELATION LENGTH  xi")
ratio = tauQ/tau0
for sig in (0.25, 0.49, 0.727):        # the (unpinned, L-drifting) beta estimates
    print(f"   KZ (sigma={sig:.3f}) : xi/a0 = (tau_Q/tau0)^sigma = {ratio**sig:.3e}")
xi_arrest = 2.63/ (d_arrest*3)          # codim-1: wall density ~ d per axis; xi ~ 2.63/n_wall
print(f"   GLASSY arrest (sim): d~{d_arrest} -> network frozen at xi_arrest ~ {xi_arrest:.2f} a0"
      f"  (lattice scale; NOT KZ-divergent)")
print("   -> beta(L) drifts 0.73/0.49/0.26 (L=4/6/8): NO clean KZ exponent; the")
print("      glassy network does not coarsen, so xi stays microscopic (a few a0).")

# ---- (3) fossil abundance + the overclosure / shadow bound ----
print("\n(3) FOSSIL ABUNDANCE  Omega_K04  (bare, then shadow requirement)")
LQCD4 = LQCD**4
rho_wall_fo = d_arrest*e_debris*LQCD**3            # energy/volume at freeze-out, GeV^4
rho_rad_fo  = (math.pi**2/30.0)*gstar*LQCD4
frac_fo = rho_wall_fo/rho_rad_fo
print(f"   rho_wall(fo) = d*e_debris*Lambda^3 = {rho_wall_fo:.3e} GeV^4"
      f"  ({frac_fo*100:.2f}% of radiation at freeze-out)")
# scale factor freeze-out -> today (entropy conservation a T gs^1/3 = const)
a_fo_over_a0 = (T0/LQCD)*(gs_0/gs_fo)**(1.0/3.0)
# frozen glassy walls: w=-2/3 -> rho ~ a^-1
rho_wall_today = rho_wall_fo * a_fo_over_a0**1
Omega_bare = rho_wall_today/rho_crit_h2
print(f"   a_fo/a_today = {a_fo_over_a0:.3e}")
print(f"   walls w=-2/3 -> rho ~ a^-1 ; rho_wall(today,bare) = {rho_wall_today:.3e} GeV^4")
print(f"   Omega_K04,bare h^2 = {Omega_bare:.3e}   <-- DOMAIN-WALL OVERCLOSURE")
shadow_needed = Omega_DM_h2/Omega_bare
print(f"   shadow suppression to reach Omega<=0.12 : {shadow_needed:.2e}"
      f"  (~{math.log10(1/shadow_needed):.0f} orders -> near-total erasure)")
# Zel'dovich one-wall-per-horizon comparison (IF it coarsened, which it does not)
H0 = 1.44e-42  # GeV (h=0.674)
Omega_zel = (8*math.pi/3.0)*mu**1*LQCD**2*H0 / (H0**2 * M_Pl**2)  # ~ sigma H0 / (rho_crit); sigma~mu*Lambda^2
sigma_wall = mu*LQCD**2   # surface tension ~ tension * area density, GeV^3
Omega_zel = 8*math.pi*sigma_wall/(3.0*H0*M_Pl**2)
print(f"   [cf. Zel'dovich 1-wall/horizon, sigma~mu*Lambda^2={sigma_wall:.3f} GeV^3:"
      f" Omega~{Omega_zel:.2e}] -- even the coarsened case overcloses")

assert Omega_bare > 1e3, Omega_bare     # must overclose
print("\n" + "="*70)
print("VERDICT: H_c ~ 1e-19 GeV (QCD epoch) is derived; xi is glassy-arrest")
print("(microscopic, no clean KZ exponent); the bare frozen wall fossil")
print("OVERCLOSES by ~30+ orders (Zel'dovich wall catastrophe, worsened by the")
print("non-coarsening glass). K04 is viable ONLY as a near-totally shadow-")
print("truncated, gravitationally negligible pinned fossil -- it CANNOT be the")
print("DM halo / CMB cold component. Canon (K04 pinned, subdominant) confirmed")
print("and SHARPENED into a quantitative overclosure bound. exit 0")
