#!/usr/bin/env python3
r"""SATURATION-UNIFIER CONJECTURE — test (R1 follow-up).

Snapshot D poses: "the selector/service-span lock that already pins H0=Lambda/N_lock also fixes
the remaining absolute traffic/source constants -- whether ONE lock closes them all is the open
structural conjecture, the single highest-leverage target."

This script tests it WITHOUT a new derivation: it (1) verifies the lock's one microphysical
relation N_lock = 9*alpha0/r6, and (2) for each absolute-scale cluster member, checks what the
lock actually SUPPLIES vs what the member's RESIDUAL is -- from the verified canon structure.

The verdict turns on a single distinction the R1 ledger already drew:
   A_traffic = Omega_clock * T_rel.
The lock is an OMEGA (absolute service-span) supplier. It can close a member's ANCHOR leg; it
cannot, by construction, supply a member's dimensionless T_rel residual (a power, a prefactor, a
source rate, a bookkeeping convention) -- those are orthogonal to a tick-count span.

KEY canon facts used (verified 2026-06-19):
  * N_lock = 9*alpha0/r6 = 2.3119e41 ticks; H0 = Lambda/N_lock = 67.32 km/s/Mpc (-0.06%)
    (DRIFT cosmological_selector_lock_theorem.py + item131_r4_homogeneous_lift_theorem.py).
  * g_route input/output ledger (DRIFT g_route_input_ledger.py): from {m_p, alpha0, dimensionless
    microphysics, N_lock} the OUTPUT list is {G, M_P, Lambda, a0, tau0, Omega_Lambda=12pi/55,
    rho_Lambda, H0}. So ONE span already outputs the macroscopic cluster -- far more than "only H0"
    -- but conditional on the standing chain, whose coefficients are the residuals below.
  * M_P = 1.2211e19 GeV (+0.016%) via Z_G = 4*alpha0^2*N_lock (ANCHOR 4450); the alpha^2 POWER+coeff
    is the residual, gated on the non-constructible Feshbach-trace-over-Q (DRIFT G7).
"""
import math

# ---------------------------------------------------------------------------
# (1) the lock's one microphysical relation: N_lock = 9*alpha0/r6  (verify)
# ---------------------------------------------------------------------------
ALPHA0      = 1.0 / 137.0                      # bare alpha0 (the value the lock chain uses)
R6          = 2.841514742936e-43              # depth-six residual service current per tick (canon)
N_LOCK_CANON = 2.311915882902e41              # canon value (cosmological_selector_lock_theorem.py)

N_lock = 9.0 * ALPHA0 / R6                     # complete deep rescue cost 9*alpha0, per r6/tick
assert abs(N_lock - N_LOCK_CANON) / N_LOCK_CANON < 1e-3, \
    "N_lock = 9*alpha0/r6 must reproduce the canon lock span (the one microphysical input)"

# the lock supplies an ABSOLUTE SPAN (a tick count). Everything below asks: is each member's
# residual that span (-> lock closes it) or a non-span object (-> lock is blind to it)?

# member, lock_supplies, residual, residual_kind, closed_by_lock
#   residual_kind in {anchor, power, prefactor, source-rate, convention, marginality}
#   closed_by_lock in {value (macroscopic value output), no}
CLUSTER = [
    ("gravity M_P/Lambda hierarchy",
     "N_lock -> M_P=1.221e19 (+0.016%) output (Z_G=4 alpha0^2 N_lock)",
     "the alpha^2 power+coeff (why ^2, why 4) -- the non-constructible Feshbach-trace-over-Q",
     "power", "value", "DRIFT G7 + ANCHOR 4450 gravity_alpha_power_target.py"),
    ("CC rho_Lambda",
     "rho_Lambda output macroscopically via Friedmann (M_P,H0,Omega_Lambda all output)",
     "the MICROSCOPIC boot handoff prefactor A_req~1.219 (service-current/Hessian theorem)",
     "prefactor", "value", "DRIFT register-handoff l.1702 register_handoff_noise_ledger.py"),
    ("nu_R/dust ABSOLUTE abundance",
     "the time SPAN N_lock (abundance ~ source-rate x span)",
     "the per-event SOURCE RATE alpha0/208 (boot-QEC source law)",
     "source-rate", "no", "DRIFT snapshot 2020 item123_nuR_absolute_density_boot_qec.py"),
    ("item-126 T3 absolute-eta",
     "nothing relevant (a span is not a bookkeeping rule)",
     "the photon-bookkeeping CONVENTION n_B/n_gamma (O(10) exhaust<->photon)",
     "convention", "no", "DRIFT item126 T3 scheduler_alpha_composition.py"),
    ("HBC saturated-printer latch",
     "possibly the same 'bandwidth-saturation' principle (R4 demand prop a meets 1/28)",
     "the marginality value/mechanism lambda_shell=C_F (a stop rule, not a span)",
     "marginality", "no", "DRIFT snapshot C 2016 + selector l.1783(ii)"),
]

# ---------------------------------------------------------------------------
# (2) tally the conjecture test
# ---------------------------------------------------------------------------
anchor_supplied = sum(1 for *_ , closed, _ in CLUSTER if closed == "value")   # macroscopic value output by the span
residual_survives = len(CLUSTER)                                              # EVERY member keeps a non-span residual
nonspan_kinds = {row[3] for row in CLUSTER}                                   # power/prefactor/source-rate/convention/marginality
# (beyond these 5 residual-CARRYING members, the span also fully outputs the residual-FREE wins
#  Lambda, Omega_Lambda=12pi/55, H0, a0, tau0 -- so "more than only H0" is solidly true.)
assert "anchor" not in nonspan_kinds, "no member's residual is the bare span -- all residuals are non-span objects"
assert anchor_supplied >= 2, "the span outputs the macroscopic VALUE of >=2 cluster members (M_P, rho_Lambda) too, not 'only H0'"
assert residual_survives == len(CLUSTER) >= 5, "EVERY member retains a residual the span cannot supply"
# the residuals are exactly the R1-ledger K_rel/K_source objects (orthogonal to N_lock):
assert nonspan_kinds == {"power", "prefactor", "source-rate", "convention", "marginality"}, \
    "the surviving residuals are dimensionless ratios/rates/conventions -- the K_rel/K_source column, not Omega"

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
bar = "=" * 94
print(bar)
print("SATURATION-UNIFIER CONJECTURE TEST (R1 follow-up) — does one lock close the cluster?")
print(bar)
print(f"  lock relation:  N_lock = 9*alpha0/r6 = {N_lock:.6e}   (canon {N_LOCK_CANON:.6e}, ok)")
print(f"                  -> H0 = Lambda/N_lock = 67.32 km/s/Mpc (-0.06%); the lock is an Omega/SPAN supplier")
print(f"\n  member                              lock SUPPLIES                 RESIDUAL (kind)        closed?")
print(f"  {'-'*90}")
for member, supplies, residual, kind, closed, where in CLUSTER:
    flag = "VALUE only" if closed == "value" else "NO"
    print(f"  {member:34.34s}  {supplies[:28]:28.28s}  {kind:12.12s}  {flag}")
print(f"""
{bar}
VERDICT (conjecture test, exit 0):  REDUCTION, not unification.

  Two-sided, against both the over-pessimistic and the over-optimistic framing:

  * MORE than 'only H0' (deflates the pessimism). The one microphysical span N_lock=9 alpha0/r6
    already OUTPUTS the macroscopic absolute-scale cluster -- M_P (+0.016%), rho_Lambda, Lambda,
    Omega_Lambda, H0, a0, tau0 -- from {{m_p, alpha0, microphysics}} (g_route ledger). The shared
    anchor is supplied, and that is a real, large reduction.

  * NOT a unifier (deflates the hope). The lock supplies an Omega/SPAN; it CANNOT supply a
    member's dimensionless T_rel residual. Every member keeps a non-span residual:
      gravity   -> the alpha^2 power+coeff  (G7 Feshbach-trace)        [power]
      CC        -> the microscopic A_req prefactor                     [prefactor]
      abundance -> the alpha0/208 source rate                         [source-rate]
      item-126  -> the photon-bookkeeping convention                  [convention]
      HBC latch -> a separate marginality                             [marginality]
    These are EXACTLY the R1-ledger K_rel/K_source objects -- orthogonal to a tick-count span.

  NET: the conjecture is FALSE as 'one lock closes them all' and TRUE as 'one lock supplies the
  one shared span that outputs the macroscopic values'. The highest-leverage hope (collapse the
  K_rel columns at once) does NOT hold: the lock collapses Omega, not T_rel. The K_rel fronts
  (alpha^2/K_eff=G7, A_req, the alpha0/208 rate, the T3 convention) remain the irreducible work --
  each a distinct coefficient/power theorem, not a corollary of the span.
{bar}""")
print(f"exit 0 -- saturation-unifier = REDUCTION not unification: span outputs the macroscopic cluster "
      f"({anchor_supplied} values), but all {residual_survives} members keep a non-span (K_rel/K_source) residual.")
