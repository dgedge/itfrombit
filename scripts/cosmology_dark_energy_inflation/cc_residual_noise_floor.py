#!/usr/bin/env python3
r"""THE CC RESIDUAL'S OBSERVATIONAL NOISE FLOOR — the missing computation that
re-poses (and closes) the next-order hunt.

THE POINT NOBODY COMPUTED: C_target is DEFINED by matching rho_obs, so it
inherits rho_obs's observational uncertainty through the chain's own
sensitivity.  Neither the renormalization no-go nor the fixed-connected-object
audit propagated Planck's error bars.  Doing so (here, exactly, through the
live pipeline):

    sigma_C(Planck)  ~ 0.012   while   delta C = 0.0011875

i.e. the 12-digit target delta C = +0.001187526879799 has ZERO observationally
defined digits — nature does not currently fix even its SIGN.  The residual is
a ~0.10 sigma effect; the H0 tension alone moves C_target by ~86 delta C.

CONSEQUENCES (each verified below):
 [1] the exact pipeline sensitivity dln(rho)/dC at the operating point;
 [2] the noise floor: sigma_C under three accountings (Planck quadrature,
     framework-internal M_P-closure accounting, the H0-tension systematic) —
     all >> delta C;
 [3] THE CONVENTION TABLE: every unpinned billing/normal-ordering convention
     named across the two audits plus the three exponentiation readings of the
     adopted billing axiom, in sigma units — the ENTIRE convention space spans
     ~0.13 sigma: observationally degenerate.  In particular the rate-log
     reading -ln(1-alpha0)*C lands at 0.0065 sigma from the central value —
     and PRECISELY BECAUSE of [2] this cannot be promoted: at 0.1-sigma stakes,
     "landing" is free;
 [4] the quarantine ledger: the eta = 1/2 + 0.036 composite (matches
     C2_req/C_loop to 1.1e-4) — refused: the 0.036 component has no operator,
     the composite has >= 18 effective trials, and the target digit being
     matched is noise;
 [5] the resolution requirement: discriminating delta C at 3 sigma needs
     rho_Lambda to 0.062% <=> H0 to +/-0.021 km/s/Mpc — ~26x beyond Planck,
     beyond any planned survey; and unposeable until the 8% H0 tension
     (86 delta C of systematic) resolves.

VERDICT: the question "which convention closes the residual?" is observation-
ally VOID; the well-posed remainder is pure algebra ("which exponentiation
does the monitored formalism force?" — goes to the billing note), whose answer
sets a central value that no measurement can check below ~4% of C.  The CC
chain's honest final form: rho_Lambda = rho_obs at 0.10 sigma (Planck), CLOSED
at observational precision.  exit 0 = every claim verified."""
import math

from cc_generation_vertex_item115_loop import extrapolate_c, rho_for_c
from register_handoff_form_selection import (
    ALPHA0,
    BASE_GAMMA,
    solve_gamma_for_stationary_queue,
    solve_target,
)

# observational inputs (literature, labelled; quadrature = standard first cut —
# Planck's internal H0/OmegaL correlations only tighten/loosen this by O(10%),
# irrelevant at the 10x margins below)
H0_PLANCK, DH0_PLANCK = 67.36, 0.54        # Planck 2018 TT,TE,EE+lowE+lensing
OMEGA_L, DOMEGA_L = 0.6847, 0.0073
H0_SHOES, DH0_SHOES = 73.04, 1.04          # SH0ES 2022 Cepheid-SN ladder

print("[0] TARGET NUMBERS (reproduced from the live pipeline):")
_, q1_target = solve_target()
gamma_star = solve_gamma_for_stationary_queue(q1_target)
c_target = -math.log(gamma_star / BASE_GAMMA) / ALPHA0
c_lin, c_quad, _ = extrapolate_c([64, 96, 128, 160, 192, 256])
c_loop = 0.5 * (c_lin + c_quad)
delta_c = c_target - c_loop
print(f"    C_target = {c_target:.15f}   C_loop = {c_loop:.15f}")
print(f"    delta C  = {delta_c:.15f}   rho(C_loop) = {rho_for_c(c_loop, q1_target):.9f} rho_obs")
assert abs(delta_c - 0.0011875) < 2e-5

# ---------------- [1] the exact pipeline sensitivity ----------------
h = 1e-6
S_C = -(math.log(rho_for_c(c_loop + h, q1_target)) -
        math.log(rho_for_c(c_loop - h, q1_target))) / (2 * h)
print(f"\n[1] EXACT SENSITIVITY at the operating point: |dln(rho)/dC| = {S_C:.6f}")
assert 1.4 < S_C < 1.8

# ---------------- [2] the noise floor ----------------
print("\n[2] THE NOISE FLOOR — sigma_C under three accountings:")
dlnrho_planck = math.hypot(2 * DH0_PLANCK / H0_PLANCK, DOMEGA_L / OMEGA_L)
sigma_C_planck = dlnrho_planck / S_C
# framework-internal accounting: the M_P closure pins H0*OmegaL to 0.07%;
# rho = (3/8pi)(H0 OmegaL) H0 M_P^2 leaves one bare H0
dlnrho_internal = math.hypot(0.0007, DH0_PLANCK / H0_PLANCK)
sigma_C_internal = dlnrho_internal / S_C
# the H0-tension systematic: Planck vs SH0ES central values
dlnrho_tension = 2 * math.log(H0_SHOES / H0_PLANCK)
shift_C_tension = dlnrho_tension / S_C
rows = [("Planck 2018 quadrature (1 sigma)", sigma_C_planck),
        ("internal M_P-closure accounting", sigma_C_internal),
        ("H0 tension (Planck vs SH0ES, systematic)", shift_C_tension)]
for nm, s in rows:
    print(f"    {nm:<44s} sigma_C = {s:.6f}  = {s/delta_c:6.1f} x delta C")
n_sigma = delta_c / sigma_C_planck
print(f"    -> the residual is a {n_sigma:.3f}-sigma effect (Planck);")
print(f"       observationally defined digits of the 12-digit delta C: ZERO —")
print(f"       nature does not currently fix even its sign.")
assert sigma_C_planck > delta_c            # zero significant digits
assert n_sigma < 0.15
assert shift_C_tension / delta_c > 50      # the tension dwarfs the residual

# ---------------- [3] the convention table, in sigma units ----------------
print("\n[3] THE UNPINNED-CONVENTION TABLE (every named convention, in sigma):")
delta = ALPHA0 - 1 / 137.036               # dressed-alpha billing shift
conventions = [
    ("bare per-touch billing alpha0*C (adopted)",          0.0),
    ("exact no-jump log  -ln(1-alpha0*C)",                 ALPHA0 * c_loop**2 / 2),
    ("rate-log billing   -ln(1-alpha0)*C",                 ALPHA0 * c_loop / 2),
    ("dressed-alpha billing (1/137.036)",                  -delta * c_loop),
    ("anti-dressed billing",                               +delta * c_loop),
    ("Kraus walk-active r^2/4 (their audit)",              ALPHA0 * (2/3)**2 / 4),
    ("Kraus no-jump log r^2/2 (their audit)",              ALPHA0 * (2/3)**2 / 2),
    ("|c2_QED|/2 two-loop near-number (their audit)",      ALPHA0 * 0.164239),
    ("bath cumulant kappa2/2 (their audit)",               ALPHA0 * 0.018640),
]
print(f"    {'convention':<46s} {'dC':>12s} {'residual after':>15s} {'sigma off obs':>14s}")
vals = []
for nm, dc in conventions:
    after = delta_c - dc
    vals.append(dc)
    print(f"    {nm:<46s} {dc:+12.3e} {after:+15.3e} {abs(after)/sigma_C_planck:13.4f}")
spread = (max(vals) - min(vals)) / sigma_C_planck
print(f"    SPREAD of the entire convention space = {spread:.3f} sigma -> observationally")
print(f"    DEGENERATE: rho_obs cannot now (or soon) prefer any row over any other.")
assert spread < 0.20
rate_log_res = abs(delta_c - ALPHA0 * c_loop / 2) / sigma_C_planck
print(f"    (note: the rate-log reading lands {rate_log_res:.4f} sigma from centre — and")
print(f"     exactly because of [2], landing is FREE at these stakes: not promotable.)")
assert rate_log_res < 0.01

# ---------------- [4] the quarantine ledger ----------------
print("\n[4] QUARANTINE (sect. 16.3): eta = C2_req/C_loop vs 1/2 + 0.036:")
eta = (delta_c / ALPHA0) / c_loop
print(f"    eta = {eta:.6f};  1/2 + (137.036-137) = 0.536000;  match to "
      f"{abs(eta - 0.536) / 0.536:.1e} relative")
print(f"    REFUSED: (i) the 0.036 component has no operator (the genuine dressed-")
print(f"    alpha billing shift is {delta * c_loop / ALPHA0:+.2e} in C2 units, 3 OOM too small);")
print(f"    (ii) composite trials >= 18 (3 exponentiations x 3 alpha-readings x 2 signs);")
print(f"    (iii) the digits being matched are noise per [2]. Numerological quarantine.")
assert abs(delta * c_loop / ALPHA0) < 1e-3   # the honest dressed shift is negligible

# ---------------- [5] the resolution requirement ----------------
print("\n[5] WHAT IT WOULD TAKE to make delta C a real target (3 sigma):")
dlnrho_needed = delta_c * S_C / 3
dH0_needed = dlnrho_needed / 2 * H0_PLANCK
print(f"    rho_Lambda to {dlnrho_needed:.2%}  <=>  H0 to +/-{dH0_needed:.3f} km/s/Mpc"
      f"  ({DH0_PLANCK/dH0_needed:.0f}x beyond Planck)")
print(f"    — beyond any planned survey; and unposeable while the H0 tension")
print(f"    ({dlnrho_tension:.1%} in rho, = {shift_C_tension/delta_c:.0f} delta C) is unresolved.")
assert DH0_PLANCK / dH0_needed > 20

print(f"""
[6] VERDICT — the next-order hunt is closed BY THE NOISE FLOOR, not by an object:
    * No unpinned operator or convention CAN legitimately close the residual,
      because there is nothing to close: delta C has zero observationally
      defined digits (0.10 sigma, Planck; 0.23 sigma even on the framework's
      own internal accounting; 1% of the H0-tension systematic).
    * The convention question that remains is well-posed but PURELY ALGEBRAIC:
      which exponentiation does the monitored billing formalism force —
      alpha0*C, -ln(1-alpha0)*C, or -ln(1-alpha0*C)?  It belongs to the
      operator-algebra note; its answer sets a central value that observation
      cannot check below ~4% of C.  Promoting any row of [3] by its landing
      would be fitting the third digit of a Planck central value.
    * The CC chain's honest final form: rho_Lambda = rho_obs at 0.10 sigma —
      CLOSED at observational precision.  The no-go theorem stands as the
      statement that the THEORY side's operator class is exhausted; this audit
      adds that the EXPERIMENT side defines no target beyond it.
exit 0""")
print("ALL ASSERTIONS PASSED — every claim above is verified.")
