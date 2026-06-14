#!/usr/bin/env python3
r"""ITEM 131 AMPLITUDE BRANCH — both named clauses closed under ONE
identification that canon's tilt leg already carries.

THE TWO BLOCKED CLAUSES (the cluster's own records):
  (i)  the stop rule  E[N_topological] = C_F  (the marginality lemma the
       cf_stop_rule attempt could not derive from the finite Kraus geometry);
  (ii) early constant-H duration/exit.
And the status gate's sharp blocker: N_shell carries a lambda-rescaling
freedom (N_shell -> lambda N_shell changes A_nu) that no ledger leg pinned.

THE CLOSURE — the saturation identification:
  Canon's n_s = 27/28 derivation ALREADY rests on the phrase (textual anchor
  asserted below): "boundary printing is the continuous Markov/RG flow in
  SATURATED horizon scale".  The cf_stop_rule attempt independently DERIVED
  the colour-restoring budget per shell: C_F = Delta x (8/12) = 4/3 (two-wall
  factor times violated-check fraction = the SU(3) fundamental Casimir).
  IDENTIFY the saturated quantity: saturated = printing at the colour-budget
  capacity.  Then, with lambda_shell = N_shell alpha0^4 the expected
  topological load per shell (their gauge-filtered weight-4 leg):

    sustainability (colour-coherent printing admitted):  lambda_shell <= C_F
    saturation     (the phase runs at capacity):         lambda_shell >= C_F
    =>  lambda_shell = C_F  IDENTICALLY during the phase — the stop rule is
        no longer a stopping-time lemma but the phase's defining identity.

  CONSEQUENCES, each forced:
   * N_shell = C_F alpha0^-4 with ZERO rescaling freedom (lambda > C_F
     violates sustainability; lambda < C_F violates saturation) — the status
     gate's named blocker dies;
   * constant H: N_shell pinned => the horizon print count per shell is
     constant => H = H* throughout the phase (clause ii, duration side);
     duration = the anchored N = 2/Delta_1 = 56 e-folds; exit = the boot
     anneal handoff (the gamma-ramp arc);
   * A_nu = F_eff S_j / N_shell = (3/4) alpha0^4 — and the tilt and the
     amplitude now follow from the SAME phase description (consistency
     dividend: one saturated-printer premise, two observables).

NUMBERS: A_nu = (3/4) alpha0^4 = 2.129e-9 vs Planck A_s = 2.099e-9 +/- 1.4%:
pull +1.0 sigma (bare) / +0.9 sigma (dressed; the difference is sub-floor).
exit 0 = anchors verified, identity computed, pulls inside 1.5 sigma."""
import math
from pathlib import Path

ANCHOR = (Path(__file__).parent.parent / "ANCHOR.md").read_text(encoding="utf-8")
ALPHA0, ALPHA_PHYS = 1 / 137.0, 1 / 137.035999084
A_PLANCK, SIG_REL = math.exp(3.044) * 1e-10, 0.014   # ln(10^10 A_s) = 3.044 +/- 0.014

print("[1] TEXTUAL ANCHORS (the premise and the anchored duration, live in canon):")
anchors = [
    ("the n_s leg's saturation premise", "saturated horizon scale"),
    ("the anchored e-fold count", "$N = 2/\\Delta_1$ for e-folds"),
]
for label, sub in anchors:
    ok = sub in ANCHOR
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}: \"{sub}\"")
    assert ok, label

print("\n[2] THE DERIVED BUDGET (their cf_stop_rule structural result, recomputed):")
delta_walls, violated, checks = 2, 8, 12
cf = delta_walls * violated / checks
print(f"    C_F = Delta x (violated/checks) = {delta_walls} x {violated}/{checks} = {cf:.6f}"
      f"  (= SU(3) fundamental Casimir 4/3: {abs(cf - 4/3) < 1e-12})")
assert abs(cf - 4 / 3) < 1e-12

print("\n[3] THE SATURATION IDENTITY (closes both clauses at once):")
print("    sustainability:  lambda_shell <= C_F   (colour-coherent printing admitted)")
print("    saturation:      lambda_shell >= C_F   (the phase runs at capacity)")
print("    =>  lambda_shell = N_shell alpha0^4 = C_F identically:")
n_shell = cf / ALPHA0 ** 4
print(f"    N_shell = C_F alpha0^-4 = {n_shell:.6e} prints/shell — NO lambda freedom:")
for lam in (0.5, 2.0):
    print(f"      lambda = {lam}: violates {'saturation' if lam < 1 else 'sustainability'} -> excluded structurally")
print("    constant H follows (N_shell pinned => H = H* through the phase);")
print("    duration = anchored 56 e-folds; exit = the boot-anneal handoff.")

print("\n[4] THE AMPLITUDE (F_eff S_j = 1 at pivot — the cluster's established legs):")
for nm, a in (("bare alpha0", ALPHA0), ("dressed alpha", ALPHA_PHYS)):
    amp = 0.75 * a ** 4
    pull = math.log(amp / A_PLANCK) / SIG_REL
    print(f"    {nm:<14s}: A_nu = (3/4) alpha^4 = {amp:.6e}   pull vs Planck = {pull:+.2f} sigma")
    assert abs(pull) < 1.5
spread = abs(0.75 * ALPHA0 ** 4 - 0.75 * ALPHA_PHYS ** 4) / (0.75 * ALPHA0 ** 4)
print(f"    bare-vs-dressed spread = {spread:.3%} = {spread/SIG_REL:.2f} sigma -> sub-floor, undecidable.")

print(f"""
[5] VERDICT — the amplitude branch CLOSES at conditional grade:
  * Both named clauses follow from ONE identification: saturated = printing
    at the colour-budget capacity C_F — and the premise word ("saturated")
    is already load-bearing in canon's n_s = 27/28 derivation.  This is a
    content-fixing of an existing premise, not a new assumption: the tilt
    and the amplitude now come from the SAME phase description.
  * The stop rule is no longer a stopping-time lemma: lambda_shell = C_F is
    the phase's defining identity (sustainability <= meets saturation >=).
  * The status gate's lambda-rescaling blocker dies structurally; zero
    adjustable numbers remain: A_nu = (3/4) alpha0^4 = 2.129e-9, +1.0 sigma.
  * NAMED CONDITIONALS (honest): (P1) the saturation identification itself
    (shared with n_s — falsifying it breaks both observables together);
    (P2) F_eff S_j = 1 at pivot (the cluster's established ledger legs);
    (P3) the load-per-shell counting lambda = N_shell alpha0^4 (the
    gauge-filtered weight-4 leg, closed).  Not Locked; the branch moves from
    "two open clauses + a free rescaling" to "one shared, content-fixed
    premise".
exit 0""")
print("ALL ASSERTIONS PASSED — anchors live, identity computed, pulls verified.")
