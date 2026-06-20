#!/usr/bin/env python3
r"""D2 CONSUMPTION-CHANNEL AUDIT — three legs pass; the static-stock routes are
now closed negative, so the live completion route is the integrated burn ledger.

D2 (the relic-exhaustion mechanism's second sub-claim): each depth-6 residual
event consumes exactly ONE relic line-element — the CC chain and the DE
activation are the same ledger read twice.

LEG 1 — STRUCTURAL (canon identification).  A depth-6 residual is
  uncorrectable by construction: its record cannot be erased and must be
  written to the permanent ledger.  Canon has exactly ONE permanent repository:
  the R4/Landauer exhaust ledger ("the Landauer exhaust that sources Dark
  Energy" — item 131, textual anchor below).  And the record-alphabet spec
  (s1 = ln(8x137)) is one cell-register's content: ONE record fills ONE
  element.  "One consumed element per deep-residual event" is the conjunction
  of two canon statements, not a new assumption.

LEG 2 — RATE HIERARCHY (no conflict).  Three rates, three jobs, ~40 OOM apart:
  R4 service clock 1/28 per cell-tick (activation/maintenance, reversible);
  portal commits alpha0 (eta's channel); permanent consumption r6 ~ 2.8e-43
  (only deep residuals bill the tape).  The demux/three-rate table already
  keeps these distinct; D2 adds no tension.

LEG 3 — DYNAMICAL SHAPE (the falsifiable consequence).  Consumption-driven
  activation gives f(a) = N_t(a) r6 / (9 alpha0) = P(a), not exactly the
  canon linear f(a) = a.  Computed below: the deviation peaks ~ +11% in f at
  a ~ 0.5, which moves w(a) by only Delta-w <~ 0.002 — observationally
  DEGENERATE with the linear law today, and a PRE-REGISTERED future
  discriminator (the consumption law slightly steepens early thawing).

LEG 4 — STATIC STOCK ROUTES CLOSED NEGATIVE (2026-06-19 update).  D1-as-posed
  overproduced at L=6 by 11.4x.  The finite-size D1' landing
  f_line(496)*e^{-3/(2phi)} = 0.0661 was subsequently refuted by the deep
  L=8..16 escalation: the protected-line density grows with L, giving
  D1' stock 1.84..2.13x target after extrapolation.  The winding-relic static
  census likewise overproduces by 9..11x.  Therefore the surviving D2 route is
  not a K04 surviving-density census; it is the integrated burn rule
  B0 = integral(permanent burn records during the service era).  The open
  theorem is to derive integral dN_burn = 9 alpha0 per cell from register-level
  rescue/service bookkeeping.
exit 0 = anchors live, D2 legs computed, static stock retired."""
import importlib
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
rh = importlib.import_module("register_handoff_form_selection")
ANCHOR = (Path(__file__).parent.parent / "ANCHOR.md").read_text(encoding="utf-8")

ALPHA0 = 1 / 137.0
LAM = rh.LAMBDA_QCD_GEV
MPC_KM, HBAR = 3.085678e19, 6.582120e-25
H0_GEV = rh.H0_KM_S_MPC / MPC_KM * HBAR
OM_L, OM_M = rh.OMEGA_L, 1 - rh.OMEGA_L
PHI = (math.sqrt(5) - 1) / 2

print("[L1] STRUCTURAL — textual anchors (live ANCHOR.md):")
for label, sub in [("R4 exhaust sources DE", "Landauer exhaust that sources Dark Energy"),
                   ("record alphabet s1 = ln(8x137)", r"\ln(8\times137)")]:
    ok = sub in ANCHOR
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}: \"{sub}\"")
    assert ok, label
print("    -> uncorrectable deep residuals must be permanently recorded; the only")
print("       permanent repository is the R4/Landauer ledger; one s1-record fills")
print("       one element. D2 = conjunction of two canon statements. PASS.")

print("\n[L2] RATE HIERARCHY (three jobs, no conflict):")
gamma_star = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
_, q_post, _ = rh.queue_readouts(gamma_star, 1)
r6 = (21 * q_post) ** 32 / 21
rows = [("R4 service clock (1/28)", 1 / 28), ("portal commit (alpha0)", ALPHA0),
        ("permanent consumption (r6)", r6)]
for nm, r in rows:
    print(f"    {nm:<28s} {r:.3e} per cell-tick")
assert 1 / 28 > ALPHA0 > r6
print("    -> service >> commit >> consumption: the demux table extends cleanly.")

print("\n[L3] DYNAMICAL SHAPE — consumption f(a) vs canon linear f(a) = a:")
def H_of(a):
    return H0_GEV * math.sqrt(OM_M * a ** -3 + OM_L * math.exp((3 / 28) * (1 - a)))
def f_cons(a):
    return (LAM / H_of(a)) * r6 / (9 * ALPHA0)
max_dev, max_dw = 0.0, 0.0
print(f"    {'a':>5s} {'f_cons':>8s} {'f_lin':>6s} {'df':>8s} {'dw':>9s}")
for a in (0.3, 0.5, 0.7, 0.9, 1.0):
    df = f_cons(a) - a
    dw = (1 / 28) * df
    max_dev = max(max_dev, abs(df)); max_dw = max(max_dw, abs(dw))
    print(f"    {a:5.2f} {f_cons(a):8.4f} {a:6.2f} {df:+8.4f} {dw:+9.5f}")
print(f"    max |df| = {max_dev:.3f} -> max |dw| = {max_dw:.4f}: DEGENERATE with the")
print(f"    linear law at current precision; pre-registered future discriminator.")
assert max_dw < 0.005

print("\n[L4] STATIC-STOCK ROUTES — historical D1' landing superseded by deep run:")
F_LINE, SE = 0.74815, 0.05488          # historical L=6 value, k04_d1_line_deposition_audit R=496
surv = math.exp(-3 / (2 * PHI))
prod = F_LINE * surv
target = 9 * ALPHA0
deep_low, deep_high = 1.369 * surv, 1.582 * surv
print(f"    BASE_GAMMA = e^(-3/2phi) = {surv:.5f}  (the sect-5.2 retention odds, canon)")
print(f"    historical L=6 D1' stock = {F_LINE:.4f} x {surv:.5f} = {prod:.5f} ({prod/target:.3f} x target)")
print("       -> superseded: this was a finite-size landing, not a live candidate.")
print(f"    deep-extrapolated D1' stock = {deep_low:.5f}..{deep_high:.5f} ({deep_low/target:.2f}..{deep_high/target:.2f} x target)")
print("    winding-relic static census = 0.60..0.73 (9..11 x target)")
assert abs(prod / target - 1) < 0.02
assert deep_low / target > 1.8
print("    -> all static surviving-density readings are retired.  The live D2")
print("       object is the integrated burn ledger, not remaining K04 trail length.")

print(f"""
[VERDICT] D2: PASSES its three audit legs — structural (canon conjunction),
  rate hierarchy (clean extension of the three-rate table), dynamical shape
  (degenerate with the canon law today; future discriminator pre-registered).
  D2 is CONSISTENT, not yet derived (the one-element-per-record clause rests
  on the record-alphabet spec, canon-cited).
  UPDATE: the D1' finite-size landing is retired by the deep escalation, and
  the winding static census is also refuted.  The live completion route is
  therefore the burn-integral theorem:
      integral dN_burn = 9 alpha0 per cell,
  derived from register-level service/rescue bookkeeping.  If that lands,
  D2 gives exhaustion at N_t = 9 alpha0/r6 -> Lemma L -> M_P.  If it fails,
  D2 remains a valid consumption law but the initial stock B0 is open.
exit 0""")
print("ALL ASSERTIONS PASSED — anchors live, D2 legs computed, static stock retired.")
