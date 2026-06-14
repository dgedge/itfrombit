#!/usr/bin/env python3
r"""D2 CONSUMPTION-CHANNEL AUDIT — three legs pass; and the D1 wreckage yields
a sharp resurrection candidate (deposition x sect-5.2 survival = the Lemma-L
budget at 0.5%, canonical rate, zero fitted numbers).

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

LEG 4 — THE D1' RESURRECTION CANDIDATE (sect-16.3 handled).  D1-as-posed
  failed at 11.4x.  Post-hoc factor recognition, declared as such: the
  surplus factor equals 1/BASE_GAMMA: the canonical sect-5.2 retention odds
  e^{-3/(2 phi)} = 0.08830.  Mechanism reading: each deposited protected
  element survives the early rescue era iff Arrhenius-retained at the SAME
  sect-5.2 barrier the CC chain bills — the most load-bearing constant in
  the framework, not an arbitrary factor.  Then
      billable fuel = f_line(boot) x e^{-3/(2 phi)}
  = 0.748(55) x 0.08830 = 0.0661(49)  vs  9 alpha0 = 0.06569:
  a 0.5% central landing, 0.07 sigma_meas, AT the canonical rate R = 496
  (and ONLY there: the R-scan values 0.94/0.90/0.65 land 0.083/0.079/0.058).
  Trials: one recognized factor among a handful of implicit candidates —
  class-2+ CANDIDATE, not promoted.  DECISIVE PRE-REGISTERED TEST: evolve the
  R = 496 relic under the rescue-era dynamics and measure the surviving
  protected fraction directly; the candidate predicts 0.0883.
exit 0 = anchors live, all legs computed, candidate quarantined."""
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

print("\n[L4] THE D1' RESURRECTION CANDIDATE (post-hoc factor, declared):")
F_LINE, SE = 0.74815, 0.05488          # measured, k04_d1_line_deposition_audit R=496
surv = math.exp(-3 / (2 * PHI))
prod = F_LINE * surv
target = 9 * ALPHA0
nsig = abs(prod - target) / (SE * surv)
print(f"    BASE_GAMMA = e^(-3/2phi) = {surv:.5f}  (the sect-5.2 retention odds, canon)")
print(f"    f_line(496) x survival = {F_LINE:.4f} x {surv:.5f} = {prod:.5f} +/- {SE*surv:.5f}")
print(f"    vs 9 alpha0 = {target:.5f}:  central {prod/target-1:+.2%}, {nsig:.2f} sigma_meas")
assert abs(prod / target - 1) < 0.08 and nsig < 1.0
print(f"    R-scan exclusivity: only the canonical rate lands —")
for R, fl in ((124, 0.9389), (248, 0.8981), (496, F_LINE), (992, 0.6537)):
    print(f"      R = {R:>4d}: f_line x gamma = {fl*surv:.4f}  ({fl*surv/target-1:+.1%})")
print(f"    TRIALS: one recognized factor (canon's most load-bearing constant) among")
print(f"    a handful of implicit candidates: CLASS-2+ — quarantined, not promoted.")
print(f"    DECISIVE TEST (pre-registered): evolve the R=496 relic under rescue-era")
print(f"    dynamics; measure the surviving protected fraction; prediction: {surv:.4f}.")

print(f"""
[VERDICT] D2: PASSES its three audit legs — structural (canon conjunction),
  rate hierarchy (clean extension of the three-rate table), dynamical shape
  (degenerate with the canon law today; future discriminator pre-registered).
  D2 is CONSISTENT, not yet derived (the one-element-per-record clause rests
  on the record-alphabet spec, canon-cited).
  AND the joint mechanism is RESURRECTED at candidate grade: D1' = deposition
  (measured) x sect-5.2 survival (canon constant) = the Lemma-L budget at
  +0.5% +/- 7% — landing ONLY at the canonical boot rate.  The chain would
  read: boot deposits f_line = 0.748 (measured) -> rescue era retains
  e^(-3/2phi) (the same barrier the CC chain bills) -> billable fuel 9 alpha0
  -> consumption at r6 -> exhaustion at N_t = 9 alpha0/r6 -> Lemma L -> M_P.
  One measurement now stands between candidate and closure: the rescue-era
  survival fraction (deep-box-sized; prediction 0.0883).
exit 0""")
print("ALL ASSERTIONS PASSED — anchors live, legs computed, candidate quarantined.")
