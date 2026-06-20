#!/usr/bin/env python3
r"""A_req — first free->conditional CONVERSION ATTEMPT under the traffic-operator programme (R1).

Target: is the CC register-handoff prefactor A_req a curvature / eigenvalue / capacity invariant
of its sector's symmetrized service-traffic form -- i.e. can it be moved from FREE to CONDITIONAL
the way w4/w6=2 (Q3 orbit incidence) and 10/27 (Moore capacity) already were?

This is an HONEST attempt, not a forced promotion. It (1) pins what A_req IS, (2) tests the
motivated in-sector traffic-Hessian candidate, (3) applies the §16.3 forking-path discipline and
the per-sector rule, and (4) reports the verdict.

Canon (verified 2026-06-19, DRIFT register-handoff noise ledger l.1702; B4 l.1987):
  * A_req = 1.219224785 multiplies the §5.2 ISOLATED single-fault odds to reach the required
    handoff density p_c (set by matching rho_Lambda -- so A_req currently consumes the observable).
  * p_c = 0.09719075 (required); p_iso = 0.081133 (code-derived isolated Boltzmann occupancy).
  * Motivated in-sector candidate: sqrt(3/2) = (3 walk-active legs)/(2 exiting vacuum legs),
    "already canonical in the gravity coefficient story" -- a fluctuation-DETERMINANT-ratio FORM.
  * §16.3 flag (canon): "source-free small rationals such as 39/32 land equally well or better",
    so the numerical match alone cannot adopt sqrt(3/2).
  * Cross-sector object: the KMS-weighted raw boundary current 1.2219346 (B4) lives on the
    HORIZON Q-graph, NOT the CC crystallisation queue.
Promotion window (canon): ~0.6% in A_req (the factor-two-in-rho tolerance, +-0.578% in p).
"""

A_REQ      = 1.219224785          # canon target (register-handoff noise ledger)
P_C        = 0.09719075          # required handoff density (rho_Lambda-set)
P_ISO      = 0.081133            # isolated §5.2 single-fault occupancy (code-derived)
WINDOW     = 0.006               # ~0.6% promotion tolerance (factor-two in rho)

def odds(p):
    return p / (1.0 - p)

# (1) what A_req IS: the handoff odds-ratio (definitional -> shows WHY it is currently free)
A_req_identity = odds(P_C) / odds(P_ISO)
assert abs(A_req_identity - A_REQ) / A_REQ < 2e-3, \
    "A_req must reproduce as odds(p_c)/odds(p_iso) -- it is a handoff odds-ratio, p_c set by rho_Lambda (so FREE)"

# (2) the candidate field
SQRT_32  = (3.0 / 2.0) ** 0.5     # in-sector traffic-HESSIAN form: sqrt(3 active / 2 vacuum legs)
R_3932   = 39.0 / 32.0            # source-free rational (NO sector motivation)
KMS_BCUR = 1.2219346             # KMS-weighted boundary current -- HORIZON sector (cross-sector)

def dev(x):
    return abs(x - A_REQ) / A_REQ

candidates = {
    "sqrt(3/2)  [in-sector determinant-ratio, MOTIVATED]": (SQRT_32, "in-sector"),
    "39/32      [source-free rational, UNMOTIVATED]":       (R_3932,  "none"),
    "KMS bndy current 1.2219 [HORIZON Q-graph]":            (KMS_BCUR, "cross-sector"),
}

in_window          = {k: v for k, (v, _) in candidates.items() if dev(v) < WINDOW}
in_sector_inwindow = [k for k, (v, sec) in candidates.items() if dev(v) < WINDOW and sec == "in-sector"]
better_unmotivated = [k for k, (v, sec) in candidates.items()
                      if dev(v) < dev(SQRT_32) and sec != "in-sector"]

# (3) the discipline:
assert dev(SQRT_32) < WINDOW, "the motivated determinant-ratio candidate sqrt(3/2) IS within the ~0.6% window"
assert len(better_unmotivated) >= 1, \
    "§16.3 forking: an unmotivated/cross-sector object fits at-least-as-well -> match cannot select sqrt(3/2)"
assert any("HORIZON" in k for k in candidates), "the closest service-current invariant is cross-sector (per-sector rule rejects it)"

# (4) verdict: NOT a clean conversion -- numerical match is blocked; needs the Hessian determinant theorem
converted = (len(in_sector_inwindow) == 1) and (len(better_unmotivated) == 0)   # would need: unique, in-sector, unbeaten
assert converted is False, "A_req is NOT promoted by candidate-matching -- the forking-path forbids it"

bar = "=" * 92
print(bar)
print("A_req — first free->conditional conversion attempt (traffic-operator programme, R1)")
print(bar)
print(f"  A_req (canon)                         = {A_REQ:.9f}")
print(f"  = odds(p_c)/odds(p_iso)               = {A_req_identity:.9f}   (p_c rho_Lambda-set -> FREE)")
print(f"\n  candidate                                              value      dev%     sector")
print(f"  {'-'*86}")
for k, (v, sec) in candidates.items():
    print(f"  {k:52.52s}  {v:.6f}  {dev(v)*100:6.3f}%  {sec}")
print(f"""
{bar}
VERDICT (A_req attempt, exit 0):  SHARPENED, not converted.

  * A_req has the RIGHT KIND of traffic candidate: sqrt(3/2) = sqrt(3 walk-active / 2 vacuum legs)
    -- a fluctuation-DETERMINANT ratio -- at {dev(SQRT_32)*100:.2f}%, INSIDE the ~0.6% promotion window.
    This is genuine progress: a free coefficient now has a motivated, right-shaped in-sector target.
  * But it is NOT promoted. §16.3 forking: an unmotivated source-free rational (39/32, {dev(R_3932)*100:.2f}%)
    fits BETTER, so the numerical match cannot select sqrt(3/2). And the closest service-current
    invariant (KMS boundary current, {dev(KMS_BCUR)*100:.2f}%) is CROSS-SECTOR (horizon Q-graph, not the CC
    crystallisation queue) -- rejected by the per-sector rule.
  * The ONLY legitimate conversion is to COMPUTE the register-handoff service-action Hessian and
    show its fluctuation-determinant ratio = 3/2 exactly (so sqrt(3/2) is FORCED, not shopped).
    Until then A_req stays FREE.

  NET: the first traffic-conversion attempt turns A_req from an arbitrary real into a motivated
  determinant-ratio TARGET sqrt(3/2), conditional on the in-sector Hessian-determinant theorem.
  Progress (right-kind candidate, in tolerance), but not a promotion -- the forking is defeated by
  computing the determinant, not by the numerical hit. Confirms the audit's 'A_req = free'.
{bar}""")
print(f"exit 0 -- A_req sharpened to sqrt(3/2) determinant-ratio target ({dev(SQRT_32)*100:.2f}%); NOT converted "
      f"(§16.3 forking vs 39/32; KMS current cross-sector); Hessian-determinant theorem is the open step.")
