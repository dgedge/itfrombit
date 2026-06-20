#!/usr/bin/env python3
r"""A_req register-handoff fluctuation determinant — PRE-REGISTERED theorem attempt (R1).

TARGET THEOREM (pre-registered, no coefficient hunting):
  From the §5.2 handoff service action ALONE, build the quadratic fluctuation Hessian at the
  handoff saddle and show
        det'(H_active) / det'(H_vacuum) = 3/2 ,   hence  A_req = sqrt(3/2).
  (prime = pseudodeterminant: drop gauge/zero modes before comparing.)

ALLOWED inputs : ONLY the §5.2 register-handoff graph + its service/action rules + its forced
                 active/vacuum split.
FORBIDDEN      : cross-sector KMS current, horizon scheduler numerology, post-hoc 39/32, fits.

OUTPUT must be exactly one of:
  CONVERTED       — det' ratio forced to 3/2 from §5.2 alone.
  REJECTED        — ratio != 3/2, or 3/2 only after a branch choice.
  UNDERSPECIFIED  — §5.2 does not define the handoff quadratic action enough to build the Hessian
                    (then: name the missing definition precisely).

This script encodes what §5.2 + the register-handoff ledger ACTUALLY provide (verified 2026-06-19,
ANCHOR §5.2 l.520, l.3896; DRIFT register-handoff l.1702) and tests whether a quadratic Hessian
and its pseudodeterminant ratio are constructible without forbidden inputs.
"""

# ---------------------------------------------------------------------------
# (A) WHAT §5.2 + THE HANDOFF LEDGER ACTUALLY PROVIDE (verified, allowed inputs)
# ---------------------------------------------------------------------------
sec52 = {
    "dynamics_kind"      : "scalar Boltzmann/dwell-time occupancy  M=exp(F/2phi)  (one effective DOF)",
    "strain_recurrence"  : "scalar in q:  level1 q1=(35/2)p^4 ; levels>=2 q_{k+1}=21 q_k^2",
    "isolated_occupancy" : 0.081133,            # p_iso (single-fault Boltzmann occupancy) — a SCALAR
    "arrest_density"     : 0.0971907500,        # p_c (required handoff density)
    "form_selected"      : False,               # DRIFT l.1702 addendum: "the handoff theorem must FIRST
                                                #   select the form (equilibrium §5.2 queue vs slaving)"
    "active_legs"        : 3,                    # §3.2 walk channels {V_em,V_strong,V_weak} = KRAUS LABELS
    "vacuum_legs"        : 2,                    # exiting vacuum/exhaust legs = KRAUS LABELS
    "legs_are"           : "Kraus LABELS in the boot pair ledger (singlet bit x 3x3, l.3896), "
                           "NOT fluctuation MODES with an inter-leg stiffness/coupling matrix",
    "interleg_coupling"  : None,                # §5.2 supplies NO quadratic stiffness matrix on the legs
}

# ---------------------------------------------------------------------------
# (B) ATTEMPT THE CONSTRUCTION, STEP BY STEP — and record where it blocks
# ---------------------------------------------------------------------------
blocks = []   # (step, what is missing)

# Step 1: a fluctuation Hessian needs a MULTI-MODE quadratic action. §5.2's handoff dynamics is
#         a SCALAR occupancy + a SCALAR recurrence -> its Hessian at the fixed point is 1x1.
sec52_quadratic_modes = 0
if sec52_quadratic_modes == 0:
    blocks.append(("multi-mode quadratic action",
                   "§5.2 handoff dynamics is SCALAR (occupancy p, recurrence in q); d^2(action)/d q^2 is "
                   "1x1 — there is no 3-vs-2 mode space to take a det' ratio over"))

# Step 2: the 3 active / 2 vacuum legs are §3.2 Kraus LABELS, not §5.2 fluctuation modes; and §5.2
#         provides NO inter-leg coupling matrix, so no H_active / H_vacuum can be assembled.
if sec52["interleg_coupling"] is None:
    blocks.append(("active/vacuum Hessian matrices",
                   "the 3+2 legs are Kraus labels (counts); §5.2 gives no stiffness/coupling among them, "
                   "so H_active (3x3) and H_vacuum (2x2) cannot be populated from §5.2 alone"))

# Step 3: even the SADDLE is ambiguous — the dynamical FORM (equilibrium queue vs per-tick slaving)
#         is unselected, so there is no unique action to expand around.
if not sec52["form_selected"]:
    blocks.append(("the saddle itself",
                   "the dynamical form (balance §5.2 queue vs slaving) is unselected (DRIFT l.1702), "
                   "so there is no unique handoff action to take a second variation of"))

# Step 4 (control): does the count-only reading FORCE 3/2 without extra structure? No.
#   The only assumption-free reading of '3 active, 2 vacuum modes' is independent modes of a common
#   stiffness k: det'_active/det'_vacuum = k^3/k^2 = k — DIMENSIONFUL, k-dependent, NOT 3/2.
#   So 3/2 does NOT follow from the leg COUNTS; it would require a specific (unsupplied) spectrum.
def count_only_ratio(k):                         # independent modes, common stiffness k
    return (k ** 3) / (k ** 2)                    # = k, not 3/2
count_forces_3_2 = all(abs(count_only_ratio(k) - 1.5) < 1e-9 for k in (0.7, 1.0, 2.0, 3.3))
assert count_forces_3_2 is False, \
    "the leg COUNTS alone do not force 3/2 (independent-mode det' ratio is k, not 3/2) — extra structure needed"

# ---------------------------------------------------------------------------
# (C) VERDICT
# ---------------------------------------------------------------------------
verdict = "UNDERSPECIFIED" if blocks else "REJECTED"   # (CONVERTED only if a forced 3/2 had been built)
assert verdict == "UNDERSPECIFIED", "with these blocks the only honest status is UNDERSPECIFIED"
assert len(blocks) >= 3, "the construction blocks at >=3 distinct missing definitions"

bar = "=" * 94
print(bar)
print("A_req HANDOFF FLUCTUATION DETERMINANT — pre-registered theorem attempt (R1)")
print(bar)
print("  TARGET:  det'(H_active)/det'(H_vacuum) = 3/2  =>  A_req = sqrt(3/2)   (from §5.2 alone)")
print(f"\n  §5.2 provides:")
print(f"    dynamics      : {sec52['dynamics_kind']}")
print(f"    recurrence    : {sec52['strain_recurrence']}")
print(f"    legs          : {sec52['active_legs']} active + {sec52['vacuum_legs']} vacuum — {sec52['legs_are']}")
print(f"    inter-leg form: {sec52['interleg_coupling']}   (no quadratic stiffness matrix)")
print(f"    form selected : {sec52['form_selected']}   (balance vs slaving unresolved)")
print(f"\n  CONSTRUCTION BLOCKS (missing definitions, precisely):")
for i, (step, missing) in enumerate(blocks, 1):
    print(f"    {i}. {step}: {missing}")
print(f"\n  control: leg COUNTS alone give det' ratio = k (k-dependent), NOT 3/2 -> 3/2 is not forced by counts")
print(f"""
{bar}
VERDICT:  {verdict}

  §5.2 does NOT define the handoff quadratic action enough to compute a fluctuation Hessian.
  The determinant theorem is blocked TWO steps upstream of any computation:
    (1) §5.2's handoff dynamics is SCALAR (occupancy + recurrence) — no multi-mode action exists;
    (2) the 3 active / 2 vacuum legs are §3.2 Kraus LABELS, with no §5.2 inter-leg stiffness, so
        H_active / H_vacuum cannot be populated;
    (3) the dynamical FORM (balance vs slaving) is itself unselected — no unique saddle.
  And the count-only reading does not rescue it: independent modes give a ratio of k, not 3/2.

  CONSEQUENCE: A_req cannot be converted free->conditional from §5.2 as it stands. This is NOT a
  rejection of sqrt(3/2) — §5.2 can currently neither force nor refute it. The precise missing
  definition is a MULTI-MODE handoff service action: pick the form, then specify the quadratic
  stiffness/coupling that promotes the 3+2 Kraus labels to fluctuation modes. That object — not a
  determinant computation — is the real open step; the determinant follows once it exists.
  A_req stays FREE; the audit classification is unchanged.
{bar}""")
print(f"exit 0 -- {verdict}: §5.2 handoff action is scalar + legs are labels + form unselected; "
      f"the missing definition is a multi-mode handoff action, named precisely. A_req stays free.")
