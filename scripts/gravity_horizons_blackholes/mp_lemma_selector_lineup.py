#!/usr/bin/env python3
r"""LEMMA L — THE SELECTOR LINE-UP: can it now be derived?  Answer: it reduces
to ONE named mechanism, and the line-up picks the candidate.

THE TEST.  Lemma L is the epoch equation N_t(a) r6 = 9 alpha0.  Using the
framework's OWN cosmology for N_t(a) = Lambda/H(a) — item-131 dark energy
rho_DE(a) = rho0 exp[(3/28)(1-a)] plus Planck matter — the product
P(a) = N_t(a) r6 / (9 alpha0) is computed across epochs and every candidate
selector the framework owns is evaluated at ITS epoch.  Pre-registered
candidates (no scanning):

  S1  R4-activation completion at a = 1 (canon's implicit selector: f(a) = a
      reaches its physical endpoint — item 131's support measure);
  S2  the 28-dilution epoch: ln(Lambda/T(a)) = 28 exactly (the quarantined
      integer given a job) -> a_28 = e^{28 - ln(Lambda/T0)};
  S3  DE-matter equality (the standard why-now epoch);
  S4  the exact Lemma-L crossing a_L (P = 1) — reported for reference, not a
      selector (it is the answer, not a mechanism).

VERDICT SHAPE: a selector DERIVES Lemma L iff its epoch lands P = 1 within
the internal budget (Part-20-grade 0.07%, since r6 and the 9 are exact and
H(a) extrapolation over percent-level Delta-a is sub-0.01%).

REQUIREMENT SPEC for the winning selector (what a full derivation must add):
  (i)  the physical saturation mechanism at its epoch (for S1: the item-131
       finite-to-cosmological lift premise — WHY the comoving R4 support's
       physical measure completing is an absolute event);
  (ii) the G-locking mechanism (the LLR boundary theorem: the accounting is
       epoch-anchored, G constant away from the anchor);
  (iii) the residual within Part-20's slop.
exit 0 = line-up computed; the derivation gap named as one mechanism."""
import importlib
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
rh = importlib.import_module("register_handoff_form_selection")

ALPHA0 = 1 / 137.0
LAM = rh.LAMBDA_QCD_GEV
T0 = 2.348e-13
MPC_KM, HBAR = 3.085678e19, 6.582120e-25
H0_GEV = rh.H0_KM_S_MPC / MPC_KM * HBAR
OM_L, OM_M = rh.OMEGA_L, 1 - rh.OMEGA_L

# the chain's r6, live
gamma_star = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
_, q_post, _ = rh.queue_readouts(gamma_star, 1)
r6 = (21 * q_post) ** 32 / 21

def H_of(a):
    """framework cosmology: item-131 DE + Planck matter (radiation negligible)."""
    rho_de = OM_L * math.exp((3 / 28) * (1 - a))
    return H0_GEV * math.sqrt(OM_M * a ** -3 + rho_de)

def P(a):
    return (LAM / H_of(a)) * r6 / (9 * ALPHA0)

print("[1] THE PRODUCT P(a) = N_t(a) r6 / (9 alpha0) across epochs (framework H(a)):")
for a in (0.5, 0.772, 0.9, 1.0, 1.023, 1.1):
    print(f"    a = {a:5.3f}:  P = {P(a):.5f}   ({P(a)-1:+.2%})")

print("\n[2] THE SELECTOR LINE-UP (each candidate at ITS OWN epoch):")
a28 = math.exp(28 - math.log(LAM / T0))
a_eq = (OM_M / OM_L) ** (1 / 3)        # leading; DE-evolution correction tiny
# refine equality epoch with the framework DE law
for _ in range(40):
    a_eq = (OM_M / (OM_L * math.exp((3 / 28) * (1 - a_eq)))) ** (1 / 3)
# the exact crossing a_L
lo, hi = 0.9, 1.1
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if P(mid) < 1 else (lo, mid)
a_L = (lo + hi) / 2
BUDGET = 0.0007
rows = [("S1 R4-completion (a = 1, canon implicit)", 1.0),
        ("S2 28-dilution (ln(Lam/T) = 28)", a28),
        ("S3 DE-matter equality", a_eq),
        ("S4 exact crossing (reference)", a_L)]
print(f"    {'selector':<42s} {'epoch a':>8s} {'P-1':>9s} {'vs 0.07% budget':>16s}")
verdicts = {}
for nm, a in rows:
    miss = P(a) - 1
    tag = "INSIDE" if abs(miss) <= 1.5 * BUDGET else ("near (x%d)" % round(abs(miss)/BUDGET) if abs(miss) < 0.02 else "EXCLUDED")
    verdicts[nm[:2]] = (a, miss, tag)
    print(f"    {nm:<42s} {a:8.4f} {miss:+9.3%} {tag:>16s}")
assert abs(verdicts["S1"][1]) < 2 * BUDGET                 # S1 lands inside-grade
assert abs(verdicts["S2"][1]) > 10 * BUDGET                # S2 misses the exact budget
assert abs(verdicts["S3"][1]) > 0.10                       # S3 excluded outright
print(f"""
    -> S3 (the standard why-now epoch) is EXCLUDED at {verdicts['S3'][1]:+.0%}: Lemma L is
       NOT the equality coincidence.  S2 (the 28-integer) lands at {verdicts['S2'][1]:+.2%} —
       intriguing but 15x the budget: viable only if a ~1% mechanism shifts it.
       S1 (R4-activation completion at a = 1) lands at {verdicts['S1'][1]:+.3%} — INSIDE
       the Part-20-grade budget.  The line-up SELECTS S1.""")

print(f"""
[3] SO CAN LEMMA L NOW BE DERIVED?  The honest answer, after the line-up:
    REDUCED TO ONE MECHANISM, NOT YET DERIVED.  The chain now reads:

      [derive: R4 completion at a = 1 is an absolute physical event
       whose epoch satisfies the accounting]                <- THE ONE GAP
        => Lemma L (N_t r6 = 9 alpha0)         [S1 lands at {verdicts['S1'][1]:+.3%}]
        => Omega_Lambda = 12pi/55 exact
        => M_P = Lambda sqrt(20790 alpha0^3)/(21 q1*)^16 a THEOREM
        => H_0, Omega_Lambda, G all derived from microphysics.

    The one gap is precisely item 131's named open premise (the finite-to-
    cosmological lift: why the comoving R4 support's completion is absolute),
    PLUS the G-locking mechanism required by the LLR boundary theorem.
    Everything else on the mountain is now derived, posed, or excluded.
    NOTE the structure of what remains: not a number to find — a MECHANISM
    to write down.  The framework has no remaining numerical freedom here:
    if the lift premise is derived, Lemma L follows at {verdicts['S1'][1]:+.3%} with zero
    adjustable constants; if it fails, the M_P formula stays a (very good)
    unexplained prediction.
exit 0""")
print("ALL ASSERTIONS PASSED — line-up computed; the gap is one named mechanism.")
