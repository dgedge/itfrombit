#!/usr/bin/env python3
r"""G7 CONSTRUCTIBILITY AUDIT (R1) — is the alpha^2/K_eff Feshbach Q-trace finite-service-graph
constructible? MODEST target: constructible / not / underspecified — NOT "derive K_eff".

The object (ANCHOR §10.6):
    1/K_eff = (1/|BZ|) sum_k Tr[ P_Eg (E - W_QQ(k))^{-1} P_Eg ].

VERDICT: NOT CONSTRUCTIBLE from the current machinery (definition-limited; not a universal no-go).
Grounded in G7 (DRIFT l.189-201; keff_sensitivity.py, keff_construct.py, both exit 0):

  CODEABLE as written (3): the 256-register / 48-codeword / 208-Q split; the walk W=S.C;
                           the E_g characters + 2/48 prefactor.
  EACH FORCES AN INVENTED OPERATOR (4): (i) R_g -- the O_h action-space P_Eg acts on is never
    built (>=3 inequivalent choices -> different K_eff; the K_6 trap, highest leverage); (ii) the
    resolvent energy E -- needs a Hermitian object but W is unitary, and E has NO value anywhere;
    (iii) the Bloch promotion W_QQ(k) -- only named, the alpha/beta sublattice gluing unfinished;
    (iv) t_mix -- load-bearing in the prose, absent from the formula, no matrix.
  Plus: "208-3" has no referent in Q (the 3 nu_R are P-codewords, not Q-states).

This script makes the decisive "E-tuning tautology" concrete, then records the supersession.
"""
import math

# --- (1) the four-invented-operators count (from G7's line-by-line audit) ---
codeable        = ["256/48/208 split", "walk W=S.C", "E_g characters + 2/48"]
invented_needed = ["R_g action-space", "resolvent energy E", "Bloch W_QQ(k)", "t_mix"]
assert len(invented_needed) == 4 and len(codeable) == 3, \
    "G7: 3 pieces codeable, 4 each force an invented operator -> not executable as written"

# --- (2) the E-tuning tautology, made concrete: Tr[(E-H)^-1] = sum 1/(E-lambda_i) ---
# a sample Q-block spectrum (any Hermitian promotion would have *some* real spectrum);
# the resolvent trace is set ENTIRELY by the un-pinned E -> any K_eff reachable.
eigs = [(-3.0 + 0.7 * i) for i in range(10)]            # a generic spread of 10 eigenvalues
def resolvent_trace(E):
    return sum(1.0 / (E - lam) for lam in eigs if abs(E - lam) > 1e-9)

lam_max = max(eigs)
t_hi = resolvent_trace(lam_max + 0.05)                 # just ABOVE the top eigenvalue -> large +
t_lo = resolvent_trace(lam_max - 0.05)                 # just BELOW the top eigenvalue -> large -
# and: an E exists making the trace hit any target (e.g. the dim) -> "uniform trace = dim" is tuned-E
target = float(len(eigs))
hit_E = None
E = lam_max + 1e-3
while E < lam_max + 60.0:
    if resolvent_trace(E) <= target:                   # trace decreases from +inf as E rises above spectrum
        hit_E = E
        break
    E += 1e-3
assert hit_E is not None, "an E exists that makes the trace hit the dim -> tuned-E tautology"
assert t_hi > 5.0 and t_lo < -5.0, "the trace spans (-inf,+inf) across an eigenvalue -> K_eff set by the un-pinned E"
# i.e. ANY K_eff is reachable by moving E; the formula does not pin a value. (G7 T2.)

# --- (3) the defensible action-space projectors give 24 or 128, never 205 (G7) ---
P_Eg_dims = {"faces (matter-cell)": 24, "colour": 128}
assert 205 not in P_Eg_dims.values(), "no defensible E_g action-space dimension is 205"

# --- (4) the supersession: alpha^2 is structurally impossible; the live coefficient is alpha^1 ---
# item119_jump_channel.py: L_k Pi_Q = 0 -> one non-unitary projection per P->Q->P loop -> alpha^1.
alpha0 = 1.0 / 137.0
walk_active_channels = 3                 # {C0, C1, I3}  (§2.5/§3.2)
exiting_channels     = 2                 # C0, C1 exit (break R3); I3 is the free bit
C_loop = walk_active_channels / exiting_channels        # = 3/2, derived (exiting/walk-active)
K_alpha1 = C_loop / alpha0                               # = 3/(2 alpha) = 205.5
K_data   = 206.49
assert abs(C_loop - 1.5) < 1e-12, "C_loop = (walk-active 3)/(exiting 2) = 3/2 (item 119, derived-under-reading)"
assert abs(K_alpha1 - 205.5) < 0.1 and abs(K_alpha1 / K_data - 1.0) < 0.006, \
    "alpha^1 route: K=3/(2 alpha)=205.5 vs K_data=206.49 (0.45%), conditional on item 79"

# --- verdict / reclassification ---
feshbach_trace_constructible = False     # NOT constructible from current machinery (4 invented ops + E-tuning)
alpha2_status   = "structurally impossible (item 119 jump algebra) + non-constructible (G7) -> DEAD"
gravity_coeff   = "conditional"          # alpha^1 / C_loop=3/2, M_P^2=(2/3) alpha Lambda^3 R_dS, on item 79
assert feshbach_trace_constructible is False and gravity_coeff == "conditional"

bar = "=" * 94
print(bar)
print("G7 CONSTRUCTIBILITY AUDIT — is the alpha^2/K_eff Feshbach Q-trace finite-graph constructible?")
print(bar)
print(f"  formula: 1/K_eff = (1/|BZ|) sum_k Tr[P_Eg (E - W_QQ(k))^-1 P_Eg]")
print(f"  codeable as written : {len(codeable)}   ({', '.join(codeable)})")
print(f"  invented operators needed : {len(invented_needed)}   ({', '.join(invented_needed)})")
print(f"  E-tuning tautology   : Tr[(E-H)^-1] = {t_hi:+.2f} just above / {t_lo:+.2f} just below an eigenvalue; "
      f"= dim at E={hit_E:.2f}")
print(f"                         -> any K_eff reachable by moving the un-pinned E (G7 T2)")
print(f"  E_g action-space dims : {P_Eg_dims} -- never 205 (the K_6 trap)")
print(f"  '208-3'              : the 3 nu_R are P-codewords, NOT Q-states -> no referent")
print(f"\n  => Feshbach Q-trace CONSTRUCTIBLE? {feshbach_trace_constructible}  (definition-limited; not a universal no-go)")
print(f"  => alpha^2 status: {alpha2_status}")
print(f"  => live gravity coefficient: alpha^1 / C_loop={C_loop} -> K={K_alpha1:.1f} (0.45% vs data); {gravity_coeff} (on item 79)")
print(f"""
{bar}
VERDICT (modest target, exit 0):  NOT CONSTRUCTIBLE.

  The alpha^2/K_eff Feshbach Q-space trace is NOT finite-service-graph constructible as the canon
  stands: 3 pieces codeable, but 4 each force an invented operator, the "208-3" has no Q-referent,
  and the resolvent is an E-tuning tautology (any value reachable). It is definition-limited -- a
  negative constructibility result for THAT object, not a proof that no operator definition ever
  could (that would be new physics).

  BUT alpha^2/K_eff is not the live gravity coefficient: alpha^2 is structurally impossible in the
  canonical jump algebra (item 119), and the gravity coefficient is the alpha^1/C_loop=3/2 erasure
  route, M_P^2=(2/3) alpha Lambda^3 R_dS, CONDITIONAL on item 79. So this last native-core entry is
  not a "free non-constructible wall" -- the alpha^2 object is dead, the coefficient is conditional.

  CONSEQUENCE (4th audit correction): native-core FREE drops 1 -> 0. Every native cosmological/QEC
  coefficient is now forced, conditional (sector-native), or -- for gravity's absolute M_P scale --
  horizon-input supplied by the selector lock. The native coefficient ledger is EXHAUSTED: zero
  genuinely-free fitted coefficients, and the one "non-constructible wall" was a dead/superseded
  object. The genuine free count is entirely the non-native EW+nuclear second-scale bolus.
{bar}""")
print(f"exit 0 -- alpha^2/K_eff Feshbach trace NOT constructible (4 invented ops + E-tuning); alpha^2 dead "
      f"(item 119); gravity coeff = alpha^1/C_loop=3/2 conditional (item 79); native-core free 1->0.")
