#!/usr/bin/env python3
r"""R1 -- THE RECORD-ACTION / MAXIMUM-CALIBER PRINCIPLE (proposition): the two forced legs.

Proposed missing definition layer. The substrate "code + walk + ledgers" lacks a *selection
rule over histories*. Add the least-biased (maximum-caliber) measure on monitored QEC service
histories, with the action read off the history's own record:

      P[gamma] ~ exp(-A_rec[gamma]),
      A_rec[gamma] = A_traffic[gamma]  -  (1/2) dS_rec[gamma]  +  i*Phi_rec[gamma]
                     \_______ ______/     \______ ______/        \____ ____/
                      time-SYMMETRIC        time-ANTISYM           COMPLEX
                      (dynamical activity)  (Crooks entropy)       (recovery holonomy)

Slogan: the record IS the action. One new object (the exponential glue); each term is a
pre-existing canon object.

This script verifies ONLY the two legs that are forced by existing canon + standard stochastic
thermodynamics, and explicitly does NOT touch the three open legs.

  FORCED / verified here
  ----------------------
  (1) RECORD-ACTION UNIT. The per-event action quantum is the already-derived record entropy
      s1 = ln(8 x 137) = 6.999 k_B  (record_alphabet_derivation.py, exit 0; address alphabet
      8 = record_content_from_syndrome.py Lemma 1; channel 137 = item 79). The principle's
      only new claim: exp(-s1) is the per-event WEIGHT -- the same record that measures entropy
      (item 144: dS/dt = alpha*log2(8) = 3*alpha) also weights probability. log2(8)=3 is the
      same 3. (We re-confirm the arithmetic + the item-144 tie; the derivation lives in the
      cited scripts, not re-done here.)

  (2) CROOKS SPLIT + KMS EQUILIBRIUM LIMIT on a minimal monitored CTMC (a 3-state ring -- the
      smallest with a cycle, so a genuine arrow is possible).
      (2a) entropy fixes ONLY the antisymmetric part: ln(P[g]/P[g~]) = dS_rec[g], and the
           time-symmetric TRAFFIC (survival/dwell factors) CANCELS in the ratio -> dwell times
           do not enter dS at all. (Crooks fixes ratios of forward/reverse, NOT absolute rates;
           this is exactly why A_traffic must be a separate term -- the user's correction.)
      (2b) the equilibrium limit IS the canonical KMS scheduler (DRIFT B1): with B1's
           half-Boltzmann rule  W_ij = A_ij exp[-beta(F_j-F_i)/2]  and zero drive, the
           stationary law is Boltzmann pi_i ~ exp(-beta F_i) and detailed balance holds; a
           nonconservative drive breaks detailed balance and forces dS_ss > 0 (the arrow).
      (2c) NUMEROLOGY LOCUS: rescaling every rate by a common clock c leaves pi, the per-jump
           entropy, and every Crooks ratio invariant -- only wall-clock time rescales. The
           absolute scale is ONE multiplier (the traffic clock); the entropy/KMS sector is
           scale-free. (Why the parameter ledger can be favourable: the free multipliers live
           in A_traffic, not in the entropy or -- if the connection is canonical -- the phase.)

  NOT verified here (status now tracked by r13_record_action_measure_status.py)
  ---------------------------------------------------------------------------
   - Phi_rec (recovery holonomy): CP phases / baryon sign / Majorana phase  -- still open.
   - the traffic CLOCK pin beyond the native selector lock -- still partial.
   - Born = closed record-action pair -- closed separately by born_closed_record_pair.py
     plus the Gleason/non-contextuality route; not re-derived in this forced-leg script.
   - GUARDRAIL (DRIFT E5 + item 131 T2): this principle does NOT license inserting s1 or any
     KMS heat factor into an AMPLITUDE (A_nu) -- that stays an illegal convention-dependent
     coefficient. s1 is a per-EVENT entropy/photon-bookkeeping unit only.
"""
import math
from itertools import product

# ============================================================================
# LEG (1) -- the record-action unit (cite + re-confirm arithmetic and item-144 tie)
# ============================================================================
ADDRESS = 8          # record_content_from_syndrome.py Lemma 1 (vertex-label alphabet, DERIVED)
CHANNEL = 137        # item 79 commit channel (136 ideal-code + 1); the integer canon reading
s1 = math.log(ADDRESS * CHANNEL)            # = ln(8 x 137): the per-event record-action quantum
bits_per_cell = math.log2(ADDRESS)          # = 3: the SAME 3 as item 144's dS/dt = 3*alpha
w_event = math.exp(-s1)                      # the principle's new reading: per-event WEIGHT

assert abs(s1 - math.log(1096)) < 1e-12
assert abs(s1 - 6.9994) < 1e-3,                "s1 = ln(8x137) must be 6.999 k_B (canon, -0.1 sigma)"
assert abs(bits_per_cell - 3.0) < 1e-12,       "log2(8)=3 ties the record content to item 144 (3 alpha)"
assert abs(s1 - (math.log(ADDRESS) + math.log(CHANNEL))) < 1e-12, \
    "action additive: record-content part ln8 + coupling part ln137"
assert abs(w_event - 1.0 / 1096.0) < 1e-12,    "per-event weight is exp(-s1) = 1/(8x137)"

# ============================================================================
# LEG (2) -- minimal monitored CTMC: a 3-state ring 0->1->2->0
#   symmetric adjacency A_ij = A_ji  (orientation-blind one-bit service = the TRAFFIC),
#   free energies F, inverse temperature beta, nonconservative drive f per forward edge.
# ============================================================================
RING = [(0, 1), (1, 2), (2, 0)]                       # the oriented cycle edges
EDGES = RING + [(j, i) for (i, j) in RING]
FWD = set(RING)

def rates(F, beta, A, f):
    """B1 half-Boltzmann rule + optional nonconservative drive f (DRIFT B1)."""
    k = {}
    for (i, j) in EDGES:
        drive = (f / 2.0) if (i, j) in FWD else (-f / 2.0)
        k[(i, j)] = A[frozenset((i, j))] * math.exp(-beta * (F[j] - F[i]) / 2.0 + drive)
    return k

def escape(k):
    lam = {s: 0.0 for s in (0, 1, 2)}
    for (i, j), v in k.items():
        lam[i] += v
    return lam

def stationary(k):
    """Matrix-tree theorem on 3 nodes: pi_r ~ sum over spanning trees rooted at r of edge-products."""
    nodes = (0, 1, 2)
    def rooted(r):
        others = [n for n in nodes if n != r]
        choices = [[(n, m) for m in nodes if m != n] for n in others]
        tot = 0.0
        for combo in product(*choices):
            succ = {e[0]: e[1] for e in combo}
            ok = True
            for n in others:                          # every node must reach r acyclically
                seen, cur = set(), n
                while cur != r:
                    if cur in seen or cur not in succ:
                        ok = False
                        break
                    seen.add(cur)
                    cur = succ[cur]
                if not ok:
                    break
            if ok:
                p = 1.0
                for e in combo:
                    p *= k[e]
                tot += p
        return tot
    w = [rooted(r) for r in nodes]
    Z = sum(w)
    return [x / Z for x in w]

def path_weight(states, dwell, k, lam, pi0):
    """P[gamma] = pi(x0) * prod_jumps k(x_i->x_{i+1}) * prod_sojourns exp(-lam[x] tau)."""
    w = pi0[states[0]]
    for idx in range(len(states) - 1):
        w *= k[(states[idx], states[idx + 1])]
    surv = 1.0
    for s, tau in zip(states, dwell):
        surv *= math.exp(-lam[s] * tau)
    return w * surv, surv

def entropy_production(states, k, pi):
    dS_sys = math.log(pi[states[0]]) - math.log(pi[states[-1]])
    dS_med = sum(math.log(k[(states[t], states[t + 1])] / k[(states[t + 1], states[t])])
                 for t in range(len(states) - 1))
    return dS_sys + dS_med

# ---- (2b) KMS EQUILIBRIUM LIMIT: zero drive -> Boltzmann stationarity + detailed balance ----
F = {0: 0.0, 1: 1.3, 2: 0.4}
beta = 0.85
A = {frozenset(e): a for e, a in {(0, 1): 1.0, (1, 2): 2.0, (2, 0): 1.5}.items()}
k_eq = rates(F, beta, A, f=0.0)
pi_eq = stationary(k_eq)
lam_eq = escape(k_eq)
Z = sum(math.exp(-beta * F[s]) for s in (0, 1, 2))
for s in (0, 1, 2):
    assert abs(pi_eq[s] - math.exp(-beta * F[s]) / Z) < 1e-12, \
        "zero-drive stationary must be the KMS Boltzmann law pi_i ~ exp(-beta F_i) (DRIFT B1)"
for (i, j) in RING:
    assert abs(pi_eq[i] * k_eq[(i, j)] - pi_eq[j] * k_eq[(j, i)]) < 1e-12, \
        "KMS detailed balance pi_i W_ij = pi_j W_ji must hold at f=0"

# ---- (2a) CROOKS SPLIT: ln(P/P~) = dS_rec, traffic cancels, dwell-independent ----
k_dr = rates(F, beta, A, f=0.9)                       # driven: a genuine arrow
lam_dr = escape(k_dr)
pi_dr = stationary(k_dr)
trajs = [([0, 1, 2, 0], [0.7, 0.3, 1.1, 0.5]),
         ([1, 2, 0, 1, 2], [0.2, 0.9, 0.4, 0.6, 0.3]),
         ([2, 1, 0], [1.0, 0.5, 0.8])]
for st, dw in trajs:
    Pf, surv_f = path_weight(st, dw, k_dr, lam_dr, pi_dr)
    Pr, surv_r = path_weight(st[::-1], dw[::-1], k_dr, lam_dr, pi_dr)
    assert abs(surv_f - surv_r) < 1e-12, \
        "TRAFFIC (survival/dwell) is time-symmetric -> cancels in the Crooks ratio"
    lhs = math.log(Pf / Pr)
    rhs = entropy_production(st, k_dr, pi_dr)
    assert abs(lhs - rhs) < 1e-9, "ln(P/P~) must equal dS_rec (entropy = the antisymmetric part only)"
    # dwell-independence: rescale all dwell times by 3.7 -> dS unchanged (timing is pure traffic)
    Pf2, _ = path_weight(st, [3.7 * x for x in dw], k_dr, lam_dr, pi_dr)
    Pr2, _ = path_weight(st[::-1], [3.7 * x for x in dw[::-1]], k_dr, lam_dr, pi_dr)
    assert abs(math.log(Pf2 / Pr2) - rhs) < 1e-9, "dS_rec must be independent of dwell times"

# equilibrium reversibility: a CLOSED loop at f=0 has zero entropy production (KMS)
assert abs(entropy_production([0, 1, 2, 0], k_eq, pi_eq)) < 1e-12, \
    "closed-loop entropy production must vanish in the KMS (f=0) limit"

# NESS arrow: the driven steady cycle has strictly positive entropy production
J = pi_dr[0] * k_dr[(0, 1)] - pi_dr[1] * k_dr[(1, 0)]               # steady cycle current
affinity = sum(math.log(k_dr[(i, j)] / k_dr[(j, i)]) for (i, j) in RING)   # = 3f
sigma_ss = J * affinity
assert abs(affinity - 3 * 0.9) < 1e-12, "ring cycle affinity must be 3f (the nonconservative drive)"
assert J > 1e-9 and sigma_ss > 1e-9, "driven NESS must have a positive current and entropy production (the arrow)"
assert abs(pi_dr[0] * k_dr[(0, 1)] - pi_dr[1] * k_dr[(1, 0)]) > 1e-6, "drive must break KMS detailed balance"

# ---- (2c) NUMEROLOGY LOCUS: a common rate rescale is a pure clock change ----
c_clock = 12.34
k_fast = {e: c_clock * v for e, v in k_dr.items()}
pi_fast = stationary(k_fast)
for s in (0, 1, 2):
    assert abs(pi_fast[s] - pi_dr[s]) < 1e-12, "stationary law invariant under a common clock rescale"
st0 = trajs[0][0]
assert abs(entropy_production(st0, k_fast, pi_fast) - entropy_production(st0, k_dr, pi_dr)) < 1e-12, \
    "entropy production invariant under a common clock rescale -> absolute scale is ONE multiplier (the traffic clock)"

# ---------------------------------- verdict ----------------------------------
bar = "=" * 88
print(bar)
print("R1  THE RECORD-ACTION / MAXIMUM-CALIBER PRINCIPLE -- two forced legs (proposition)")
print(bar)
print(f"  (1) record-action unit  s1 = ln(8 x 137)      = {s1:.5f} k_B   (canon, -0.1 sigma)")
print(f"      record content      log2(8)               = {bits_per_cell:.3f} bits  (= item 144's 3 in 3*alpha)")
print(f"      per-event WEIGHT    exp(-s1)              = {w_event:.6e}  (= 1/1096)")
print(f"  (2b) KMS limit (f=0)    pi = Boltzmann(beta)  detailed balance OK   (DRIFT B1 form)")
print(f"  (2a) Crooks split       ln(P/P~) = dS_rec     traffic cancels; dwell-independent")
print(f"  (2c) clock rescale x{c_clock}  pi, dS invariant      -> absolute scale = 1 multiplier")
print(f"       driven NESS        affinity = 3f = {affinity:.3f}  sigma_ss = {sigma_ss:.4f} > 0  (arrow)")
print(f"""
{bar}
VERDICT (R1 forced legs, exit 0):
  The record-action measure's two forced legs hold: (1) its per-event quantum is the
  already-derived s1 = ln(8x137) (the entropy unit re-read as a weight), and (2) on a minimal
  monitored CTMC, entropy production fixes ONLY the time-antisymmetric part (traffic cancels in
  the Crooks ratio, dwell-independent), the equilibrium limit IS the canonical KMS scheduler
  (DRIFT B1), and a common clock rescale shows the absolute scale is a single traffic
  multiplier. This is why the parameter ledger can be favourable: the free multipliers are
  concentrated in A_traffic.

  TIER: forced-leg audit. FORCED + verified: the unit, the Crooks antisymmetry split, the
  KMS equilibrium limit, the scale-locus. Current R13 status is tracked by
  r13_record_action_measure_status.py: Born and bare alpha0 are no longer open R13 legs;
  the live make-or-break legs are the recovery holonomy Phi_rec (CP/sign) and the
  traffic-clock universalisation beyond the native selector lock. GUARDRAIL: per DRIFT E5
  + item 131 T2, this does NOT license inserting s1 or any KMS heat factor into an
  amplitude -- s1 is a per-event entropy/photon-bookkeeping unit only.
{bar}""")
print("exit 0 -- R1 forced legs verified (record-action unit; Crooks split + KMS limit + scale-locus).")
