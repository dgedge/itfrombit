#!/usr/bin/env python3
r"""THE d <-> p MAP — the exact cell-failure law, the correlation-class theorem, and
what the measured K04 plateau (d ~ 0.30) actually implies for the boot recurrence.

[1] THE EXACT CELL FAILURE LAW (from the strain-decoder theorem, re-verified here):
    a crystallised cell with k bit-faults under full strain readout:
       k <= 3 : corrected exactly (coboundary injective)        -> f(k) = 0
       k  = 4 : complement confusion, tie-break                 -> f(4) = 1/2
       k >= 5 : min-weight decoding applies the COMPLEMENT (weight 8-k <= 3),
                completing the CPT flip with certainty          -> f(k) = 1
    So q1 = E[f(k)] over the fault-placement ensemble. The familiar (35/2)p^4 is just
    the leading term of the exact i.i.d. law.

[2] THE CORRELATION-CLASS THEOREM: at fixed global fault fraction d, q1 is NOT a
    function of d — it spans [0, ~d] depending on placement:
       capped-spread (<= 3 per cell, possible iff d <= 3/8)  -> q1 = 0
       i.i.d.                                                -> binomial law
       supra-cell clusters (whole cells faulted)             -> q1 = d
    d alone determines NOTHING; the correlation class decides everything.

[3] THE MEASURED K04 DEFECTS ARE SUPRA-CELL (whole non-cube components), so IF the
    wrecked material enters the code layer: q1 ~ d ~ 0.30 -> 21 q1 >> 1 -> the
    concatenation DIVERGES: that reading is dead. The living reading: registers
    instantiate only on crystallised cells (the wreckage is substrate-dark — canon's
    own R4/dark-sector narrative); K04 then measures the YIELD, not p; the boot p is
    the BIT-layer noise on good cells — a different, unmeasured observable.

[4] REQUIREMENT REFINEMENT: with the exact law, the rho_Lambda closure requirement
    moves from p = 0.1088 (the 17.5 p^4 truncation) to the exact-law inversion.
Self-asserting; exit 0 = every number verified."""
import itertools as it
import math
import numpy as np

# ---------------- [1] re-verify the failure law from the coboundary classes ----------------
pts = list(range(8))
EDGES = [(u, v) for u in pts for v in pts if u < v and bin(u ^ v).count("1") == 1]
def syn(S):
    v = [1 if i in S else 0 for i in range(8)]
    return tuple((v[a] + v[b]) % 2 for a, b in EDGES)
by_syn = {}
for w in range(9):
    for S in it.combinations(range(8), w):
        by_syn.setdefault(syn(set(S)), []).append(frozenset(S))
def fail_prob(k):
    """P(logical fault | k faults, min-weight strain decoding), exact."""
    tot = fail = 0.0
    for S in it.combinations(range(8), k):
        S = frozenset(S)
        cls = by_syn[syn(set(S))]
        lightest = min(len(T) for T in cls)
        mins = [T for T in cls if len(T) == lightest]
        # residual after correcting with T: S xor T in {empty, V}; fault iff residual = V
        p_fault = sum(1.0 for T in mins if T != S) / len(mins)
        tot += 1; fail += p_fault
    return fail / tot
F = {k: fail_prob(k) for k in range(9)}
print("[1] EXACT CELL FAILURE LAW f(k) (re-verified from the coboundary classes):")
print("    " + "  ".join(f"f({k})={F[k]:.2f}" for k in range(9)))
assert all(F[k] == 0.0 for k in range(4)) and F[4] == 0.5 and all(F[k] == 1.0 for k in range(5, 9))

# ---------------- [2] the correlation-class theorem ----------------
def q1_iid(p):
    return sum(math.comb(8, k) * p**k * (1 - p)**(8 - k) * F[k] for k in range(9))
d = 0.30
print(f"\n[2] CORRELATION-CLASS THEOREM at fixed d = {d}:")
print(f"    capped-spread (<=3/cell; possible since d <= 3/8): q1 = 0   (all corrected)")
print(f"    i.i.d.:                                            q1 = {q1_iid(d):.4f}")
print(f"    supra-cell clusters (whole cells):                 q1 = {d:.4f}")
print(f"    -> q1 spans [0, {d}]: d alone determines NOTHING; the placement class does.")
assert d <= 3 / 8 and 0 < q1_iid(d) < d

# ---------------- [3] the two readings of the measured plateau ----------------
thr = 1 / 21
print(f"\n[3] THE MEASURED PLATEAU (d ~ 0.30, supra-cell by construction — whole non-cube")
print(f"    components): IF the wreckage enters the code layer, q1 ~ 0.30 vs the pair-")
print(f"    recursion threshold 1/21 = {thr:.4f}: 21 q1 = {21*d:.1f} >> 1 — DIVERGES. Even the")
print(f"    i.i.d. re-reading at d = 0.30 gives q1 = {q1_iid(d):.3f} > {thr:.4f}: dead either way.")
print(f"    THE LIVING READING: registers instantiate only on CRYSTALLISED cells; the")
print(f"    wrecked fraction is substrate-dark (canon's own R4/dark-sector narrative —")
print(f"    flagged observation, NOT adopted: the measured excluded fraction ~0.30 sits")
print(f"    in the same decade as the dark fraction; class-1, vertex-vs-energy categories")
print(f"    differ). K04 then measures the YIELD (~0.70), not p; the boot p is the")
print(f"    bit-layer noise on good cells — unmeasured by the K04 sweep.")
assert 21 * d > 1 and q1_iid(d) > thr

# ---------------- [4] the requirement, refined by the exact law ----------------
alpha0, LAM = 1 / 137.0, 0.332
M_P = 1.220890e19
H0 = 67.36 / 3.085678e19 * 6.582120e-25
rho_obs = 3 * 0.6847 * H0 * H0 * M_P * M_P / (8 * math.pi)
q_top = rho_obs / (alpha0 * LAM ** 4)
lg_q1req = (math.log(21) + math.log(q_top)) / 2 ** 5 - math.log(21)
q1_req = math.exp(lg_q1req)
lo, hi = 0.01, 0.3
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if q1_iid(mid) < q1_req else (lo, mid)
p_req_exact = (lo + hi) / 2
p_req_old = (q1_req / 17.5) ** 0.25
print(f"\n[4] REQUIREMENT REFINEMENT (chain: q1_req = {q1_req:.3e} at ell = 6):")
print(f"    truncated law (35/2)p^4:  p_req = {p_req_old:.4f}  (the registered 0.1088)")
print(f"    EXACT law (k=4..8 terms): p_req = {p_req_exact:.4f}")
print(f"    The exact law moves the requirement {100*(p_req_exact/p_req_old-1):+.1f}% — AWAY from the")
print(f"    historical '0.110': the right-number-wrong-role coincidence WEAKENS (honest).")
assert abs(p_req_old - 0.1088) < 0.0005
assert 0.090 < p_req_exact < 0.105

print(f"""
VERDICT — the 'd <-> p map' out, executed:
  * The map is now EXACT at the cell level (the failure law f(k) = 0/0/0/0/.5/1/1/1/1,
    re-verified), and it is a CORRELATION-CLASS map, not a function: q1 in [0, d].
  * The out as originally hoped — 'correlated defects at d = 0.30 might correspond to
    a small effective i.i.d. p' — FAILS in both directions that keep the wreckage in
    the code layer: supra-cell gives q1 = 0.30, i.i.d. re-reading gives 0.126; both
    exceed the 1/21 threshold; the chain diverges. NO placement of d = 0.30 INTO the
    code layer rescues the closure except the capped-spread class, which the measured
    whole-component defect structure explicitly is NOT.
  * What survives is the CATEGORY SPLIT: crystallisation YIELD (what K04 measures,
    ~0.70 good) vs bit-layer NOISE on good cells (what the recurrence's p actually
    is — unmeasured). The rho_Lambda requirement re-targets the bit layer:
    p(boot) = {p_req_exact:.4f} (exact-law refinement of the registered 0.1088), and the
    eta-derived T8 bound (p < 1.0e-3/tick at equilibrium) becomes its late-time
    complement: the framework now owes a NOISE-ANNEALING TRAJECTORY p(t) from
    ~{p_req_exact:.2f} at boot to < 1e-3 at equilibrium — a named, falsifiable object.
  * Flagged observation (class-1, NOT adopted): the excluded fraction ~0.30 vs the
    dark fraction — canon's own R4 narrative owns the qualitative story; categories
    (vertex vs energy, boot vs today) differ. exit 0""")
print("ALL ASSERTIONS PASSED — every number above is verified.")
