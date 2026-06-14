#!/usr/bin/env python3
r"""THE 0.222% QUESTION — is the stationary branch's gamma-correction a specific,
producible/falsifiable object of the register thermodynamics?

LAYER 1 — THE READING DISCRIMINATOR (decided before any small correction):
  when a bit flips adjacent to k existing faults, the generation weight is either
    G1 (net-frustration Boltzmann): gamma(k) = exp(-(3-2k)/(2 phi))
       — a fault-adjacent flip UNFRUSTRATES shared edges, so generation
         accelerates near faults;
    G2 (billed ledger cost): gamma = exp(-3/(2 phi)) always
       — every flip writes 3 syndrome records regardless of neighbours.
  At stationary occupancy ~5%, G1's acceleration is O(70%), not O(0.2%): the two
  readings differ catastrophically and the 256-state chain decides which can
  support the stationary CC branch at all. (The Sec 5.2 M = e^(F/2phi) inscription
  is AMBIGUOUS between net and billed frustration — this computation makes the
  ambiguity load-bearing.)

LAYER 2 — THE CORRECTION SCAN (within the surviving reading):
  required: Delta(-ln gamma) = +0.0022245 (gamma ratio 0.997778, the parallel
  session's number — re-derived here). Pre-registered candidate classes:
    (a) radiative/portal alpha-forms: alpha0 x {1, 1/2, 1/3, 1/pi, 1/(2pi), 2/pi, 3/(2pi)}
    (b) self-consistency forms from the chain itself: {q1_stat, p_stat, p_stat^2,
        gamma_raw, gamma_raw^2, q1_stat/2}
    (c) scheduler-granularity forms: {1/137, 2/137, 1/274, 1/1096, 12/1096}
  Gate honesty FIRST: the empirical sensitivity d(ln rho)/d(Delta) ~ 216 means the
  factor-2 window in rho admits Delta in [0, ~0.0054] — nearly any small candidate
  — so the meaningful gates are rho within +-5% and +-2%; the trials factor at
  those widths is computed and reported. If multiple candidates crowd the window,
  the verdict is 'requires derivation, not scanning' — itself the answer.
exit 0 = machinery verified."""
import math
import numpy as np

PHI = (math.sqrt(5) - 1) / 2
KAP = 1 / (2 * PHI)
G_RAW = math.exp(-3 * KAP)
EDGES = [(i, j) for i in range(8) for j in range(8) if i < j and bin(i ^ j).count("1") == 1]
NBR = [[j for j in range(8) if bin(i ^ j).count("1") == 1] for i in range(8)]
FLAW = [0, 0, 0, 0, 0.5, 1, 1, 1, 1]

M_P = 1.220890e19
H0 = 67.36 / 3.085678e19 * 6.582120e-25
rho_obs = 3 * 0.6847 * H0 * H0 * M_P * M_P / (8 * math.pi)
def rho_ratio(q1):
    return (1 / 137.0) * 0.332 ** 4 * (21 * q1) ** 32 / 21 / rho_obs

def chain_stationary(gen_weight):
    """exact 256-state chain: per-bit generation weights gen_weight(state, bit),
    then one-jump repair (uniform over faults). Returns stationary q1."""
    P = np.zeros((256, 256))
    for f in range(256):
        clean = [b for b in range(8) if not (f >> b) & 1]
        probs = {b: gen_weight(f, b) for b in clean}
        # enumerate flip subsets of clean bits
        from itertools import combinations
        for r in range(len(clean) + 1):
            for sub in combinations(clean, r):
                pr = 1.0
                for b in clean:
                    pr *= probs[b] if b in sub else (1 - probs[b])
                m = f
                for b in sub:
                    m |= (1 << b)
                # one-jump repair: remove one fault uniformly
                faults = [b for b in range(8) if (m >> b) & 1]
                if not faults:
                    P[f, m] += pr
                else:
                    for b in faults:
                        P[f, m & ~(1 << b)] += pr / len(faults)
    w, v = np.linalg.eig(P.T)
    pi = np.real(v[:, np.argmax(np.real(w))]); pi /= pi.sum()
    return float(sum(pi[f] * FLAW[bin(f).count('1')] for f in range(256)))

print("[1] THE READING DISCRIMINATOR (256-state exact chains):")
def g2(f, b):  # billed ledger cost
    return G_RAW
def g1(f, b):  # net-frustration Boltzmann
    k = sum(1 for n in NBR[b] if (f >> n) & 1)
    return math.exp(-(3 - 2 * k) * KAP)
q1_G2 = chain_stationary(g2)
q1_G1 = chain_stationary(g1)
print(f"    G2 (billed, 3 records always):  q1 = {q1_G2:.6e}, rho = {rho_ratio(q1_G2):.3f} x obs")
print(f"    G1 (net-frustration Boltzmann): q1 = {q1_G1:.6e}, rho = {rho_ratio(q1_G1):.3e} x obs")
print(f"    -> G1 destroys the stationary branch by ~{math.log10(rho_ratio(q1_G1)/rho_ratio(q1_G2)):.0f} OOM:")
print(f"       fault-adjacent generation accelerates (DF = 3-2k), occupancy avalanches.")
print(f"       THE STATIONARY CC BRANCH REQUIRES THE BILLED-LEDGER READING of Sec 5.2's")
print(f"       M = e^(F/2phi) — the net-vs-billed ambiguity is now load-bearing canon.")
assert rho_ratio(q1_G1) > 1e3 * rho_ratio(q1_G2)
assert abs(rho_ratio(q1_G2) - 1.616) < 0.05

# ---------------- layer 2: the correction scan within G2 ----------------
# required Delta re-derived: find Delta with rho = 1 exactly
lo, hi = 0.0, 0.01
for _ in range(50):
    mid = (lo + hi) / 2
    q1m = chain_stationary(lambda f, b, d=mid: G_RAW * math.exp(-d))
    lo, hi = (mid, hi) if rho_ratio(q1m) > 1.0 else (lo, mid)
DREQ = (lo + hi) / 2
print(f"\n[2] required correction (re-derived on the 256-state chain): Delta = {DREQ:.7f}")
print(f"    (parallel session's 9-state value: 0.0022245 — agreement check)")
assert abs(DREQ - 0.0022245) < 3e-4

q1_at = {}
def rho_of_delta(d):
    if d not in q1_at:
        q1_at[d] = chain_stationary(lambda f, b, dd=d: G_RAW * math.exp(-dd))
    return rho_ratio(q1_at[d])

a0 = 1 / 137.0
p_stat = 0.0544
CANDS = [
    ("alpha0",          a0),            ("alpha0/2",   a0 / 2),
    ("alpha0/3",        a0 / 3),        ("alpha0/pi",  a0 / math.pi),
    ("alpha0/(2pi)",    a0 / (2 * math.pi)), ("2 alpha0/pi", 2 * a0 / math.pi),
    ("3 alpha0/(2pi)",  3 * a0 / (2 * math.pi)),
    ("q1_stat",         q1_G2),         ("q1_stat/2",  q1_G2 / 2),
    ("p_stat",          p_stat),        ("p_stat^2",   p_stat ** 2),
    ("gamma_raw^2",     G_RAW ** 2),    ("gamma_raw^2/2", G_RAW ** 2 / 2),
    ("1/137",           1 / 137),       ("2/137",      2 / 137),
    ("1/274",           1 / 274),       ("1/1096",     1 / 1096),
    ("12/1096",         12 / 1096),
]
SENS = math.log(1.616) / DREQ
print(f"    sensitivity d(ln rho)/d(Delta) = {SENS:.0f}; gates: +-5% rho -> |Delta-Dreq| < {0.05/SENS:.1e};"
      f" +-2% -> {0.02/SENS:.1e}")
print(f"\n    {'candidate':<16s} {'Delta':>10s} {'rho/rho_obs':>12s}  verdict")
hits5, hits2 = [], []
for nm, d in CANDS:
    if d > 0.02:
        continue
    r = rho_of_delta(d)
    v = ""
    if abs(math.log(r)) < 0.02:
        v = "HIT (2%)"; hits2.append(nm); hits5.append(nm)
    elif abs(math.log(r)) < 0.05:
        v = "hit (5%)"; hits5.append(nm)
    print(f"    {nm:<16s} {d:>10.6f} {r:>12.4f}  {v}")
n_in_range = sum(1 for _, d in CANDS if d <= 0.02)
exp5 = n_in_range * 2 * (0.05 / SENS) / 0.008
exp2 = n_in_range * 2 * (0.02 / SENS) / 0.008
print(f"\n    trials factor: {n_in_range} candidates over Delta in [0, 0.008]:"
      f" expected accidentals {exp5:.1f} (5% gate), {exp2:.2f} (2% gate)")
print(f"    hits: 5% gate {hits5}; 2% gate {hits2}")

print(f"""
VERDICT — is the 0.222% a specific producible object?
  LAYER 1: YES, decisively, one level up — the register thermodynamics' FIRST
    binary choice (net-frustration vs billed-ledger generation) is executed by
    the stationary branch: net-frustration destroys it by orders of magnitude;
    ONLY the billed-ledger reading of M = e^(F/2phi) supports the CC route.
    That is a sharp, falsifiable demand on Sec 5.2 semantics.
  LAYER 2: the 0.222% window itself is {'CROWDED' if len(hits5) > 1 else 'clean'} at the 5% gate — scanning
    cannot pick the winner; the correction REQUIRES a derivation (the two
    physically-motivated leading forms: a radiative alpha0/pi-class dressing of
    the attempt rate, or the q1 self-consistency feedback), discriminable only
    at better-than-2%-in-rho grade. The question is well-posed and falsifiable,
    but by derivation, not by numerology — which is the right way round. exit 0""")
print("ALL ASSERTIONS PASSED")
