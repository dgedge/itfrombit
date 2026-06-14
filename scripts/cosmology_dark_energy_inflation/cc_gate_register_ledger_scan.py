#!/usr/bin/env python3
r"""THE CC GATE — structure-first scan of canon-citable register-handoff ledgers
against p* = p_req +- 0.0006 (the user's gate; anything outside is falsified).

THE TWO LEDGER FORMS (each with its own target exponent):
  BALANCE  (fault creation vs thermal removal, detailed balance):
      p* = 1/(1+e^x),  x_target = ln((1-p_req)/p_req)
  SLAVING  (per-tick full correction resets faults; standing p = per-tick creation
      — THE CANON-FAVOURED FORM, derived in noise_annealing_trajectory.py):
      p* = e^(-x),     x_target = -ln(p_req)

PRE-REGISTERED CANDIDATE CLASSES (cited; no post-hoc additions):
  CLASS K — energy forms x = w_s * kappa: syndrome weight w_s in {1,2,3,4,6}
      (single bit trips 3 edges; pairs trip 4 or 6; sub-edge units 1,2),
      kappa in the canonical constant set:
        1/(2 phi)  [Sec 5.2 Group II inscription]
        phi, 1/phi, phi^2, 1/phi^2   [Group II family]
        ln 2       [bit entropy]
        1, 1/2, 3/4 [rational service weights]
  CLASS C — count-logarithm forms x = ln(C)/m, m in {1,2,3,4}, C in the PRIMARY
      canon counts only: {8, 12, 14, 21, 28, 35, 48, 64, 70, 137, 256, 1096}
      (alphabet, edges, hyperplanes, pair coefficient, channels, complement
      classes, codewords, register dim, weight-4 sets, fine structure, states,
      8x137 = the ADOPTED record alphabet of the s1 closure).
  CLASS R — bare rationals p = a/b, b <= 144: enumerated ONLY to measure the
      coincidence density at this tolerance (Sec 16.3 discipline): the class is
      reported and DISMISSED as non-discriminating.
Trials-factor accounting is printed per class. exit 0 = scan verified."""
import math

# ---------------- the requirement, recomputed from the chain ----------------
def f_k(k): return 0.0 if k <= 3 else (0.5 if k == 4 else 1.0)
def q1_of(p): return sum(math.comb(8, k) * p**k * (1 - p)**(8 - k) * f_k(k) for k in range(9))
M_P = 1.220890e19
H0 = 67.36 / 3.085678e19 * 6.582120e-25
rho_obs = 3 * 0.6847 * H0 * H0 * M_P * M_P / (8 * math.pi)
q_top = rho_obs / ((1 / 137.0) * 0.332 ** 4)
q1_req = math.exp((math.log(21) + math.log(q_top)) / 32 - math.log(21))
lo, hi = 0.01, 0.3
for _ in range(80):
    mid = (lo + hi) / 2
    lo, hi = (mid, hi) if q1_of(mid) < q1_req else (lo, mid)
p_req = (lo + hi) / 2
TOL = 0.0006
x_bal = math.log((1 - p_req) / p_req)
x_slv = -math.log(p_req)
# gate width in x for each form: dx = dp/(p(1-p)) (balance), dp/p (slaving)
dx_bal = TOL / (p_req * (1 - p_req))
dx_slv = TOL / p_req
print(f"[0] p_req = {p_req:.6f} +- {TOL}  (DF/T balance target {x_bal:.5f} +- {dx_bal:.5f};"
      f" slaving target {x_slv:.5f} +- {dx_slv:.5f})")
assert abs(p_req - 0.0972) < 0.0005

PHI = (math.sqrt(5) - 1) / 2
KAPPAS = [("1/(2phi) [5.2 GrpII]", 1 / (2 * PHI)), ("phi", PHI), ("1/phi", 1 / PHI),
          ("phi^2", PHI ** 2), ("1/phi^2", 1 / PHI ** 2), ("ln2", math.log(2)),
          ("1", 1.0), ("1/2", 0.5), ("3/4", 0.75)]
WS = [1, 2, 3, 4, 6]
COUNTS = [8, 12, 14, 21, 28, 35, 48, 64, 70, 137, 256, 1096]
MS = [1, 2, 3, 4]

def rho_factor(p):
    """rho(p)/rho_obs through the exact law (the 120-power sensitivity)."""
    return (q1_of(p) / q1_req) ** 32

def scan(cands, form):
    x_t, dx = (x_bal, dx_bal) if form == "balance" else (x_slv, dx_slv)
    hits = []
    for nm, x in cands:
        if abs(x - x_t) < dx:
            p = 1 / (1 + math.exp(x)) if form == "balance" else math.exp(-x)
            hits.append((nm, x, p, rho_factor(p)))
    return hits

# CLASS K
candK = []
seen = set()
for w in WS:
    for nm, k in KAPPAS:
        x = w * k
        key = round(x, 6)
        if key in seen or not 0.5 < x < 4.5:
            continue
        seen.add(key)
        candK.append((f"{w} x {nm}", x))
# CLASS C
candC = []
for C in COUNTS:
    for m in MS:
        x = math.log(C) / m
        key = round(x, 6)
        if key in seen or not 0.5 < x < 4.5:
            continue
        seen.add(key)
        candC.append((f"ln({C})/{m}", x))

print(f"\n[1] CLASS K (energy forms): {len(candK)} candidates over x in (0.5, 4.5)")
print(f"    CLASS C (count-log forms): {len(candC)} candidates")
RANGE = 4.0
for form in ("balance", "slaving"):
    dx = dx_bal if form == "balance" else dx_slv
    expK = len(candK) * 2 * dx / RANGE
    expC = len(candC) * 2 * dx / RANGE
    print(f"\n  -- {form.upper()} form (expected accidentals: K {expK:.2f}, C {expC:.2f}) --")
    for label, cands in (("K", candK), ("C", candC)):
        for nm, x, p, rf in scan(cands, form):
            print(f"    HIT [{label}] {nm}: x = {x:.5f}, p* = {p:.6f}, "
                  f"rho_Lambda = {rf:.3f} x observed")
    if not (scan(candK, form) + scan(candC, form)):
        print(f"    no hits")

# CLASS R — coincidence density only
nR = 0
for b in range(2, 145):
    for a in range(1, b):
        if abs(a / b - p_req) < TOL:
            nR += 1
print(f"\n[2] CLASS R (bare rationals, b <= 144): {nR} rationals inside the gate")
print(f"    -> the rational class is NON-DISCRIMINATING at this tolerance (16.3):")
print(f"       any p-rational 'hit' (e.g. 7/72 = {7/72:.6f}) carries no evidential weight")
print(f"       without a derived ledger producing it. Class dismissed.")
assert nR > 3

print(f"""
[3] VERDICT PROTOCOL (pre-registered): a structured-class hit is a CANDIDATE
    (class-2), not a derivation, until the register-handoff theorem produces its
    exponent from the 8-bit/12-edge service ledger. The trials factors above
    quantify exactly how seriously to take any single hit. exit 0""")
print("ALL ASSERTIONS PASSED — scan verified.")
