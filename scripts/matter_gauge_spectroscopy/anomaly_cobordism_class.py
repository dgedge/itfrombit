#!/usr/bin/env python3
"""
TQO Target A — the trivial cobordism class, COMPUTED comprehensively (not cited).

Derives the full SM irrep content of one generation from the P3 SO(10)-spinor weight
map (so10_spinor_isomorphism.py), then computes EVERY anomaly that defines the
relevant cobordism/anomaly class:

  A1  U(1)_Y^3
  A2  grav^2-U(1)_Y
  A3  [SU(2)_L]^2-U(1)_Y
  A4  [SU(3)_c]^2-U(1)_Y
  A5  [SU(3)_c]^3
  A6  Witten SU(2)_L global (# doublets mod 2)
  A7  global Z_16 (Wang-Wen 2018): #Weyl mod 16

All vanish => the [8,4,4] boundary content is the TRIVIAL cobordism class. This
strengthens ANCHOR 2.11 + smg_code_projection.py by computing the COMPLETE set from
the explicit P3 charges rather than citing A1..A6.

KEY CONSEQUENCE (printed): trivial cobordism class == TRIVIAL bulk SPT == ZERO anomaly
inflow. So 'anomaly inflow nullifies the measure obstruction' is SELF-CONTRADICTORY for
this content: a trivial class has no inflow to nullify anything. The SM is anomaly-free
in 4D; the bulk's role is regularisation (domain-wall), and the residual obstruction
(decouple the mirror gauge-invariantly) is SMG (dynamical, fails for the framework's CSS
realisation) or the overlap measure (Lüscher, open) -- the SAME wall, unmoved.

exact Fraction arithmetic; numpy-free; self-asserting.
"""
import itertools
from fractions import Fraction as F

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
H = F(1, 2)


def R1(c): return not (c[G0] == 1 and c[G1] == 1)
def R2(c): return c[W] == c[CHI]
def R3(c):
    return (c[C0], c[C1]) == (0, 0) if c[LQ] == 0 else (c[C0], c[C1]) != (0, 0)


# --- P3 weight map -> left-Weyl SO(10) hypercharge per codeword --------------
TRIPLET_L = [(+H, -H, -H), (-H, +H, -H), (-H, -H, +H)]
SINGLET_L = (+H, +H, +H)
QC = {(0, 1): 0, (1, 0): 1, (1, 1): 2}


def colour_w(c, left):
    w = SINGLET_L if c[LQ] == 0 else TRIPLET_L[QC[(c[C0], c[C1])]]
    return w if left else tuple(-x for x in w)


def weak_w(c, left):
    if left:
        s = +H if c[I3] == 0 else -H
        return (s, s)
    t3r = -H if c[I3] == 0 else +H
    return (t3r, -t3r)


def YT(c):
    left = (c[CHI] == 0)
    e1, e2, e3 = colour_w(c, left)
    e4, e5 = weak_w(c, left)
    BL = F(-2, 3) * (e1 + e2 + e3)
    T3R = (e4 - e5) / 2
    return BL / 2 + T3R                     # Y (left-Weyl SO(10) hypercharge)


# --- derive the SM irrep content of one generation ---------------------------
gen = [c for c in itertools.product([0, 1], repeat=8)
       if (c[G0], c[G1]) == (0, 0) and R1(c) and R2(c) and R3(c)]
assert len(gen) == 16


def colour_rep(c):
    if c[LQ] == 0: return "1"
    return "3" if c[CHI] == 0 else "3bar"


def su2(c): return "2" if c[CHI] == 0 else "1"


DIM3 = {"1": 1, "3": 3, "3bar": 3}
groups = {}
for c in gen:
    key = (colour_rep(c), su2(c), YT(c))
    groups[key] = groups.get(key, 0) + 1

# convert codeword counts -> rep instances (count / rep dimension)
reps = []
for (cr, w2, Y), n in groups.items():
    dim = DIM3[cr] * (2 if w2 == "2" else 1)
    assert n % dim == 0, f"group {(cr,w2,Y)} count {n} not a multiple of dim {dim}"
    reps.append((cr, w2, Y, n // dim))

print("Derived SM irrep content of one generation (from P3 weights):")
print(f"  {'SU(3)':>5} {'SU(2)':>5} {'Y':>6} {'inst':>4}")
for cr, w2, Y, inst in sorted(reps, key=lambda r: (r[0], r[1], r[2])):
    print(f"  {cr:>5} {w2:>5} {str(Y):>6} {inst:>4}")
assert sorted((cr, w2, Y, i) for cr, w2, Y, i in reps) == sorted([
    ("1", "1", F(0), 1), ("1", "1", F(1), 1), ("1", "2", F(-1, 2), 1),
    ("3", "2", F(1, 6), 1), ("3bar", "1", F(-2, 3), 1), ("3bar", "1", F(1, 3), 1)])
print("  => exactly Q + u^c + d^c + L + e^c + nu^c  (the SO(10) 16). ✓")

# --- anomaly coefficients ----------------------------------------------------
T_FUND = F(1, 2)                              # Dynkin index of SU(N) fundamental
A_SU3 = {"1": 0, "3": 1, "3bar": -1}         # cubic anomaly coefficient


def d3(cr): return DIM3[cr]
def d2(w2): return 2 if w2 == "2" else 1


A1 = sum(inst * d3(cr) * d2(w2) * Y**3 for cr, w2, Y, inst in reps)            # U(1)^3
A2 = sum(inst * d3(cr) * d2(w2) * Y for cr, w2, Y, inst in reps)               # grav^2 U(1)
A3 = sum(inst * d3(cr) * T_FUND * Y for cr, w2, Y, inst in reps if w2 == "2")  # SU(2)^2 U(1)
A4 = sum(inst * d2(w2) * T_FUND * Y for cr, w2, Y, inst in reps if cr != "1")  # SU(3)^2 U(1)
A5 = sum(inst * d2(w2) * A_SU3[cr] for cr, w2, Y, inst in reps)                # SU(3)^3
A6 = sum(inst * d3(cr) for cr, w2, Y, inst in reps if w2 == "2") % 2           # Witten
A7 = sum(inst * d3(cr) * d2(w2) for cr, w2, Y, inst in reps) % 16              # Z_16

print("\nAnomaly coefficients (one generation, all left-Weyl):")
for name, val, exp in [("A1  U(1)_Y^3", A1, 0), ("A2  grav^2-U(1)_Y", A2, 0),
                       ("A3  [SU(2)_L]^2-U(1)_Y", A3, 0), ("A4  [SU(3)_c]^2-U(1)_Y", A4, 0),
                       ("A5  [SU(3)_c]^3", A5, 0), ("A6  Witten SU(2)_L (mod 2)", A6, 0),
                       ("A7  global Z_16 (#Weyl mod 16)", A7, 0)]:
    print(f"   {name:<28} = {str(val):>4}   [expect {exp}]")
assert A1 == 0 and A2 == 0 and A3 == 0 and A4 == 0 and A5 == 0 and A6 == 0 and A7 == 0
print("\n  ALL anomalies vanish => the [8,4,4] boundary content is the TRIVIAL cobordism class.")

print("""
WHAT THIS DOES AND DOES NOT DO
==============================
DONE (Target A, comprehensively & computed): every gauge/gravitational/global anomaly
  of the framework's 16 vanishes. The boundary fermion content is the trivial cobordism
  class. (Strengthens 2.11 + smg_code_projection: full set computed from P3 charges.)

DOES NOT close the chiral frontier -- and the reason is the SPT argument's own logic:
  * trivial cobordism class  ==  TRIVIAL bulk SPT  ==  ZERO anomaly inflow.
  * 'anomaly inflow nullifies the measure obstruction' needs a NON-trivial bulk (nonzero
    inflow). That is mutually exclusive with a trivial class. For anomaly-free content
    (the SM) there is NOTHING to inflow -- the SM is anomaly-free in 4D and needs no bulk
    for its consistency.
  * So the bulk's role here is REGULARISATION (domain-wall), not anomaly cancellation.
    Domain-wall == overlap; Lüscher's measure problem is FORMULATED in exactly this
    construction, not bypassed by it. The residual obstruction is decoupling the
    doubler/mirror GAUGE-INVARIANTLY: SMG (dynamical -- the framework's CSS realisation
    collapses, DRIFT K4) or the overlap measure (non-abelian, open).
  * Anomaly matching is regularisation-independent, so the obstruction is invariant under
    the continuous<->discrete<->SPT relabelling. Target A (precondition) is met; the wall
    (dynamical mirror-decoupling) is unmoved.
ALL ASSERTS PASSED.
""")
