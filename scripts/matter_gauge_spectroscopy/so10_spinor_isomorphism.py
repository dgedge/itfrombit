#!/usr/bin/env python3
"""
TQO Protocol 3 — the SO(10) spinor weight-lattice isomorphism (EXACT, no hardcoding).

Tests whether the framework's 16 codewords/generation map BIJECTIVELY onto the SO(10)
chiral spinor 16 = { (+-1/2)^5 : even number of -1/2 } (the even-minus weight set),
with the SM electric charge recovered from the 5 Cartan weights by the STANDARD
Pati-Salam SU(4)xSU(2)_L xSU(2)_R formula -- not by fitting.

Basis is the committed code (smg_code_projection.py, ANCHOR 2.1/2.2/2.8):
  bits c=(G0,G1,LQ,C0,C1,I3,CHI,W); one gen = (G0,G1)=(0,0); R1&R2&R3; charge() = 2.8.
SO(10) 16 = (4,2,1) [left, CHI=0] (+) (4bar,1,2) [right, CHI=1].
  color (e1,e2,e3): 4 = even-minus 3-tuple (lepton=singlet (+++)),
                    4bar = odd-minus (lepton singlet (---), = conjugate).
  weak  (e4,e5):    left  T3R=0 -> (s,s);  right T3L=0 -> (T3R,-T3R).
Charges from weights (standard): B-L = -2/3 (e1+e2+e3); T3L=(e4+e5)/2; T3R=(e4-e5)/2;
  Y = T3R + (B-L)/2;  Q_weight = T3L + Y.
Check: Q_weight == +charge(c) for left (CHI=0); == -charge(c) for right (CHI=1)
  (right SO(10) states are the conjugates e^c,u^c,d^c,nu^c; charge() ignores CHI).

numpy-free; exact Fraction arithmetic; self-asserting.
"""
import itertools
from fractions import Fraction as F

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
H = F(1, 2)


def R1(c): return not (c[G0] == 1 and c[G1] == 1)
def R2(c): return c[W] == c[CHI]
def R3(c):
    return (c[C0], c[C1]) == (0, 0) if c[LQ] == 0 else (c[C0], c[C1]) != (0, 0)


# committed 2.8 charge (independent of CHI/W -> a Dirac-style electric charge)
def Zf(c): return 1 if c[I3] == 0 else -1
def sumZc(c): return -3 if (c[C0], c[C1]) == (0, 0) else -1
def charge(c): return H * Zf(c) + F(1, 3) * sumZc(c) + H


# fixed colour-triplet labelling for the 3 quark colours (any fixed bijection -> same
# SU(3) triplet as a SET; the order is a basis choice, not physics)
QUARK_COLOURS = {(0, 1): 0, (1, 0): 1, (1, 1): 2}
TRIPLET_LEFT = [(+H, -H, -H), (-H, +H, -H), (-H, -H, +H)]   # the 3 non-singlet even-minus
SINGLET_LEFT = (+H, +H, +H)


def colour_weight(c, left):
    if c[LQ] == 0:                                  # lepton = SU(3) singlet (4th colour)
        w = SINGLET_LEFT
    else:
        w = TRIPLET_LEFT[QUARK_COLOURS[(c[C0], c[C1])]]
    return w if left else tuple(-x for x in w)      # 4bar = conjugate (odd-minus)


def weak_weight(c, left):
    if left:                                        # (4,2,1): T3R=0, T3L=+-1/2 from I3
        s = +H if c[I3] == 0 else -H
        return (s, s)
    # (4bar,1,2): T3L=0, T3R from conjugated I3 (u^c,nu^c have T3R=-1/2; d^c,e^c +1/2)
    t3r = -H if c[I3] == 0 else +H
    return (t3r, -t3r)


def weight5(c):
    left = (c[CHI] == 0)
    return colour_weight(c, left) + weak_weight(c, left)


def charges_from_weight(w):
    e1, e2, e3, e4, e5 = w
    BL = F(-2, 3) * (e1 + e2 + e3)
    T3L = (e4 + e5) / 2
    T3R = (e4 - e5) / 2
    Y = T3R + BL / 2
    return BL, T3L, T3R, Y, T3L + Y


# ---- build the generation and run the checks --------------------------------
gen = [c for c in itertools.product([0, 1], repeat=8)
       if (c[G0], c[G1]) == (0, 0) and R1(c) and R2(c) and R3(c)]
assert len(gen) == 16, f"expected 16 codewords, got {len(gen)}"

# the target: SO(10) chiral spinor 16 = even number of -1/2 over 5 slots
spinor16 = {w for w in itertools.product([+H, -H], repeat=5)
            if sum(1 for x in w if x == -H) % 2 == 0}
assert len(spinor16) == 16

weights = [weight5(c) for c in gen]
print("=" * 78)
print("SO(10) spinor weight-lattice isomorphism test")
print("=" * 78)
print(f"  framework codewords/gen : {len(gen)}")
print(f"  SO(10) spinor-16 target : {len(spinor16)} weights (even # of -1/2)")

# (1) bijection onto the spinor 16
assert len(set(weights)) == 16, "weights collide -> not a bijection"
assert set(weights) == spinor16, "image != SO(10) spinor 16"
print("\n  (1) BIJECTION: the 16 codewords map onto EXACTLY the even-minus spinor-16 set. PASS")

# (2) standard-formula electric charge recovered (sign = chirality/conjugation)
print("\n  (2) electric charge from the 5 Cartan weights (standard PS formula):")
print(f"      {'codeword':>9} {'spec':>4} {'CHI':>3}  {'weight(e1..e5)':>26}  {'B-L':>5} {'Y':>6} {'Q_w':>5} {'2.8 Q':>6}")
all_ok = True
for c in sorted(gen, key=lambda x: (x[CHI], x[LQ], x[I3], x[C0], x[C1])):
    w = weight5(c)
    BL, T3L, T3R, Y, Qw = charges_from_weight(w)
    left = (c[CHI] == 0)
    expected = charge(c) if left else -charge(c)
    spec = (("nu" if c[I3] == 0 else "e") if c[LQ] == 0 else ("u" if c[I3] == 0 else "d"))
    ok = (Qw == expected)
    all_ok &= ok
    wt = "(" + ",".join(f"{x}" for x in w) + ")"
    print(f"      {''.join(map(str,c)):>9} {spec:>4} {c[CHI]:>3}  {wt:>26}  {str(BL):>5} {str(Y):>6} {str(Qw):>5} {str(charge(c)):>6}  {'ok' if ok else 'XX'}")
assert all_ok, "electric charge mismatch under the standard SO(10) formula"
print("\n      Q_weight == +Q(2.8) for left (CHI=0), == -Q(2.8) for right (CHI=1, conjugates). PASS")

# (3) B-L spectrum matches the SM (sanity)
bl_left = {(("nu" if c[I3]==0 else "e") if c[LQ]==0 else ("u" if c[I3]==0 else "d")):
           charges_from_weight(weight5(c))[0] for c in gen if c[CHI]==0}
print(f"\n  (3) B-L (left fields): {{ {', '.join(f'{k}:{v}' for k,v in sorted(bl_left.items()))} }}  (lepton -1, quark +1/3). PASS")
assert bl_left == {"nu": F(-1), "e": F(-1), "u": F(1,3), "d": F(1,3)}

print("""
RESULT: the framework's 16 codewords ARE the SO(10) chiral spinor 16 -- an EXACT
weight-lattice isomorphism, with the SM electric charge reproduced from the 5 Cartan
weights by the standard Pati-Salam formula (no fitting).

SCOPE (honest):
 * This is the ANOMALY-FREE PRECONDITION (necessary), strengthening item 116 with the
   explicit Cartan/weight structure. It is NOT a bypass of Lüscher: cobordism-trivial
   is necessary, not sufficient; the dynamical SMG mirror-gapping stays open/negative
   (DRIFT K4). The Wen-Wang theorem's *applicability* is confirmed; its *realisation*
   is not delivered by a static weight map.
 * The isomorphism is STATIC. The Z2xZ2-vs-Z3 colour issue (ANCHOR L209,
   matter_coupled_proxy.py) lives at the GAUGING/dynamics level (the SU(3) centre /
   triality acting on the qubits), which a weight-lattice bijection does not touch.
   The 3 quark colours form a proper SU(3) triplet AS A SET here; whether the gauge
   *dynamics* realises the Z3 centre is the separate, still-open question.
ALL ASSERTS PASSED.
""")
