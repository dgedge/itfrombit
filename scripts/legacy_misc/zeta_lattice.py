#!/usr/bin/env python3
"""
zeta(3) from a discrete LATTICE iterated integral — proof-of-route for the
refined a_e obstruction (§15 item 138 / ae_test.py).

Context. ae_test.py shows §5.4's combinatorial COUNT (S4 permutations over beta1
homology classes -> -24/7, a rational; geometric series -> 21) cannot reach the
QED 2-loop coefficient c2 = 197/144 + pi^2/12 - (pi^2/2)ln2 + (3/4)zeta(3), which
contains the transcendental period zeta(3). The refined diagnosis: the obstruction
is SHORTCUT-SPECIFIC, not substrate-fundamental — §5.4 replaced an iterated
integral with a count, and zeta(3) is precisely what the iterated dx/x integral
PRODUCES. This script demonstrates the constructive half: a discrete lattice
evaluation of the iterated integral converges to zeta(3) natively, with no period
inserted by hand.

What is and is NOT shown:
  SHOWN  — the iterated-integral machinery class reaches the QED period zeta(3)
           on a discrete grid (O(h^2) convergence, verified by Richardson).
  NOT    — that the TCH substrate's specific electron vertex graph F2(0) yields
           the QED coefficient. That is the real open computation (build the
           D_TCH vertex, integrate it). This is proof-of-route only.

Identity used: zeta(3) = (1/2) \int_0^1 ln^2(x)/(1-x) dx   [standard period form;
the ln^2 x supplies two dx/x integrations, 1/(1-x) = sum x^n the third].
Substitution x = e^{-t} gives the smooth, singularity-free integrand
  (1/2) \int_0^infty t^2 e^{-t}/(1 - e^{-t}) dt
which a plain uniform-lattice trapezoidal sum evaluates cleanly.

Needs only mpmath. Nothing fitted.
"""
import mpmath as mp
mp.mp.dps = 30
z3 = mp.zeta(3)

print("target  zeta(3) =", mp.nstr(z3, 20))

# --- 1. continuous identity, adaptive (tanh-sinh) quadrature on [0,1] ---
I_ad = mp.mpf(1)/2 * mp.quad(lambda x: mp.log(x)**2 / (1 - x), [0, 1])
print("\n[1] adaptive quadrature  (1/2) int_0^1 ln^2 x/(1-x) dx")
print("    =", mp.nstr(I_ad, 20), "  matches zeta(3):", mp.almosteq(I_ad, z3),
      "  err =", mp.nstr(abs(I_ad - z3), 3))

# --- 2. genuine UNIFORM LATTICE (trapezoidal in t after x=e^{-t}) ---
def integrand(t):
    return t**2 * mp.e**(-t) / (1 - mp.e**(-t))   # smooth; ->0 at t=0 and t->inf

def lattice(h, T=80):
    N = int(T / h)
    return mp.mpf(1)/2 * h * sum(integrand(k * h) for k in range(1, N + 1))  # f(0)=0

print("\n[2] uniform lattice in t (spacing h, trapezoidal); a real discrete grid:")
for h in [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.02')]:
    Ih = lattice(h)
    print(f"    h={float(h):5.2f}:  I={mp.nstr(Ih,14)}  err={mp.nstr(abs(Ih-z3),3)}")

# --- 3. Richardson rate check: O(h^2) trapezoidal => error ratio -> 4 ---
print("\n[3] Richardson convergence rate (ratio -> 4 confirms O(h^2)):")
prev = None
for h in [mp.mpf('0.2'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('0.025'), mp.mpf('0.0125')]:
    e = abs(lattice(h) - z3)
    r = "" if prev is None else "ratio=%.2f" % float(prev / e)
    print(f"    h={float(h):7.4f}  err={mp.nstr(e,3)}  {r}")
    prev = e

# --- 4. contrast: §5.4's COUNT is a rational, cannot limit to zeta(3) ---
S = mp.nsum(lambda k: (-1)**(k-1) * 24 * mp.mpf(1)/7**(k-1), [1, mp.inf])
print("\n[4] contrast — §5.4 COUNT: sum_k (-1)^{k-1} 24 (1/7)^{k-1} = %s (= 21, RATIONAL)"
      % mp.nstr(S, 8))
print("    |21 - zeta(3)| = %s  — a convergent rational series cannot reach the period."
      % mp.nstr(abs(S - z3), 4))

print("""
VERDICT: a discrete LATTICE iterated integral converges to zeta(3) (O(h^2),
ratio 4.00), with no trig/period inserted by hand — the period emerges from the
iterated dx/x structure. §5.4's combinatorial count discards that structure and
lands on a rational. So the a_e obstruction is SHORTCUT-SPECIFIC, not
substrate-fundamental. Proof-of-route only: the TCH F2(0) vertex integral
(build it, integrate it) remains the open computation that would actually yield a_e.
""")
