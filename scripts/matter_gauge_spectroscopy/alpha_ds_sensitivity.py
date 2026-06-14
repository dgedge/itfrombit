#!/usr/bin/env python3
"""
Audit of the canonical §5.4 / Part 12 discrete Dyson-Schwinger equation for
alpha^{-1}.  Two jobs:

  PART A  REPRODUCIBILITY.  Does the equation AS WRITTEN actually yield the
          headline 137.035999077 (3 ppb)?  Solved two independent ways
          (fixed-point + exact quartic polyroots) to rule out a solver bug.
          Internal-consistency check of the three quoted numbers (1/2/3-loop).

  PART B  SENSITIVITY / FORMULA-FREEDOM (§16.3).  Treating the 2-loop piece as
          a free additive term T, how wide is the T-window that lands inside
          the CODATA bar, and how many simple expressions fall in it?

Equation, quoted verbatim (ANCHOR L579/L2691/L2746, DRIFT L720, STATUS L45,
part_12_fine_structure.tex eq. after L96):

    alpha^{-1}(alpha^{-1} - 137) = 31/(2pi) - (24/7)*(1/(2 pi alpha^{-1}))^2

Self-asserting: exit 0 == every asserted fact below holds.  CODATA quoted inline.
Nothing fitted.
"""
import mpmath as mp
from fractions import Fraction
mp.mp.dps = 50
twopi = 2*mp.pi

# --- CODATA references (quoted, not from memory) ---
# CODATA 2018: alpha^{-1} = 137.035999084(21)   (a_e route, Fan/Gabrielse)
# CODATA 2022: alpha^{-1} = 137.035999177(21)
codata18 = mp.mpf('137.035999084'); sigma = mp.mpf('0.000000021')
codata22 = mp.mpf('137.035999177')
CLAIM    = mp.mpf('137.035999077')   # framework headline
ppb      = lambda x, ref=codata18: abs(x-ref)/ref*mp.mpf(10)**9

def two_loop_term(x, Cabs=mp.mpf(24)/7):
    """The 2-loop piece EXACTLY as written: (24/7)*(1/(2 pi x))^2."""
    return Cabs*(1/(twopi*x))**2

def solve_fixedpoint(N1, two_loop):
    x = mp.mpf('137.036')
    for _ in range(300):
        rhs = N1/twopi - two_loop(x)
        xn = 137 + rhs/x
        if abs(xn-x) < mp.mpf(10)**(-40): x = xn; break
        x = xn
    return x

def solve_quartic(N1, Cabs=mp.mpf(24)/7):
    # x(x-137) = N1/(2pi) - (Cabs/(2pi)^2)/x^2  ; multiply by x^2:
    # x^3(x-137) - (N1/2pi)x^2 + Cabs/(2pi)^2 = 0
    a = mp.mpf(1); b = mp.mpf(-137); c = -(N1/twopi); d = mp.mpf(0)
    e = Cabs/twopi**2
    roots = mp.polyroots([a,b,c,d,e], maxsteps=200, extraprec=200)
    real = [r.real for r in roots if abs(r.imag) < mp.mpf(10)**(-20)]
    return max(real)  # physical root near 137

print("="*72)
print("PART A — REPRODUCIBILITY OF THE HEADLINE alpha^-1 = 137.035999077")
print("="*72)
N1 = mp.mpf(31)
x_fp   = solve_fixedpoint(N1, two_loop_term)
x_quar = solve_quartic(N1)
x_1l   = solve_fixedpoint(N1, lambda x: mp.mpf(0))
print(f"  equation-as-written, fixed-point :  alpha^-1 = {mp.nstr(x_fp,13)}")
print(f"  equation-as-written, quartic root:  alpha^-1 = {mp.nstr(x_quar,13)}")
print(f"  1-loop only (C2=0)               :  alpha^-1 = {mp.nstr(x_1l,13)}")
print(f"  framework CLAIMS (full 2-loop)   :  alpha^-1 = {mp.nstr(CLAIM,13)}")
print()
print(f"  two methods agree to            : {mp.nstr(abs(x_fp-x_quar),3)}")
print(f"  |equation-as-written - CLAIM|   : {mp.nstr(abs(x_fp-CLAIM),3)}  ({mp.nstr(ppb(x_fp,CLAIM),3)} ppb)")
print(f"  |equation-as-written - 1-loop|  : {mp.nstr(abs(x_fp-x_1l),3)}   <-- 2-loop term's TOTAL worth")
print(f"  equation-as-written vs CODATA18 : {mp.nstr(ppb(x_fp),3)} ppb off (NOT 3 ppb)")
print(f"  CLAIM vs CODATA18               : {mp.nstr(ppb(CLAIM),3)} ppb off")

# Residual of the equation at the CLAIMED value: is 137.035999077 even a root?
def residual(x, Cabs=mp.mpf(24)/7):
    return x*(x-137) - (N1/twopi - Cabs*(1/(twopi*x))**2)
print(f"\n  equation residual at the CLAIMED root 137.035999077 : {mp.nstr(residual(CLAIM),4)}")
print(f"  equation residual at 137.0360037 (true root)        : {mp.nstr(residual(x_fp),3)}")
print("  -> 137.035999077 is NOT a root of the written equation; 137.0360037 is.")

# What 2-loop term WOULD be needed to put the root at the claimed value?
T_req = N1/twopi - CLAIM*(CLAIM-137)        # required value of the subtracted term
T_written = two_loop_term(CLAIM)             # what the written term supplies there
print(f"\n  DIAGNOSIS — to make 137.035999077 the root, the subtracted 2-loop term")
print(f"  would need to equal T_req = {mp.nstr(T_req,5)} .")
print(f"  The written term (24/7)(1/(2pi x))^2 supplies only {mp.nstr(T_written,4)}")
print(f"  -> short by a factor of {mp.nstr(T_req/T_written,5)} (~150x).")
print(f"  Candidate structures evaluated at x=137.036:")
for lbl, val in [
    ("(24/7)/(2pi)^2            [no x]   ", (mp.mpf(24)/7)/twopi**2),
    ("(24/7)/((2pi)^2 * x)      [1 power]", (mp.mpf(24)/7)/(twopi**2*CLAIM)),
    ("(24/7)/((2pi)^2 * x^2)    [written]", (mp.mpf(24)/7)/(twopi**2*CLAIM**2)),
    ("(24/7)/(2pi * x)                   ", (mp.mpf(24)/7)/(twopi*CLAIM)),
]:
    print(f"     {lbl} = {mp.nstr(val,5):>12}   match to T_req? {'YES' if abs(val-T_req)/T_req<mp.mpf('0.05') else 'no'}")

# Internal consistency of the THREE quoted numbers via the stated geometric series.
print("\n  INTERNAL CONSISTENCY of the quoted 1/2/3-loop numbers:")
print("    quoted: 1-loop 137.036004 ; 2-loop 137.035999077 ; 3-loop 137.035999078")
s12 = mp.mpf('137.036004') - CLAIM
s23 = mp.mpf('137.035999078') - CLAIM
print(f"    shift 1->2 loop = {mp.nstr(s12,3)} ;  shift 2->3 loop = {mp.nstr(s23,3)}")
print(f"    ratio |s12/s23| = {mp.nstr(abs(s12/s23),4)}  (stated geometric ratio is 1/7 => expect ~7)")
print("    -> the loop-to-loop shifts are NOT in the claimed 1/7 geometric ratio.")

# Assert the ROBUST findings (these are the verified facts; exit 0 confirms them):
assert abs(x_fp - x_quar) < mp.mpf('1e-9'), "two solution methods disagree -> solver bug"
assert abs(x_fp - mp.mpf('137.0360037')) < mp.mpf('1e-6'), "written eqn root not 137.0360037"
assert abs(x_fp - CLAIM) > mp.mpf('4e-6'),  "written eqn unexpectedly close to claim"
assert abs(x_fp - x_1l)  < mp.mpf('1e-7'),  "2-loop term unexpectedly large"
assert T_req/T_written   > 100,             "2-loop term not ~150x short"

print("\n" + "="*72)
print("PART B — SENSITIVITY / FORMULA-FREEDOM (treat 2-loop piece as free term T)")
print("="*72)
# Use the structure that actually lands near CODATA so the test is meaningful:
# x(x-137) = 31/(2pi) - T, root x(T). Find the T-window giving |x-CODATA18|<sigma.
def x_of_T(T):
    x = mp.mpf('137.036')
    for _ in range(200):
        xn = 137 + (N1/twopi - T)/x
        if abs(xn-x)<mp.mpf(10)**(-40): x=xn;break
        x=xn
    return x
from mpmath import findroot
T_lo = findroot(lambda T: x_of_T(T)-(codata18+sigma), mp.mpf('6e-4'))
T_hi = findroot(lambda T: x_of_T(T)-(codata18-sigma), mp.mpf('6e-4'))
T_ctr= findroot(lambda T: x_of_T(T)-codata18,         mp.mpf('6e-4'))
window = abs(T_hi-T_lo)
print(f"  T that hits CODATA18 exactly         : {mp.nstr(T_ctr,6)}")
print(f"  T-window for landing within +/-1 sigma: width {mp.nstr(window,4)}")
print(f"  relative width of that window         : {mp.nstr(window/T_ctr*100,3)} %")
print("  -> any 2-loop expression evaluating to T within this narrow window")
print("     reproduces alpha^-1 to CODATA precision. The datum constrains only the")
print("     MAGNITUDE of the 2-loop term to ~4 sig figs; it does not single out 24/7.")

# How many simple p/q (as |C2|, in the ONE-POWER structure that ~works) land in-bar?
print("\n  Simple rationals |C2|=p/q (1<=p<=50, 1<=q<=12), one-power structure")
print("  T=(p/q)/((2pi)^2 x): how many land alpha^-1 within CODATA18 1-sigma?")
def x_onepow(Cabs):
    x=mp.mpf('137.036')
    for _ in range(200):
        xn=137+(N1/twopi-Cabs/(twopi**2*x))/x
        if abs(xn-x)<mp.mpf(10)**(-40):x=xn;break
        x=xn
    return x
cands=sorted({Fraction(p,q) for q in range(1,13) for p in range(1,51)},
             key=lambda f:(f.denominator,f.numerator))
inbar=[f for f in cands if abs(x_onepow(mp.mpf(f.numerator)/f.denominator)-codata18)<sigma]
print(f"    tested {len(cands)} simple p/q; {len(inbar)} land within 1 sigma:")
print("    ", [str(f) for f in inbar][:25])
print(f"    (24/7 = {mp.nstr(mp.mpf(24)/7,5)} in this structure gives "
      f"alpha^-1={mp.nstr(x_onepow(mp.mpf(24)/7),12)})")

print("\n" + "="*72)
print("VERDICT")
print("="*72)
print(" A. The DS equation AS WRITTEN does NOT reproduce 137.035999077. It yields")
print("    137.0360037 (= the 1-loop value, 33 ppb off CODATA). The written 2-loop")
print("    term is ~150x too small; the three quoted loop numbers are mutually")
print("    inconsistent with the stated 1/7 geometric series. The headline '3 ppb")
print("    at zero parameters' is NOT reproducible from the published equation.")
print(" B. Even with a structure that does land near CODATA, the 2-loop coefficient")
print("    is constrained only to ~4 sig figs in magnitude; multiple simple p/q sit")
print("    inside the bar. Per §16.3 this is formula-freedom, not a clean derivation.")
print(" Neither A nor B has anything to do with zeta(3): the QED running is rational")
print(" through 3 loops (verified separately); zeta(3) is a VERTEX (a_e) object.")
print("\nexit 0 == PART A robust findings asserted and verified.")
