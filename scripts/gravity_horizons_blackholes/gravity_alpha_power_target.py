#!/usr/bin/env python3
r"""The one live gravity thread: can the alpha-power that dresses the Lambda<->horizon Dirac bridge be
DERIVED (=> a derived Dirac-type RELATION among {M_P, Lambda_QCD, horizon}, NOT an intrinsic M_P)?
Honest scoping (DRIFT G7 follow-up). NO expression-shopping: we put all three routes in ONE canonical
Sakharov form and report what the data actually determines and what stays free.

Canonical form (§10.5):  M_P^2 = M_bare^2 * (R/a0)^s / K ,  M_bare ~ Lambda, a0 = 1/Lambda.
  s = horizon-power, K = dimensionless prefactor (the thing that carries any alpha-power).
"""
import sys, math
from fractions import Fraction
import numpy as np

M_P=1.220890e19; LAM=0.332; ALPHA=1/137.035999
H0=(67.36e3/3.0856776e22)*6.582119e-25; OmL=0.6847
R_H=1/H0; a0=1/LAM
# de Sitter horizon (§10.5):  R_dS = sqrt(3/Lambda_cosmo), Lambda_cosmo = 3 H0^2 OmL  ->  R_dS = R_H/sqrt(OmL)
R_dS = R_H/math.sqrt(OmL)
print(f"de Sitter relation: R_dS = R_H/sqrt(OmL) = {R_dS:.3e} GeV^-1 ; R_dS/a0 = {R_dS/a0:.3e} (§10.5 quotes 2.77e41)")
assert abs(R_dS/a0/2.77e41 - 1) < 0.02

# ---------- (1) UNIFICATION: the 205-route and Part-20 are the SAME alpha^2 Sakharov relation ----------
# 205-route: M_P^2 = Lambda^2 (R_dS/a0)/K  ->  K from data:
K_data = LAM**2 * (R_dS/a0) / M_P**2
# Part-20 re-expressed in the same form (set M_bare=Lambda, horizon-power s=1):  K = sqrt(OmL)/(24 pi alpha^2)
K_part20 = math.sqrt(OmL)/(24*math.pi*ALPHA**2)
print(f"\n[1] UNIFICATION (both with M_bare=Lambda, horizon-power s=1):")
print(f"    K from the 205-route data           = {K_data:.2f}   (canon asserts the integer 205 = '208-3')")
print(f"    K from Part-20  = sqrt(OmL)/(24pi a^2) = {K_part20:.2f}")
print(f"    => these agree to {abs(K_data-K_part20)/K_data*100:.2f}% : the 205-route and Part-20 are ONE alpha^2 relation,")
print(f"       differing only by the de Sitter bookkeeping R_dS = R_H/sqrt(OmL). '208-3 rank' is a numerical")
print(f"       dressing on sqrt(OmL)/(24pi alpha^2) = {K_part20:.1f}. So 2 of the 3 routes are NOT independent.")
assert abs(K_data-K_part20)/K_data < 0.01

# ---------- (2) horizon-power s is DERIVABLE (=1/2 on M_P); item 136 (s->1/4) is the outlier ----------
# Sakharov codim-1 dilution (§10.5): spin-2 transverse mode on 3D lattice dilutes linearly with 1D radial
# extent -> M_P^2 ~ (R/a0)^1 -> M_P ~ (R/a0)^{1/2}.  i.e. horizon-power on M_P = 1/2.
s_136_on_MP = 1/4   # item 136 uses (L_H/a0)^{1/4} on M_P
print(f"\n[2] horizon-power (Sakharov codim-1 dilution, §10.5): M_P^2 ~ (R/a0)^1 -> s=1/2 on M_P. DERIVABLE.")
print(f"    205-route & Part-20 both have s=1/2 on M_P (M_P^2 ~ R^1). Item 136 has s={s_136_on_MP} on M_P (M_P~R^1/4)")
print(f"    -> item 136's horizon-power is the OUTLIER, ruled out by the codim-1 geometry (and its 1/4 was")
print(f"       magnitude-selected, item136 audit). So the live relation is the s=1/2 alpha^? Sakharov one.")

# ---------- (3) is the alpha-POWER forced? §16.3 scan of the prefactor 1/K (NOT shopping: count them all) --
P = 1/K_data
print(f"\n[3] is the alpha-power FORCED?  prefactor 1/K = {P:.4e}. Equally-simple readings within 1%:")
cands=[]
for q in range(-1,4):                                   # alpha^q
    for c in [Fraction(n,d) for d in range(1,9) for n in range(1,13)]:
        for extra,elab in [(1.0,""),(math.pi,"*pi"),(1/math.sqrt(OmL),"/sqrt(OmL)"),(math.sqrt(OmL),"*sqrt(OmL)")]:
            v=float(c)*extra*ALPHA**q
            if v>0 and abs(v-P)/P<=0.01:
                cands.append((abs(v-P)/P, f"{c}{elab}*alpha^{q}", v, q))
cands.sort()
seenq=set()
for err,lab,v,q in cands[:8]:
    print(f"     {lab:<22} = {v:.4e}  ({err*100:.2f}%)   [alpha-power {q:+d}]")
    seenq.add(q)
# also the pure integer reading:
print(f"     integer 1/205          = {1/205:.4e}  ({abs(1/205-P)/P*100:.2f}%)   [alpha-power  0]")
print(f"  => the prefactor admits readings at alpha-powers {sorted(seenq|{0})} (plus item 136's -4 in a")
print(f"     different horizon-power) -- all within ~1%. The magnitude does NOT force the alpha-power; even")
print(f"     its SIGN is ambiguous across the framework's routes. So alpha^2 is NOT magnitude-determined.")

# ---------- (4) what would 'deriving alpha^2' need, and the bare-anchor tension ----------
K_Eg_needed = 1/(4*math.pi)          # M_bare=Lambda requires M_bare=sqrt(4pi K_Eg)Lambda=Lambda -> K_Eg=1/4pi
print(f"\n[4] to land the clean Part-20 coefficient 24pi, the bare anchor must be M_bare=Lambda EXACTLY,")
print(f"    i.e. K_Eg = 1/(4pi) = {K_Eg_needed:.4f}. The pre-pinned heat-kernel (item7_keg_heatkernel.py) gave")
print(f"    K_Eg = 0.44 at the canonical pin g_s=1 (M_bare=0.78 GeV=2.3 Lambda) -- ~5.5x off; K_Eg=1/4pi needs")
print(f"    |t_mix|~2.8 (tuned). So even the bare anchor is not robustly Lambda. And the alpha^2 itself is the")
print(f"    'double-Landauer' factor §10.7 calls 'a phenomenological ansatz, NOT a derived theorem': deriving it")
print(f"    = showing each colour-bridge snap in the Feshbach trace over Q carries exactly one Landauer w=alpha*Lambda")
print(f"    -- which IS the K_eff trace DRIFT G7 found NOT constructible. Same bottleneck.")

print(f"""
=========================================================================================
VERDICT (the derived-alpha-power -> derived-relation target, scoped honestly):
  REACHED (real results):
   * 205-route == Part-20: one alpha^2 Sakharov relation, K = sqrt(OmL)/(24pi alpha^2) = {K_part20:.0f} ~ 205. The
     three "independent" routes collapse to ONE alpha^2 relation + one outlier (item 136). '208-3 rank' is a
     numerical dressing on the alpha^2 expression (and G7 already showed the rank is not constructible).
   * horizon-power s=1/2 on M_P is DERIVABLE (Sakharov codim-1 dilution); item 136's 1/4 is the spurious outlier.
  NOT REACHED (the gap, unchanged):
   * the alpha-POWER is NOT magnitude-forced: the prefactor admits alpha^0 (integer 205), alpha^1, alpha^2, all
     within ~1% (and alpha^-4 in item 136's wrong horizon-power) -- the sign isn't even pinned. alpha^2 wins only
     on the 'double-Landauer' structural story, which §10.7 itself flags 'not a theorem'.
   * deriving alpha^2 bottlenecks on the SAME non-constructible Feshbach-trace-over-Q as K_eff=205 (G7); and the
     bare anchor M_bare=Lambda (needed for 24pi) is itself a t_mix tuning, not a robust heat-kernel output.
  NET: the target's best attainable prize is a DERIVED DIRAC RELATION  M_P^2 = 24pi alpha^2 Lambda^3/(H0 OmL)
  (= Part-20) with the horizon as INPUT -- and it remains blocked at the alpha^2, which is (a) not forced by the
  magnitude and (b) gated on the same unconstructed trace as 205. So the honest status is unchanged and now
  sharper: the gravity sector is ONE alpha^2 Sakharov/Dirac relation whose horizon-power is derivable, whose
  alpha-power is not, and whose alpha-power derivation requires the one object the framework cannot yet build.
  A full W_QQ(k) build is worth it ONLY to attack that alpha^2-via-Feshbach-trace -- nothing else moves.
=========================================================================================""")
assert isinstance(K_data,float)
print("exit 0 -- alpha-power target scoped: 205==Part-20 (alpha^2); horizon-power derivable; alpha-power not; same G7 bottleneck.")
