#!/usr/bin/env python3
"""
decuplet_curvature_dgg.py
=========================
Reproducibility + SIGN-LOCK NO-GO for ANCHOR §15 item 141 (decuplet curvature residual).

Item 141 anchors: the empirical J=3/2 decuplet mass steps DECREASE with strangeness
(152.6 -> 148.8 -> 139 MeV, curvature ~ -14 MeV), while DGG static one-gluon-exchange
hyperfine predicts the OPPOSITE (increasing steps, curvature ~ +10 MeV). This script
verifies that and proves the stronger statement the prose only asserted:

  the static-hyperfine curvature  step3 - step1  =  (A/2)(1/m_u - 1/m_s)^2
  is a PERFECT SQUARE -> >= 0 for ANY positive (A, m_u, m_s).

So the negative empirical curvature is NOT a parameter-tuning miss; it is a SIGN
OBSTRUCTION -- unreachable by any static-hyperfine parameters. The ~14-MeV residual
is therefore genuine second-order physics (configuration mixing / relativistic), and
the framework's equal-spacing (zero-curvature) leading order is the better description.

Honest scope: <S_i.S_j> = +1/4 is pure J=3/2 angular-momentum algebra (NOT framework-
specific -- any spin-1/2 triple in the symmetric state gives it); (m_u,m_s,A) are DGG-
fitted constituent scales (item 142: in the baryon sector the framework IS the CQM).
This script confirms the framework INHERITS DGG's wrong-sign curvature -- a consistency
check, not new physics.

numpy only. exit 0 == every quoted number verified.
"""
import numpy as np

fails = []
def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond: fails.append(msg)

# ---- 1. <S_i.S_j> = +1/4 for the decuplet (J=3/2), two independent routes ----
sx = 0.5*np.array([[0,1],[1,0]], complex)
sy = 0.5*np.array([[0,-1j],[1j,0]])
sz = 0.5*np.array([[1,0],[0,-1]], complex)
def op(s, i):                                   # embed 1-spin op on factor i of 3
    o = [np.eye(2, dtype=complex)]*3; o[i] = s
    return np.kron(np.kron(o[0], o[1]), o[2])
def SiSj(i, j):
    return op(sx,i)@op(sx,j) + op(sy,i)@op(sy,j) + op(sz,i)@op(sz,j)

print("[1] decuplet spin algebra")
psi = np.zeros(8, complex); psi[0] = 1.0                     # |uuu>, Sz=3/2 (stretched)
W   = np.zeros(8, complex); W[1]=W[2]=W[4] = 1/np.sqrt(3)    # symmetric Sz=1/2 member
for st, nm in [(psi, "Sz=3/2 stretched"), (W, "Sz=1/2 symmetric")]:
    vals = [float(np.real(st.conj()@SiSj(i,j)@st)) for i,j in [(0,1),(0,2),(1,2)]]
    check(all(abs(v-0.25) < 1e-12 for v in vals),
          f"<S_i.S_j> = +1/4 on all pairs ({nm}): {[round(v,3) for v in vals]}")
Stot = [sum(op(s,i) for i in range(3)) for s in (sx,sy,sz)]
S2 = sum(S@S for S in Stot)
check(abs(np.real(psi.conj()@S2@psi) - 15/4) < 1e-12, "S^2 = 15/4 (J=3/2) cross-check")

# ---- 2. DGG decuplet masses with the anchored constituent scales ----
mu, ms = 362.0, 535.0
Amu2 = 146.5/0.75                              # A/m_u^2 fixed by Delta hyperfine = (3/4)(A/m_u^2) = kappa
A = Amu2*mu*mu
print(f"\n[2] DGG static hyperfine  (m_u={mu:.0f}, m_s={ms:.0f} MeV, A/m_u^2={Amu2:.1f} MeV = kappa)")
def mass(ns):
    m = [mu]*(3-ns) + [ms]*ns
    return sum(m) + sum((A*0.25)/(m[i]*m[j]) for i,j in [(0,1),(0,2),(1,2)])
names = ["Delta ", "Sigma*", "Xi*   ", "Omega "]
M = [mass(n) for n in range(4)]
for nm, m in zip(names, M): print(f"      {nm} = {m:8.2f} MeV")
steps = [M[i+1]-M[i] for i in range(3)]
curv = steps[2] - steps[0]
print(f"      steps = {steps[0]:.2f} -> {steps[1]:.2f} -> {steps[2]:.2f} MeV ;  curvature = {curv:+.2f} MeV")
check(abs(M[0]-1232.47) < 0.05 and abs(M[1]-1373.90) < 0.05 and
      abs(M[2]-1520.43) < 0.05 and abs(M[3]-1672.06) < 0.05, "absolute masses reproduce the writeup")
check(abs(steps[0]-141.42) < 0.05 and abs(steps[1]-146.53) < 0.05 and abs(steps[2]-151.63) < 0.05,
      "steps reproduce 141.42 -> 146.53 -> 151.63")
check(abs(curv-10.21) < 0.05, "theoretical curvature = +10.21 MeV")

# ---- 3. SIGN-LOCK NO-GO: curvature = (A/2)(1/m_u - 1/m_s)^2  >= 0 always ----
print("\n[3] sign-lock no-go")
curv_closed = (A/2)*(1/mu - 1/ms)**2
check(abs(curv - curv_closed) < 1e-6,
      f"curvature == (A/2)(1/m_u - 1/m_s)^2 = {curv_closed:+.4f} MeV  (closed form)")
# deterministic grid scan: every positive (A, m_u, m_s) gives curvature >= 0
worst = min((np.float64(a)/2)*(1/u - 1/s)**2
            for a in (1e5, 1e6, 5e6, 2e7)
            for u in (250., 300., 362., 450.)
            for s in (400., 535., 700., 250.))     # includes m_s < m_u
check(worst >= 0.0, f"grid scan over (A,m_u,m_s): min curvature = {worst:.4f} >= 0 (sign-locked)")
print("      => the empirical NEGATIVE curvature is unreachable by any static-hyperfine params:")
print("         a sign obstruction (second-order physics needed), not a parameter-tuning gap.")

# ---- 4. empirical PDG + equal-spacing beats DGG ----
print("\n[4] empirical (PDG, isospin-averaged) vs leading-order models")
pdg = [1232.0, 1384.6, 1533.4, 1672.45]
es = [pdg[i+1]-pdg[i] for i in range(3)]
ecurv = es[2]-es[0]
print(f"      PDG steps = {es[0]:.1f} -> {es[1]:.1f} -> {es[2]:.1f} MeV ;  empirical curvature = {ecurv:+.1f} MeV")
res_equal = abs(ecurv - 0.0)        # equal-spacing model: zero curvature
res_dgg   = abs(ecurv - curv)       # DGG model: +10.21 curvature
print(f"      |empirical - equal-spacing(0)|   = {res_equal:.1f} MeV")
print(f"      |empirical - DGG(+{curv:.1f})|        = {res_dgg:.1f} MeV")
check(ecurv < 0, "empirical decuplet curvature is NEGATIVE")
check(res_equal < res_dgg, "equal-spacing (zero curvature) beats DGG on the curvature direction")

print("\n" + "="*70)
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]
    import sys; sys.exit(1)
print("RESULT: exit 0 -- decuplet curvature verified; static DGG hyperfine is")
print("sign-locked to POSITIVE curvature (perfect square), so the empirical")
print("negative curvature is a second-order-physics residual, confirming §15 item 141.")
print("="*70)
