#!/usr/bin/env python3
"""
SU(3) baryon mass-relation test of the framework's baryon machinery.
Consolidated 2026-05-30; CORRECTED 2026-05-30 after reading §15 item 129.

The framework spreads its baryon-mass content across THREE complementary
mechanisms, not one. An honest test must hit each on its own scale:
  §9.10 boxed w-formula  -> isospin/EM splittings (few MeV)
  §15 item 129 (C_8 Nyquist fractions) -> octet STRANGENESS ladder (tens-hundreds MeV)
  §15 item 130 (O_h irrep dim D)       -> J^P octet-vs-decuplet split (~293 MeV, N->Delta)

Framework inputs (verbatim from canon):
  §2.8 quark charges : u=+2/3, d=-1/3, s=-1/3 ; §5.9 active-bit I3 : u=0, d=1, s=1
  §9.10 : dm = w[ d(sum I3) - (7/8)d(Qtot^2) + 4 d(sum_{i<j} q_iq_j) ], w=alpha*Lambda=2.4227 MeV
  §129  : M = 939 + (k/8) Lambda_QCD ; Lambda_QCD=332 ; k = 4(Lam,I=0), 6(Sig,I=1), 9(Xi)
  §130  : M = M0 + kappa*D ; D = O_h irrep dim (octet 2, decuplet 4)

FINDINGS (corrected):
  PASS  Coleman-Glashow (p-n)+(Xi0-Xi-)-(Sig+-Sig-) = 0 identically  [Tier-A; now canon item 140]
  PASS  n-p, Xi, Sigma span via §9.10  at <0.1 MeV
  PASS  octet strangeness centres Lam/Sig/Xi via §129 at 0.2-0.9%
  PARTIAL Sigma-Lambda IS captured by §129 (6/8 vs 4/8 fractions at same n_s) = 83 vs 77 MeV (~7%).
        This is the I(I+1) content -- present in canon as an isospin-keyed fraction
        (Lambda_QCD/8 = 41.5 MeV per Casimir unit), hand-assigned (motivated, not forced).
        [Earlier consolidation WRONGLY reported "framework gives Sigma-Lambda=0, FAIL" -- that
         tested only the §9.10 isospin formula and ignored §129. Corrected here.]
  WEAK  Sigma0 EM curvature (a pure isospin quantity, §9.10's own scale): 3.5x too large;
        §9.10's "stepwise" paragraph EM = (exp - w) is a retrofit.
  FIT   §130 is a 2-point fit (N,Delta); within-octet J^P degenerate by design (strangeness
        lives in §129, not §130) -- so "all octet equal" under §130 is correct, not a failure.

  GENUINE GAPS (the user's two points, both verified against canon):
   (1) DECUPLET spectrum: §129 stops at the octet; Sigma*,Xi*,Omega are NOT anchored. This
       is the real over-simplification. A per-n_s increment recovers equal-spacing (Omega- ~1%);
       the increment's SCALE P_s = strange constituent mass stays Target-B-blocked.
   (2) PRESENTATION: three mechanisms (§129 fractions / §130 irrep-dim / §9.10 EM) carry content
       a single GMO-style M = a + bY + c[I(I+1)-Y^2/4] expresses in one formula. Cosmetic, not structural.

Needs only mpmath + numpy. Nothing fitted except where labelled.
"""
from fractions import Fraction as F
import numpy as np
import mpmath as mp

w = F(24227, 10000)
Q = {'u': F(2,3), 'd': F(-1,3), 's': F(-1,3)}
I3 = {'u': 0, 'd': 1, 's': 1}
C = {'p':'uud','n':'udd','Sp':'uus','S0':'uds','Sm':'dds','Xi0':'uss','Xim':'dss','Lam':'uds'}
PDG = {'p':938.272,'n':939.565,'Sp':1189.37,'S0':1192.64,'Sm':1197.45,
       'Xi0':1314.86,'Xim':1321.71,'Lam':1115.68}
def Qt(b):  return sum(Q[q] for q in C[b])
def SI(b):  return sum(I3[q] for q in C[b])
def Sqq(b):
    q=[Q[x] for x in C[b]]; return q[0]*q[1]+q[0]*q[2]+q[1]*q[2]
def dm(b1,b2):
    return float(w*((SI(b1)-SI(b2)) - F(7,8)*(Qt(b1)**2-Qt(b2)**2) + 4*(Sqq(b1)-Sqq(b2))))
def H(s): print("\n"+"="*70+"\n"+s+"\n"+"="*70)

H("1. Isospin/EM splittings - §9.10 boxed w-formula vs PDG (MeV)")
for b1,b2 in [('n','p'),('Sm','Sp'),('Sm','S0'),('Xim','Xi0')]:
    f=dm(b1,b2); e=PDG[b1]-PDG[b2]
    print(f"  {b1}-{b2:3s}: framework={f:+.3f}  exp={e:+.3f}  dev={abs(f-e):.2f}")

H("2. Coleman-Glashow: (p-n)+(Xi0-Xi-)-(Sig+ -Sig-) = 0 ?  [canon item 140]")
def cg(fn): return (fn('p')-fn('n'))+(fn('Xi0')-fn('Xim'))-(fn('Sp')-fn('Sm'))
print(f"  sum I3 term = {cg(SI)} ; Qtot^2 term = {cg(lambda b: Qt(b)**2)} ; sum q_iq_j = {cg(Sqq)}")
full=dm('p','n')+dm('Xi0','Xim')-dm('Sp','Sm')
print(f"  => full CG = {full:+.6f} MeV ; identically zero: {abs(full)<1e-9}   [PASS, Tier-A]")

H("3. Sigma0 EM curvature - §9.10's OWN (isospin) scale; stepwise paragraph is a retrofit")
s1,s2=dm('S0','Sp'),dm('Sm','S0'); e1,e2=PDG['S0']-PDG['Sp'],PDG['Sm']-PDG['S0']
print(f"  step Sig0-Sig+: formula EM={s1-2.4227:+.2f}  paragraph EM=+0.85  (=exp-w={e1-2.4227:+.2f})")
print(f"  step Sig- -Sig0: formula EM={s2-2.4227:+.2f}  paragraph EM=+2.39  (=exp-w={e2-2.4227:+.2f})")
print("  => paragraph EM = (experiment - w): a RETROFIT, not the formula's EM.")
curv_f=dm('Sm','Sp')-2*dm('S0','Sp'); curv_e=PDG['Sp']+PDG['Sm']-2*PDG['S0']
print(f"  frame-indep curvature Sig+ +Sig- -2Sig0: formula={curv_f:+.2f} PDG={curv_e:+.2f}"
      f" ({curv_f/curv_e:.1f}x)  [WEAK -- §9.10 EM only; the strangeness sector is §129's job]")

H("4. Octet STRANGENESS ladder - §15 item 129 (C_8 Nyquist fractions). Captures Sigma-Lambda.")
L=332.0; MN=939.0
lad={'Lam':(F(4,8),1,PDG['Lam']),'Sig':(F(6,8),1,(PDG['Sp']+PDG['S0']+PDG['Sm'])/3),
     'Xi':(F(9,8),2,(PDG['Xi0']+PDG['Xim'])/2)}
for b,(k,ns,pdg) in lad.items():
    M=MN+float(k)*L
    print(f"  {b:4s}(I={'0' if b=='Lam' else ('1' if b=='Sig' else '1/2')},n_s={ns}): "
          f"{int(k*8)}/8 -> {M:.1f}  PDG={pdg:.1f}  resid={M-pdg:+.1f} ({(M-pdg)/pdg*100:+.2f}%)")
print("  (all residuals same sign => small common systematic in the 939/332 base, not random scatter)")
SmL=(float(F(6,8))-float(F(4,8)))*L; SmL_e=PDG['S0']-PDG['Lam']
print(f"  Sigma-Lambda = (6/8 - 4/8)*Lambda_QCD = {SmL:.1f}  vs PDG {SmL_e:.2f}  ({(SmL-SmL_e)/SmL_e*100:+.1f}%)")
print(f"  6/8(I=1) vs 4/8(I=0) at the SAME n_s = an isospin-keyed mass term: functionally the")
print(f"  I(I+1) Casimir with coeff Lambda_QCD/8 = {L/8:.1f} MeV/unit (empirical {SmL_e/2:.1f}).")
print(f"  => the I(I+1) structure IS in canon -- as a hand-assigned fraction (k=4 logical dim;")
print(f"     3 axes x 2 bits = 6; 3^2 cross-term = 9), motivated but not forced. Postdiction, ~7% on the split.")

H("5. §130 (O_h irrep dim) and §129 (strangeness) are COMPLEMENTARY axes, not rivals")
M0,kap=np.linalg.solve(np.array([[1,2],[1,4]]),np.array([939.0,1232.0]))
print(f"  §130 J^P split: M=M0+kappa*D from (N,Delta) -> M0={M0:.0f}, kappa={kap:.1f}")
print(f"    2 eqns / 2 unknowns = exact fit; defends via kappa~M_pi(140), M0~2*Lambda_QCD(664).")
print(f"  Within-octet, all baryons share D=2 -> §130 gives them equal. CORRECT: the strangeness")
print(f"    splitting is §129's axis (4/8,6/8,9/8), not §130's. 'all octet equal' is not a failure.")
print(f"  N->Delta (J^P, §130, {1232-939:.0f} MeV) _|_ N->Lam->Sig->Xi (strangeness, §129).")

H("6. GENUINE GAPS: (1) decuplet ladder unanchored  (2) three mechanisms vs one GMO formula")
dec={'Delta':1232.0,'Sigma*':1384.6,'Xi*':1531.8,'Omega':1672.45}; nsd={'Delta':0,'Sigma*':1,'Xi*':2,'Omega':3}
Ps=dec['Sigma*']-dec['Delta']
print(f"  (1) §129 stops at the OCTET. Decuplet Sigma*,Xi*,Omega NOT anchored. Per-n_s (Ps={Ps:.1f}")
print(f"      from Delta,Sigma* only) recovers the rest:")
for bn in ['Xi*','Omega']:
    pr=dec['Delta']+nsd[bn]*Ps
    print(f"        {bn:6s}: pred={pr:.1f} PDG={dec[bn]:.1f} dev={pr-dec[bn]:+.1f} ({(pr-dec[bn])/dec[bn]*100:+.2f}%)")
print(f"      Ps SCALE = strange constituent mass = Target-B-blocked (decuplet step {Ps:.0f} vs octet N->Lam step 166;")
mp.mp.dps=20; Lm=mp.mpf('332'); phi=(1+mp.sqrt(5))/2
fac=[(mp.mpf(1),''),(mp.sqrt(2),'r2'),(1/mp.sqrt(2),'/r2'),(phi,'phi'),(1/phi,'/phi')]
hits=sum(1 for p in range(1,13) for q in range(1,13) for f,_ in fac if 140<=Lm*mp.mpf(p)/q*f<=190)
print(f"      neither derived; {hits} simple Lambda-ratios crowd [140,190] MeV so a 'match' is unearned, S16.3).")
oct_={'N':939.0,'Lambda':1115.7,'Sigma':1193.2,'Xi':1318.3}; qn={'N':(1,.5),'Lambda':(0,0),'Sigma':(0,1),'Xi':(-1,.5)}
M=[[1,Y,I*(I+1)-Y*Y/4] for (Y,I) in (qn[b] for b in ['N','Lambda','Sigma','Xi'])]
yv=[oct_[b] for b in ['N','Lambda','Sigma','Xi']]
coef,_,_,_=np.linalg.lstsq(np.array(M),np.array(yv),rcond=None)
gl=2*(oct_['N']+oct_['Xi']); gr=3*oct_['Lambda']+oct_['Sigma']
print(f"  (2) One GMO M=a+bY+c[I(I+1)-Y^2/4] reproduces the same octet content §129 splits across")
print(f"      fractions: GMO 2(N+Xi)={gl:.1f} vs 3Lam+Sig={gr:.1f} (dev {abs(gl-gr)/gl*100:.2f}%).")
print(f"      Unifying §129/§130/§9.10 under one Casimir formula is cleaner -- presentational, not a bug.")

H("NET")
print("""Corrected from the earlier consolidation. The framework's baryon content is SPREAD
across §9.10 (isospin/EM + Coleman-Glashow identity, now item 140), §129 (octet
strangeness ladder incl. Sigma-Lambda at ~7%, centres <1%), and §130 (N/Delta J^P).
Tested per-mechanism, it captures more than a one-formula probe suggests -- the
within-octet structure is NOT missing (my earlier "Sigma-Lambda=0 FAIL" was wrong:
it ignored §129). Two genuine gaps remain: (1) the DECUPLET strangeness ladder
(Sigma*,Xi*,Omega) is unanchored -- the real over-simplification; a per-n_s increment
recovers it but its SCALE P_s stays Target-B-blocked; (2) three mechanisms carry what
one GMO-style Casimir formula would -- presentational, not structural.""")
