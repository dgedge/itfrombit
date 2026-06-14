#!/usr/bin/env python3
"""
Deeper look at the J^P=3/2+ baryon DECUPLET ladder -- the gap flagged in
su3_baryon_test.py probe 6(1). EXPLORATORY (uncommitted).

Question: can the framework's ALREADY-ANCHORED constants reach the decuplet
spectrum (Delta, Sigma*, Xi*, Omega-) without a new free parameter?

Canon pieces in play:
  item 129 : octet strangeness ladder in units of Lambda_QCD/8; N->Lam step = 4/8 = 166 MeV
  item 130 : M = M0 + kappa*D, kappa = 146.5 MeV ~ M_pi, fixed by (N, Delta) only
  PDG decuplet (isospin-averaged): Delta 1232.0, Sigma* 1384.6, Xi* 1533.4, Omega 1672.45

This script REVISES the earlier su3_baryon_test verdict ("P_s Target-B-blocked, a free
strange mass, Omega off ~1%"). That used a 2-point free fit P_s=152.6 from (Delta,Sigma*).
Here we test the better hypothesis: the decuplet strange step IS item-130's kappa.

Nothing fitted except the explicit GMO cross-check at the end (4 params / 8 baryons).
"""
import numpy as np

# ---- PDG isospin-averaged masses (MeV) ----
DEC = {'Delta':1232.0,'Sigma*':1384.6,'Xi*':1533.4,'Omega':1672.45}
NS  = {'Delta':0,'Sigma*':1,'Xi*':2,'Omega':3}     # number of strange quarks
I   = {'Delta':1.5,'Sigma*':1.0,'Xi*':0.5,'Omega':0.0}  # isospin -> note I=(3-n_s)/2
Ybar= {'Delta':1,'Sigma*':0,'Xi*':-1,'Omega':-2}   # hypercharge = 1 - n_s
OCT = {'N':939.0,'Lambda':1115.7,'Sigma':1193.2,'Xi':1318.3}
L   = 332.0            # Lambda_QCD (canon §1.4)
KAP = 146.5            # item-130 kappa = (M_Delta - M_N)/2, fixed by (N,Delta) ONLY
MPI = 139.57           # charged pion (item-130 identifies kappa ~ M_pi)
def H(s): print("\n"+"="*72+"\n"+s+"\n"+"="*72)

H("1. The raw decuplet: is it equal-spaced? (Gell-Mann's Omega- prediction)")
ks=list(DEC); sp=[DEC[ks[i+1]]-DEC[ks[i]] for i in range(3)]
print(f"  masses : {[DEC[k] for k in ks]}")
print(f"  spacings: {[f'{s:.1f}' for s in sp]}  mean={sum(sp)/3:.2f}  (decreasing => slight curvature)")
print(f"  => equal-spacing holds to +-{(max(sp)-min(sp))/2:.0f} MeV; curvature ~{sp[0]-sp[2]:.0f} MeV over 3 steps")

H("2. WHY equal-spacing is structural here: I is slaved to n_s (I=(3-n_s)/2)")
print("   On the decuplet the GMO Casimir [I(I+1)-Y^2/4] is itself LINEAR in n_s:")
for k in ks:
    cas=I[k]*(I[k]+1)-Ybar[k]**2/4
    print(f"     {k:7s}: n_s={NS[k]} I={I[k]} Y={Ybar[k]:+d}  I(I+1)-Y^2/4 = {cas:+.2f}")
print("   3.50, 2.00, 0.50, -1.00 -> constant step -1.5 => any Casimir term collapses to")
print("   linear-in-n_s. So a single strange ladder REPRODUCES equal-spacing (unlike the")
print("   octet, where I is NOT slaved to n_s and 6/8-vs-4/8 contamination appears).")

H("3. KEY TEST: decuplet strange step = item-130 kappa? (zero new parameters)")
print(f"   item-130 kappa = {KAP} (fixed by N,Delta).  decuplet mean spacing = {(DEC['Omega']-DEC['Delta'])/3:.2f}")
print(f"   agreement: {abs(KAP-(DEC['Omega']-DEC['Delta'])/3)/KAP*100:.2f}%   (vs octet N->Lam step 166 -- see probe 5)")
print(f"   M(n_s) = M_Delta + kappa*n_s  [parameter-free; kappa NOT fit to decuplet]:")
for k in ks:
    pr=DEC['Delta']+KAP*NS[k]; print(f"     {k:7s}: pred={pr:7.1f}  PDG={DEC[k]:7.1f}  dev={pr-DEC[k]:+5.1f} ({(pr-DEC[k])/DEC[k]*100:+.2f}%)")
print("   vs the earlier 2-point free fit P_s=152.6 (Delta,Sigma*): Omega 1689.8 (+1.0%).")
print("   => item-130's kappa fits the decuplet BETTER than the free P_s did, and is NOT free.")

H("4. Unified form: M_decuplet = M0 + kappa*(D + n_s),  D=4")
M0=646.0
for k in ks:
    pr=M0+KAP*(4+NS[k]); print(f"     {k:7s}: M0+kappa*(4+{NS[k]}) = {pr:7.1f}  PDG={DEC[k]:7.1f}  dev={pr-DEC[k]:+5.1f}")
print(f"   J^P tension (D) and strange step (n_s) carried by the SAME kappa~M_pi -> one ladder.")

H("5. The octet step (166-177) != decuplet step (147): a spin-dependent signal")
print(f"   octet  N->Lambda  = {OCT['Lambda']-OCT['N']:.1f} MeV empirical  (item-129 4/8 = 166)")
print(f"   decuplet per-n_s  = {(DEC['Omega']-DEC['Delta'])/3:.1f} MeV  (= kappa)")
print(f"   difference {OCT['Lambda']-OCT['N']-(DEC['Omega']-DEC['Delta'])/3:.0f} MeV: DIRECTION matches QCD -- the octet strange step")
print("   carries spin-dependent (chromomagnetic) structure that the spin-aligned decuplet")
print("   does not, so the decuplet step is the cleaner constituent increment. Magnitude")
print("   attribution is qualitative, NOT a derived hyperfine value.")

H("6. Does one GMO formula M=a+b*Y+c*[I(I+1)-Y^2/4] cover BOTH multiplets? (4 params, 8 data)")
OCTET={'N','Lambda','Sigma','Xi'}
rows=[]; ys=[]
data=[('N',939.0,0.5,1),('Lambda',1115.7,0,0),('Sigma',1193.2,1,0),('Xi',1318.3,0.5,-1),
      ('Delta',1232.0,1.5,1),('Sigma*',1384.6,1.0,0),('Xi*',1533.4,0.5,-1),('Omega',1672.45,0,-2)]
for name,mass,ii,yy in data:
    D=2 if name in OCTET else 4
    rows.append([1,D,yy,ii*(ii+1)-yy*yy/4]); ys.append(mass)
Aarr=np.array(rows,float); yv=np.array(ys)
coef,res,_,_=np.linalg.lstsq(Aarr,yv,rcond=None); pred=Aarr@coef
print(f"   fit M = a + d*Dim + b*Y + c*[I(I+1)-Y^2/4]:  a={coef[0]:.1f} d={coef[1]:.1f} b={coef[2]:.1f} c={coef[3]:.1f}")
for (n,_,_,_),p,o in zip(data,pred,yv): print(f"     {n:7s}: pred={p:7.1f} PDG={o:7.1f} dev={p-o:+5.1f}")
print(f"   max|dev| over all 8 = {np.max(np.abs(pred-yv)):.1f} MeV.")
print(f"   CAVEAT: d (irrep-dim slope) = {coef[1]:.0f} != kappa={KAP}. Dim and the Casimir are")
print(f"   correlated across the two multiplets, so the global fit SPLITS the N->Delta gap")
print(f"   (2d+3c={2*coef[1]+3*coef[3]:.0f} ~ 293) between them -- d is NOT independently kappa.")
print(f"   So: one 4-param GMO covers all 8 baryons to <={np.max(np.abs(pred-yv)):.0f} MeV, but its parameters")
print(f"   do NOT map cleanly onto framework constants. The clean result is the n_s ladder (probe 3).")

H("7. Honest checks: is kappa~147 'unearned', and what is NOT derived")
import mpmath as mp; mp.mp.dps=20; phi=(1+mp.sqrt(5))/2
fac=[(mp.mpf(1),'1'),(mp.sqrt(2),'r2'),(1/mp.sqrt(2),'/r2'),(phi,'phi'),(1/phi,'/phi')]
hits=[(p,q,t) for p in range(1,13) for q in range(1,13) for f,t in fac if 140<=float(mp.mpf(p)/q*f*L)<=153]
print(f"   crowding: {len(hits)} simple Lambda_QCD*(p/q)*{{1,r2,phi}} land in [140,153] MeV.")
print("   BUT kappa is NOT chosen from that band -- it is item-130's (M_Delta-M_N)/2, fixed by")
print("   independent data. So 'decuplet step = kappa' is a 1-number POSTDICTION (0.2% on the")
print("   mean), not a free pick. That is materially stronger than the earlier 'free P_s' framing.")
print(f"   Still NOT derived: (a) the ~{sp[0]-sp[2]:.0f} MeV curvature (linear kappa misses it, ~0.4% residual);")
print(f"   (b) that P_s == kappa is an IDENTIFICATION (real QCD: constituent-Dm_s vs hyperfine differ),")
print(f"   plausible (both ~M_pi) but unproven; (c) kappa~M_pi itself is item-130/§9.9's claim,")
print(f"      and kappa=146.5 vs M_pi=139.6 is ~5% loose -- the decuplet wants kappa, not bare M_pi.")

H("VERDICT (revises su3_baryon_test probe 6(1))")
print("""The decuplet gap is NARROWER than the earlier consolidation claimed, but not
closed. Extending the ladder needs zero NEW constants IF one grants the identification
P_s = kappa (item-130's J^P tension): M = M_Delta + kappa*n_s then predicts the whole
decuplet to <=0.55% and Omega- to 0.06%, with the MEAN strange step matching kappa to
0.22%. kappa is fixed by (N,Delta) alone, so this is a genuine 1-number postdiction,
materially stronger than the earlier '2-point free P_s, Omega off 1%'. Equal-spacing
is structural here: on the decuplet I is slaved to n_s, so the Casimir collapses to
linear and one strange ladder reproduces Gell-Mann's rule.

WHAT IS NOT CLOSED (Tier: Proposition, not Locked):
 - The identification P_s = kappa is UNPROVEN. In standard QCD the J^P (N-Delta
   hyperfine) splitting and the constituent strange increment are DIFFERENT operators;
   their near-equality (146.5 vs 146.8) may be partly coincidental -- the [140,153] MeV
   band holds ~12 simple Lambda-ratios (S16.3), so a few-% match is not by itself strong.
 - The ~14 MeV curvature (steps 152.6 -> 148.8 -> 139.0) is unmodeled by linear kappa.
 - The scale question doesn't vanish, it RECURSES: 'why kappa ~ M_pi ~ Lambda_QCD/2?' is
   item-130 / §9.9 GMOR's open problem, not a new free strange mass. (kappa=146.5 vs bare
   M_pi=139.6 is itself 5% loose -- the decuplet prefers kappa, not the pion.)

NET: revise the su3 probe-6(1) line from 'P_s Target-B-blocked, free strange mass' to
'decuplet reachable parameter-free via P_s=kappa; the residual gap is justifying that
identification + the curvature, both Proposition-tier.'""")
