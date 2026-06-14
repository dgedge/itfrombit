#!/usr/bin/env python3
"""
Specialness audit for the proposed first-generation quark MASS suppression by a
colour-trace factor 1/3 (tch_yukawa_hierarchy_probe.py + tch_quark_node_trace_proof.py).

Distinct from ANCHOR item 96, which is the quark Koide PHASES (delta_u=2/27, delta_d=1/9).
This is the proposed mass-MAGNITUDE m_light -> m_light/3. Question: is 1/3 derived, or
chosen to fit? Asserts on every quoted number (evidence-not-assertion).

Findings (all asserted below):
  (1) the colour-resolved operator gives m_node/3 for WHATEVER generation is tagged rank-one;
      light is picked because it fits (up 7.4% / down 9.8%) -- middle/heavy are far worse.
      => the 1/3 is an input, not a derivation (same lesson as the gauge-dressing tautology).
  (2) the down sector's BEST operator description is R=sqrt(12/5), NO trace, 1.2%; the trace
      HARMS it (->17.4%). down needs no colour trace.
  (3) continuous sweep: up prefers ~1/7 (2.3%), not 1/3 (7.4%); down is s-floored at 9.8%
      independent of the d-factor. The 10% log-ratio bar passes a ~x4-9 wide band.
  (4) any rank-one colour projector has Tr/3 = 1/3 (the Boltzmann weights are irrelevant).
"""
import math, random
ORD={"up":["u","c","t"],"down":["d","s","b"]}
M={"up":{"u":2.16e-3,"c":1.27,"t":172.76},"down":{"d":4.67e-3,"s":93.4e-3,"b":4.18}}
LIGHT,MIDDLE,HEAVY=1,2,0
def lam(R,d,n): return 1+R*math.cos(d+2*math.pi*n/3)
def masses(R,d): return [lam(R,d,n)**2 for n in range(3)]
def obs(sec):
    h=ORD[sec][-1]; return {n:math.log10(M[sec][n]/M[sec][h]) for n in ORD[sec]}
def predict(R,d,sec,factor=1.0,which=LIGHT):
    mv=masses(R,d); em=[mv[n]*(factor if n==which else 1.0) for n in range(3)]
    nm=ORD[sec]; a={nm[0]:em[LIGHT],nm[1]:em[MIDDLE],nm[2]:em[HEAVY]}; h=nm[-1]
    return {n:math.log10(a[n]/a[h]) for n in nm}
def maxrel(sec,p):
    o=obs(sec); return max(abs(p[n]-o[n])/abs(o[n]) for n in o if abs(o[n])>1e-12)
def hd(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)

UP=(math.sqrt(3),2/27); DN=(math.sqrt(2),1/9)

hd("(1) which generation is rank-one is FREE (the 1/3 is an input, not derived)")
res={}
for sec,(R,d) in [("up",UP),("down",DN)]:
    row={}
    for which,lab in [(LIGHT,"light"),(MIDDLE,"middle"),(HEAVY,"heavy")]:
        e=maxrel(sec,predict(R,d,sec,1/3,which)); row[lab]=e
        print(f"  {sec:4s} rank-one on {lab:6s}: max rel = {e:5.1%}")
    res[sec]=row
assert abs(res["up"]["light"]-0.074)<0.003 and res["up"]["middle"]>0.15 and res["up"]["heavy"]>0.20
assert abs(res["down"]["light"]-0.098)<0.003 and res["down"]["middle"]>0.15 and res["down"]["heavy"]>0.30
print("  -> light gives 7.4%/9.8%; middle&heavy far worse -> light chosen to fit, not derived.")

hd("(2) branch audit: the down sector's best description uses NO trace")
tab=[("up   R=sqrt3     ","up",UP),("down R=sqrt2     ","down",DN),
     ("down R=sqrt(12/5)","down",(math.sqrt(12/5),1/9))]
vals={}
for lab,sec,(R,d) in tab:
    nt=maxrel(sec,predict(R,d,sec,1.0)); tr=maxrel(sec,predict(R,d,sec,1/3))
    vals[lab.strip()]=(nt,tr); print(f"  {lab}:  no-trace {nt:5.1%} | with d/3 {tr:5.1%}")
nt_b,tr_b=vals["down R=sqrt(12/5)"]
assert nt_b<0.02 and tr_b>0.15   # 1.2% no-trace, 17.4% with-trace
print("  -> down R=sqrt(12/5) NO-trace 1.2% beats down R=sqrt2 + d/3 9.8%; trace HARMS it.")

hd("(3) continuous sweep: up prefers ~1/7; down is s-floored, indifferent to the d-factor")
for sec,(R,d),first in [("up",UP,"u"),("down",DN,"d")]:
    base=predict(R,d,sec,1.0); o=obs(sec)
    fstar=10**(o[first]-base[first])   # factor that exactly nails the first-gen ratio
    e_star=maxrel(sec,predict(R,d,sec,fstar))
    e_third=maxrel(sec,predict(R,d,sec,1/3))
    print(f"  {sec:4s}: f*(exact {first}) = {fstar:.3f} (~1/{1/fstar:.1f}), sector max {e_star:.1%}"
          f"  | at 1/3: {e_third:.1%}")
    # passing band on the 10% bar
    lo=hi=None; f=0.02
    while f<=1.5:
        if maxrel(sec,predict(R,d,sec,f))<=0.10:
            lo=lo or f; hi=f
        f+=0.002
    print(f"        10%-bar passing band [{lo:.3f},{hi:.3f}] (width x{hi/lo:.1f})")
# up minimizer ~1/7, NOT 1/3
up_fstar=10**(obs("up")["u"]-predict(*UP,"up",1.0)["u"])
assert 0.12<up_fstar<0.17           # ~0.145 ~ 1/7, not 1/3
# down s-floored: every f in [1/3..1/7] gives the same 9.8% max
floor=[round(maxrel("down",predict(*DN,"down",1/k)),3) for k in (3,4,5,6,7)]
assert floor==[0.098]*5
print(f"  down max rel for d-factor 1/3..1/7 = {floor} -> s-floored, d-factor irrelevant.")

hd("(4) Tr(rank-one)/3 = 1/3 is rank-one triviality (Boltzmann weights irrelevant)")
random.seed(0)
for k in range(4):
    v=[complex(random.gauss(0,1),random.gauss(0,1)) for _ in range(3)]
    nrm=math.sqrt(sum(abs(x)**2 for x in v)); v=[x/nrm for x in v]
    t=sum(abs(v[i])**2 for i in range(3))/3
    assert abs(t-1/3)<1e-12
print("  4 random rank-one colour projectors: Tr/3 = 1/3 exactly -> no colour content.")

print("\nALL ASSERTS PASSED. The 1/3 mass suppression is a target-to-derive, not a result:")
print("  it is a fit-chosen input; down needs no trace; up prefers 1/7.")
