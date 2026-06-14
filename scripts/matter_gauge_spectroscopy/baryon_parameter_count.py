#!/usr/bin/env python3
"""
How many genuinely-free scales does the baryon octet+decuplet have, and which does
the framework supply? EXPLORATORY (uncommitted). Follows the DGG seam from
baryon_step_operators.py.

Thesis under test: the framework presents THREE baryon mechanisms --
  item 129  (C_8 Nyquist fractions {4,6,9}/8 * Lambda)  -> octet strangeness
  item 130  (O_h irrep dim, kappa)                      -> N-Delta J^P split
  §9.10     (w = alpha*Lambda)                          -> EM/isospin
-- but in QCD the whole spin-averaged spectrum collapses to just THREE physical
scales {m_u, m_s, A_hyperfine} (De Rujula-Georgi-Glashow). We fit those 3 to all 8
isospin-averaged masses, map them onto framework constants, and show items 129 & 130
SHARE one operator (A) -- the genuine redundancy.

DGG: M = sum_i m_i + A * sum_{i<j} <S_i.S_j>/(m_i m_j).
Spin sums (ground-state, symmetric space => spin-flavour symmetric):
  octet J=1/2: N,Lam <SS>net via diquark; decuplet J=3/2: all pairs +1/4.
Nothing fitted beyond the 3 DGG scales; every framework mapping is a CHECK, not a fit.
"""
import numpy as np, mpmath as mp
mp.mp.dps=20
def H(s): print("\n"+"="*74+"\n"+s+"\n"+"="*74)

# isospin-averaged PDG masses (MeV), EM removed by averaging
M = {'N':938.92,'Lam':1115.68,'Sig':1193.15,'Xi':1318.29,
     'Del':1232.0,'Sig*':1384.6,'Xi*':1533.4,'Om':1672.45}
NL = {'N':3,'Lam':2,'Sig':2,'Xi':1,'Del':3,'Sig*':2,'Xi*':1,'Om':0}   # light quarks
NS = {'N':0,'Lam':1,'Sig':1,'Xi':2,'Del':0,'Sig*':1,'Xi*':2,'Om':3}   # strange quarks
ORD=['N','Lam','Sig','Xi','Del','Sig*','Xi*','Om']
L=332.0; MPI=139.57; KAP=146.5

def hyper(b,u,s):
    """<sum S_i.S_j/(m_i m_j)> for baryon b given constituent masses u,s."""
    uu,us,ss=1/u**2,1/(u*s),1/s**2
    return {  # coefficients * (1/4 base unit folded in)
      'N':(-0.75)*uu, 'Lam':(-0.75)*uu,
      'Sig':(0.25)*uu-1.0*us, 'Xi':(0.25)*ss-1.0*us,
      'Del':(0.75)*uu, 'Sig*':(0.25)*uu+0.5*us,
      'Xi*':0.5*us+0.25*ss, 'Om':(0.75)*ss }[b]

H("1. Global 3-parameter DGG fit (m_u, m_s, A) to all 8 isospin-averaged masses")
best=None
for u in np.arange(340.,390.,0.5):
    for s in np.arange(490.,600.,0.5):
        G=np.array([hyper(b,u,s) for b in ORD])
        base=np.array([NL[b]*u+NS[b]*s for b in ORD])
        y=np.array([M[b] for b in ORD])
        A=float(((y-base)@G)/(G@G))                 # linear LSQ for A
        pred=base+A*G; rms=float(np.sqrt(np.mean((pred-y)**2)))
        if best is None or rms<best[0]: best=(rms,u,s,A,pred)
rms,u,s,A,pred=best
print(f"   best fit: m_u={u:.1f}  m_s={s:.1f}  A={A:.3e}  A/m_u^2={A/u**2:.1f}  RMS={rms:.1f} MeV over 8 baryons")
print(f"   {'baryon':7s}{'pred':>9s}{'PDG':>9s}{'dev':>7s}{'%':>7s}")
for b,p in zip(ORD,pred):
    print(f"   {b:7s}{p:9.1f}{M[b]:9.1f}{p-M[b]:+7.1f}{(p-M[b])/M[b]*100:+7.2f}")
print(f"   => 3 scales reproduce 8 masses; residuals concentrate in the OCTET (the spin-")
print(f"      dependent sector). 8 data - 3 params = 5 genuine predictions.")

H("2. UNIFICATION: are N-Delta and Sigma-Lambda the SAME hyperfine A? (items 130 vs 129)")
Au=A/u**2; r=u/s
A_ND=(M['Del']-M['N'])/1.5                      # Delta-N = A/u^2 * 3/2
A_SL=(M['Sig']-M['Lam'])/(1-r)                  # Sigma-Lam = A/u^2 *(1-r)
print(f"   N-Delta   gives A/u^2 = (M_Del-M_N)/1.5      = {A_ND:.1f}")
print(f"   Sigma-Lam gives A/u^2 = (M_Sig-M_Lam)/(1-r)  = {A_SL:.1f}   (r=m_u/m_s={r:.3f})")
print(f"   ratio {A_SL/A_ND:.2f}: the SAME operator, ~{abs(A_SL/A_ND-1)*100:.0f}% sloppy across octet/decuplet")
print(f"   (the known DGG m_s tension: octet wants m_s~{u/(1-(M['Sig']-M['Lam'])/A_ND):.0f}, decuplet/Omega wants ~{s:.0f}).")
print(f"   KEY: item-130 kappa = (3/4)(A/u^2) = {0.75*A_ND:.1f}  ->  kappa IS the hyperfine A, up to 3/4.")
print(f"        item-129's 6/8-vs-4/8 isospin split (Sigma-Lambda) is the SAME A acting in the octet.")
print(f"   => items 129 and 130 are NOT two mechanisms; they are ONE chromomagnetic operator A,")
print(f"      fragmented across two anchors. Merging them removes a redundancy (the GMO move).")

H("3. Map the 3 DGG scales onto framework constants -- which are native, which free?")
print(f"   (a) m_u = {u:.1f}  vs  Lambda_QCD = {L:.0f}   -> ratio {u/L:.3f} ({(u/L-1)*100:+.0f}%).")
print(f"       Constituent m_u ~ chiral-breaking scale ~ Lambda_QCD is QCD folklore; framework-")
print(f"       native at the ~9%% level (the framework uses M_N=2sqrt2*Lambda directly, bypassing m_u).")
print(f"   (b) A (hyperfine) <-> kappa: kappa=(3/4)(A/u^2). EXACT when A/u^2 from (M_Del-M_N):")
print(f"       (3/4)*195.4 = {0.75*195.4:.1f} = kappa={KAP}. The GLOBAL-fit A/u^2={Au:.1f} gives")
print(f"       (3/4)*{Au:.1f}={0.75*Au:.1f} -- the 1%% gap IS the octet/decuplet tension from sec.2.")
print(f"       And kappa is empirical (item-130 fits N,Delta), 'kappa~M_pi' Proposition")
print(f"       ({(KAP/MPI-1)*100:+.0f}%% loose). So A is native ONLY as far as kappa is.")
print(f"   (c) m_s = {s:.1f}: the genuinely NEW strange scale. m_s-m_u={s-u:.0f}. NOT fixed by Lambda")
print(f"       or kappa. This is the one Target-B continuous parameter of the baryon sector.")

H("4. Does item-129 'derive' m_s, or just discretize it? (the {4,6,9}/8 claim)")
print(f"   item-129 adds integer*Lambda/8 to M_N=939; unit Lambda/8 = {L/8:.1f} MeV.")
print(f"   What integers does the DATA demand? round((M_b - M_N)/(Lambda/8)):")
for b in ['Lam','Sig','Xi']:
    real=(M[b]-M['N']); k=real/(L/8)
    print(f"     {b:4s}: (M-{M['N']:.0f})={real:6.1f}  /41.5 = {k:5.2f}  -> nearest int {round(k)}  (item-129 uses {{4,6,9}}[{['Lam','Sig','Xi'].index(b)}])")
print(f"   => {{4,6,9}} are round(data/(Lambda/8)) to ~3%%. item-129 is a DISCRETIZATION of the")
print(f"      measured splittings in units of Lambda/8, NOT a 0-parameter derivation: the unit")
print(f"      is native (Lambda/8) but the integers are READ OFF the data. Info ~ 1 strange scale.")

H("5. Numerology kill: is m_s = Lambda*phi (or any 'nice' form) earned?")
phi=(1+mp.sqrt(5))/2
cands=[('Lambda*phi',float(L*phi)),('Lambda*sqrt(5/2)',float(L*mp.sqrt(mp.mpf(5)/2))),
       ('Lambda*8/5',L*8/5),('3*Lambda/2-... n/a',0)]
ms_oct = u/(1-(M['Sig']-M['Lam'])/A_ND)
print(f"   m_s is only defined to ~{abs(ms_oct-s)/s*100:.0f}%%: decuplet(Omega)->{s:.0f}, octet(Sig-Lam)->{ms_oct:.0f}.")
for nm,v in cands[:3]:
    print(f"     {nm:18s}= {v:6.1f}  (|dev| from decuplet m_s = {abs(v-s):.0f} = {abs(v-s)/s*100:.1f}%%)")
print(f"   Lambda*phi=537 sits 0.4%% from the decuplet m_s -- BUT m_s itself is only ~12%% defined,")
print(f"   so 0.4%% is spurious precision. With phi in the alphabet (S16.3), this is the SAME")
print(f"   numerology trap flagged for P_s. VERDICT: m_s=Lambda*phi NOT earned. m_s stays free.")

H("6. The classic 3-quark win that the framework CANNOT yet take: mu_p/mu_n = -3/2")
print(f"   SU(6) spin-flavour wavefunctions give mu_p/mu_n = -3/2 parameter-free; exp = {2.793/-1.913:.3f}")
print(f"   (2.7%% from -1.5). This needs the baryon SPIN-FLAVOUR wavefunction -- exactly the §15")
print(f"   item-138 spin gap. So the single cleanest constituent-model prediction is GATED behind")
print(f"   the same missing piece that blocks the within-multiplet hierarchy. One root cause.")

H("NET -- the seam, quantified")
print(f"""The baryon octet+decuplet is a THREE-scale system {{m_u={u:.0f}, m_s={s:.0f}, A/u^2={Au:.0f}}};
those 3 reproduce all 8 spin-averaged masses to RMS {rms:.0f} MeV (octet-dominated).
Mapping onto the framework:
  * m_u  <-> Lambda_QCD     : native at ~9% (constituent mass ~ chiral scale).
  * A    <-> kappa (item130): A/u^2 = (4/3)kappa EXACTLY; native only as far as kappa
             is (empirical, 'kappa~M_pi' Proposition, 5% loose).
  * m_s                     : the ONE genuinely free strange scale. Target-B. item-129
             does not derive it -- it DISCRETIZES it as round(split/(Lambda/8))={{4,6,9}}.
Two structural results:
  (i)  items 129 and 130 SHARE the hyperfine operator A -- they are one mechanism
       fragmented, not two. (Sigma-Lambda and N-Delta are both A.) Merge -> 1 GMO formula.
  (ii) the baryon sector's irreducible free content is ~ONE strange scale; everything
       else is Lambda + kappa. The framework re-expresses QCD's 3-scale economy in
       substrate language -- neither beating nor beaten by DGG -- with one redundancy to fix.
The cleanest forward win (mu_p/mu_n=-3/2) and the curvature and the within-multiplet
hierarchy are ALL gated by the same missing spin-flavour wavefunction (item 138).""")
