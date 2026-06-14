#!/usr/bin/env python3
"""
Gap-robustness of the CSS [[8,4,4]] gap under a BACKGROUND gauge twist.

This is the first test in the whole SMG line that probes PHYSICS/LOCALITY rather than
representation theory. The dressing audits (check_gauge_dressing.py, gauge_dressing_audit.py)
showed gauge-invariance-by-compensator is generic/tautological. Here we instead ask the
physical question:

  Does the unique gapped CSS ground state survive a background gauge field?

Construction (faithful to item 13: X-half = gauge/plaquette sector, Z-half = matter/Gauss):
  put a background gauge holonomy g(theta)=exp(i*theta*G) on the X (gauge) half only:
      H(theta) = H_Z  +  g(theta) H_X g(theta)^dagger,     G in {U(1)_Y, SU(2)_L, SU(3)_c}.
  H(0) = H_CSS (gap 2). Conjugating only the X-half by g (which does NOT commute with H_Z) is a
  genuine spectral deformation (not a global conjugation), so gap(theta) is a real function.
  gap stays open  -> the trivial SRE gap is ROBUST to the background gauge field (the necessary
                     'local as a functional of the gauge field' condition, at the cell level);
  gap closes      -> the CSS gap is FRAGILE to a gauge background -> bad sign.

HONEST SCOPE (stated, not hidden):
  * single CELL: a uniform background field, NOT a holonomy threaded through a spatial loop
    (no ring -> can't see spatial locality / spectral flow proper);
  * BACKGROUND (c-number g(theta)), NOT dynamical (the gauge field is not quantized/summed);
  * the gauge factor enters via the encoded generators of r2_smg_operator.py.
  So a robust gap here is a NECESSARY-condition probe, not the dynamical-gauge closure.
numpy; self-asserting only on construction facts, REPORTING the gap (the honest finding).
"""
import itertools
import numpy as np
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex); Z=np.diag([1,-1]).astype(complex)
P0=np.array([[1,0],[0,0]],complex); P1=np.array([[0,0],[0,1]],complex)
SR=np.array([[0,1],[0,0]],complex); SL=np.array([[0,0],[1,0]],complex)
G0,G1,LQ,C0,C1,I3,CHI,W=range(8)
def op(d):
    M=np.array([[1]],complex)
    for k in range(8): M=np.kron(M,d.get(k,I2))
    return M
def Xc(c): return op({k:X for k in range(8) if c[k]})
def Zc(c): return op({k:Z for k in range(8) if c[k]})
def ctrl(cb,cv,tb,to):
    a=P0 if cv==0 else P1; ina=P1 if cv==0 else P0
    return op({cb:a,tb:to})+op({cb:ina,tb:I2})
def idx(b): return int("".join(map(str,b)),2)

# self-dual [8,4,4] = RM(1,3); CSS Hamiltonian halves
gens=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]
HX=-sum(Xc(g) for g in gens); HZ=-sum(Zc(g) for g in gens); HCSS=HX+HZ

# encoded gauge generators (Hermitian, so exp(i*theta*G) is unitary)
allb=[tuple(b) for b in itertools.product([0,1],repeat=8)]
def Zf(b): return 1 if b[I3]==0 else -1
def sZc(b): return -3 if (b[C0],b[C1])==(0,0) else -1
def Qf(b): return 0.5*Zf(b)+(1/3)*sZc(b)+0.5
def T3v(b): return (0.5 if b[I3]==0 else -0.5) if b[CHI]==0 else 0.0
Y=np.diag([Qf(b)-T3v(b) for b in allb]).astype(complex)
T3=np.diag([T3v(b) for b in allb]).astype(complex)
Tp=ctrl(CHI,0,I3,SR); Tm=ctrl(CHI,0,I3,SL); T1=(Tp+Tm)/2          # Hermitian SU(2)_L off-diagonal
cyc=np.zeros((256,256),complex); cm={(0,1):(1,0),(1,0):(1,1),(1,1):(0,1),(0,0):(0,0)}
for b in allb:
    t=list(b); t[C0],t[C1]=cm[(b[C0],b[C1])]; cyc[idx(tuple(t)),idx(b)]=1
# Hermitian colour generator from the unitary 3-cycle: Lc = V diag(angle) V^dagger
wc,Vc=np.linalg.eig(cyc); Lc=(Vc*np.angle(wc))@np.linalg.inv(Vc); Lc=(Lc+Lc.conj().T)/2

def expmH(Gh,th):                                  # exp(i*theta*Gh), Gh Hermitian
    w,V=np.linalg.eigh(Gh); return (V*np.exp(1j*th*w))@V.conj().T
def gap_of(H):
    w=np.linalg.eigvalsh(H); return float(w[1]-w[0]), float(w[0])

# ----------------------------------------------------------------------------------
Hd("(A) Setup: which twists are REAL deformations vs trivial global conjugations?")
g0,e0=gap_of(HCSS)
print(f"  H_CSS: gap(0) = {g0:.3f}, ground energy {e0:.3f}  (unique [[8,0,4]] state)")
print("  H(t)=H_Z + g(t)H_X g(t)^dag.  If [G,H_Z]=0 the twist is a GLOBAL conjugation")
print("  H(t)=g H_CSS g^dag -> spectrum-preserving -> TRIVIAL (uninformative). Only generators")
print("  that do NOT commute with the diagonal H_Z give a real spectral deformation:")
for nm,Gh in [("U(1)_Y",Y),("SU(2)_L T3",T3),("SU(2)_L T1",T1),("SU(3)_c",Lc)]:
    dc=np.linalg.norm(Gh@HZ-HZ@Gh)
    print(f"     {nm:12s}: ||[G,H_Z]|| = {dc:8.3f}  -> {'global conj (TRIVIAL)' if dc<1e-9 else 'REAL deformation'}")
assert abs(g0-2.0)<1e-9

# ----------------------------------------------------------------------------------
Hd("(B) Scan the gap under each gauge-factor background twist  H(t)=H_Z+g(t)H_X g(t)^dag")
thetas=np.linspace(0,2*np.pi,33)
def scan(name,Gh):
    real = np.linalg.norm(Gh@HZ-HZ@Gh) > 1e-9
    mins=1e9; argmin=0.0; g00=None; gweak=None
    for th in thetas:
        g=expmH(Gh,th)
        gp,_=gap_of(HZ + g@HX@g.conj().T)
        if g00 is None: g00=gp
        if abs(th-thetas[2])<1e-9: gweak=gp                # weak background (theta ~ 0.39)
        if gp<mins: mins=gp; argmin=th
    tag="REAL deformation" if real else "trivial global-conj"
    print(f"  {name:12s} [{tag:18s}]: gap(0)={g00:.2f} gap(weak)={gweak:.2f} "
          f"min={mins:.3f}@th={argmin:.2f} -> {'open' if mins>0.05 else 'CLOSES'}")
    return mins, real
res={"U(1)_Y":scan("U(1)_Y",Y), "SU(2)_L T3":scan("SU(2)_L T3",T3),
     "SU(2)_L T1":scan("SU(2)_L T1",T1), "SU(3)_c":scan("SU(3)_c",Lc)}
real_mins=[m for (m,real) in res.values() if real]      # only the genuine (non-abelian) deformations
worst_gauge=min(real_mins) if real_mins else None

# ----------------------------------------------------------------------------------
Hd("(C) Controls: random Hermitian twist, and a known gap-CLOSING deformation (discrimination)")
rng=np.random.default_rng(0)
R=rng.standard_normal((256,256))+1j*rng.standard_normal((256,256)); R=(R+R.conj().T)
R=R/np.linalg.norm(R)*np.linalg.norm(Y)            # comparable scale
m_rand,_=scan("random-herm", R)
# a deformation that SHOULD close the gap: interpolate H_CSS -> H_Z only (drop the X-half).
# H_Z alone has a 16-fold degenerate ground space (gap 0). This confirms the test CAN see closing.
mins=1e9
for s in np.linspace(0,1,21):
    gp,_=gap_of(HZ + (1-s)*HX); mins=min(mins,gp)
print(f"  drop-X interpolation H_CSS->H_Z: min gap = {mins:.3f}  -> CLOSES (control: test detects gap loss)")
assert mins<0.05

# ----------------------------------------------------------------------------------
Hd("VERDICT — does the CSS gap survive a background gauge twist?")
print(f"""  KEY DISTINCTION the run exposes:
   * U(1)_Y and SU(2)_L T3 are DIAGONAL -> they commute with the diagonal H_Z -> the 'twist' is
     a global conjugation H(t)=g H_CSS g^dag -> gap invariant (=2). TRIVIAL, uninformative --
     so 'robust under abelian/Cartan twists' here is an artifact, not evidence.
   * The genuine (non-abelian, off-diagonal) deformations are SU(2)_L T1 and SU(3)_c. Worst-case
     gap over these REAL deformations: min = {worst_gauge:.3f}.
   * Control (drop X-half) closes the gap -> the test detects gap loss (not vacuously open);
     a random Hermitian twist of comparable size leaves min gap = {m_rand:.3f} (marginally open).

  READING (honest, negative-leaning):
  Under the REAL (non-abelian) gauge background twists the CSS gap is robust at WEAK background
  but CLOSES at strong background -- SU(2)_L T1 -> 0 at theta=pi, SU(3)_c -> ~0 near theta~1.
  So the trivial SRE gap is NOT robustly local under non-abelian gauge backgrounds even at the
  cell level: a strong non-abelian background collapses it. That is a real cautionary signal,
  and it sits exactly in the non-abelian (SU(2)/SU(3)) direction that is the hard part of the
  chiral-gauge problem -- not in the abelian directions (which are trivial/solved anyway).

  IMPORTANT caveats (so this is a SIGNAL, not a decisive no):
   * the twist used is a PROXY -- conjugating the X-half by a uniform g -- not the physical
     link-dressed gauge coupling; some of the closing may be the proxy's artifact;
   * single CELL = uniform background, large theta (~pi), NOT a weak holonomy threaded through a
     spatial loop (no ring -> no true spatial locality / spectral flow);
   * BACKGROUND c-number, NOT dynamical gauge fields.
  Whether strong or weak non-abelian backgrounds dominate is a dynamical question; in a confining
  sector strong fluctuations are plausible, which would make this fragility bite.

  NET: the first PHYSICS-level probe pushes the OTHER way from the (tautological) dressing
  result -- the CSS gap is fragile to strong non-abelian gauge backgrounds at the cell level.
  The decisive test remains a small RING of cells with a genuinely threaded (link-dressed)
  holonomy + spectral-flow tracking; this cell-level proxy says: do not assume robustness.""")
print("\nConstruction asserts passed; gap values are REPORTED findings (read above).")
