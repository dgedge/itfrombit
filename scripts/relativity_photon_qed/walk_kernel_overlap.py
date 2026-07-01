#!/usr/bin/env python3
"""
Item 93/143 — the FINAL step: feed the GENUINE TCH walk kernel (item 97 / W=S.C, the
8-body-diagonal shift) through the domain-wall->overlap reduction of
bulk_domainwall_overlap.py, and count physical species. Does the walk lift to ONE
chiral fermion (continuum-Locked) or not? EXPLORATORY, honest, self-asserting (numpy).

Walk kernel (reproduced from gw_nogo_bodydiagonal.py, ANCHOR Sec 3.1/3.5):
  kinetic  D_kin(k) = i * sum_i V_i(k) alpha_i,  alpha_i = sigma_x^chi (x) sigma_i
  V(k) = 8(sin kx cos ky cos kz, cos kx sin ky cos kz, cos kx cos ky sin kz)
  body-diagonal Wilson scalar (the WALK-SOURCED doubler-lifter, finite-range):
  W_body(k) = sum_{8 corners}(1 - cos(k.R)) = 8(1 - cos kx cos ky cos kz)
  D_W^walk(k) = D_kin + (M0 + r W_body) I    ;   H_W = gamma5 D_W^walk ;  D_ov = 1 + g5 sgn(H_W)

Honesty: this asks whether the walk's OWN body-diagonal Wilson term lifts all doublers to
one species. gw_nogo found the body-diagonal NAIVE operator has nodal LINES ('enhanced
doubling'); here we test the Wilson-completed operator at the overlap level. The number of
surviving species is the verdict, computed not assumed.
"""
import numpy as np, itertools
def kron(a,b): return np.kron(a,b)
I2=np.eye(2,dtype=complex); I4=np.eye(4,dtype=complex)
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
def anti(A,B): return A@B+B@A
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
a1=kron(sx,sx); a2=kron(sx,sy); a3=kron(sx,sz); alpha=[a1,a2,a3]
g5=kron(sy,I2)

def Vbody(k):
    kx,ky,kz=k
    return 8*np.array([np.sin(kx)*np.cos(ky)*np.cos(kz),
                       np.cos(kx)*np.sin(ky)*np.cos(kz),
                       np.cos(kx)*np.cos(ky)*np.sin(kz)])
def Wbody(k):
    kx,ky,kz=k
    return 8*(1-np.cos(kx)*np.cos(ky)*np.cos(kz))
def Wface(k):                       # EXTERNAL (NOT walk-sourced) face-direction Wilson term
    return sum(1-np.cos(ki) for ki in k)

def DW_walk(k,M0=-1.0,r=1.0,wilson="body"):
    V=Vbody(k)
    kin=1j*sum(V[i]*alpha[i] for i in range(3))
    W=Wbody(k) if wilson=="body" else Wface(k)
    return kin + (M0 + r*W)*I4
def HW(k,**kw): return g5@DW_walk(k,**kw)
def sgnH(M):
    w,Vv=np.linalg.eigh(M); return (Vv*np.sign(w))@Vv.conj().T
def Dov(k,**kw): return I4 + g5@sgnH(HW(k,**kw))

Hd("(A) walk kernel is genuine: D_kin = i sum V_i alpha_i (body-diagonal), H_W Hermitian")
rng=np.random.default_rng(0)
pts=[rng.uniform(-np.pi,np.pi,3) for _ in range(300)]
herm=max(np.linalg.norm(HW(p)-HW(p).conj().T) for p in pts)
chiral=max(np.linalg.norm(anti(g5,1j*sum(Vbody(p)[i]*alpha[i] for i in range(3)))) for p in pts)
print(f"  H_W = g5 D_W^walk Hermitian: residual {herm:.1e}")
print(f"  naive kinetic is exactly chiral {{g5,D_kin}}=0: residual {chiral:.1e}  (massless walk Dirac)")
assert herm<1e-10 and chiral<1e-10

Hd("(B) overlap from the walk kernel is still GW + g5-Hermitian (reduction machinery is fine)")
gw =max(np.linalg.norm(anti(g5,Dov(p))-Dov(p)@g5@Dov(p)) for p in pts)
gh =max(np.linalg.norm(g5@Dov(p)@g5-Dov(p).conj().T) for p in pts)
print(f"  max ||{{g5,D_ov}}-D_ov g5 D_ov|| : {gw:.1e}   (GW, exact)")
print(f"  max ||g5 D_ov g5 - D_ov^dag||    : {gh:.1e}   (g5-Hermitian, exact)")
assert gw<1e-9 and gh<1e-9

Hd("(C) SPECIES COUNT — does the WALK-SOURCED (body-diagonal) Wilson term give ONE species?")
corners=list(itertools.product([0.0,np.pi],repeat=3))
def species_corners(wilson,M0=-1.0):
    n=0; lst=[]
    for p in corners:
        if np.min(np.abs(np.linalg.eigvals(Dov(np.array(p),M0=M0,wilson=wilson))))<1e-8:
            n+=1; lst.append(tuple(int(round(x/np.pi)) for x in p))
    return n,lst
print("  body-diagonal Wilson scalar at the 8 corners (n = #pi-components):")
for p in corners:
    n=sum(1 for x in p if abs(x-np.pi)<1e-9)
    print(f"     k/pi={tuple(int(round(x/np.pi)) for x in p)}  n={n}  W_body={Wbody(np.array(p)):4.0f}  "
          f"{'(NOT lifted: W_body=0)' if abs(Wbody(np.array(p)))<1e-9 else '(lifted)'}")
nb,lb=species_corners("body")
print(f"\n  WALK-SOURCED body-diagonal Wilson, M0=-1: massless corners = {nb}  at {lb}")
# also scan a body-diagonal nodal line kx=ky=pi/2 to be sure it is lifted
linemin=min(np.min(np.abs(np.linalg.eigvals(Dov(np.array([np.pi/2,np.pi/2,kz]),wilson="body")))) for kz in np.linspace(-np.pi,np.pi,41))
print(f"  along the naive nodal line kx=ky=pi/2: min |D_ov eig| = {linemin:.2f}  (lifted, gapped)")
nf,lf=species_corners("face")
print(f"\n  EXTERNAL face-direction Wilson (NOT walk-sourced), M0=-1: massless corners = {nf}  at {lf}")
print("\n  -> the body-diagonal nodal LINES are lifted by W_body, BUT the even corners")
print("     (n even: k=0 and the three (pi,pi,0)-type) have W_body=0 and survive.")
print(f"  VERDICT: walk-sourced kernel -> {nb} species ; a (non-walk) face Wilson term -> {nf} species.")
# PIN the honest-negative result so exit-0 certifies it (assert-on-every-quoted-number):
assert nb==4 and nf==1 and linemin>1.0   # walk's own kernel keeps 4 species (NOT 1); only the
#   non-walk face-Wilson term gives 1; the naive nodal line IS lifted (gapped) -> the closing
#   fails specifically at the 4 even corners, not for lack of nodal-line lifting.

Hd("(D) the reduction (domain-wall) machinery still converges on the walk kernel")
def eps_Ls(M,Ls):
    w,Vv=np.linalg.eigh(M); T=(1-w)/(1+w); e=(1-T**Ls)/(1+T**Ls); return (Vv*e)@Vv.conj().T
def Dov_Ls(p,Ls,**kw): return I4+g5@eps_Ls(HW(p,**kw),Ls)
tp=[rng.uniform(-np.pi,np.pi,3) for _ in range(40)]
print(f"  {'L_s':>4} | max ||D_ov^Ls - D_ov||  (walk kernel, body Wilson)")
prev=None
for Ls in [8,16,24,32]:
    e=max(np.linalg.norm(Dov_Ls(p,Ls)-Dov(p)) for p in tp)
    print(f"  {Ls:>4} | {e:.3e}" + (f"   ({prev/e:.1f}x)" if prev else ""))
    prev=e
print("  -> the bulk-domain-wall reduction converges regardless of WHICH kernel: the machinery")
print("     is kernel-agnostic. What it converges TO is the multi-species walk overlap.")

Hd("VERDICT — does the walk kernel move item 93/143 to continuum-Locked?")
locked = (nb==1)
lines_no = [
 "NO. The genuine TCH walk kernel does NOT reduce to a single chiral fermion.",
 "",
 f" * The body-diagonal Wilson term W_body=8(1-cos kx cos ky cos kz) IS walk-sourced and",
 f"   finite-range, and it DOES lift the naive nodal LINES (min|D_ov| on kx=ky=pi/2 = {linemin:.2f}).",
 f" * BUT it vanishes on the EVEN corners (n=0,2): k=0 plus (pi,pi,0),(pi,0,pi),(0,pi,pi).",
 f"   So the walk-sourced overlap keeps {nb} degenerate species, not 1 -- the 'enhanced",
 "   doubling' of the body-diagonal geometry (gw_nogo) made concrete at the overlap level.",
 f" * A face-direction Wilson term gives the single species ({nf}), but that term is NOT",
 "   sourced by the body-diagonal shift S -- it is the EXTERNAL kernel gw_nogo named.",
 "",
 "So: the domain-wall->overlap REDUCTION is kernel-agnostic and closed (bulk_domainwall_",
 "overlap.py); but the WALK's own kernel + its own Wilson term yields a 4-species overlap.",
 "Selecting one species still requires an external (non-walk) doubler-lifting term.",
 "",
 "CONSEQUENCE for item 93/143: NOT continuum-Locked. The residual is now pinned EXACTLY:",
 "a single-species, walk-SOURCED doubler-lifting term. gw_nogo proved no finite-range one",
 "exists (Poincare-Hopf); this script shows the natural body-diagonal Wilson term gives 4",
 "species, not 1. Continuum closure stays at 'overlap (exp-local), external-Wilson' tier.",
 "The bare-substrate theta_UV=0 (item 93) and arg det M=0 (item 143) are untouched and stand;",
 "only the continuum doubler-lifting remains external. HONEST: the final step does not close.",
]
lines_yes = ["UNEXPECTED: the walk kernel gives a single species -- investigate (would be a real advance)."]
for ln in (lines_yes if locked else lines_no):
    print(ln)
print("\nALL ASSERTS PASSED (computed claims verified; verdict is the honest species count).")
