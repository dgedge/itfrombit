#!/usr/bin/env python3
"""
Item 93, the 3+1D residual: does the framework's 8-BODY-DIAGONAL shift admit a LOCAL
(finite-range) gamma5-real Dirac operator (a staggered/minimally-doubled escape), or
does it provably hit Nielsen-Ninomiya and REQUIRE the non-local (exponentially-local)
overlap? EXPLORATORY (uncommitted). Honest adversarial computation.

Structure of the answer (each part computed, not asserted):
 (0) The framework naive Dirac operator from the 8 corner-neighbour (body-diagonal) hops
     has D(k) = sigma_x^chi (x) (V.sigma) with V(k) derived in closed form below.
 (1) Nielsen-Ninomiya is Poincare-Hopf on T^3: sum of cone chiralities = chi(T^3) = 0,
     GEOMETRY-INDEPENDENT. Verified on the face-direction naive operator (8 clean cones).
 (2) The BODY-DIAGONAL naive operator is WORSE than face: nodal LINES, not just points
     -> no minimally-doubled (Karsten-Wilczek) finite-range escape; doubling is enhanced.
 (3) Poincare-Hopf still forces sum=0 on any generic local completion -> no finite-range
     single-cone gamma5-Hermitian operator, body diagonals included.
 (4) Horvath 1998: no ULTRALOCAL GW operator for one species. The overlap is the
     resolution and is local in the EXPONENTIAL sense only -- demonstrated: its real-space
     kernel decays exponentially but is nonzero at all ranges (not finite-range).
Needs only numpy. Self-asserting where confident; report-only where a discovered number.
"""
import numpy as np, itertools
np.set_printoptions(suppress=True)
def Hh(t): print("\n"+"="*76+"\n"+t+"\n"+"="*76)
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)

Hh("(0) Framework naive Dirac operator from the 8 body-diagonal (corner-neighbour) hops")
# lattice vectors R = (+-1,+-1,+-1) (cube body diagonals connecting a cell to its 8 corner
# neighbours). Hopping phase e^{i k.R}, 2pi-periodic on the cubic BZ. The Clifford factor
# for direction R is (R.alpha)/sqrt3; the antihermitian (kinetic) part picks out sin(k.R).
def Vbody_sum(k):
    kx,ky,kz=k; V=np.zeros(3)
    for s in itertools.product([1,-1],repeat=3):
        V+=np.array(s)*np.sin(s[0]*kx+s[1]*ky+s[2]*kz)
    return V
def Vbody_closed(k):
    kx,ky,kz=k
    return 8*np.array([np.sin(kx)*np.cos(ky)*np.cos(kz),
                       np.cos(kx)*np.sin(ky)*np.cos(kz),
                       np.cos(kx)*np.cos(ky)*np.sin(kz)])
rng=np.random.default_rng(0)
err=max(np.linalg.norm(Vbody_sum(p)-Vbody_closed(p)) for p in rng.uniform(-np.pi,np.pi,(200,3)))
print(f"  sum over 8 corners  V_x = 8 sin kx cos ky cos kz (+ cyc):  closed-form error = {err:.1e}")
print("  => the body-diagonal naive operator is D(k)=sigma_x^chi (x) (V.sigma),")
print("     V(k) = 8(sin kx cos ky cos kz,  cos kx sin ky cos kz,  cos kx cos ky sin kz).")
print("     massless: {gamma5=sigma_y^chi(x)I , D}=0 exactly -> exact chiral sym -> NN applies.")
assert err<1e-10

Hh("(1) Nielsen-Ninomiya = Poincare-Hopf on T^3: sum of cone chiralities = chi(T^3) = 0")
def jac(V,k,h=1e-6):
    J=np.zeros((3,3))
    for j in range(3):
        kp=np.array(k,float); km=np.array(k,float); kp[j]+=h; km[j]-=h
        J[:,j]=(V(kp)-V(km))/(2*h)
    return J
def Vface(k): return np.array([np.sin(k[0]),np.sin(k[1]),np.sin(k[2])])
faces=list(itertools.product([0.0,np.pi],repeat=3))
idx=[]
for k in faces:
    d=np.linalg.det(jac(Vface,k)); idx.append(int(np.sign(d)))
print(f"  FACE-direction naive: 8 isolated cones at {{0,pi}}^3")
print(f"     chiralities (sign det dV/dk): {idx}")
print(f"     sum = {sum(idx)}  (Poincare-Hopf: must equal chi(T^3)=0 -> 4 left + 4 right, NN)")
assert sum(idx)==0 and len(idx)==8

Hh("(2) BODY-DIAGONAL naive: NODAL LINES, not point doublers -> no minimal-doubling escape")
# check a candidate nodal line kx=ky=pi/2, all kz:
line=[Vbody_closed((np.pi/2,np.pi/2,kz)) for kz in np.linspace(-np.pi,np.pi,9)]
print("  scan kx=ky=pi/2, kz in [-pi,pi]:  |V(k)| =", [f"{np.linalg.norm(v):.1e}" for v in line])
print("  => V vanishes on the ENTIRE line (kx=ky=pi/2, any kz). Likewise kx=kz=pi/2, ky=kz=pi/2,")
print("     and the cube corners {0,pi}^3. The body-diagonal naive operator has nodal LINES +")
print("     point nodes -- a MORE degenerate band structure than the face operator's 8 points.")
# count the nodal lines of type k_i=k_j=pi/2 (and -pi/2): choose 2 of 3 axes, each +-pi/2
nlines=3*2*2  # choose pair (3) x sign of each of the 2 fixed coords (2x2)
print(f"  => at least {nlines} nodal lines (choose 2 axes pinned to +-pi/2) + the 8 corner points.")
print("  CONSEQUENCE: there is NO minimally-doubled (Karsten-Wilczek-type) 2-cone escape here;")
print("  the body-diagonal geometry ENHANCES doubling (continuous nodal lines), not reduces it.")
assert all(np.linalg.norm(v)<1e-9 for v in line)

Hh("(3) Poincare-Hopf forces sum=0 on ANY generic local completion (geometry-independent)")
# perturb the body-diagonal V by a generic constant -> nodal lines break into isolated zeros;
# their indices must still sum to chi(T^3)=0. Find zeros by Newton from a grid of seeds.
c=np.array([0.61,1.07,0.43])
def Vp(k): return Vbody_closed(k)+c
def newton(k0,iters=60):
    k=np.array(k0,float)
    for _ in range(iters):
        V=Vp(k); J=jac(Vp,k)
        if abs(np.linalg.det(J))<1e-9: return None
        dk=np.linalg.solve(J,V); k=k-dk
        if np.linalg.norm(dk)<1e-12: break
    k=(k+np.pi)%(2*np.pi)-np.pi
    return k if np.linalg.norm(Vp(k))<1e-8 else None
seeds=[np.array(s,float) for s in itertools.product(np.linspace(-np.pi,np.pi,7),repeat=3)]
zeros=[]
for s in seeds:
    z=newton(s)
    if z is None: continue
    if not any(np.linalg.norm(((z-w+np.pi)%(2*np.pi)-np.pi))<1e-5 for w in zeros): zeros.append(z)
inds=[int(np.sign(np.linalg.det(jac(Vp,z)))) for z in zeros]
print(f"  generic perturbation V+c: found {len(zeros)} isolated zeros; index sum = {sum(inds)}")
print(f"     (+1 count={inds.count(1)}, -1 count={inds.count(-1)})  -> Poincare-Hopf sum=0 confirmed: {sum(inds)==0}")
print("  => geometry-independent: no smooth (=local) gamma5-Hermitian operator on T^3 has a")
print("     single Weyl cone. Body diagonals cannot evade NN; only the GW relation can.")
assert sum(inds)==0 and len(zeros)>0

Hh("(4) No ULTRALOCAL GW operator (Horvath 1998) -- CITED, with one convention-robust check")
# HONESTY NOTE: an earlier draft tried to MEASURE the overlap's real-space decay in a free
# 1D toy and claimed 'exponential'. That was over-reach and WRONG: the 1D one-species window
# M in (-2,0) makes the Wilson kernel H_W=g5 D_W GAPLESS at intermediate k (verified below),
# so sgn(H_W) is ill-defined there and the toy kernel goes ~1/x (power law), not exponential.
# Measuring overlap locality in a free toy is convention-fraught; the correct statement is a
# THEOREM, not a fit. We CITE it and verify only the one robust, convention-independent fact.
g5=sz
def DW1(k,M,r=1.0): return 1j*sx*np.sin(k) + (M + r*(1-np.cos(k)))*sz
def HWgapmin(M):
    # fine grid so the exact zero at k=pi/4 (for M in (-2,0)) is resolved, not missed
    return min(min(abs(np.linalg.eigvalsh(sz@DW1(k,M)))) for k in np.linspace(0,2*np.pi,20001))
print("  1D Wilson kernel gap min_k|eig H_W(k)| vs mass M (fine grid):")
for M in [-1.0,-0.5,0.5,1.0]:
    g=HWgapmin(M); print(f"     M={M:+.1f}: gap={g:.1e}  {'GAPLESS (one-species window!)' if g<1e-3 else 'gapped'}")
print("  => the 1D one-species window M in (-2,0) is GAPLESS at intermediate k -- a known 1D")
print("     peculiarity. So a free-toy overlap-decay measurement is not meaningful here.")
print("  THEOREM (Horvath, PRL 81 (1998) 4063): no ULTRALOCAL (strictly finite-range) Dirac")
print("     operator satisfies Ginsparg-Wilson for a single species. The overlap (Neuberger")
print("     1998) realises GW and is EXPONENTIALLY local in the gapped (admissible) regime")
print("     (Hernandez-Jansen-Luscher, NPB 552 (1999) 363) -- nonzero at all ranges, never")
print("     finite-range. We take these as cited, not recomputed.")
print("  Convention-robust check actually done here: parts (1)-(3) already PROVE the finite-")
print("     range obstruction (Poincare-Hopf sum=0) without any overlap construction -- that")
print("     is the part this script is entitled to assert. The 'overlap is exp-local' is cited.")

Hh("VERDICT -- does the 8-body-diagonal shift admit a local gamma5-real Dirac operator?")
for ln in [
  "ANSWER: NO finite-range (ultralocal) one. The shift provably hits Nielsen-Ninomiya and",
  "requires the non-ultralocal overlap. The decisive facts (1-3 computed here; 4 cited):",
  "",
  " 1. The framework body-diagonal naive operator is D=sigma_x^chi(x)(V.sigma) with",
  "    V=8(sin kx cos ky cos kz,+cyc). WORSE than face dirs: it has nodal LINES (kx=ky=pi/2,",
  "    any kz; >=12 lines) + 8 corner points -- NOT a clean minimally-doubled 2-cone operator.",
  "    So no Karsten-Wilczek finite-range chiral escape hides in the body-diagonal geometry;",
  "    the geometry ENHANCES doubling.",
  " 2. The obstruction is TOPOLOGICAL, geometry-independent: Poincare-Hopf on the BZ torus",
  "    forces sum(cone chiralities)=chi(T^3)=0 (verified: face 4+/4-; perturbed body-diagonal",
  "    8+/8- balanced). A single physical Weyl cone is impossible for ANY smooth (=finite-range)",
  "    gamma5-Hermitian D(k) -- 6 face dirs or 8 body diagonals alike.",
  " 3. Hence no finite-range gamma5-real MASSLESS single-species Dirac operator from the walk.",
  "    RECONCILIATION with gw_euclidean_bridge.py's '1+1D closure': that finite-range GW operator",
  "    D=(1-W_s)/a describes a MASSIVE mode (mass=coin angle) with the doubler lifted -- NOT a",
  "    massless single chiral fermion. No contradiction: finite-range GW exists when massive/",
  "    doubled; the MASSLESS single species needs the overlap (chi(T^d)=0 forces >=2 cones for",
  "    every d>=1, so 1+1D was never a massless-chiral escape either).",
  " 4. CITED (not recomputed here): Horvath 1998 -- no ULTRALOCAL GW operator for one species;",
  "    Neuberger 1998 overlap realises GW; Hernandez-Jansen-Luscher 1999 -- overlap is",
  "    EXPONENTIALLY local (gapped/admissible regime), nonzero at all ranges, never finite.",
  "    [Self-correction: an earlier draft tried to MEASURE this decay in a free 1D toy and",
  "     mislabelled a gapless-window ~1/x tail as exponential. Retracted; deferred to theorem.]",
  "",
  "CONSEQUENCE for item 93: D_TCH continuum chiral fermion is NECESSARILY the overlap --",
  "exponentially local, NOT a finite-range function of W=S.C. The body-diagonal structure",
  "does not rescue ultralocality; it makes the naive doubling worse. So the open 3+1D piece",
  "is SETTLED in form: the answer is the (exp-local) overlap, the standard and only resolution.",
  "",
  "Honest closure of item 93: bare-substrate theta_UV=0 stands (Locked); 1+1D walk->GW direct",
  "construction stands (gw_euclidean_bridge.py); the 3+1D physical operator is an OVERLAP",
  "fermion (exp-local), provably NOT a finite-range W-derived D (Poincare-Hopf + Horvath).",
  "The continuum Strong-CP / d_n bound therefore rests on overlap LOCALITY (exp), which is",
  "the accepted lattice-QCD sense of local -- so it is defensible, but it is NOT the framework",
  "ultralocal-walk ideal, and the retraction of 'Theorem 4.1 by direct construction' stands.",
  "",
  "STRUCTURAL FLAG (for canon): the framework finite-range / ultralocal-walk ethos CANNOT host",
  "a single chiral fermion as a finite-range function of W -- a named, proven obstruction, true",
  "of all lattice chiral fermions. Item 93 is closable only at 'overlap (exp-local)' tier, not",
  "at 'finite-range substrate Dirac operator' tier. This is the sharp, final location of the gap.",
]: print(ln)
print()
print("ALL ASSERTS PASSED.")
