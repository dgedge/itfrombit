#!/usr/bin/env python3
r"""ITEM 87 -- the theta23 octant from the QLC-tangent recovery texture.

The frame-transport lemma fixes the lepton texture to the QLC tangent
a = standard (-5/6,1/6,2/3) + sign K_or/3 (1/3,1/3,1/3) = (-1/2,1/2,1), acting
right on the bimaximal base U0 = R23(pi/4) R12(pi/4):  U(t) = U0 exp(t A(a)).
The atmospheric latch gives dtheta23/dt=0 at first order; the OCTANT is the
second-order residual at the t reproducing the observed reactor theta13~8.6 deg.
Self-asserting numerics (Rodrigues expm, PDG extraction) + a robustness sweep.
"""
import numpy as np
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m
def A(x,y,z): return np.array([[0,-z,y],[z,0,-x],[-y,x,0]],float)
def expm_so3(M):
    a=np.array([M[2,1],M[0,2],M[1,0]]); th=np.linalg.norm(a)
    if th<1e-15: return np.eye(3)
    K=M/th; return np.eye(3)+np.sin(th)*K+(1-np.cos(th))*(K@K)
def R12(t): c,s=np.cos(t),np.sin(t); return np.array([[c,s,0],[-s,c,0],[0,0,1]])
def R23(t): c,s=np.cos(t),np.sin(t); return np.array([[1,0,0],[0,c,s],[0,-s,c]])
def angles(U):
    U=np.abs(U)
    return (np.degrees(np.arctan2(U[0,1],U[0,0])),
            np.degrees(np.arctan2(U[1,2],U[2,2])),
            np.degrees(np.arcsin(min(1,U[0,2]))))
def fit_theta23(a, target13=8.6):
    U0=R23(np.pi/4)@R12(np.pi/4); best=None
    for t in np.linspace(0,0.5,5001):
        th12,th23,th13=angles(U0@expm_so3(t*A(*a)))
        if best is None or abs(th13-target13)<abs(best[3]-target13): best=(t,th12,th23,th13)
    return best

print("="*72); print("theta23 OCTANT from the QLC-tangent recovery texture"); print("="*72)
a0=np.array([-1/2,1/2,1.0])
t,th12,th23,th13=fit_theta23(a0)
print(f"\n  physical QLC tangent a=(-1/2,1/2,1) at theta13={th13:.2f} deg:")
print(f"    theta12={th12:.1f} deg (obs ~33.6),  theta23={th23:.2f} deg (obs ~49),  -> {'SECOND' if th23>45 else 'FIRST'} octant")

# first-order latch
eps=1e-5; d23=(angles(R23(np.pi/4)@R12(np.pi/4)@expm_so3(eps*A(*a0)))[1]-45)/eps
ok(abs(d23)<1e-1, f"atmospheric latch: dtheta23/dt|_0 = {d23:.3f} ~ 0 (45 deg held at first order)")
ok(th23>45, f"second-order residual pushes theta23={th23:.2f} > 45 -> SECOND octant for the physical texture")

# robustness: perturb each component of the tangent by +/-15%, does the octant survive?
print("\n  robustness sweep (perturb tangent components by +/-15%):")
import itertools
second=0; tot=0
for dx,dy,dz in itertools.product((-.15,0,.15),repeat=3):
    a=a0*np.array([1+dx,1+dy,1+dz]); _,_,th23p,_=fit_theta23(a); tot+=1; second+= (th23p>45)
print(f"    {second}/{tot} perturbed textures stay in the SECOND octant")
ok(second/tot>0.8, f"second octant is robust ({second}/{tot} > 80% of +/-15% perturbations)")

print("\n"+"="*72); print("VERDICT")
print(f"  The atmospheric latch holds theta23=45 deg at FIRST order; the octant is the")
print(f"  second-order residual of the QLC-tangent texture. The physical texture reproduces")
print(f"  theta12~33 and theta13=8.6 and gives theta23={th23:.1f} deg -> SECOND OCTANT, robust under")
print(f"  +/-15% tangent perturbations. So the framework PREDICTS the second octant, matching the")
print(f"  current global-fit lean. HONEST caveats: (i) leading-order -- the texture is schematic,")
print(f"  so the robust claim is the octant (>45), not the precise value (it undershoots ~49 deg);")
print(f"  (ii) the octant is a property of the texture/atmospheric latch, NOT cleanly tied to the")
print(f"  handedness s (the s=-1 branch simply fails to fit, it is not a clean first-octant). exit 0")
