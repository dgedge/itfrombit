#!/usr/bin/env python3
r"""The actual "K6" — item 102's 3D velocity-unification, done with the SYMMETRY-CORRECT k.p coupling.
(Self-correction: a naive 6-face exp(i k.dr) model reduces to K6={5,-1x5} at Gamma but gives ZERO first-
order splitting and wrong velocity ratios -- it fails to realise the T_1u<->E_g coupling. The correct
first-order k.p within the -1 manifold is fixed by O_h selection rules and reproduces item 102 exactly.)

Setup: at Gamma the 6x6 Bloch H = complete-graph K6 adjacency J-I, eigenvalues {5 (A_1g), -1x5}; the -1
manifold = T_1u (+) E_g. First-order degenerate k.p couples T_1u<->E_g (the ONLY parity-allowed channel:
T_1u (x) k(T_1u) (x) E_g contains A_1g, since T_1u (x) E_g = T_1u (+) T_2u). The coupling tensor is the
E_g-symmetric bilinear of (k, p):  u ~ (2 k_z p_z - k_x p_x - k_y p_y)/sqrt6,  v ~ (k_x p_x - k_y p_y)/sqrt2.
This is a DERIVATION from O_h, not a fit; the single overall scale lambda is set to 1.
"""
import numpy as np

# ---- (1) Gamma point: complete-graph K6 = J-I, {5, -1x5} ----
J = np.ones((6,6)) - np.eye(6)
w0 = np.sort(np.linalg.eigvalsh(J))
print(f"[Gamma] K6 adjacency J-I eigenvalues = {np.round(w0,4)}  => A_1g(+5), T_1u(+)E_g (-1, x5). Confirmed.")
assert np.allclose(w0, [-1,-1,-1,-1,-1,5])

# ---- (2) the symmetry-fixed first-order k.p on the -1 manifold, basis (px,py,pz | u,v) ----
def Hkp(k):
    kx,ky,kz = k
    H = np.zeros((5,5))
    # u (index 3) <-> p_j
    H[3,0]=-kx/np.sqrt(6); H[3,1]=-ky/np.sqrt(6); H[3,2]=2*kz/np.sqrt(6)
    # v (index 4) <-> p_j
    H[4,0]= kx/np.sqrt(2); H[4,1]=-ky/np.sqrt(2); H[4,2]=0.0
    H = H + H.T                                   # Hermitian (real symmetric)
    return H

def I4(n):
    n=np.array(n,float); n/=np.linalg.norm(n); x,y,z=n
    return x*x*y*y + y*y*z*z + z*z*x*x

print(f"\n[velocities] eigenvalues of the k.p matrix along n_hat (|eig| = bare group velocity; linear in k):")
print(f"  {'dir':>6} {'eigs (=±v, 0...)':>34} {'v_max':>9} {'item102 v_g=sqrt(1/3±(1/3)sqrt(1-3I4))':>40}")
ok=True
for lbl,n in [('[100]',(1,0,0)),('[110]',(1,1,0)),('[111]',(1,1,1))]:
    nh=np.array(n,float); nh/=np.linalg.norm(nh)
    eigs=np.linalg.eigvalsh(Hkp(nh))             # eigenvalues at unit |k| => velocities
    vmax=np.abs(eigs).max()
    v_plus=np.sqrt(1/3 + (1/3)*np.sqrt(max(0,1-3*I4(n))))
    ok &= abs(vmax - v_plus) < 1e-9
    print(f"  {lbl:>6} {np.array2string(np.round(eigs,4)):>34} {vmax:>9.4f} {f'v+={v_plus:.4f}':>40}")
print(f"  targets: sqrt(2/3)={np.sqrt(2/3):.4f}, 1/sqrt2={1/np.sqrt(2):.4f}, 1/sqrt3={1/np.sqrt(3):.4f}")
assert ok, "the k.p velocities must equal item 102's v_+ branch"
print("  => the symmetry-fixed k.p REPRODUCES item 102's velocity formula EXACTLY (v_[100]=sqrt(2/3) etc.).")

# ---- (3) linear vs quadratic: eigenvalues scale LINEARLY with |k| => finite v_g => MARGINAL ----
ks=np.array([0.01,0.02,0.04])
n=np.array([1,0,0],float)
spread=np.array([np.ptp(np.linalg.eigvalsh(Hkp(kk*n))) for kk in ks])
ratio=spread[2]/spread[0]            # 4 if linear (k:0.01->0.04 is x4), 16 if quadratic
print(f"\n[order] splitting along [100] at |k|={ks}: {np.round(spread,5)}")
print(f"        split(4k)/split(k) = {ratio:.2f}  (4 = LINEAR/marginal ; 16 = quadratic/irrelevant)")
linear = abs(ratio-4) < abs(ratio-16)
print(f"        => dispersion is {'LINEAR -> the anisotropy is MARGINAL (a direction-dependent speed of light)' if linear else 'QUADRATIC -> irrelevant'}")
dv0=np.sqrt(2/3)-1/np.sqrt(3)

print(f"""
=========================================================================================
VERDICT (item-102 K6, 3D velocity-unification -- the object the "K6" request meant):
  * Gamma = complete-graph K6 {{5,-1x5}} = A_1g (+) T_1u (+) E_g: confirmed.
  * The first-order k.p is the O_h-FORCED T_1u<->E_g coupling; it reproduces item 102's velocity formula
    v_g(n)=sqrt(1/3 ± (1/3)sqrt(1-3 I4)) EXACTLY: v_[100]=sqrt(2/3), v_[110]=1/sqrt2, v_[111]=1/sqrt3.
  * The splitting is LINEAR in |k| (finite direction-dependent group velocity) -> the bare anisotropy is
    MARGINAL (dimension-4), NOT an O(k^4) irrelevant operator.
  RESOLUTION of the item-102 internal contradiction: in favour of (a) FIRST-ORDER / MARGINAL. The item's
  'Status update' Wilsonian closure -- treating the anisotropy as a dim-6 O(k^4) irrelevant operator that
  flows away by pure power-counting -- is INCONSISTENT with its own first-order velocity formula and is
  WRONG. A marginal velocity anisotropy is NOT removed by power-counting; its IR fate (converge to one c,
  or Lifshitz-diverge) is decided by the INTERACTING RG.
  => Item 115 / the Chadha-Nielsen loop IS the right question, and the sign is genuinely OPEN (per K6).
     The correct 3D bare input is Delta v_0 = v_[100]-v_[111] = {dv0:.3f} (the full 41% anisotropy), NOT the
     2D-8x8 value 0.073. Photon side (3D SC Laplacian 6-2*sum cos, §7.3) is already pinned.
  NEXT (disciplined): graphene known-case harness for the loop machinery, THEN the 3D T_1u<->E_g loop with
  this k.p vertex against the SC photon -- to compute the sign of c_4.8.8.
=========================================================================================""")
assert linear
print("exit 0 -- correct k.p built; item-102 velocity formula reproduced; anisotropy is MARGINAL (loop needed, sign open).")
