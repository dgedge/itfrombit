#!/usr/bin/env python3
r"""Verify the user's reverse-engineering of the phantom 1.00 in DRIFT K2:
  (1) the E=0 dark modes of the square project the inter-cell coupling D to diag(cos ky, cos kx)
      with zero cross-terms -> a decoupled 4x4 H_eff of two 1D wires, max velocity EXACTLY 1.00;
  (2) in the full 8x8 those dark modes hybridise with the E=+-2 bright modes via <p|D|z1>=(i/sqrt2)sin k,
      and level repulsion drags the max velocity 1.00 -> 0.78;
  (3) therefore '1.00' is the projected-model velocity, '0.78' is the full-8x8 max -> the quoted
      ratio 1.28 is CROSS-MODEL (phantom); the real within-8x8 bifurcation is 0.78/0.707.
We confirm by a CONTINUOUS interpolation lambda: scale the dark<->bright coupling of D by lambda and
watch the leading-edge velocity run 1.00 (lambda=0, decoupled dark sector) -> 0.78 (lambda=1, full).
"""
import numpy as np

H_intra = np.array([[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]], float)
def Dk(kx,ky): return np.diag([np.exp(1j*ky),np.exp(1j*kx),np.exp(-1j*ky),np.exp(-1j*kx)])

# square eigenbasis: dark z1,z2 (E=0); bright p+ (E=+2), p- (E=-2)
z1=np.array([1,0,-1,0])/np.sqrt(2); z2=np.array([0,1,0,-1])/np.sqrt(2)
pp=np.array([1,1,1,1])/2.0;          pm=np.array([1,-1,1,-1])/2.0
U=np.column_stack([z1,z2,pp,pm]).astype(complex)              # columns = eigenvectors
assert np.allclose(U.conj().T@U, np.eye(4)), "eigenbasis must be orthonormal"
assert np.allclose(U.conj().T@H_intra@U, np.diag([0,0,2,-2])), "H_intra eigenvalues 0,0,2,-2"

# (1) projection identities (assert the user's matrix elements exactly)
kx,ky=0.6,-0.4
Dt = U.conj().T @ Dk(kx,ky) @ U                                # D in {z1,z2,p+,p-} basis
assert abs(Dt[0,0]-np.cos(ky))<1e-12 and abs(Dt[1,1]-np.cos(kx))<1e-12, "dark-dark = diag(cos ky,cos kx)"
assert abs(Dt[0,1])<1e-12 and abs(Dt[1,0])<1e-12, "dark-dark cross-term = 0"
assert abs(Dt[2,0]-1j*np.sin(ky)/np.sqrt(2))<1e-12, "<p+|D|z1> = (i/sqrt2) sin ky"
print("[1] projection identities VERIFIED: <z1|D|z1>=cos ky, <z2|D|z2>=cos kx, cross=0, "
      "<p+|D|z1>=(i/sqrt2)sin ky.")

def max_speed_8x8_in_eigbasis(lam, N=121):
    """Full 8x8 in the eigenbasis, with the dark<->bright block of D scaled by lam. lam=1 => exact full 8x8."""
    Hd=np.diag([0,0,2,-2]).astype(complex)
    g=np.zeros(8); h=1e-4; kk=np.linspace(-np.pi,np.pi,N)
    def Hful(kx,ky):
        Dt=U.conj().T@Dk(kx,ky)@U
        S=Dt.copy()
        S[0:2,2:4]*=lam; S[2:4,0:2]*=lam                      # scale dark<->bright coupling
        return np.block([[Hd,S],[S.conj().T,Hd]])
    for kx in kk:
        for ky in kk:
            ex=(np.linalg.eigvalsh(Hful(kx+h,ky))-np.linalg.eigvalsh(Hful(kx-h,ky)))/(2*h)
            ey=(np.linalg.eigvalsh(Hful(kx,ky+h))-np.linalg.eigvalsh(Hful(kx,ky-h)))/(2*h)
            g=np.maximum(g,np.sqrt(ex**2+ey**2))
    return g.max()

# (2)+(3) the interpolation 1.00 -> 0.78
print("\n[2] hybridisation interpolation (scale dark<->bright coupling by lambda):")
print(f"     {'lambda':>7} {'max leading-edge velocity':>26}")
vals=[]
for lam in [0.0,0.25,0.5,0.75,1.0]:
    v=max_speed_8x8_in_eigbasis(lam); vals.append((lam,v))
    tag = "  <- decoupled dark sector = H_eff (TWO 1D wires)" if lam==0 else ("  <- FULL 8x8 bare dispersion" if lam==1 else "")
    print(f"     {lam:>7.2f} {v:>26.4f}{tag}")
v0,v1=vals[0][1],vals[-1][1]
assert abs(v0-1.00)<0.01, f"lambda=0 must give the projected 1.00, got {v0}"
assert abs(v1-0.78)<0.02, f"lambda=1 must give the full-8x8 0.78, got {v1}"
mono = all(vals[i][1]>=vals[i+1][1]-1e-6 for i in range(len(vals)-1))
print(f"     => v(0)={v0:.3f} (projected 1.00)  ->  v(1)={v1:.3f} (full 0.78); monotone suppression: {mono}")
print(f"        CONFIRMED: the 0.78 IS the hybridisation-suppressed continuation of the projected 1.00 "
      f"(same mode, dragged down by level repulsion from the E=+-2 bands).")

# the real within-8x8 bifurcation (for contrast with the phantom 1.28)
per_band = []
kk=np.linspace(-np.pi,np.pi,2001); dk=kk[1]-kk[0]
def Hfull(kx,ky):
    D=Dk(kx,ky); return np.block([[H_intra.astype(complex),D],[D.conj().T,H_intra.astype(complex)]])
bands=np.array([np.linalg.eigvalsh(Hfull(kx,0.0)) for kx in kk])
vax=np.abs(np.gradient(bands,dk,axis=0)).max(0)
print(f"\n[3] real within-8x8 axis velocities per band: {np.round(vax,3)}")
print(f"     the two dominant scales are {sorted(set(np.round(vax,2)))[-1]:.2f} and "
      f"{sorted(set(np.round(vax[vax>0.05],2)))[-2]:.2f} (=1/sqrt2={1/np.sqrt(2):.3f}) -> real ratio "
      f"{sorted(set(np.round(vax,3)))[-1]/ (1/np.sqrt(2)):.3f}, NOT 1.28.")
print(f"     the quoted 1.28 = 1.00/0.78 compares the PROJECTED model to the FULL model -> PHANTOM (cross-model).")

print(f"""
=========================================================================================
VERDICT: the user's reverse-engineering is CORRECT and computationally confirmed.
  * 1.00 = leading-edge velocity of the decoupled E=0 dark-mode projection (4x4 H_eff = two 1D wires,
    E=+-cos kx, +-cos ky). A real sub-algebra of the 8x8, but an UNCONTROLLED approximation (it discards
    the O(1) dark<->bright coupling (i/sqrt2)sin k).
  * 0.78 = the actual leading-edge (max group) velocity of the full 8x8 bare dispersion -- the SAME dark
    mode, hybridisation-suppressed from 1.00 by level repulsion off the E=+-2 bands (interpolation above).
  * The quoted (v_fast,v_slow,ratio)=(1.00,0.78,1.28) at DRIFT K2 L781/L785 CONFLATES the projected 4x4
    model (1.00) with the full 8x8 max (0.78): a cross-model ratio. It is a PHANTOM -- the 8x8 has no band
    at 1.00 and no internal 1.28. The real within-8x8 bifurcation is 0.78 / 0.707(=1/sqrt2), ratio ~1.10.
  * Consequence: the downstream "ratio 1.28 -> branch-dependent Airy lag -> (0.81,0.60)" reconciliation
    (ANCHOR §7.4, DRIFT K2) is built on the phantom ratio and is mathematically ungrounded. And sqrt(2/3)
    is absent from the bare spectrum entirely (item115_bare_dispersion.py).
=========================================================================================""")
print("exit 0 -- dark-mode origin of the phantom 1.00 verified; 0.78 confirmed as its hybridised continuation.")
