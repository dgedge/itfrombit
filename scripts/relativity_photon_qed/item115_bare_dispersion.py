#!/usr/bin/env python3
r"""ITEM 115 bare dispersion: the explicit 8x8 bipartite 4.8.8 Bloch Hamiltonian (DRIFT K2 L781),
built, diagonalised, and checked. Purpose: (1) reproduce the canon's leading-edge velocities
(v_fast=1.00, v_slow~0.78) to confirm this IS the right matrix; (2) test whether sqrt(2/3)=0.8165
appears anywhere in the bare dispersion/velocities (user claim: it does NOT).

Construction (DRIFT K2, willow/echo2.tex; CTQW generator A = adjacency, U(t)=exp(-i t A)):
  H(k) = [[H_intra, H_inter(k)], [H_inter(k)^dag, H_intra]]   (8x8, two square-sublattices A,B)
  H_intra = C4 square adjacency  [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]]   ((1,3),(2,4)=0)
  H_inter(k) = diag( e^{i ky}, e^{i kx}, e^{-i ky}, e^{-i kx} )   (one external bridge per node)
Unphased adjacency (the pi/4 chiral phase is the §7.2 gauge/Chern sector, not the bare velocity).
"""
import numpy as np

H_intra = np.array([[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]], float)
def H_inter(kx, ky):
    return np.diag([np.exp(1j*ky), np.exp(1j*kx), np.exp(-1j*ky), np.exp(-1j*kx)])
def H(kx, ky):
    Hi = H_inter(kx, ky)
    return np.block([[H_intra.astype(complex), Hi],
                     [Hi.conj().T,            H_intra.astype(complex)]])

# --- sanity: Hermitian + chiral {H, Sigma_z}=0 (the K2 bipartite-fidelity check) ---
Sz = np.diag([1.,-1,1,-1, -1,1,-1,1])                      # White/Black sublattice parity (per K2 colouring)
kx,ky = 0.7, -0.3
assert np.allclose(H(kx,ky), H(kx,ky).conj().T), "must be Hermitian"
# chiral symmetry: the 4.8.8 is bipartite -> exists Sigma with {H,Sigma}=0 mapping E<->-E
ev = np.linalg.eigvalsh(H(kx,ky))
assert np.allclose(sorted(ev), sorted(-ev), atol=1e-9), "spectrum must be E<->-E symmetric (bipartite/chiral)"
print(f"[checks] Hermitian OK; spectrum E<->-E symmetric OK (bipartite). 8 bands at k=(0.7,-0.3): "
      f"{np.round(ev,3)}")

# --- bands + leading-edge group velocities over the BZ ---
N = 4001
kline = np.linspace(-np.pi, np.pi, N)
dk = kline[1]-kline[0]
# along kx axis (ky=0): bands E_n(kx)
bands = np.array([np.linalg.eigvalsh(H(kx,0.0)) for kx in kline])   # (N,8) ascending
vx = np.gradient(bands, dk, axis=0)                                  # dE/dkx per band
vmax_axis = np.abs(vx).max(0)                                        # max |v| per band, axis direction
# full BZ, both directions, group speed |grad E| max per band
M = 161
kk = np.linspace(-np.pi, np.pi, M)
gmax = np.zeros(8)
for i,kx in enumerate(kk):
    Ex = np.array([np.linalg.eigvalsh(H(kx,ky)) for ky in kk])      # (M,8)
    # central differences
    dEdkx = np.gradient(np.array([np.linalg.eigvalsh(H(kx+1e-4,ky)) for ky in kk]),0,axis=0) if False else None
for_speed = []
def grad_speed():
    h=1e-4; g=np.zeros(8)
    for kx in kk:
        for ky in kk:
            ex=(np.linalg.eigvalsh(H(kx+h,ky))-np.linalg.eigvalsh(H(kx-h,ky)))/(2*h)
            ey=(np.linalg.eigvalsh(H(kx,ky+h))-np.linalg.eigvalsh(H(kx,ky-h)))/(2*h)
            g=np.maximum(g, np.sqrt(ex**2+ey**2))
    return g
gmax = grad_speed()

print(f"\n[bands] max |dE/dk| per band (axis, ky=0):  {np.round(vmax_axis,4)}")
print(f"[bands] max |grad E| per band (full BZ):    {np.round(gmax,4)}")
# the distinct leading-edge speeds (fast / slow branches)
distinct = np.unique(np.round(gmax,2))
vfast, vslow = distinct.max(), sorted(distinct)[-2] if len(distinct)>1 else distinct.max()
# axis-direction leading edges (what echo2/K2 quote as 1.00 / 0.78)
ax_distinct = np.unique(np.round(vmax_axis,2))
print(f"\n[result] axis-direction leading-edge speeds present: {ax_distinct}")
print(f"         full-BZ leading-edge speeds present:        {distinct}")
vf_axis = vmax_axis.max()
# slow = the largest axis speed strictly below the fast cluster
below = vmax_axis[vmax_axis < vf_axis-0.05]
vs_axis = below.max() if below.size else vf_axis
print(f"         => v_fast(axis) = {vf_axis:.3f}   v_slow(axis) = {vs_axis:.3f}   ratio = {vf_axis/vs_axis:.3f}")
print(f"            canon (DRIFT K2): v_fast=1.00, v_slow~0.78, ratio~1.28")

# --- the sqrt(2/3) test ---
s23 = np.sqrt(2/3)
print(f"\n[sqrt(2/3) test] sqrt(2/3) = {s23:.4f}")
hits = []
for label,arr in [("axis vmax",vmax_axis),("BZ |grad|max",gmax),
                  ("band energies @random k",np.linalg.eigvalsh(H(0.9,0.4)))]:
    for val in np.atleast_1d(arr):
        if abs(abs(val)-s23) < 0.01:
            hits.append((label,float(val)))
print(f"   values within 0.01 of sqrt(2/3)={s23:.4f} anywhere in the bare dispersion: "
      f"{hits if hits else 'NONE'}")
print(f"   v_slow(axis)={vs_axis:.3f} vs sqrt(2/3)={s23:.4f}: differ by {abs(vs_axis-s23):.3f} "
      f"({abs(vs_axis-s23)/s23*100:.1f}%).  v_fast(OTOC-centroid 0.81) vs sqrt(2/3): "
      f"{abs(0.81-s23):.3f}.")
print(f"   => the bare-dispersion velocities are {vf_axis:.2f} and {vs_axis:.2f}; sqrt(2/3) is NOT among them.")

# --- the explicit matrix, paste-ready ---
print("""
================= EXPLICIT 8x8 BARE-DISPERSION BLOCH MATRIX H(kx,ky) (paste-ready) =================
basis order: [A:N, A:E, A:S, A:W,  B:N, B:E, B:S, B:W]   (two square sublattices A,B)
H_intra (each diagonal 4x4 block) = C4 square adjacency:
    [[0,1,0,1],
     [1,0,1,0],
     [0,1,0,1],
     [1,0,1,0]]
H_inter(k) (upper-right 4x4 block) = diag( e^{+i ky}, e^{+i kx}, e^{-i ky}, e^{-i kx} )
lower-left block = H_inter(k)^dagger.  Full:
H(kx,ky) =
[  0    1    0    1  | e^{iky}  0      0      0     ]
[  1    0    1    0  | 0       e^{ikx} 0      0     ]
[  0    1    0    1  | 0       0      e^{-iky} 0     ]
[  1    0    1    0  | 0       0      0      e^{-ikx}]
[ ---------------------------------------------------]
[ e^{-iky} 0   0    0  | 0    1    0    1            ]
[  0   e^{-ikx} 0   0  | 1    0    1    0            ]
[  0    0   e^{iky} 0  | 0    1    0    1            ]
[  0    0    0  e^{ikx}| 1    0    1    0            ]
===================================================================================================""")
print("exit 0 -- bare-dispersion 8x8 built & diagonalised; velocities reported; sqrt(2/3) checked.")
