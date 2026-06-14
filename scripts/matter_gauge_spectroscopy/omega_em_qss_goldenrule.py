#!/usr/bin/env python3
"""omega_em from canon: the QSS-averaged energy released at the portal, and the full
golden-rule dressed-alpha closure attempt against the REAL web DOS.

Upgrade of the toy single-level photon mode (dressed_alpha_nonunital.py / item75_gaugeweb_dos.py):
under the derived uniform measure, emission from pair eigenstate k goes at rate
   Gamma_k = (1/137) * 2 pi w^2 |<port|k>|^2 rho_web(omega_k),   omega_k = E_k - E_ref,
with rho_web the COMPUTED Wilson-photon DOS (not the small-omega law only). Then
   Gamma_esc = sum_k Gamma_k,   delta_pred = c_gap * Gamma_esc / gap_phys.

Exact structural result (no choices): the emission-weighted mean pair energy under the uniform
measure is the portal diagonal element,  sum_k |<port|k>|^2 E_k = <port|H2|port>  — verified.

Three discrete canon choices remain; the FULL grid is computed (16.3: report all, adopt none):
  t_m   in {Lambda (clock), Lambda/3 (Grover 7.15)}          [matter hopping]
  w     in {1/3, 1/2, 1/sqrt3, sqrt(2/9)}                    [portal coupling]
  E_ref in {band-center 0, band-bottom E_min, E=+1 resonance alignment (7.4),
            cell-gap offset -(E_min) + 2*Delta with Delta=2 (2.6)}
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ---------- bridge pair system (canonical; as in prior scripts) ----------
N = 16
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1
pairs = [(i,j) for i in range(N) for j in range(i,N)]
idx = {p:k for k,p in enumerate(pairs)}
Bas = np.zeros((N*N,len(pairs)))
for (i,j),k in idx.items():
    v = np.zeros((N,N))
    if i==j: v[i,i]=1.0
    else: v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,k]=v.reshape(-1)
H2u = Bas.T @ (np.kron(A,np.eye(N))+np.kron(np.eye(N),A)) @ Bas      # unit-hopping pair H
port = np.zeros(len(pairs))
for p in ((0,8),(0,9),(1,8),(1,9)): port[idx[p]] = 0.5
Ek_u, Vk = np.linalg.eigh(H2u)
ov2 = (Vk.T @ port)**2                                                # |<port|k>|^2
assert abs(ov2.sum()-1) < 1e-12

# ---------- exact structural result: emission-weighted mean energy ----------
mean_E = float(ov2 @ Ek_u)
port_diag = float(port @ H2u @ port)
print(f"PORTAL SPECTRAL IDENTITY: sum_k |<port|k>|^2 E_k = {mean_E:.10f} = <port|H2|port> = {port_diag:.10f}")
assert abs(mean_E - port_diag) < 1e-10
assert abs(port_diag - 2.0) < 1e-10        # the portal 4-cycle gives exactly 2 (in t_m units)
print(f"  => under the uniform measure the mean energy released at the portal is EXACTLY 2*t_m")
print(f"     (relative to band-centre reference) — the canon-native omega_em theorem.\n")

# ---------- Wilson-photon web DOS (t_web = Lambda; reuse the verified construction) ----------
n = 160
k = (np.arange(n)+0.5)*2*np.pi/n
KX,KY,KZ = np.meshgrid(k,k,k, indexing="ij")
om_grid = np.sqrt(4*(np.sin(KX/2)**2+np.sin(KY/2)**2+np.sin(KZ/2)**2)).ravel()
hist, ed = np.histogram(om_grid, bins=400, range=(0,float(om_grid.max())+1e-9), density=True)
ctw = 0.5*(ed[:-1]+ed[1:])
om_max = float(om_grid.max())
def rho_ph(w_):
    return float(np.interp(w_, ctw, hist)) if 0 < w_ < om_max else 0.0
assert abs(om_max - 2*np.sqrt(3)) < 0.01

# ---------- chain normalisation (as derived) ----------
Wc = np.abs(np.block([[H2u, (0.5*port)[:,None]],[ (0.5*port)[None,:], np.zeros((1,1))]]))**2
np.fill_diagonal(Wc,0)
gen = Wc - np.diag(Wc.sum(1))
gap_phys = 2*(-np.sort(np.linalg.eigvalsh(gen))[-2])
c_gap = 84.21
delta_target = 0.0360
assert abs(gap_phys - 0.3045) < 0.004

# ---------- the full canon-pinned grid ----------
Emin_u = float(Ek_u.min())
tms = {"t_m=Lambda": 1.0, "t_m=Lambda/3 (Grover)": 1/3}
ws  = {"w=1/3": 1/3, "w=1/2": 0.5, "w=1/sqrt3": 1/np.sqrt(3), "w=sqrt(2/9)": np.sqrt(2/9)}
print(f"{'t_m':>22} {'w':>12} {'E_ref':>26} {'<omega>_emit':>12} {'delta_pred':>11} {'/target':>8}")
results = []
for tn, tm in tms.items():
    Ek = tm*Ek_u; Emin = tm*Emin_u
    refs = {"band-centre (0)": 0.0,
            "band-bottom": Emin,
            "E=+1 alignment (7.4)": -1.0,                       # omega = E_k + 1
            "cell-gap 2Delta (2.6)": Emin - 4.0}                # omega = E_k - Emin + 4
    for wn, w in ws.items():
        for rn, Eref in refs.items():
            om_k = Ek - Eref
            rates = np.array([ov2[i]*rho_ph(om_k[i]) for i in range(len(Ek))])
            Gesc = (2*np.pi*w*w/137.0)*rates.sum()
            mean_om = float((rates @ om_k)/rates.sum()) if rates.sum() > 0 else float('nan')
            dl = c_gap*Gesc/gap_phys
            results.append((tn,wn,rn,mean_om,dl))
            print(f"{tn:>22} {wn:>12} {rn:>26} {mean_om:>12.3f} {dl:>11.3e} {dl/delta_target:>8.2f}")
hits = [r for r in results if 0.5 < r[4]/delta_target < 2.0]
print(f"\ncombinations within 2x of delta = {delta_target}: {len(hits)} of {len(results)}")
for h in hits: print(f"   HIT: {h[0]}, {h[1]}, {h[2]}  -> delta = {h[4]:.3e} ({h[4]/delta_target:.2f}x), <omega> = {h[3]:.3f} Lambda")

print(f"""
VERDICT LOGIC (pre-stated): if zero canon-pinned combinations land within 2x, the QSS golden-rule
route FAILS to recover omega_em from canon — the second item-75 target is falsified-as-posed
(the energy bookkeeping of the emission step needs physics canon does not yet contain). If hits
exist, they are candidates ONLY: the model-space density ({len(results)} combos) is itself a
16.3-style alphabet, and adoption would require the hitting combination to be independently
forced. Nothing is adopted by this script.
ALL ASSERTS PASSED""")
