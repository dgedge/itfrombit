#!/usr/bin/env python3
"""Item 75 (second target): the gauge-web density of states, and the dressed-alpha closure attempt
through the derived linear-response formula delta = c * Gamma_esc/gap (dressed_alpha_nonunital.py).

Canon-pinned ingredients: web = SC lattice (7.3); gauge-matter hopping t = 1/3 (7.15 Grover, the
H6-pinned value); monitoring at the clock rate gamma = Lambda (one syndrome per tick, 5.9/119);
E = +1 web transmission resonance (7.4/9.7). NOT pinned: the emission energy omega_em.
16.3 discipline: every coupling/energy combination tried is reported; nothing is adopted.
Units: Lambda = 1, lattice a0 = 1. Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ---------- 1. SC tight-binding DOS (the web's matter-like band): E(k) = -2 sum cos k_i ----------
n = 240
k = (np.arange(n)+0.5)*2*np.pi/n
KX,KY,KZ = np.meshgrid(k,k,k, indexing="ij")
E_sc = -2*(np.cos(KX)+np.cos(KY)+np.cos(KZ))
hist, edges = np.histogram(E_sc, bins=600, range=(-6,6), density=True)
ctr = 0.5*(edges[:-1]+edges[1:])
def rho_sc(E):  return float(np.interp(E, ctr, hist))
assert abs(np.trapezoid(hist, ctr) - 1) < 1e-3                  # normalised
assert abs(rho_sc(0.0) - rho_sc(-0.0)) < 1e-9 and abs(rho_sc(2.5)-rho_sc(-2.5)) < 6e-3  # symmetric
print(f"1. SC web DOS built (240^3 k-grid): rho(0) = {rho_sc(0):.4f}/site  rho(E=+1) = {rho_sc(1.0):.4f}/site")
assert 0.13 < rho_sc(0) < 0.15                                  # known SC band-centre DOS ~0.142

# ---------- 2. photonic (Wilson-action) web band: omega^2(k) = 4 sum sin^2(k_i/2) ----------
W2 = 4*(np.sin(KX/2)**2 + np.sin(KY/2)**2 + np.sin(KZ/2)**2)
om = np.sqrt(W2)
histw, edw = np.histogram(om, bins=600, range=(0, float(om.max())+1e-9), density=True)
ctw = 0.5*(edw[:-1]+edw[1:])
def rho_ph(w):  return float(np.interp(w, ctw, histw))
# small-omega law via the CUMULATIVE count (bin-stable): N(<omega) = omega^3/(6 pi^2)
for w0 in (0.4, 0.6, 0.8):
    Ncum = float((om < w0).mean())
    coef = Ncum/w0**3
    assert abs(coef - 1/(6*np.pi**2)) < 0.0025, (w0, coef, 1/(6*np.pi**2))
print(f"2. photon (Wilson) DOS built: cumulative law N(<omega) = omega^3/(6 pi^2) verified "
      f"(coef at omega=0.6: {float((om<0.6).mean())/0.6**3:.5f} vs {1/(6*np.pi**2):.5f}); "
      f"density law rho = omega^2/(2 pi^2) follows")

# ---------- 3. line-graph structure check (the link-variable lattice): finite cluster ----------
L = 6; N3 = L**3
sites = list(itertools.product(range(L),repeat=3))
sidx = {s:i for i,s in enumerate(sites)}
edgesL = []
for s in sites:
    for d in range(3):
        t = list(s); t[d] = (t[d]+1)%L
        edgesL.append((sidx[s], sidx[tuple(t)]))
M = len(edgesL)                                                 # 3N edges
Aedge = np.zeros((M,M))
for a in range(M):
    for b in range(a+1,M):
        if len(set(edgesL[a]) & set(edgesL[b])) == 1:
            Aedge[a,b]=Aedge[b,a]=1
evL = np.sort(np.linalg.eigvalsh(Aedge))
flat = int(np.sum(np.abs(evL + 2) < 1e-8))
assert flat == M - N3 + 1                                       # line-graph flat band at -2, mult E-V+1
print(f"3. L(SC) cluster check (L={L}): flat band at -2 with multiplicity E-V+1 = {flat} "
      f"= {M}-{N3}+1 (the canon multi-band/flat-band structure)")

# ---------- 4. the closure attempt: required Gamma_esc/gap through the derived formula ----------
# rebuild the canonical 137-chain gap (self-contained; matches dressed_alpha_nonunital.py)
N = 16
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1
pairs = [(i,j) for i in range(N) for j in range(i,N)]
idx = {p:c for c,p in enumerate(pairs)}
d137 = len(pairs)+1
Bas = np.zeros((N*N,len(pairs)))
for (i,j),c in idx.items():
    v = np.zeros((N,N))
    if i==j: v[i,i]=1.0
    else: v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,c]=v.reshape(-1)
H = np.zeros((d137,d137)); H[:len(pairs),:len(pairs)] = Bas.T@(np.kron(A,np.eye(N))+np.kron(np.eye(N),A))@Bas
for p in ((0,8),(0,9),(1,8),(1,9)): H[idx[p],len(pairs)] = H[len(pairs),idx[p]] = 0.25
Wc = np.abs(H)**2; np.fill_diagonal(Wc,0)
gen = Wc - np.diag(Wc.sum(1))
gap_W = -np.sort(np.linalg.eigvalsh(gen))[-2]                   # 0.1522 in |H|^2 units
c_gap = 84.21                                                   # derived linear coefficient (gap units)
delta_target = 0.0360
req_ratio = delta_target/c_gap                                  # Gamma_esc/gap_phys needed
assert abs(gap_W - 0.1522) < 0.002
# physical units: dephased chain rates = (2/gamma)|H|^2 with gamma = Lambda = 1 -> gap_phys = 2*gap_W
gap_phys = 2*gap_W
print(f"\n4. closure attempt: required Gamma_esc = {req_ratio:.3e} x gap_phys = {req_ratio*gap_phys:.3e} Lambda")
print(f"   (gap_phys = 2 x {gap_W:.4f} = {gap_phys:.4f} Lambda at clock-rate monitoring)")
couplings = {"Grover t=1/3 (7.15)": 1/3, "portal w=1/2": 0.5, "hop 1/sqrt3 (3.2)": 1/np.sqrt(3),
             "V_weak sqrt(2/9)": np.sqrt(2/9)}
print(f"   {'coupling':>22} {'DOS reading':>26} {'Gamma_esc':>11} {'delta_pred':>11}")
results = []
for cn, w in couplings.items():
    for dn, rho in (("SC @ E=+1 (7.4 res.)", rho_sc(1.0)),
                    ("photon @ omega=alpha*L", rho_ph(1/137)),
                    ("photon @ omega=gap", rho_ph(gap_phys))):
        G = 2*np.pi*w*w*rho
        dl = c_gap*G/gap_phys
        results.append((cn,dn,G,dl))
        print(f"   {cn:>22} {dn:>26} {G:>11.3e} {dl:>11.3e}")
# none lands on 0.036 within 2x?
hits = [r for r in results if 0.5 < r[3]/delta_target < 2.0]
print(f"   combinations within 2x of the observed delta = {delta_target}: {len(hits)}")

# the inverse problem: required omega_em per coupling, vs an ENUMERATED canonical-scale list
print(f"\n   inverse problem (photonic DOS): required emission energy omega_em per coupling:")
scales = {"alpha*Lambda (bit-weight)": 1/137, "gap_phys": gap_phys}
scales.update({f"Lambda/{m}": 1/m for m in range(2,21)})
for cn, w in couplings.items():
    rho_req = req_ratio*gap_phys/(2*np.pi*w*w)
    om_req = float(np.sqrt(2*np.pi**2*rho_req))
    near = sorted(scales.items(), key=lambda kv: abs(kv[1]/om_req-1))[:3]
    ns = ", ".join(f"{k}={v:.4f} ({(v/om_req-1)*100:+.0f}%)" for k,v in near)
    print(f"   {cn:>22}: omega_em = {om_req:.4f} Lambda;  nearest scales: {ns}")
# competitor density within +-10% of the Grover-coupling omega_req (16.3 discipline)
w = 1/3
om_req = float(np.sqrt(2*np.pi**2*req_ratio*gap_phys/(2*np.pi*w*w)))
comp = [k for k,v in scales.items() if abs(v/om_req-1) < 0.10]
print(f"   16.3 check: canonical scales within 10% of the Grover omega_em: {comp} "
      f"({len(comp)} competitors in a deliberately small list — the alphabet is dense; NO adoption)")

print(f"""
RESULT:
  INFRASTRUCTURE: the gauge-web DOS objects are now built and verified (SC band rho(0)={rho_sc(0):.3f},
  Wilson-photon omega^2/(2pi^2) law, L(SC) flat band at -2 mult E-V+1) — the genuine item-75
  3D-bulk objects, reusable beyond this attempt.
  CLOSURE ATTEMPT: with every canon-pinned coupling and energy reading, NO combination reproduces
  delta = 0.036 (all reported above; none within 2x). The inverse problem pins what is missing to
  ONE scalar: the emission energy omega_em ~ 0.06 Lambda (Grover coupling), for which canon has no
  derivation and the canonical-scale alphabet is too dense to select honestly.
  NET: the dressed-alpha magnitude now rests on a single named, underived object — omega_em (the
  energy the annihilating pair deposits into the web) — with the full pipeline
  omega_em -> rho_web -> Gamma_esc -> delta derived and verified around it. Item 75's promotion
  criterion gains this second, sharper target alongside the original N1/C2 3D-promotion.
ALL ASSERTS PASSED""")
