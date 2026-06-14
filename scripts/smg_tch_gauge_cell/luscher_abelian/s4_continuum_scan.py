#!/usr/bin/env python3
"""S4 — continuum-limit (a->0) scan of the abelian anomaly-cancellation residual.
Berry-phase machinery of S3 with a BAND-LIMITED NON-SEPARABLE smooth U(1) background
(fixed low wavevectors, L-independent coeffs => same physical field sampled finer => O(a^2)
artifacts ~1/L^2). Memory-light sparse build of H_W (one dense materialisation for eigh);
P_hat_- basis = positive-H_W eigvecs (one eigh per config). Reports spread + cancel per L."""
import sys, itertools, numpy as np, scipy.sparse as sp
np.seterr(divide="ignore", invalid="ignore")
try:
    import cupy as xp                         # GPU: ~16x faster FP64 eigh on the box's 3090
    GPU = True
except Exception:
    xp = np; GPU = False                      # CPU fallback (Mac)
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
amp_base = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
amp_dlt  = float(sys.argv[3]) if len(sys.argv) > 3 else 0.06
amp_dlt *= (4.0/L)**2          # scale perturbation ~1/L^2: Grassmann-loop Sum(theta^2) ~ dim ~ L^4,
                               # so this holds |loop| (conditioning) constant across L. arg is the L=4 value.
M = 1.0
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.array([[1,0],[0,-1]],complex)]
I2 = np.eye(2,dtype=complex); Zr = np.zeros((2,2),complex)
g = [np.block([[Zr,-1j*s[k]],[1j*s[k],Zr]]) for k in range(3)] + [np.block([[I2,Zr],[Zr,-I2]])]
g5 = g[0]@g[1]@g[2]@g[3]
gA = [g5@(np.eye(4)-g[mu]) for mu in range(4)]              # g5(I-gamma_mu) for forward hop
gB = [g5@(np.eye(4)+g[mu]) for mu in range(4)]              # g5(I+gamma_mu) for backward hop
sites = list(itertools.product(range(L), repeat=4)); N = len(sites); idx = {c:i for i,c in enumerate(sites)}
coord = np.array(sites, float)
KS = []; seen = set()                                       # fixed non-separable low wavevectors
for k in itertools.product((-1,0,1), repeat=4):
    if not any(k) or tuple(-np.array(k)) in seen: continue
    seen.add(k); KS.append(np.array(k, float))
def shift_sp(mu):                                          # forward shift permutation (sparse)
    rows = [idx[tuple((np.array(c)+np.eye(4,dtype=int)[mu]) % L)] for c in sites]
    return sp.csr_matrix((np.ones(N), (rows, range(N))), shape=(N, N), dtype=complex)
SHIFTP = {mu: shift_sp(mu) for mu in range(4)}
G5c = sp.kron(sp.identity(N, dtype=complex), g5, format="csr")
def smooth(seed, amp):
    r = np.random.default_rng(seed)
    c = r.uniform(-1, 1, (4, len(KS))); phi = r.uniform(0, 2*np.pi, (4, len(KS)))
    th = []
    for mu in range(4):
        f = np.zeros(N)
        for ik, k in enumerate(KS):
            f += c[mu, ik]*np.cos(2*np.pi*(coord @ k)/L + phi[mu, ik])
        th.append(amp*f/np.sqrt(len(KS)))
    return th
A0 = smooth(1, amp_base); d1 = smooth(2, amp_dlt); d2 = smooth(3, amp_dlt)
configs = [A0, [A0[mu]+d1[mu] for mu in range(4)], [A0[mu]+d2[mu] for mu in range(4)]]
def subspace(theta):
    HW = (4-M)*G5c.copy()
    for mu in range(4):
        Tp = SHIFTP[mu] @ sp.diags(np.exp(1j*theta[mu]))   # scale forward link by phase
        Tm = Tp.conj().T
        HW = HW - 0.5*(sp.kron(Tp, gA[mu], format="csr") + sp.kron(Tm, gB[mu], format="csr"))
    HWd = HW.toarray()
    if GPU: HWd = xp.asarray(HWd)
    w, V = xp.linalg.eigh(HWd)                             # one dense eigh (GPU if available)
    pos = w > 0
    index = int(round(0.5*(int(xp.count_nonzero(~pos)) - int(xp.count_nonzero(pos)))))
    W = V[:, pos].copy()                                   # P_hat_- basis; copy so V can be freed
    mh = float(xp.min(xp.abs(w)))
    if GPU:
        del HWd, w, V, pos
        xp.get_default_memory_pool().free_all_blocks()     # keep GPU memory bounded for L>=8
    return W, index, mh
def berry(Y):
    Ws, ixs, mhs = [], [], []
    for cf in configs:
        W, ix, mh = subspace([Y*cf[mu] for mu in range(4)]); Ws.append(W); ixs.append(ix); mhs.append(mh)
    if len(set(ixs)) != 1 or len(set(W.shape[1] for W in Ws)) != 1:
        return None, ixs, min(mhs), 0.0
    ph, la = 1.0+0j, 0.0
    for a, b in ((0,1),(1,2),(2,0)):
        sg, l = xp.linalg.slogdet(Ws[a].conj().T@Ws[b]); ph *= complex(sg); la += float(l)
    if GPU: xp.get_default_memory_pool().free_all_blocks()
    return float(-np.angle(ph)), ixs[0], min(mhs), float(np.exp(la))
charges = [(1/6,6),(-2/3,3),(1/3,3),(-1/2,2),(1.0,1),(0.0,1)]
print(f"S4 L={L} (dim {N*4}) amp_base={amp_base} amp_dlt={amp_dlt} |KS|={len(KS)} backend={'GPU' if GPU else 'CPU'}", flush=True)
tot, tot_abs, ratios = 0.0, 0.0, []
for Y, mult in charges:
    b, ix, mh, mod = berry(Y)
    if b is None: print(f"  Y={Y:+.4f} x{mult}: SKIP index/rank {ix}", flush=True); continue
    if abs(Y) > 1e-12: ratios.append(b/Y**3)
    tot += mult*b; tot_abs += mult*abs(b)
    print(f"  Y={Y:+.4f} x{mult}: Berry={b:+.3e} Berry/Y^3={(b/Y**3 if abs(Y)>1e-12 else float('nan')):+.3e} "
          f"idx={ix} min|H_W|={mh:.3f} |loop|={mod:.3f}", flush=True)
rr = np.array(ratios)
spread = rr.std()/abs(rr.mean()); cancel = abs(tot)/tot_abs
print(f"RESULT L={L}: spread={spread:.4e}  cancel={cancel:.4e}  scale={tot_abs:.3e}  (a~1/L={1/L:.4f})", flush=True)
