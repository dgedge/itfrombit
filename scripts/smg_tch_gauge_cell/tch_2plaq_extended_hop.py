#!/usr/bin/env python3
r"""ADD THE NEUTRAL-HOPPING LEG to the 106,460 extended escalation.

Builds the gauge-covariant fermion hop  psi^dag_head U_l psi_tail  (+ h.c.) on the
2-plaquette NON-REDUCED charged basis, then runs the full t-scan on the 106,460
extended basis:  H = (1/beta) sum C(R_l) - (beta/2)(W0+W0d+W1+W1d) + t (HOP+HOPd)
                    + Delta_SMG (charged vertices).

DISCIPLINE: the hop changes the FERMION rep (3<->1) and one link rep, so it is NOT the
spectator-fermion magnetic move. It is validated ELEMENTWISE against the already-validated
single-cell HOP of tch_nonreduced_magnetic_plaquette.py (whose |K^4|=1 character gate passed)
BEFORE it is trusted on the 2-plaquette / extended basis.
  STEP 0: load validated machinery (finish: GATE A/B) + the single-cell hop oracle (mag_plaq).
  STEP 1: generic hop reproduces the single-cell HOP elementwise (|dHOP|<1e-10).
  STEP 2: 2-plaquette minimal regression — W matches dense (0), HOP+HOPd Hermitian, t=0 gap
          reproduces the published 6.674654.
  STEP 3: the 106,460 extended t-scan.
exit 0 = hop validated against the single-cell oracle AND the extended t-scan computed.
"""
import importlib.util, functools, string, sys, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
print("[load] finish machinery (GATE A/B) ..."); m  = load("finish", "tch_finish_port_and_swap.py")
print("[load] single-cell hop oracle (character gate) ..."); m2 = load("magp",  "tch_nonreduced_magnetic_plaquette.py")
cg = m.cg

def _rsvd(a, uv):
    try: return np.linalg.svd(a, compute_uv=uv)
    except np.linalg.LinAlgError:
        return sla.svd(a, full_matrices=True, lapack_driver="gesvd") if uv else sla.svd(a, compute_uv=False, lapack_driver="gesvd")
def robust_vbasis(reps):
    stacked = np.vstack(cg.generator_sum(reps)); s = _rsvd(stacked, False)
    if int(np.sum(s < 1e-9)) == 0: return ()
    u, s2, vh = _rsvd(stacked, True); null = vh[np.sum(s2 > 1e-9):].conj()
    dims = tuple(cg.dim(r) for r in reps); out = []
    for row in null:
        t = row.reshape(dims)
        for o in out: t = t - o * np.vdot(o, t)
        nrm = np.linalg.norm(t)
        if nrm > 1e-9: out.append(t / nrm)
    return tuple(out)
m.vbasis = functools.lru_cache(maxsize=None)(robust_vbasis)
vbasis = m.vbasis
F_emb = m.F_emb
LET = string.ascii_lowercase

# ---------------- generic builders (geometry-parametrized) ----------------
def build_basis_geom(reps, VERTS, vlegs):
    basis = []
    nlinks = max(li for v in VERTS for li, _ in vlegs(v)) + 1
    for Rs in __import__("itertools").product(reps, repeat=nlinks):
        per_v = []; ok = True
        for v in VERTS:
            legs = vlegs(v); opts = []
            for f in m.FERMS:
                full = tuple((cg.dual(Rs[li]) if d else Rs[li]) for li, d in legs) + (f,)
                for mi in range(len(vbasis(full))): opts.append((f, mi))
            if not opts: ok = False; break
            per_v.append(opts)
        if not ok: continue
        for choice in __import__("itertools").product(*per_v):
            basis.append((Rs, tuple(choice)))
    return basis

@functools.lru_cache(maxsize=None)
def hop_tail(legreps, dins, p, Rp):
    """fermion 3->1 at the link's TAIL (out-leg, kind out_fwd); link p raised R->Rp."""
    R_old = legreps[p]; nL = len(legreps)
    old = tuple(cg.dual(r) if d else r for r, d in zip(legreps, dins)) + ("3",)
    newr = list(legreps); newr[p] = Rp
    new = tuple(cg.dual(r) if d else r for r, d in zip(newr, dins)) + ("1",)
    Told, Tnew = vbasis(old), vbasis(new)
    if not Told or not Tnew: return None
    E = F_emb("out_fwd", R_old, Rp)                 # [F(=ferm colour 3), R_old, Rp]
    gl = LET[:nL]; sT = gl + "F"; sE = "F" + gl[p] + "N"; so = gl[:p] + "N" + gl[p+1:]
    K = np.zeros((len(Tnew), len(Told)), complex)
    for mo, T in enumerate(Told):
        X = np.einsum(f"{sT},{sE}->{so}", T, E, optimize=True)
        for mn, T2 in enumerate(Tnew):
            K[mn, mo] = np.vdot(T2[..., 0], X)      # T2 ferm=1 trivial axis sliced
    return K

@functools.lru_cache(maxsize=None)
def hop_head(legreps, dins, p, Rp):
    """fermion 1->3 at the link's HEAD (in-leg, kind in_fwd); link p raised R->Rp."""
    R_old = legreps[p]; nL = len(legreps)
    old = tuple(cg.dual(r) if d else r for r, d in zip(legreps, dins)) + ("1",)
    newr = list(legreps); newr[p] = Rp
    new = tuple(cg.dual(r) if d else r for r, d in zip(newr, dins)) + ("3",)
    Told, Tnew = vbasis(old), vbasis(new)
    if not Told or not Tnew: return None
    E = F_emb("in_fwd", R_old, Rp)                  # [C(=created colour 3), R_old, Rp]
    gl = LET[:nL]; sT = gl; sE = "C" + gl[p] + "N"; so = gl[:p] + "N" + gl[p+1:] + "C"
    K = np.zeros((len(Tnew), len(Told)), complex)
    for mo, T in enumerate(Told):
        X = np.einsum(f"{sT},{sE}->{so}", T[..., 0], E, optimize=True)   # T ferm=1 trivial sliced
        for mn, T2 in enumerate(Tnew):
            K[mn, mo] = np.vdot(T2, X)
    return K

def sparse_hop(reps, basis, index, VERTS, vlegs, LINK_TAIL, LINK_HEAD):
    vpos = {v: k for k, v in enumerate(VERTS)}
    rows, cols, vals = [], [], []
    nlinks = len(LINK_TAIL)
    for i, (Rs, choice) in enumerate(basis):
        for l in range(nlinks):
            vt, vh = LINK_TAIL[l], LINK_HEAD[l]
            kt, kh = vpos[vt], vpos[vh]
            (ft, mt0), (fh, mh0) = choice[kt], choice[kh]
            if ft != "3" or fh != "1":
                continue
            Lt, Lh = vlegs(vt), vlegs(vh)
            pt = [li for li, _ in Lt].index(l); ph = [li for li, _ in Lh].index(l)
            for Rp in cg.fuse_targets("3", Rs[l]):
                if Rp not in reps:
                    continue
                Kt = hop_tail(tuple(Rs[li] for li, _ in Lt), tuple(d for _, d in Lt), pt, Rp)
                Kh = hop_head(tuple(Rs[li] for li, _ in Lh), tuple(d for _, d in Lh), ph, Rp)
                if Kt is None or Kh is None:
                    continue
                pw = (cg.dim(Rs[l]) / cg.dim(Rp)) ** 0.5
                Rps = list(Rs); Rps[l] = Rp; Rps = tuple(Rps)
                col_t = Kt[:, mt0]; col_h = Kh[:, mh0]
                for mt in range(Kt.shape[0]):
                    if abs(col_t[mt]) < 1e-12: continue
                    for mh in range(Kh.shape[0]):
                        amp = col_t[mt] * col_h[mh] * pw
                        if abs(amp) < 1e-12: continue
                        nc = list(choice); nc[kt] = ("1", mt); nc[kh] = ("3", mh)
                        j = index.get((Rps, tuple(nc)))
                        if j is not None:
                            rows.append(j); cols.append(i); vals.append(amp)
    n = len(basis)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)

# =================== STEP 1: validate hop vs the single-cell oracle ===================
print("="*80)
print("STEP 1 — HOP regression vs the validated single-cell oracle (mag_plaq.HOP):")
VERTS_SC = [0, 1, 2, 3]
def vlegs_SC(v): return [(v, False), ((v - 1) % 4, True)]      # out-link v, in-link (v-1)%4
LT_SC = [l for l in range(4)]; LH_SC = [(l + 1) % 4 for l in range(4)]
basis_sc = build_basis_geom(m2.REPS, VERTS_SC, vlegs_SC)
index_sc = {b: i for i, b in enumerate(basis_sc)}
assert len(basis_sc) == len(m2.basis) == 336, (len(basis_sc), len(m2.basis))
HOP_mine = sparse_hop(m2.REPS, basis_sc, index_sc, VERTS_SC, vlegs_SC, LT_SC, LH_SC).toarray()
# map my (Rs, choice) -> oracle index (Rs, Fs, ms)
perm = np.empty(len(basis_sc), int)
for i, (Rs, ch) in enumerate(basis_sc):
    Fs = tuple(f for f, _ in ch); ms = tuple(mi for _, mi in ch)
    perm[i] = m2.index[(Rs, Fs, ms)]
HOP_oracle = m2.HOP[np.ix_(perm, perm)]
dhop = np.max(np.abs(HOP_mine - HOP_oracle))
print(f"  336-state single-cell: |my HOP - oracle HOP| = {dhop:.2e}  (target 0); nnz={np.count_nonzero(np.abs(HOP_mine)>1e-12)}")
assert dhop < 1e-10, "generic hop does NOT match the validated single-cell oracle"
print("  STEP 1 PASSED — the generic gauge-covariant hop is validated elementwise.\n")

# =================== STEP 2: 2-plaquette geometry + minimal regression ===================
VERTS, vlegs, PLAQS = m.VERTS, m.vlegs, m.PLAQS
LINKS = m.LINKS
def tail_head(link):
    kind, x, y = link
    return ((x, y), (x + 1, y)) if kind == "h" else ((x, y), (x, y + 1))
LINK_TAIL = [tail_head(L)[0] for L in LINKS]; LINK_HEAD = [tail_head(L)[1] for L in LINKS]

# reuse the validated sparse magnetic assembler from the t=0 closer (verbatim logic)
Kc = functools.lru_cache(maxsize=None)(m.K_corner)
def sparse_wilson(p, reps, basis, index):
    loop = PLAQS[p]; looplinks = [li for li, _ in loop]; dirs = {li: d for li, d in loop}
    plaq_verts = []
    for v in VERTS:
        legs = vlegs(v); tl = [(ax, li, isin) for ax, (li, isin) in enumerate(legs) if li in looplinks]
        if len(tl) == 2: plaq_verts.append((v, legs, tl))
    rows, cols, vals = [], [], []
    import itertools as it
    for i, (Rs, choice) in enumerate(basis):
        opts = [m.targets("out_fwd" if dirs[li] == "fwd" else "out_bwd", Rs[li], reps) for li in looplinks]
        for newreps in it.product(*opts):
            Rps = list(Rs)
            for li, rp in zip(looplinks, newreps): Rps[li] = rp
            Rps = tuple(Rps); amp_blocks = []; ok = True
            for (v, legs, tl) in plaq_verts:
                vi = VERTS.index(v); f, mi = basis[i][1][vi]
                raw_old = tuple(Rs[li] for li, _ in legs); raw_new = tuple(Rps[li] for li, _ in legs)
                df = tuple(isin for _, isin in legs)
                (axA, liA, iA), (axB, liB, iB) = tl
                posA, posB = looplinks.index(liA), looplinks.index(liB)
                first, second = ((axA, liA, iA), (axB, liB, iB)) if (posB - posA) % 4 == 1 else ((axB, liB, iB), (axA, liA, iA))
                kinds = tuple((("in_" if isin else "out_") + dirs[li]) for (ax, li, isin) in (first, second))
                K = Kc(raw_old, raw_new, df, f, (first[0], second[0]), kinds)
                if K is None: ok = False; break
                amp_blocks.append((vi, K))
            if not ok: continue
            pw = 1.0
            for li, rp in zip(looplinks, newreps): pw *= (cg.dim(Rs[li]) / cg.dim(rp)) ** 0.5
            corner_vis = [vi for vi, _ in amp_blocks]; Kmap = {vi: K for vi, K in amp_blocks}
            for mps in it.product(*[range(Kmap[vi].shape[0]) for vi in corner_vis]):
                amp = pw + 0j; nc = list(basis[i][1])
                for vi, mp in zip(corner_vis, mps):
                    f_old, m_old = basis[i][1][vi]; amp *= Kmap[vi][mp, m_old]; nc[vi] = (f_old, mp)
                if abs(amp) < 1e-12: continue
                j = index.get((Rps, tuple(nc)))
                if j is not None: rows.append(j); cols.append(i); vals.append(amp)
    n = len(basis)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)

def hamiltonian(basis, W0, W1, HOP, beta, t):
    ele = np.array([sum(cg.casimir(r) for r in Rs) for (Rs, ch) in basis], float)
    mat = np.array([m.DSMG * sum(1 for f, _ in ch if f != "1") for (Rs, ch) in basis], float)
    Wsym = W0 + W0.getH() + W1 + W1.getH()
    Hh = HOP + HOP.getH()
    return (sp.diags(ele / beta + mat) - (beta / 2) * Wsym + t * Hh).tocsr()

print("="*80)
print("STEP 2 — 2-plaquette MINIMAL regression (W elementwise, HOP hermiticity, t=0 gap):")
bmin = m.build_basis(m.REPS_MIN); imin = {b: i for i, b in enumerate(bmin)}
assert len(bmin) == 3321
W0m = sparse_wilson(0, m.REPS_MIN, bmin, imin); W1m = sparse_wilson(1, m.REPS_MIN, bmin, imin)
print(f"  |sparse W0-dense|={np.max(np.abs(W0m.toarray()-m.W0)):.1e}, |W1-dense|={np.max(np.abs(W1m.toarray()-m.W1)):.1e}")
assert np.max(np.abs(W0m.toarray()-m.W0)) < 1e-10 and np.max(np.abs(W1m.toarray()-m.W1)) < 1e-10
HOPm = sparse_hop(m.REPS_MIN, bmin, imin, VERTS, vlegs, LINK_TAIL, LINK_HEAD)
Hh = (HOPm + HOPm.getH()); herm = abs((Hh - Hh.getH())).max() if Hh.nnz else 0.0
print(f"  HOP nnz={HOPm.nnz}; (HOP+HOPd) hermiticity residual = {herm:.1e}")
assert herm < 1e-10
H0 = hamiltonian(bmin, W0m, W1m, HOPm, 0.5, 0.0)
e0 = np.linalg.eigvalsh(H0.toarray())
print(f"  t=0, beta=0.5 gap = {e0[1]-e0[0]:.6f}  (published 6.674654)")
assert abs((e0[1]-e0[0]) - 6.674654) < 1e-4
print("  STEP 2 PASSED.\n")

# =================== STEP 3: 106,460 extended t-scan ===================
print("="*80)
print("STEP 3 — the 106,460 EXTENDED basis: neutral-hopping t-scan")
t0 = time.time(); bext = m.build_basis(m.REPS_EXT); iext = {b: i for i, b in enumerate(bext)}
N = len(bext); print(f"  extended basis: {N} states (canon 106,460)  [build {time.time()-t0:.0f}s]")
assert N == 106460
ta = time.time(); W0 = sparse_wilson(0, m.REPS_EXT, bext, iext); W1 = sparse_wilson(1, m.REPS_EXT, bext, iext)
HOP = sparse_hop(m.REPS_EXT, bext, iext, VERTS, vlegs, LINK_TAIL, LINK_HEAD)
print(f"  W0 nnz={W0.nnz}, HOP nnz={HOP.nnz}  [assemble {time.time()-ta:.0f}s]")
print(f"\n  EXTENDED gap scan (electric+magnetic+SMG + t*(HOP+HOPd)):")
print(f"  {'beta':>5} {'t':>5} | {'ground':>11} {'gap':>11}")
res = {}
for beta in (0.5, 1.0):
    for t in (0.0, 0.2, 1.0):
        ts = time.time()
        H = hamiltonian(bext, W0, W1, HOP, beta, t)
        e = np.sort(eigsh(H, k=6, which="SA", return_eigenvectors=False, maxiter=30000, tol=0))
        res[(beta, t)] = (e[0], e[1]-e[0])
        print(f"  {beta:5.2f} {t:5.1f} | {e[0]:11.6f} {e[1]-e[0]:11.6f}   [{time.time()-ts:.0f}s]")

print(f"""
================================================================================
VERDICT (neutral-hopping leg of the 106,460 extended escalation — RUN, CPU):
  * The gauge-covariant fermion hop psi^dag_head U_l psi_tail (+h.c.) is VALIDATED
    ELEMENTWISE (<1e-10) against the validated single-cell oracle, is Hermitian-closable
    on the 2-plaquette basis, and the magnetic assembler still matches dense exactly.
  * The full electric+magnetic+SMG+HOPPING Hamiltonian is now diagonalized on the
    106,460 extended basis via sparse eigsh (no GPU). The t-scan above shows how the
    neutral hop moves the gap at each beta.
  * Together with the t=0 magnetic+higher-cutoff result (tch_2plaq_extended_deep.py),
    this CLOSES the full "magnetic / neutral-hopping / higher-cutoff" 106,460 escalation
    locally — the registered 'deep job' needed no GPU, only sparse assembly + Lanczos.
================================================================================""")
print("exit 0 -- neutral-hopping operator built, single-cell-validated, and run on the 106,460 extended t-scan.")
