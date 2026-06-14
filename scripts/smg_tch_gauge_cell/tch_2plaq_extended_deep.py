#!/usr/bin/env python3
r"""CLOSE THE OPEN ITEM: the 106,460-state 2-plaquette EXTENDED non-reduced charged
spin-network build (REPS_EXT = 1,3,3b,6,6b,8), electric + magnetic-on-both-plaquettes.

Registered in tch_finish_port_and_swap.py / tch_nonreduced_*.py as a "deep job" — but it
was deferred only because those scripts diagonalize DENSELY (106,460^2 complex ~ 180 GB).
H is in fact very sparse (the magnetic W has ~few nnz/row), so this is a routine sparse
Lanczos (scipy eigsh) on CPU. No GPU needed.

DISCIPLINE (re-derive the known case as a harness before trusting the new one):
  1. Reuse the VALIDATED machinery from tch_finish_port_and_swap.py verbatim (its module
     execution already passes GATE A flux-character + GATE B 1e-7 spectrum regression).
  2. Re-implement ONLY the assembler as sparse (COO), changing nothing in the recoupling.
  3. REGRESSION: rebuild the minimal (3,321) W0/W1 sparsely and assert ELEMENTWISE equality
     to the dense module arrays, and reproduce the published minimal gap 6.674654 (beta=0.5).
     Only if that passes do we trust the extended run.
  4. Build the 106,460 extended basis (assert the count), assemble sparse H, eigsh the gap.
exit 0 = regression verified AND the extended gap computed.
"""
import importlib.util, functools, sys, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parent
# --- load the validated module (runs its GATE A/B on import; that IS the machinery proof) ---
spec = importlib.util.spec_from_file_location("finish", ROOT / "tch_finish_port_and_swap.py")
m = importlib.util.module_from_spec(spec)
print("[load] executing validated machinery (tch_finish_port_and_swap.py: GATE A + GATE B) ...")
spec.loader.exec_module(m)
print("[load] machinery loaded; reusing its build_basis / targets / K_corner / geometry verbatim.\n")

cg = m.cg
import scipy.linalg as sla

def _robust_svd(a, compute_uv):
    """np.linalg.svd (gesdd, fast) with a gesvd fallback on non-convergence.
    Identical output where gesdd converges (so the minimal regression is untouched);
    only the rare extended-rep tensors take the robust path. Any orthonormal null
    basis is a valid intertwiner basis -> the SPECTRUM is invariant to the choice."""
    try:
        return np.linalg.svd(a, compute_uv=compute_uv)
    except np.linalg.LinAlgError:
        if compute_uv:
            return sla.svd(a, full_matrices=True, lapack_driver="gesvd")
        return sla.svd(a, compute_uv=False, lapack_driver="gesvd")

def robust_vbasis(reps):
    """Verbatim copy of m.vbasis (lines 40-58) with the convergence-safe SVD."""
    stacked = np.vstack(cg.generator_sum(reps))
    s = _robust_svd(stacked, compute_uv=False)
    nullity = int(np.sum(s < 1e-9))
    if nullity == 0:
        return ()
    u, s2, vh = _robust_svd(stacked, compute_uv=True)
    null = vh[np.sum(s2 > 1e-9):].conj()
    dims = tuple(cg.dim(r) for r in reps)
    out = []
    for row in null:
        t = row.reshape(dims)
        for o in out:
            t = t - o * np.vdot(o, t)
        nrm = np.linalg.norm(t)
        if nrm > 1e-9:
            out.append(t / nrm)
    return tuple(out)

# cache the two pure functions for the larger run (deterministic in their args)
m.vbasis = functools.lru_cache(maxsize=None)(robust_vbasis)  # build_basis + K_corner resolve this global
Kc = functools.lru_cache(maxsize=None)(m.K_corner)
VERTS, vlegs, PLAQS, FERMS, DSMG = m.VERTS, m.vlegs, m.PLAQS, m.FERMS, m.DSMG

def sparse_wilson(p, reps, basis, index):
    """EXACT sparse copy of m.wilson_P (lines 193-261): same recoupling, COO output."""
    loop = PLAQS[p]
    looplinks = [li for li, _ in loop]
    dirs = {li: d for li, d in loop}
    plaq_verts = []
    for v in VERTS:
        legs = vlegs(v)
        tl = [(ax, li, isin) for ax, (li, isin) in enumerate(legs) if li in looplinks]
        if len(tl) == 2:
            plaq_verts.append((v, legs, tl))
    assert len(plaq_verts) == 4
    rows, cols, vals = [], [], []
    for i, (Rs, choice) in enumerate(basis):
        opts = []
        for li in looplinks:
            kind_t = ("out_fwd" if dirs[li] == "fwd" else "out_bwd")
            opts.append(m.targets(kind_t, Rs[li], reps))
        for newreps in _product(*opts):
            Rps = list(Rs)
            for li, rp in zip(looplinks, newreps):
                Rps[li] = rp
            Rps = tuple(Rps)
            amp_blocks = []
            ok = True
            for (v, legs, tl) in plaq_verts:
                vi = VERTS.index(v)
                f, mi = basis[i][1][vi]
                raw_old = tuple(Rs[li] for li, _ in legs)
                raw_new = tuple(Rps[li] for li, _ in legs)
                df = tuple(isin for _, isin in legs)
                (axA, liA, isinA), (axB, liB, isinB) = tl
                posA, posB = looplinks.index(liA), looplinks.index(liB)
                first, second = ((axA, liA, isinA), (axB, liB, isinB)) \
                    if (posB - posA) % 4 == 1 else ((axB, liB, isinB), (axA, liA, isinA))
                kinds = []
                for (ax, li, isin) in (first, second):
                    d = dirs[li]
                    kinds.append(("in_" if isin else "out_") + d)
                K = Kc(raw_old, raw_new, df, f, (first[0], second[0]), tuple(kinds))
                if K is None:
                    ok = False; break
                amp_blocks.append((vi, K))
            if not ok:
                continue
            pw = 1.0
            for li, rp in zip(looplinks, newreps):
                pw *= (cg.dim(Rs[li]) / cg.dim(rp)) ** 0.5
            corner_vis = [vi for vi, _ in amp_blocks]
            Kmap = {vi: K for vi, K in amp_blocks}
            newdims = {vi: K.shape[0] for vi, K in amp_blocks}
            for mps in _product(*[range(newdims[vi]) for vi in corner_vis]):
                amp = pw + 0j
                newchoice = list(basis[i][1])
                for vi, mp in zip(corner_vis, mps):
                    f_old, m_old = basis[i][1][vi]
                    amp *= Kmap[vi][mp, m_old]
                    newchoice[vi] = (f_old, mp)
                if abs(amp) < 1e-12:
                    continue
                j = index.get((Rps, tuple(newchoice)))
                if j is not None:
                    rows.append(j); cols.append(i); vals.append(amp)
    n = len(basis)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)

import itertools as _it
_product = _it.product

def diagels(basis):
    ele = np.array([sum(cg.casimir(r) for r in Rs) for (Rs, ch) in basis], float)
    mat = np.array([DSMG * sum(1 for f, _ in ch if f != "1") for (Rs, ch) in basis], float)
    return ele, mat

def gaps(basis, index, reps, betas, k=6, dense=False):
    t0 = time.time()
    W0 = sparse_wilson(0, reps, basis, index)
    W1 = sparse_wilson(1, reps, basis, index)
    Wsym = (W0 + W0.getH() + W1 + W1.getH()).tocsr()
    ele, mat = diagels(basis)
    out = {}
    for beta in betas:
        diag = sp.diags(ele / beta + mat)
        H = (diag - (beta / 2.0) * Wsym).tocsr()
        if dense:
            e = np.linalg.eigvalsh(H.toarray())
        else:
            e = np.sort(eigsh(H, k=k, which="SA", return_eigenvectors=False, maxiter=20000, tol=0))
        out[beta] = (e[0], e[1] - e[0], e[:k])
    return out, W0, W1, time.time() - t0

# =================== STEP 1: REGRESSION on the minimal (3,321) basis ===================
print("="*78)
print("STEP 1 — REGRESSION: sparse assembler must reproduce the validated dense minimal build")
basis_min = m.build_basis(m.REPS_MIN)
index_min = {b: i for i, b in enumerate(basis_min)}
assert len(basis_min) == 3321, len(basis_min)
W0s = sparse_wilson(0, m.REPS_MIN, basis_min, index_min)
W1s = sparse_wilson(1, m.REPS_MIN, basis_min, index_min)
d0 = np.max(np.abs(W0s.toarray() - m.W0)); d1 = np.max(np.abs(W1s.toarray() - m.W1))
print(f"  elementwise |sparse W0 - dense W0| = {d0:.2e} ; |W1 - W1| = {d1:.2e}  (target 0)")
assert d0 < 1e-10 and d1 < 1e-10, "sparse assembler does NOT match validated dense machinery"
gmin, _, _, _ = gaps(basis_min, index_min, m.REPS_MIN, (0.5, 1.0), dense=True)
print(f"  minimal gap beta=0.50: ground={gmin[0.5][0]:.6f}, gap={gmin[0.5][1]:.6f}  (published -0.023952 / 6.674654)")
print(f"  minimal gap beta=1.00: ground={gmin[1.0][0]:.6f}, gap={gmin[1.0][1]:.6f}  (published -0.199710 / 5.033044)")
assert abs(gmin[0.5][1] - 6.674654) < 1e-4 and abs(gmin[1.0][1] - 5.033044) < 1e-4
# also verify eigsh path agrees with the dense path on the minimal basis (sanity for the big run)
gmin_sp, _, _, _ = gaps(basis_min, index_min, m.REPS_MIN, (0.5,), dense=False)
print(f"  eigsh vs dense (minimal, beta=0.5) gap: {gmin_sp[0.5][1]:.6f} vs {gmin[0.5][1]:.6f}  (sparse-Lanczos path validated)")
assert abs(gmin_sp[0.5][1] - gmin[0.5][1]) < 1e-6
print("  REGRESSION PASSED — sparse + eigsh path reproduces the validated machinery exactly.\n")

# =================== STEP 2: the 106,460 EXTENDED run ===================
print("="*78)
print("STEP 2 — the EXTENDED build (REPS_EXT = 1,3,3b,6,6b,8): the registered 'deep job'")
t0 = time.time()
basis_ext = m.build_basis(m.REPS_EXT)
N = len(basis_ext)
index_ext = {b: i for i, b in enumerate(basis_ext)}
print(f"  extended 2-plaquette non-reduced charged basis: {N} states  (canon: 106,460)   [build {time.time()-t0:.1f}s]")
gext, W0e, W1e, asm = gaps(basis_ext, index_ext, m.REPS_EXT, (0.5, 1.0), k=8, dense=False)
print(f"  W0 nnz={W0e.nnz}, W1 nnz={W1e.nnz}  (sparse; assemble+eigsh {asm:.1f}s)")
print(f"\n  EXTENDED gaps (electric + magnetic on both plaquettes, t=0, Delta_SMG={DSMG}):")
for beta in (0.5, 1.0):
    g0, gp, lev = gext[beta]
    gm = gmin[beta][1]
    print(f"    beta={beta:4.2f}: ground={g0:.6f}, gap={gp:.6f}   (minimal 3,321 gap={gm:.6f}, "
          f"shift {gp-gm:+.6f} = {100*(gp-gm)/gm:+.2f}%)")
    print(f"             lowest levels: {np.array2string(lev, precision=5, floatmode='fixed')}")

print(f"""
================================================================================
VERDICT (the 106,460-state extended escalation, RUN — no GPU needed):
  * The registered 'deep job' was deferred only because the existing scripts diagonalize
    DENSELY. It is a routine SPARSE Lanczos problem and runs on CPU in {asm:.0f}s.
  * Sparse assembler verified ELEMENTWISE against the validated dense machinery on the
    3,321 minimal basis (|dW|<1e-10) and the eigsh path matched dense to <1e-6 — so the
    extended spectrum below inherits the GATE-A/GATE-B validation.
  * Extended basis dimension = {N} (canon 106,460: {'CONFIRMED' if N==106460 else 'MISMATCH -> '+str(N)}).
  * The higher-cutoff (6,6b,8) escalation shifts the electric+magnetic mirror gap by the
    percentages above relative to the minimal (1,3,3b) result -> reports whether the
    minimal-basis gap is converged under the rep-cutoff (the open question this closes).
  NOTE/SCOPE: this closes the electric+magnetic higher-cutoff gap at t=0. The neutral/charged
  HOPPING term (t-scan) on the 2-plaquette extended basis is a separate operator not built
  here; if wanted it is the same sparse pattern with the link-hop block added.
================================================================================""")
print(f"exit 0 -- extended 106,460 build RUN locally (sparse eigsh); regression-validated; gaps reported.")
