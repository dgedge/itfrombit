#!/usr/bin/env python3
r"""Item 48 -- KPM local DOS of the rho 'P4 string' embedded in the 3D substrate continuum.

PURPOSE.  Yesterday's exact-diagonalisation (item48_feshbach_substrate.py, L<=12, N<=1728)
found the phi-mode dissolves into a band-wide SMEAR at the substrate gap Delta/t~=1.78. That
verdict could in principle be a finite-volume artefact: a genuinely NARROW Fano resonance can
fall between the coarse levels of a small box and be missed. The Kernel Polynomial Method
(KPM) resolves the *intrinsic* line shape directly in the thermodynamic limit -- no matrix
inversion, only SpMV -- so we can settle whether the rho is a sharp pole near ~1.52 (framework
survives) or a >GeV smear (light-vector sector needs rebuild).

WHAT THIS IS NOT.  This does not invent a number. It computes the path-projected LDOS of an
explicitly-defined Hamiltonian and reports peak / FWHM / background. Two modelling inputs remain
UNPINNED (flagged, not hidden):
  * the bulk graph identity -- ANCHOR pins it two ways: the macroscopic gauge web is the
    *line graph of the SC lattice* (sec 7.3, "...the line graph of the SC lattice"; band [-2,10],
    z=10), while yesterday's exact-diag used the SC lattice itself as a PROXY (band [-6,6], z=6).
    We run BOTH (--graph) so the verdict's robustness to this choice is visible, not assumed.
  * the leak amplitude t_leak = t_path = 1 (unweighted induced P4 in the bulk). A weaker leak
    makes survival easier; this embedding cannot tune it independently (same caveat as yesterday).

CORRECTION TO THE PROSE SPEC.  The spec said "E_max/min ~ +-6t, a~7, b~0". That is the band of
the SC lattice, NOT of the line graph it tells you to build (N=3L^3). The line graph L(SC) has
band [-2,10]; rescaling it with a=7,b=0 puts eigenvalues at +10/7>1 and the Chebyshev recurrence
DIVERGES. Bounds are therefore estimated per-graph by Lanczos/power iteration with a safety pad.

SELF-ASSERTING.  --validate runs a small lattice and asserts (i) KPM moments == exact
<j|T_n(Htilde)|j> from full eigh to <1e-9, (ii) reconstructed LDOS integrates to 1, (iii) the
spectrum sits strictly inside the rescaled (-1,1). Exit 0 => the KPM engine is verified; the
physics verdict is then REPORTED for the requested (graph, L, Delta, moments).

Backends: --backend cpu (numpy/scipy) or gpu (cupy/cupyx). Same code path; the 4 path-site LDOS
vectors are carried as 4 columns of one dense block so a single SpMM amortises each matrix read.
"""
import argparse
import sys
import time

import numpy as np
import scipy.sparse as ssp


# ----------------------------------------------------------------------------- geometry
def build_sc(L):
    """Simple-cubic L^3, periodic. Returns (CSR float64 adjacency, path node indices).
    Band [-6,6], z=6. This is yesterday's PROXY bulk."""
    N = L ** 3

    def lin(x, y, z):
        return (x % L) * L * L + (y % L) * L + (z % L)

    x, y, z = np.meshgrid(np.arange(L), np.arange(L), np.arange(L), indexing="ij")
    src = lin(x, y, z).ravel()
    rows = []
    cols = []
    for dx, dy, dz in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:
        dst = lin(x + dx, y + dy, z + dz).ravel()
        rows.append(src); cols.append(dst)      # +direction
        rows.append(dst); cols.append(src)      # -direction (Hermitian)
    rows = np.concatenate(rows).astype(np.int32)
    cols = np.concatenate(cols).astype(np.int32)
    data = np.ones(rows.size, dtype=np.float64)
    A = ssp.csr_matrix((data, (rows, cols)), shape=(N, N))
    A.sum_duplicates()
    cen = L // 2
    path = [lin(cen + i, cen, cen) for i in range(4)]   # 4 collinear sites = induced P4
    return A, path


def build_linegraph(L):
    """Line graph L(SC) of the simple-cubic L^3 (periodic): nodes = SC edges, N=3L^3.
    Two edges adjacent iff they share an SC vertex => each degree-6 vertex contributes a K6
    clique on its 6 incident edges. Band [-2,10], z=10. This is ANCHOR sec 7.3's macroscopic
    gauge web (line graph of the SC lattice). Returns (CSR float64 adjacency, path indices)."""
    L3 = L ** 3

    def node(x, y, z, d):              # edge (vertex (x,y,z) -> +e_d), d in {0,1,2}
        return 3 * ((x % L) * L * L + (y % L) * L + (z % L)) + d

    x, y, z = np.meshgrid(np.arange(L), np.arange(L), np.arange(L), indexing="ij")
    x = x.ravel(); y = y.ravel(); z = z.ravel()
    # 6 edges incident to each vertex v: 3 outgoing (v,d) + 3 incoming (v-e_d, d)
    inc = np.stack([
        node(x, y, z, 0), node(x, y, z, 1), node(x, y, z, 2),
        node(x - 1, y, z, 0), node(x, y - 1, z, 1), node(x, y, z - 1, 2),
    ], axis=1).astype(np.int32)        # (L3, 6)
    rows = []
    cols = []
    for i in range(6):
        for j in range(i + 1, 6):      # 15 unordered pairs -> K6 clique edges
            rows.append(inc[:, i]); cols.append(inc[:, j])
            rows.append(inc[:, j]); cols.append(inc[:, i])
    rows = np.concatenate(rows)
    cols = np.concatenate(cols)
    data = np.ones(rows.size, dtype=np.float64)
    A = ssp.csr_matrix((data, (rows, cols)), shape=(3 * L3, 3 * L3))
    A.sum_duplicates()                 # collapse the (rare) shared-vertex multiplicities to 1
    A.data[:] = 1.0                    # adjacency, not multigraph
    cen = L // 2
    # 4 consecutive collinear x-edges share endpoints pairwise => induced P4
    path = [node(cen + i, cen, cen, 0) for i in range(4)]
    return A, path


# cube-graph Q3 edges on labels 0..7 (Hamming-distance-1); label 0 = shared apex/centre
_Q3_EDGES = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6),
             (4, 5), (4, 6), (3, 7), (5, 7), (6, 7)]


def build_matter(L, t_leak, sublattice=False, gauge="sc"):
    r"""'matter' substrate -- the honest isolated-matter / silver-ratio-leak model.

    The rho is NOT a bulk node here. It is an internally-coherent matter string (P4, hopping 1,
    leading eigenvalue phi) that couples to the gauge continuum ONLY through weak 'triangular-face'
    channels at the silver-ratio amplitude t_leak = sqrt(2)-1. This is the configuration the earlier
    runs did NOT test (they used t_leak = t_path = 1, i.e. the string WAS bulk).

      * Gauge continuum  : gauge='sc' -> SC web of truncated-cube centres (L^3, band [-6,6], the PROXY);
                           gauge='linegraph' -> the ACTUAL L(SC) gauge web of sec 7.3 (3L^3, band [-2,10],
                           with the flat band at -2). The continuum DOS at Delta+phi differs between the
                           two, so the resonance WIDTH does -- gauge='linegraph' tests the proxy artifact.
      * rho-host         : 4 P4 nodes (internal hopping 1; on-site Delta added later via add_defect),
                           each coupled to one central gauge node by t_leak (4 triangular-face channels).
      * sublattice (opt) : ANCHOR's matter cells -- per cubic cell, 3 orthogonal bipyramids as Q3 cube
                           graphs (8 nodes), the three SHARING ONE per-cell centre node (line 61) and
                           NO inter-cell matter bonds (line 51: one bipyramid per tiling vertex -> cells
                           share nothing). Each bipyramid couples to its cell's gauge node by ONE
                           t_leak channel. These cells are spectators: being isolated, they can dress
                           the rho-host only at order t_leak^4 -- the flag lets us verify that numerically.

    Returns (CSR adjacency, rho-host path indices). Band has no clean analytic form -> BAND['matter']=None.
    """
    # gauge continuum (the bulk the rho leaks into) -- reuse the validated builders
    Agauge, gpath = build_sc(L) if gauge == "sc" else build_linegraph(L)
    Ng = Agauge.shape[0]
    rho = [Ng + i for i in range(4)]
    Ag = Agauge.tocoo()
    rows = [Ag.row.astype(np.int64)]; cols = [Ag.col.astype(np.int64)]; data = [Ag.data.astype(np.float64)]

    def bond(i, j, t):
        rows.append(np.array([i, j], np.int64)); cols.append(np.array([j, i], np.int64))
        data.append(np.array([t, t], np.float64))

    # rho-host P4 (separate coherent chain) coupled to the 4 central gauge nodes by silver-ratio channels
    for i in range(3):
        bond(rho[i], rho[i + 1], 1.0)                       # P4 internal hopping -> leading eigenvalue phi
    for i in range(4):
        bond(rho[i], int(gpath[i]), t_leak)                 # triangular-face channel, silver ratio

    nxt = Ng + 4
    if sublattice:
        if gauge != "sc":
            raise ValueError("--sublattice is implemented only for the SC gauge continuum")
        def g(x, y, z):
            return ((x % L) * L + (y % L)) * L + (z % L)
        cells = np.arange(Ng)
        cx, cy, cz = np.unravel_index(cells, (L, L, L))
        centre = nxt + cells                                # one centre node per cell
        bip0 = nxt + Ng                                     # first bipyramid node index
        for c in range(Ng):
            gc = g(cx[c], cy[c], cz[c])
            ctr = int(centre[c])
            bond(ctr, gc, t_leak)                            # centre couples to gauge
            for d in range(3):
                lab = {0: ctr}
                for L_ in range(1, 8):
                    lab[L_] = bip0 + ((c * 3 + d) * 7 + (L_ - 1))
                for a, b in _Q3_EDGES:
                    bond(lab[a], lab[b], 1.0)                # Q3 cube-graph internal hopping
                bond(lab[1], gc, t_leak)                     # ONE triangular-face channel per bipyramid
        N = nxt + Ng + 21 * Ng
    else:
        N = nxt

    rows = np.concatenate(rows); cols = np.concatenate(cols); data = np.concatenate(data)
    A = ssp.csr_matrix((data, (rows, cols)), shape=(N, N))
    A.sum_duplicates()
    return A, rho


def add_defect(A, path, delta):
    """H = A + delta * (on-site projector on the 4 P4 path nodes). Returns CSR."""
    N = A.shape[0]
    diag = np.zeros(N, dtype=np.float64)
    diag[path] = delta
    return (A + ssp.diags(diag)).tocsr()


def build_graph(graph, L, t_leak=None, sublattice=False, gauge="sc"):
    if graph == "sc":
        return build_sc(L)
    if graph == "linegraph":
        return build_linegraph(L)
    if graph == "matter":
        return build_matter(L, t_leak, sublattice, gauge)
    raise ValueError(graph)


# band used by extract() for the in-band / bound-state test (matter: the gauge continuum it leaks into)
EXTRACT_BAND = {"sc": (-6.0, 6.0), "linegraph": (-2.0, 10.0), "matter": (-6.0, 6.0)}


# ----------------------------------------------------------------------------- spectral bounds
# analytic clean-bulk band edges (defect-free), per graph. None => bound both edges by power iteration
# (matter mode: gauge continuum + discrete matter levels + leak broadening has no simple analytic band).
BAND = {"sc": (-6.0, 6.0), "linegraph": (-2.0, 10.0), "matter": None}


def _power_top(H, shift, iters, seed):
    """Rayleigh quotient of the dominant eigenvector of (H + shift*I) -> lambda_max(H) (if shift makes
    the top dominant) ; caller picks shift so the wanted edge is the largest-magnitude one."""
    N = H.shape[0]
    rng = np.random.default_rng(seed)
    v = rng.standard_normal(N); v /= np.linalg.norm(v)
    for _ in range(iters):
        w = H @ v + shift * v
        v = w / np.linalg.norm(w)
    return float(v @ (H @ v))


def spectral_bounds(H, band, delta, pad=0.02, iters=200, seed=0):
    """(Emin, Emax) WITHOUT eigsh (which never converges on 1e7-1e8-nnz matrices, and is hopeless on
    the line graph's hugely degenerate E=-2 flat band).
      * band given: emin = analytic clean lower edge (Weyl: H = A + Delta*P, P PSD, Delta>0 =>
        lambda_min(H) >= lambda_min(A) = band_lo, rigorous & tight); emax = shifted power iteration.
      * band None (matter mode): bound BOTH edges by shifted power iteration -- emax from (H + s+ I),
        emin from -(- H + s- I) i.e. the top of (sI - H) gives lambda_min(H).
    Returns (emin_pad, emax_pad, emin_raw, emax_raw)."""
    if band is not None:
        lo, _ = band
        emax = _power_top(H, -lo + delta + 2.0, iters, seed)   # shift -> (H+sI) all positive, top dominant
        emin = lo
    else:
        absmax = float(np.asarray(np.abs(H).sum(axis=1)).ravel().max()) + delta   # Gershgorin radius, safe shift
        emax = _power_top(H, absmax, iters, seed)               # top of (H + |.|I) -> lambda_max
        emin = -_power_top(-H, absmax, iters, seed + 1)         # top of (-H + |.|I) -> -lambda_min => lambda_min
    assert emax > emin, f"degenerate bracket [{emin},{emax}]"
    span = emax - emin
    return emin - pad * span, emax + pad * span, emin, emax


# ----------------------------------------------------------------------------- KPM core
def jackson(Nc):
    n = np.arange(Nc)
    q = np.pi / (Nc + 1)
    return ((Nc - n + 1) * np.cos(q * n) + np.sin(q * n) / np.tan(q)) / (Nc + 1)


def kpm_local_moments(H, path, a, b, Nc, backend="cpu", device=0, report=True, start=None):
    r"""Chebyshev moments mu_n = (1/K) sum_k <c_k| T_n((H-b)/a) |c_k> over K unit columns c_k.
    start=None -> the K=len(path) unit site vectors (per-site path LDOS). start=(N,K) host array of
    unit-norm columns -> e.g. random vectors for the stochastic per-site bulk DOS (1/N)Tr T_n.
    Carries the columns as a dense (N,K) block; one SpMM per step. Doubling => ~Nc/2 SpMMs."""
    N = H.shape[0]
    K = len(path) if start is None else start.shape[1]
    if backend == "gpu":
        import cupy as xp
        import cupyx.scipy.sparse as xsp
        xp.cuda.Device(device).use()
        Hs = xsp.csr_matrix(H.astype(np.float64))
        inv_a = np.float64(1.0 / a)
        b_ = np.float64(b)
        sync = lambda: xp.cuda.Stream.null.synchronize()
    else:
        import numpy as xp
        Hs = H
        inv_a = 1.0 / a
        b_ = b
        sync = lambda: None

    # Htilde @ X  without forming Htilde: (H@X - b X)/a
    def Ht(X):
        return (Hs @ X - b_ * X) * inv_a

    if start is None:
        R0 = xp.zeros((N, K), dtype=xp.float64)
        cols = xp.asarray(np.array(path, dtype=np.int64))
        R0[cols, xp.arange(K)] = 1.0
    else:
        R0 = xp.asarray(start.astype(np.float64))
    a0 = R0
    a1 = Ht(a0)

    mu = xp.zeros(Nc, dtype=xp.float64)
    # mu_n = <a0 | a_n>, taken column-wise then averaged over the K path sites
    def coldot(P, Q):
        return (P * Q).sum(axis=0)
    mu[0] = coldot(a0, a0).mean()
    if Nc > 1:
        mu[1] = coldot(a0, a1).mean()
    # doubling: mu_{2k} = 2<a_k|a_k> - mu_0 ; mu_{2k+1} = 2<a_{k+1}|a_k> - mu_1
    half = (Nc + 1) // 2
    am1, a_k = a0, a1
    t0 = time.time()
    for k in range(1, half):
        a_next = 2.0 * Ht(a_k) - am1
        am1, a_k = a_k, a_next          # now a_k = alpha_{k+1}? -> track indices explicitly below
        # after this step a_k holds alpha_{k+1}; am1 holds alpha_k
        idx2 = 2 * k
        idx2p = 2 * k + 1
        if idx2 < Nc:
            mu[idx2] = 2.0 * coldot(am1, am1).mean() - mu[0]
        if idx2p < Nc:
            mu[idx2p] = 2.0 * coldot(a_k, am1).mean() - mu[1]
        if report and backend == "gpu" and (k % max(1, half // 10) == 0):
            sync()
            print(f"    [kpm] {2 * k}/{Nc} moments  {time.time() - t0:6.1f}s", flush=True)
    sync()
    mu_host = mu.get() if backend == "gpu" else np.asarray(mu)
    return mu_host


def reconstruct(mu, a, b, npts=20000, chunk=4000):
    """rho(E) on a fine grid from Jackson-damped moments. Returns (E, rho) in physical E units.
    Chunked over the energy grid so the npts x Nc Chebyshev matrix never blows up memory."""
    Nc = mu.size
    n = np.arange(Nc)
    gm = jackson(Nc) * mu
    xs = np.linspace(-1, 1, npts + 2)[1:-1]              # avoid the +-1 singularity
    rho_x = np.empty_like(xs)
    for i in range(0, xs.size, chunk):
        xc = xs[i:i + chunk]
        T = np.cos(np.outer(np.arccos(xc), n))           # T_n(x) = cos(n arccos x)
        rho_x[i:i + chunk] = gm[0] + 2.0 * (T[:, 1:] @ gm[1:])
    rho_x /= (np.pi * np.sqrt(1 - xs ** 2))
    return a * xs + b, rho_x / a                          # Jacobian; integral over E stays 1


# ----------------------------------------------------------------------------- features
def _fwhm_around(E, rho, ipk):
    half = rho[ipk] / 2.0
    li = ipk
    while li > 0 and rho[li] > half:
        li -= 1
    ri = ipk
    while ri < len(rho) - 1 and rho[ri] > half:
        ri += 1
    return E[li], E[ri], E[ri] - E[li]


def _smooth(E, rho, sig):
    """Gaussian-smooth rho(E) with std sig (energy units). For robust peak/FWHM finding."""
    dE = E[1] - E[0]
    half = int(max(1, round(4 * sig / dE)))
    k = np.exp(-0.5 * (np.arange(-half, half + 1) * dE / sig) ** 2)
    k /= k.sum()
    return np.convolve(rho, k, mode="same")


def extract(E, rho, delta, phi, t_mev, band, raw_emin, raw_emax, gamma_mev=149.0):
    """Discriminate 'sharp resonance' from 'dissolved smear' for the path LDOS rho(E).

    Rigorous diagnostics (analytic band edges; no fragile DOS auto-detection, no edge artifact):
      * bound gaps -- raw largest/smallest eigenvalue beyond the clean band => a split-off bound state.
      * conc_phys  -- fraction of the string's weight inside a physical-width (Gamma=149 MeV) window
                      at the would-be rho level Delta+phi, vs a flat smear over the band (enhancement).
      * lineshape  -- peak & FWHM of rho smoothed at the physical resolution, near Delta+phi.
    centroid/sigma are reported but flagged trivial (mu_1=H_jj=Delta; mu_2-mu_1^2 = sum|H_jk|^2 = degree).
    """
    Elo, Ehi = band                                      # analytic clean-band edges, defect-free
    rho = np.clip(rho, 0, None)
    norm = float(np.trapezoid(rho, E))
    w = rho / norm
    centroid = float(np.trapezoid(w * E, E))
    sigma = float(np.trapezoid(w * (E - centroid) ** 2, E)) ** 0.5
    E_res = delta + phi
    gw = gamma_mev / t_mev
    out = dict(norm=norm, centroid=centroid, centroid_rel=centroid - delta, sigma=sigma,
               sigma_mev=sigma * t_mev, E_res=E_res, gw=gw, band_lo=Elo, band_hi=Ehi)

    # bound state = an eigenvalue clearly beyond the clean band (the defect pushes UP for Delta>0)
    out["gap_above"] = raw_emax - Ehi                    # >0 and O(0.1+) => split-off state above the band
    out["gap_below"] = Elo - raw_emin
    out["E_res_in_band"] = Elo < E_res < Ehi

    # concentration at the would-be rho level vs a flat (fully smeared) band of the same width
    win = (E > E_res - gw / 2) & (E < E_res + gw / 2)
    flat = gw / (Ehi - Elo)
    out["conc_phys"] = float(np.trapezoid(rho[win], E[win]) / norm)
    out["conc_over_flat"] = out["conc_phys"] / flat if flat > 0 else float("nan")

    # lineshape of the phi-level: smooth only lightly (~25 MeV floor, NOT the 149 MeV physical width --
    # else every feature reads as >=149 MeV and a surviving narrow rho is masked), and lock the search
    # to a tight window around Delta+phi so we track the phi-level, not the Delta+1/phi sub-level.
    smooth_mev = 25.0
    rs = _smooth(E, w, (smooth_mev / t_mev) / 2.355)
    wmask = (E > E_res - 0.45) & (E < E_res + 0.45)
    idxs = np.where(wmask)[0]
    ip = idxs[int(np.argmax(rs[idxs]))]
    _, _, fw = _fwhm_around(E, rs, ip)
    fw_deconv = (max(fw * t_mev, smooth_mev) ** 2 - smooth_mev ** 2) ** 0.5   # remove smoothing in quadrature
    out.update(Epk=float(E[ip]), Epk_rel=float(E[ip] - delta), pk_fwhm=fw,
               pk_fwhm_mev=fw * t_mev, pk_fwhm_deconv_mev=fw_deconv, smooth_mev=smooth_mev,
               rho_at_Eres=float(rs[np.argmin(np.abs(E - E_res))] * norm))   # un-normalised LDOS height at Delta+phi
    return out


def make_plot(E, rho, rho_bulk, f, delta, phi, fname):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"# (no plot: {e})")
        return
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(E, rho, lw=1.1, label="rho_path (P4 string sites)")
    if rho_bulk is not None:
        rb = np.clip(rho_bulk, 0, None); rb /= np.trapezoid(rb, E)
        ax.plot(E, rb, lw=0.9, alpha=0.6, label="rho_bulk (defect-free per-site)")
    for x, c, lab in [(delta + phi, "r", "Delta+phi (would-be rho)"),
                      (delta + 1.52, "orange", "Delta+1.52 (Schur)"),
                      (delta, "gray", "Delta (on-site)")]:
        ax.axvline(x, color=c, ls="--", lw=0.8, label=lab)
    ax.set_xlabel("E / t"); ax.set_ylabel("path LDOS")
    ax.set_title(f"KPM path LDOS  gap_above={f['gap_above']:+.2f}t  conc@D+phi={f['conc_phys']*100:.1f}% "
                 f"({f['conc_over_flat']:.1f}x flat)  FWHM={f['pk_fwhm_mev']:.0f} MeV")
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(fname, dpi=110)
    print(f"# saved plot {fname}")


# ----------------------------------------------------------------------------- validation
def validate(graph, backend="cpu", L=6, delta=1.78, Nc=400, t_mev=332.0, t_leak=0.4142135623730951,
             sublattice=False, gauge="sc"):
    """Re-derive a KNOWN case: small lattice where full eigh is exact. Assert the KPM recurrence
    (on the REQUESTED backend) reproduces the exact spectral moments and the reconstructed LDOS
    integrates to 1. Exit 0 here == the KPM engine is correct, independent of the physics question."""
    Lg = 4 if graph == "matter" and sublattice else L      # sublattice cell count grows fast; keep tiny
    print(f"[validate] graph={graph} gauge={gauge} backend={backend} L={Lg} delta={delta} Nc={Nc} t_leak={t_leak:.4f} sub={sublattice}")
    A, path = build_graph(graph, Lg, t_leak, sublattice, gauge)
    H = add_defect(A, path, delta)
    N = H.shape[0]
    emin, emax, _, _ = spectral_bounds(H, BAND[graph], delta)
    a = (emax - emin) / 2.0
    b = (emax + emin) / 2.0

    # exact reference via dense eigh
    Hd = H.toarray()
    evals, evecs = np.linalg.eigh(Hd)
    assert evals.min() > emin - 1e-9 and evals.max() < emax + 1e-9, "bounds must bracket spectrum"
    xk = (evals - b) / a
    assert np.all(np.abs(xk) < 1.0), f"rescaled spectrum must be inside (-1,1): max|x|={np.abs(xk).max()}"
    K = len(path)
    wj = (np.abs(evecs[path, :]) ** 2).sum(0) / K          # path-averaged spectral weight per eigenstate
    n = np.arange(Nc)
    Tnk = np.cos(np.outer(np.arccos(xk), n))               # T_n(x_k)
    mu_exact = (wj[:, None] * Tnk).sum(0)                  # <path|T_n|path>/K, exact

    mu_kpm = kpm_local_moments(H, path, a, b, Nc, backend=backend, report=False)
    err = float(np.max(np.abs(mu_kpm - mu_exact)))
    print(f"[validate] N={N}  backend={backend}  max|mu_KPM - mu_exact| = {err:.2e}")
    assert err < 1e-9, f"KPM moments disagree with exact eigh: {err}"

    E, rho = reconstruct(mu_kpm, a, b, npts=4000)
    integral = float(np.trapezoid(rho, E))
    print(f"[validate] integral of reconstructed path-LDOS = {integral:.6f}  (must be ~1)")
    assert abs(integral - 1.0) < 5e-3, f"LDOS not normalised: {integral}"
    assert mu_kpm[0] > 0.999 and mu_kpm[0] < 1.001, f"mu_0 must be 1: {mu_kpm[0]}"
    print("[validate] PASS: KPM engine reproduces exact moments and a normalised LDOS.")
    return True


# ----------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graph", choices=["sc", "linegraph", "matter"], default="linegraph")
    ap.add_argument("--L", type=int, default=64)
    ap.add_argument("--delta", type=float, default=1.78, help="on-site confinement gap on the P4 path (Delta/t)")
    ap.add_argument("--t-leak", type=float, default=2 ** 0.5 - 1,
                    help="matter mode: gauge-matter triangular-face channel hopping (default silver ratio sqrt2-1)")
    ap.add_argument("--sublattice", action="store_true",
                    help="matter mode: add the full 3-bipyramid/cell matter sublattice (spectator-invariance check)")
    ap.add_argument("--gauge", choices=["sc", "linegraph"], default="sc",
                    help="matter mode: the gauge continuum the rho leaks into -- 'sc' [-6,6] proxy or the real L(SC) [-2,10]")
    ap.add_argument("--moments", type=int, default=8000)
    ap.add_argument("--backend", choices=["auto", "cpu", "gpu"], default="auto")
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--t-mev", type=float, default=332.0, help="t ~= Lambda_QCD in MeV (item48 uses 332)")
    ap.add_argument("--egrid", type=int, default=20000)
    ap.add_argument("--validate", action="store_true", help="run the self-test and exit")
    ap.add_argument("--no-validate", action="store_true", help="skip the self-test before the production run")
    args = ap.parse_args()

    phi = (1 + 5 ** 0.5) / 2

    if args.backend == "auto":
        try:
            import cupy  # noqa: F401
            backend = "gpu"
        except Exception:
            backend = "cpu"
    else:
        backend = args.backend
    print(f"# backend={backend}  graph={args.graph}  L={args.L}  delta={args.delta}  moments={args.moments}")

    if args.validate:
        validate(args.graph, backend, t_leak=args.t_leak, sublattice=args.sublattice, gauge=args.gauge)
        print("VALIDATION OK (exit 0).")
        return 0

    if not args.no_validate:
        validate(args.graph, backend, t_leak=args.t_leak, sublattice=args.sublattice, gauge=args.gauge)

    # ---- production ----
    if args.graph == "matter":
        print(f"# matter mode: gauge continuum={args.gauge} (sc=[-6,6] proxy, linegraph=real L(SC) [-2,10]); "
              f"t_leak={args.t_leak:.4f} (silver ratio sqrt2-1={2**0.5-1:.4f}); sublattice={args.sublattice}")
    t0 = time.time()
    A, path = build_graph(args.graph, args.L, args.t_leak, args.sublattice, args.gauge)
    H = add_defect(A, path, args.delta)
    N = H.shape[0]
    nnz = H.nnz
    print(f"# H built: N={N:,}  nnz={nnz:,}  ({time.time()-t0:.1f}s)  est CSR fp64+i32 = {nnz*12/1e9:.2f} GB")
    emin, emax, raw_emin, raw_emax = spectral_bounds(H, BAND[args.graph], args.delta)
    a = (emax - emin) / 2.0
    b = (emax + emin) / 2.0
    print(f"# spectral bounds [{emin:.3f},{emax:.3f}]  a={a:.4f} b={b:.4f}  "
          f"phi={phi:.4f} in-band={emin < phi < emax}  Delta+phi={args.delta+phi:.4f} in-band={emin < args.delta+phi < emax}")
    assert a > 0, f"non-positive scale a={a}"
    # fail-fast: rescaled spectral radius must be < 1 or the Chebyshev recurrence diverges to NaN
    rng = np.random.default_rng(0); v = rng.standard_normal(N); v /= np.linalg.norm(v)
    rad = 0.0
    for _ in range(40):
        w = (H @ v - b * v) / a; rad = float(np.linalg.norm(w)); v = w / rad
    print(f"# rescaled spectral radius est = {rad:.4f}  (must be < 1)")
    assert rad < 1.0 + 1e-6, f"rescaled spectral radius {rad:.4f} >= 1 -> bounds too tight; widen pad"

    tk = time.time()
    mu = kpm_local_moments(H, path, a, b, args.moments, backend=backend, device=args.device)
    print(f"# {args.moments} path moments in {time.time()-tk:.1f}s "
          f"({(time.time()-tk)/args.moments*1e3:.2f} ms/moment)")
    if backend == "gpu":                                  # free GPU pool before the bulk run
        import cupy as _cp
        _cp.get_default_memory_pool().free_all_blocks()

    # defect-free bulk per-site DOS (stochastic trace, Kb random unit vectors) -- the continuum reference.
    # For matter mode the continuum the rho leaks into is the chosen GAUGE web (sc or the real L(SC)).
    A_bulk = (build_sc(args.L)[0] if args.gauge == "sc" else build_linegraph(args.L)[0]) \
        if args.graph == "matter" else A
    Nb = A_bulk.shape[0]
    Kb = 16 if Nb < 2_000_000 else 4
    rng = np.random.default_rng(1)
    rad = np.sign(rng.standard_normal((Nb, Kb))).astype(np.float64)
    rad /= np.linalg.norm(rad, axis=0, keepdims=True)
    tb = time.time()
    mu_bulk = kpm_local_moments(A_bulk, [0, 1, 2, 3], a, b, args.moments, backend=backend,
                                device=args.device, start=rad, report=False)
    print(f"# {args.moments} bulk moments in {time.time()-tb:.1f}s")

    E, rho = reconstruct(mu, a, b, npts=args.egrid)
    _, rho_bulk = reconstruct(mu_bulk, a, b, npts=args.egrid)
    xband = EXTRACT_BAND[args.gauge] if args.graph == "matter" else EXTRACT_BAND[args.graph]
    f = extract(E, rho, args.delta, phi, args.t_mev, xband, raw_emin, raw_emax)
    eta = np.pi * a / args.moments
    print("\n=== path-projected LDOS features (graph=%s L=%d Delta=%.3f Nc=%d) ===" %
          (args.graph, args.L, args.delta, args.moments))
    print(f"  KPM resolution eta ~ pi a/Nc        = {eta:.4f} t = {eta*args.t_mev:.1f} MeV "
          f"(<< 149 MeV: a real narrow rho WOULD be resolved)")
    print(f"  clean band [{f['band_lo']:.1f},{f['band_hi']:.1f}] t  would-be rho Delta+phi = {f['E_res']:.3f} t "
          f"({f['E_res']*args.t_mev:.0f} MeV) -> {'INSIDE -> resonance' if f['E_res_in_band'] else 'OUTSIDE band'}")
    print(f"  [trivial moments] centroid {f['centroid']:.3f} (rel {f['centroid_rel']:+.3f}=mu1-Delta), "
          f"sigma {f['sigma']:.3f}t=sqrt(deg) -> total spread {2.355*f['sigma_mev']:.0f} MeV")
    print("  -- physics discriminators --")
    print(f"  (1) bound state above band? raw lam_max-band_top = {f['gap_above']:+.3f} t  "
          f"(needs O(0.1+) gap for a split-off state; ~0 => none)")
    print(f"  (2) weight in physical {args.t_mev*f['gw']:.0f}-MeV window at Delta+phi = {f['conc_phys']*100:.2f}%  "
          f"(= {f['conc_over_flat']:.2f}x a flat smear; ~1x = NO enhancement)")
    print(f"  (3) phi-level peak at Delta+phi      = {f['Epk']:.3f} t  (rel to Delta: {f['Epk_rel']:+.3f}; cf phi={phi:.3f}, 1.52)")
    print(f"      FWHM (smoothed {f['smooth_mev']:.0f} MeV)      = {f['pk_fwhm_mev']:.0f} MeV ; deconvolved ~{f['pk_fwhm_deconv_mev']:.0f} MeV  (physical rho 149)")

    # verdict (REPORTED, not asserted). Survival = the phi-level keeps a substantial fraction of the
    # string weight as a NARROW (<~300 MeV) resonance at Delta+phi, OR a split-off bound state.
    survives = (f["gap_above"] > 0.15) or (f["conc_over_flat"] > 3.0 and f["pk_fwhm_deconv_mev"] < 300
                                           and abs(f["Epk_rel"] - phi) < 0.45)
    print("\nVERDICT (reported):")
    if survives:
        print(f"  SURVIVES: the phi-level keeps {f['conc_phys']*100:.0f}% of the string weight ({f['conc_over_flat']:.1f}x a flat")
        print(f"  smear) as a resonance at rel {f['Epk_rel']:+.2f} (cf phi={phi:.3f}) of width ~{f['pk_fwhm_deconv_mev']:.0f} MeV")
        print("  (empirical 149). A distinct vector meson IS present -- the embedded-mode dissolution was an")
        print("  artifact of t_leak=1; at the silver-ratio gauge-matter coupling the rho is a finite-width object.")
    else:
        print(f"  DISSOLVED: gap {f['gap_above']:+.2f}t, only {f['conc_phys']*100:.1f}% of the string weight in a 149-MeV")
        print(f"  window at Delta+phi ({f['conc_over_flat']:.1f}x flat), feature ~{f['pk_fwhm_deconv_mev']:.0f} MeV wide. No distinct narrow rho.")
    tag = f"_{args.gauge}_tl{args.t_leak:.3f}{'_sub' if args.sublattice else ''}" if args.graph == "matter" else ""
    out = f"kpm_{args.graph}_L{args.L}_D{args.delta}_Nc{args.moments}{tag}"
    np.savez(out + ".npz", E=E, rho=rho, rho_bulk=rho_bulk, mu=mu, mu_bulk=mu_bulk, a=a, b=b,
             **{k: v for k, v in f.items()})
    make_plot(E, rho, rho_bulk, f, args.delta, phi, out + ".png")
    print(f"# saved {out}.npz")
    return 0


if __name__ == "__main__":
    sys.exit(main())
