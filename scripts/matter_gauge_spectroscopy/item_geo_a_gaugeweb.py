#!/usr/bin/env python
"""
item_geo_a_gaugeweb.py -- Experiment (A): entanglement geometry of the substrate GAUGE WEB.

Canon (ANCHOR sec 7.3): the gauge web is a compact-U(1) Wilson field on the cubic dual lattice with
kernel K(k) = 6 - 2(cos kx + cos ky + cos kz). In its free/Gaussian form this is an exactly solvable
lattice scalar -- so we reach large L by the correlation-matrix method and finally give the curvature
meter the dynamic range ED could not.

Tested on the ACTUAL canon kernel (massless, IR-regulated) vs a massive contrast:
  [1] RT kinematics: AREA LAW -- S(block) ~ boundary area (l^2), not volume (l^3).
  [2] EMERGENT GEOMETRY: d(r) = -log I(i:j). d ~ log r => scale-free (AdS-like emergent/RG bulk);
      d ~ r => flat. HONEST FRAMING: this is the emergent RG/holographic geometry reconstructed from
      entanglement, NOT the physical flat-Z^3 space the substrate lives on -- the two differ, exactly
      as a flat-space CFT carries a hyperbolic AdS bulk.
  [3] ISOTROPY of MI across lattice directions (reconstruction sanity).

Free real scalar H = 1/2 sum pi^2 + 1/2 phi^T K phi. Ground state X=<phi phi>=K^{-1/2}/2,
P=<pi pi>=K^{1/2}/2. Region entropy = symplectic eigenvalues of X_A P_A (Casini-Huerta / Peschel).
"""
import json
import numpy as np


def correlators(L, m2):
    k = 2 * np.pi * np.arange(L) / L
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    Kk = m2 + 6 - 2 * (np.cos(KX) + np.cos(KY) + np.cos(KZ))
    Xr = np.fft.ifftn(0.5 * Kk ** -0.5).real
    Pr = np.fft.ifftn(0.5 * Kk ** 0.5).real
    return Xr, Pr


def sympl_entropy(XA, PA):
    ev = np.linalg.eigvals(XA @ PA)
    nu = np.sqrt(np.clip(ev.real, 0, None))
    nu = np.clip(nu, 0.5 + 1e-12, None)
    return float(np.sum((nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)))


def block_entropy(Xr, Pr, ell, L):
    sites = [(x, y, z) for x in range(ell) for y in range(ell) for z in range(ell)]
    n = len(sites)
    XA = np.empty((n, n)); PA = np.empty((n, n))
    for a, (x1, y1, z1) in enumerate(sites):
        for b, (x2, y2, z2) in enumerate(sites):
            XA[a, b] = Xr[(x1 - x2) % L, (y1 - y2) % L, (z1 - z2) % L]
            PA[a, b] = Pr[(x1 - x2) % L, (y1 - y2) % L, (z1 - z2) % L]
    return sympl_entropy(XA, PA)


def site_entropy(Xr, Pr):
    return sympl_entropy(np.array([[Xr[0, 0, 0]]]), np.array([[Pr[0, 0, 0]]]))


def mi_pair(Xr, Pr, dr, L, S1):
    dx, dy, dz = dr
    X0, Xd = Xr[0, 0, 0], Xr[dx % L, dy % L, dz % L]
    P0, Pd = Pr[0, 0, 0], Pr[dx % L, dy % L, dz % L]
    S2 = sympl_entropy(np.array([[X0, Xd], [Xd, X0]]), np.array([[P0, Pd], [Pd, P0]]))
    return max(2 * S1 - S2, 1e-15)


def linfit_r2(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    yh = A @ c
    ss = np.sum((y - yh) ** 2); st = np.sum((y - y.mean()) ** 2)
    return float(c[0]), float(1 - ss / st if st > 0 else 0)


def run(L, m2, tag):
    Xr, Pr = correlators(L, m2)
    S1 = site_entropy(Xr, Pr)

    ells = list(range(2, 9))
    Sb = [block_entropy(Xr, Pr, e, L) for e in ells]
    _, r2_area = linfit_r2([e ** 2 for e in ells], Sb)
    _, r2_vol = linfit_r2([e ** 3 for e in ells], Sb)

    rs = list(range(1, L // 2))
    Ir = [mi_pair(Xr, Pr, (r, 0, 0), L, S1) for r in rs]
    keep = [(r, I) for r, I in zip(rs, Ir) if I > 1e-12]
    rk = [r for r, _ in keep]; dk = [-np.log(I) for _, I in keep]
    s_log, r2_log = linfit_r2(np.log(rk), dk)
    s_lin, r2_lin = linfit_r2(rk, dk)
    enough = len(keep) >= 4
    geom = ("SCALE-FREE  d~log r  => AdS-like emergent bulk" if (enough and r2_log > r2_lin)
            else "FLAT  d~r" if enough else f"gapped: only {len(keep)} pts (flat by default)")

    iso = {"axis_(3,0,0)|r|=3.00": mi_pair(Xr, Pr, (3, 0, 0), L, S1),
           "face_(2,2,0)|r|=2.83": mi_pair(Xr, Pr, (2, 2, 0), L, S1),
           "body_(2,2,1)|r|=3.00": mi_pair(Xr, Pr, (2, 2, 1), L, S1)}

    return {"L": L, "m2": m2, "tag": tag, "S_site": S1,
            "area_R2": r2_area, "vol_R2": r2_vol, "Sblocks": dict(zip(ells, Sb)),
            "geom_logr_R2": r2_log, "geom_r_R2": r2_lin, "geom_logr_slope": s_log,
            "n_pts": len(keep), "geom_verdict": geom,
            "MI_r": {r: float(I) for r, I in zip(rs, Ir)}, "isotropy": iso}


def main():
    L = 32
    out = {}
    for m2, tag in [(1e-4, "massless_canon_kernel"), (1.0, "massive_contrast")]:
        res = run(L, m2, tag)
        out[tag] = res
        print(f"\n=== {tag}  (L={L}, m^2={m2}) ===")
        print(f"  S(single site) = {res['S_site']:.4f}")
        print(f"  [1] AREA LAW : S(block) vs l^2 R2={res['area_R2']:.4f}  vs l^3 R2={res['vol_R2']:.4f}"
              f"  -> {'AREA' if res['area_R2'] > res['vol_R2'] else 'VOLUME'} law")
        print(f"      S(l) for l=2..8: " + " ".join(f"{v:.2f}" for v in res['Sblocks'].values()))
        print(f"  [2] EMERGENT GEOM: d=-logI vs log r R2={res['geom_logr_R2']:.4f}  vs r R2={res['geom_r_R2']:.4f}"
              f"  ({res['n_pts']} pts, log-slope={res['geom_logr_slope']:.3f})")
        print(f"      -> {res['geom_verdict']}")
        iso = res['isotropy']
        print("  [3] ISOTROPY MI at |r|~3: " + "  ".join(f"{k}={v:.3e}" for k, v in iso.items()))
    with open("item_geo_a_gaugeweb_results.json", "w") as f:
        json.dump(out, f, indent=2, default=float)
    print("\nsaved -> item_geo_a_gaugeweb_results.json")


if __name__ == "__main__":
    main()
