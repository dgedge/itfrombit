#!/usr/bin/env python3
r"""SMG/TCH continuum lift -- phase 2: finite-VOLUME axis (generic Z3 lattice).

Generalises the one-plaquette cell to a multi-plaquette strip so the mirror gap can
be tracked vs volume L (number of plaquettes). Same local mirror-Fock block, Z3 open
links, Gauss-law basis, electric + magnetic + gauge-covariant hopping, matrix-free
Lanczos. Geometries: L=1 (4 sites, cross-check vs phase 1) and L=2 strip (6 sites, 7
links, 2 plaquettes). The mirror gap is the gap to the first state with matter
excitation > 0.5 above the ground. Goal: show the gap does not close as L grows.

numpy only; writes JSON.
"""
import importlib.util, json, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("block", HERE / "css_spatial_mirror_fock_gauge_gap.py")
block = importlib.util.module_from_spec(spec); spec.loader.exec_module(block)

# (sites, edges as (tail,head), plaquettes as (link-indices, signs))
GEOM = {
    1: dict(sites=4,
            edges=[(0, 1), (1, 2), (3, 2), (0, 3)],
            plaqs=[((0, 1, 2, 3), (1, 1, -1, -1))]),
    2: dict(sites=6,
            edges=[(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)],
            plaqs=[((0, 5, 2, 4), (1, 1, -1, -1)), ((1, 6, 3, 5), (1, 1, -1, -1))]),
}

def electric_energy(flux): return 0.0 if flux % 3 == 0 else 1.0

def gauss_ok(charges, fluxes, edges, sites):
    for v in range(sites):
        val = charges[v]
        for link, (tail, head) in enumerate(edges):
            if tail == v: val += fluxes[link]
            if head == v: val -= fluxes[link]
        if val % 3 != 0: return False
    return True

def build_basis(charges, edges, sites):
    ld = len(charges); nl = len(edges)
    basis = []
    for config in np.ndindex(*(ld for _ in range(sites))):
        ct = tuple(int(charges[s]) for s in config)
        if sum(ct) % 3 != 0: continue
        for fluxes in np.ndindex(*([3] * nl)):
            if gauss_ok(ct, fluxes, edges, sites):
                basis.append((tuple(config), tuple(int(f) for f in fluxes)))
    return basis, {it: i for i, it in enumerate(basis)}

def transition_tables(annihilators):
    ann, create = [], []
    for mode in annihilators:
        ab, cb = [], []
        for ch, M in enumerate(mode):
            asrc, csrc = [], []
            D = M.conj().T
            for s in range(M.shape[1]):
                asrc.append([(int(t), M[t, s]) for t in np.nonzero(np.abs(M[:, s]) > 1e-12)[0]])
                csrc.append([(int(t), D[t, s]) for t in np.nonzero(np.abs(D[:, s]) > 1e-12)[0]])
            ab.append(asrc); cb.append(csrc)
        ann.append(ab); create.append(cb)
    return ann, create

class LatticeH:
    def __init__(self, L, beta, hopping, sp):
        g = GEOM[L]; self.edges = g["edges"]; self.plaqs = g["plaqs"]; self.sites = g["sites"]
        self.beta = beta; self.hopping = hopping; self.g2 = 1.0 / beta
        self.energies, self.charges, self.ann_mats, _ = block.local_mirror_dressed_block(sp)
        self.basis, self.index = build_basis(self.charges, self.edges, self.sites)
        self.dim = len(self.basis)
        self.ann, self.create = transition_tables(self.ann_mats)
        self.diag = np.zeros(self.dim); self.mdiag = np.zeros(self.dim)
        for p, (cfg, fl) in enumerate(self.basis):
            mt = sum(self.energies[s] for s in cfg)
            self.mdiag[p] = mt
            self.diag[p] = mt + self.g2 * sum(electric_energy(f) for f in fl)
        self.trans = self._build()

    def _pshift(self, fluxes, plaq, amount):
        out = list(fluxes); links, signs = plaq
        for lk, sg in zip(links, signs): out[lk] = (out[lk] + amount * sg) % 3
        return tuple(out)

    def _fshift(self, fluxes, link, amount):
        out = list(fluxes); out[link] = (out[link] + amount) % 3; return tuple(out)

    def _build(self):
        T = []
        for col, (cfg, fl) in enumerate(self.basis):
            ct = []
            for plaq in self.plaqs:
                for amt in (1, -1):
                    r = self.index.get((cfg, self._pshift(fl, plaq, amt)))
                    if r is not None: ct.append((r, -(self.beta / 2)))
            if abs(self.hopping) > 1e-15:
                cl = list(cfg)
                for link, (tail, head) in enumerate(self.edges):
                    ot, oh = cl[tail], cl[head]
                    for mode in range(len(self.ann)):
                        for ch in range(3):
                            for nt, at in self.create[mode][ch][ot]:
                                for nh, ah in self.ann[mode][ch][oh]:
                                    nc = list(cl); nc[tail] = nt; nc[head] = nh
                                    r = self.index.get((tuple(nc), self._fshift(fl, link, ch)))
                                    if r is not None: ct.append((r, -self.hopping * at * ah))
                            for nt, at in self.ann[mode][ch][ot]:
                                for nh, ah in self.create[mode][ch][oh]:
                                    nc = list(cl); nc[tail] = nt; nc[head] = nh
                                    r = self.index.get((tuple(nc), self._fshift(fl, link, -ch)))
                                    if r is not None: ct.append((r, -self.hopping * at * ah))
            T.append(ct)
        return T

    def matvec(self, v):
        out = self.diag * v
        for col, ct in enumerate(self.trans):
            a = v[col]
            if abs(a) < 1e-15: continue
            for r, c in ct: out[r] += c * a
        return out

    def matter(self, v): return float(np.real(np.vdot(v, self.mdiag * v)))

def lanczos(model, n_iter=110, n_eigs=16, seed=143):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=model.dim) + 1j * rng.normal(size=model.dim); q /= np.linalg.norm(q)
    qprev = np.zeros_like(q); bprev = 0.0; Q = []; al = []; be = []
    for _ in range(min(n_iter, model.dim)):
        Q.append(q); v = model.matvec(q); a = float(np.real(np.vdot(q, v)))
        v = v - a * q - bprev * qprev
        for bv in Q: v -= np.vdot(bv, v) * bv
        b = float(np.linalg.norm(v)); al.append(a)
        if b < 1e-12: break
        be.append(b); qprev = q; q = v / b; bprev = b
    k = len(al); Tm = np.diag(al)
    for i in range(k - 1): Tm[i, i + 1] = be[i]; Tm[i + 1, i] = be[i]
    ev, ec = np.linalg.eigh(Tm); R = np.column_stack(Q) @ ec[:, :n_eigs]
    return ev[:n_eigs], R

def mirror_gap(model):
    ev, R = lanczos(model)
    g0 = model.matter(R[:, 0]); mir = float("inf")
    for i in range(1, R.shape[1]):
        if model.matter(R[:, i]) - g0 > 0.5: mir = float(ev[i] - ev[0]); break
    return float(ev[1] - ev[0]), mir

def main():
    results = []
    plan = [(1, 2), (1, 3), (2, 2), (2, 3)]   # (L, sp)
    for L, sp in plan:
        for beta in (0.5, 1.0, 2.0):
            for t in (0.0, 1.0):
                t0 = time.time()
                m = LatticeH(L, beta, t, sp)
                full, mir = mirror_gap(m); dt = time.time() - t0
                row = dict(L=L, sp=sp, beta=beta, t=t, dim=int(m.dim),
                           full_gap=full, mirror_gap=(None if mir == float("inf") else mir), sec=dt)
                results.append(row)
                print(f"  L={L} sp={sp} beta={beta:<4} t={t:<4} dim={m.dim:<8} "
                      f"full={full:.4f} mirror={'--' if mir==float('inf') else f'{mir:.4f}'} ({dt:.1f}s)", flush=True)
                json.dump(results, open(HERE / "smg_volume_results.json", "w"), indent=2)
    # volume comparison at matched (sp,beta,t)
    print("\n  mirror gap vs volume L (matched sp,beta,t):")
    for sp in (2, 3):
        for beta in (0.5, 1.0, 2.0):
            for t in (0.0, 1.0):
                row = {r["L"]: r["mirror_gap"] for r in results if r["sp"] == sp and r["beta"] == beta and r["t"] == t}
                if 1 in row and 2 in row and row[1] and row[2]:
                    print(f"    sp={sp} beta={beta:<4} t={t:<4}: L1={row[1]:.3f}  L2={row[2]:.3f}  (ratio {row[2]/row[1]:.3f})")
    print("\nDONE.")

if __name__ == "__main__":
    main()
