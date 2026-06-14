#!/usr/bin/env python3
"""
Approach B gold-standard test: genuine multi-cell fermion-Fock ED of the SU(4) Pati-Salam
SMG mirror gap vs spatial hopping. NO hardcoding -- real creation/annihilation operators,
CSS X/Z stabilisers, sparse eigsh.

Guardrails (asserted): t=0 (decoupled cells) -> gap 2.0 (full Fock AND mirror sector);
the hopping genuinely competes (||[H_hop, X-stab]|| != 0).

Result: the SMG mirror gap COLLAPSES under hopping for SU(4) (2.0 -> ~0.19 at t=4 -> ~0.085 at
t=8), no better than SU(3) (-> ~0.09 -> ~0.04); at intermediate t SU(4) collapses faster. This
is full-Fock CONVERGED (no n/q truncation) for a 2-cell system -> the n/q->inf limit at that
size; it matches the n/q=6 strip collapse.

Scope note: a *proper Z4-gauged* strip (Gauss law) needs the F2-register <-> Z4-colour-charge
map, which is the OPEN section-2.8 problem -> not cleanly definable yet. The gauge-covariance
question that a gauged strip would test is settled operator-exactly in su4_smg_real_test.py
(X-stab fails the colour-centre commutation for SU(4), = SU(3)).

Needs numpy + scipy.
"""
import numpy as np
import scipy.sparse as sp
from scipy.sparse import identity, csr_matrix
from scipy.sparse.linalg import eigsh

NB = 8
G0, G1, LQ, C0, C1, I3, CHI, W = range(NB)
GEN = [[1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]]
MATTER = [G0, G1, LQ, C0, C1, I3]   # hopped modes; CHI,W preserved -> stays in mirror sector


def fermion_ops(M):
    DIM = 1 << M
    bits = [tuple((s >> (M - 1 - b)) & 1 for b in range(M)) for s in range(DIM)]
    index = {b: i for i, b in enumerate(bits)}
    c = []; g = []; P = []
    for k in range(M):
        ra, ca, va = [], [], []
        rg, cg, vg = [], [], []
        dP = np.ones(DIM)
        for s, b in enumerate(bits):
            sign = (-1) ** sum(b[:k]); t = list(b); t[k] ^= 1; tt = index[tuple(t)]
            rg.append(tt); cg.append(s); vg.append(sign)
            if b[k] == 1:
                ra.append(tt); ca.append(s); va.append(sign)
            dP[s] = 1 - 2 * b[k]
        c.append(csr_matrix((va, (ra, ca)), shape=(DIM, DIM)))
        g.append(csr_matrix((vg, (rg, cg)), shape=(DIM, DIM)))
        P.append(csr_matrix((dP, (range(DIM), range(DIM))), shape=(DIM, DIM)))
    cd = [x.getH() for x in c]
    return c, cd, g, P, bits


def cell_H(g, P, cell):
    off = cell * NB; H = None
    for row in GEN:
        sup = [off + b for b in range(NB) if row[b]]
        X = identity(g[0].shape[0], format='csr'); Z = identity(P[0].shape[0], format='csr')
        for b in sup: X = X @ g[b]
        for b in sup: Z = Z @ P[b]
        H = (X + Z) if H is None else H + X + Z
    return -H


def hopping(c, cd, A, B):
    H = None
    for b in MATTER:
        a, bb = A * NB + b, B * NB + b
        term = cd[a] @ c[bb]; term = term + term.getH()
        H = term if H is None else H + term
    return H


def lowgap(H, k=6):
    ev = np.sort(eigsh(H, k=k, which='SA', return_eigenvectors=False).real); lv = []
    for v in ev:
        if not lv or v - lv[-1] > 1e-6: lv.append(v)
    return lv[1] - lv[0] if len(lv) > 1 else float('nan')


def main():
    ncells = 2; M = NB * ncells
    print(f"building genuine fermion ops on M={M} modes (dim={1 << M}) ...", flush=True)
    c, cd, g, P, bits = fermion_ops(M)
    Hc = sum(cell_H(g, P, cell) for cell in range(ncells))
    g0 = lowgap(Hc); print(f"GUARDRAIL t=0 full Fock gap = {g0:.4f} (expect 2.0)")
    assert abs(g0 - 2.0) < 1e-6

    Hh = hopping(c, cd, 0, 1)
    xsA = cell_H(g, P, 0)  # not the stabiliser itself, but check competition via a stabiliser
    # competition check against one X-stabiliser:
    Xs = identity(1 << M, format='csr')
    for b in [b for b in range(NB) if GEN[1][b]]:
        Xs = Xs @ g[b]
    comp = abs((Hh @ Xs - Xs @ Hh)).max()
    print(f"[check] ||[H_hop, X-stab]|| = {comp:.3f} (nonzero -> competes)")
    assert comp > 1e-6

    mir = [i for i, b in enumerate(bits) if b[CHI] != b[W] and b[NB + CHI] != b[NB + W]]
    Pm = csr_matrix((np.ones(len(mir)), (range(len(mir)), mir)), shape=(len(mir), 1 << M))
    R = lambda Op: (Pm @ Op @ Pm.getH()).tocsr()
    Hc_m, Hh_m = R(Hc), R(Hh)
    mb = [bits[i] for i in mir]
    p3 = np.array([0.0 if ((x[C0], x[C1]) == (0, 0) or (x[NB + C0], x[NB + C1]) == (0, 0)) else 1.0 for x in mb])
    P3 = sp.diags(p3)
    Hh_m_su3 = (P3 @ Hh_m @ P3).tocsr()

    gm0 = lowgap(Hc_m); print(f"GUARDRAIL t=0 mirror gap = {gm0:.4f} (expect 2.0)")
    assert abs(gm0 - 2.0) < 1e-6

    print(f"\n  {'t':>5} {'SU(4) gap':>10} {'SU(3) gap':>10}")
    res = {}
    for t in [0.0, 1.0, 2.0, 4.0, 8.0]:
        g4 = lowgap(Hc_m + t * Hh_m); g3 = lowgap(Hc_m + t * Hh_m_su3)
        res[t] = (g4, g3); print(f"  {t:>5.1f} {g4:>10.4f} {g3:>10.4f}")
    assert res[4.0][0] < 0.5 and res[4.0][1] < 0.5   # both collapse
    print("\nRESULT: SU(4) SMG mirror gap COLLAPSES under hopping (2.0 -> ~0.19 at t=4),")
    print("no better than SU(3). Pati-Salam route does NOT close the frontier.")
    print("ALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
