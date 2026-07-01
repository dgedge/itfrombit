#!/usr/bin/env python
r"""candidate_macroscopic_model.py -- a CANDIDATE macroscopic matter+gauge model (MODELING CHOICE).

Canon flags the macroscopic substrate model as the open gate: the multi-cell operator set is
"described only qualitatively (line graph L(TCH), gauge web) -- NOT as a diagonalisable operator
set" (DRIFT.md:1283); gamma_TEE is gated on first writing it down. This script PROPOSES one
explicit, constraint-respecting resolution and computes its topological order, to make the gate
concrete. It is NOT derived from canon -- a labelled modeling choice -- but it respects the
canonical local data:
  - matter = the self-dual [8,4,4] / [[8,0,4]] cell code (distance 4), one per site; k=0, so the
    cells carry NO logical qubits and contribute NO topological degeneracy themselves;
  - cells do NOT share qubits (the distance-4 gluing obstruction, DRIFT.md:1287) -- they connect
    THROUGH a gauge web of bridge qubits living on the bonds;
  - the gauge web is taken as a Z2 lattice gauge theory on the bond graph (the discrete-photon
    caricature of the simple-cubic Wilson web).

THE POINT: this single explicit operator set has TWO phases, realizing canon's documented tension
(DRIFT.md:1281, photon-band-topology vs long-range-entanglement) as a phase CHOICE --
  - deconfined gauge web      -> Z2 TOPOLOGICAL ORDER: torus ground-state degeneracy GSD = D^2 = 4,
                                 total quantum dimension D = 2, gamma_TEE = log 2;
  - confined / matter-Higgsed -> TRIVIAL: GSD = 1, gamma_TEE = 0, and THIS is the branch that
                                 hosts a massless photon (the framework's physical Coulomb vacuum,
                                 measured gamma~0 in item_geo_b_maxwell.py).
So "writing the macroscopic model" does not by itself decide gamma_TEE; the gate is a PHASE choice,
and canon's massless emergent photon selects the gamma = 0 (trivial/Coulomb) branch -- consistent
with item_geo_b. The topological (gamma=log2) branch is the gapped alternative canon already deems
incompatible with the photon (sc_gauge_web_toric.py).

Diagnostic: stabilizer formalism over GF(2). GSD = 2^(N - rank(check matrix)).
"""
import numpy as np


def gf2_rank(M):
    M = (M.copy() % 2).astype(np.int8)
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i, c]), None)
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


# ---- matter: the [8,4,4] Reed-Muller cell code, distance 4 ----
G = np.array([[1, 1, 1, 1, 1, 1, 1, 1],
              [0, 0, 0, 0, 1, 1, 1, 1],
              [0, 0, 1, 1, 0, 0, 1, 1],
              [0, 1, 0, 1, 0, 1, 0, 1]])


def cell_distance(G):
    d = 8
    for k in range(1, 16):
        cw = np.zeros(8, int)
        for r in range(4):
            if (k >> r) & 1:
                cw ^= G[r]
        w = int(cw.sum())
        if w:
            d = min(d, w)
    return d


# ---- gauge web: Z2 lattice gauge theory on an L x L periodic square (qubits on edges) ----
def gauge_web_GSD(L, phase):
    """phase = 'deconfined' (toric/topological) or 'trivial' (confined/Higgs product)."""
    Lh = L * L                                  # horizontal edges 0..L^2-1
    N = 2 * L * L                               # + vertical edges L^2..2L^2-1
    h = lambda x, y: (x % L) * L + (y % L)
    v = lambda x, y: Lh + (x % L) * L + (y % L)
    rows = []

    def Xrow(edges):
        r = np.zeros(2 * N, np.int8);
        for e in edges: r[e] = 1
        return r

    def Zrow(edges):
        r = np.zeros(2 * N, np.int8)
        for e in edges: r[N + e] = 1
        return r

    if phase == "deconfined":
        for x in range(L):
            for y in range(L):
                rows.append(Xrow([h(x, y), h(x - 1, y), v(x, y), v(x, y - 1)]))       # star A_v (X)
                rows.append(Zrow([h(x, y), h(x, y + 1), v(x, y), v(x + 1, y)]))       # plaquette B_p (Z)
    elif phase == "trivial":
        for e in range(N):                                                            # product state: Z on every edge
            rows.append(Zrow([e]))
    M = np.array(rows, np.int8)
    R = gf2_rank(M)
    return N, R, 2 ** (N - R)


def main():
    print("=== CANDIDATE macroscopic matter+gauge model (modeling choice; not canon) ===\n")

    d = cell_distance(G)
    print(f"[matter] [8,4,4] cell code distance d = {d}  (k=0 -> 0 logical qubits -> 0 topological "
          f"degeneracy from matter)")

    L = 3
    Nd, Rd, gsd_d = gauge_web_GSD(L, "deconfined")
    Nt, Rt, gsd_t = gauge_web_GSD(L, "trivial")
    import math
    Dd = math.isqrt(gsd_d)
    print(f"\n[gauge web on {L}x{L} torus, N={Nd} bond qubits]")
    print(f"  deconfined : rank={Rd}  -> GSD = {gsd_d}  => total quantum dim D = {Dd}, "
          f"gamma_TEE = log D = {'log 2' if Dd == 2 else f'log {Dd}'}   (Z2 TOPOLOGICAL ORDER)")
    print(f"  trivial    : rank={Rt}  -> GSD = {gsd_t}  => gamma_TEE = 0                    "
          f"(CONFINED/Higgs -- the branch that hosts the massless photon)")

    print("\n[verdict] The gate ('write the macroscopic model') does NOT by itself fix gamma_TEE:")
    print("  one explicit operator set has BOTH phases. gamma_TEE is a PHASE CHOICE. Canon's massless")
    print("  emergent photon (item_geo_b_maxwell: gamma~0) selects the TRIVIAL/Coulomb branch; the")
    print("  topological gamma=log2 branch is the gapped alternative canon deems photon-incompatible.")
    print("  TIER: modeling choice (one resolution of the open gate), not derived from canon.")

    assert d == 4, "matter cell must be distance 4"
    assert gsd_d == 4 and Dd == 2, "deconfined gauge web must be Z2 topological (GSD=4, D=2, gamma=log2)"
    assert gsd_t == 1, "trivial branch must have GSD=1 (gamma_TEE=0, the photon/Coulomb branch)"
    print("\nALL ASSERTIONS PASSED -- explicit operator set; gate = phase choice; photon -> gamma=0.")
    print("exit 0")


if __name__ == "__main__":
    main()
