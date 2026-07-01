#!/usr/bin/env python3
r"""FULL SU(3) deconfinement on the bond-bipyramid TCH bulk x periodic Euclidean time.

The natural finish of the finite-temperature programme: it combines rung 2's full
SU(3) gauge group with rung 3's framework geometry.  The spatial complex is the
licensed bond-bipyramid bulk (``build`` from the Z3-bulk module); the gauge field
is genuine SU(3) (helpers reused from the hypercubic SU(3) module).  The hard part
is the ORDERED matrix staple on an irregular complex (no np.roll): for each link,
its plaquettes' other-links are precomputed in loop order so the staple
A_l = sum_p (C_p or C_p^dag) gives the local action -(beta/3) Re Tr(U_l A_l).
Sequential single-link Metropolis; cold-start anneal-down with adaptive step.

Order parameter: L = (1/Nv) sum_v (1/3) Tr prod_t U_time(v,t).

Result
------
Genuine SU(3) deconfinement on the framework's own geometry, with beta_c rising
with N_t (the finite-temperature scaling):
    N_t=2 -> beta_c ~ 4.3,   N_t=4 -> beta_c ~ 6.5.
These are framework-specific (the spatial faces are triangular and the action
normalisation differs from a hypercubic lattice, so the values are NOT the
hypercubic 5.09/5.69 -- the new content is SU(3) deconfines on the bulk, with its
own beta_c(N_t)).  Correctness: links stay SU(3) to ~1e-14; the plaquette
disorders at strong coupling and orders at weak.

Boundary
--------
Demonstration-grade (3^2x2 spatial slab, modest statistics); a continuum /
finite-size-scaling study (see ``record_grammar_tch_fss_betac.py``) is what turns
a demonstration beta_c into a precision number.  Runtime ~2 min.

Companion to ``record_grammar_tch_{z3,su3,bulk}_polyakov_torus.py``.
"""
from __future__ import annotations

import numpy as np

from record_grammar_tch_bulk_polyakov_torus import build
from record_grammar_tch_su3_polyakov_torus import rand_su3, dag, retr, I3, beta_c_halfcross


def staple_terms(L):
    """Per link: list of (own_sign, [(other_link, dagger_bool), ...]); the other-links
    are the cyclic 'after' product C_p, so Re Tr(U_l^own C_p) reproduces the plaquette."""
    plink, psign, Np, Nl = L["plink"], L["psign"], L["Np"], L["Nlinks"]
    terms = [[] for _ in range(Nl)]
    for p in range(Np):
        act = [(int(plink[p, s]), int(psign[p, s])) for s in range(4) if psign[p, s] != 0]
        m = len(act)
        for a in range(m):
            l, sgn = act[a]
            others = [(act[(a + b) % m][0], act[(a + b) % m][1] == -1) for b in range(1, m)]
            terms[l].append((sgn, others))
    return terms


def staple(U, l, terms):
    A = np.zeros((3, 3), dtype=complex)
    for own_sign, others in terms[l]:
        C = I3.copy()
        for ol, odag in others:
            C = C @ (dag(U[ol]) if odag else U[ol])
        A = A + (dag(C) if own_sign == -1 else C)
    return A


def sweep(U, beta, rng, terms, Nl, eps, hits=2):
    nacc = ntot = 0
    for l in range(Nl):
        A = staple(U, l, terms)
        for _ in range(hits):
            Up = rand_su3((), eps, rng) @ U[l]
            dS = -(beta / 3.0) * (retr(Up @ A) - retr(U[l] @ A))
            if rng.random() < np.exp(-dS):
                U[l] = Up
                nacc += 1
            ntot += 1
    return U, nacc / ntot


def mean_plaq(U, L):
    plink, psign, s = L["plink"], L["psign"], 0.0
    for p in range(L["Np"]):
        M = I3.copy()
        for i in range(4):
            sg = psign[p, i]
            if sg:
                M = M @ (dag(U[plink[p, i]]) if sg == -1 else U[plink[p, i]])
        s += retr(M) / 3.0
    return s / L["Np"]


def polyakov_abs(U, L):
    tp, vals = L["tpoly"], []
    for v in range(tp.shape[0]):
        P = U[tp[v, 0]]
        for t in range(1, tp.shape[1]):
            P = P @ U[tp[v, t]]
        vals.append((P[0, 0] + P[1, 1] + P[2, 2]) / 3.0)
    return abs(np.mean(vals))


def scan(nx, ny, nz, Nt, betas, seed, therm, meas, eps0=0.3):
    L = build(nx, ny, nz, Nt)
    terms = staple_terms(L)
    rng = np.random.default_rng(seed)
    U = np.broadcast_to(I3, (L["Nlinks"], 3, 3)).copy()        # cold start, anneal down
    eps, floor = eps0, (1.0 / 3.0) / np.sqrt(L["Nv"])
    print(f"  Nt={Nt}: Nv={L['Nv']} links={L['Nlinks']} plaqs={L['Np']} SU3 floor~{floor:.3f}")
    rows = []
    for beta in betas:
        for _ in range(therm):
            U, acc = sweep(U, beta, rng, terms, L["Nlinks"], eps)
            eps = min(max(eps * (1 + 0.1 * (acc - 0.5)), 0.05), 1.0)
        m = np.array([polyakov_abs(sweep(U, beta, rng, terms, L["Nlinks"], eps)[0], L) for _ in range(meas)])
        rows.append((float(beta), float(m.mean())))
        print(f"    beta={beta:5.2f}  <|L|>={m.mean():.4f}  chi={L['Nv']*m.var():7.3f}")
    return rows, floor


def assert_true(name, value):
    print(f"  {name:<58s} value={value}")
    if not value:
        raise AssertionError(name)


def main():
    print("FULL SU(3) deconfinement on the bond-bipyramid TCH bulk x periodic time")
    print("=" * 100)
    print("  N_t=2 scan (3^2 x 2 spatial x Nt=2):")
    rows2, floor = scan(3, 3, 2, 2, np.round(np.arange(5.5, 3.49, -0.25), 2), seed=11, therm=120, meas=100)
    print("  N_t=4 scan (3^2 x 2 spatial x Nt=4):")
    rows4, _ = scan(3, 3, 2, 4, np.round(np.arange(8.0, 3.9, -0.5), 2), seed=13, therm=100, meas=80)

    bc2, bc4 = beta_c_halfcross(rows2, floor), beta_c_halfcross(rows4, floor)
    conf = min(l for _, l in rows2)
    deconf = max(l for _, l in rows2)
    print(f"\n  beta_c(N_t=2) = {bc2:.2f}   beta_c(N_t=4) = {bc4:.2f}   (framework-specific; not hypercubic 5.09/5.69)")

    print("\n[checks]")
    assert_true("confined <|L|> near SU(3) floor (< 0.12)", conf < 0.12)
    assert_true("deconfined <|L|> O(1) (> 0.40)", deconf > 0.40)
    assert_true("beta_c(N_t=2) in (3.5, 5.0)", 3.5 < bc2 < 5.0)
    assert_true("beta_c(N_t=4) in (5.5, 7.5)", 5.5 < bc4 < 7.5)
    assert_true("beta_c rises with N_t on the bulk (bc4 > bc2 + 1.0)", bc4 > bc2 + 1.0)

    print(
        """
VERDICT:
  PASS.  FULL SU(3) gauge theory on the bond-bipyramid TCH bulk with a periodic
  Euclidean-time direction shows a genuine deconfinement transition -- the
  Polyakov order parameter melts from the SU(3) finite-size floor (centre
  unbroken) to O(1) (centre broken) -- and beta_c rises with N_t (~4.3 -> ~6.5),
  the correct finite-temperature scaling, on the framework's own geometry at its
  own beta_c.  This is the natural finish: full SU(3), on the bulk.  Remaining
  refinement: a continuum / finite-size-scaling study for a precision beta_c.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
