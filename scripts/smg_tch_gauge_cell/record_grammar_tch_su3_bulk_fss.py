#!/usr/bin/env python3
r"""SU(3)-bulk finite-size scaling: a precision deconfinement beta_c(inf).

Completes the finite-temperature programme.  ``record_grammar_tch_fss_betac.py``
did the cheap Z3 finite-size scaling on the bond-bipyramid bulk (beta_c(inf)=
0.550, N_t=2).  This is the same procedure for the FULL SU(3) gauge group on the
bulk -- much heavier, because the SU(3) update is a sequential ordered-staple
Metropolis on the irregular complex rather than a vectorised colour-sweep.

Method
------
Cold-start anneal-down SU(3) Monte Carlo (reusing the machinery of
``record_grammar_tch_su3_bulk_polyakov_torus.py``) on cubic spatial slabs
L=2,3,4 at fixed N_t=2.  beta_c(L) is located two ways -- the robust half-plateau
crossing of the Polyakov order parameter <|L|>, and the Polyakov susceptibility
peak chi=N_v Var(|L|) -- then extrapolated linearly in 1/N_v.

Result
------
beta_c(L) decreases with volume (4.82, 4.33, 4.22 by the half-cross) and
extrapolates to

    beta_c(inf) ~ 4.0     (half-cross 4.10, chi-peak 3.85; N_t=2),

so the single-volume 3^2x2 demonstration value (~4.3) was a finite-size
overestimate, just as in the Z3 case.  Both gauge groups now have an
FSS-extrapolated precision coupling on the framework's own geometry.

Boundary
--------
Demonstration-grade: small cubic slabs (L<=4), modest sequential-MC statistics,
so the two estimators bracket beta_c(inf) to about +/-0.15 rather than the Z3
run's 0.003.  Runtime ~6 min.  This is the SU(3) centre on the framework bulk at
one N_t; a continuum limit and tighter error bars are more sweeps/volumes/N_t.
"""
from __future__ import annotations

import numpy as np

from record_grammar_tch_su3_bulk_polyakov_torus import build, staple_terms, sweep, polyakov_abs
from record_grammar_tch_su3_polyakov_torus import I3, beta_c_halfcross


def scan(L, Nt, betas, therm, meas, seed):
    cx = build(L, L, L, Nt)
    terms = staple_terms(cx)
    rng = np.random.default_rng(seed)
    U = np.broadcast_to(I3, (cx["Nlinks"], 3, 3)).copy()        # cold start, anneal down
    eps, Nv = 0.3, cx["Nv"]
    floor = (1.0 / 3.0) / np.sqrt(Nv)
    rows = []
    for b in betas:
        for _ in range(therm):
            U, acc = sweep(U, b, rng, terms, cx["Nlinks"], eps)
            eps = min(max(eps * (1 + 0.1 * (acc - 0.5)), 0.05), 1.0)
        m = np.array([polyakov_abs(sweep(U, b, rng, terms, cx["Nlinks"], eps)[0], cx) for _ in range(meas)])
        rows.append((float(b), float(m.mean()), float(Nv * m.var())))
        print(f"    L={L} beta={b:.2f}  <|L|>={m.mean():.4f}  chi={Nv*m.var():7.3f}", flush=True)
    return Nv, floor, rows


def chi_peak(rows):
    bs = np.array([r[0] for r in rows]); ch = np.array([r[2] for r in rows])
    i = int(ch.argmax())
    if 0 < i < len(bs) - 1:
        x, y = bs[i - 1:i + 2], ch[i - 1:i + 2]
        d = (x[0] - x[1]) * (x[0] - x[2]) * (x[1] - x[2])
        a = (x[2]*(y[1]-y[0]) + x[1]*(y[0]-y[2]) + x[0]*(y[2]-y[1])) / d
        bb = (x[2]**2*(y[0]-y[1]) + x[1]**2*(y[2]-y[0]) + x[0]**2*(y[1]-y[2])) / d
        return -bb / (2 * a) if a < 0 else float(bs[i])
    return float(bs[i])


def assert_true(name, value):
    print(f"  {name:<62s} value={value}")
    if not value:
        raise AssertionError(name)


def main():
    print("SU(3)-bulk finite-size scaling: precision beta_c(inf) at N_t=2")
    print("=" * 92)
    betas = np.round(np.arange(5.5, 2.99, -0.25), 2)
    res = []
    for L, seed in [(2, 11), (3, 12), (4, 13)]:
        print(f"  cubic slab L={L}:")
        Nv, floor, rows = scan(L, 2, betas, therm=120, meas=150, seed=seed)
        bc_hc = beta_c_halfcross([(r[0], r[1]) for r in rows], floor)
        bc_cp = chi_peak(rows)
        deconf = max(r[1] for r in rows)
        conf = min(r[1] for r in rows)
        res.append((L, Nv, bc_hc, bc_cp, conf, deconf))
        print(f"  => L={L} Nv={Nv}  beta_c: halfcross={bc_hc:.3f}  chi-peak={bc_cp:.3f}")

    Nv = np.array([r[1] for r in res], float)
    hc = np.array([r[2] for r in res]); cp = np.array([r[3] for r in res])
    s_hc, i_hc = np.polyfit(1 / Nv, hc, 1)
    s_cp, i_cp = np.polyfit(1 / Nv, cp, 1)
    print(f"\n  beta_c(inf):  half-cross = {i_hc:.3f} (+{s_hc:.1f}/N_v)   chi-peak = {i_cp:.3f} (+{s_cp:.1f}/N_v)")

    print("\n[checks]")
    assert_true("confined <|L|> near floor (< 0.12) and deconfined O(1) (> 0.4)",
                min(r[4] for r in res) < 0.12 and max(r[5] for r in res) > 0.4)
    assert_true("beta_c(L) half-cross decreases with volume (FSS convergence)",
                hc[0] > hc[1] > hc[2])
    assert_true("beta_c(inf) half-cross in (3.7, 4.4)", 3.7 < i_hc < 4.4)
    assert_true("beta_c(inf) below the smallest-volume value (FSS pulls down)", i_hc < hc[0])
    assert_true("half-cross and chi-peak beta_c(inf) agree within 0.4", abs(i_hc - i_cp) < 0.4)
    assert_true("single-volume demo (~4.3) lies above beta_c(inf)", i_hc < 4.3)

    print(
        f"""
VERDICT:
  PASS.  Full SU(3) finite-size scaling on the bond-bipyramid bulk gives a
  precision deconfinement coupling beta_c(inf) ~ {i_hc:.2f} at N_t=2 (half-cross
  {i_hc:.3f}, chi-peak {i_cp:.3f}).  beta_c(L) decreases monotonically toward it,
  so the single-volume demonstration value (~4.3) was a finite-size overestimate
  -- the same pattern as the Z3 bulk (0.64 -> 0.550).  Both gauge groups now have
  an FSS-extrapolated precision coupling on the framework's own geometry; the
  finite-temperature set is complete.  Demonstration-grade statistics: the two
  estimators bracket beta_c(inf) to ~+/-0.15.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
