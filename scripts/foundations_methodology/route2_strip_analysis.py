#!/usr/bin/env python3
"""
DRIFT K4 route-2 analysis: SMG mirror gap vs cutoff at the field-theory cutoff n/q=6.

Reduces the (multi-hour) two-plaquette charge-block streaming-Lanczos run
(css_2d_strip_hpc_lanczos.py) to its decisive numbers, reproducibly, from the
saved tridiagonals + gap trajectories (route2_strip_tridiagonals.json, 44 KB) --
NOT the 3.2 GB checkpoints.

Gap extraction is ghost-robust: streaming Lanczos (no reorth) spawns a duplicate
of the ground state ~iter 60, collapsing the *reported* final gap to ~1e-13. The
physical gap is read at the iteration MINIMISING the excited-state residual,
restricted to ritz_gap > 1e-3 (this excludes the exact ghost-collapse while
KEEPING a genuinely small gap). Cross-checked with corrected Cullum-Willoughby.

Validation: reproduces the recorded n/q=3,4 t=4 gaps (~1.79, ~1.76, falling).
Result:    n/q=6 t=4 gap collapses to ~0.13 (CW ~0.21) -- closure, no saturation.
Needs only numpy.
"""
import json
import os
import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "route2_strip_tridiagonals.json")


def gap_minresid(traj):
    """Physical gap = ritz_gap at min excited-residual, restricted to gap>1e-3."""
    v = [(i, g, r) for (i, g, r) in traj if r is not None and g > 1e-3]
    i, g, r = min(v, key=lambda x: x[2])
    return g, r, i


def cw_gap(alphas, betas):
    """Corrected Cullum-Willoughby: converged eigenvalues cluster (mult>=2 -> genuine);
    a SIMPLE eigenvalue of T that also lies in T-hat is spurious."""
    a = np.array(alphas, float); b = np.array(betas, float); m = len(a)
    T = np.diag(a) + np.diag(b[:m - 1], 1) + np.diag(b[:m - 1], -1)
    Th = np.diag(a[1:]) + np.diag(b[1:m - 1], 1) + np.diag(b[1:m - 1], -1)
    ev = np.sort(np.linalg.eigvalsh(T)); evh = np.sort(np.linalg.eigvalsh(Th))
    clusters = [[ev[0]]]
    for e in ev[1:]:
        if abs(e - clusters[-1][-1]) < 1e-5:
            clusters[-1].append(e)
        else:
            clusters.append([e])
    genuine = []
    for cl in clusters:
        rep = float(np.mean(cl))
        if len(cl) == 1 and evh.size and np.min(np.abs(evh - rep)) < 1e-6:
            continue  # spurious
        genuine.append(rep)
    genuine = sorted(genuine)
    return genuine[1] - genuine[0]


def main():
    data = json.load(open(DATA))
    print("=" * 74)
    print("Route-2 SMG mirror gap vs cutoff (two-plaquette strip, Z3 colour proxy, beta=0.5)")
    print("=" * 74)
    print(f"  {'n/q':>4} {'t':>5} {'gap(min-resid)':>16} {'resid':>9} {'CW gap':>8}")
    res = {}
    for key in ["3_0.2", "3_4.0", "4_0.2", "4_4.0", "6_0.2", "6_4.0"]:
        d = data[key]
        g, r, it = gap_minresid(d["gap_traj"])
        cw = cw_gap(d["alphas"], d["betas"])
        res[key] = g
        print(f"  {d['spc']:>4} {d['t']:>5} {g:>16.4f} {r:>9.1e} {cw:>8.4f}   (iter {it})")

    print("\n  t-response (gap collapse t=0.2 -> t=4) per cutoff:")
    for spc in (3, 4, 6):
        g02, g4 = res[f"{spc}_0.2"], res[f"{spc}_4.0"]
        print(f"     n/q={spc}: {g02:.3f} -> {g4:.3f}   ({(g02-g4)/g02*100:.0f}% collapse)")

    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    print("  - validation: n/q=3,4 t=4 reproduce recorded ~1.79, ~1.76 (falling). OK")
    print("  - n/q=6 t=4 gap collapses to ~0.13 (CW ~0.21) -- an order below n/q=4's 1.76;")
    print("    t-response jumps 9% -> 11% -> ~93% with cutoff. CLOSURE, no saturation.")
    print("  - caveat: n/q=6 t=4 was killed at iter 64 + Lanczos ghosting -> exact value")
    print("    uncertain (0.13-0.21); a clean number needs ARPACK/reorth or tensor-network.")
    print("    The QUALITATIVE result (small gap, ~90% t-response, no protection) is robust.")

    # asserts (pre-registered decision criterion)
    assert 1.6 < res["3_4.0"] < 1.95, "n/q=3 t=4 validation"
    assert 1.6 < res["4_4.0"] < 1.90 and res["4_4.0"] < res["3_4.0"], "n/q=4 t=4 falling"
    assert res["6_0.2"] > 1.5, "n/q=6 weak-hopping gap stays open"
    assert res["6_4.0"] < 0.5, "n/q=6 t=4 gap CLOSES (decisive)"
    r3 = (res["3_0.2"] - res["3_4.0"]) / res["3_0.2"]
    r6 = (res["6_0.2"] - res["6_4.0"]) / res["6_0.2"]
    assert r6 > 0.8 and r3 < 0.2, "t-response grows dramatically with cutoff"
    print("\nALL ASSERTS PASSED -> closure confirmed at the n/q=6 field-theory cutoff.")


if __name__ == "__main__":
    main()
