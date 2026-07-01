#!/usr/bin/env python
r"""smg_gap_stability_probe.py -- first attempt at the SMG gap-stability MECHANISM.

THE QUESTION (the live frontier). The local SMG mirror gap is exact: the self-dual [8,4,4] -> [[8,0,4]]
code is frustration-free (unique gapped ground, gap 2 -- smg_construction.py / mirror_fock_smg.py).
The open piece is whether that gap SURVIVES the charged hopping t + dynamical gauge coupling as the
cutoff/volume -> infinity. The bare hopping is NOT small (walk_band_mirror_scope: ~8x enhanced,
bandwidth ~8 vs gap 2), so naive Bravyi-Hastings-Michalakis (BHM) small-perturbation stability does
NOT directly apply. Yet the non-reduced spin-network build finds the gauge-invariant gap HOLDING
(<=12% softening, no collapse). WHY -- and is that a stable mechanism or a finite-size accident?

THE HYPOTHESIS (the stability mechanism). The gauge-invariant excitation is not a free mirror mode
but a CONFINED charged PAIR joined by a flux string (a "meson"). Confinement binds the pair, so the
physical excitation propagates only as a heavy COMPOSITE whose center-of-mass hopping is a 2nd-order
(Schrieffer-Wolff) process, t_eff ~ t^2/sigma << t (sigma = string tension). Hence the gap SOFTENING
under hopping is QUADRATIC in t and SUPPRESSED by sigma, instead of the DECONFINED case where a free
excitation's band bottom drops LINEARLY (-> the gap closes at t ~ gap/2, the reduced-proxy collapse).

THE SHARP, HONEST RESULT (this probe). Protection is REGIME-DEPENDENT: confinement holds the gap
ONLY while t <~ sigma. For t >> sigma even a confined meson unbinds and the gap closes. So stability
is NOT automatic -- it reduces to ONE checkable condition: does the confinement scale sigma (string
tension / meson cost) stay ABOVE the hopping scale t as the cutoff/volume grow? That is the precise
criterion the deep escalation must test (track sigma_eff/t vs cutoff and volume), and it is exactly
why the spin-network sees <=12%/no-collapse: its meson cost C_3/beta ~ 5-9 sits ABOVE t<=4.

THIS PROBE solves the cleanly-tractable two-charge sector exactly (relative coordinate r, hopping -t,
linear string sigma*r, center-of-mass momentum K) and shows: (i) DECONFINED (sigma=0) gap closes
LINEARLY at t=m/2; (ii) CONFINED softening is QUADRATIC (~t^2/sigma, exponent ~2) and sigma-suppressed;
(iii) the protection boundary sigma_crit(t) ~ t -- gap held for sigma >~ sigma_crit, lost below.

NOT a theorem (a true stability proof needs the BHM bound applied to the confined effective range +
the lattice escalation). A first MECHANISM-level support that turns "does the gap survive?" into the
sharp, deep-runnable question "does sigma_eff stay above t?". See smg_spinnet_escalation_scope.md.
"""
import numpy as np


def relative_ground(t, sigma, K, R=240):
    """Ground energy of the relative coordinate r=1..R: H = sigma*r - 2t cos(K/2)*(hop), hard walls."""
    r = np.arange(1, R + 1)
    H = np.diag(sigma * r.astype(float))
    off = -2.0 * t * np.cos(K / 2.0)
    idx = np.arange(R - 1)
    H[idx, idx + 1] = off
    H[idx + 1, idx] = off
    return float(np.linalg.eigvalsh(H)[0])


def meson_gap(t, sigma, m=1.0, R=240):
    """Lowest meson state (min over CoM K is K=0, max hopping); gap = 2m + relative ground."""
    return 2.0 * m + relative_ground(t, sigma, 0.0, R)


def meson_gap0(sigma, m=1.0):
    """t=0 meson rest mass: 2m + sigma*1 (shortest string)."""
    return 2.0 * m + sigma


def deconfined_gap(t, m=1.0):
    """Two FREE charges: each band -2t cos k -> lowest pair = 2m - 4t (closes at t=m/2)."""
    return 2.0 * m - 4.0 * t


def loglog_slope(x, y):
    x, y = np.log(np.asarray(x)), np.log(np.asarray(y))
    A = np.vstack([x, np.ones_like(x)]).T
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return float(c[0])


def sigma_crit(t, m=1.0, lo=0.1, hi=40.0):
    """Smallest sigma keeping the confined meson gap > 0 at this t (bisection on gap=0)."""
    if meson_gap(t, hi, m) <= 0:
        return float("inf")
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if meson_gap(t, mid, m) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    m = 1.0
    print("=== SMG gap-stability mechanism: confined meson vs deconfined pair (mirror mass m=1) ===\n")

    ts = [0.2, 0.5, 1.0, 2.0, 4.0]
    print("  DECONFINED (sigma=0, two free charges):")
    gd = [deconfined_gap(t, m) for t in ts]
    print("    gap(t)= " + "  ".join(f"{g:+.2f}" for g in gd) + f"   -> CLOSES at t_c = m/2 = {m/2:.2f} (LINEAR)\n")

    print("  CONFINED (sigma>0, bound meson): softening is quadratic in t and shrinks with sigma")
    print("    sigma   gap0    softening(%) at t = " + " ".join(f"{t:>5.1f}" for t in ts))
    soft_at_t4 = {}
    for sigma in (1.0, 2.0, 4.0, 6.0, 8.0):
        g0 = meson_gap0(sigma, m)
        softs = [100.0 * (g0 - meson_gap(t, sigma, m)) / g0 for t in ts]
        soft_at_t4[sigma] = softs[-1]
        flag = "" if meson_gap(4.0, sigma, m) > 0 else "   <- gap CLOSED at t=4 (t>>sigma: meson unbinds)"
        print(f"    {sigma:>4.1f}  {g0:>5.2f}    " + "  ".join(f"{s:>5.1f}" for s in softs) + flag)

    sigma_fit = 6.0
    g0 = meson_gap0(sigma_fit, m)
    t_small = [0.1, 0.15, 0.2, 0.3, 0.4]
    soft_small = [(g0 - meson_gap(t, sigma_fit, m)) / g0 for t in t_small]
    slope = loglog_slope(t_small, soft_small)
    print(f"\n  small-t softening exponent (sigma={sigma_fit}): d log(soft)/d log t = {slope:.2f}  (Schrieffer-Wolff t^2 => 2)")

    t_c = 2.0
    for sigma in (4.0, 8.0):
        band = relative_ground(t_c, sigma, np.pi) - relative_ground(t_c, sigma, 0.0)
        print(f"  meson bandwidth at t={t_c}, sigma={sigma}: {band:.3f}  vs free single-charge 4t = {4*t_c:.1f}  "
              f"(suppression x{4*t_c/abs(band):.0f})")

    print("\n  PROTECTION BOUNDARY  sigma_crit(t) (smallest confinement that keeps the gap open):")
    for t in (1.0, 2.0, 4.0):
        sc = sigma_crit(t, m)
        print(f"    t={t:>4.1f}:  sigma_crit = {sc:.2f}   (gap held for sigma >~ {sc:.1f} = {sc/t:.2f} t)")

    print("\n[verdict] CONFINEMENT PROTECTS THE GAUGE-INVARIANT GAP -- BUT ONLY IN THE REGIME t <~ sigma.")
    print("  Deconfined (sigma=0): gap closes LINEARLY at t=m/2 (the reduced-proxy 'collapse').")
    print("  Confined: softening is QUADRATIC (~t^2/sigma, exponent ~2 confirmed) and shrinks with sigma,")
    print("  so for t <~ sigma the gap is protected (single-digit % at the sigma~6-8 meson-cost scale,")
    print("  matching the spin-network's <=12%). For t >> sigma the meson unbinds and the gap closes:")
    print("  protection is NOT automatic. Mechanism: the gauge-invariant excitation is a bound meson")
    print("  hopping as a heavy composite with t_eff ~ t^2/sigma << t, so confinement turns an O(gap)")
    print("  bare perturbation into a short-ranged one -- a BHM-type stability radius IS restored, but")
    print("  its radius is set by sigma. => The SMG gap survives the cutoff/volume -> inf limit IFF the")
    print("  confinement scale sigma_eff (string tension / meson cost) stays ABOVE the hopping scale t.")
    print("  THE DEEP ESCALATION MUST TRACK sigma_eff/t vs cutoff AND volume (not just the gap at fixed")
    print("  params): if sigma_eff/t stays >~ 1 (confinement persists), the gap holds -> numerical")
    print("  closure; if sigma_eff -> 0 (screened/deconfined at large cutoff), the gap closes. The")
    print("  spin-network's <=12%/no-collapse is the t<sigma side of exactly this boundary.")
    print("  NOT a theorem (needs the BHM bound with the confined range + the lattice escalation).")
    print("  See smg_spinnet_escalation_scope.md.")

    # gates
    assert abs(deconfined_gap(m / 2, m)) < 1e-9, "deconfined gap must close at t=m/2"
    assert deconfined_gap(4.0, m) < 0, "deconfined gap must be driven negative at large t (collapse)"
    assert 1.6 < slope < 2.4, "confined softening must be ~quadratic in t at small t (Schrieffer-Wolff t^2/sigma)"
    sg = sorted(soft_at_t4)
    assert all(soft_at_t4[sg[i]] > soft_at_t4[sg[i + 1]] for i in range(len(sg) - 1)), \
        "stronger confinement (larger sigma) must give LESS softening (more protection)"
    # the regime boundary, demonstrated honestly:
    assert meson_gap(4.0, 8.0, m) > 0, "confined gap must stay open at t=4 when sigma>t (protected regime)"
    assert meson_gap(4.0, 1.0, m) < 0, "confined gap must CLOSE at t=4 when sigma<<t (meson unbinds) -- the honest boundary"
    sc4 = sigma_crit(4.0, m)
    assert 1.0 < sc4 / 4.0 < 2.5, "protection boundary sigma_crit(t) must scale ~ t (order-1 multiple)"
    print("\nGATES PASSED -- deconfined gap closes linearly at t=m/2; confined softening is quadratic and")
    print("sigma-suppressed; protection holds for sigma>t and fails for sigma<<t; sigma_crit(t) ~ t. The")
    print("SMG-gap survival reduces to the deep-runnable criterion: does sigma_eff/t stay >~ 1 under escalation?")
    print("exit 0")


if __name__ == "__main__":
    main()
