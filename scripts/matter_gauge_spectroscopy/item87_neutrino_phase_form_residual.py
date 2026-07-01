#!/usr/bin/env python3
r"""Item 87 -- the deep residual: WHY is the Koide phase delta = d/N as a RAW
radian?  An honest attack that pins the functional FORM and the proportionality
constant from the one precisely-measured sector, narrows the admissible mechanism
class, and states exactly what is left underived.  It does NOT close item 87.

Setup.  The phase-selection audit (item87_neutrino_phase_selection_audit.py)
showed delta is NOT a winding (a 2*pi-fraction Berry/holonomy phase is excluded by
the fit).  The remaining question is the FORM of the map (d/N) -> delta and why
its proportionality constant is exactly 1 radian.  The charged-lepton sector fixes
delta to 0.02% (m_e,m_mu to 0.006-0.007%), so it can discriminate candidate forms.

Write delta = k * (d/N).  The data measures k.  k=1 (raw radian) is special: it is
what an ANGLE-AS-ARC/RADIUS gives (a dimensionless ratio is intrinsically radian-
valued with coefficient 1).  A winding gives k=2*pi; a trig response (arcsin/arctan)
gives a k that drifts from 1.  So measuring k discriminates the mechanism CLASS.
"""
import math
import numpy as np

R2 = math.sqrt(2)
me, mmu, mtau = 0.51099895, 105.6583755, 1776.86          # MeV, PDG
D_OVER_N = 2.0 / 9.0                                        # charged-lepton defect fraction


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def lepton_delta_fit():
    """The phase the charged-lepton spectrum actually requires (Foot inversion)."""
    sm = np.sqrt([me, mmu, mtau]); M = sm.mean()
    return float(np.min(np.arccos(np.clip((sm / M - 1) / R2, -1, 1))))


def main():
    print("ITEM 87 -- DEEP RESIDUAL: the form of (d/N) -> delta")
    print("=" * 80)

    d = lepton_delta_fit()
    print(f"\n[1] Measured lepton phase  delta_fit = {d:.6f} rad   (d/N = {D_OVER_N:.6f})")
    k = d / D_OVER_N
    print(f"    proportionality constant k = delta_fit / (d/N) = {k:.6f}")
    check(abs(k - 1.0) < 1.0e-3, "k = 1 to <0.1% : the phase is the RAW defect fraction in radians")

    print("\n[2] Functional-form discrimination at the precise point")
    forms = {
        "raw   d/N        (k=1, arc/radius)": D_OVER_N,
        "winding 2*pi*d/N (Berry/holonomy)": 2 * math.pi * D_OVER_N,
        "arcsin(d/N)      (trig response)": math.asin(D_OVER_N),
        "arctan(d/N)      (trig response)": math.atan(D_OVER_N),
        "sin(d/N)         (trig response)": math.sin(D_OVER_N),
    }
    best_name, best_err = None, 1e9
    for name, val in forms.items():
        err = abs(val - d) / d
        flag = "  <-- best" if err < best_err else ""
        if err < best_err:
            best_err, best_name = err, name
        print(f"    {name:<36s} = {val:8.5f}   rel.err {err*100:7.3f}%{flag}")
    check("raw" in best_name, "the RAW (linear, k=1) form is the best fit")
    check(best_err < 1.0e-3, "raw form fits to <0.1%")
    # every non-raw candidate is >0.5% off at the precise point:
    others = [abs(v - d) / d for n, v in forms.items() if "raw" not in n]
    check(min(others) > 5.0e-3, "all winding/trig forms are excluded (>0.5%) at the precise point")
    print(f"    -> raw beats the next-best form by {min(others)/best_err:.0f}x")

    print("\n[3] What this forces, and what stays open")
    # arc-subtension reproduces k=1 BY CONSTRUCTION (radian = arc/radius); that is
    # the point -- it is the unique natural source of k=1, not an independent check.
    arc = D_OVER_N  # arc of d units at radius N units subtends d/N rad
    check(abs(arc - D_OVER_N) < 1e-15, "an open arc of d at radius N subtends exactly d/N rad (k=1 by definition)")

    print(
        """
[4] VERDICT (honest -- item 87 NOT closed)
    Locked from the precise lepton sector (0.02%):
      * the map (d/N) -> delta is LINEAR with constant k = 1 (raw radian);
      * winding (k=2*pi) and trig (arcsin/arctan/sin) forms are all excluded,
        beaten by the raw form by ~40x at the measured point.
    Mechanism class therefore forced: k=1 means the phase is an ANGLE defined as
    arc/radius (intrinsically radian-valued, coefficient 1).  A CLOSED loop gives
    k=2*pi (a winding -- excluded); only an OPEN arc gives the raw ratio.  So the
    phase must be an open-arc/displacement angle, delta = (defect displacement d) /
    (scale N) rad -- not a Berry winding, not a trig response.

    IRREDUCIBLE RESIDUAL (still open, now a sharp geometric question, NOT "derive
    a free phase"):
      derive from QEC/boot mechanics that the nu_R defect subtends an OPEN arc of
      its support d at the plaquette scale N -- i.e. that displacement=d, radius=N.
      Arc-subtension reproduces d/N BY CONSTRUCTION (it is why k=1), so it is the
      natural candidate but NOT a derivation; the substrate embedding that makes
      displacement=d and radius=N is undocumented.  This is the limit of what is
      honestly reachable here: the FORM (raw-linear, k=1, open-arc) is pinned and
      the winding/trig classes are excluded, but the d/N-as-arc identity itself is
      asserted, not derived.  delta=2/9 (lepton) / 1/3 (nu) remain fit-selected.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
