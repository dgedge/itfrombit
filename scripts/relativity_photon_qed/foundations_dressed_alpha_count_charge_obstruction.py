#!/usr/bin/env python3
r"""foundations_dressed_alpha_count_charge_obstruction.py

Dressed-alpha (137.035999) "push further" attempt — does the G7 resolution transfer? Result: NO, and there
is a PRINCIPLED reason, which downgrades dressed-alpha from "open, keep trying" to "bounded by a principled
obstruction." Self-asserting; exit 0. (Pure arithmetic; runs anywhere, kept on ~/bin/py13_7 for the venv.)

The state (dressed_alpha_*; ANCHOR l.612): bare 1/alpha_0=137=T(16)+1 is DERIVED (a flat service COUNT via
the monitored equipartition, Evans-Frigerio, item 79/R14). The dressing delta = 1/alpha - 137 = 0.035999
needs, written one-loop-style delta = N1 * alpha/(2pi), a kernel N1 ~= 31. Canon's decisive negatives:
N1=31 is a COUNT (2*16-1 = dim SO(10)-spinor - 1), the photon self-energy is CHARGE-WEIGHTED
(sum Q_f^2 N_c = 8 Dirac / ~12 web-dressed); "no loop order takes 8 to 31"; the sector-billing map and a
pinned Ward/Kubo normalisation both fail (theorem-shaped negatives).

THE G7 PARALLEL (and why it does NOT rescue dressed-alpha):
  In G7 the same clash appeared -- the framework's native billing is a COUNT/monitor (unital), while the
  standard observable is CHARGE-WEIGHTED/golden-rule. There it RESOLVED: golden-rule was a CATEGORY ERROR
  (the coherent walk is code-preserving, makes no errors), so the count/monitor billing is the whole story
  and C_loop=3/2 firmed up. The reason that worked: the G7 observable (M_P) is itself UN-CLOSABLE (rank-1 +
  T7), so there is no external charge-weighted observable to contradict the count billing.

  Dressed-alpha is the MIRROR case. Here the physical observable 1/alpha = 137.035999 is MEASURED (electron
  g-2, atomic recoil) -- it IS a charge-weighted QED vertex/self-energy quantity. So the charge-weighting is
  NOT a category error here; it is the genuine observable. The count near-match (N1=31) therefore CANNOT be
  promoted to a derivation: it would have to OVERTURN the charge-weighting of a measured QED observable,
  which no framework principle supplies. The very feature that let G7 resolve (no external charge-weighted
  observable) is ABSENT here (1/alpha is exactly that observable).

So the obstruction is principled, not a missing trick: dressed-alpha is the ONE place the framework's
count-billing meets a directly-measured charge-weighted observable head-on. The 31-count gives a striking
~tens-of-ppb numerical near-match, but it is coincidence-grade by construction (wrong operator class), and
the framework's own flat-count proxy (the monitor occupation) gives 0.9956*31 ~= 30.86, not 31 -- a residual
that is itself the count-vs-charge gap. Bare 1/alpha_0=137 (the flat count) stays the clean win.
"""
import math


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    A0_INV, A_INV = 137.0, 137.035999084
    A0, A = 1 / A0_INV, 1 / A_INV
    delta = A_INV - A0_INV
    TPI = 2 * math.pi

    print("=== the dressing as a one-loop-style kernel N1 ===")
    N1_bare = delta * TPI / A0
    N1_sc = delta * TPI / A          # self-consistent (dressed alpha in the loop)
    print(f"    delta = 1/alpha - 137 = {delta:.6f}")
    ok(abs(N1_bare - 31) < 0.1, f"N1 = delta*2pi/alpha_0 = {N1_bare:.3f} ~= 31 (the kernel the shift demands)")
    ok(abs(N1_sc - 31) < 0.1, f"N1(self-consistent) = {N1_sc:.3f} ~= 31")

    print("\n=== two candidate kernels: a COUNT vs a CHARGE-WEIGHTING ===")
    count_kernel = 2 * 16 - 1                                   # 2*16-1 = dim(SO(10) Dirac spinor) - 1
    charge_weight = 3 * (1 * 1.0**2 + 3 * (2/3)**2 + 3 * (1/3)**2)   # sum Q_f^2 N_c over 3 SM generations
    ok(count_kernel == 31, f"COUNT kernel 2*16-1 = {count_kernel} (Dirac-doubled matter modes minus the photon)")
    ok(abs(charge_weight - 8) < 1e-9, f"CHARGE-WEIGHT kernel sum Q_f^2 N_c (3 gen) = {charge_weight:.0f} (web-dressed ~12); 'no loop order takes 8->31'")

    def inv_alpha(N):
        return A0_INV + N * A / TPI                              # self-consistent one-loop-style
    print("\n=== predictions ===")
    ia_count = inv_alpha(31); ia_cw = inv_alpha(charge_weight)
    ppb_count = abs(ia_count / A_INV - 1) * 1e9
    print(f"    COUNT (N=31):          1/alpha = {ia_count:.6f}  (obs {A_INV:.6f}; {ppb_count:.0f} ppb) -- striking, but wrong operator class")
    print(f"    CHARGE-WEIGHT (N=8):   1/alpha = {ia_cw:.6f}  -> shift {100*(ia_cw-A0_INV)/delta:.0f}% of observed (undershoots; ~12 web-dressed ~ 39%)")
    ok(ppb_count < 200, f"the N1=31 COUNT lands the physical value to {ppb_count:.0f} ppb -- a real near-match (canon: a mode-count FIT, not a derivation)")
    ok((ia_cw - A0_INV) / delta < 0.4, "the genuine CHARGE-WEIGHTED kernel (8, ~12 dressed) UNDERSHOOTS the shift (<40%)")

    print("\n=== the framework's own flat-count proxy (monitor occupation) ===")
    monitor = 0.9956 * 31                                        # the monitor near-hit, l.612
    print(f"    monitor occupation = 0.9956 * 31 = {monitor:.2f} (NOT 31) -> 1/alpha = {inv_alpha(monitor):.6f}")
    ok(abs(monitor - 31) > 0.1, "even the framework's flat-count proxy gives ~30.86, not 31: the 0.44% gap IS the count-vs-charge residual")

    print("\n[verdict] dressed-alpha 'push':")
    print("  - The G7 resolution does NOT transfer. In G7 the charge-weighted (golden-rule) object was a")
    print("    CATEGORY ERROR because the observable (M_P) is itself un-closable (T7) -- so the count/monitor")
    print("    billing was the whole story. Dressed-alpha is the MIRROR: 1/alpha=137.036 is a MEASURED,")
    print("    genuinely charge-weighted QED observable (g-2). So the charge-weighting here is REAL, not a")
    print("    category error, and the N1=31 COUNT near-match (~tens of ppb) CANNOT be promoted -- it would")
    print("    have to overturn the charge-weighting of a measured observable, which no principle supplies.")
    print("  - PRINCIPLED OBSTRUCTION (the push's actual result): dressed-alpha is the one place the framework's")
    print("    count-billing meets a directly-measured charge-weighted observable head-on. That is WHY it is")
    print("    stuck -- not a missing trick. It is bounded, coincidence-grade on the count side, undershooting")
    print("    on the (correct) charge-weighted side; bare 1/alpha_0=137 (the flat count) stays the clean win.")
    print("  - So: NOT pushable to closure; the contribution is the principled diagnosis (count-vs-charge clash,")
    print("    the exact mirror of G7) that downgrades it from 'open, keep trying' to 'bounded by obstruction'.")
    print("  exit 0")


if __name__ == "__main__":
    main()
