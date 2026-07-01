#!/usr/bin/env python3
r"""foundations_koide_feshbach_pin.py

Pinning the lepton Koide scale eps from the EXPLICIT §5.8 two-channel Feshbach, with the C8 standing-wave
energy as the bare nu_R channel. (Final step of the mu_l = m_N/3 thread.)

§5.8 Feshbach: H_eff = [[0, xi],[xi*, eps]] -- open channel = the charged lepton (diagonal 0), closed
channel = the forbidden nu_R pseudocodeword (diagonal eps, the bare channel energy), coupling xi. Canon:
"the k=0 pole reads m = eps", so the lepton mass scale tracks the bare nu_R channel energy eps. Pinning
mu_l == pinning eps.

TWO routes for the bare nu_R channel energy:
  ROUTE A (bare nu_R = the single-(colourless-)fermion CONSTITUENT/dynamical mass = M_N/3):
    eps_A = M_N/3 = (2 sqrt2/3) Lambda = 313 MeV -> m = eps_A = mu_l (313.84, 0.26%), TUNING-FREE.
    Clean, but takes the bare channel = the constituent mass M_N/3 (the nucleon-as-3-constituents relation).
  ROUTE B (bare nu_R = one C8 transverse mode sqrt2*Lambda, + Feshbach dressing):
    requires xi ~ 1.49 Lambda to dress down to mu_l -> NOT parameter-free (xi must be supplied). RULED OUT
    as a parameter-free pin.

RESULT: the Feshbach m=eps pins mu_l = eps = the constituent mass = M_N/3 (Route A, tuning-free); the
single-mode-plus-dressing alternative (Route B) does not close without an independent xi. The crux is
therefore REDUCED to a single, standard-physics-shared residual: derive 'constituent mass = M_N/3'
intrinsically from the C8 standing wave -- M_N = (lambda_1+lambda_7)Lambda = 2 sqrt2 Lambda comes from the
TWO degenerate transverse modes, while the /3 is 3 constituents/colours, and '2 modes != 3 constituents'
is not yet reconciled inside the C8. So eps is pinned to the constituent mass; the C8-intrinsic /3 is the
last (narrow) step. NOT fully closed -- an honest reduction, not a theorem.

Self-asserting (on the routes' verdicts); exit 0.
"""
import numpy as np
from scipy.optimize import brentq


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def feshbach(eps, xi):
    return np.linalg.eigvalsh(np.array([[0.0, xi], [xi, eps]]))


def main():
    Lam = 332.0; MN = 2 * np.sqrt(2) * Lam
    muL = ((np.sqrt(0.51099895) + np.sqrt(105.6583755) + np.sqrt(1776.86)) / 3) ** 2
    print("=== §5.8 Feshbach H_eff = [[0, xi],[xi*, eps]] ; m = eps (k=0 pole) ===")
    ev0 = feshbach(150.0, 0.0)
    ok(abs(ev0[0]) < 1e-9 and abs(ev0[1] - 150.0) < 1e-9,
       "xi->0 eigenvalues {0 (lepton), eps (nu_R)} -> the lepton scale tracks the bare channel eps")

    # Route A
    epsA = MN / 3
    print(f"\n[Route A] bare nu_R = constituent mass M_N/3 = (2 sqrt2/3) Lambda = {epsA:.2f} MeV")
    print(f"          m = eps_A ; mu_l(meas) = {muL:.2f} ({100*(muL/epsA-1):+.2f}%)")
    ok(abs(muL / epsA - 1) < 0.004, "Route A: m = eps = M_N/3 to <0.4%, TUNING-FREE (assumes bare = constituent mass)")

    # Route B
    epsB = np.sqrt(2) * Lam
    xi = brentq(lambda x: abs(feshbach(epsB, x)[0]) - muL, 1, 3000)
    print(f"\n[Route B] bare nu_R = one C8 mode sqrt2*Lambda = {epsB:.1f} MeV ; tune xi for m = mu_l")
    print(f"          xi needed = {xi:.1f} MeV = {xi/Lam:.2f} Lambda")
    ok(xi / Lam > 1.0, "Route B: needs xi ~ 1.5 Lambda (large, must be supplied) -> NOT a parameter-free pin (ruled out)")

    print("\n[verdict] Feshbach pin of mu_l:")
    print("  - m = eps pins the lepton scale to the bare nu_R channel energy; the tuning-free value is")
    print("    eps = M_N/3 (Route A) = the substrate single-fermion constituent/dynamical mass = 313 MeV (0.26%).")
    print("  - the single-mode+dressing alternative (Route B) needs xi ~ 1.5 Lambda -> not parameter-free.")
    print("  - so mu_l = eps = the constituent mass is the pinned result; the crux REDUCES to the narrow,")
    print("    standard-physics-shared residual 'constituent = M_N/3' from the C8 2-transverse-mode standing")
    print("    wave (2 modes vs 3 constituents not yet reconciled inside the C8). REDUCED, not fully closed. exit 0")


if __name__ == "__main__":
    main()
