#!/usr/bin/env python3
r"""Non-circular route to the dressed-alpha magnitude (+0.036): does the field-content 1-loop give it?

The web/vertex route is circular (it uses t=alpha0 to dress alpha0). But the Part-12 Schwinger-Dyson
form is NOT circular:

    alpha^-1 (alpha^-1 - 137) = N1/(2*pi)        =>   delta = alpha^-1 - 137 ~= N1 * alpha0 / (2*pi),

where alpha0 = 1/137 is the SEPARATELY-DERIVED bare value (R14) and N1 is a pure, alpha-independent
1-loop KERNEL coefficient. So the dressing is (bare coupling) x (pure loop number) -- standard
perturbation theory, not circular. The only non-circular question is: does the kernel N1 come from the
substrate's CHARGED FIELD CONTENT (the photon self-energy = charge-weighted fermion loop), or is N1=31
a fit (the K5 caveat: the count != the genuine integral)?

This script computes the genuine, non-circular field-content kernels and compares to the N1=31 that
reproduces +0.036:
  - Sigma Q_f^2 N_c  (Dirac fermions, the QED vacuum-polarization sum, 3 generations) = 8
  - Sum Q^2          (the SO(10) 16 Weyl count, the framework's anomaly invariant, 3 generations) = 16
The verdict: a non-circular route gives the RIGHT ORDER but UNDERSHOOTS +0.036 by ~2-4x; the N1=31 that
hits it is NOT the charge-weighted coefficient (=> a fit, K5). Circularity is resolved; the kernel gap
is the open object.
"""
from __future__ import annotations

import math

ALPHA0 = 1.0 / 137.0
DELTA_PHYS = 0.036                       # measured dressed shift: 137.035999 - 137


def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c:
        raise AssertionError(m)


def main():
    print("NON-CIRCULAR dressed-alpha magnitude: field-content 1-loop kernel vs the physical +0.036")
    print("=" * 98)

    # ---- the non-circular structure: delta = N1 * alpha0 / (2 pi); the kernel N1 the physical value needs ----
    N1_target = DELTA_PHYS * 2 * math.pi / ALPHA0      # solve delta = N1 alpha0/2pi  (=> N1 = delta*2pi*137)
    print("\n[1] the non-circular structure (Part-12 Schwinger-Dyson):")
    print("    delta = N1 * alpha0 / (2 pi),  alpha0 = 1/137 the DERIVED bare value, N1 an alpha-independent kernel")
    print(f"    the physical delta = {DELTA_PHYS} needs kernel  N1 = delta*2pi/alpha0 = {N1_target:.2f}  (Part-12's N1=31)")
    print("    NB this is NOT circular: bare coupling x pure loop number. (The web/vertex route's t=alpha0 was.)")
    check(abs(N1_target - 31) < 0.5, "the physical dressing needs kernel N1 ~= 31 (Part-12)")

    # ---- the genuine, non-circular kernels from the charged field content ----
    print("\n[2] genuine non-circular kernels from the substrate's CHARGED field content (no alpha used):")
    # QED vacuum-polarization sum over Dirac fermions: Sigma_f Q_f^2 * N_c, 3 generations
    leptons = 3 * 1 * (1.0) ** 2                        # e, mu, tau: Q=-1, N_c=1
    up = 3 * 3 * (2.0 / 3) ** 2                         # u, c, t:   Q=+2/3, N_c=3
    down = 3 * 3 * (1.0 / 3) ** 2                       # d, s, b:   Q=-1/3, N_c=3
    sumQ2_dirac = leptons + up + down                  # = 8
    sumQ2_weyl = 2 * sumQ2_dirac                        # SO(10) 16 = L+R Weyl => 16 (framework anomaly invariant)
    print(f"    Sigma Q_f^2 N_c (Dirac, 3 gen)  = {leptons:.0f}(lep) + {up:.0f}(up) + {down:.0f}(down) = {sumQ2_dirac:.0f}")
    print(f"    Sum Q^2 (Weyl / SO(10) 16, 3 gen) = 2 x {sumQ2_dirac:.0f} = {sumQ2_weyl:.0f}   (= the framework's anomaly Sum Q^2)")
    check(abs(sumQ2_dirac - 8) < 1e-9 and abs(sumQ2_weyl - 16) < 1e-9,
          "field-content charge sums are 8 (Dirac vacuum-pol) and 16 (Weyl/SO(10) 16) — both alpha-free")

    # ---- the delta each genuine kernel predicts, vs the physical +0.036 ----
    print("\n[3] delta from each GENUINE (non-circular) kernel, vs physical +0.036:")
    print(f"    {'kernel':>26} {'value':>6} {'delta = kernel*alpha0/2pi':>26} {'vs 0.036':>10}")
    for name, k in (("Sigma Q_f^2 (Dirac)", sumQ2_dirac), ("Sum Q^2 (Weyl 16)", sumQ2_weyl), ("Part-12 N1 (fit)", 31)):
        dl = k * ALPHA0 / (2 * math.pi)
        print(f"    {name:>26} {k:>6.0f} {dl:>26.4f} {dl/DELTA_PHYS:>9.2f}x")
    d_dirac = sumQ2_dirac * ALPHA0 / (2 * math.pi)
    d_weyl = sumQ2_weyl * ALPHA0 / (2 * math.pi)
    check(d_dirac < DELTA_PHYS and d_weyl < DELTA_PHYS, "BOTH genuine charge-weighted kernels UNDERSHOOT +0.036")
    check(d_dirac > DELTA_PHYS / 5 and d_weyl > DELTA_PHYS / 3, "...but to the RIGHT ORDER (within ~2-4x)")

    # ---- is N1=31 a clean field-content count? (the K5 'count vs integral' question) ----
    print("\n[4] does N1=31 have a non-circular field-content derivation?")
    print(f"    target N1 = {N1_target:.1f}.  Genuine charge-weighted sums: Dirac 8, Weyl 16.")
    print(f"    31 / 8  = {31/8:.2f}  (~4x the Dirac sum);   31 / 16 = {31/16:.2f}  (~2x the Weyl sum)")
    # candidate raw mode counts (NOT charge-weighted -- the wrong physics for a photon self-energy):
    charged_weyl_3gen = 2 * (3 * (1 + 3 + 3))           # per gen: e+3u+3d charged Weyl = 7, x2 (anti), x3 gen = 42
    print(f"    raw charged-Weyl mode count (3 gen, +anti) = {charged_weyl_3gen} != 31, so 31 is not that count either;")
    print(f"    '2*16-1=31' conflates the 16 RECORD words with the fermion loop (a category slip).")
    check(charged_weyl_3gen != 31, "31 is not the raw charged-fermion mode count (42); nor the charge-weighted sum (8/16)")
    print("    => N1=31 matches no clean charge-weighted OR raw field-content count -> it is a FIT (K5 confirmed).")

    print(f"""
{"=" * 98}
VERDICT (exit 0):  the CIRCULARITY is resolved, but a non-circular route does NOT reach +0.036 -- it
gives the right ORDER and undershoots by ~2-4x; the exact value's kernel (N1=31) is a fit.

  CIRCULARITY RESOLVED. The dressed shift's structure delta = N1*alpha0/(2pi) is NOT circular: alpha0
  is the separately-derived bare value (R14) and N1 is a pure alpha-independent loop kernel. The
  circularity I flagged earlier was specific to the web/vertex route (t=alpha0). So "the magnitude is
  set by the vertex ~alpha, circular" is too pessimistic -- the perturbative structure is fine.

  THE GENUINE NON-CIRCULAR KERNEL UNDERSHOOTS. The photon self-energy is a CHARGE-WEIGHTED fermion
  loop, so the field-content kernel is Sigma Q_f^2 N_c = 8 (Dirac, 3 gen) or Sum Q^2 = 16 (Weyl). These
  give delta = 0.009 and 0.019 -- the RIGHT ORDER (delta ~ a few x alpha0) but UNDERSHOOTING the
  physical +0.036 by ~4x and ~2x. So a genuinely non-circular 1-loop reaches the right ballpark, not
  the number.

  N1=31 IS A FIT. The kernel that reproduces +0.036 is N1~31, which is ~4x the Dirac sum / ~2x the Weyl
  sum, and is neither the charge-weighted coefficient nor a clean raw mode count (the charged-Weyl
  count is 42; '2x16-1' mis-uses the record-16). So 31 has no non-circular field-content derivation --
  exactly the K5 'count != integral' caveat, now sharpened to a specific ~2-4x gap.

  NET: the dressed-alpha magnitude is non-circular IN PRINCIPLE and right ORDER from the field content,
  but the exact +0.036 needs a kernel (~31) that the genuine charge-weighted loop does not supply (it
  gives ~8-16). The open object is no longer "circularity" but this ~2-4x kernel gap -- plausibly
  higher-loop / chiral / colour structure, not pinned. 137.036 stays a fit (Part-12) until the kernel
  is derived charge-weighted.
{"=" * 98}""")
    print(f"exit 0 -- non-circular structure CONFIRMED (delta=N1*alpha0/2pi, bare x pure kernel); genuine "
          f"charge-weighted kernel = 8-16 -> delta 0.009-0.019 (right order, undershoots 0.036 by 2-4x); "
          f"N1=31 is a fit (no charge-weighted/raw-count derivation). Magnitude open = the kernel gap, not circularity.")


if __name__ == "__main__":
    main()
