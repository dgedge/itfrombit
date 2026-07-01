#!/usr/bin/env python3
r"""Does 2-loop + chiral/colour structure enhance the charge-weighted self-energy kernel 8->~31?

dressed_alpha_noncircular_attempt.py left a precise target: the physical dressed shift +0.036 needs a
kernel N1~31 (in delta=N1*alpha0/2pi), while the genuine non-circular 1-loop charge-weighted kernel is
Sigma Q_f^2 N_c = 8 (Dirac, 3 gen) or 16 (Weyl, chiral-doubled). The open hypothesis: do 2-loop QED +
QCD (colour) + chiral corrections multiply 8 up to ~31 (a 2-4x enhancement)?

This is decidable. The dressing is alpha-SUPPRESSED (0.036 ~ 5*alpha0), so the loop expansion is
PERTURBATIVE; a 2-4x enhancement would require the 2-loop term to be 100-300% of the 1-loop, i.e. the
expansion to diverge. The genuine perturbative corrections are tiny:
  - 2-loop QED:  ~ (alpha/pi) * (Sigma Q_f^4 / Sigma Q_f^2)  ~ 0.1-0.2%
  - QCD on quark loops:  (1 + alpha_s/pi + ...)  ~ +7-16% on the quark part
  - chiral:  Dirac->Weyl is the factor 2 already in the 8->16 step (not an extra enhancement).
So the maximal genuine charge-weighted kernel is ~16-18, giving delta ~0.020 -- still 1.8-3.6x short of
the N1=31 the fit needs. The 2-4x gap is NOT a higher-order self-energy effect.
"""
from __future__ import annotations

import math

ALPHA0 = 1.0 / 137.0
TWO_PI = 2 * math.pi
DELTA_PHYS = 0.036
N1_TARGET = DELTA_PHYS * TWO_PI / ALPHA0          # = 31


def kernel_to_delta(K):
    return K * ALPHA0 / TWO_PI


def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c:
        raise AssertionError(m)


def main():
    print("2-LOOP / chiral / colour enhancement of the charge-weighted self-energy kernel — decisive")
    print("=" * 98)
    print(f"\n  target: physical delta = {DELTA_PHYS} needs kernel N1 = {N1_TARGET:.1f}; 1-loop Dirac kernel = 8")

    # ---- charge sums (3 generations) ----
    Q = {"lep": (-1.0, 1, 3), "up": (2/3, 3, 3), "down": (-1/3, 3, 3)}   # (charge, N_c, n_gen)
    sumQ2 = sum(q * q * nc * ng for q, nc, ng in Q.values())            # Dirac = 8
    sumQ4 = sum(q ** 4 * nc * ng for q, nc, ng in Q.values())          # for the 2-loop coefficient
    sumQ2_quark = sum(q * q * nc * ng for k, (q, nc, ng) in Q.items() if k != "lep")  # quark part = 5
    print(f"\n[1] 1-loop charge-weighted kernels (non-circular, from field content):")
    print(f"    Sigma Q_f^2 N_c (Dirac) = {sumQ2:.2f}   Sigma Q_f^4 N_c = {sumQ4:.2f}   quark part Sigma Q^2 = {sumQ2_quark:.2f}")
    print(f"    Dirac 1-loop -> delta = {kernel_to_delta(sumQ2):.4f};  Weyl (x2 chiral) 16 -> delta = {kernel_to_delta(2*sumQ2):.4f}")
    check(abs(sumQ2 - 8) < 1e-9, "Dirac 1-loop charge-weighted kernel = 8")

    # ---- 2-loop QED enhancement: relative size ~ (alpha/pi) * (Sigma Q^4 / Sigma Q^2) ----
    two_loop_rel = (ALPHA0 / math.pi) * (sumQ4 / sumQ2)                 # ~ O(alpha/pi)
    print(f"\n[2] 2-loop QED enhancement (alpha-suppressed): ~ (alpha/pi)*(SumQ^4/SumQ^2) = {two_loop_rel:.5f}"
          f"  = {two_loop_rel*100:.3f}%")
    check(two_loop_rel < 0.01, "2-loop QED correction is < 1% (alpha-suppressed) — nowhere near a 2-4x enhancement")

    # ---- QCD on quark loops: (1 + alpha_s/pi) at a low scale; enhances only the quark part ----
    for alpha_s in (0.30, 0.50):
        qcd_factor = 1 + alpha_s / math.pi                              # leading R-ratio correction
        K_qcd = (sumQ2 - sumQ2_quark) + sumQ2_quark * qcd_factor        # leptons unchanged, quarks enhanced
        K_full = K_qcd * (1 + two_loop_rel)                             # + 2-loop QED on everything
        K_weyl = 2 * (sumQ2 - sumQ2_quark) + 2 * sumQ2_quark * qcd_factor  # Weyl (chiral) version
        print(f"\n[3] alpha_s={alpha_s}: QCD factor (1+a_s/pi)={qcd_factor:.3f}")
        print(f"    Dirac kernel 8 -> {K_full:.2f}  (delta={kernel_to_delta(K_full):.4f});  "
              f"Weyl kernel 16 -> {K_weyl:.2f}  (delta={kernel_to_delta(K_weyl):.4f})")
        print(f"    vs target N1=31 (delta=0.036): short by {N1_TARGET/K_full:.2f}x (Dirac), {N1_TARGET/K_weyl:.2f}x (Weyl)")

    # use the generous case (alpha_s=0.5, Weyl) as the MAXIMAL genuine kernel
    qcd = 1 + 0.5 / math.pi
    K_max = 2 * (sumQ2 - sumQ2_quark) + 2 * sumQ2_quark * qcd
    K_max *= (1 + two_loop_rel)
    print(f"\n[4] MAXIMAL genuine charge-weighted kernel (Weyl + QCD@a_s=0.5 + 2-loop) = {K_max:.2f}")
    print(f"    -> delta = {kernel_to_delta(K_max):.4f}  (target 0.036);  still short by {N1_TARGET/K_max:.2f}x")
    check(K_max < 20, "even the maximal genuine charge-weighted kernel is ~18 << 31 — the gap is NOT closed by loops")
    # the decisive logical point: a 2-4x enhancement needs the 2-loop ~100-300% of 1-loop
    needed_2loop_rel = N1_TARGET / sumQ2 - 1                            # fractional enhancement 8->31
    print(f"\n[5] to take the 1-loop 8 up to 31 needs a {needed_2loop_rel*100:.0f}% enhancement, i.e. the 'higher-order'")
    print(f"    term must be ~{needed_2loop_rel:.1f}x the 1-loop — a DIVERGENT expansion, contradicting the dressing")
    print(f"    being alpha-suppressed (0.036 ~ 5*alpha0, manifestly perturbative). So it cannot be loops.")
    check(needed_2loop_rel > 2.0, "the required enhancement (>200%) is non-perturbative — impossible for an alpha-suppressed dressing")

    print(f"""
{"=" * 98}
VERDICT (exit 0):  DECISIVE NO. The genuine 2-loop + chiral + colour charge-weighted self-energy does
NOT multiply the 1-loop kernel up to ~31. Perturbative corrections are ~10-20%, not 2-4x.

  THE NUMBERS. 1-loop charge-weighted kernel = 8 (Dirac) / 16 (Weyl, chiral-doubled). 2-loop QED is
  ~{two_loop_rel*100:.2f}% (alpha-suppressed). QCD on the quark loops adds ~7-16% (the standard 1+alpha_s/pi). The
  MAXIMAL genuine kernel (Weyl + QCD@strong-coupling + 2-loop) is ~{K_max:.0f}, giving delta ~ {kernel_to_delta(K_max):.3f} --
  still {N1_TARGET/K_max:.1f}x short of the physical 0.036 (N1=31). Chiral structure is the 8->16 factor 2,
  already counted; it does not give a further enhancement.

  THE LOGIC. Reaching 31 from 8 needs a ~{needed_2loop_rel*100:.0f}% enhancement -- the higher-order term would have to be
  ~{needed_2loop_rel:.0f}x the 1-loop, i.e. a DIVERGENT loop expansion. But the dressing is alpha-suppressed
  (0.036 ~ 5*alpha0), so the expansion is manifestly perturbative and every higher order is ~alpha/pi ~
  0.2% smaller. A 2-4x enhancement is therefore logically excluded, not merely unobserved.

  CONSEQUENCE. The 2-4x kernel gap is NOT a higher-order charge-weighted self-energy effect. N1=31 is
  ~2x the maximal genuine charge-weighted kernel (~16-18) and is therefore confirmed a FIT, not the
  substrate's physical vacuum polarisation at any loop order. The dressed-alpha magnitude +0.036 is NOT
  derivable from the charge-weighted self-energy; the honest status is settled: 137.036 is matched by
  the Part-12 count (a fit), the genuine non-circular self-energy gives only the right ORDER
  (delta ~ 0.01-0.02), and no perturbative route closes the remaining ~2x. The magnitude is a genuine,
  bounded open residual -- the dressing is the one EM-observable number the substrate does not derive.
{"=" * 98}""")
    print(f"exit 0 -- DECISIVE NO: 2-loop QED ~{two_loop_rel*100:.2f}%, QCD ~10%, max genuine kernel ~{K_max:.0f} (delta~{kernel_to_delta(K_max):.3f}); "
          f"reaching 31 needs ~{needed_2loop_rel*100:.0f}% (divergent), excluded for an alpha-suppressed dressing. N1=31 is a fit; +0.036 not loop-derivable.")


if __name__ == "__main__":
    main()
