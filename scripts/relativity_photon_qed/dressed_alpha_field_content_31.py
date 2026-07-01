#!/usr/bin/env python3
r"""The last door: is N1 = 31 a charge-weighted field-content sum, or a different invariant?

The Part-12 fit needs the self-energy kernel N1 = delta*2pi/alpha0 = 0.036*2pi*137 ~= 31. The genuine
photon self-energy is CHARGE-WEIGHTED: Pi ~ sum_f Q_f^2 N_c (each fermion enters as Q^2 via the
photon-fermion vertex). For the SM that sum is 8 (3-gen Dirac) / 16 (Weyl-counted). DRIFT 1429 already
showed no loop order takes 8 to 31. This script asks the remaining field-content question directly:

  Can ANY physical charged content give a charge-weighted sum = 31?  And if 31 is not charge-weighted,
  what invariant IS it?

The answers are decisive:
  [1] the SM charge-weighted self-energy kernel is 8 (Dirac) / 16 (Weyl); + W/charged-Higgs ~ 11-18;
  [2] reaching 31 as a genuine charge-weighted sum needs ~12 Dirac generations (3 is fixed) -- excluded;
  [3] 31 is a COUNT, not a charge-weighting: 31 = 2^5 - 1 = 2*16 - 1 = dim(SO(10) Dirac spinor) - 1;
  [4] the substrate's OWN current-weighted kernel (the Peierls-J spectral sum) is ~12 (physical,
      web-dressed) -- charge-weighted-like (8-16), NOT 31.
So no physical field content supplies a charge-weighted 31; the N1=31 fit uses a MODE COUNT where the
self-energy physics requires a charge-weighted sum -- a category mismatch, and that mismatch IS the gap.
"""
from __future__ import annotations

import numpy as np

import dressed_alpha_bridge_web_open_system as bw
import bridge_web_lindblad_keldysh_poles as poles


ALPHA0 = 1.0 / 137.0
TWO_PI = 2 * np.pi
DELTA_PHYS = 0.036
N1_TARGET = DELTA_PHYS * TWO_PI / ALPHA0          # ~= 31


def charge_weighted_sum(generations: float, weyl: bool = False) -> float:
    """sum_f Q_f^2 N_c over the SM charged fermions, x generations. Dirac unless weyl=True (x2)."""
    per_gen = 1 * (1.0) ** 2          # charged lepton e: Q=-1, N_c=1
    per_gen += 3 * (2.0 / 3) ** 2     # up-type:   Q=+2/3, N_c=3
    per_gen += 3 * (1.0 / 3) ** 2     # down-type: Q=-1/3, N_c=3
    s = per_gen * generations         # = 8/3 per gen Dirac
    return 2 * s if weyl else s


def substrate_current_kernel():
    """N1_sub from the substrate's Peierls current-current spectral sum (bare + web-dressed)."""
    H2, _pairs, idx, _ = bw.build_pair_system()
    t_m = 1.0 / 3.0
    evals = t_m * np.linalg.eigvalsh(H2)
    evals_u, evecs = np.linalg.eigh(H2)
    omega = t_m * evals_u - float((t_m * evals_u).min())
    J = bw.current_portal(idx, bw.ETA_PIN)
    amps = evecs.T @ J
    groups = poles.eigen_groups(t_m * evals_u)
    s_eta, _ = bw.photon_form_factor(n_grid=240)
    bare = web = 0.0
    for g in groups:
        w = float(np.sum(amps[g] ** 2))
        og = float(np.mean(omega[g]))
        if og < 1e-6:
            continue
        bare += w / og
        web += w * s_eta(og, bw.ETA_PIN) / og
    # delta = 137 * gv2 * disp ; N1 = delta*2pi/alpha0
    gv2 = ALPHA0 * bw.G_PORTAL ** 2
    n1_bare = (137.0 * gv2 * bare) * TWO_PI / ALPHA0
    n1_web = (137.0 * gv2 * web) * TWO_PI / ALPHA0
    return n1_bare, n1_web


def main() -> None:
    print("THE LAST DOOR — is N1=31 a charge-weighted field-content sum?")
    print("=" * 92)
    print(f"\n  target kernel N1 = delta*2pi/alpha0 = {N1_TARGET:.2f}  (Part-12 fit needs ~31)")

    # [1] the physical charge-weighted self-energy kernel
    s3_dirac = charge_weighted_sum(3)
    s3_weyl = charge_weighted_sum(3, weyl=True)
    print("\n[1] charge-weighted photon self-energy  sum_f Q_f^2 N_c  (the PHYSICAL kernel):")
    print(f"    per generation (Dirac)        = {charge_weighted_sum(1):.4f}  (= 8/3)")
    print(f"    3 generations (Dirac)         = {s3_dirac:.2f}   <- the physical 1-loop kernel")
    print(f"    3 generations (Weyl-counted)  = {s3_weyl:.2f}")
    # + W boson (Q=+-1) and a charged Higgs (Q=+-1), as generous extra charged content
    s_plus_wh = s3_dirac + 2 * 1.0 ** 2 + 1.0 ** 2     # naive Q^2 for W+- and H+-
    print(f"    3 gen + W^+- + H^+- (generous) = {s_plus_wh:.2f}   (still nowhere near 31)")

    # [2] what content would a genuine charge-weighted 31 require?
    gens_for_31 = N1_TARGET / charge_weighted_sum(1)
    print("\n[2] generations needed for a GENUINE charge-weighted sum = 31:")
    print(f"    31 / (8/3 per gen) = {gens_for_31:.1f} Dirac generations  (the substrate fixes 3) -> EXCLUDED")
    print(f"    even Weyl-counted: 31 / (16/3) = {N1_TARGET/charge_weighted_sum(1, weyl=True):.1f} gen -> still excluded")

    # [3] 31 is a COUNT, not a charge-weighting
    print("\n[3] 31 as a MODE COUNT (not charge-weighted):")
    print(f"    2^5 - 1                       = {2**5 - 1}")
    print(f"    2 x 16 - 1 (Weyl count x2 -1) = {2*16 - 1}")
    print(f"    dim(SO(10) Dirac spinor) - 1  = {2**5} - 1 = {2**5 - 1}   (32 = 16 + 16bar)")
    print("    -> 31 is the SO(10) Dirac-spinor dimension minus a singlet: a COUNT of modes, with NO")
    print("       charge weighting. The photon self-energy weights each mode by Q^2, which 31 does not.")

    # [4] the substrate's own current-weighted kernel
    n1_bare, n1_web = substrate_current_kernel()
    print("\n[4] the substrate's OWN current-weighted kernel (Peierls-J spectral sum):")
    print(f"    N1_sub (bare matter loop)     = {n1_bare:.1f}   (overshoots; omits the web DOS)")
    print(f"    N1_sub (web-dressed, physical)= {n1_web:.1f}   <- charge-weighted-like (cf. 8-16), NOT 31")

    # [5] verdict
    print("\n" + "=" * 92)
    print("VERDICT — the last door is closed:")
    print(f"  The physical photon self-energy is CHARGE-WEIGHTED: sum_f Q_f^2 N_c = {s3_dirac:.0f} (3-gen Dirac),")
    print(f"  {s3_weyl:.0f} (Weyl). Reaching a genuine charge-weighted 31 needs ~{gens_for_31:.0f} generations -- excluded")
    print("  (the substrate fixes 3). 31 is NOT a charge-weighted sum; it is a MODE COUNT (2^5-1 =")
    print("  2*16-1 = the SO(10) Dirac-spinor dimension minus a singlet). The substrate's own physical")
    print(f"  current-weighted kernel ({n1_web:.0f}) is charge-weighted-like (8-16), not 31 -- confirming it.")
    print("  So the field content does NOT supply a charge-weighted 31: the N1=31 fit reads a COUNT")
    print("  where the self-energy physics requires a CHARGE-WEIGHTED sum, and that category mismatch")
    print("  (count 31 vs charge-weighted 8-16) IS the 2-4x gap. No field-content door reaches 0.036;")
    print("  137.036 stays the Part-12 fit, now closed from the field-content side too.")


if __name__ == "__main__":
    main()
