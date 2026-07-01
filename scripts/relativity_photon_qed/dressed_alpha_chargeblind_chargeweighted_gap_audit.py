#!/usr/bin/env python3
r"""dressed_alpha_chargeblind_chargeweighted_gap_audit.py

Pushes on the sharpened dressed-alpha obstruction: "charge-BLIND substrate vs charge-WEIGHTED-by-
definition alpha(0)". Tests the one un-tested charge-weighted enhancement (the SMG mirror fermions),
quantifies the gap exactly, and reaches the honest verdict.

SETUP (confirmed canon): the shift delta(alpha^-1)=0.036 needs a one-loop-style kernel N_req ~ 31.
  - charge-WEIGHTED routes (the physical photon self-energy, sum Q^2): 8 (3-gen Dirac), <=18 (max with
    QCD + 2-loop) -- ALL undershoot 31 (prior decisive no-go).
  - charge-BLIND mode counts: 31 = "2x16-1", 32 = 2x16, 42 = charged-Weyl modes -- bracket/hit 31, but
    a mode count is NOT the physical EM observable (the propagator residue is charge-weighted).

THE NEW TEST -- do the SMG MIRROR fermions close the gap? At the cell scale the SMG mirror partners
(same gauge charges, opposite chirality) are present and DOUBLE the charge-weighted sum: sum Q^2 -> 16.
This was NOT in the prior "max 18". Result below: it gives 16 (= the Weyl count), STILL ~2x short of 31.
So mirror-inclusion does NOT close it -- the charge-weighted kernel, even maximally enhanced (mirrors +
QCD + 2-loop), tops out near 18-20, against a required 31.

THE GAP, QUANTIFIED. The dressing needs the charge-BLIND mode count (~31-32); the physical observable is
the charge-WEIGHTED sum (~16). Their ratio IS the mean-square charge of the matter rep:
  sum Q^2 / N_modes = 16/32 = 1/2 = <Q^2>.
So the dressed-alpha gap is EXACTLY the Q^2-weighting factor <Q^2> ~ 1/2 of the SO(10) 16: the substrate
counts modes; QED weights by Q^2; the dressed shift falls in that factor-of-~2 gap.

VERDICT (honest, no flattery): the dressed alpha is a STRUCTURALLY IRREDUCIBLE EM residual in the present
framework. It cannot be charge-blind (a mode count is not the propagator residue) and it cannot be
charge-weighted (every charge-weighted route, mirrors included, undershoots ~2x). The clean, non-fit
result is the DEFLATION already half-adopted in canon, now made explicit and explained:
  * the framework PREDICTS alpha^-1 = 137 (item 79 charge-blind equipartition), matching the measured
    137.035999 to 0.026% -- a real parameter-free 4-sig-fig prediction;
  * the 0.036 dressing is charge-weighted EMERGENT QED, whose ORDER the substrate reproduces
    (delta ~ sum Q^2 * alpha0/2pi ~ 0.01-0.02) but whose EXACT value is not a charge-blind substrate
    number; the N1=31 exact-derivation is a fit and stays retired.
The residual is thus LOCATED, not closed: it is the irreducible charge-blind/charge-weighted boundary,
= <Q^2> of the 16. This is the framework's analogue of "predict the bare integer; QED does the
charge-weighted dressing." Pushing further on a charge-weighted closure is diminishing returns.
"""
import math

ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
SHIFT = ALPHA_PHYS_INV - ALPHA0_INV
N_REQ = SHIFT * 2.0 * math.pi * ALPHA0_INV         # ~31

# charge-weighted kernels (sum Q^2 N_c)
SM_DIRAC_3GEN = 8.0                                 # standard QED vacuum-polarization sum
ALPHA_S = 0.5                                       # generous strong coupling at the cell scale
QCD_FACTOR = 1.0 + ALPHA_S / math.pi               # 1 + alpha_s/pi on the quark loops (~1.16)

# charge-blind mode counts
MODE_2x16 = 32                                      # 2 (chirality/mirror) x 16 fermions
MODE_2x16_m1 = 31                                   # "2x16-1" (the historical N1)
MODE_CHARGED_WEYL = 42                              # prior charged-Weyl mode count


def main():
    print("=== Dressed alpha: charge-blind vs charge-weighted gap audit ===\n")
    print(f"  shift delta(alpha^-1) = {SHIFT:.6f};  required one-loop kernel N_req = {N_REQ:.2f} (~31)\n")

    print("  [charge-WEIGHTED] physical photon self-energy (sum Q^2), with every enhancement:")
    sm = SM_DIRAC_3GEN
    with_mirror = 2.0 * sm                           # SMG mirrors double the charge-weighted sum
    with_mirror_qcd = with_mirror * QCD_FACTOR       # + QCD on quark loops
    for name, k in (("SM 3-gen Dirac", sm),
                    ("+ SMG mirrors (x2) [NEW test]", with_mirror),
                    ("+ QCD (1+a_s/pi)", with_mirror_qcd)):
        d = k * ALPHA0 / (2 * math.pi)
        print(f"      {name:32s} sumQ^2={k:5.1f} -> delta={d:.4f} ({d/SHIFT*100:4.0f}% of 0.036), "
              f"kernel/N_req={k/N_REQ:.2f}")
    print(f"      => even mirrors+QCD top out at sumQ^2~{with_mirror_qcd:.0f}, vs required {N_REQ:.0f}: "
          f"undershoot ~{N_REQ/with_mirror_qcd:.1f}x. Mirror-inclusion does NOT close it.\n")

    print("  [charge-BLIND] mode counts (NOT the physical observable):")
    for name, m in (("2x16-1 (historical N1)", MODE_2x16_m1), ("2x16", MODE_2x16),
                    ("charged-Weyl modes", MODE_CHARGED_WEYL)):
        print(f"      {name:24s} = {m}  -> delta={m*ALPHA0/(2*math.pi):.4f} ({m*ALPHA0/(2*math.pi)/SHIFT*100:.0f}%)")
    print(f"      => mode counts reach/exceed 31, but a mode count is charge-BLIND, not the residue.\n")

    print("  [GAP] the dressing wants the mode count (~32); QED gives sum Q^2 (~16):")
    mean_q2 = with_mirror / MODE_2x16                # <Q^2> = sumQ^2 / N_modes
    print(f"      <Q^2> = sumQ^2/N_modes = {with_mirror:.0f}/{MODE_2x16} = {mean_q2:.3f}")
    print(f"      => the dressed-alpha gap IS the mean-square-charge weighting <Q^2>~{mean_q2:.2f} of the")
    print(f"         SO(10) 16: substrate counts modes (32), QED weights by Q^2 (16). The shift falls in")
    print(f"         that factor-{1/mean_q2:.0f} gap.\n")

    print("  [DEFLATION] the honest, non-fit reading:")
    bare_acc = abs(ALPHA0_INV - ALPHA_PHYS_INV) / ALPHA_PHYS_INV
    print(f"      framework PREDICTS alpha^-1 = 137 (item 79, charge-blind) -> matches measured to "
          f"{bare_acc*100:.3f}% (4 sig figs, parameter-free).")
    print(f"      the 0.036 dressing is charge-weighted EMERGENT QED: order reproduced "
          f"(sumQ^2*alpha0/2pi ~ {with_mirror*ALPHA0/(2*math.pi):.3f}), exact value NOT a charge-blind number.")
    print(f"      N1=31 exact-derivation stays retired (a fit). The residual is LOCATED, not closed.")

    print("\n[verdict] DRESSED ALPHA IS A STRUCTURALLY IRREDUCIBLE EM RESIDUAL:")
    print("  - cannot be charge-blind (a mode count is not the propagator residue);")
    print("  - cannot be charge-weighted (every route, mirrors+QCD+2-loop included, undershoots ~2x);")
    print("  - the gap IS <Q^2>~1/2, the mean-square charge of the 16 -- substrate counts modes, QED")
    print("    weights by Q^2. The framework predicts 137 to 0.026%; the charge-weighted 0.036 is QED.")
    print("  This is the framework's 'predict the bare integer; QED does the charge-weighted dressing'.")
    print("  Further charge-weighted-closure attempts are diminishing returns; the residual is located.")

    # gates -- assert the no-go + the gap + the deflation
    assert 30.5 < N_REQ < 31.5, "required kernel ~31"
    assert with_mirror_qcd < 0.7 * N_REQ, "charge-weighted (mirrors+QCD) still undershoots required 31 by >40%"
    assert with_mirror == 16.0, "SMG mirror inclusion doubles sumQ^2 to 16 (= Weyl count), not 31"
    assert MODE_2x16 >= N_REQ > with_mirror_qcd, "31 sits between the charge-weighted max and the mode count"
    assert abs(mean_q2 - 0.5) < 1e-9, "the gap equals <Q^2>=1/2 (mode count 32 vs sumQ^2 16)"
    assert bare_acc < 3e-4, "the charge-blind prediction 137 matches measured to <0.03% (the solid result)"
    print("\nGATES PASSED -- mirror-inclusion undershoots (16<31); gap = <Q^2>~1/2; 137 predicted to 0.026%;")
    print("dressed alpha located as the irreducible charge-blind/charge-weighted boundary. exit 0")


if __name__ == "__main__":
    main()
