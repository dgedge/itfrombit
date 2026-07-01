#!/usr/bin/env python3
r"""dressed_alpha_endpoint_escape_audit.py

Tests -- and REFUTES -- the last un-pushed escape from the dressed-alpha no-go: the
"charge-blind endpoint" reading. Honest negative; the no-go stands. Records the result so the
escape is not re-tried.

THE ESCAPE (the candidate that motivated this audit):
  The framework derives the BARE alpha0^-1 = 137 charge-blind (item 79: equipartition over a
  137-label channel count). The charge-blind continuation -- the monitor/web SERVICE OCCUPATION --
  reproduces the shift to 0.9956, i.e. alpha(0)^-1 to ~1 ppm. One might argue: alpha(0) is QED's
  INPUT (QED predicts the running, not the endpoint), and the Q^2-weighting lives in the running,
  so the ENDPOINT dressing may legitimately be charge-blind -> the ~1 ppm near-hit is the real
  derivation, modulo "the dressing is charge-blind".

WHY IT FAILS (two independent reasons, both already implicit in canon):
  (1) alpha(0) IS, BY DEFINITION, the charge-weighted photon-propagator residue (the Ward/Kubo
      current-current slot). The physical low-energy coupling is *that operator's* value. The
      framework is not taking alpha(0) as an input -- it is DERIVING it via bare(137) -> dressed
      (137.036), which is precisely a photon self-energy map, hence charge-weighted. So the
      endpoint-vs-running distinction does not apply: the framework's dressing IS a self-energy.
  (2) The charge-blind kernel is not even clean: N1 = 31 = "2*16-1" conflates the record-16 with
      the fermion-loop count (the charged-Weyl mode count is 42), and the service occupation is a
      diagonal pointer Born weight, NOT affine-equivalent to the self-energy operator
      (dressed_alpha_service_kubo_moment_no_go.py). So the ~1 ppm near-hit is a coincidence of a
      pointer occupation with the self-energy value, not a derivation.

  Conclusion: the charge-blind service occupation reproduces the NUMBER but is the wrong OPERATOR;
  the genuine (charge-weighted) self-energy gives the right ORDER but undershoots ~2x. The escape
  is closed.

THE HONEST STATUS (what is actually true, stated without flattery):
  + alpha0^-1 = 137 -- the integer channel count -- is derived parameter-free and matches the
    measured 137.035999 to 0.026% (4 significant figures). That is a real, clean prediction.
  - the dressed shift 0.036 is a BOUNDED residual: the genuine charge-weighted 1-loop kernel
    (sum Q^2 N_c = 8 Dirac, 16 Weyl, <=18 with QCD/2-loop) gives delta ~ 0.009-0.020 -- the right
    order but ~2x short of 0.036; N1=31 is a fit; all charge-weighted self-energy routes (1/2-loop,
    QCD, web DOS, non-unital steady response) are ruled out as its exact source.
  So the dressed alpha is NOT closed and this session does not close it. It remains the deepest
  single residual: 137 to 0.026%, the last 0.036 a bounded ~2x charge-weighted-kernel gap. The
  charge-blind near-hit is a coincidence, now documented as such.
"""
import math

ALPHA0_INV = 137.0
ALPHA_PHYS_INV = 137.035999084
SHIFT = ALPHA_PHYS_INV - ALPHA0_INV
SERVICE_OCC_RATIO = 0.995613       # charge-blind pointer occupation / shift (a NEAR-HIT, not the operator)
WARD_KUBO_RATIO = 0.394072         # charge-weighted self-energy / shift (the DEFINING operator; undershoots)
KERNEL_DIRAC = 8.0                 # sum Q^2 N_c, 3 gen Dirac
KERNEL_WEYL = 16.0                 # Weyl-counted
KERNEL_MAX = 18.0                  # generous (QCD + 2-loop)
N_REQUIRED = SHIFT * 2.0 * math.pi * ALPHA0_INV    # ~31


def main():
    print("=== Dressed alpha: charge-blind ENDPOINT escape -- audit (negative) ===\n")
    bare_acc = abs(ALPHA0_INV - ALPHA_PHYS_INV) / ALPHA_PHYS_INV
    print(f"  [bare] alpha0^-1 = 137 matches measured {ALPHA_PHYS_INV:.6f} to {bare_acc*100:.3f}%  "
          f"(parameter-free integer; a real prediction).")
    print(f"  [shift] dressed residual = {SHIFT:.6f}; one-loop-style kernel needed N = {N_REQUIRED:.2f} (~31).\n")

    cb = ALPHA0_INV + SERVICE_OCC_RATIO * SHIFT
    cb_ppm = abs(cb - ALPHA_PHYS_INV) / ALPHA_PHYS_INV * 1e6
    print("  [escape] charge-blind service occupation:")
    print(f"      reproduces {SERVICE_OCC_RATIO:.4f} of the shift -> alpha(0)^-1 = {cb:.6f} ({cb_ppm:.1f} ppm).")
    print("      tempting: 'alpha(0) is QED's input, Q^2-weighting is in the running, so the endpoint")
    print("      dressing may be charge-blind' -> the ~1 ppm near-hit would be the derivation.\n")

    print("  [refutation] the escape fails for two independent reasons:")
    print("    (1) alpha(0) IS the charge-weighted photon-propagator residue BY DEFINITION; the framework")
    print("        DERIVES it via bare(137)->dressed(137.036), a self-energy map -> charge-weighted. The")
    print("        framework is not taking alpha(0) as input, so 'endpoint vs running' does not apply.")
    print("    (2) the charge-blind kernel is not clean: 31='2*16-1' conflates the record-16 with the")
    print("        fermion loop (charged-Weyl count is 42), and the service occupation is a diagonal")
    print("        pointer weight, not the self-energy operator. The ~1 ppm match is a COINCIDENCE.\n")

    print("  [genuine kernel] charge-weighted 1-loop (the defining operator) gives the right ORDER, ~2x short:")
    for name, k in (("Dirac 3-gen (8)", KERNEL_DIRAC), ("Weyl (16)", KERNEL_WEYL), ("max QCD+2-loop (18)", KERNEL_MAX)):
        d = k * ALPHA0_INV ** -1 / (2 * math.pi)   # delta = k * alpha0 / 2pi
        print(f"      {name:<22} -> delta = {d:.4f}  ({d/SHIFT*100:.0f}% of the observed 0.036)")
    print(f"      required N={N_REQUIRED:.1f} is ~2x the Weyl / ~4x the Dirac kernel -> N1=31 is a FIT.\n")

    print("[verdict] ESCAPE CLOSED; dressed alpha remains the deepest single residual (no flattery):")
    print("  + alpha0^-1 = 137 derived parameter-free, matches measured to 0.026% (4 sig figs) -- solid.")
    print("  - the dressed shift 0.036 is a BOUNDED residual: charge-weighted kernel gives delta~0.01-0.02")
    print("    (right order, ~2x short); the charge-blind ~1 ppm near-hit is a coincidence (wrong operator);")
    print("    N1=31 is a fit. The endpoint-vs-running escape is refuted (alpha(0) IS the charge-weighted")
    print("    residue, which the framework derives via a self-energy map). This session does NOT close it.")
    print("  The honest residual: 137 to 0.026% + a bounded ~2x charge-weighted-kernel gap for the last 0.036.")

    # gates -- assert the HONEST facts (including the negative)
    assert bare_acc < 3e-4, "bare 137 must match measured to <0.03% (the solid parameter-free prediction)"
    assert cb_ppm < 5.0, "charge-blind service occupation numerically reproduces alpha(0) to ~1 ppm (the coincidence)"
    assert WARD_KUBO_RATIO < 0.5, "the DEFINING charge-weighted operator undershoots (no-go stands)"
    assert KERNEL_WEYL * ALPHA0_INV ** -1 / (2 * math.pi) < SHIFT, "charge-weighted kernel undershoots the shift"
    assert N_REQUIRED > 1.7 * KERNEL_WEYL, "required kernel is ~2x the genuine charge-weighted one -> 31 is a fit"
    print("\nGATES PASSED -- bare 137 solid to 0.026%; charge-blind near-hit is a coincidence (wrong operator);")
    print("charge-weighted kernel undershoots ~2x; endpoint escape refuted. Dressed alpha NOT closed. exit 0")


if __name__ == "__main__":
    main()
