#!/usr/bin/env python3
r"""Dressed fine-structure constant status consolidation.

This script records the current sharp status of

    alpha_0^{-1} = 137        versus        alpha_phys^{-1}=137.035999...

after the R14 bare-rate theorem and the dressed-alpha Ward/Kubo audits.

It deliberately does not try to re-fit the physical value.  Its job is to keep
the three logically different claims separated:

  1. Bare service weight 1/137: derived by the monitored record-pair alphabet.
  2. Monitor/web occupation near-hit: internally selected and continuum stable,
     but a monitored service-label occupation.
  3. Physical low-energy QED coupling: Ward/Kubo/Peierls billed by the
     current-current / photon self-energy slot, which gives a charge-weighted
     one-loop-sized kernel, not the mode-count N1=31.

Current conclusion: dressed alpha remains a bounded electromagnetic residual.
The Part-12 N1=31 Dyson-Schwinger line is retained only as a historical
mode-count fit unless a new sector-billing theorem maps the physical EM vertex
residue to the monitored service occupation.
"""

from __future__ import annotations

import math


ALPHA0_INV = 137.0
ALPHA_PHYS_INV = 137.035999084
DELTA_INV = ALPHA_PHYS_INV - ALPHA0_INV
ALPHA0 = 1.0 / ALPHA0_INV

# The physical inverse-alpha shift written as a one-loop-style kernel
# delta = N_eff * alpha0 / (2*pi).
N_REQUIRED = DELTA_INV * 2.0 * math.pi / ALPHA0

# Existing audit outputs.  These are not fitted here; they are the stable
# status numbers from the named scripts below.
MONITOR_OCCUPATION_RATIO = 0.995613
WARD_WEB_SELF_ENERGY_RATIO = 0.394072
WARD_BARE_SELF_ENERGY_RATIO = 2.694240
MAX_CHARGE_WEIGHTED_KERNEL = 17.62

# Standard charge-weighted EM self-energy kernels.
DIRAC_THREE_GEN_KERNEL = 8.0
WEYL_THREE_GEN_KERNEL = 16.0


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> None:
    print("DRESSED ALPHA STATUS CONSOLIDATION")
    print("=" * 88)
    print(f"bare service alpha0^-1       = {ALPHA0_INV:.0f}")
    print(f"physical alpha^-1            = {ALPHA_PHYS_INV:.9f}")
    print(f"inverse-alpha shift          = {DELTA_INV:.9f}")
    print(f"kernel required by shift     = {N_REQUIRED:.3f}")

    print("\n[1] Bare service rate")
    print("  R14 derives one firing projector in a 137-label monitored record-pair")
    print("  alphabet.  This closes alpha0=1/137 at the bare service layer.")
    check(abs(ALPHA0_INV - 137.0) < 1e-12, "bare service alphabet is 137")

    print("\n[2] Monitor/web occupation route")
    print(f"  monitored Peierls-link service occupation / observed shift = {MONITOR_OCCUPATION_RATIO:.6f}")
    print("  This is internally selected and continuum stable, but it is a pointer")
    print("  occupation, not yet the physical Ward/Kubo EM vertex residue.")
    check(0.99 < MONITOR_OCCUPATION_RATIO < 1.01, "monitor occupation gives the near-hit")

    print("\n[3] Physical Ward/Kubo self-energy route")
    print(f"  web-dressed Ward/Kubo self-energy / observed shift = {WARD_WEB_SELF_ENERGY_RATIO:.6f}")
    print(f"  bare Ward/Kubo self-energy / observed shift        = {WARD_BARE_SELF_ENERGY_RATIO:.6f}")
    print("  The self-energy slot is the QED-correct observable, but it gives a")
    print("  one-loop-sized undershoot after web dressing, not the near-hit.")
    check(0.3 < WARD_WEB_SELF_ENERGY_RATIO < 0.5, "physical web-dressed slot undershoots")
    check(abs(WARD_WEB_SELF_ENERGY_RATIO - MONITOR_OCCUPATION_RATIO) > 0.4, "service and self-energy slots are distinct")

    print("\n[4] N1=31 field-content test")
    print(f"  required N1                         = {N_REQUIRED:.2f}")
    print(f"  charge-weighted 3-generation Dirac  = {DIRAC_THREE_GEN_KERNEL:.2f}")
    print(f"  charge-weighted Weyl-counted max    = {WEYL_THREE_GEN_KERNEL:.2f}")
    print(f"  generous max with QCD/two-loop      = {MAX_CHARGE_WEIGHTED_KERNEL:.2f}")
    print("  N1=31 is a mode count (2*16-1), not a charge-weighted EM kernel.")
    check(30.0 < N_REQUIRED < 32.0, "observed shift requires the historical N1~31 kernel")
    check(MAX_CHARGE_WEIGHTED_KERNEL < 0.60 * N_REQUIRED, "charge-weighted kernels cannot reach N1=31")

    print("\nVERDICT")
    print("  Bare 137 is derived.  The dressed shift is not currently derived.")
    print("  The monitor/web route is the best internal near-hit, but the physical")
    print("  low-energy electromagnetic coupling is Ward/Kubo billed through the")
    print("  self-energy slot, and the honest charge-weighted kernel undershoots.")
    print("  Therefore alpha_phys^{-1}=137.035999... remains the one native")
    print("  electromagnetic number not produced by the substrate.  Closing it")
    print("  requires a new sector-billing theorem or a pinned Ward/Kubo")
    print("  normalisation/subtraction, not another mode-count fit.")
    print("exit 0 -- dressed alpha remains a bounded EM residual; N1=31 is not promoted.")


if __name__ == "__main__":
    main()
