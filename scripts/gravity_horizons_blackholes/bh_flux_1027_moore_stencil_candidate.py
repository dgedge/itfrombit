#!/usr/bin/env python3
r"""Black-hole flux coefficient: 10/27 Moore-stencil candidate.

The absolute Hawking flux gate reduces the coefficient to an attempt-rate
theorem:

    Gamma0 = (10/27) alpha0 Lambda_QCD

lands at 0.997 of the Stefan-Hawking coefficient once the beta-one shell and
Schwarzschild escape-cone scaling are used.  The open question is whether
10/27 is an actual service-geometry count or just a numerical near-hit.

This audit supplies a concrete, falsifiable candidate:

  * the local causal service alphabet is the 3 x 3 x 3 Moore stencil;
  * outgoing service across a horizon shell selects the 3 x 3 outward face;
  * the horizon latch contributes one local stay/commit slot.

Then

    C_flux = (9 outward face slots + 1 latch slot) / 27 = 10/27.

This is NOT a derivation.  It is a candidate operator statement: prove that the
localized horizon scheduler uses the 27-slot Moore service alphabet and that
outward emission is exactly "outward face plus latch", or reject the 10/27
coefficient and return to the general greybody/attempt-rate theorem.
"""

from __future__ import annotations

import itertools
import math


PHI = (math.sqrt(5.0) - 1.0) / 2.0
EPS_F = 1.0 / (2.0 * PHI)
ALPHA0 = 1.0 / 137.0
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def source_ladder(beta_eff: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta_eff * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return dict(sorted((f, v / z) for f, v in raw.items()))


def mean_f(beta_eff: float) -> float:
    src = source_ladder(beta_eff)
    return sum(f * p for f, p in src.items())


def required_gamma_for_stefan(beta_eff: float) -> float:
    """Gamma0/Lambda required for P=M_P^4/(15360 pi M^2)."""

    rho_lambda = beta_eff / (2.0 * math.pi * EPS_F)
    target_c = 1.0 / (15360.0 * math.pi)
    denom = (27.0 * math.pi / 32.0) * EPS_F * mean_f(beta_eff) * rho_lambda**4
    return target_c / denom


def model_over_target(gamma: float, beta_eff: float = 1.0) -> float:
    return gamma / required_gamma_for_stefan(beta_eff)


def moore_stencil() -> list[tuple[int, int, int]]:
    return list(itertools.product((-1, 0, 1), repeat=3))


def main() -> None:
    print("BLACK-HOLE 10/27 FLUX STENCIL CANDIDATE")
    print("=" * 88)

    print("[1] Local service alphabet")
    stencil = moore_stencil()
    outward = [step for step in stencil if step[2] == 1]
    latch = [(0, 0, 0)]
    outward_plus_latch = sorted(set(outward + latch))
    pref = len(outward_plus_latch) / len(stencil)
    print(f"    Moore stencil slots        = {len(stencil)}")
    print(f"    outward face slots         = {len(outward)}")
    print(f"    latch/stay slots           = {len(latch)}")
    print(f"    outward + latch slots      = {len(outward_plus_latch)}")
    print(f"    C_flux                     = {pref:.12f} = 10/27")
    check(len(stencil) == 27, "3x3x3 local service alphabet has 27 slots")
    check(len(outward) == 9, "one oriented horizon face has 9 slots")
    check(len(outward_plus_latch) == 10, "outward face plus latch has 10 slots")
    check(abs(pref - 10.0 / 27.0) < 1.0e-15, "candidate coefficient is exactly 10/27")

    print("\n[2] Flux consequence at the beta-one shell")
    gamma_req = required_gamma_for_stefan(1.0)
    gamma_face = (9.0 / 27.0) * ALPHA0
    gamma_face_latch = (10.0 / 27.0) * ALPHA0
    gamma_non_inward = (18.0 / 27.0) * ALPHA0
    print(f"    Gamma_req/Lambda           = {gamma_req:.12e}")
    print(f"    outward face only          = {gamma_face:.12e}, P/P_SB={model_over_target(gamma_face):.9f}")
    print(f"    outward face + latch       = {gamma_face_latch:.12e}, P/P_SB={model_over_target(gamma_face_latch):.9f}")
    print(f"    non-inward control         = {gamma_non_inward:.12e}, P/P_SB={model_over_target(gamma_non_inward):.9f}")
    check(abs(model_over_target(gamma_face_latch) - 1.0) < 4.0e-3, "outward+latch candidate lands within 0.4 percent")
    check(abs(model_over_target(gamma_face) - 1.0) > 0.05, "face-only control misses at order ten percent")
    check(model_over_target(gamma_non_inward) > 1.7, "counting tangential slots as escaping strongly over-emits")

    print(
        """
[3] VERDICT
    The 10/27 coefficient is upgraded from a naked near-hit to a precise
    geometric candidate:

        10/27 = (outward Moore face + local horizon latch) / full 27-slot
                local service alphabet.

    This candidate script alone does not lock the Hawking flux coefficient.
    The follow-up transfer-theorem script supplies the conditional proof.  The
    theorem boundary is now sharper:

      1. prove the horizon service scheduler uses the 27-slot local
         causal/Moore alphabet;
      2. prove emitted service current is exactly outward face plus latch;
      3. include spin/partial-wave greybody transfer and rerun the coefficient
         gate if the target is not the pure Stefan-Hawking scalar coefficient.

ALL ASSERTIONS PASSED -- 10/27 is a concrete flux-stencil candidate; see the
transfer-theorem script for the conditional proof."""
    )


if __name__ == "__main__":
    main()
