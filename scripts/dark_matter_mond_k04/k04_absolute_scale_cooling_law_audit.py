#!/usr/bin/env python3
r"""K04 absolute-scale / cooling-law audit.

The K04 residual has three pieces that are easy to mix up:

  1. the dimensionless local ratio w4/w6;
  2. the absolute energy bridge w6/Lambda;
  3. the physical cooling law that turns a boot history into a Kibble-Zurek
     defect density.

This script consolidates the current canon status.  It is deliberately a
status resolver, not a new large simulation:

  * w4/w6 = 2 is conditionally derived by Q3 loop-orbit edge-ledger
    normal-ordering.
  * w6/Lambda cannot be derived by K04 Metropolis dynamics alone, because
    (w4,w6,lambda,T) -> s(w4,w6,lambda,T) is an exact symmetry of every
    acceptance probability.
  * gamma = 0.995 is a provenance-strong simulation protocol.  Its closeness
    to 2^{-alpha0} is real, but the K04 entropy-spike gate failed; the
    physical cosmological cooling law is instead the boundary-printer dilution
    dF/dt = -n H F, with residual inputs H_c, mu, and w6/Lambda.

Exit 0 means the arithmetic behind that split is internally consistent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction


ALPHA0 = 1.0 / 137.035999084
GAMMA_PAPER = 0.995
T_START_OVER_W6 = 6.0
T_END_OVER_W6 = 0.5
LAMBDA_QCD_GEV = 0.332
DEBRIS_EXCESS_PER_VERTEX_W6 = 2.17

ENTROPY_TARGET = 2.508170
ENTROPY_LATENT_BEST = 0.870


@dataclass(frozen=True)
class Anchor:
    name: str
    t_over_w6: float
    status: str

    @property
    def w6_over_lambda(self) -> float:
        return 1.0 / self.t_over_w6

    @property
    def debris_over_lambda(self) -> float:
        return DEBRIS_EXCESS_PER_VERTEX_W6 * self.w6_over_lambda


def metropolis_accept(delta_e: float, temperature: float) -> float:
    return min(1.0, math.exp(-delta_e / temperature))


def assert_metropolis_scale_invariance() -> None:
    for delta_e, temperature, scale in [
        (1.0, 6.0, 11.0),
        (4.0, 2.8, 0.37),
        (16.0, 3.15, 29.0),
    ]:
        base = metropolis_accept(delta_e, temperature)
        scaled = metropolis_accept(scale * delta_e, scale * temperature)
        assert abs(base - scaled) < 1e-15


def k04_weight_ratio() -> tuple[Fraction, Fraction, Fraction]:
    """Return (w4, w6, w4/w6) under the orbit normal-ordering rule."""

    c4_weight = Fraction(1, 2)  # one orbit; each edge is used by 2 square loops
    c6_antipodal = Fraction(1, 2)  # 4 loops; each edge incidence 2
    c6_adjacent = Fraction(1, 6)  # 12 loops; each edge incidence 6
    c6_weight = (4 * c6_antipodal + 12 * c6_adjacent) / 16
    ratio = c4_weight / c6_weight
    assert c6_weight == Fraction(1, 4)
    assert ratio == 2
    return c4_weight, c6_weight, ratio


def gamma_to_ramp(gamma: float) -> float:
    return math.log(T_START_OVER_W6 / T_END_OVER_W6) / -math.log(gamma)


def printer_cooling_rate(n: int = 3, h_rate: float = 1.0) -> float:
    """Check dF/dt = -n H F over one short finite difference."""

    f4 = math.exp(-n * h_rate * 0.4)
    f5 = math.exp(-n * h_rate * 0.5)
    measured = -(math.log(f5) - math.log(f4)) / 0.1
    assert abs(measured - n * h_rate) < 1e-12
    return measured


def main() -> None:
    assert_metropolis_scale_invariance()
    w4, w6, ratio = k04_weight_ratio()
    measured_printer_rate = printer_cooling_rate()

    anchors = [
        Anchor("ramp start", 6.0, "conditional bridge: T_start = Lambda"),
        Anchor("transition low", 2.8, "alternative postulate: T_c = Lambda"),
        Anchor("transition high", 3.0, "alternative postulate: T_c = Lambda"),
        Anchor("latent crossing", 3.15, "failed entropy-gamma gate; not selected"),
    ]
    selected = anchors[0]
    assert abs(selected.w6_over_lambda - (1.0 / 6.0)) < 1e-15

    gamma_alpha = math.exp(-ALPHA0 * math.log(2.0))
    r_paper = gamma_to_ramp(GAMMA_PAPER)
    r_alpha = gamma_to_ramp(gamma_alpha)
    entropy_miss = ENTROPY_LATENT_BEST / ENTROPY_TARGET - 1.0
    assert abs(r_alpha / r_paper - 1.0) < 0.02
    assert entropy_miss < -0.5

    print("[0] K04 absolute-scale / cooling-law audit")
    print("    Metropolis scale invariance: PASS")
    print("    status: K04 temperatures are ratios; absolute w6 needs a bridge.")
    print()

    print("[1] w4/w6")
    print(f"    C4 edge-normal-ordered weight w4 = {w4}")
    print(f"    C6 coarse-grained orbit weight w6 = {w6}")
    print(f"    w4/w6 = {ratio}")
    print("    tier: conditionally derived under local edge-ledger normal ordering.")
    print()

    print("[2] w6/Lambda candidate anchors")
    print(
        f"    {'anchor':<18s} {'T/w6':>8s} {'w6/Lambda':>12s} "
        f"{'e_D/Lambda':>12s}  status"
    )
    for anchor in anchors:
        chosen = " SELECTED" if anchor == selected else ""
        print(
            f"    {anchor.name:<18s} {anchor.t_over_w6:8.3f} "
            f"{anchor.w6_over_lambda:12.6f} {anchor.debris_over_lambda:12.6f}  "
            f"{anchor.status}{chosen}"
        )
    print(
        f"    selected w6 = {selected.w6_over_lambda * LAMBDA_QCD_GEV:.6f} GeV "
        f"for Lambda = {LAMBDA_QCD_GEV:.3f} GeV"
    )
    print("    tier: conditional bridge, not a K04-internal derivation.")
    print()

    print("[3] gamma / cooling")
    print(f"    protocol gamma = {GAMMA_PAPER:.6f} -> R = {r_paper:.3f} sweeps")
    print(f"    alpha-bit gamma = exp(-alpha0 ln2) = {gamma_alpha:.9f} -> R = {r_alpha:.3f}")
    print(
        f"    entropy-spike gate: target {ENTROPY_TARGET:.6f}, "
        f"latent estimator {ENTROPY_LATENT_BEST:.6f}, miss {100*entropy_miss:.1f}%"
    )
    print("    verdict: alpha-bit proximity is not a derivation of the protocol gamma.")
    print()

    print("[4] physical cosmological cooling law")
    print("    boundary printer mints fresh S=0 cells and dilutes frustration:")
    print("        dF/dt = -n H F, so |d ln F/dt| = nH")
    print(f"    finite-difference check at n=3 gives {measured_printer_rate:.3f} H")
    print("    KZ quench time is tau_Q = 1/(n H_c), not a free Metropolis gamma.")
    print("    residuals: crystallisation epoch H_c, KZ exponent mu, and w6/Lambda.")
    print()

    print("[5] canon verdict")
    print("    CLOSED/conditional: w4/w6 = 2 under edge-ledger normal ordering.")
    print("    NOT K04-derived: w6/Lambda, by exact Metropolis scale symmetry.")
    print("    SUPERSEDED target: deriving gamma=0.995 from the K04 entropy spike.")
    print("    LIVE cooling law: boundary-printer dilution; gamma is a simulation proxy.")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
