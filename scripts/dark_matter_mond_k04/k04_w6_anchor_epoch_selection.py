#!/usr/bin/env python3
"""Select the K04 w6 anchor epoch under the current canon premises.

The K04 Metropolis dynamics determines only temperatures in units of w6.  The
absolute bridge has to identify one physical epoch temperature with the canon
scale Lambda.  This script audits the live candidates:

  * ramp start: the external bath/protocol temperature T_start = 6 w6;
  * transition/arrest: emergent K04 response temperatures near 2.8-3.0 w6;
  * latent crossing: the branch-resolved entropy-spike crossing near 3.15 w6.

Under the present canon premise "the boot bath is priced at Lambda", the
external ramp-start temperature is the selected anchor.  Transition anchoring is
not numerically forbidden, but it is a different physical premise:
"crystallisation itself occurs at Lambda."  K04 does not derive that premise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


LAMBDA_QCD_GEV = 0.332
W6_CARRIED_BAND = (0.05, 0.27)

T_RAMP_START = 6.0
T_TRANSITION_LOW = 2.8
T_TRANSITION_HIGH = 3.0
T_ARREST_ENTRY = 2.99
T_LATENT_CROSSING = 3.15

DEBRIS_EXCESS_PER_VERTEX_W6 = 2.17


@dataclass(frozen=True)
class AnchorCandidate:
    name: str
    t_over_w6: float
    external_bath: bool
    reading: str

    @property
    def w6_over_lambda(self) -> float:
        return 1.0 / self.t_over_w6

    @property
    def e_debris_over_lambda(self) -> float:
        return DEBRIS_EXCESS_PER_VERTEX_W6 * self.w6_over_lambda

    @property
    def in_carried_band(self) -> bool:
        lo, hi = W6_CARRIED_BAND
        return lo <= self.w6_over_lambda <= hi


def metropolis_accept(delta_e: float, temperature: float) -> float:
    return min(1.0, math.exp(-delta_e / temperature))


def verify_scale_invariance() -> None:
    for delta_e, temperature, scale in [
        (1.0, 6.0, 7.0),
        (4.0, 2.8, 13.0),
        (16.0, 3.15, 0.37),
    ]:
        base = metropolis_accept(delta_e, temperature)
        scaled = metropolis_accept(scale * delta_e, scale * temperature)
        assert abs(base - scaled) < 1e-15


def candidates() -> list[AnchorCandidate]:
    return [
        AnchorCandidate(
            "ramp start",
            T_RAMP_START,
            True,
            "external boot bath/protocol temperature",
        ),
        AnchorCandidate(
            "C_v transition low",
            T_TRANSITION_LOW,
            False,
            "emergent response; anchoring it is a new T_c=Lambda premise",
        ),
        AnchorCandidate(
            "C_v transition high",
            T_TRANSITION_HIGH,
            False,
            "emergent response; anchoring it is a new T_c=Lambda premise",
        ),
        AnchorCandidate(
            "arrest entry",
            T_ARREST_ENTRY,
            False,
            "emergent response; from the superseded gamma-bridge attempt",
        ),
        AnchorCandidate(
            "latent crossing",
            T_LATENT_CROSSING,
            False,
            "branch-resolved entropy-spike crossing; gamma route failed",
        ),
    ]


def select_anchor(rows: list[AnchorCandidate]) -> AnchorCandidate:
    external = [row for row in rows if row.external_bath]
    assert len(external) == 1, external
    selected = external[0]
    assert selected.name == "ramp start"
    assert selected.in_carried_band
    return selected


def main() -> None:
    verify_scale_invariance()
    rows = candidates()
    selected = select_anchor(rows)

    transition = [row for row in rows if row.name.startswith("C_v transition")]
    assert all(not row.external_bath for row in transition)
    assert all(not row.in_carried_band for row in transition)

    print("[0] K04 w6 anchor-epoch selection")
    print("    Metropolis scale invariance: PASS")
    print("    K04 can determine T_anchor/w6, not w6/Lambda by itself.")
    print()

    print("[1] candidate epochs")
    print(
        f"    {'epoch':<20s} {'T/w6':>8s} {'w6/Lambda':>12s} "
        f"{'e_D/Lambda':>12s} {'band':>8s}  reading"
    )
    for row in rows:
        band = "inside" if row.in_carried_band else "outside"
        chosen = " SELECTED" if row == selected else ""
        print(
            f"    {row.name:<20s} {row.t_over_w6:8.3f} "
            f"{row.w6_over_lambda:12.6f} {row.e_debris_over_lambda:12.6f} "
            f"{band:>8s}  {row.reading}{chosen}"
        )
    print()

    w6_gev = selected.w6_over_lambda * LAMBDA_QCD_GEV
    print("[2] selected bridge under current canon")
    print("    premise: the boot bath/ramp-start temperature is Lambda.")
    print(f"    selected anchor: {selected.name}, T_start = {T_RAMP_START:g} w6")
    print(f"    w6/Lambda = 1/6 = {selected.w6_over_lambda:.6f}")
    print(f"    w6 = {w6_gev:.6f} GeV for Lambda = {LAMBDA_QCD_GEV:.3f} GeV")
    print(
        "    debris mass scale: e_D ~= "
        f"{DEBRIS_EXCESS_PER_VERTEX_W6:.2f} w6 = "
        f"{selected.e_debris_over_lambda:.3f} Lambda"
    )
    print()

    print("[3] verdict")
    print("    The discrete choice is conditionally closed: ramp-start is selected")
    print("    by the canon boot-bath premise.  Transition anchoring remains a")
    print("    named alternative postulate, not a K04-derived scale.")
    print("    This does not derive the boot cooling law gamma; the entropy gate")
    print("    for gamma=exp(-alpha0 ln2) remains failed.")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
