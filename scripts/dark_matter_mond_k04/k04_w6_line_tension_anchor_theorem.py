#!/usr/bin/env python3
r"""K04 w6 anchor epoch: line-tension service theorem.

This audit separates the no-go from the positive bridge.

No-go:
  K04 Metropolis dynamics is exactly scale invariant.  Response observables
  such as the heat-capacity transition locate only T/w6.  They cannot determine
  w6/Lambda_QCD by themselves.

Positive bridge:
  Once the Q3 edge-ledger normal ordering gives w4/w6 = 2, the native K04
  line-service tension is

      mu = w4 + 4 w6 = 6 w6.

  The simulation ramp starts at T_start = 6 w6.  Therefore the ramp start is
  not just an arbitrary high-temperature endpoint: it is exactly the native
  line-service quantum.  If the strong-sector anchor Lambda_QCD prices that
  line-service quantum, then

      mu = T_start = Lambda_QCD,     w6 = Lambda_QCD / 6.

  Transition anchoring remains a different physical postulate.  It assigns the
  heat-capacity response temperature to Lambda_QCD and consequently makes the
  native line-service quantum mu about two Lambda_QCD, introducing an extra
  dimensionless offset not supplied by the K04 service ledger.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction


LAMBDA_QCD_GEV = 0.332
W4_OVER_W6 = Fraction(2, 1)
MU_OVER_W6 = W4_OVER_W6 + 4
SIGMA_WALL_OVER_W6 = 2 * W4_OVER_W6 + 12
T_START_OVER_W6 = Fraction(6, 1)
DEBRIS_EXCESS_PER_VERTEX_W6 = 2.17


@dataclass(frozen=True)
class Candidate:
    name: str
    t_over_w6: float
    role: str

    @property
    def w6_over_lambda_if_temperature_anchored(self) -> float:
        return 1.0 / self.t_over_w6

    @property
    def mu_over_lambda_if_temperature_anchored(self) -> float:
        return float(MU_OVER_W6) / self.t_over_w6

    @property
    def sigma_over_lambda_if_temperature_anchored(self) -> float:
        return float(SIGMA_WALL_OVER_W6) / self.t_over_w6

    @property
    def debris_over_lambda_if_temperature_anchored(self) -> float:
        return DEBRIS_EXCESS_PER_VERTEX_W6 / self.t_over_w6


def metropolis_accept(delta_e: float, temperature: float) -> float:
    return min(1.0, math.exp(-delta_e / temperature))


def assert_scale_invariance() -> None:
    for delta_e, temperature, scale in [
        (1.0, 6.0, 9.0),
        (4.0, 2.8, 0.41),
        (16.0, 3.15, 23.0),
    ]:
        base = metropolis_accept(delta_e, temperature)
        scaled = metropolis_accept(scale * delta_e, scale * temperature)
        assert abs(base - scaled) < 1e-15


def candidates() -> list[Candidate]:
    return [
        Candidate("ramp start", 6.0, "external bath endpoint and line-service quantum"),
        Candidate("C_v transition low", 2.8, "emergent response temperature"),
        Candidate("C_v transition high", 3.0, "emergent response temperature"),
        Candidate("arrest entry", 2.99, "emergent arrest response"),
        Candidate("latent crossing", 3.15, "entropy-spike crossing; gamma route failed"),
    ]


def select_by_line_service(rows: list[Candidate]) -> Candidate:
    selected = [
        row
        for row in rows
        if abs(row.t_over_w6 - float(MU_OVER_W6)) < 1e-12
    ]
    assert len(selected) == 1
    assert selected[0].name == "ramp start"
    return selected[0]


def main() -> None:
    assert_scale_invariance()
    rows = candidates()
    selected = select_by_line_service(rows)

    assert W4_OVER_W6 == 2
    assert MU_OVER_W6 == T_START_OVER_W6
    assert SIGMA_WALL_OVER_W6 == 16
    assert abs(selected.mu_over_lambda_if_temperature_anchored - 1.0) < 1e-12
    transition_rows = [row for row in rows if "transition" in row.name]
    assert all(row.mu_over_lambda_if_temperature_anchored > 1.8 for row in transition_rows)

    w6_gev = selected.w6_over_lambda_if_temperature_anchored * LAMBDA_QCD_GEV
    debris_gev = selected.debris_over_lambda_if_temperature_anchored * LAMBDA_QCD_GEV
    sigma_gev = selected.sigma_over_lambda_if_temperature_anchored * LAMBDA_QCD_GEV

    print("[0] K04 w6 line-tension anchor theorem")
    print("    Metropolis response scale invariance: PASS")
    print("    Therefore K04 transition data alone cannot set w6/Lambda_QCD.")
    print()

    print("[1] Native K04 service quantum")
    print(f"    edge-ledger input: w4/w6 = {W4_OVER_W6}")
    print(f"    line-service tension: mu = w4 + 4 w6 = {MU_OVER_W6} w6")
    print(f"    frustrated wall cut tension: sigma = 2 w4 + 12 w6 = {SIGMA_WALL_OVER_W6} w6")
    print(f"    ramp start: T_start = {T_START_OVER_W6} w6")
    print("    identity: T_start = mu")
    print()

    print("[2] Candidate temperature anchors if each were set equal to Lambda_QCD")
    print(
        f"    {'candidate':<22s} {'T/w6':>7s} {'w6/Lambda':>11s} "
        f"{'mu/Lambda':>10s} {'sigma/Lambda':>13s} {'e_D/Lambda':>12s}  role"
    )
    for row in rows:
        marker = " SELECTED" if row == selected else ""
        print(
            f"    {row.name:<22s} {row.t_over_w6:7.3f} "
            f"{row.w6_over_lambda_if_temperature_anchored:11.6f} "
            f"{row.mu_over_lambda_if_temperature_anchored:10.6f} "
            f"{row.sigma_over_lambda_if_temperature_anchored:13.6f} "
            f"{row.debris_over_lambda_if_temperature_anchored:12.6f}  "
            f"{row.role}{marker}"
        )
    print()

    print("[3] Absolute bridge under the line-service reading")
    print("    premise discharged to: the strong-sector anchor prices the native")
    print("    K04 line-service quantum, mu = Lambda_QCD.")
    print(f"    w6/Lambda_QCD = 1/6 = {selected.w6_over_lambda_if_temperature_anchored:.6f}")
    print(f"    w6 = {w6_gev:.6f} GeV for Lambda_QCD = {LAMBDA_QCD_GEV:.3f} GeV")
    print(f"    e_D = {DEBRIS_EXCESS_PER_VERTEX_W6:.2f} w6 = {debris_gev:.6f} GeV")
    print(f"    sigma_wall = 16 w6 = {sigma_gev:.6f} GeV per cut cell")
    print()

    print("[4] Verdict")
    print("    Ramp-start anchoring is selected once the anchor is a line-service")
    print("    tension.  Transition anchoring is not selected by the present K04")
    print("    service algebra; it is a new postulate that the response temperature,")
    print("    rather than the native service tension, is the object priced by Lambda_QCD.")
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
