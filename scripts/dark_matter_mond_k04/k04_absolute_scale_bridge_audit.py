#!/usr/bin/env python3
"""Audit the two K04 absolute-scale bridges.

The embedded K04/debris runs determine dimensionless lattice quantities:
cycle-count energies in units of w6, defect fractions, and cooling sweeps.
This script checks which parts are fixed by the K04 ledger itself and which
remain external scale identifications.
"""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

ALPHA0 = 1.0 / 137.035999084
LAMBDA_QCD_GEV = 0.332
T_CMB_GEV = 2.348e-13
GAMMA_PAPER = 0.995
T_START = 6.0
T_END = 0.5


def metropolis_accept(delta_e: float, temperature: float) -> float:
    return min(1.0, math.exp(-delta_e / temperature))


def verify_scale_invariance() -> None:
    """Metropolis dynamics cannot determine an absolute energy scale."""

    examples = [
        (1.0, 2.0, 0.5),
        (12.0, 3.0, 7.0),
        (208.3, 6.0, 113.0),
    ]
    for delta_e, temperature, scale in examples:
        base = metropolis_accept(delta_e, temperature)
        scaled = metropolis_accept(scale * delta_e, scale * temperature)
        assert abs(base - scaled) < 1e-15, (delta_e, temperature, scale, base, scaled)


def verify_gamma_source() -> str:
    tex = (ROOT / "crystallize" / "crystallisation.tex").read_text(encoding="utf-8")
    required = [
        "geometric cooling schedule",
        r"\gamma^k",
        r"\gamma = 0.995",
        r"approximately 80\% of uphill moves",
    ]
    missing = [needle for needle in required if needle not in tex]
    assert not missing, missing
    return "crystallize/crystallisation.tex"


def main() -> None:
    verify_scale_invariance()
    source = verify_gamma_source()

    r_paper = math.log(T_START / T_END) / (-math.log(GAMMA_PAPER))
    sweeps_per_efold = 1.0 / (-math.log(GAMMA_PAPER))

    gamma_alpha_bit = math.exp(-ALPHA0 * math.log(2.0))
    r_alpha_bit = math.log(T_START / T_END) / (ALPHA0 * math.log(2.0))

    gamma_28 = math.exp(-1.0 / 28.0)
    gamma_56 = math.exp(-1.0 / 56.0)
    r_28 = 28.0 * math.log(T_START / T_END)
    r_56 = 56.0 * math.log(T_START / T_END)

    w6_floor = 0.05 * LAMBDA_QCD_GEV
    w6_ceiling = 0.27 * LAMBDA_QCD_GEV
    arrhenius_floor = 3.0 * w6_floor / T_CMB_GEV
    arrhenius_ceiling = 3.0 * w6_ceiling / T_CMB_GEV

    print("[0] K04 absolute-scale bridge audit")
    print("    Metropolis scale invariance: PASS")
    print("    If (w4,w6,lambda,T) -> s*(w4,w6,lambda,T),")
    print("    every acceptance probability is unchanged.")
    print()

    print("[1] w6 <-> Lambda")
    print("    K04 fixes sign/order data: w4 > w6 > 0 and lambda >> 1.")
    print("    It does not fix an absolute GeV scale for w6.")
    print(f"    carried band: w6/Lambda in [0.05, 0.27]")
    print(
        "    late thermal barrier exponent 3 w6/T0 spans "
        f"[{arrhenius_floor:.3e}, {arrhenius_ceiling:.3e}]"
    )
    print("    so mobility conclusions are insensitive to the exact band.")
    print()

    print("[2] gamma source and candidates")
    print(f"    paper source: {source}")
    print(f"    gamma_paper = {GAMMA_PAPER:.9f}")
    print(f"    R_paper = ln(12)/[-ln(gamma)] = {r_paper:.3f} sweeps")
    print(f"    sweeps per e-fold = {sweeps_per_efold:.3f}")
    print()
    print("    candidate alpha-bit cooling:")
    print(f"      gamma = exp(-alpha0 ln2) = {gamma_alpha_bit:.9f}")
    print(f"      R = ln(12)/(alpha0 ln2) = {r_alpha_bit:.3f}")
    print(
        "      relative R shift from paper = "
        f"{(r_alpha_bit / r_paper - 1.0):+.3%}"
    )
    print()
    print("    expansion-clock branches:")
    print(f"      28-clock gamma = exp(-1/28) = {gamma_28:.9f}, R = {r_28:.3f}")
    print(f"      56-clock gamma = exp(-1/56) = {gamma_56:.9f}, R = {r_56:.3f}")
    print()

    print("[3] verdict")
    print("    w6 scale: underived; K04 has a scale-invariance obstruction.")
    print("    gamma=0.995: canon-pinned protocol parameter, not a boot law.")
    print("    alpha-bit cooling is a sharp theorem target, not a closure.")


if __name__ == "__main__":
    main()
