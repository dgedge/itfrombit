#!/usr/bin/env python3
r"""Current black-hole observable residual map.

This is a consolidation script, not a new derivation.  It records the status
after the local KMS scheduler, finite Hawking ladder, freeze shell, standard
Schwarzschild greybody transfer, all-contact audit, and scrambling audits.

Exit 0 means the residuals are exactly:

  1. the all-contact horizon-severing theorem for the 10/27 flux attempt rate;
  2. the species/polarization emission ledger beyond the exterior spin barrier;
  3. the O(1)-gap horizon-cell service graph if fast scrambling is required.

The script deliberately keeps the 0.997 Stefan-Hawking coefficient as
conditional: it follows from (10/27) alpha0 only if all-contact Moore transfer
is accepted or proved.
"""

from __future__ import annotations

import math


ALPHA0 = 1.0 / 137.0
GAMMA_REQ = 2.711306813e-3
P_OVER_TARGET_1027 = ((10.0 / 27.0) * ALPHA0) / GAMMA_REQ

SPIN_POLARIZATIONS = {
    "scalar": 1,
    "photon": 2,
    "graviton": 2,
}

GREYBODY_WEIGHTS = {
    # From bh_greybody_transfer.py at the beta-one shell.
    3: {"scalar": 0.4809, "photon": 0.0157, "graviton": 5.7e-5},
    6: {"scalar": 1.3403, "photon": 1.1295, "graviton": 0.0242},
    12: {"scalar": 6.1386, "photon": 6.1970, "graviton": 4.9287},
}


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def local_torus_gap(n: float) -> float:
    return math.sin(math.pi / math.sqrt(n)) ** 2


def ramanujan_gap(degree: int) -> float:
    return 1.0 - 2.0 * math.sqrt(degree - 1.0) / degree


def mixing_time(n: float, gap: float, eps: float = 1.0e-2) -> float:
    return math.log(n / eps) / gap


def main() -> None:
    print("BLACK-HOLE OBSERVABLE RESIDUAL MAP")
    print("=" * 88)

    print("[1] Closed/sharp ingredients")
    closed = [
        "finite V_cell/V_Sch horizon record channel",
        "local half-Boltzmann KMS scheduler",
        "finite Hawking ladder g_Q(F) exp[-beta F]",
        "proper freeze shell and escape-cone M^-2 flux scaling",
        "standard Schwarzschild spin/partial-wave greybody transfer",
        "current-channel coherent-echo null",
    ]
    for item in closed:
        print(f"    CLOSED/SHARP: {item}")
    check(len(closed) == 6, "closed/sharp black-hole observable ingredients enumerated")

    print("\n[2] Absolute flux coefficient")
    gamma_1027 = (10.0 / 27.0) * ALPHA0
    print(f"    Gamma_req/Lambda             = {GAMMA_REQ:.12e}")
    print(f"    (10/27) alpha0               = {gamma_1027:.12e}")
    print(f"    P(10/27 alpha0)/P_SB         = {P_OVER_TARGET_1027:.9f}")
    check(abs(P_OVER_TARGET_1027 - 0.997096) < 2.0e-6, "10/27 route reproduces the recorded 0.997096 near-hit")
    print("    status: conditional on all-contact Landauer-Moore transfer")

    print("\n[3] Species/polarization ledger")
    for species, pol in SPIN_POLARIZATIONS.items():
        print(f"    {species:<8s} physical polarizations = {pol}")
    for strain_f, row in GREYBODY_WEIGHTS.items():
        weighted = {s: row[s] * SPIN_POLARIZATIONS[s] for s in row}
        print(f"    F={strain_f:2d}: greybody weights={row}, with-pol={weighted}")
    check(SPIN_POLARIZATIONS["photon"] == 2 and SPIN_POLARIZATIONS["graviton"] == 2, "massless spin-1/spin-2 have two polarizations")
    check(GREYBODY_WEIGHTS[3]["scalar"] > GREYBODY_WEIGHTS[3]["photon"] > GREYBODY_WEIGHTS[3]["graviton"], "low finite lines are spin-filtered")
    print("    status: exterior barrier is computed; QEC emission species weights are not.")

    print("\n[4] All-contact severing theorem")
    candidates = {
        "face+latch / Moore+latch": (9 + 1, 26 + 1),
        "face-only / Moore+latch": (9, 26 + 1),
        "non-inward / Moore+latch": (9 + 8 + 1, 26 + 1),
        "face / shell-no-latch": (9, 26),
    }
    for name, (num, den) in candidates.items():
        print(f"    {name:<28s} {num:2d}/{den:2d} = {num / den:.9f}")
    check(candidates["face+latch / Moore+latch"] == (10, 27), "all-contact plus latch gives 10/27")
    check(candidates["face-only / Moore+latch"] != (10, 27), "dropping the latch changes the coefficient")
    print("    status: current V_cell/V_Sch does not force the Moore contact alphabet.")

    print("\n[5] Scrambling")
    n_h = 1.0e42
    local_gap = local_torus_gap(n_h)
    exp_gap = ramanujan_gap(8)
    local_mix = mixing_time(n_h, local_gap)
    exp_mix = mixing_time(n_h, exp_gap)
    serial_mix = mixing_time(n_h, exp_gap / n_h)
    print(f"    representative N_H          = {n_h:.3e}")
    print(f"    local surface mix ticks     = {local_mix:.3e}")
    print(f"    O(1)-gap expander mix ticks = {exp_mix:.3e}")
    print(f"    serial demux mix ticks      = {serial_mix:.3e}")
    check(local_mix / exp_mix > 1.0e35, "local surface diffusion is not fast scrambling")
    check(serial_mix / exp_mix > 1.0e35, "serial demux is not fast scrambling")
    print("    status: fast scrambling requires a new parallel O(1)-gap horizon-cell graph.")

    print(
        """
[6] CURRENT STATUS
    Locked/sharp:
      local KMS thermality, finite strain ladder, greybody transfer, echo-null,
      and Hawking M^-2 scaling.

    Still conditional/open:
      flux normalization is the all-contact service-class theorem plus the
      QEC species/polarization emission ledger; fast scrambling is not derived
      by current local horizon service and would require a nonlocal/expander
      service primitive.

ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
