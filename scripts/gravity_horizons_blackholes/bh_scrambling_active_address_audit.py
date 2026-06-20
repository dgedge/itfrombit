#!/usr/bin/env python3
r"""Black-hole fast-scrambling graph: active address is not yet an expander.

Target
------
The local black-hole sector now has a finite horizon channel, a local
Schwarzschild/KMS scheduler, and an active-address demux principle from the
register handoff work.  The remaining scrambling question is whether those
ingredients already imply black-hole fast scrambling,

    t_* ~ beta log S,

or whether an additional nonlocal horizon-cell service graph is still required.

Result
------
The existing canon does not yet derive fast scrambling.  It derives only:

  * in-cell active-address repair: a register service operation chooses the
    tracked active syndrome address inside a cell;
  * local shell composition: V_Sch(M) is a direct sum over orthogonal horizon
    cells;
  * local one-bit KMS jumps on the 208-state invalid horizon subspace.

None of these supplies a nonlocal adjacency on the set of horizon-cell labels.
The graph alternatives are decisive:

  * local two-dimensional surface diffusion has gap O(1/N_H) and mixes in
    O(N_H log N_H);
  * a globally serial active-address queue also has gap O(1/N_H);
  * an in-cell demux has O(1) internal gap but support radius zero on the
    horizon surface;
  * only a parallel O(1)-gap expander/active-address graph on horizon-cell
    labels gives O(log N_H) mixing.

Therefore the sharp theorem target is not "derive active demux" again; that
exists.  It is: derive a parallel O(1)-gap service-current graph on horizon-cell
labels from the localized-mass QEC scheduler.
"""

from __future__ import annotations

import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LAMBDA_QCD_GEV = 0.332
M_PLANCK_GEV = 1.220890e19
GEV_PER_KG = 5.60959e26
M_SUN_GEV = 1.98892e30 * GEV_PER_KG


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def horizon_radius(m_gev: float) -> float:
    return 2.0 * m_gev / M_PLANCK_GEV**2


def horizon_cells(m_solar: float) -> float:
    m_gev = m_solar * M_SUN_GEV
    rs = horizon_radius(m_gev)
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2


def local_torus_gap(n: float) -> float:
    return math.sin(math.pi / math.sqrt(n)) ** 2


def ramanujan_gap(degree: int) -> float:
    return 1.0 - 2.0 * math.sqrt(degree - 1.0) / degree


def mixing_time(n: float, gap: float, eps: float = 1.0e-2) -> float:
    return math.log(n / eps) / gap


def source_contains(relpath: str, *needles: str) -> bool:
    text = (ROOT / relpath).read_text(encoding="utf-8")
    return all(needle in text for needle in needles)


def main() -> None:
    print("BLACK-HOLE SCRAMBLING ACTIVE-ADDRESS AUDIT")
    print("=" * 86)

    print("\n[1] What existing canon actually supplies")
    check(
        source_contains("python_code/register_handoff_demux_decision_audit.py", "active-address demux"),
        "register handoff supplies active-address demux inside the service ledger",
    )
    check(
        source_contains(
            "python_code/bh_schwarzschild_channel_derivation.py",
            "V_Sch(M) is the direct-sum isometry",
            "Orthogonal shell labels",
        ),
        "V_Sch is a direct sum over orthogonal shell-cell labels",
    )
    check(
        source_contains(
            "python_code/bh_kms_scheduler_derivation.py",
            "symmetric one-bit QEC service",
            "thermal/GNS",
        ),
        "local KMS scheduler acts on one-bit horizon-register neighbours",
    )
    print("    -> these are local/register statements; none is an inter-cell expander.")

    print("\n[2] Four graph readings")
    exp_gap = ramanujan_gap(8)
    print(f"    8-regular expander gap lower scale = {exp_gap:.6f}")
    check(0.3 < exp_gap < 0.4, "constant-degree expander has O(1) gap")

    for m_solar in (3.0, 30.0, 4.3e6):
        n = horizon_cells(m_solar)
        local_gap = local_torus_gap(n)
        serial_gap = exp_gap / n
        local_mix = mixing_time(n, local_gap)
        serial_mix = mixing_time(n, serial_gap)
        exp_mix = mixing_time(n, exp_gap)
        print(f"\n    M={m_solar:g} Msun: N_H={n:.3e}, log N_H={math.log(n):.2f}")
        print(f"      local surface diffusion gap={local_gap:.3e}, mix={local_mix:.3e}")
        print(f"      in-cell demux internal gap=O(1), horizon support radius=0")
        print(f"      serial active-address gap={serial_gap:.3e}, mix={serial_mix:.3e}")
        print(f"      parallel expander gap={exp_gap:.3e}, mix={exp_mix:.3e}")
        check(local_mix / exp_mix > 1.0e35, "local surface diffusion fails fast scrambling")
        check(serial_mix / exp_mix > 1.0e35, "serial active-address service fails fast scrambling")

    print("\n[3] Decision table")
    rows = [
        ("local surface graph", "derived geometry", "FAIL", "gap O(1/N_H)"),
        ("in-cell active demux", "derived register recovery", "FAIL", "no horizon-cell spreading"),
        ("serial global demux", "possible queue reading", "FAIL", "gap O(1/N_H)"),
        ("parallel expander on cell labels", "not derived", "PASS if supplied", "gap O(1)"),
    ]
    for name, provenance, verdict, reason in rows:
        print(f"    {name:<32s} {provenance:<28s} {verdict:<16s} {reason}")

    print(
        r"""
[4] VERDICT
  Fast scrambling is NOT closed by the current active-demux result.

  The active demux selected by the cosmological-constant and R4 ledgers is an
  in-cell service-address theorem: it says the decoder services the tracked
  active register address rather than a blind channel.  It does not define an
  O(1)-gap graph on the N_H horizon-cell labels.  Since V_Sch is currently a
  direct sum over orthogonal shell cells, the local channel factorises until a
  separate surface service-current graph is supplied.

  Closed/refuted:
    local surface diffusion and serial global demux cannot give fast scrambling.

  Remaining exact theorem:
    derive a parallel O(1)-gap active-address/expander graph on horizon-cell
    labels from the localized-mass QEC scheduler, or accept that the present
    horizon-QEC channel does not by itself prove fast scrambling.
ALL ASSERTIONS PASSED -- scrambling graph remains a precise expander theorem target."""
    )


if __name__ == "__main__":
    main()
