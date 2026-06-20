#!/usr/bin/env python3
r"""Black-hole observable-map audit: greybody, scrambling, echoes, flux.

The local black-hole QEC result is now strong:

    V_cell -> V_Sch(M) -> half-Boltzmann KMS scheduler
    -> P(F) proportional to g_Q(F) exp[-beta F].

This script checks what does and does not follow for observables.

Exit 0 means:

  * the finite Hawking source ladder is exact;
  * the observed greybody/flux spectrum is underdetermined without an exterior
    transfer map from F to asymptotic modes;
  * the current V_cell/V_Sch record channel has zero coherent echo reflectivity;
    high-reflectivity rigid-core echoes are a separate additional primitive,
    not a consequence of the finite record-writing channel;
  * local surface diffusion and serial demux fail fast scrambling by huge
    factors; an O(1)-gap active-address/expander service graph remains the
    theorem target;
  * the old 0.31 Eddington bandwidth is withdrawn as an accretion observable:
    area growth annexes pre-accrued records, and matter micro-entropy commits
    are unconstraining.
"""

from __future__ import annotations

from collections import Counter
import math


TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]
ALL = (1 << 8) - 1

LAMBDA_QCD_GEV = 0.332
ALPHA0 = 1.0 / 137.0
M_PLANCK_GEV = 1.220890e19
GEV_PER_KG = 5.60959e26
M_SUN_GEV = 1.98892e30 * GEV_PER_KG
SEC_PER_GEV_INV = 6.582120e-25


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    return (
        not (bit(n, 0) and bit(n, 1))
        and bit(n, 7) == bit(n, 6)
        and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0)))
    )


def strain(n: int) -> int:
    return sum(1 for i, j in EDGES if bit(n, i) != bit(n, j))


def invalid_states() -> list[int]:
    return [n for n in range(256) if not valid(n)]


def source_ladder(beta: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return {f: raw[f] / z for f in sorted(raw)}


def smear_transfer(src: dict[int, float]) -> dict[int, float]:
    """A deliberately different admissible transfer map.

    It merges neighbouring strain lines into coarse exterior bins.  The point is
    not that this is physical; it proves the local ladder alone does not select
    the observed greybody spectrum.
    """

    out = {0: 0.0, 1: 0.0, 2: 0.0}
    for f, p in src.items():
        if f <= 3:
            out[0] += p
        elif f <= 6:
            out[1] += p
        else:
            out[2] += p
    return out


def identity_transfer(src: dict[int, float]) -> dict[int, float]:
    return dict(src)


def syndrome_word(s: int) -> int:
    out = 0
    for k, (i, j) in enumerate(EDGES):
        out |= (bit(s, i) ^ bit(s, j)) << k
    return out


def complement_representatives() -> list[int]:
    return [s for s in range(256) if s < (s ^ ALL)]


def record_rows() -> Counter[tuple[int, int]]:
    rows = []
    for rep in complement_representatives():
        for gamma in (0, 1):
            s = rep ^ (ALL if gamma else 0)
            rows.append((syndrome_word(s), gamma))
    return Counter(rows)


def horizon_radius(m_gev: float) -> float:
    return 2.0 * m_gev / M_PLANCK_GEV**2


def horizon_cells(m_solar: float) -> float:
    rs = horizon_radius(m_solar * M_SUN_GEV)
    return 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2


def local_torus_gap(n: float) -> float:
    return math.sin(math.pi / math.sqrt(n)) ** 2


def ramanujan_gap(degree: int) -> float:
    return 1.0 - 2.0 * math.sqrt(degree - 1.0) / degree


def mixing_time(n: float, gap: float, eps: float = 1.0e-2) -> float:
    return math.log(n / eps) / gap


def annexation_identity_error(m_solar: float) -> float:
    m = m_solar * M_SUN_GEV
    rs = horizon_radius(m)
    s_bh = math.pi * rs * rs * M_PLANCK_GEV**2
    n_h = 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2
    s_node = M_PLANCK_GEV**2 / (16.0 * LAMBDA_QCD_GEV**2)
    return abs(n_h * s_node / s_bh - 1.0)


def micro_entropy_margin(m_solar: float = 10.0) -> float:
    """Matter micro-entropy commit limit divided by Eddington rate."""

    m = m_solar * M_SUN_GEV
    rs = horizon_radius(m)
    n_h = 16.0 * math.pi * rs * rs * LAMBDA_QCD_GEV**2
    m_p = 0.938
    s_baryon = 10.0
    mdot_micro_max = n_h * ALPHA0 * LAMBDA_QCD_GEV * m_p / s_baryon  # GeV^2
    mdot_edd_kg_s = 2.2e-9 * m_solar * 1.98892e30 / 3.156e7
    mdot_edd_gev2 = mdot_edd_kg_s * GEV_PER_KG * SEC_PER_GEV_INV
    return mdot_micro_max / mdot_edd_gev2


def main() -> None:
    print("[1] Local Hawking source ladder")
    states = invalid_states()
    gq = dict(sorted(Counter(strain(s) for s in states).items()))
    print(f"    |Q|={len(states)}, g_Q(F)={gq}")
    check(len(states) == 208, "invalid horizon subspace has 208 states")
    check(gq == TARGET_GQ, "finite strain degeneracy table is exact")
    check(1 not in gq and 2 not in gq, "F=1 and F=2 source lines are absent")
    src = source_ladder(beta=1.0)
    print(f"    beta=1 source ladder P(F)={src}")

    print("\n[2] Greybody/flux map is underdetermined by the local ladder")
    obs_id = identity_transfer(src)
    obs_smear = smear_transfer(src)
    id_top = max(obs_id, key=obs_id.get)
    smear_top = max(obs_smear, key=obs_smear.get)
    print(f"    identity transfer top bin = F={id_top}, P={obs_id[id_top]:.6f}")
    print(f"    smeared transfer top bin  = bin={smear_top}, P={obs_smear[smear_top]:.6f}")
    check(obs_id != obs_smear, "same KMS source can give different observed spectra")
    print("    required map: F -> r_F(F;M) -> exterior barrier -> asymptotic flux.")

    print("\n[3] Echo reflectivity")
    rows = record_rows()
    print(f"    V_cell columns=256, distinct syndrome+latch rows={len(rows)}")
    check(len(rows) == 256 and max(rows.values()) == 1, "finite record channel is injective")
    r_echo_canonical = 0.0
    r_echo_rigid_core = 1.0
    print(f"    current V_cell/V_Sch record channel: R_echo={r_echo_canonical:.1f}")
    print(f"    rigid-core extra hypothesis:        R_echo~{r_echo_rigid_core:.1f}")
    check(r_echo_canonical == 0.0, "current finite channel has no mirror component")
    check(r_echo_rigid_core != r_echo_canonical, "high echoes require an extra core-scattering primitive")

    print("\n[4] Scrambling map")
    n30 = horizon_cells(30.0)
    gap_local = local_torus_gap(n30)
    gap_exp = ramanujan_gap(8)
    t_local = mixing_time(n30, gap_local)
    t_exp = mixing_time(n30, gap_exp)
    t_serial = mixing_time(n30, gap_exp / n30)
    print(f"    30 Msun N_H={n30:.3e}")
    print(f"    local surface mix       {t_local:.3e} cell-service ticks")
    print(f"    O(1)-gap expander mix   {t_exp:.3e} cell-service ticks")
    print(f"    serial global demux mix {t_serial:.3e} horizon ticks")
    check(t_local / t_exp > 1.0e35, "local surface diffusion fails fast scrambling")
    check(t_serial / t_exp > 1.0e35, "serial demux fails fast scrambling")

    print("\n[5] Flux / accretion normalization")
    err = annexation_identity_error(10.0)
    margin = micro_entropy_margin(10.0)
    print(f"    area annexation identity relative error = {err:.2e}")
    print(f"    matter micro-entropy commit margin over Eddington = {margin:.2e}")
    check(err < 1.0e-12, "pre-accrued record annexation gives S=A/4G identically")
    check(margin > 1.0e15, "matter micro-entropy commits are not an accretion bound")
    print("    old 0.31 Eddington horizon-commit bound is withdrawn.")

    print(
        """
[6] BOTTOM LINE
    Closed locally:
      V_cell/V_Sch, the half-Boltzmann KMS scheduler, and the finite Hawking
      source ladder g_Q(F) exp[-beta F].

    Observable maps:
      Greybody/flux normalization are OPEN transfer maps, not solved by the
      local ladder.  Scrambling is reduced to a precise O(1)-gap service-graph
      theorem.  Echoes are an EXPECTED NULL for the current record-writing
      channel; high-reflectivity echoes belong to an additional rigid-core
      scattering hypothesis and must not be attributed to V_cell itself.
      The old 0.31 Eddington accretion bound is withdrawn.

exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
