#!/usr/bin/env python3
"""
TCH first-generation quark node suppression probe.

Question:
  Can the u,d -> (u,d)/3 correction found by the Yukawa hierarchy probe be
  derived from framework structure rather than searched?

What this script proves:
  A rank-one local quark-color node has normalized color trace 1/3:

      tau(P_color) = Tr_color(P_color) / N_c = 1/3.

  Equivalently, the Z3 color-orbit twirl of one color projector is I_3/3.
  By Schur's lemma the full SU(3) twirl of a fundamental rank-one projector
  gives the same scalar operator I_3/3. This is a trace/twirl theorem, not a
  claim that a single quark color vector is an SU(3) singlet.

  Important audit status: this is a normalized-rank identity. Any rank-one
  projector in a 3D color space has trace 1 and therefore normalized trace
  1/3. The TCH-specific work would be to derive why the first-generation
  quark mass node has that rank-one block. This file does not derive that.

Framework inputs used:
  1. The [8,4,4] code has one lepton column plus three quark-color columns.
  2. The Pati-Salam stabilizer identity equates a lepton-column defect with a
     tri-color quark composite on physical states.
  3. ANCHOR item 96 already places quark Koide phases in a color-coherent
     defect-fraction structure: delta_u = 2/27, delta_d = 1/9.

Conditional step:
  If the first-generation quark node is the rank-one color footprint of the
  TCH mass operator before the color trace, then its scalar mass footprint is
  suppressed by 1/3 relative to the tri-color Pati-Salam support.

What this script does not prove:
  It does not construct the full TCH mass operator. The load-bearing remaining
  lemma is that the first-generation u,d nodes are exactly those rank-one color
  footprints, and that the chosen down-sector R does not already include the
  same trace factor. The documented 256-state walk operator was checked in
  tch_256_walk_mass_rank_audit.py and does not produce this rank pattern.
"""

from __future__ import annotations

import math
from fractions import Fraction


N_C = 3

MASSES_GEV = {
    "up": {"u": 2.16e-3, "c": 1.27, "t": 172.76},
    "down": {"d": 4.67e-3, "s": 93.4e-3, "b": 4.18},
    "lepton": {"e": 0.510998950e-3, "mu": 105.6583755e-3, "tau": 1.77686},
}

SECTOR_ORDER = {
    "up": ["u", "c", "t"],
    "down": ["d", "s", "b"],
    "lepton": ["e", "mu", "tau"],
}


def hd(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def zero_matrix(n: int) -> list[list[Fraction]]:
    return [[Fraction(0) for _ in range(n)] for _ in range(n)]


def eye_matrix(n: int) -> list[list[Fraction]]:
    out = zero_matrix(n)
    for i in range(n):
        out[i][i] = Fraction(1)
    return out


def projector(index: int, n: int = N_C) -> list[list[Fraction]]:
    out = zero_matrix(n)
    out[index][index] = Fraction(1)
    return out


def coherent_orbit_projector(n: int = N_C) -> list[list[Fraction]]:
    return [[Fraction(1, n) for _ in range(n)] for _ in range(n)]


def cycle_conjugate(matrix: list[list[Fraction]], shift: int) -> list[list[Fraction]]:
    n = len(matrix)
    out = zero_matrix(n)
    for i in range(n):
        for j in range(n):
            out[(i + shift) % n][(j + shift) % n] = matrix[i][j]
    return out


def add_matrices(mats: list[list[list[Fraction]]]) -> list[list[Fraction]]:
    n = len(mats[0])
    out = zero_matrix(n)
    for mat in mats:
        for i in range(n):
            for j in range(n):
                out[i][j] += mat[i][j]
    return out


def scale_matrix(matrix: list[list[Fraction]], scale: Fraction) -> list[list[Fraction]]:
    return [[scale * value for value in row] for row in matrix]


def trace(matrix: list[list[Fraction]]) -> Fraction:
    return sum(matrix[i][i] for i in range(len(matrix)))


def rank_diagonal_projector(matrix: list[list[Fraction]]) -> int:
    """Rank for the diagonal projectors used in the color-trace proof."""
    return sum(1 for i in range(len(matrix)) if matrix[i][i] != 0)


def format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def color_trace_theorem() -> Fraction:
    hd("A. Exact Normalized-Rank Identity")

    p_r = projector(0)
    p_g = projector(1)
    p_b = projector(2)
    identity = eye_matrix(N_C)
    tri_color = add_matrices([p_r, p_g, p_b])

    normalized_single_trace = trace(p_r) / N_C
    normalized_tri_trace = trace(tri_color) / N_C
    single_over_tri = normalized_single_trace / normalized_tri_trace

    z3_twirl = scale_matrix(
        add_matrices([cycle_conjugate(p_r, shift) for shift in range(N_C)]),
        Fraction(1, N_C),
    )
    expected_twirl = scale_matrix(identity, Fraction(1, N_C))
    coherent = coherent_orbit_projector()

    assert tri_color == identity
    assert z3_twirl == expected_twirl
    assert trace(coherent) == 1
    assert all(coherent[i][i] == Fraction(1, N_C) for i in range(N_C))
    assert single_over_tri == Fraction(1, 3)

    print(f"  color dimension N_c                         = {N_C}")
    print(f"  rank(P_single_color)                        = {rank_diagonal_projector(p_r)}")
    print(f"  rank(P_tri_color = P_r+P_g+P_b)             = {rank_diagonal_projector(tri_color)}")
    print(
        "  normalized trace tau(P_single_color)        = "
        f"{format_fraction(normalized_single_trace)}"
    )
    print(
        "  normalized trace tau(P_tri_color)           = "
        f"{format_fraction(normalized_tri_trace)}"
    )
    print(
        "  relative first-node color footprint         = "
        f"{format_fraction(single_over_tri)}"
    )
    print("  Z3 orbit twirl of P_single_color            = I_3/3 exactly")
    print("  coherent color-orbit projector diagonal     = (1/3, 1/3, 1/3)")
    print(
        "  full SU(3) reading                          = same scalar by Schur twirl, "
        "not a one-quark singlet vector"
    )
    print("  audit status                                = rank-one input, not node derivation")
    return single_over_tri


def koide_shape_values(r_value: float, delta: float) -> list[float]:
    return [
        (1 + r_value * math.cos(delta + 2 * math.pi * index / 3)) ** 2
        for index in range(3)
    ]


def koide_predicted_logs(
    sector: str,
    r_value: float,
    delta: float,
    assignment: tuple[int, int, int] = (1, 2, 0),
) -> dict[str, float]:
    names = SECTOR_ORDER[sector]
    values = koide_shape_values(r_value, delta)
    assigned = {name: math.log10(values[index]) for name, index in zip(names, assignment)}
    heavy = names[-1]
    return {name: assigned[name] - assigned[heavy] for name in names}


def corrected_predicted_logs(
    base: dict[str, float],
    correction: dict[str, float],
) -> dict[str, float]:
    return {
        particle: log_value + math.log10(correction.get(particle, 1.0))
        for particle, log_value in base.items()
    }


def observed_log_ratios(sector: str) -> dict[str, float]:
    names = SECTOR_ORDER[sector]
    heavy = names[-1]
    return {
        name: math.log10(MASSES_GEV[sector][name] / MASSES_GEV[sector][heavy])
        for name in names
    }


def ratio_errors(
    observed: dict[str, float],
    predicted: dict[str, float],
) -> dict[str, tuple[float, float]]:
    out = {}
    for name, obs in observed.items():
        if abs(obs) < 1e-12:
            continue
        err = predicted[name] - obs
        out[name] = (err, abs(err) / abs(obs))
    return out


def max_rel_for_sector(sector: str, predicted: dict[str, float]) -> float:
    return max(rel for _, rel in ratio_errors(observed_log_ratios(sector), predicted).values())


def print_ratio_table(title: str, sector: str, predicted: dict[str, float]) -> float:
    observed = observed_log_ratios(sector)
    errors = ratio_errors(observed, predicted)
    print(f"\n{title} / {sector}")
    print("  particle  obs log10(m/m_heavy)  pred       abs err    rel err")
    for name in SECTOR_ORDER[sector][:-1]:
        err, rel = errors[name]
        print(
            f"  {name:<8s} "
            f"{observed[name]:>10.4f}             "
            f"{predicted[name]:>8.4f}  "
            f"{abs(err):>8.4f}  "
            f"{rel:>7.1%}"
        )
    max_rel = max(rel for _, rel in errors.values())
    verdict = "PASS" if max_rel <= 0.10 else "FAIL"
    print(f"  -> {verdict}: max relative log-ratio error = {max_rel:.1%}")
    return max_rel


def hierarchy_consequence(trace_factor: Fraction) -> None:
    hd("B. Hierarchy Consequence With Inserted Trace Factor")
    factor = float(trace_factor)

    lepton = koide_predicted_logs("lepton", math.sqrt(2), 2 / 9)
    up_base = koide_predicted_logs("up", math.sqrt(3), 2 / 27)
    up_trace = corrected_predicted_logs(up_base, {"u": factor})
    down_two_path_base = koide_predicted_logs("down", math.sqrt(2), 1 / 9)
    down_two_path_trace = corrected_predicted_logs(down_two_path_base, {"d": factor})

    print(
        "  No denominator search is run here, but the factor is applied by choosing "
        f"the first-generation quark nodes to carry rank one: N_c^-1 = {format_fraction(trace_factor)}."
    )
    lepton_err = print_ratio_table("charged lepton: R=sqrt2, delta=2/9", "lepton", lepton)
    up_no_trace_err = print_ratio_table("up: R=sqrt3, delta=2/27, no node trace", "up", up_base)
    up_trace_err = print_ratio_table(
        "up: R=sqrt3, delta=2/27, u-node color trace 1/3",
        "up",
        up_trace,
    )
    down_no_trace_err = print_ratio_table(
        "down branch A: R=sqrt2, delta=1/9, no node trace",
        "down",
        down_two_path_base,
    )
    down_trace_err = print_ratio_table(
        "down branch A: R=sqrt2, delta=1/9, d-node color trace 1/3",
        "down",
        down_two_path_trace,
    )

    combined = max(lepton_err, up_trace_err, down_trace_err)
    print(
        "\n  combined score on branch A after derived trace factor: "
        f"{combined:.1%} max relative log-ratio error"
    )
    print(
        "  branch A passes the 10% shape bar only if the TCH mass-node lemma is accepted."
    )
    print(
        "  improvement from trace: up "
        f"{up_no_trace_err:.1%} -> {up_trace_err:.1%}, down "
        f"{down_no_trace_err:.1%} -> {down_trace_err:.1%}"
    )


def down_branch_audit(trace_factor: Fraction) -> None:
    hd("C. Down-Sector Double-Counting Audit")
    print(
        "  The 1/3 trace cannot be blindly multiplied into every historical down "
        "ansatz. If an R value was fitted after color dressing, it may already "
        "contain part of this normalization."
    )
    factor = float(trace_factor)
    options = [
        ("two-path discrete baseline", math.sqrt(2)),
        ("Part 01 near-fit sqrt(12/5)", math.sqrt(12 / 5)),
        ("three-color path sqrt3", math.sqrt(3)),
    ]
    print("\n  down ansatz                    no trace     with d-node 1/3")
    for label, r_value in options:
        base = koide_predicted_logs("down", r_value, 1 / 9)
        traced = corrected_predicted_logs(base, {"d": factor})
        print(
            f"  {label:<28s} "
            f"{max_rel_for_sector('down', base):>8.1%}   "
            f"{max_rel_for_sector('down', traced):>12.1%}"
        )
    print(
        "\n  Reading: the color trace closes the down hierarchy only on the "
        "unfitted two-path R=sqrt2 branch. Keeping the Part 01 fitted "
        "R~sqrt(12/5) and also applying 1/3 overcorrects d. Since the no-trace "
        "R~sqrt(12/5) branch is much better, the down sector does not provide "
        "positive evidence for a separate trace factor."
    )


def verdict() -> None:
    hd("D. Verdict")
    print(
        "  Proved: any rank-one color block in a 3D color space has normalized "
        "trace 1/3. This is exact but generic."
    )
    print(
        "  Not proved: applying it specifically to first-generation u,d. The "
        "documented 256-state walk operator does not produce C_light rank-one "
        "or C_middle=C_heavy=I_3."
    )
    print(
        "  Consequence: with that extra lemma and the two-path down branch, the "
        "charged hierarchy shape can be made to pass the 10% bar. Without that "
        "lemma it remains an ansatz selected by the fit."
    )
    print(
        "  Down-sector warning: the cleaner R=sqrt(12/5) no-trace branch beats "
        "the traced branch, so the trace is not supported there."
    )


def main() -> None:
    trace_factor = color_trace_theorem()
    hierarchy_consequence(trace_factor)
    down_branch_audit(trace_factor)
    verdict()


if __name__ == "__main__":
    main()
