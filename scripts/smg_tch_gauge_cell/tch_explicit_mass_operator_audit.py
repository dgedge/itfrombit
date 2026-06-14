#!/usr/bin/env python3
"""
Explicit color-resolved TCH Koide mass-operator audit.

Purpose:
  Construct the smallest explicit mass operator that can test the proposed
  first-generation quark node suppression.

Operator built here:
  K(R,delta) = I + (R/2) * (exp(i delta) S + exp(-i delta) S^dagger)

  K is the documented three-generation circulant amplitude operator. Its
  eigenvalues are

      lambda_n = 1 + R cos(delta + 2 pi n / 3), n = 0,1,2.

  The mass operator is K^2. The physical ordering used by the existing
  hierarchy probe is n=0 heavy, n=1 light, n=2 middle.

Color-resolved candidate:
  M = sum_n lambda_n^2 P_n (x) C_n

  with C_0 = C_2 = I_color and, for the proposed first-generation quark node,
  C_1 = rho_color where Tr(rho_color) = 1. The normalized color trace is

      tau_color(C_1) = Tr(C_1) / 3 = 1/3.

Therefore partial color tracing gives m_light -> m_light / 3. This is an
operator identity conditional on inserting C_1 as a rank-one block. It is not
a derivation of C_1; placing the rank-one block on light/middle/heavy is an
input unless a separate 256-state Yukawa operator produces that pattern.

What this can decide:
  * If the quark mass operator already supplies this rank-one light-node color
    block, the 1/3 suppression is real at operator level.
  * If a sector's fitted R already contains the color normalization, applying
    this trace again double-counts it. The down-sector audit exposes exactly
    that branch ambiguity.

What this cannot decide alone:
  The repository does not contain a full color-resolved Yukawa mass operator
  derived from the 256-state walk. The CKM paper gives a color-projected weak
  transition operator, not the Yukawa mass operator. So this script constructs
  the minimal candidate operator and reports the remaining canonical lemma.

Follow-up:
  tch_256_walk_mass_rank_audit.py implements the documented 256-state walk
  operator directly and finds that this candidate rank pattern is not present
  in that operator. Treat this file as the conditional ansatz audit, not the
  canonical 256-state derivation.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass


N_GEN = 3
N_COLOR = 3
LIGHT_INDEX = 1
MIDDLE_INDEX = 2
HEAVY_INDEX = 0
NODE_LABELS = {
    HEAVY_INDEX: "heavy n=0",
    LIGHT_INDEX: "light n=1",
    MIDDLE_INDEX: "middle n=2",
}

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


Matrix = list[list[complex]]


def hd(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def zeros(n: int, m: int) -> Matrix:
    return [[0j for _ in range(m)] for _ in range(n)]


def identity(n: int) -> Matrix:
    out = zeros(n, n)
    for i in range(n):
        out[i][i] = 1.0 + 0j
    return out


def matadd(a: Matrix, b: Matrix) -> Matrix:
    return [[a[i][j] + b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def matscale(z: complex, a: Matrix) -> Matrix:
    return [[z * value for value in row] for row in a]


def matmul(a: Matrix, b: Matrix) -> Matrix:
    n = len(a)
    p = len(b)
    m = len(b[0])
    out = zeros(n, m)
    for i in range(n):
        for k in range(p):
            if abs(a[i][k]) < 1e-15:
                continue
            aik = a[i][k]
            for j in range(m):
                out[i][j] += aik * b[k][j]
    return out


def kron(a: Matrix, b: Matrix) -> Matrix:
    out = zeros(len(a) * len(b), len(a[0]) * len(b[0]))
    for i in range(len(a)):
        for j in range(len(a[0])):
            for k in range(len(b)):
                for l in range(len(b[0])):
                    out[i * len(b) + k][j * len(b[0]) + l] = a[i][j] * b[k][l]
    return out


def trace(a: Matrix) -> complex:
    return sum(a[i][i] for i in range(len(a)))


def projector_from_vector(v: list[complex]) -> Matrix:
    return [[v[i] * v[j].conjugate() for j in range(len(v))] for i in range(len(v))]


def gen_eigenvector(n: int) -> list[complex]:
    return [
        cmath.exp(-2j * math.pi * n * j / N_GEN) / math.sqrt(N_GEN)
        for j in range(N_GEN)
    ]


def gen_projector(n: int) -> Matrix:
    return projector_from_vector(gen_eigenvector(n))


def color_projector(index: int = 0) -> Matrix:
    out = zeros(N_COLOR, N_COLOR)
    out[index][index] = 1.0 + 0j
    return out


def boltzmann_color_projector(delta: float = 2 / 9) -> Matrix:
    # R=(0,1), G=(1,0), B=(1,1) have color-bit Hamming weights 1,1,2.
    raw = [math.exp(-delta), math.exp(-delta), math.exp(-2 * delta)]
    norm = math.sqrt(sum(value * value for value in raw))
    vector = [value / norm for value in raw]
    return projector_from_vector([complex(value) for value in vector])


def color_trace_factor(c: Matrix) -> float:
    return float((trace(c) / N_COLOR).real)


def circulant_amplitude(r_value: float, delta: float) -> Matrix:
    k = identity(N_GEN)
    for source in range(N_GEN):
        k[(source + 1) % N_GEN][source] += (r_value / 2) * cmath.exp(1j * delta)
        k[(source - 1) % N_GEN][source] += (r_value / 2) * cmath.exp(-1j * delta)
    return k


def lambda_value(r_value: float, delta: float, n: int) -> float:
    return 1 + r_value * math.cos(delta + 2 * math.pi * n / N_GEN)


def mass_values(r_value: float, delta: float) -> list[float]:
    return [lambda_value(r_value, delta, n) ** 2 for n in range(N_GEN)]


def color_resolved_mass_operator(
    r_value: float,
    delta: float,
    rank_one_color_block: Matrix,
    rank_one_index: int = LIGHT_INDEX,
) -> Matrix:
    values = mass_values(r_value, delta)
    color_identity = identity(N_COLOR)
    out = zeros(N_GEN * N_COLOR, N_GEN * N_COLOR)
    for n, value in enumerate(values):
        color_block = rank_one_color_block if n == rank_one_index else color_identity
        out = matadd(out, matscale(value, kron(gen_projector(n), color_block)))
    return out


def partial_color_trace(m: Matrix) -> Matrix:
    out = zeros(N_GEN, N_GEN)
    for a in range(N_GEN):
        for b in range(N_GEN):
            acc = 0j
            for c in range(N_COLOR):
                acc += m[a * N_COLOR + c][b * N_COLOR + c]
            out[a][b] = acc / N_COLOR
    return out


def rayleigh(a: Matrix, v: list[complex]) -> complex:
    av = [sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(v))]
    return sum(v[i].conjugate() * av[i] for i in range(len(v)))


def max_abs_matrix(a: Matrix) -> float:
    return max(abs(value) for row in a for value in row)


def matrix_sub(a: Matrix, b: Matrix) -> Matrix:
    return [[a[i][j] - b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def observed_log_ratios(sector: str) -> dict[str, float]:
    names = SECTOR_ORDER[sector]
    heavy = names[-1]
    return {
        name: math.log10(MASSES_GEV[sector][name] / MASSES_GEV[sector][heavy])
        for name in names
    }


def predicted_logs_from_operator(
    r_value: float,
    delta: float,
    rank_one_color_block: Matrix,
    sector: str,
    rank_one_index: int = LIGHT_INDEX,
) -> dict[str, float]:
    traced = partial_color_trace(
        color_resolved_mass_operator(
            r_value,
            delta,
            rank_one_color_block,
            rank_one_index=rank_one_index,
        )
    )
    eigen_masses = [rayleigh(traced, gen_eigenvector(n)).real for n in range(N_GEN)]
    # Existing Koide assignment: n=1 light, n=2 middle, n=0 heavy.
    names = SECTOR_ORDER[sector]
    assigned = {
        names[0]: eigen_masses[LIGHT_INDEX],
        names[1]: eigen_masses[MIDDLE_INDEX],
        names[2]: eigen_masses[HEAVY_INDEX],
    }
    heavy = names[-1]
    return {name: math.log10(assigned[name] / assigned[heavy]) for name in names}


def ratio_errors(sector: str, predicted: dict[str, float]) -> dict[str, tuple[float, float]]:
    observed = observed_log_ratios(sector)
    out = {}
    for name, obs in observed.items():
        if abs(obs) < 1e-12:
            continue
        err = predicted[name] - obs
        out[name] = (err, abs(err) / abs(obs))
    return out


def max_rel_error(sector: str, predicted: dict[str, float]) -> float:
    return max(rel for _, rel in ratio_errors(sector, predicted).values())


def print_table(title: str, sector: str, predicted: dict[str, float]) -> float:
    observed = observed_log_ratios(sector)
    errors = ratio_errors(sector, predicted)
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


@dataclass(frozen=True)
class Branch:
    name: str
    sector: str
    r_value: float
    delta: float
    source: str


def operator_sanity() -> None:
    hd("A. Circulant Mass Operator Sanity Check")
    r_value = math.sqrt(2)
    delta = 2 / 9
    k = circulant_amplitude(r_value, delta)
    k2 = matmul(k, k)
    spectral = zeros(N_GEN, N_GEN)
    for n, value in enumerate(mass_values(r_value, delta)):
        spectral = matadd(spectral, matscale(value, gen_projector(n)))
    residual = max_abs_matrix(matrix_sub(k2, spectral))
    print("  K(R,delta) = I + R/2 (exp(i delta) S + exp(-i delta) S^dagger)")
    for n in range(N_GEN):
        lam = rayleigh(k, gen_eigenvector(n)).real
        target = lambda_value(r_value, delta, n)
        print(f"  n={n}: eigenvalue {lam:.12f}, target {target:.12f}")
    print(f"  ||K^2 - sum_n lambda_n^2 P_n||_max = {residual:.3e}")
    assert residual < 1e-12


def color_node_sanity() -> None:
    hd("B. Rank-One Color Trace Is a Normalized-Rank Identity")
    color_identity = identity(N_COLOR)
    single = color_projector(0)
    boltz = boltzmann_color_projector()
    print(f"  tau_color(I_3)                  = {color_trace_factor(color_identity):.12f}")
    print(f"  tau_color(|r><r|)               = {color_trace_factor(single):.12f}")
    print(f"  tau_color(|w_Boltzmann><w|)     = {color_trace_factor(boltz):.12f}")
    print("  rank-one condition               = Tr(C_node)=1 instead of Tr(I_3)=3")
    print("  operator consequence             = whichever tagged node gets multiplied by 1/3")
    print("  status                           = tautological until the node choice is derived")
    assert abs(color_trace_factor(single) - 1 / 3) < 1e-12
    assert abs(color_trace_factor(boltz) - 1 / 3) < 1e-12


def rank_one_assignment_probe() -> None:
    hd("C. Rank-One Assignment Probe")
    rank_one_color = color_projector(0)
    print(
        "  The same 9x9 ansatz gives a 1/3 factor to whichever generation is "
        "tagged rank-one. The light-node choice is therefore not derived here."
    )
    print("\n  rank-one block on        up max rel    down A max rel")
    for node in [LIGHT_INDEX, MIDDLE_INDEX, HEAVY_INDEX]:
        up_pred = predicted_logs_from_operator(
            math.sqrt(3),
            2 / 27,
            rank_one_color,
            "up",
            rank_one_index=node,
        )
        down_pred = predicted_logs_from_operator(
            math.sqrt(2),
            1 / 9,
            rank_one_color,
            "down",
            rank_one_index=node,
        )
        print(
            f"  {NODE_LABELS[node]:<22s} "
            f"{max_rel_error('up', up_pred):>10.1%}   "
            f"{max_rel_error('down', down_pred):>13.1%}"
        )
    print(
        "  -> The fit selects the light node. This file supplies no independent "
        "reason that light, rather than middle or heavy, carries rank one."
    )


def branch_audit() -> None:
    hd("D. Branch Audit: Does the Trace Beat No-Trace Alternatives?")
    identity_color = identity(N_COLOR)
    rank_one_color = color_projector(0)
    branches = [
        Branch("up bare color paths", "up", math.sqrt(3), 2 / 27, "ANCHOR item 96 / Part 01"),
        Branch("down A two-path baseline", "down", math.sqrt(2), 1 / 9, "unfitted two-path branch"),
        Branch("down B color-dressed R", "down", math.sqrt(12 / 5), 1 / 9, "Part 01 near-fit R"),
    ]

    print(
        "  Each branch is evaluated twice from the same 9x9 operator: "
        "C_light=I_3 (no node trace) and C_light=rank-one (node trace)."
    )
    print("\n  branch                       no trace    rank-one node trace   reading")
    for branch in branches:
        no_trace = predicted_logs_from_operator(
            branch.r_value,
            branch.delta,
            identity_color,
            branch.sector,
        )
        traced = predicted_logs_from_operator(
            branch.r_value,
            branch.delta,
            rank_one_color,
            branch.sector,
        )
        no_err = max_rel_error(branch.sector, no_trace)
        tr_err = max_rel_error(branch.sector, traced)
        if tr_err <= 0.10 and no_err > 0.10:
            reading = "trace supplies missing factor"
        elif no_err <= 0.10 and tr_err > 0.10:
            reading = "trace disfavored / double-counts"
        elif no_err <= 0.10 and tr_err <= 0.10:
            reading = "ambiguous"
        else:
            reading = "neither branch passes"
        print(f"  {branch.name:<28s} {no_err:>8.1%} {tr_err:>18.1%}   {reading}")

    up_traced = predicted_logs_from_operator(math.sqrt(3), 2 / 27, rank_one_color, "up")
    down_a_traced = predicted_logs_from_operator(math.sqrt(2), 1 / 9, rank_one_color, "down")
    down_b_notr = predicted_logs_from_operator(math.sqrt(12 / 5), 1 / 9, identity_color, "down")
    print_table("up, R=sqrt3, delta=2/27, rank-one light color node", "up", up_traced)
    print_table("down A, R=sqrt2, delta=1/9, rank-one light color node", "down", down_a_traced)
    print_table("down B, R=sqrt(12/5), delta=1/9, no node trace", "down", down_b_notr)


def lemma_status() -> None:
    hd("E. Lemma Status")
    print(
        "  Conditional ansatz result:\n"
        "    M = sum_n lambda_n^2 P_n x C_n with C_light rank-one gives\n"
        "    partial_color_trace(M) = same Koide mass operator but with\n"
        "    m_light -> m_light/3. This is an operator identity after the\n"
        "    rank-one light block has been inserted, not a derivation of that block."
    )
    print(
        "\n  What the branch audit says:\n"
        "    Up sector: R=sqrt3 needs the rank-one node trace to pass the 10% bar.\n"
        "    Down branch A: R=sqrt2 also needs the rank-one node trace and then passes.\n"
        "    Down branch B: R=sqrt(12/5) passes without the trace and fails with it,\n"
        "    so the best down-sector description argues against adding a separate trace."
    )
    print(
        "\n  Remaining canonical proof obligation:\n"
        "    derive C_light=rank-one and C_middle=C_heavy=I_3 from an actual\n"
        "    color-resolved Yukawa operator. The documented 256-state walk operator\n"
        "    was checked in tch_256_walk_mass_rank_audit.py and does not produce\n"
        "    this pattern."
    )


def main() -> None:
    operator_sanity()
    color_node_sanity()
    rank_one_assignment_probe()
    branch_audit()
    lemma_status()


if __name__ == "__main__":
    main()
