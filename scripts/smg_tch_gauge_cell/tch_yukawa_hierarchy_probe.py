#!/usr/bin/env python3
"""
TCH charged-fermion Yukawa hierarchy probe.

Goal:
  Test whether existing TCH invariants can explain the charged-fermion Yukawa
  hierarchy, rather than merely restating y_f proportional to m_f.

Quality bar:
  Match log10 mass/Yukawa ratios within each charged sector at ~10% precision
  without fitted continuous shape parameters. An overall sector scale is allowed
  only when explicitly testing shape, because ratios are the hierarchy question.

Existing framework ingredients tested:
  1. The literal Q3 edge-frustration F used in narrow_higgs:
        M(c) ~ exp(F(c)/(2 phi)).
  2. The CKM/shell-routing generation dressing:
        Gen 1 local, Gen 2 shell 4*12, Gen 3 shell 9*18.
  3. Discrete Koide-sector shapes:
        m_n = mu * (1 + R cos(delta + 2 pi n/3))^2
     with anchored candidate phases such as lepton delta=2/9, up delta=2/27,
     down delta=1/9.

The script is an audit, not a fit paper. It reports what works, what fails,
and how many continuous shape parameters would be needed to repair a failure.
"""

import itertools
import math

import numpy as np


PHI = (1 + 5**0.5) / 2
Q3_EDGES = [
    (0, 1),
    (0, 2),
    (0, 4),
    (1, 3),
    (1, 5),
    (2, 3),
    (2, 6),
    (3, 7),
    (4, 5),
    (4, 6),
    (5, 7),
    (6, 7),
]

# Representative low-scale masses. Ratios are what matter here; y_f ratios are
# identical to mass ratios because y_f = sqrt(2) m_f / v.
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

GEN_BITS = {
    1: (0, 0),
    2: (0, 1),
    3: (1, 0),
}

GEN_BY_PARTICLE = {
    "u": 1,
    "d": 1,
    "e": 1,
    "c": 2,
    "s": 2,
    "mu": 2,
    "t": 3,
    "b": 3,
    "tau": 3,
}

COLOR_BITS = [(0, 1), (1, 0), (1, 1)]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def q3_frustration(bits):
    return sum(bits[left] ^ bits[right] for left, right in Q3_EDGES)


def codeword(gen, kind, color=(0, 1), chirality=(0, 0)):
    """Return bits in canonical order G0,G1,LQ,C0,C1,I3,chi,W."""
    g0, g1 = GEN_BITS[gen]
    chi, weak = chirality
    if kind == "e":
        return (g0, g1, 0, 0, 0, 1, chi, weak)
    if kind == "u":
        c0, c1 = color
        return (g0, g1, 1, c0, c1, 0, chi, weak)
    if kind == "d":
        c0, c1 = color
        return (g0, g1, 1, c0, c1, 1, chi, weak)
    raise ValueError(kind)


def averaged_frustration(particle):
    gen = GEN_BY_PARTICLE[particle]
    if particle in ["e", "mu", "tau"]:
        return q3_frustration(codeword(gen, "e"))
    kind = "u" if particle in ["u", "c", "t"] else "d"
    return sum(q3_frustration(codeword(gen, kind, color)) for color in COLOR_BITS) / 3


def observed_log_ratios(sector):
    names = SECTOR_ORDER[sector]
    heavy = names[-1]
    return {name: math.log10(MASSES_GEV[sector][name] / MASSES_GEV[sector][heavy]) for name in names}


def shape_log_ratios_from_scores(names, scores):
    heavy_score = scores[names[-1]]
    return {name: (scores[name] - heavy_score) / math.log(10) for name in names}


def ratio_errors(observed, predicted):
    out = {}
    for name, obs in observed.items():
        if abs(obs) < 1e-12:
            continue
        err = predicted[name] - obs
        out[name] = (err, abs(err) / abs(obs))
    return out


def print_ratio_table(title, sector, predicted):
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


def canonical_friction_scores(include_shell=True):
    scores = {}
    shell = {
        1: 1.0,
        2: 4 * 12,
        3: 9 * 18,
    }
    for sector, names in SECTOR_ORDER.items():
        for name in names:
            f_value = averaged_frustration(name)
            score = f_value / (2 * PHI)
            if include_shell:
                score += math.log(shell[GEN_BY_PARTICLE[name]])
            scores[name] = score
    return scores


def koide_shape_values(r_value, delta):
    return [
        (1 + r_value * math.cos(delta + 2 * math.pi * index / 3)) ** 2
        for index in range(3)
    ]


def koide_predicted_logs(sector, r_value, delta, assignment=(1, 2, 0)):
    """Default assignment: light=n1, middle=n2, heavy=n0."""
    names = SECTOR_ORDER[sector]
    values = koide_shape_values(r_value, delta)
    assigned = {name: math.log10(values[index]) for name, index in zip(names, assignment)}
    heavy = names[-1]
    return {name: assigned[name] - assigned[heavy] for name in names}


def best_one_slope_generation_model():
    """Fit only a universal slope plus one offset per sector as a diagnostic."""
    # Generation tier candidates built from the two-bit generation labels.
    tiers = {
        "linear": {1: 0, 2: 1, 3: 2},
        "hamming_to_gen3": {
            gen: sum(a ^ b for a, b in zip(GEN_BITS[gen], GEN_BITS[3]))
            for gen in [1, 2, 3]
        },
        "shell_log": {1: math.log(1), 2: math.log(4 * 12), 3: math.log(9 * 18)},
    }
    rows = []
    for name, tier in tiers.items():
        xs = []
        ys = []
        sector_cols = []
        particles = []
        for sector, names in SECTOR_ORDER.items():
            for particle in names:
                xs.append(tier[GEN_BY_PARTICLE[particle]])
                ys.append(math.log10(MASSES_GEV[sector][particle]))
                sector_cols.append(sector)
                particles.append(particle)

        sectors = list(SECTOR_ORDER)
        # Design: universal slope + sector offsets.
        design = []
        for x, sector in zip(xs, sector_cols):
            row = [x]
            row.extend(1.0 if sector == s else 0.0 for s in sectors)
            design.append(row)
        beta, *_ = np.linalg.lstsq(np.array(design), np.array(ys), rcond=None)
        pred = np.array(design) @ beta
        rms = float(np.sqrt(np.mean((pred - ys) ** 2)))
        max_abs = float(np.max(np.abs(pred - ys)))
        rows.append((rms, max_abs, name, beta[0], dict(zip(particles, pred - ys))))
    return sorted(rows)


def empirical_summary():
    hd("A. Empirical Charged-Fermion Hierarchy")
    for sector, names in SECTOR_ORDER.items():
        print(f"  {sector}:")
        heavy = names[-1]
        for name in names:
            ratio = MASSES_GEV[sector][heavy] / MASSES_GEV[sector][name]
            print(f"    {heavy}/{name:<3s} = {ratio:>10.3g}, log10 = {math.log10(ratio):.4f}")


def test_canonical_friction():
    hd("B. Existing TCH F + Shell-Routing Formula")
    print("  F values from canonical 8-bit register on Q3 edges:")
    for sector, names in SECTOR_ORDER.items():
        print("  " + sector + ": " + ", ".join(f"{name}:F={averaged_frustration(name):.3g}" for name in names))

    scores = canonical_friction_scores(include_shell=True)
    max_errors = []
    for sector, names in SECTOR_ORDER.items():
        predicted = shape_log_ratios_from_scores(names, scores)
        max_errors.append(print_ratio_table("M ~ exp(F/(2phi)) * shell", sector, predicted))
    print(f"\n  overall max relative log-ratio error = {max(max_errors):.1%}")
    return max(max_errors)


def test_koide_candidates():
    hd("C. Discrete Koide-Shape Candidates")
    candidates = {
        "lepton anchored: R=sqrt2, delta=2/9": ("lepton", math.sqrt(2), 2 / 9),
        "up candidate: R=sqrt3, delta=2/27": ("up", math.sqrt(3), 2 / 27),
        "down candidate: R=sqrt3, delta=1/9": ("down", math.sqrt(3), 1 / 9),
        "down alt: R=sqrt2, delta=1/9": ("down", math.sqrt(2), 1 / 9),
    }
    max_errors = {}
    for title, (sector, r_value, delta) in candidates.items():
        predicted = koide_predicted_logs(sector, r_value, delta)
        max_errors[title] = print_ratio_table(title, sector, predicted)
    return max_errors


def diagnostic_fit():
    hd("D. Best Universal Generation-Tier Diagnostic Fit")
    print(
        "  This diagnostic allows one fitted universal slope plus one offset per "
        "sector. If even this fails, generation tiers alone cannot explain the "
        "9 Yukawas."
    )
    print("  tier model          rms dex   max abs dex   fitted slope")
    for rms, max_abs, name, slope, residuals in best_one_slope_generation_model():
        print(f"  {name:<18s} {rms:<9.4f} {max_abs:<12.4f} {slope:+.4f}")
    print("  -> This is a fit diagnostic only, not a parameter-free TCH result.")


def corrected_predicted_logs(base, correction):
    corrected = {}
    for particle, log_value in base.items():
        corrected[particle] = log_value + math.log10(correction.get(particle, 1.0))
    return corrected


def max_rel_for_sector(sector, predicted):
    return max(rel for _, rel in ratio_errors(observed_log_ratios(sector), predicted).values())


def integer_node_correction_search():
    hd("E. Exploratory Integer Node-Correction Search")
    print(
        "  This is deliberately labelled exploratory. It searches small integer "
        "denominators applied only to the first-generation quark node states "
        "(u,d), because the existing quark notes already identify the lightest "
        "quarks as spectral-node sensitive. A hit here is a target to derive, "
        "not a derivation."
    )
    base_up = koide_predicted_logs("up", math.sqrt(3), 2 / 27)
    base_down = koide_predicted_logs("down", math.sqrt(2), 1 / 9)
    trials = []
    denominators = range(2, 13)
    for du in denominators:
        for dd in denominators:
            pred_up = corrected_predicted_logs(base_up, {"u": 1 / du})
            pred_down = corrected_predicted_logs(base_down, {"d": 1 / dd})
            up_err = max_rel_for_sector("up", pred_up)
            down_err = max_rel_for_sector("down", pred_down)
            trials.append((max(up_err, down_err), up_err, down_err, du, dd, pred_up, pred_down))
    trials.sort(key=lambda row: row[0])
    print(f"  searched {len(trials)} pairs: d_u,d_d in [{min(denominators)}, {max(denominators)}]")
    print("  best candidates:")
    print("    d_u  d_d   max rel err   up rel err   down rel err")
    for max_err, up_err, down_err, du, dd, _, _ in trials[:5]:
        print(f"    {du:<4d} {dd:<5d} {max_err:<13.1%} {up_err:<12.1%} {down_err:<.1%}")

    best = trials[0]
    _, _, _, du, dd, pred_up, pred_down = best
    print_ratio_table(f"up with node correction u -> u/{du}", "up", pred_up)
    print_ratio_table(f"down with node correction d -> d/{dd}", "down", pred_down)
    print(
        "  Candidate target: derive the integer suppressions "
        f"u-node 1/{du} and d-node 1/{dd} from colour/bipartite trace "
        "rather than fitting them."
    )
    return best


def main():
    empirical_summary()
    f_error = test_canonical_friction()
    koide_errors = test_koide_candidates()
    diagnostic_fit()
    best_node = integer_node_correction_search()

    hd("Verdict")
    lepton_ok = koide_errors["lepton anchored: R=sqrt2, delta=2/9"] <= 0.10
    quark_ok = all(
        error <= 0.10
        for key, error in koide_errors.items()
        if key.startswith("up") or key.startswith("down candidate")
    )
    print(
        "  Charged leptons: "
        + ("PASS" if lepton_ok else "FAIL")
        + " under the existing discrete Koide shape."
    )
    print(
        "  Quarks: "
        + ("PASS" if quark_ok else "FAIL")
        + " under the current no-continuous-shape candidates."
    )
    print(
        "  Canonical F+shell friction: "
        + ("PASS" if f_error <= 0.10 else "FAIL")
        + " for the all-sector 10% log-ratio bar."
    )
    print(
        "\n  Honest conclusion: TCH has a real charged-lepton hierarchy handle, "
        "but the current Hamming/friction/shell invariants do not yet derive "
        "all 9 charged Yukawas. The missing piece is a quark-sector colour/QCD "
        "dressing invariant, not another fitted slope."
    )
    print(
        "  New concrete target from this run: the smallest useful correction is "
        f"an integer node suppression u->u/{best_node[3]}, d->d/{best_node[4]}. "
        "That pattern passes the ratio bar, but it is currently a searched "
        "postulate over 121 integer pairs, not a TCH theorem."
    )
    corrected_max = max(
        koide_errors["lepton anchored: R=sqrt2, delta=2/9"],
        best_node[1],
        best_node[2],
    )
    print(
        f"  Combined hierarchy-shape score after that discrete correction: "
        f"max relative log-ratio error = {corrected_max:.1%}."
    )
    if not lepton_ok or quark_ok or f_error <= 0.10:
        raise SystemExit("unexpected verdict changed; inspect assumptions")


if __name__ == "__main__":
    main()
