#!/usr/bin/env python3
r"""Strict thermodynamic recovery protocol for ANCHOR item 126.

The old item-126 claim treated

    eta = (3/14) alpha^4 + (1/3) alpha^5

as a parameter-free baryogenesis derivation. The earlier audits already showed the
problem: alpha^4 is a real ideal-code fault scale, but the coefficient was imported
from an ideal [8,4,4] count rather than derived from the actual 48-state matter
set, and the alpha^5 term is a fit-improving correction.

This script turns that into a strict recovery ledger using ANCHOR section 16.4:
separate the recoverable scale, the conditional ideal-code branch fraction, and
the still-missing microscopic event-current theorem.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product


ALPHA = 1 / 137.036
ETA_OBS = 6.12e-10
ETA_ERR = 0.04e-10

# Standard generator for the ideal extended Hamming [8,4,4] code, with bit i
# packed into integer bit i to match the item126 colour-count audit.
GENERATOR = (
    (1, 0, 0, 0, 0, 1, 1, 1),
    (0, 1, 0, 0, 1, 0, 1, 1),
    (0, 0, 1, 0, 1, 1, 0, 1),
    (0, 0, 0, 1, 1, 1, 1, 0),
)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def weight(n: int) -> int:
    return n.bit_count()


def colour(n: int) -> tuple[int, int]:
    return bit(n, 3), bit(n, 4)


def pack(bits: tuple[int, ...]) -> int:
    return sum(b << i for i, b in enumerate(bits))


def xor_rows(selector: tuple[int, ...]) -> int:
    out = [0] * 8
    for take, row in zip(selector, GENERATOR):
        if take:
            out = [(a ^ b) for a, b in zip(out, row)]
    return pack(tuple(out))


def ideal_844_codewords() -> list[int]:
    return [xor_rows(selector) for selector in product((0, 1), repeat=4)]


def in_matter_48(n: int) -> bool:
    """Canonical 48-state R1/R2/R3 matter-set rules used by the audits."""

    if bit(n, 0) and bit(n, 1):
        return False
    if bit(n, 7) != bit(n, 6):
        return False

    c = colour(n)
    if bit(n, 2) == 0 and c != (0, 0):
        return False
    if bit(n, 2) == 1 and c == (0, 0):
        return False

    return True


def pct_delta(value: float, target: float) -> float:
    return 100.0 * (value - target) / target


def sigma(value: float, target: float = ETA_OBS, err: float = ETA_ERR) -> float:
    return abs(value - target) / err


def rational_hits(required: float, tolerance: float, qmax: int = 29) -> list[Fraction]:
    hits: set[Fraction] = set()
    for q in range(1, qmax + 1):
        for p in range(1, q):
            r = Fraction(p, q)
            if abs(float(r) - required) / required <= tolerance:
                hits.add(r)
    return sorted(hits, key=float)


def nearest_ratio(candidates: dict[str, Fraction], target: float) -> tuple[str, Fraction, float]:
    name, ratio = min(candidates.items(), key=lambda kv: abs(float(kv[1]) - target))
    return name, ratio, abs(float(ratio) - target) / target


def main() -> None:
    a4 = ALPHA**4
    a5 = ALPHA**5
    eta_leading = float(Fraction(3, 14)) * a4
    eta_full = eta_leading + float(Fraction(1, 3)) * a5
    required_coeff = ETA_OBS / a4
    leading_tolerance = abs(eta_leading - ETA_OBS) / ETA_OBS

    print("ITEM 126 STRICT THERMODYNAMIC RECOVERY PROTOCOL")
    print("=" * 72)
    print("[1] numerical target")
    print(f"alpha^4                       = {a4:.6e}")
    print(
        f"(3/14) alpha^4                = {eta_leading:.6e}  "
        f"delta={pct_delta(eta_leading, ETA_OBS):+.2f}%  "
        f"sigma={sigma(eta_leading):.2f}"
    )
    print(
        f"(3/14) alpha^4 + (1/3)alpha^5 = {eta_full:.6e}  "
        f"delta={pct_delta(eta_full, ETA_OBS):+.2f}%  "
        f"sigma={sigma(eta_full):.2f}"
    )
    print(f"required coefficient eta_obs/alpha^4 = {required_coeff:.6f}")

    assert abs(eta_leading - 6.076e-10) < 2e-13
    assert abs(eta_full - 6.145e-10) < 3e-13

    ideal = ideal_844_codewords()
    ideal_weights = Counter(weight(n) for n in ideal)
    ideal_w4 = [n for n in ideal if weight(n) == 4]
    ideal_colour_counts = Counter(colour(n) for n in ideal_w4)
    ideal_colourless = [n for n in ideal_w4 if colour(n) == (0, 0)]

    print("\n[2] ideal [8,4,4] code facts")
    print(f"codeword count        = {len(ideal)}")
    print(f"weight enumerator     = {dict(sorted(ideal_weights.items()))}")
    print(f"weight-4 colour split = {dict(sorted(ideal_colour_counts.items()))}")
    print(f"colourless/weight-4   = {len(ideal_colourless)}/{len(ideal_w4)}")

    assert len(ideal) == 16
    assert ideal_weights == Counter({0: 1, 4: 14, 8: 1})
    assert len(ideal_colourless) == 3

    hits = rational_hits(required_coeff, leading_tolerance)
    print("\n[3] coefficient uniqueness at the leading tolerance")
    print(f"tolerance from leading match = {100 * leading_tolerance:.3f}%")
    print("simple p/q hits, q <= 29    = " + ", ".join(str(h) for h in hits))
    print("interpretation              = 3/14 is a real ideal-code coincidence")
    assert hits == [Fraction(3, 14)]

    matter = [n for n in range(256) if in_matter_48(n)]
    matter_weights = Counter(weight(n) for n in matter)
    min_nonzero_weight = min(weight(n) for n in matter if n)
    leptons = [n for n in matter if bit(n, 2) == 0]
    quarks = [n for n in matter if bit(n, 2) == 1]
    sterile_colourless = [
        n for n in matter if bit(n, 2) == 0 and colour(n) == (0, 0) and bit(n, 6) == 1 and bit(n, 7) == 1
    ]
    colourless_w4 = [n for n in matter if weight(n) == 4 and colour(n) == (0, 0)]

    print("\n[4] actual 48-state matter set")
    print(f"state count           = {len(matter)}")
    print(f"weight distribution   = {dict(sorted(matter_weights.items()))}")
    print(f"minimum nonzero weight= {min_nonzero_weight}")
    print(f"weight-4 count        = {matter_weights[4]}")

    assert len(matter) == 48
    assert min_nonzero_weight == 1
    assert matter_weights[4] == 11

    matter_ratios = {
        "sterile colourless / 48": Fraction(len(sterile_colourless), len(matter)),
        "leptons / 48": Fraction(len(leptons), len(matter)),
        "quarks / 48": Fraction(len(quarks), len(matter)),
        "colourless weight-4 / weight-4": Fraction(len(colourless_w4), matter_weights[4]),
        "leptons / quarks": Fraction(len(leptons), len(quarks)),
    }
    print("\n[5] principled 48-set ratios")
    for name, ratio in matter_ratios.items():
        miss = abs(float(ratio) - required_coeff) / required_coeff
        print(f"{name:<34} = {str(ratio):>5} = {float(ratio):.6f}  miss={100 * miss:.1f}%")

    best_name, best_ratio, best_miss = nearest_ratio(matter_ratios, required_coeff)
    print(f"nearest 48-set ratio          = {best_name} = {best_ratio}, miss={100 * best_miss:.1f}%")
    assert best_miss > 0.10

    print("\n[6] section 16.4 thermodynamic gates")
    gates = (
        ("T1 event algebra", "CONDITIONAL", "ideal [8,4,4] count explicit; actual 48-state baryon Kraus map open"),
        ("T2 Landauer status", "CONDITIONAL", "alpha is used as fault amplitude/rate; no thermal equality coefficient is derived"),
        ("T3 mean current", "OPEN", "no absolute baryon-producing service current or freeze-out event density"),
        ("T4 covariance/Fano", "OPEN", "no fluctuation ledger for survival asymmetry"),
        ("T5 correlation volume", "OPEN", "no mode/cell/horizon coarse-graining volume for eta"),
        ("T6 observable map", "OPEN", "branching fraction to baryon-to-photon ratio is not derived"),
        ("T7 scale accounting", "PASS", "dimensionless eta; no hidden H0 or horizon scale in this coefficient audit"),
        ("T8 alternatives/nulls", "PASS", "48-set null and alpha^5 fit-improver are explicitly separated"),
        ("T9 promotion statement", "PASS", "scale, coefficient, correction, and open theorem are split below"),
    )
    for gate, status, reason in gates:
        print(f"{gate:<22} {status:<12} {reason}")

    print("\n[7] recovery verdict")
    print("RECOVERED:")
    print("  * eta ~ alpha^4 as an ideal logical-fault scale.")
    print("  * 3/14 as a real, fairly unique ideal [8,4,4] colourless/weight-4 count.")
    print("NOT RECOVERED:")
    print("  * a parameter-free baryogenesis theorem on the actual 48-state matter set.")
    print("  * the old 0.6 sigma headline, because it uses the extra +(1/3) alpha^5 fit improver.")
    print("  * the event-current map from QEC/Landauer ledger to eta = n_B/n_gamma.")
    print("NEXT THEOREM TARGET:")
    print("  Derive a baryon-producing Kraus/current ledger on the actual matter state space,")
    print("  or prove a valid three-generation [8,4,4] lift whose 3/14 branch fraction maps")
    print("  to baryon number and to the photon exhaust ledger.")
    print("\nexit 0 -- item126 recovered as an alpha^4 scale plus conditional ideal-code proposition, not Locked.")


if __name__ == "__main__":
    main()
