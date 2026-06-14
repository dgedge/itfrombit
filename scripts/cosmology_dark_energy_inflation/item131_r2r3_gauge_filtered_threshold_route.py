#!/usr/bin/env python3
"""ITEM 131 / 42: gauge-filtered route to N_print alpha^4 = 4/3.

The raw R2/R3 threshold audit found a real blocker: Boolean R2/R3 violations
occur at weight 1, so R2/R3 alone cannot force an alpha^4 threshold.

This audit asks the narrower question that survived that failure:

    If "R2/R3 activation current" means the post-decoder, gauge/strain-filtered
    topology-changing current, can the actual ledger geometry force the alpha^4
    leg of

        N_print(H_*) alpha_0^4 = C_F,  C_F = 4/3 ?

Answer:
    It gives a plausible derivation route for the alpha^4 power, but not a full
    derivation of the threshold equation.

Closed by existing geometry:
    * the 12-edge strain decoder uniquely corrects all weight <= 3 faults;
    * true topology change is therefore the weight-4 logical/portal sector;
    * the licensed weight-4 commit budget is alpha_0^4 per photon/exhaust unit.

Still open:
    * prove item-42 activation samples the post-decoder topology-changing current
      rather than the raw R2/R3 Boolean violation current;
    * prove the stopping load is the unbroken SU(3) fundamental Casimir C_F=4/3;
    * prove the stop rule E[N_topological commits per printed shell] = C_F.
"""

from __future__ import annotations

import itertools
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
ALPHA0 = 1.0 / 137.0
C_F = Fraction(4, 3)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def bits3(x: int) -> tuple[int, int, int]:
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


def hyperplanes() -> list[frozenset[int]]:
    out: set[frozenset[int]] = set()
    for a in range(1, 8):
        avec = bits3(a)
        for b in (0, 1):
            out.add(
                frozenset(
                    p
                    for p in range(8)
                    if sum(u * v for u, v in zip(bits3(p), avec)) % 2 == b
                )
            )
    return sorted(out, key=lambda h: tuple(sorted(h)))


POINTS = tuple(range(8))
EDGES = tuple(
    (u, v)
    for u in POINTS
    for v in POINTS
    if u < v and bin(u ^ v).count("1") == 1
)


def strain_syndrome(support: frozenset[int]) -> tuple[int, ...]:
    return tuple(((u in support) + (v in support)) % 2 for u, v in EDGES)


def code_syndrome(support: frozenset[int]) -> tuple[int, ...]:
    checks = [
        frozenset(POINTS),
        frozenset(p for p in POINTS if bits3(p)[0] == 1),
        frozenset(p for p in POINTS if bits3(p)[1] == 1),
        frozenset(p for p in POINTS if bits3(p)[2] == 1),
    ]
    return tuple(len(support & check_word) % 2 for check_word in checks)


def source_contains(path: str, needle: str) -> bool:
    return needle in (ROOT / path).read_text(encoding="utf-8")


def main() -> None:
    print("ITEM 131 / 42 GAUGE-FILTERED R2/R3 THRESHOLD ROUTE")

    print("\n[1] Strain filtering removes lower-weight activation candidates")
    by_syndrome: dict[tuple[int, ...], list[frozenset[int]]] = {}
    for weight in range(0, 9):
        for support in itertools.combinations(POINTS, weight):
            fs = frozenset(support)
            by_syndrome.setdefault(strain_syndrome(fs), []).append(fs)

    low_weight_unique = True
    for syndrome, supports in by_syndrome.items():
        lows = [s for s in supports if len(s) <= 3]
        if lows:
            low_weight_unique = low_weight_unique and len(lows) == 1
        else:
            low_weight_unique = low_weight_unique and True
    check(low_weight_unique, "every weight <= 3 fault has a unique 12-edge strain syndrome")

    weight4_pairs = []
    all8 = frozenset(POINTS)
    weight4_complements = True
    for supports in by_syndrome.values():
        fours = [s for s in supports if len(s) == 4]
        if fours:
            weight4_complements = weight4_complements and len(fours) == 2 and fours[0] == all8 - fours[1]
            weight4_pairs.append(tuple(fours))
    check(weight4_complements, "weight-4 strain ambiguity is complement pairing")
    check(len(weight4_pairs) == 35, "70 weight-4 subsets reduce to 35 complement ambiguity classes")
    print("  Interpretation: raw one-, two-, and three-bit R2/R3 violations are")
    print("  detectable strain events.  They may heat or strain the ledger, but")
    print("  they are not the first uncorrected topology-changing current.")

    print("\n[2] The topology-changing sector is weight-4")
    hyps = hyperplanes()
    check(len(hyps) == 14 and all(len(h) == 4 for h in hyps), "[8,4,4] code has 14 weight-4 hyperplane logical supports")
    silent = [h for h in hyps if code_syndrome(h) == (0, 0, 0, 0)]
    check(len(silent) == 14, "all hyperplane supports are code-syndrome silent")
    tripped_counts = sorted({sum(strain_syndrome(h)) for h in hyps})
    check(tripped_counts == [4, 8, 12], "weight-4 logicals are strain/frustration visible")
    check(
        source_contains("python_code/feshbach_alpha2_target.py", "real TOPOLOGY CHANGE is a weight-4 logical fault"),
        "canon audit already separates topology change from weight-1 stiffness faults",
    )
    check(
        source_contains("python_code/scheduler_alpha_composition.py", "commits per photon = alpha_0^4"),
        "licensed weight-4 commit ledger gives alpha_0^4 per photon/exhaust unit",
    )
    print("  Consequence: the alpha^4 leg can be grounded if the item-42")
    print("  activation current is defined after the strain/QEC filter, i.e. as")
    print("  topology-changing logical commits rather than raw Boolean violations.")

    print("\n[3] The 4/3 leg is real but not yet selected by the threshold ledger")
    c_su3 = Fraction(3 * 3 - 1, 2 * 3)
    check(c_su3 == C_F, "R3 colour triplet gives SU(3) fundamental C_F=4/3")
    check(
        source_contains("ANCHOR.md", "Quark-sector defect-fraction scaling")
        and source_contains("ANCHOR.md", "OPEN. The"),
        "canon still flags SU(3) coherent-superposition scaling as open elsewhere",
    )
    print("  Interpretation: C_F=4/3 is available as the unbroken colour")
    print("  restoring-load invariant.  The prior audit already showed the")
    print("  broken Pati-Salam activation loads are 3/4, 1/2, or 1-class, not")
    print("  4/3.  Therefore the threshold must prove it couples to the colour")
    print("  restoring load, not merely to the broken R2/R3 coset.")

    print("\n[4] What the route would derive, and where it still stops")
    n_print_required = float(C_F) * ALPHA0 ** -4
    print(f"  If the three remaining clauses hold: N_print = C_F alpha^-4 = {n_print_required:.9e}")
    print(f"  Equivalently: N_print alpha^4 = {float(C_F):.12f}")
    gates = [
        (
            "alpha^4 power",
            "CONDITIONAL ROUTE",
            "closed if activation means post-decoder topology-changing current",
        ),
        (
            "4/3 coefficient",
            "PARTIAL",
            "R3 supplies C_F, but threshold must select colour-restoring load",
        ),
        (
            "stopping rule",
            "OPEN",
            "no ledger law yet says E[N_topological commits per shell] = C_F",
        ),
    ]
    for name, status, note in gates:
        print(f"  {name:16s} {status:18s} {note}")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  There is a way forward, but it is not a finished derivation.")
    print("  The failed raw-R2/R3 route should be replaced by a gauge/strain-filtered")
    print("  route:")
    print("")
    print("    raw lower-weight R2/R3 faults -> detected/corrected strain current;")
    print("    first uncorrected topology-changing current -> weight-4 logical sector;")
    print("    licensed topology-changing commits -> alpha_0^4 per printed-shell unit;")
    print("    colour restoring load -> candidate C_F=4/3;")
    print("    exit condition -> still-open equality E[N_topological]=C_F.")
    print("")
    print("  Thus the alpha^4 power is nearly derivable from existing QEC/Kraus")
    print("  geometry once the activation-current definition is made post-decoder.")
    print("  The coefficient and the equality remain the live theorem, not an")
    print("  accounting identity.  No promotion until those clauses are proved.")
    print("=" * 104)
    print("exit 0 -- route narrowed; alpha^4 rescue conditional, C_F/stopping rule still open.")


if __name__ == "__main__":
    main()
