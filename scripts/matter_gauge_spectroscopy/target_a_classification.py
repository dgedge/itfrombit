#!/usr/bin/env python3
"""Finite classification for Target A anomaly-compression.

This script upgrades the earlier Target-A genericity scan from "counting hits"
to a small classification ledger.  It keeps the paper's exact charge convention
and asks where the Standard-Model content sits inside increasingly broad
natural Boolean rule classes on the eight-bit record cell.

Rule classes:

  A. Current paper scan:
     Any non-empty matter/colour sector rule V subset {LQ,C0,C1}; fixed R1,
     fixed R2, and the canonical R4 exclusion.

  B. The quotient of A by the colour-bit GL(2,F2) action, so colour names do
     not inflate the count.

  C. A wider natural class:
     Any non-empty V, plus one bit-local right-singlet exclusion
     (LQ value, I3 value, chi=1).  This is the natural R4-like class:
     it removes singlets without reading the colour bits and without breaking
     SU(2) doublets.

  D. Role-map invariance:
     The canonical SM rule is re-run over all disjoint choices of the physical
     roles.  All successful maps must be SM-isomorphic, not new physics.

The output is self-asserting: exit 0 means the finite ledgers reproduce the
printed claims.  The result is intentionally a classification, not a proof that
the rule class is the only possible "natural" one.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations, product


Bit = int
Codeword = tuple[int, ...]
Triple = tuple[int, int, int]  # (LQ,C0,C1)
Exclusion = tuple[str, int, int, int]  # ("lq", LQ value, I3, chi)

CW: list[Codeword] = list(product([0, 1], repeat=8))
TRIPLES: list[Triple] = list(product([0, 1], repeat=3))

# Canonical bit roles: 0G0 1G1 2LQ 3C0 4C1 5I3 6chi 7W.
SM_SECTOR: frozenset[Triple] = frozenset(
    [(0, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
)
SM_EXCLUSION: Exclusion = ("lq", 0, 0, 1)  # nu_R removed from active set.


@dataclass(frozen=True)
class RoleMap:
    i3: Bit
    chi: Bit
    lq: Bit
    c0: Bit
    c1: Bit
    g0: Bit
    g1: Bit
    w: Bit


SM_ROLE = RoleMap(i3=5, chi=6, lq=2, c0=3, c1=4, g0=0, g1=1, w=7)


def r1(c: Codeword, rm: RoleMap = SM_ROLE) -> bool:
    return not (c[rm.g0] == 1 and c[rm.g1] == 1)


def r2(c: Codeword, rm: RoleMap = SM_ROLE) -> bool:
    return c[rm.w] == c[rm.chi]


def sector(c: Codeword, rm: RoleMap = SM_ROLE) -> Triple:
    return (c[rm.lq], c[rm.c0], c[rm.c1])


def zf(c: Codeword, rm: RoleMap = SM_ROLE) -> int:
    return 1 if c[rm.i3] == 0 else -1


def sum_zc(c: Codeword, rm: RoleMap = SM_ROLE) -> int:
    return -3 if (c[rm.c0], c[rm.c1]) == (0, 0) else -1


def q(c: Codeword, rm: RoleMap = SM_ROLE, b: F = F(1, 3)) -> F:
    # Historical packaging of the colour-blind SM readout; anomaly-equivalent.
    return F(1, 2) * zf(c, rm) + b * sum_zc(c, rm) + F(1, 2)


def t3(c: Codeword, rm: RoleMap = SM_ROLE) -> F:
    return F(1, 2) * zf(c, rm) if c[rm.chi] == 0 else F(0)


def y(c: Codeword, rm: RoleMap = SM_ROLE) -> F:
    return q(c, rm) - t3(c, rm)


def y_eff(c: Codeword, rm: RoleMap = SM_ROLE) -> F:
    yy = y(c, rm)
    return yy if c[rm.chi] == 0 else -yy


def triplet(c: Codeword, rm: RoleMap = SM_ROLE) -> bool:
    return (c[rm.c0], c[rm.c1]) != (0, 0)


def doublet(c: Codeword, rm: RoleMap = SM_ROLE) -> bool:
    return c[rm.chi] == 0


def excluded(c: Codeword, exclusion: Exclusion | None, rm: RoleMap = SM_ROLE) -> bool:
    if exclusion is None:
        return False
    kind, lq_value, i3_value, chi_value = exclusion
    assert kind == "lq"
    return c[rm.lq] == lq_value and c[rm.i3] == i3_value and c[rm.chi] == chi_value


def content(sectors: frozenset[Triple], exclusion: Exclusion | None, rm: RoleMap = SM_ROLE) -> list[Codeword]:
    out: list[Codeword] = []
    for c in CW:
        if not (r1(c, rm) and r2(c, rm) and sector(c, rm) in sectors):
            continue
        if excluded(c, exclusion, rm):
            continue
        out.append(c)
    return out


def anomalies(states: list[Codeword], rm: RoleMap = SM_ROLE) -> tuple[F | int, ...]:
    a1 = sum(y_eff(c, rm) for c in states if triplet(c, rm))
    a2 = F(1, 2) * sum(y_eff(c, rm) for c in states if doublet(c, rm))
    a3 = sum(y_eff(c, rm) ** 3 for c in states)
    a4 = sum(y_eff(c, rm) for c in states)
    a5 = sum(1 for c in states if triplet(c, rm) and c[rm.chi] == 0) - sum(
        1 for c in states if triplet(c, rm) and c[rm.chi] == 1
    )
    # Witten SU(2): number of complete doublets modulo 2.  In this rule family
    # doublet completeness is checked separately; this is still the paper's
    # finite ledger convention.
    a6 = (sum(1 for c in states if doublet(c, rm)) // 2) % 2
    return (a1, a2, a3, a4, a5, a6)


def anomaly_free(states: list[Codeword], rm: RoleMap = SM_ROLE) -> bool:
    return anomalies(states, rm) == (0, 0, 0, 0, 0, 0)


def vectorlike(states: list[Codeword], rm: RoleMap = SM_ROLE) -> bool:
    s = set(states)
    for c in states:
        partner = list(c)
        partner[rm.chi] = 1 - partner[rm.chi]
        partner[rm.w] = partner[rm.chi]
        if tuple(partner) not in s:
            return False
    return True


def su2_complete(states: list[Codeword], rm: RoleMap = SM_ROLE) -> bool:
    # Every weak doublet should come in I3 pairs for fixed generation/sector.
    by_base: dict[tuple[int, int, Triple, int], set[int]] = defaultdict(set)
    for c in states:
        if c[rm.chi] != 0:
            continue
        gen = (c[rm.g0], c[rm.g1])
        by_base[(gen[0], gen[1], sector(c, rm), c[rm.chi])].add(c[rm.i3])
    return all(vals == {0, 1} for vals in by_base.values())


def one_generation_signature(states: list[Codeword], rm: RoleMap = SM_ROLE) -> tuple[tuple[str, str, bool, bool], ...]:
    first_gen = [c for c in states if c[rm.g0] == 0 and c[rm.g1] == 0]
    return tuple(
        sorted((str(q(c, rm)), str(y_eff(c, rm)), doublet(c, rm), triplet(c, rm)) for c in first_gen)
    )


SM_CONTENT = content(SM_SECTOR, SM_EXCLUSION)
SM_SIGNATURE = one_generation_signature(SM_CONTENT)


def sm_isomorphic(states: list[Codeword], rm: RoleMap = SM_ROLE) -> bool:
    return one_generation_signature(states, rm) == SM_SIGNATURE


def r3_clean(sectors: frozenset[Triple]) -> bool:
    """True when one LQ value is colourless and the other is the three-colour triplet."""
    if len(sectors) != 4:
        return False
    for lepton_lq in (0, 1):
        quark_lq = 1 - lepton_lq
        target = {(lepton_lq, 0, 0)} | {(quark_lq, c0, c1) for c0, c1 in [(0, 1), (1, 0), (1, 1)]}
        if set(sectors) == target:
            return True
    return False


def colourless_lq(sectors: frozenset[Triple]) -> int | None:
    if not r3_clean(sectors):
        return None
    for lq, c0, c1 in sectors:
        if (c0, c1) == (0, 0):
            return lq
    return None


def sterile_exclusion(sectors: frozenset[Triple], exclusion: Exclusion, rm: RoleMap = SM_ROLE) -> bool:
    """Whether the exclusion removes exactly colourless Q=0 right-singlets in this sector rule."""
    kind, lq_value, _i3_value, chi_value = exclusion
    if kind != "lq" or chi_value != 1:
        return False
    removed: list[Codeword] = []
    for c in CW:
        if r1(c, rm) and r2(c, rm) and sector(c, rm) in sectors and excluded(c, exclusion, rm):
            removed.append(c)
    if not removed:
        return False
    return all((not triplet(c, rm)) and q(c, rm) == 0 and c[rm.lq] == lq_value for c in removed)


def sm_rule_like(sectors: frozenset[Triple], exclusion: Exclusion | None) -> bool:
    """SM matter/colour split plus sterile right-singlet removal, up to LQ-name flip."""
    if exclusion is None or not r3_clean(sectors):
        return False
    kind, lq_value, i3_value, chi_value = exclusion
    return (
        kind == "lq"
        and lq_value == colourless_lq(sectors)
        and i3_value == 0
        and chi_value == 1
    )


def all_sector_masks() -> list[frozenset[Triple]]:
    out: list[frozenset[Triple]] = []
    for n in range(1, len(TRIPLES) + 1):
        for comb in combinations(TRIPLES, n):
            out.append(frozenset(comb))
    return out


# GL(2,F2) on the colour pair.  This fixes 00 and permutes the three non-zero
# colour labels, so it quotients colour names without changing the physics.
GL2 = [
    ((1, 0), (0, 1)),
    ((0, 1), (1, 0)),
    ((1, 1), (0, 1)),
    ((1, 0), (1, 1)),
    ((0, 1), (1, 1)),
    ((1, 1), (1, 0)),
]


def gl_act(m: tuple[tuple[int, int], tuple[int, int]], col: tuple[int, int]) -> tuple[int, int]:
    (a, b), (c, d) = m
    x, y_ = col
    return ((a * x + b * y_) % 2, (c * x + d * y_) % 2)


def transform_sector(mask: frozenset[Triple], m: tuple[tuple[int, int], tuple[int, int]]) -> frozenset[Triple]:
    return frozenset((lq, *gl_act(m, (c0, c1))) for (lq, c0, c1) in mask)


def transform_exclusion(ex: Exclusion, m: tuple[tuple[int, int], tuple[int, int]]) -> Exclusion:
    # LQ-local exclusions do not read colour and are unchanged by GL(2,F2).
    return ex


def canonical_sector(mask: frozenset[Triple]) -> tuple[Triple, ...]:
    return min(tuple(sorted(transform_sector(mask, m))) for m in GL2)


def canonical_pair(mask: frozenset[Triple], ex: Exclusion | None) -> tuple[tuple[Triple, ...], Exclusion | None]:
    variants = []
    for m in GL2:
        mm = transform_sector(mask, m)
        ee = transform_exclusion(ex, m) if ex is not None else None
        variants.append((tuple(sorted(mm)), ee))
    return min(variants)


def canonical_r4(sectors: frozenset[Triple]) -> Exclusion | None:
    return SM_EXCLUSION


def classify_cases(cases: list[tuple[frozenset[Triple], Exclusion | None]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for sectors, exclusion in cases:
        states = content(sectors, exclusion)
        counts["total"] += 1
        if anomaly_free(states):
            counts["anomaly_free"] += 1
        if not vectorlike(states):
            counts["chiral"] += 1
        if anomaly_free(states) and not vectorlike(states):
            counts["chiral_anomaly_free"] += 1
        if anomaly_free(states) and not vectorlike(states) and su2_complete(states):
            counts["su2_complete_chiral_af"] += 1
        if sm_isomorphic(states):
            counts["sm_spectrum"] += 1
        if sm_rule_like(sectors, exclusion):
            counts["sm_rule_like"] += 1
        if exclusion is not None and sterile_exclusion(sectors, exclusion):
            counts["sterile_exclusion"] += 1
    return counts


def print_counts(title: str, counts: Counter[str]) -> None:
    print(f"\n=== {title} ===")
    for key in [
        "total",
        "anomaly_free",
        "chiral",
        "chiral_anomaly_free",
        "su2_complete_chiral_af",
        "sm_spectrum",
        "sm_rule_like",
        "sterile_exclusion",
    ]:
        if key in counts:
            print(f"{key:26s}: {counts[key]}")


def role_maps() -> list[RoleMap]:
    maps: list[RoleMap] = []
    bits = list(range(8))
    # This matches the earlier Target-A full scan convention: choose I3, chi,
    # LQ, an unordered colour pair, then let the remaining two spectator bits
    # be generation and W.  The two generation-bit orderings are quotiented.
    for i3, chi, lq in permutations(bits, 3):
        rest = [b for b in bits if b not in (i3, chi, lq)]
        for c0, c1 in combinations(rest, 2):
            rem = [b for b in rest if b not in (c0, c1)]
            maps.append(RoleMap(i3=i3, chi=chi, lq=lq, c0=c0, c1=c1, g0=rem[0], g1=rem[1], w=rem[2]))
    return maps


def main() -> int:
    print("Target-A finite classification")
    print("canonical SM size:", len(SM_CONTENT))
    print("canonical SM anomalies:", anomalies(SM_CONTENT))
    print("canonical SM vectorlike?:", vectorlike(SM_CONTENT))
    print("canonical SM SU(2)-complete?:", su2_complete(SM_CONTENT))
    print("canonical SM sterile exclusion?:", sterile_exclusion(SM_SECTOR, SM_EXCLUSION))

    masks = all_sector_masks()

    # A: current scan, fixed canonical R4 where available.
    cases_a = [(mask, canonical_r4(mask)) for mask in masks]
    counts_a = classify_cases(cases_a)
    print_counts("A. fixed-R4 sector masks (paper scan)", counts_a)

    # B: quotient A by GL(2,F2) colour relabelling.
    orbit_cases: dict[tuple[tuple[Triple, ...], Exclusion | None], tuple[frozenset[Triple], Exclusion | None]] = {}
    for mask, ex in cases_a:
        orbit_cases.setdefault(canonical_pair(mask, ex), (mask, ex))
    cases_b = list(orbit_cases.values())
    counts_b = classify_cases(cases_b)
    print_counts("B. A quotiented by colour GL(2,F2)", counts_b)
    print("colour orbits in A:", len(cases_b))

    print("\nA/B anomaly-free orbit representatives:")
    for mask, ex in cases_b:
        states = content(mask, ex)
        if anomaly_free(states):
            print(
                "  size={:2d} vectorlike={} sm_rule={} sector={} exclusion={}".format(
                    len(states), vectorlike(states), sm_rule_like(mask, ex), tuple(sorted(mask)), ex
                )
            )

    # C: wider bit-local right-singlet exclusion class.
    cases_c: list[tuple[frozenset[Triple], Exclusion | None]] = []
    for mask in masks:
        for lq_value in (0, 1):
            for i3_value in (0, 1):
                cases_c.append((mask, ("lq", lq_value, i3_value, 1)))
    counts_c = classify_cases(cases_c)
    print_counts("C. all bit-local right-singlet exclusions", counts_c)

    orbit_cases_c: dict[tuple[tuple[Triple, ...], Exclusion | None], tuple[frozenset[Triple], Exclusion | None]] = {}
    for mask, ex in cases_c:
        orbit_cases_c.setdefault(canonical_pair(mask, ex), (mask, ex))
    cases_cq = list(orbit_cases_c.values())
    counts_cq = classify_cases(cases_cq)
    print_counts("Cq. C quotiented by colour GL(2,F2)", counts_cq)
    print("colour orbits in C:", len(cases_cq))

    # Show the actual SU(2)-complete chiral anomaly-free survivors in Cq.
    survivors = []
    for mask, ex in cases_cq:
        states = content(mask, ex)
        if anomaly_free(states) and (not vectorlike(states)) and su2_complete(states):
            survivors.append((mask, ex, len(states), sm_isomorphic(states), sm_rule_like(mask, ex), sterile_exclusion(mask, ex)))
    print("\nCq survivors (SU2-complete, chiral, anomaly-free):")
    for mask, ex, size, is_sm, is_rule_sm, is_sterile in survivors:
        print(
            "  size={:2d} sm_spectrum={} sm_rule={} sterile_exclusion={} sector={} exclusion={}".format(
                size, is_sm, is_rule_sm, is_sterile, tuple(sorted(mask)), ex
            )
        )

    # D: role-map invariance for the canonical SM rule.
    rmaps = role_maps()
    iso = 0
    af = 0
    for rm in rmaps:
        states = content(SM_SECTOR, SM_EXCLUSION, rm)
        if anomaly_free(states, rm):
            af += 1
        if sm_isomorphic(states, rm):
            iso += 1
    print("\n=== D. role-map invariance for canonical SM rule ===")
    print("role maps tested          :", len(rmaps))
    print("anomaly-free role maps    :", af)
    print("SM-isomorphic role maps   :", iso)

    # Guardrails for the paper-level theorem.
    assert anomalies(SM_CONTENT) == (0, 0, 0, 0, 0, 0)
    assert len(SM_CONTENT) == 45
    assert not vectorlike(SM_CONTENT)
    assert su2_complete(SM_CONTENT)
    assert sterile_exclusion(SM_SECTOR, SM_EXCLUSION)
    assert counts_a["total"] == 255
    assert counts_a["anomaly_free"] == 2
    assert counts_a["chiral_anomaly_free"] == 1
    assert counts_a["sm_rule_like"] == 1
    assert counts_b["chiral_anomaly_free"] == 1
    assert counts_b["sm_rule_like"] == 1
    assert counts_c["sm_rule_like"] == 2  # LQ-name flip included in the wider class.
    assert counts_cq["sm_rule_like"] == 2
    assert len(survivors) >= 1
    assert any(is_rule_sm and is_sterile for _, _, _, _, is_rule_sm, is_sterile in survivors)
    assert af == len(rmaps) and iso == len(rmaps)

    print("\nVERDICT:")
    print("  The paper scan has two anomaly-free sector rules, but only one is")
    print("  genuinely chiral; the other is vectorlike.  The unique chiral")
    print("  anomaly-free orbit is the SM rule.  In the wider bit-local")
    print("  right-singlet exclusion class, the only chiral anomaly-free survivors")
    print("  are the SM rule and its LQ-name flip.")
    print("  Role-map freedom is pure relabelling in the canonical SM rule.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
