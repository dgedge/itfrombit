#!/usr/bin/env python3
"""ITEM 131 / 42: can actual R2/R3 geometry derive N_print alpha^4 = 4/3?

Question
--------
The early-printer threshold audit reduced the H_* problem to

    N_print(H_*) * alpha_0^4 = 4/3.

Can that equation be derived from the actual R2/R3 Kraus/current geometry?

Verdict
-------
Not yet.  The pieces are real, but current canon does not connect them.

What the actual R2/R3 geometry DOES derive:
  * R3 exposes the SU(3) colour triplet, whose fundamental quadratic Casimir is
    C_F = 4/3 under standard generator normalisation.
  * The portal-licensed [8,4,4] logical sector supplies alpha^4 per weight-4
    logical commit.

What it does NOT derive:
  * R2/R3 Boolean enforcement itself is not quartic.  On the actual 48-state
    register, one-bit flips can violate R2 or R3.  The alpha^4 power therefore
    comes from importing the weight-4 logical-Kraus sector, not from R2/R3 alone.
  * The R2/R3 Pati-Salam broken-coset quadratic loads are not 4/3.  The SU(4)
    broken-coset average is 3/4, the right-SU(2) broken load is 1/2 on the
    right sector (1/4 averaged over the spinor), and the full PS broken average
    over the 16 is 1.  The 4/3 is the unbroken SU(3) quark-triplet Casimir.
  * No current conservation / stopping rule in canon says "activate when the
    expected quartic logical commits per printed shell equals C_F."

So the live derivation target becomes two clauses:

  (A) logical-activation lemma:
      R2/R3 threshold current is the weight-4 portal logical current, not the
      lower-weight Boolean violation current.

  (B) Casimir-load criterion:
      the stopping load is the SU(3) fundamental C_F, not a broken-coset PS
      Casimir or another nearby simple invariant.
"""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction as F


G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def bits8(n: int) -> tuple[int, ...]:
    return tuple((n >> (7 - i)) & 1 for i in range(8))


def idx(c: tuple[int, ...]) -> int:
    return int("".join(str(x) for x in c), 2)


def flip(c: tuple[int, ...], support: set[int] | frozenset[int]) -> tuple[int, ...]:
    out = list(c)
    for bit in support:
        out[bit] ^= 1
    return tuple(out)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] and c[G1])


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    return (c[C0], c[C1]) == (0, 0) if c[LQ] == 0 else (c[C0], c[C1]) != (0, 0)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def violation_tag(c: tuple[int, ...]) -> str:
    tags = []
    if not r1(c):
        tags.append("R1")
    if not r2(c):
        tags.append("R2")
    if not r3(c):
        tags.append("R3")
    return "+".join(tags) if tags else "P"


def ag3_bits(x: int) -> tuple[int, int, int]:
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


def hyperplanes() -> list[frozenset[int]]:
    out = set()
    for a in range(1, 8):
        avec = ag3_bits(a)
        for b in (0, 1):
            out.add(
                frozenset(
                    p
                    for p in range(8)
                    if sum(u * v for u, v in zip(ag3_bits(p), avec)) % 2 == b
                )
            )
    return sorted(out, key=lambda h: tuple(sorted(h)))


def su_n_fundamental_casimir(n: int) -> F:
    return F(n * n - 1, 2 * n)


def main() -> None:
    print("ITEM 131 / 42 R2/R3 THRESHOLD LEDGER DERIVATION ATTEMPT")

    physical = [bits8(n) for n in range(256) if valid_r123(bits8(n))]
    one_gen = [c for c in physical if (c[G0], c[G1]) == (0, 0)]
    check(len(physical) == 48, "actual R1/R2/R3 physical register has 48 states")
    check(len(one_gen) == 16, "one generation has the SO(10)/Pati-Salam 16 states")

    print("\n[1] Boolean R2/R3 enforcement is lower-weight, not quartic")
    min_bad_weight = 9
    bad_by_weight: dict[int, Counter[str]] = {}
    for c in physical:
        for weight in range(1, 5):
            for support in itertools.combinations(range(8), weight):
                tag = violation_tag(flip(c, set(support)))
                if "R2" in tag or "R3" in tag:
                    min_bad_weight = min(min_bad_weight, weight)
                    bad_by_weight.setdefault(weight, Counter())[tag] += 1
            if min_bad_weight < 9:
                break
    print(f"  minimum R2/R3-violating Boolean fault weight = {min_bad_weight}")
    for weight in sorted(bad_by_weight):
        print(f"  weight {weight}: {dict(bad_by_weight[weight])}")
    check(min_bad_weight == 1, "single-bit flips can violate R2 or R3")
    print("  Consequence: R2/R3 alone does not force alpha^4.  The alpha^4 power")
    print("  must come from the separate [8,4,4] weight-4 logical-Kraus sector.")

    print("\n[2] Weight-4 logical sector interaction with R2/R3")
    hyps = hyperplanes()
    check(len(hyps) == 14 and all(len(h) == 4 for h in hyps), "[8,4,4] logical sector has 14 weight-4 hyperplanes")
    per_channel = []
    globally_preserves = 0
    for h in hyps:
        tags = Counter(violation_tag(flip(c, h)) for c in physical)
        if set(tags) == {"P"}:
            globally_preserves += 1
        per_channel.append((h, tags))
    print(f"  channels preserving all 48 R1/R2/R3 states = {globally_preserves}/14")
    for h, tags in per_channel:
        support = "{" + ",".join(NAMES[i] for i in sorted(h)) + "}"
        touches_colour = bool(h & {C0, C1})
        print(f"  {support:28s} touches_colour={touches_colour!s:5s}  {dict(tags)}")
    check(globally_preserves == 0, "no weight-4 hyperplane is simply the R2/R3 activation map on all states")
    check(any(not (h & {C0, C1}) for h in hyps), "three colour-silent hyperplanes exist, but they are a filtered subset")
    print("  Consequence: the weight-4 logical current is compatible with R2/R3 bookkeeping,")
    print("  but the actual R2/R3 action table does not by itself select the whole threshold current.")

    print("\n[3] Casimir accounting from the actual PS/SM representation")
    c_su3 = su_n_fundamental_casimir(3)
    c_su4 = su_n_fundamental_casimir(4)
    # SU(4) -> SU(3) x U(1)_{B-L}; T15 = diag(1,1,1,-3)/(2 sqrt(6)).
    u1_quark = F(1, 24)
    u1_lepton = F(3, 8)
    su4_broken_quark = c_su4 - c_su3 - u1_quark
    su4_broken_lepton = c_su4 - u1_lepton
    su4_broken_average = (3 * su4_broken_quark + su4_broken_lepton) / 4
    # SU(2)_R -> U(1)_{T3R}; right doublet broken load = C2(2)-T3R^2.
    c_su2 = su_n_fundamental_casimir(2)
    su2r_broken_right = c_su2 - F(1, 4)
    su2r_broken_spinor_average = su2r_broken_right / 2
    ps_broken_spinor_average = su4_broken_average + su2r_broken_spinor_average
    q_left_sm_casimir = c_su3 + c_su2 + F(1, 36)
    rows = [
        ("SU(3) fundamental C_F", c_su3),
        ("SU(4) fundamental C2", c_su4),
        ("SU(4)/SU(3)xU1 broken, quark", su4_broken_quark),
        ("SU(4)/SU(3)xU1 broken, lepton", su4_broken_lepton),
        ("SU(4) broken average over 4", su4_broken_average),
        ("SU(2)_R broken on right doublet", su2r_broken_right),
        ("SU(2)_R broken average over 16", su2r_broken_spinor_average),
        ("full PS broken average over 16", ps_broken_spinor_average),
        ("left quark SM Casimir C3+C2+Y^2", q_left_sm_casimir),
    ]
    for label, value in rows:
        print(f"  {label:38s} = {value} = {float(value):.9f}")
    check(c_su3 == F(4, 3), "R3 colour triplet gives SU(3) C_F=4/3")
    check(ps_broken_spinor_average == F(1, 1), "actual PS broken-coset average is 1, not 4/3")
    check(su4_broken_average == F(3, 4), "R3 broken SU(4) coset average is 3/4, not 4/3")
    print("  Consequence: C_F=4/3 is a real R3/SU(3) representation invariant,")
    print("  but it is not the quadratic load of the broken R2/R3 coset itself.")

    print("\n[4] Gate verdict")
    gates = [
        ("alpha^4 power", "OPEN", "requires logical-activation lemma: R2/R3 threshold samples weight-4 portal logical commits rather than lower-weight Boolean faults"),
        ("4/3 coefficient", "PARTIAL", "R3 supplies SU(3) fundamental C_F=4/3, but broken-coset loads differ"),
        ("threshold rule", "OPEN", "no current law says activation exits at expected quartic logical load = C_F"),
        ("scale bridge", "CONDITIONAL", "if the three clauses above hold, N_print=S_dS/ln2 gives H*=1.199e15 GeV"),
    ]
    for name, status, note in gates:
        print(f"  {name:18s} {status:12s} {note}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  I do not see a full derivation of N_print alpha_0^4 = 4/3 from the")
    print("  actual R2/R3 Kraus/current geometry yet.  The route is now precise:")
    print("")
    print("    1. derive a logical-activation lemma forcing the R2/R3 threshold")
    print("       current into the portal-licensed [8,4,4] weight-4 sector;")
    print("    2. derive a Casimir-load criterion selecting the R3 quark-triplet")
    print("       SU(3) C_F=4/3, rather than the PS broken-coset loads;")
    print("    3. derive the stopping rule E[N_logical]=C_F per printed shell.")
    print("")
    print("  Current canon proves the ingredients separately, but not their product.")
    print("  Without (1), R2/R3 admits one-bit violation currents and alpha^4 is")
    print("  imported.  Without (2), the actual broken-coset geometry gives 3/4,")
    print("  1/2, or 1-class loads, not 4/3.  Without (3), the threshold equation")
    print("  is still the alpha^4 amplitude candidate rewritten as a printer-exit law.")
    print("=" * 108)
    print("exit 0 -- R2/R3 route sharpened; no independent derivation of N_print alpha^4=4/3.")


if __name__ == "__main__":
    main()
