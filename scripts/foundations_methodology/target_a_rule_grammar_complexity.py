#!/usr/bin/env python3
"""Small-grammar complexity audit for Target A sector rules.

This script makes the "simplest natural rule" claim executable.  It does not
try to prove absolute Kolmogorov minimality.  Instead it fixes a small Boolean
grammar over the local sector records and computes the exact minimal AST size
for every sector mask on (LQ,C0,C1).

Two grammars are compared:

  bit-literal:
      atoms are equality predicates on LQ, C0, C1, plus TRUE.

  record-local:
      the bit-literal atoms plus primitive colour-plane predicates C=00 and
      C!=00.  These are allowed because the framework treats the two colour
      bits as one local record plane with a singlet plus a triplet.

Cost convention: every atom and every Boolean connective (NOT, AND, OR) costs
one node.  The output records where the anomaly-free and chiral anomaly-free
sector rules sit by this description length.

Exit 0 means the finite grammar audit reproduces the printed claims.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import inf

import target_a_classification as target_a


Triple = target_a.Triple
Mask = frozenset[Triple]
UNIVERSE: Mask = frozenset(target_a.TRIPLES)
SM_SECTOR: Mask = target_a.SM_SECTOR


@dataclass(frozen=True)
class Expr:
    cost: int
    text: str


def atom_masks(kind: str) -> dict[Mask, Expr]:
    atoms: list[tuple[str, Mask]] = [
        ("TRUE", UNIVERSE),
        ("LQ=0", frozenset(t for t in UNIVERSE if t[0] == 0)),
        ("LQ=1", frozenset(t for t in UNIVERSE if t[0] == 1)),
        ("C0=0", frozenset(t for t in UNIVERSE if t[1] == 0)),
        ("C0=1", frozenset(t for t in UNIVERSE if t[1] == 1)),
        ("C1=0", frozenset(t for t in UNIVERSE if t[2] == 0)),
        ("C1=1", frozenset(t for t in UNIVERSE if t[2] == 1)),
    ]
    if kind == "record-local":
        atoms.extend(
            [
                ("C=00", frozenset(t for t in UNIVERSE if (t[1], t[2]) == (0, 0))),
                ("C!=00", frozenset(t for t in UNIVERSE if (t[1], t[2]) != (0, 0))),
            ]
        )
    elif kind != "bit-literal":
        raise ValueError(f"unknown grammar kind: {kind}")

    best: dict[Mask, Expr] = {}
    for text, mask in atoms:
        update(best, mask, Expr(1, text))
    return best


def update(best: dict[Mask, Expr], mask: Mask, expr: Expr) -> bool:
    old = best.get(mask)
    if old is None or (expr.cost, len(expr.text), expr.text) < (old.cost, len(old.text), old.text):
        best[mask] = expr
        return True
    return False


def synthesize(kind: str, max_cost: int = 21) -> dict[Mask, Expr]:
    """Exact finite synthesis up to max_cost over the chosen grammar."""

    best = atom_masks(kind)
    changed = True
    while changed:
        changed = False
        items = list(best.items())

        for mask, expr in items:
            c = UNIVERSE - mask
            new = Expr(expr.cost + 1, f"not({expr.text})")
            if new.cost <= max_cost:
                changed |= update(best, c, new)

        items = list(best.items())
        for m1, e1 in items:
            for m2, e2 in items:
                # Canonicalize textual order to suppress irrelevant tie noise.
                if e2.text < e1.text:
                    continue
                cost = e1.cost + e2.cost + 1
                if cost > max_cost:
                    continue
                changed |= update(best, m1 & m2, Expr(cost, f"({e1.text} and {e2.text})"))
                changed |= update(best, m1 | m2, Expr(cost, f"({e1.text} or {e2.text})"))

        if len(best) == 256 and all(expr.cost < max_cost for expr in best.values()):
            # A fixed point below max_cost has all masks; one more pass would
            # only seek tie improvements that do not affect the audit.
            pass

    if len(best) != 256:
        missing = 256 - len(best)
        raise RuntimeError(f"{kind} grammar failed to synthesize {missing} masks <= {max_cost}")
    return best


def fixed_r4_cases() -> list[tuple[Mask, target_a.Exclusion | None, list[target_a.Codeword]]]:
    cases = []
    for mask in target_a.all_sector_masks():
        exclusion = target_a.canonical_r4(mask)
        states = target_a.content(mask, exclusion)
        cases.append((mask, exclusion, states))
    return cases


def wide_r4_cases() -> list[tuple[Mask, target_a.Exclusion, list[target_a.Codeword]]]:
    cases = []
    exclusions = [("lq", lq, i3, 1) for lq in (0, 1) for i3 in (0, 1)]
    for mask in target_a.all_sector_masks():
        for exclusion in exclusions:
            states = target_a.content(mask, exclusion)
            cases.append((mask, exclusion, states))
    return cases


def case_flags(states: list[target_a.Codeword]) -> tuple[bool, bool, bool]:
    af = target_a.anomaly_free(states)
    chiral = not target_a.vectorlike(states)
    su2 = target_a.su2_complete(states)
    return af, chiral, su2


def short(mask: Mask) -> str:
    return "{" + ",".join("".join(map(str, t)) for t in sorted(mask)) + "}"


def summarize_grammar(kind: str, best: dict[Mask, Expr]) -> None:
    sm = best[SM_SECTOR]
    fixed = fixed_r4_cases()
    wide = wide_r4_cases()

    fixed_af = []
    fixed_chiral_af = []
    for mask, exclusion, states in fixed:
        af, chiral, su2 = case_flags(states)
        if af:
            fixed_af.append((mask, exclusion, states))
        if af and chiral and su2:
            fixed_chiral_af.append((mask, exclusion, states))

    wide_chiral_af = []
    for mask, exclusion, states in wide:
        af, chiral, su2 = case_flags(states)
        if af and chiral and su2:
            wide_chiral_af.append((mask, exclusion, states))

    cheaper_masks = [m for m, e in best.items() if m and e.cost < sm.cost]
    equal_masks = [m for m, e in best.items() if m and e.cost == sm.cost]

    print(f"\n=== {kind} grammar ===")
    print(f"SM sector cost      : {sm.cost}")
    print(f"SM sector formula   : {sm.text}")
    print(f"nonempty masks cheaper than SM: {len(cheaper_masks)}")
    print(f"nonempty masks equal to SM    : {len(equal_masks)}")

    print("\nfixed-R4 anomaly-free sector rules:")
    for mask, _exclusion, states in sorted(fixed_af, key=lambda row: (best[row[0]].cost, short(row[0]))):
        af, chiral, su2 = case_flags(states)
        print(
            f"  cost={best[mask].cost:2d} chiral={chiral!s:5s} su2={su2!s:5s} "
            f"sm_rule={target_a.sm_rule_like(mask, target_a.SM_EXCLUSION)!s:5s} mask={short(mask)}"
        )

    print("\nwide-R4 chiral SU2 anomaly-free survivors:")
    for mask, exclusion, _states in sorted(wide_chiral_af, key=lambda row: (best[row[0]].cost, short(row[0]), row[1])):
        print(
            f"  cost={best[mask].cost:2d} sm_rule={target_a.sm_rule_like(mask, exclusion)!s:5s} "
            f"exclusion={exclusion} mask={short(mask)}"
        )

    assert len(fixed_af) == 2
    assert len(fixed_chiral_af) == 1
    assert fixed_chiral_af[0][0] == SM_SECTOR
    assert len(wide_chiral_af) == 2
    assert all(target_a.sm_rule_like(mask, exclusion) for mask, exclusion, _states in wide_chiral_af)
    assert all(best[mask].cost == sm.cost for mask, _exclusion, _states in wide_chiral_af)
    assert any((target_a.anomaly_free(states) and target_a.vectorlike(states) and best[mask].cost < sm.cost) for mask, _exclusion, states in fixed)


def main() -> int:
    print("Target-A small Boolean grammar complexity audit")
    print("cost = atom/connective AST node count; grammar is finite over sector triples")

    for kind in ("bit-literal", "record-local"):
        best = synthesize(kind)
        summarize_grammar(kind, best)

    print("\nVERDICT:")
    print("  Complexity alone does not select the SM: a simpler vectorlike")
    print("  anomaly-free sector exists.  But once chirality, SU(2) completeness,")
    print("  and anomaly cancellation are imposed inside the natural local rule")
    print("  class, the SM sector is the unique fixed-R4 survivor.  In the wider")
    print("  bit-local sterile-exclusion class, the only survivors have the same")
    print("  minimal grammar cost and are the SM rule plus its LQ-name flip.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
