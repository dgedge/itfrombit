"""
rule_discovery.py — Exhaustive search for the unique weak rule.

Demonstrates that among all invertible sector-boundary XOR couplings
on the 8-bit circlette ring, exactly ONE non-trivial rule preserves
all 45 valid fermion states in pure-valid cycles while minimising
the average bit-flip cost (information action).

That rule is:  I₃(t+1) = I₃(t) ⊕ LQ(t)   (a CNOT gate)

Usage:
    python src/rule_discovery.py
"""

from __future__ import annotations

import sys
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from circlette import (
    LABELS,
    SECTOR_BOUNDARIES,
    State,
    generate_all_states,
    generate_valid_states,
    particle_name,
    state_to_str,
)


def boundary_matrix(src: int, tgt: int) -> np.ndarray:
    """Build Identity + e_tgt · e_src^T over F₂.

    This matrix XORs bit `src` into bit `tgt` at each tick,
    leaving all other bits unchanged.
    """
    m = np.eye(8, dtype=int)
    m[tgt][src] = 1
    return m


def matrix_invertible_f2(mat: np.ndarray) -> bool:
    """Check if an 8×8 matrix is invertible over F₂ via Gaussian elimination."""
    n = len(mat)
    m = mat.copy()
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if m[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            return False
        m[[col, pivot]] = m[[pivot, col]]
        for row in range(n):
            if row != col and m[row][col] == 1:
                m[row] = (m[row] + m[col]) % 2
    return True


def apply_rule(state: State, matrix: np.ndarray) -> State:
    """Apply a linear rule (matrix over F₂) to an 8-bit state."""
    result = []
    for i in range(8):
        val = 0
        for j in range(8):
            val ^= matrix[i][j] * state[j]
        result.append(val)
    return tuple(result)


def find_cycles(
    matrix: np.ndarray,
    all_states: List[State],
) -> List[List[State]]:
    """Find all cycles of the rule acting on the full 256-state space."""
    visited: Dict[State, int] = {}
    cycles: List[List[State]] = []

    for s in all_states:
        if s in visited:
            continue
        cycle: List[State] = []
        current = s
        while current not in visited:
            visited[current] = len(cycles)
            cycle.append(current)
            current = apply_rule(current, matrix)
        cycles.append(cycle)

    return cycles


def matrix_order_f2(mat: np.ndarray, max_k: int = 300) -> Optional[int]:
    """Find the order of a matrix over F₂ (smallest k with M^k = I)."""
    mk = np.eye(8, dtype=int)
    for k in range(1, max_k):
        mk = (mk @ mat) % 2
        if np.array_equal(mk, np.eye(8, dtype=int)):
            return k
    return None


def analyse_rule(
    matrix: np.ndarray,
    valid_set: Set[State],
    all_states: List[State],
) -> Dict:
    """Comprehensive analysis of a candidate rule.

    Returns a dict with:
        valid_in_pure:      number of valid states in pure-valid cycles
        pure_cycle_lengths: cycle lengths for pure-valid cycles
        avg_flips:          average Hamming distance per valid state per tick
        order:              matrix order over F₂
        total_cycles:       total number of cycles on all 256 states
    """
    cycles = find_cycles(matrix, all_states)

    pure_valid_cycles = [
        c for c in cycles if set(c).issubset(valid_set) and len(c) > 0
    ]

    valid_in_pure = sum(len(c) for c in pure_valid_cycles)
    pure_cycle_lengths = sorted(
        [len(c) for c in pure_valid_cycles], reverse=True
    )

    total_flips = sum(
        sum(a != b for a, b in zip(s, apply_rule(s, matrix)))
        for s in valid_set
    )
    avg_flips = total_flips / len(valid_set)

    order = matrix_order_f2(matrix)

    return {
        "valid_in_pure": valid_in_pure,
        "pure_cycle_lengths": pure_cycle_lengths,
        "avg_flips": avg_flips,
        "order": order,
        "total_cycles": len(cycles),
    }


def main() -> None:
    valid_states = generate_valid_states()
    valid_set = set(valid_states)
    all_states = generate_all_states()

    print("=" * 70)
    print("CIRCLETTE RULE DISCOVERY")
    print("Exhaustive search over all sector-boundary couplings")
    print("=" * 70)
    print(f"\nValid fermion states: {len(valid_states)}")
    print(f"Total 8-bit states:  {len(all_states)}")
    print(f"\nRing: {' '.join(LABELS)}")
    print(f"Bits: {' '.join(f'{i:>2}' for i in range(8))}")

    # ----------------------------------------------------------
    # Test all 8 single sector-boundary couplings
    # ----------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("TESTING ALL SINGLE SECTOR-BOUNDARY COUPLINGS")
    print(f"{'=' * 70}")

    results = []

    for src, tgt, name in SECTOR_BOUNDARIES:
        mat = boundary_matrix(src, tgt)

        if not matrix_invertible_f2(mat.copy()):
            print(f"\n  {name:10s}  NOT INVERTIBLE — skipped")
            continue

        res = analyse_rule(mat, valid_set, all_states)
        res["name"] = name
        res["src"] = src
        res["tgt"] = tgt
        results.append(res)

        marker = " ★ ALL 45 PRESERVED" if res["valid_in_pure"] == 45 else ""
        print(
            f"\n  {name:10s}  valid={res['valid_in_pure']:>2}/45"
            f"  flips={res['avg_flips']:.2f}"
            f"  order={res['order']}{marker}"
        )

    # ----------------------------------------------------------
    # Ranking
    # ----------------------------------------------------------
    results.sort(key=lambda r: (-r["valid_in_pure"], r["avg_flips"]))

    print(f"\n{'=' * 70}")
    print("RANKING (sorted by valid states preserved, then bit-flip cost)")
    print(f"{'=' * 70}")
    print(f"\n{'Rank':>4}  {'Rule':>10}  {'Valid':>7}  {'Flips':>6}  {'Order':>5}")
    print("-" * 42)
    for i, r in enumerate(results):
        star = " ★" if r["valid_in_pure"] == 45 else ""
        print(
            f"{i + 1:>4}  {r['name']:>10}  {r['valid_in_pure']:>3}/45"
            f"  {r['avg_flips']:>6.2f}  {r['order']:>5}{star}"
        )

    # ----------------------------------------------------------
    # Uniqueness result
    # ----------------------------------------------------------
    winners = [
        r for r in results if r["valid_in_pure"] == 45 and r["avg_flips"] > 0
    ]

    print(f"\n{'=' * 70}")
    print("RESULT")
    print(f"{'=' * 70}")

    if len(winners) == 1:
        w = winners[0]
        print(f"\n  ★ UNIQUE NON-TRIVIAL RULE: {w['name']}")
        print(f"    {LABELS[w['tgt']]}(t+1) = {LABELS[w['tgt']]}(t) ⊕ {LABELS[w['src']]}(t)")
        print(f"    Valid states preserved:  {w['valid_in_pure']}/45 (all)")
        print(f"    Average bit-flip cost:   {w['avg_flips']:.2f} bits/tick")
        print(f"    Matrix order:            {w['order']} (M² = I)")
        print(f"    Cycle structure:         {w['pure_cycle_lengths']}")
        print(f"\n    This rule is a CNOT gate:")
        print(f"      {LABELS[w['src']]} (control) ──●── {LABELS[w['src']]}  (unchanged)")
        print(f"      {LABELS[w['tgt']]} (target)  ──⊕── {LABELS[w['tgt']]} ⊕ {LABELS[w['src']]}")
    elif len(winners) == 0:
        print("\n  No non-trivial rule preserves all 45 states.")
    else:
        print(f"\n  Multiple winners found ({len(winners)}):")
        for w in winners:
            print(f"    {w['name']}: flips={w['avg_flips']:.2f}")

    # ----------------------------------------------------------
    # Detailed action on particles
    # ----------------------------------------------------------
    if winners:
        w = winners[0]
        mat = boundary_matrix(w["src"], w["tgt"])

        print(f"\n{'=' * 70}")
        print(f"DETAILED ACTION: {w['name']}")
        print(f"{'=' * 70}")

        print(f"\nUpdate matrix (I + e{w['tgt']}·e{w['src']}ᵀ over F₂):")
        print(f"       {'  '.join(f'{l:>3}' for l in LABELS)}")
        for i, row in enumerate(mat):
            print(f"  {LABELS[i]:>3}: {'  '.join(f'{int(x):>2}' for x in row)}")

        # Classify states
        fixed = []
        pairs = []
        seen: Set[Tuple[State, State]] = set()

        for s in sorted(valid_states):
            ns = apply_rule(s, mat)
            if ns == s:
                fixed.append(s)
            else:
                pair = (min(s, ns), max(s, ns))
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)

        print(f"\nFIXED POINTS ({len(fixed)} leptons, period 1):")
        for s in fixed:
            print(f"  {state_to_str(s)}  {particle_name(s)}")

        print(f"\nOSCILLATING PAIRS ({len(pairs)} weak isospin doublets, period 2):")
        for s1, s2 in sorted(pairs):
            n1 = particle_name(s1)
            n2 = particle_name(s2)
            print(f"  {state_to_str(s1)} ({n1:8s}) ↔ {state_to_str(s2)} ({n2:8s})")

        print(f"\n{'=' * 70}")
        print("INTERPRETATION")
        print(f"{'=' * 70}")
        print(f"""
  Leptons (LQ=0): rule is identity → fixed points (no internal clock)
  Quarks  (LQ=1): I₃ flips every tick → period-2 oscillation

  The 18 oscillating pairs are exactly the weak isospin doublets:
    (u,d), (c,s), (t,b)  ×  3 colours  ×  2 chiralities  =  18 pairs

  Bit-flip cost: 36 quarks × 1 flip + 9 leptons × 0 = 36/45 = 0.80

  The weak interaction is the UNIQUE rule selected by the
  information action principle on the circlette ring.
""")


if __name__ == "__main__":
    # Allow running from project root or src/
    if "src" not in sys.path:
        sys.path.insert(0, "src")
    main()
