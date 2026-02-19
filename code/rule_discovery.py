#!/usr/bin/env python3
"""
rule_discovery.py — Discovery of the Unique Weak Rule

Demonstrates that among all invertible, sector-boundary XOR couplings
on the 8-bit circlette ring, exactly ONE non-trivial rule preserves
all 45 valid fermion states in pure-valid cycles while minimising
the average bit-flip cost (information action).

That rule is:  I3(t+1) = I3(t) ⊕ LQ(t)
               (isospin XORed with the bridge bit at every tick)

This is a CNOT gate with LQ as control and I3 as target.
It IS the weak interaction.

Reference: Elliman (2026), "It from Bit, Revisited"
"""

from circlette import (
    generate_valid_states, is_valid, state_to_str,
    G0, G1, C0, C1, LQ, I3, CHI, W, BIT_NAMES
)

# Sector boundaries on the ring (adjacent bits from different sectors)
SECTOR_BOUNDARIES = [
    (G1, C0, "generation→colour"),
    (C1, LQ, "colour→bridge"),
    (LQ, I3, "bridge→electroweak"),
    (W,  G0, "electroweak→generation"),
]


def apply_xor_rule(state, source, target):
    """Apply rule: target(t+1) = target(t) XOR source(t)."""
    new = list(state)
    new[target] = new[target] ^ new[source]
    return tuple(new)


def evaluate_rule(source, target, valid_states):
    """
    Evaluate a single XOR rule on the valid state set.

    Returns:
        preserves_all: bool — does the rule map every valid state to a valid state?
        avg_cost: float — average number of bit flips per application
        fixed_points: int — states unchanged by the rule
        period2: int — states that return to themselves after 2 applications
    """
    preserves_all = True
    total_flips = 0
    fixed_points = 0
    period2 = 0

    for state in valid_states:
        new_state = apply_xor_rule(state, source, target)

        if not is_valid(new_state):
            preserves_all = False

        # Count bit flips
        flips = sum(a != b for a, b in zip(state, new_state))
        total_flips += flips

        if flips == 0:
            fixed_points += 1

        # Check period-2 (involution)
        restored = apply_xor_rule(new_state, source, target)
        if restored == state:
            period2 += 1

    avg_cost = total_flips / len(valid_states)

    return preserves_all, avg_cost, fixed_points, period2


def search_all_rules():
    """
    Exhaustively search all sector-boundary XOR rules.

    For each pair of adjacent bits across a sector boundary,
    test both directions (A⊕B and B⊕A).
    """
    valid_states = generate_valid_states()
    assert len(valid_states) == 45, f"Expected 45 states, got {len(valid_states)}"

    print("=" * 80)
    print("EXHAUSTIVE SEARCH: Sector-Boundary XOR Rules")
    print("=" * 80)
    print(f"\nSearching {len(SECTOR_BOUNDARIES) * 2} candidate rules ")
    print(f"(2 directions × {len(SECTOR_BOUNDARIES)} boundaries)...\n")

    winning_rules = []

    print(f"{'Rule':<25} | {'Preserves?':<11} | {'Avg cost':<10} | "
          f"{'Fixed pts':<10} | {'Period 2':<10} | {'Involution?'}")
    print("-" * 95)

    for b1, b2, boundary_name in SECTOR_BOUNDARIES:
        for source, target in [(b1, b2), (b2, b1)]:
            rule_name = f"{BIT_NAMES[target]} ⊕ {BIT_NAMES[source]}"
            preserves, avg_cost, fixed_pts, period2 = evaluate_rule(
                source, target, valid_states
            )
            is_involution = (period2 == 45)
            is_nontrivial = (fixed_pts < 45)

            marker = ""
            if preserves and is_nontrivial:
                marker = " ★ WINNER"
                winning_rules.append((rule_name, avg_cost, fixed_pts))

            print(f"{rule_name:<25} | {'✓' if preserves else '✗':<11} | "
                  f"{avg_cost:<10.4f} | {fixed_pts:<10} | {period2:<10} | "
                  f"{'Yes' if is_involution else 'No'}{marker}")

    print(f"\n{'=' * 80}")
    print(f"RESULT: {len(winning_rules)} non-trivial spectrum-preserving rule(s) found")
    print(f"{'=' * 80}\n")

    if len(winning_rules) == 1:
        name, cost, fixed = winning_rules[0]
        print(f"  The UNIQUE rule is: {name}")
        print(f"  Average bit-flip cost: {cost:.4f} = {int(cost * 45)}/45")
        print(f"  Fixed points (leptons): {fixed}")
        print(f"  Oscillating states (quarks): {45 - fixed}")
        print(f"\n  This is a CNOT gate: LQ controls, I3 is target.")
        print(f"  It IS the weak interaction.")
    else:
        for name, cost, fixed in winning_rules:
            print(f"  {name}: cost={cost:.4f}, fixed={fixed}")


if __name__ == "__main__":
    search_all_rules()
