#!/usr/bin/env python3
"""
circlette.py — Core encoding, constraints, and state generation.

Defines the 8-bit ring structure, the four local constraints (R1–R4),
and generates the 45 valid matter fermion states of the Standard Model.

Canonical ring layout (from the paper):

    Index:   0   1   2   3   4   5   6   7
    Bit:    G0  G1  C0  C1  LQ  I3  χ   W
    Sector: [generation] [colour] [bridge] [electroweak]

This ordering achieves maximum constraint locality: all four rules
span windows of at most 3 adjacent bits on the ring.

Reference: Elliman (2026), "It from Bit, Revisited"
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Tuple

# Type alias for an 8-bit state
State = Tuple[int, ...]

# ============================================================
# Canonical bit positions on the ring
# ============================================================
G0, G1, C0, C1, LQ, I3, CHI, W = range(8)

BIT_NAMES = ["G0", "G1", "C0", "C1", "LQ", "I3", "χ", "W"]

SECTOR_MAP = {
    "generation":  (G0, G1),
    "colour":      (C0, C1),
    "bridge":      (LQ,),
    "electroweak": (I3, CHI, W),
}


# ============================================================
# Constraint checking
# ============================================================

def check_constraints(bits: State) -> Tuple[bool, bool, bool, bool]:
    """
    Apply the four local constraints and return (r1, r2, r3, r4).

    R1: G0,G1 ≠ (1,1)             — three generations only
    R2: χ = W                      — chirality gate
    R3: C1=1 ⟹ LQ=1              — colour requires quark identity
        (equivalently: leptons colourless, quarks carry colour)
    R4: LQ=0 ∧ I3=0 ⟹ χ=0       — no right-handed neutrinos
    """
    g0, g1, c0, c1, lq, i3, chi, w = bits

    r1 = not (g0 == 1 and g1 == 1)
    r2 = (chi == w)

    # R3: leptons must be colourless; quarks must carry colour
    if lq == 0:
        r3 = (c0 == 0 and c1 == 0)
    else:
        r3 = not (c0 == 0 and c1 == 0)

    r4 = not (lq == 0 and i3 == 0 and chi == 1)

    return (r1, r2, r3, r4)


def is_valid(bits: State) -> bool:
    """Return True if all four constraints are satisfied."""
    return all(check_constraints(bits))


# ============================================================
# State generation
# ============================================================

def generate_all_states() -> List[State]:
    """Generate all 256 possible 8-bit states."""
    return list(product([0, 1], repeat=8))


def generate_valid_states() -> List[State]:
    """Generate the 45 valid matter fermion states."""
    return [s for s in generate_all_states() if is_valid(s)]


def generate_leptons() -> List[State]:
    """Generate the 9 lepton states (LQ=0, fixed points under CNOT)."""
    return [s for s in generate_valid_states() if s[LQ] == 0]


def generate_quarks() -> List[State]:
    """Generate the 36 quark states (LQ=1, oscillators under CNOT)."""
    return [s for s in generate_valid_states() if s[LQ] == 1]


def generate_sterile_neutrinos() -> List[State]:
    """
    Generate the 3 sterile neutrino candidates.
    These are states that fail ONLY R4 (all other constraints pass).
    """
    result = []
    for s in generate_all_states():
        r1, r2, r3, r4 = check_constraints(s)
        if r1 and r2 and r3 and not r4:
            result.append(s)
    return result


# ============================================================
# CNOT update rule
# ============================================================

def apply_cnot(bits: State) -> State:
    """
    Apply the unique update rule: I3(t+1) = I3(t) ⊕ LQ(t).

    This is a CNOT gate with LQ as control and I3 as target.
    It is the weak interaction, derived from the information action
    principle as the minimum-cost non-trivial dynamics preserving
    the 45-state spectrum.
    """
    lst = list(bits)
    lst[I3] = lst[I3] ^ lst[LQ]
    return tuple(lst)


def get_cnot_pairs() -> List[Tuple[State, State]]:
    """
    Return the 18 quark doublet pairs linked by the CNOT rule.
    Each pair consists of two quark states related by I3 ↔ (1-I3).
    """
    quarks = generate_quarks()
    pairs = []
    seen = set()
    for s in quarks:
        partner = apply_cnot(s)
        pair = (min(s, partner), max(s, partner))
        if pair not in seen:
            seen.add(pair)
            pairs.append(pair)
    return pairs


# ============================================================
# Physical properties
# ============================================================

GENERATION_MAP = {(0, 0): 1, (0, 1): 2, (1, 0): 3}
COLOUR_MAP = {(0, 1): "green", (1, 0): "red", (1, 1): "blue"}


def get_generation(bits: State) -> int:
    """Return the generation number (1, 2, or 3)."""
    return GENERATION_MAP[(bits[G0], bits[G1])]


def get_colour(bits: State) -> str:
    """Return the colour label for quarks, or 'colourless' for leptons."""
    if bits[LQ] == 0:
        return "colourless"
    return COLOUR_MAP.get((bits[C0], bits[C1]), "unknown")


def get_chirality(bits: State) -> str:
    """Return 'L' (left-handed) or 'R' (right-handed)."""
    return "R" if bits[CHI] == 1 else "L"


def get_electric_charge(bits: State) -> float:
    """
    Compute electric charge from the ring bits.
    Q = f(LQ, I3) following Standard Model assignments.
    """
    lq, i3, chi = bits[LQ], bits[I3], bits[CHI]
    if lq == 0:  # Lepton
        if i3 == 0:
            return 0.0   # neutrino (L-handed; R-handed excluded by R4)
        else:
            return -1.0   # charged lepton
    else:  # Quark
        if chi == 0:  # Left-handed
            if i3 == 0:
                return 2.0 / 3.0   # up-type
            else:
                return -1.0 / 3.0  # down-type
        else:  # Right-handed
            if i3 == 0:
                return 2.0 / 3.0   # up-type
            else:
                return -1.0 / 3.0  # down-type


def state_to_str(bits: State) -> str:
    """Format a state as a binary string."""
    return "".join(map(str, bits))


# ============================================================
# Summary
# ============================================================

def print_spectrum():
    """Print the full 45-state fermion spectrum."""
    valid = generate_valid_states()
    leptons = [s for s in valid if s[LQ] == 0]
    quarks = [s for s in valid if s[LQ] == 1]

    print(f"Total valid states: {len(valid)} / 256")
    print(f"  Leptons (LQ=0, fixed points): {len(leptons)}")
    print(f"  Quarks  (LQ=1, oscillators):  {len(quarks)}")
    print(f"  Sterile neutrino candidates:  {len(generate_sterile_neutrinos())}")
    print()

    print("Leptons:")
    for s in leptons:
        gen = get_generation(s)
        chi = get_chirality(s)
        q = get_electric_charge(s)
        print(f"  Gen {gen} {chi:>1s} Q={q:+5.2f}  {state_to_str(s)}")

    print("\nQuark doublets (CNOT pairs):")
    for s1, s2 in get_cnot_pairs():
        gen = get_generation(s1)
        col = get_colour(s1)
        chi = get_chirality(s1)
        q1 = get_electric_charge(s1)
        q2 = get_electric_charge(s2)
        print(f"  Gen {gen} {col:>5s} {chi:>1s}  "
              f"I3={s1[I3]}(Q={q1:+.2f}) ↔ I3={s2[I3]}(Q={q2:+.2f})  "
              f"{state_to_str(s1)} ↔ {state_to_str(s2)}")


if __name__ == "__main__":
    print_spectrum()
