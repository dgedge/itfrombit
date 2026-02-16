#!/usr/bin/env python3
"""
test_circlette.py — Unit tests for the circlette framework.

Run with: python -m pytest tests/test_circlette.py -v
Or:       python tests/test_circlette.py
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from circlette import (
    generate_valid_states, generate_leptons, generate_quarks,
    generate_sterile_neutrinos, get_cnot_pairs, apply_cnot,
    is_valid, check_constraints, get_electric_charge,
    get_generation, get_colour, get_chirality,
    G0, G1, C0, C1, LQ, I3, CHI, W
)


def test_state_counts():
    assert len(generate_valid_states()) == 45
    assert len(generate_leptons()) == 9
    assert len(generate_quarks()) == 36
    assert len(generate_sterile_neutrinos()) == 3
    assert len(get_cnot_pairs()) == 18


def test_cnot_is_involution():
    for s in generate_valid_states():
        assert apply_cnot(apply_cnot(s)) == s, f"CNOT not involution for {s}"


def test_cnot_preserves_spectrum():
    for s in generate_valid_states():
        assert is_valid(apply_cnot(s)), f"CNOT maps {s} to invalid state"


def test_leptons_are_fixed_points():
    for s in generate_leptons():
        assert apply_cnot(s) == s, f"Lepton {s} not fixed under CNOT"


def test_quarks_oscillate():
    for s in generate_quarks():
        assert apply_cnot(s) != s, f"Quark {s} is fixed (should oscillate)"


def test_bit_flip_cost():
    valid = generate_valid_states()
    flips = sum(1 for s in valid if apply_cnot(s) != s)
    assert flips == 36
    assert abs(flips / len(valid) - 0.8) < 1e-10


def test_lepton_charges():
    charges = {get_electric_charge(s) for s in generate_leptons()}
    assert charges == {0.0, -1.0}


def test_quark_charges():
    charges = {get_electric_charge(s) for s in generate_quarks()}
    assert charges == {2.0 / 3.0, -1.0 / 3.0}


def test_generations():
    for gen in [1, 2, 3]:
        gen_states = [s for s in generate_valid_states()
                      if get_generation(s) == gen]
        assert len(gen_states) == 15, f"Gen {gen} has {len(gen_states)} states"


def test_colour_structure():
    quarks = generate_quarks()
    colours = {get_colour(s) for s in quarks}
    assert colours == {"red", "green", "blue"}

    for s in generate_leptons():
        assert get_colour(s) == "colourless"


def test_sterile_neutrinos_fail_only_r4():
    for s in generate_sterile_neutrinos():
        r1, r2, r3, r4 = check_constraints(s)
        assert r1 and r2 and r3 and not r4


def test_doublet_pairs_per_generation():
    pairs = get_cnot_pairs()
    for gen in [1, 2, 3]:
        gen_pairs = [p for p in pairs if get_generation(p[0]) == gen]
        assert len(gen_pairs) == 6  # 3 colours × 2 chiralities


def test_cnot_rule_is_unique():
    """Verify no other sector-boundary XOR rule preserves the spectrum."""
    valid = generate_valid_states()
    BOUNDARIES = [(G1, C0), (C0, G1), (C1, LQ), (LQ, C1),
                  (LQ, I3), (I3, LQ), (W, G0), (G0, W)]

    nontrivial_preserving = []
    for source, target in BOUNDARIES:
        preserves = True
        n_fixed = 0
        for s in valid:
            new = list(s)
            new[target] = new[target] ^ new[source]
            new = tuple(new)
            if not is_valid(new):
                preserves = False
                break
            if new == s:
                n_fixed += 1

        if preserves and n_fixed < 45:
            nontrivial_preserving.append((source, target))

    assert len(nontrivial_preserving) == 1
    assert nontrivial_preserving[0] == (LQ, I3)


if __name__ == "__main__":
    # Simple test runner without pytest
    import traceback
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = failed = 0
    for test in tests:
        try:
            test()
            print(f"  ✓ {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"  ✗ {test.__name__}: {e}")
            traceback.print_exc()
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
