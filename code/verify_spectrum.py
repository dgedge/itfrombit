#!/usr/bin/env python3
"""
verify_spectrum.py — Verify all numerical claims in the paper.

Checks:
  - 45 valid states (9 leptons + 36 quarks)
  - 3 sterile neutrino candidates (R4-only failures)
  - 18 CNOT quark doublet pairs
  - Bit-flip cost = 36/45 = 0.800
  - All leptons are fixed points under CNOT
  - CNOT is an involution (M² = I)
  - Charge assignments match Standard Model

Reference: Elliman (2026), "It from Bit, Revisited"
"""

from circlette import (
    generate_valid_states, generate_leptons, generate_quarks,
    generate_sterile_neutrinos, get_cnot_pairs, apply_cnot,
    get_electric_charge, get_generation, get_colour, get_chirality,
    state_to_str, LQ, I3
)


def verify_state_counts():
    """Verify the 45-state spectrum decomposition."""
    valid = generate_valid_states()
    leptons = generate_leptons()
    quarks = generate_quarks()
    sterile = generate_sterile_neutrinos()
    pairs = get_cnot_pairs()

    print("=" * 60)
    print("STATE COUNT VERIFICATION")
    print("=" * 60)

    checks = [
        ("Total valid states", len(valid), 45),
        ("Leptons (LQ=0)", len(leptons), 9),
        ("Quarks (LQ=1)", len(quarks), 36),
        ("Sterile neutrinos", len(sterile), 3),
        ("CNOT doublet pairs", len(pairs), 18),
    ]

    all_pass = True
    for name, actual, expected in checks:
        ok = actual == expected
        all_pass = all_pass and ok
        status = "✓" if ok else "✗"
        print(f"  {status} {name}: {actual} (expected {expected})")

    return all_pass


def verify_cnot_properties():
    """Verify CNOT rule properties."""
    valid = generate_valid_states()
    leptons = generate_leptons()
    quarks = generate_quarks()

    print("\n" + "=" * 60)
    print("CNOT RULE VERIFICATION")
    print("=" * 60)

    # Leptons are fixed points
    lepton_fixed = all(apply_cnot(s) == s for s in leptons)
    print(f"  {'✓' if lepton_fixed else '✗'} All leptons are fixed points")

    # Quarks oscillate (I3 flips)
    quark_flips = all(apply_cnot(s) != s for s in quarks)
    print(f"  {'✓' if quark_flips else '✗'} All quarks flip I3")

    # CNOT is an involution (M² = I)
    involution = all(apply_cnot(apply_cnot(s)) == s for s in valid)
    print(f"  {'✓' if involution else '✗'} CNOT is an involution (M²=I)")

    # CNOT preserves the 45-state spectrum
    from circlette import is_valid
    preserves = all(is_valid(apply_cnot(s)) for s in valid)
    print(f"  {'✓' if preserves else '✗'} CNOT preserves all 45 states")

    # Bit-flip cost = 36/45
    flips = sum(1 for s in valid if apply_cnot(s) != s)
    cost = flips / len(valid)
    cost_ok = abs(cost - 36 / 45) < 1e-10
    print(f"  {'✓' if cost_ok else '✗'} Bit-flip cost: {flips}/{len(valid)} = {cost:.4f}")

    return lepton_fixed and quark_flips and involution and preserves and cost_ok


def verify_charges():
    """Verify electric charge assignments match the Standard Model."""
    valid = generate_valid_states()

    print("\n" + "=" * 60)
    print("CHARGE VERIFICATION")
    print("=" * 60)

    # Expected charges
    expected_lepton_charges = {0.0, -1.0}  # neutrinos and charged leptons
    expected_quark_charges = {2.0 / 3.0, -1.0 / 3.0}  # up-type and down-type

    lepton_charges = set()
    quark_charges = set()
    for s in valid:
        q = get_electric_charge(s)
        if s[LQ] == 0:
            lepton_charges.add(q)
        else:
            quark_charges.add(q)

    lep_ok = lepton_charges == expected_lepton_charges
    quark_ok = quark_charges == expected_quark_charges
    print(f"  {'✓' if lep_ok else '✗'} Lepton charges: {sorted(lepton_charges)}")
    print(f"  {'✓' if quark_ok else '✗'} Quark charges: {sorted(quark_charges)}")

    return lep_ok and quark_ok


def verify_generation_structure():
    """Verify generation assignments and CNOT pair structure."""
    pairs = get_cnot_pairs()

    print("\n" + "=" * 60)
    print("GENERATION AND DOUBLET STRUCTURE")
    print("=" * 60)

    gen_counts = {1: 0, 2: 0, 3: 0}
    for s1, s2 in pairs:
        gen = get_generation(s1)
        gen_counts[gen] += 1

    # Each generation should have 6 pairs (3 colours × 2 chiralities)
    gen_ok = all(v == 6 for v in gen_counts.values())
    for g, c in sorted(gen_counts.items()):
        print(f"  {'✓' if c == 6 else '✗'} Generation {g}: {c} doublet pairs (expected 6)")

    return gen_ok


if __name__ == "__main__":
    results = [
        verify_state_counts(),
        verify_cnot_properties(),
        verify_charges(),
        verify_generation_structure(),
    ]

    print("\n" + "=" * 60)
    if all(results):
        print("ALL VERIFICATIONS PASSED ✓")
    else:
        print("SOME VERIFICATIONS FAILED ✗")
    print("=" * 60)
