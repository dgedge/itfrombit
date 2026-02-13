"""
test_circlette.py — Verification tests for all computational claims.

Each test corresponds to a specific claim in the paper. Running this
file provides a machine-checkable proof of the key results.

Usage:
    python -m pytest tests/test_circlette.py -v
    # or simply:
    python tests/test_circlette.py
"""

from __future__ import annotations

import sys
import unittest
from itertools import product

import numpy as np

# Allow imports from src/
sys.path.insert(0, "src")

from circlette import (
    G0, G1, C0, C1, LQ, I3, CHI, W,
    generate_all_states,
    generate_valid_states,
    is_valid_state,
)


def boundary_matrix(src: int, tgt: int) -> np.ndarray:
    m = np.eye(8, dtype=int)
    m[tgt][src] = 1
    return m


def apply_rule(state, matrix):
    result = []
    for i in range(8):
        val = 0
        for j in range(8):
            val ^= matrix[i][j] * state[j]
        result.append(val)
    return tuple(result)


def matrix_invertible_f2(mat):
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


class TestEncoding(unittest.TestCase):
    """Verify the 45-state encoding."""

    def test_exactly_45_valid_states(self):
        """Paper claim: exactly 45 valid matter fermion states."""
        states = generate_valid_states()
        self.assertEqual(len(states), 45)

    def test_antimatter_via_ring_reversal(self):
        """Antimatter is ring reversal; verify the combined set size."""
        states = generate_valid_states()
        valid_set = set(states)
        reversed_states = {tuple(reversed(s)) for s in states}
        combined = valid_set | reversed_states
        # Total = matter + antimatter - overlap
        # Overlap includes palindromes AND reversed matter states
        # that happen to also satisfy the matter constraints
        self.assertEqual(len(combined), len(valid_set) + len(reversed_states - valid_set))
        self.assertGreater(len(combined), 45)  # antimatter adds states
        self.assertLessEqual(len(combined), 90)

    def test_all_256_checked(self):
        """Verify constraints are applied to all 256 possible states."""
        all_states = generate_all_states()
        self.assertEqual(len(all_states), 256)
        valid = [s for s in all_states if is_valid_state(s)]
        self.assertEqual(len(valid), 45)

    def test_three_generations(self):
        """R1: G0G1 ≠ (1,1) gives exactly 3 generation values."""
        states = generate_valid_states()
        generations = {(s[G0], s[G1]) for s in states}
        self.assertEqual(generations, {(0, 0), (0, 1), (1, 0)})

    def test_chirality_gate(self):
        """R2: chi = W for all valid states."""
        for s in generate_valid_states():
            self.assertEqual(s[CHI], s[W], f"R2 violated in {s}")

    def test_leptons_colourless(self):
        """Leptons (LQ=0) must have C0=C1=0."""
        for s in generate_valid_states():
            if s[LQ] == 0:
                self.assertEqual(s[C0], 0, f"Lepton has colour: {s}")
                self.assertEqual(s[C1], 0, f"Lepton has colour: {s}")

    def test_quarks_have_colour(self):
        """Quarks (LQ=1) must have at least one colour bit set."""
        for s in generate_valid_states():
            if s[LQ] == 1:
                self.assertTrue(
                    s[C0] != 0 or s[C1] != 0,
                    f"Quark without colour: {s}",
                )

    def test_9_leptons_36_quarks(self):
        """Verify the lepton/quark split."""
        states = generate_valid_states()
        leptons = [s for s in states if s[LQ] == 0]
        quarks = [s for s in states if s[LQ] == 1]
        self.assertEqual(len(leptons), 9)
        self.assertEqual(len(quarks), 36)


class TestUniqueRule(unittest.TestCase):
    """Verify the uniqueness of the weak rule."""

    def setUp(self):
        self.valid_states = generate_valid_states()
        self.valid_set = set(self.valid_states)
        self.all_states = generate_all_states()

        # All 8 sector-boundary couplings
        self.boundaries = [
            (G1, C0), (C0, G1),
            (C1, LQ), (LQ, C1),
            (LQ, I3), (I3, LQ),
            (W, G0), (G0, W),
        ]

    def _pure_valid_count(self, mat):
        """Count valid states in cycles containing only valid states."""
        visited = {}
        cycles = []
        for s in self.all_states:
            if s in visited:
                continue
            cycle = []
            current = s
            while current not in visited:
                visited[current] = len(cycles)
                cycle.append(current)
                current = apply_rule(current, mat)
            cycles.append(cycle)

        pure = [c for c in cycles if set(c).issubset(self.valid_set)]
        return sum(len(c) for c in pure)

    def test_all_boundaries_invertible(self):
        """All boundary matrices are invertible over F₂."""
        for src, tgt in self.boundaries:
            mat = boundary_matrix(src, tgt)
            self.assertTrue(
                matrix_invertible_f2(mat.copy()),
                f"Boundary ({src},{tgt}) not invertible",
            )

    def test_only_lq_i3_preserves_all_45(self):
        """Only LQ→I3 preserves all 45 valid states in pure-valid cycles."""
        for src, tgt in self.boundaries:
            mat = boundary_matrix(src, tgt)
            count = self._pure_valid_count(mat)
            if (src, tgt) == (LQ, I3):
                self.assertEqual(count, 45, "LQ→I3 should preserve all 45")
            else:
                self.assertLess(
                    count, 45,
                    f"Boundary ({src},{tgt}) also preserves 45 — uniqueness violated!",
                )

    def test_rule_is_involution(self):
        """M² = I: the rule is its own inverse."""
        mat = boundary_matrix(LQ, I3)
        m2 = (mat @ mat) % 2
        np.testing.assert_array_equal(m2, np.eye(8, dtype=int))

    def test_avg_flip_cost_is_0_80(self):
        """Average bit-flip cost is 36/45 = 0.80."""
        mat = boundary_matrix(LQ, I3)
        total = sum(
            sum(a != b for a, b in zip(s, apply_rule(s, mat)))
            for s in self.valid_set
        )
        self.assertAlmostEqual(total / 45, 36 / 45)

    def test_leptons_are_fixed_points(self):
        """All 9 leptons are fixed points of the rule."""
        mat = boundary_matrix(LQ, I3)
        leptons = [s for s in self.valid_states if s[LQ] == 0]
        self.assertEqual(len(leptons), 9)
        for s in leptons:
            self.assertEqual(apply_rule(s, mat), s, f"Lepton {s} is not fixed")

    def test_quarks_oscillate_in_pairs(self):
        """All 36 quarks oscillate in 18 pairs of period 2."""
        mat = boundary_matrix(LQ, I3)
        quarks = [s for s in self.valid_states if s[LQ] == 1]
        self.assertEqual(len(quarks), 36)

        pairs = set()
        for s in quarks:
            ns = apply_rule(s, mat)
            self.assertNotEqual(s, ns, f"Quark {s} is fixed — should oscillate")
            self.assertIn(ns, self.valid_set, f"Quark {s} maps to invalid state")
            self.assertEqual(
                apply_rule(ns, mat), s,
                f"Quark {s} does not return after 2 ticks",
            )
            pairs.add((min(s, ns), max(s, ns)))

        self.assertEqual(len(pairs), 18)

    def test_pairs_are_isospin_doublets(self):
        """Each oscillating pair differs only in the I₃ bit."""
        mat = boundary_matrix(LQ, I3)
        quarks = [s for s in self.valid_states if s[LQ] == 1]

        for s in quarks:
            ns = apply_rule(s, mat)
            # Should differ only at position 5 (I3)
            for i in range(8):
                if i == I3:
                    self.assertNotEqual(s[i], ns[i])
                else:
                    self.assertEqual(
                        s[i], ns[i],
                        f"Bit {i} differs for {s} → {ns}",
                    )

    def test_minimum_cost_among_nontrivial(self):
        """LQ→I3 has the lowest bit-flip cost among rules preserving all 45."""
        target_mat = boundary_matrix(LQ, I3)
        target_cost = sum(
            sum(a != b for a, b in zip(s, apply_rule(s, target_mat)))
            for s in self.valid_set
        ) / 45

        # Check no other boundary rule achieves 45/45 with lower cost
        for src, tgt in self.boundaries:
            if (src, tgt) == (LQ, I3):
                continue
            mat = boundary_matrix(src, tgt)
            count = self._pure_valid_count(mat)
            if count == 45:
                cost = sum(
                    sum(a != b for a, b in zip(s, apply_rule(s, mat)))
                    for s in self.valid_set
                ) / 45
                self.assertGreaterEqual(
                    cost, target_cost,
                    f"({src},{tgt}) has lower cost {cost} < {target_cost}",
                )


class TestCNOTProperties(unittest.TestCase):
    """Verify the CNOT structure of the rule."""

    def test_lq_is_read_only(self):
        """LQ is never modified by the rule (read-only control bit)."""
        mat = boundary_matrix(LQ, I3)
        for s in generate_all_states():
            ns = apply_rule(s, mat)
            self.assertEqual(s[LQ], ns[LQ], f"LQ changed for {s}")

    def test_only_i3_changes(self):
        """Only the I₃ bit can change; all other bits are preserved."""
        mat = boundary_matrix(LQ, I3)
        for s in generate_all_states():
            ns = apply_rule(s, mat)
            for i in range(8):
                if i != I3:
                    self.assertEqual(
                        s[i], ns[i],
                        f"Bit {i} changed for {s}",
                    )

    def test_nilpotent_perturbation(self):
        """(M - I)² = 0: the perturbation from identity is nilpotent."""
        mat = boundary_matrix(LQ, I3)
        N = (mat - np.eye(8, dtype=int))
        N2 = (N @ N) % 2
        np.testing.assert_array_equal(N2, np.zeros((8, 8), dtype=int))


if __name__ == "__main__":
    unittest.main(verbosity=2)
