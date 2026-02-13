"""
analysis.py — Deep analysis of the unique weak rule.

Explores the physical consequences of I₃(t+1) = I₃(t) ⊕ LQ(t):
  - CNOT gate structure and universal entanglement
  - Time generation from the rule
  - Computational Lagrangian factorisation
  - Quark/lepton mass ratios
  - Charge as computational permission

Usage:
    python src/analysis.py
"""

from __future__ import annotations

import sys
from typing import Dict, List, Set, Tuple

import numpy as np

from circlette import (
    LABELS,
    PARTICLE_CATALOGUE,
    State,
    generate_all_states,
    generate_valid_states,
    particle_name,
    state_to_str,
)


def apply_rule(state: State, matrix: np.ndarray) -> State:
    """Apply a linear rule (matrix over F₂) to an 8-bit state."""
    result = []
    for i in range(8):
        val = 0
        for j in range(8):
            val ^= matrix[i][j] * state[j]
        result.append(val)
    return tuple(result)


def main() -> None:
    valid_states = generate_valid_states()
    valid_set = set(valid_states)

    # The unique rule: I₃(t+1) = I₃(t) ⊕ LQ(t)
    mat = np.eye(8, dtype=int)
    mat[5][4] = 1  # I3 ← I3 ⊕ LQ

    # ==============================================================
    print("=" * 70)
    print("THE CNOT STRUCTURE")
    print("=" * 70)
    print("""
The rule M = I + e₅·e₄ᵀ is a Controlled-NOT gate over F₂:

    LQ (control) ──●── LQ  (unchanged, read-only)
                   │
    I₃ (target)  ──⊕── I₃ ⊕ LQ

The CNOT is the universal entangling gate of quantum computation.
Its appearance as the fundamental dynamics of the circlette ring
means the weak interaction is the universe's entangling operation.

Ring decomposition under the rule:
  Control register:  LQ          (1 bit, read-only)
  Target register:   I₃          (1 bit, conditionally flipped)
  Environment:       G₀G₁C₀C₁χW (6 bits, spectators)
""")

    # ==============================================================
    print("=" * 70)
    print("THE 18 WEAK ISOSPIN DOUBLETS")
    print("=" * 70)

    pairs: List[Tuple[State, State]] = []
    fixed: List[State] = []
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

    print(f"\nFixed points (leptons):      {len(fixed)}")
    print(f"Oscillating pairs (quarks):  {len(pairs)}")

    for s1, s2 in sorted(pairs):
        n1 = particle_name(s1)
        n2 = particle_name(s2)
        print(f"  {n1:8s} ↔ {n2:8s}   {state_to_str(s1)} ↔ {state_to_str(s2)}")

    # ==============================================================
    print(f"\n{'=' * 70}")
    print("TIME AS A COMPUTATIONAL PRODUCT")
    print("=" * 70)
    print("""
Before the rule acts, a circlette is a static codeword — no
distinction between tick t and tick t+1.

The rule CREATES the temporal distinction:
  Leptons (LQ=0): rule is identity → no internal change → no clock
  Quarks  (LQ=1): I₃ flips → binary alternation IS the clock

Consequences:
  • Neutrinos (no internal clock) → propagate at c
  • Quarks (internal clock, period 2) → propagate slower than c
  • At high velocity, lattice updates consumed by spatial translation,
    fewer available for internal oscillation → clock slows

The Lorentz factor γ = total_updates / internal_updates_available

When v → c: all updates → propagation, clock stops, γ → ∞

Time dilation is a computational bandwidth constraint.
""")

    # ==============================================================
    print("=" * 70)
    print("THE COMPUTATIONAL LAGRANGIAN")
    print("=" * 70)
    print("""
The total information action factorises:

    S_I = S_external + S_internal

  S_external: Fisher arc-length of lattice propagation
              → geodesic equation (gravity)

  S_internal: bit-flip cost of the rule
              = 36/45 = 0.80 bits per tick
              → isospin oscillation (weak force)

  Cross-terms: coupling between propagation and internal dynamics
              → W boson mass

This is the "Computational Lagrangian" of the weak sector.
The 0.80 bits/tick is the computational vacuum energy —
the minimum power consumption of a 45-state code with
non-trivial dynamics.
""")

    # ==============================================================
    print("=" * 70)
    print("CHARGE AS COMPUTATIONAL PERMISSION")
    print("=" * 70)
    print("""
The bridge bit LQ is never modified by the rule.
It is a READ-ONLY CLASSICAL CONTROL LINE.

  Quarks  (LQ=1): permission to oscillate → "computationally active"
  Leptons (LQ=0): no permission           → "computationally quiet"

Charge quantisation follows: you cannot have half a control line.
The value is 0 or 1 — no continuous interpolation.

This extends to all charges:
  Colour charge (C₀,C₁):  determines strong interaction participation
  Electric charge:         derived from I₃ + LQ (Gell-Mann–Nishijima)
  All charges are PERMISSIONS — read-only bits governing which
  logic gates can execute on a given pattern.
""")

    # ==============================================================
    print("=" * 70)
    print("QUARK/LEPTON MASS RATIOS")
    print("=" * 70)

    mass_data = {
        1: {"lepton": 0.511, "up": 2.2, "down": 4.7,
            "l_name": "e⁻", "u_name": "u", "d_name": "d"},
        2: {"lepton": 105.66, "up": 1270.0, "down": 95.0,
            "l_name": "μ⁻", "u_name": "c", "d_name": "s"},
        3: {"lepton": 1776.9, "up": 173100.0, "down": 4180.0,
            "l_name": "τ⁻", "u_name": "t", "d_name": "b"},
    }

    print(f"\n{'Gen':>3}  {'Up/Lepton':>12}  {'Down/Lepton':>12}  {'Up/Down':>10}")
    print("-" * 42)
    for gen, d in mass_data.items():
        u_l = d["up"] / d["lepton"]
        d_l = d["down"] / d["lepton"]
        u_d = d["up"] / d["down"]
        print(
            f"  {gen}  {d['u_name']}/{d['l_name']} = {u_l:>7.1f}"
            f"  {d['d_name']}/{d['l_name']} = {d_l:>7.1f}"
            f"  {d['u_name']}/{d['d_name']} = {u_d:>5.1f}"
        )

    print("""
The internal rule gives a BINARY mass split:
  Quarks  → period 2 (oscillation mass contribution)
  Leptons → period 1 (no oscillation mass)

But all quarks oscillate at the SAME internal rate.
Mass DIFFERENCES within quarks come from the 6 environment bits
(generation, colour, chirality) via their effect on 2D lattice
propagation — the external part of the Computational Lagrangian.
""")

    # ==============================================================
    print("=" * 70)
    print("EMERGENT SU(2) GAUGE SYMMETRY")
    print("=" * 70)
    print("""
In quantum computation:
  Single-qubit rotations + CNOT = universal gate set

Our framework provides:
  CNOT  → the unique weak rule (internal ring dynamics)
  U(θ)  → lattice-mediated rotations (external Fisher curvature)

Together: a complete universal gate set.

SU(2) gauge symmetry is the CONTINUOUS LIMIT of this discrete
CNOT-based logic circuit operating over macroscopic timescales.

The gauge coupling g₂ is the probability per tick that a
lattice-mediated rotation reaches the target bit.
Its running with energy reflects the lattice's finite bandwidth.
""")

    # ==============================================================
    print("=" * 70)
    print("MATRIX PROPERTIES")
    print("=" * 70)

    print(f"\nM² = I?  {np.array_equal((mat @ mat) % 2, np.eye(8, dtype=int))}")
    print(f"det(M) mod 2 = {int(round(np.linalg.det(mat))) % 2}")

    eigenvalues = np.linalg.eigvals(mat.astype(float))
    print(f"Eigenvalues (real): {sorted(eigenvalues.real)}")

    N = (mat - np.eye(8, dtype=int))
    N2 = (N @ N) % 2
    print(f"(M-I)² mod 2 = 0?  {np.max(N2) == 0}  (nilpotent perturbation)")


if __name__ == "__main__":
    if "src" not in sys.path:
        sys.path.insert(0, "src")
    main()
