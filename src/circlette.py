"""
circlette.py — Core encoding and constraint system for the circlette framework.

Defines the 8-bit ring structure, the four local constraints (R1-R4),
and generates the 45 valid matter fermion states of the Standard Model.

Ring layout:
    Index:   0   1   2   3   4   5   6   7
    Bit:    G0  G1  C0  C1  LQ  I3  chi  W
    Sector: [generation] [colour] [bridge] [electroweak]
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Dict, List, Set, Tuple

# Type alias for an 8-bit state
State = Tuple[int, ...]

# Bit positions on the ring
G0, G1, C0, C1, LQ, I3, CHI, W = range(8)

LABELS = ["G0", "G1", "C0", "C1", "LQ", "I3", "χ", "W"]

# Sector boundaries (source_bit, target_bit, name)
SECTOR_BOUNDARIES = [
    (G1, C0, "G1→C0"),
    (C0, G1, "C0→G1"),
    (C1, LQ, "C1→LQ"),
    (LQ, C1, "LQ→C1"),
    (LQ, I3, "LQ→I3"),
    (I3, LQ, "I3→LQ"),
    (W, G0, "W→G0"),
    (G0, W, "G0→W"),
]


@dataclass(frozen=True)
class Particle:
    """A named particle state on the circlette ring."""

    name: str
    state: State
    particle_type: str  # "lepton" or "quark"
    generation: int  # 1, 2, or 3
    chirality: str  # "L" or "R"
    charge: float  # electric charge
    mass_mev: float  # mass in MeV


def is_valid_state(bits: State) -> bool:
    """Apply the four local constraints to test whether a state is valid.

    R1: G0,G1 ≠ (1,1)               — three generations only
    R2: chi = W                      — chirality gate
    R3: C1=1 requires LQ=1           — colour requires quark identity
    R4: LQ=0, I3=0 forbids chi=1    — no right-handed neutrinos

    Plus: leptons (LQ=0) must be colourless
          quarks (LQ=1) must carry colour
    """
    g0, g1, c0, c1, lq, i3, chi, w = bits

    if g0 == 1 and g1 == 1:
        return False  # R1
    if chi != w:
        return False  # R2
    if c1 == 1 and lq == 0:
        return False  # R3
    if lq == 0 and i3 == 0 and chi == 1:
        return False  # R4
    if lq == 0 and (c0 != 0 or c1 != 0):
        return False  # leptons colourless
    if lq == 1 and c0 == 0 and c1 == 0:
        return False  # quarks carry colour

    return True


def generate_valid_states() -> List[State]:
    """Generate all 45 valid matter fermion states."""
    return [bits for bits in product([0, 1], repeat=8) if is_valid_state(bits)]


def generate_all_states() -> List[State]:
    """Generate all 256 possible 8-bit states."""
    return list(product([0, 1], repeat=8))


def state_to_str(state: State) -> str:
    """Format a state as a binary string."""
    return "".join(map(str, state))


def _build_particle_catalogue() -> Dict[State, Particle]:
    """Build the complete particle name catalogue."""
    gen_label = {(0, 0): 1, (0, 1): 2, (1, 0): 3}
    colour_label = {(1, 0): "r", (0, 1): "g", (1, 1): "b"}

    # Lepton names by generation
    lepton_names = {
        1: {"nu": "νe", "charged": "e⁻"},
        2: {"nu": "νμ", "charged": "μ⁻"},
        3: {"nu": "ντ", "charged": "τ⁻"},
    }

    # Quark names by generation and isospin
    quark_names = {
        1: {0: "u", 1: "d"},
        2: {0: "c", 1: "s"},
        3: {0: "t", 1: "b"},
    }

    # Lepton masses (MeV)
    lepton_masses = {
        "νe": 0.0, "e⁻": 0.511,
        "νμ": 0.0, "μ⁻": 105.66,
        "ντ": 0.0, "τ⁻": 1776.9,
    }

    # Quark masses (MeV)
    quark_masses = {
        "u": 2.2, "d": 4.7, "c": 1270.0, "s": 95.0,
        "t": 173100.0, "b": 4180.0,
    }

    # Electric charges
    lepton_charges = {"νe": 0, "νμ": 0, "ντ": 0, "e⁻": -1, "μ⁻": -1, "τ⁻": -1}
    quark_charges = {"u": 2 / 3, "c": 2 / 3, "t": 2 / 3, "d": -1 / 3, "s": -1 / 3, "b": -1 / 3}

    catalogue: Dict[State, Particle] = {}

    for state in generate_valid_states():
        g0, g1, c0, c1, lq, i3, chi, w = state
        gen = gen_label[(g0, g1)]
        chir = "L" if chi == 0 else "R"

        if lq == 0:
            # Lepton
            if i3 == 0:
                base = lepton_names[gen]["nu"]
            else:
                base = lepton_names[gen]["charged"]
            name = f"{base}_{chir}"
            charge = lepton_charges[base]
            mass = lepton_masses[base]
            catalogue[state] = Particle(name, state, "lepton", gen, chir, charge, mass)
        else:
            # Quark
            col = colour_label.get((c0, c1), "?")
            base = quark_names[gen][i3]
            name = f"{base}_{col}_{chir}"
            charge = quark_charges[base]
            mass = quark_masses[base]
            catalogue[state] = Particle(name, state, "quark", gen, chir, charge, mass)

    return catalogue


# Module-level catalogue
PARTICLE_CATALOGUE = _build_particle_catalogue()


def particle_name(state: State) -> str:
    """Look up the particle name for a valid state."""
    p = PARTICLE_CATALOGUE.get(state)
    return p.name if p else f"[{state_to_str(state)}]"
