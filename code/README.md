# It from Bit, Revisited

**Computational verification code for the circlette framework.**

> *The lattice does not obey quantum mechanics. Quantum mechanics obeys the lattice.*

## Overview

We propose that the Standard Model fermion spectrum is the set of valid codewords of a quantum error-correcting code on a holographic surface. Each fermion is an oriented 8-bit ring — a *circlette* — where four local constraints select exactly 45 valid matter states from 256 possibilities. The ring partitions into generation, colour, and electroweak sectors connected by a bridge bit, mirroring SU(3) × SU(2) × U(1).

The paper's central results are:

1. The **information action principle selects a unique update rule**: a CNOT gate at the bridge-isospin boundary. This is the weak interaction, derived from first principles.

2. The **Dirac equation emerges** as the continuum limit of a discrete quantum walk whose coin operator is the CNOT gate. Mass is the CNOT execution frequency. The imaginary unit *i* is forced by the unitarity requirement of a reversible Boolean swap.

3. The **Schrödinger equation** follows as the non-relativistic limit. Leptons (LQ=0) bypass the CNOT entirely and propagate as massless Weyl fermions at *c*.

```
Ring layout:  G₀ G₁ │ C₀ C₁ │ LQ │ I₃ χ W
Sectors:      [gen]    [col]  [br]  [electroweak]

The unique rule:  I₃(t+1) = I₃(t) ⊕ LQ(t)

    LQ (control) ──●── LQ  (unchanged)
                   │
    I₃ (target)  ──⊕── I₃ ⊕ LQ
```

## Repository Structure

```
├── README.md
├── LICENSE
├── requirements.txt
├── src/
│   ├── circlette.py           # Core: encoding, constraints, state generation
│   ├── rule_discovery.py      # Rule search and uniqueness proof
│   ├── wave_emergence.py      # Dirac/Schrödinger from CNOT lattice walk
│   └── verify_spectrum.py     # Full verification of paper's claims
├── tests/
│   └── test_circlette.py      # Unit tests (pytest or standalone)
└── working_notes/
    ├── cheshire_cat_bridge_bit.md
    ├── measurement_and_retrocausality.md
    └── born_rule_from_lattice.md
```

## Quick Start

```bash
pip install numpy matplotlib
python src/circlette.py          # Print the full 45-state spectrum
python src/rule_discovery.py     # Prove the update rule is unique
python src/wave_emergence.py     # Derive wave equations from CNOT
python src/verify_spectrum.py    # Run full verification suite
python tests/test_circlette.py   # Run unit tests
```

## Key Results

| Result | Value |
|--------|-------|
| Valid matter states | 45 / 256 |
| Valid incl. antimatter | 90 / 256 |
| Non-trivial spectrum-preserving rules | **1** (unique) |
| The unique rule | I₃ ⊕ LQ (CNOT gate) |
| Fixed points (leptons) | 9 |
| Oscillating pairs (quarks) | 18 (36 states) |
| Average bit-flip cost | 36/45 = 0.80 bits/tick |
| Rule order | 2 (involution: M² = I) |
| Sterile neutrino candidates | 3 (R4-only failures) |
| Dirac equation overlap | 0.986 (Bhattacharyya) |

## The Four Constraints

| Rule | Constraint | Window | Physical meaning |
|------|-----------|--------|------------------|
| R1 | G₀G₁ ≠ 11 | 2 bits | Three generations only |
| R2 | χ = W | 2 bits | Chirality gate |
| R3 | LQ=0 ⟹ C₀=C₁=0 | 3 bits | Colour requires quark identity |
| R4 | LQ=0 ∧ I₃=0 ⟹ χ=0 | 3 bits | No right-handed neutrinos |

## The Wave Equation Derivation Chain

```
CNOT (Boolean NOT on I₃)
  ↓  reversibility requires unitarity
Complex phase (i emerges)
  ↓  couple internal toggle to spatial propagation
Discrete quantum walk
  ↓  continuum limit (Δt, Δx → 0)
1+1D Dirac equation  [mass term mc²σₓ IS the CNOT]
  ↓  non-relativistic limit
Schrödinger equation
```

## Citation

```bibtex
@article{elliman2026itfrombit,
  title   = {It from Bit, Revisited: Standard Model Fermions
             as Codewords of a Holographic Error-Correcting Code},
  author  = {Elliman, D.G.},
  year    = {2026},
  note    = {Preprint}
}
```

## Author

D.G. Elliman — [Neuro-Symbolic Ltd](https://neuro-symbolic.co.uk)

## License

MIT License — see [LICENSE](LICENSE) for details.
