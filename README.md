# It from Bit, Revisited

**Standard Model Fermions as an Error-Correcting Code on a Holographic Surface**

This repository contains the computational tools and analysis supporting the paper *"It from Bit, Revisited: An Information-Geometric Framework for the Standard Model"*.

## Overview

We propose that the Standard Model fermion spectrum is the set of valid codewords of a quantum error-correcting code on a holographic surface. Each fermion is an oriented 8-bit ring — a *circlette* — where four local constraints select exactly 45 valid matter states from 256 possibilities. The ring partitions into generation, colour, and electroweak sectors connected by a bridge bit, mirroring SU(3) × SU(2) × U(1).

The paper's central result is that the **information action principle selects a unique update rule**: a CNOT gate at the bridge-isospin boundary. This is the weak interaction, derived from first principles as the minimum-cost non-trivial dynamics preserving the 45-state spectrum.

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
├── paper/
│   ├── it-from-bit-unified.tex    # LaTeX source
│   ├── references.bib             # Bibliography
│   └── it-from-bit-unified.pdf    # Compiled paper
├── src/
│   ├── circlette.py               # Core: encoding, constraints, state generation
│   ├── rule_discovery.py          # Rule search and uniqueness proof
│   └── analysis.py                # Detailed analysis of the winning rule
├── tests/
│   └── test_circlette.py          # Verification tests
└── requirements.txt
```

## Quick Start

```bash
pip install numpy
python src/rule_discovery.py
```

## Key Results

| Result | Value |
|--------|-------|
| Valid matter states | 45 / 256 |
| Valid incl. antimatter | 90 / 256 |
| Non-trivial rules preserving all 45 states | **1** |
| The unique rule | I₃ ⊕ LQ (CNOT gate) |
| Fixed points (leptons) | 9 |
| Oscillating pairs (quarks) | 18 (36 states) |
| Average bit-flip cost | 36/45 = 0.80 bits/tick |
| Rule order | 2 (involution: M² = I) |

## The Four Constraints

| Rule | Constraint | Physical meaning |
|------|-----------|------------------|
| R1 | G₀G₁ ≠ 11 | Three generations only |
| R2 | χ = W | Chirality gate |
| R3 | C₁=1 ⟹ LQ=1 | Colour requires quark identity |
| R4 | LQ=0 ∧ I₃=0 ⟹ χ=0 | No right-handed neutrinos |

Plus: leptons must be colourless; quarks must carry colour.

## Citation

```bibtex
@article{gooding2026itfrombit,
  title   = {It from Bit, Revisited: An Information-Geometric Framework
             for the Standard Model},
  author  = {Gooding, David},
  year    = {2026},
  note    = {Preprint}
}
```

## Author

David Gooding — [Neuro-Symbolic Ltd](https://neuro-symbolic.co.uk)

## License

MIT License — see [LICENSE](LICENSE) for details.
