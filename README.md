# The Holographic Circlette

**Unifying the Standard Model, Gravity, and Cosmology via Error-Correcting Codes on a Fisher-Information Lattice**

We propose that the Standard Model fermion spectrum corresponds to the set of valid codewords of an 8-bit quantum error-correcting code on a holographic lattice. Four local constraints select exactly 45 valid matter states from 256 possibilities. A unique update rule — a CNOT gate at the bridge-isospin boundary — emerges from the information action principle as the weak interaction.

From this foundation the framework derives: the Dirac equation as the continuum limit of a CNOT-driven quantum walk; gravity as curvature of the Fisher information metric; special relativity as a bandwidth constraint on the computational substrate; a dynamic cosmological constant matching DESI DR2 observations; and a resolution of the black hole information paradox.

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
├── code/                                  # All computational verification code
│   ├── circlette.py                       # Core: encoding, constraints, state generation
│   ├── rule_discovery.py                  # Rule search and uniqueness proof
│   ├── wave_emergence.py                  # Dirac/Schrödinger from CNOT lattice walk
│   ├── verify_spectrum.py                 # Full verification of paper's numerical claims
│   ├── CKM_matrix_evaluation.py           # CKM matrix derivation from the circlette lattice
│   ├── weinberg_corrected.py              # Weinberg angle / electroweak mixing computation
│   ├── dirac_matrices.py                  # Circlette Dirac algebra with anticommutation checks
│   ├── naive_square.py                    # Fermion doubling on the standard square lattice
│   ├── wilson_fermions.py                 # Wilson mass resolution (standard fix)
│   ├── lattice_488.py                     # Physical Hamiltonian on the 4.8.8 lattice
│   ├── band_structure.py                  # Comparative band structure plots
│   ├── dispersion_scan.py                 # Full Brillouin zone scans
│   ├── run_all.py                         # Run the complete analysis
│   ├── tests/
│   │   └── test_circlette.py              # Unit tests
│   ├── requirements.txt
│   └── LICENSE                            # MIT
│
├── paper/
│   ├── working docs/                      # LaTeX source, figures & research notes
│   │   ├── holographic-circlette-v2.tex   # Main paper LaTeX source
│   │   ├── circlette-lattice-v2.tex       # Lattice paper LaTeX source
│   │   ├── holographic-circlette-references.bib
│   │   ├── circlette-lattice-refs.bib
│   │   ├── circlette-lattice-diagram.png  # Lattice diagram figure
│   │   ├── circlette-lattice-diagram.pdf
│   │   ├── bands_comparison.png           # Band structure comparison figure
│   │   ├── dispersion_surface.png         # Dispersion surface figure
│   │   ├── blackholes-companion.pdf
│   │   ├── holographic-circlette-summary.md
│   │   ├── cheshire_cat_bridge_bit.md
│   │   ├── dynamic_lambda.md
│   │   ├── The Schwinger effect.md
│   │   ├── zero-point-energy-notes_3.md
│   │   ├── David_Bohm_and_the_Circlette_Framework.md
│   │   ├── The Universe is a Giant Magic Screen.md
│   │   ├── Roadmap to Gauge Interactions and the Lamb Shift
│   │   └── Advancements in Quantum Information Systems...md
│   └── working_notes/
│       ├── born_rule_from_lattice.md      # Born rule from lattice configuration counting
│       ├── cheshire_cat_bridge_bit.md     # Quantum Cheshire Cat and the bridge bit
│       └── measurement_and_retrocausality.md  # Lattice resolution of the measurement problem
│
├── it-from-bit-final.pdf                  # Compiled main paper
├── circlette-lattice.pdf                  # Compiled lattice paper
└── LICENSE                                # MIT
```

## Quick Start

```bash
cd code
pip install -r requirements.txt

python circlette.py              # Print the full 45-state spectrum
python rule_discovery.py         # Prove the update rule is unique
python wave_emergence.py         # Derive wave equations from CNOT
python verify_spectrum.py        # Run full verification suite
python -m pytest tests/          # Run unit tests
```

### Electroweak and CKM computations

```bash
python CKM_matrix_evaluation.py  # Derive the CKM matrix
python weinberg_corrected.py     # Compute the Weinberg angle
```

### Fermion doubling analysis (4.8.8 lattice)

```bash
python run_all.py                # Run the complete lattice analysis
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

| Rule | Constraint | Physical meaning |
|------|-----------|------------------|
| R1 | G₀G₁ ≠ 11 | Three generations only |
| R2 | χ = W | Chirality gate |
| R3 | LQ=0 ⟹ C₀=C₁=0 | Colour requires quark identity |
| R4 | LQ=0 ∧ I₃=0 ⟹ χ=0 | No right-handed neutrinos |

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
