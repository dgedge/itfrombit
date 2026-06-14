# It from Bit Canon Source Repository

This repository collects the source code used to audit, reproduce, and stress-test
the finite-QEC / It-from-Bit canon.  The current code migration is organized by
the seven-paper structure:

- overview and foundations/methodology,
- matter, gauge structure, and spectroscopy,
- cosmology, dark energy, and inflation,
- dark matter, MOND, and K04 debris,
- gravity, horizons, and black holes,
- relativity, photon identity, and causal-set QED.

For the papers and public write-ups, start here:

[David Elliman on ResearchGate](https://www.researchgate.net/profile/David-Elliman-2/research)

The working canon lives in `ANCHOR.md` and `DRIFT.md` in the main canon tree.  This
repo is the public source-code companion: it should contain runnable scripts,
small fixtures, topic READMEs, and manifests, but not bulky run output.

## Book Reference

The companion book manuscript is not stored in this source-code repository.  Use
the public listing instead:

[A Model of John Wheeler's It from Bit: Eight Easy Pieces](https://www.amazon.co.uk/Model-John-Wheelers-Bit-Eight-Bit/dp/1919558853/ref=tmm_hrd_swatch_0?_encoding=UTF8&dib_tag=se&dib=eyJ2IjoiMSJ9.b6Y4yUDpPou9XoHeADqxNTfg-O7J2RbPwJE2dHHC9OE.UfVGiibvXu48ZTc1jm6-6LpP3zsq1r9yDZ0mexPqnHw&qid=1781426196&sr=1-2)

## Current Code Layout

| Path | Purpose |
|---|---|
| [`scripts/`](scripts/) | Canon-current reproducibility scripts, grouped by topic. |
| [`docs/script_manifest.csv`](docs/script_manifest.csv) | Machine-readable inventory of migrated scripts, source paths, descriptions, and run commands. |
| [`tools/import_canon_scripts.py`](tools/import_canon_scripts.py) | Repeatable importer from the working canon tree into this repo. |
| [`data/`](data/) | Small, named input fixtures only. Large datasets stay external. |
| [`results/`](results/) | Small, named output summaries only. Bulk sweeps and run directories stay external. |

The old `code/` directory has been removed because it was stale.  Current source
code should be imported from `octahedrons/python_code` into `scripts/` using
`tools/import_canon_scripts.py`.

## Running Scripts

Most scripts are self-contained numerical or audit checks.  Topic READMEs give a
direct command for every migrated script.  Example:

```bash
PYTHONPATH=scripts/cosmology_dark_energy_inflation:$PYTHONPATH \
  python scripts/cosmology_dark_energy_inflation/cc_generation_vertex_item115_loop.py
```

Some scripts require heavier optional packages such as SciPy, matplotlib, TeNPy,
CuPy, Qiskit, or GPU-specific libraries.  The baseline environment is ordinary
Python with NumPy/SciPy; heavyweight scripts document their own requirements in
their docstrings or topic README row.

## Curation Standard

Every source script should have:

1. a short purpose statement;
2. a canon or paper reference where applicable;
3. inputs and outputs;
4. a run command from the repository root;
5. an honest status: derived, conditional, exploratory, retired, or imported.

The first migration pass generated descriptions from existing docstrings where
available, and from filenames where needed.  Rows marked "Purpose inferred from
filename" are explicit targets for a later docstring pass.

## What Not To Commit

Do not commit large numerical sweeps, sparse matrices, `.npz`/`.npy` artifacts,
bulk JSONL outputs, local run directories, or generated paper build products.
Small named summary tables may be committed under `results/` when they are
referenced by a README.

## Historical README

# Finding the papers

There are series of 21 papers covering the findings of this work. Most of the standard model can be derived rather than fitted to stat using this approach. Look under Elliman as an author on Zenodo. 

If anyone can endorse me on arXiv I would put the papers there for a wider readership.  I was a full Professor at Nottingham and led a successful research group, but moved a high-technology company  where my work was classified.  I no longer have a university email address.

## The Holographic Circlette

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
