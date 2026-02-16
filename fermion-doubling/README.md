# Fermion Doubling on the 4.8.8 Lattice

**Companion code to: "The Holographic Circlette" (Shorten, 2026)**

## Summary

This repository contains the numerical analysis of fermion doubling on the 4.8.8 truncated square lattice, in the context of the circlette framework where Standard Model fermions emerge from an 8-bit error-correcting code on a 2D holographic lattice.

### Key Results

1. **The naive square lattice has 4 Dirac points** (1 physical + 3 doublers) — the standard Nielsen-Ninomiya fermion doubling problem.

2. **The Wilson mass (added by hand) gaps all 3 doublers** via the momentum-dependent term M(k) = r(2 − cos kₓ − cos kᵧ), proportional to the mass matrix β = σ_z^χ ⊗ I.

3. **The 4.8.8 NNN term does NOT gap the doublers.** The physical two-hop process through interstitial squares generates a term proportional to α₁α₂ = i(I^χ ⊗ σ_z^{I₃}) × sin(kₓ)sin(kᵧ), which vanishes at all doubler locations. This rigorously disproves the geometric Wilson mass hypothesis.

4. **The CNOT coin operator resolves fermion doubling at the dynamical level.** The circlette's discrete Z₂ chirality constraint (R2: χ ⊕ W = 0) breaks the continuous U(1)_A chiral symmetry that the Nielsen-Ninomiya theorem requires. The doublers are not gapped — they never exist, because the topological winding argument has no jurisdiction over a discrete binary state space.

### The "Hardware vs Software" Separation

- **Hardware (4.8.8 topology):** Provides n=8 bandwidth matching, clean separation of 4 kinematic channels from 4 gauge channels, and the natural plaquette structure for Wilson loops.
- **Software (CNOT dynamics):** Resolves fermion doubling through discrete chiral mixing at each clock tick — the algorithmic equivalent of a Wilson mass, built into the fundamental update rule rather than added by hand.

## Files

| File | Description |
|------|-------------|
| `dirac_matrices.py` | Circlette Dirac algebra: α₁, α₂, β in the χ ⊗ I₃ space, with anticommutation verification |
| `naive_square.py` | Fermion doubling on the standard square lattice (baseline) |
| `wilson_fermions.py` | Wilson mass resolution (standard fix, added by hand) |
| `lattice_488.py` | Physical Hamiltonian on the 4.8.8 lattice with derived NNN coupling |
| `band_structure.py` | Comparative band structure plots along Γ→X→M→Γ |
| `dispersion_scan.py` | Full Brillouin zone scans with Dirac point counting |
| `run_all.py` | Run the complete analysis and generate all figures |

## Requirements

```
numpy
scipy
matplotlib
```

## Usage

```bash
# Run the complete analysis
python run_all.py

# Or run individual components
python naive_square.py
python wilson_fermions.py
python lattice_488.py
python band_structure.py
```

## The Dirac Algebra

The circlette's 4-component coin space (χ, I₃) gives the 2D Dirac matrices:

```
α₁ = σ_x^χ ⊗ σ_x^{I₃}    (hop in x-direction)
α₂ = σ_x^χ ⊗ σ_y^{I₃}    (hop in y-direction)
β  = σ_z^χ ⊗ I^{I₃}       (mass / rest energy)
```

These satisfy the full Clifford algebra Cl(2,1): {αᵢ, αⱼ} = 2δᵢⱼ, {αᵢ, β} = 0, β² = I.

The key finding is that the two-hop NNN product α₂·α₁ = i(I^χ ⊗ σ_z^{I₃}) acts on the **isospin** bit, not the **chirality** bit. Since Wilson mass requires β = σ_z^χ (chirality), the geometric mechanism fails.

## Citation

```bibtex
@article{shorten2026circlette,
  title={The Holographic Circlette: Standard Model Fermions from an 8-Bit Error-Correcting Code},
  author={Shorten, David J.},
  year={2026},
  institution={Neuro-Symbolic Ltd}
}
```

## Licence

MIT
