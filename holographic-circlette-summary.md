# The Holographic Circlette

**Standard Model Fermions, Quantum Mechanics, and the Fine Structure Constant from an 8-Bit Error-Correcting Code**

David Elliman,· Neuro-Symbolic Ltd · February 2026

---

## Premise

The paper proposes that fundamental physics emerges from a single computational structure: an 8-bit binary codeword (the *circlette*) on a 2D holographic lattice at the Planck scale. Four local parity constraints (R1–R4) select exactly **45 valid states** from 256 possible, reproducing the complete Standard Model fermion spectrum: 3 generations of quarks and leptons with correct charges, chiralities, and colour assignments. A unique spectrum-preserving bit operation — the CNOT gate copying the chirality bit to the weak bit — is identified as the **weak interaction**.

## Exact Derivations (not fitted — derived from the code)

- **3+1D Dirac equation.** The CNOT oscillation on the isospin bit I₃, embedded in the chirality–isospin tensor product space, yields a 4-component quantum walk whose continuum limit is the *exact* 3+1D Dirac equation. All 10 anticommutation relations of Cl(3,1) are computationally verified.

- **Three spatial dimensions from two.** The three Pauli generators of SU(2)\_I₃ produce three independent momentum components on a 2D lattice surface. The third spatial dimension emerges algebraically, not geometrically.

- **3+1D Schrödinger equation.** The standard non-relativistic limit, via the Pauli identity, gives the exact Schrödinger equation for a spin-½ particle.

- **Complex phase.** The imaginary unit *i* is forced by unitarity of the reversible boolean swap.

- **Dynamic dark energy.** The vacuum occupation fraction Φ = 45/256 generates a dynamic equation of state w(z) that crosses −1 at z ≈ 0.41 and fits DESI DR2 data (χ²/dof = 0.41).

## Gauge Theory and the Fine Structure Constant

Extending the framework to include U(1) gauge fields on lattice links (Wilson, 1974), we identify the QED interaction vertex as a **mandatory one-tick code violation**: the Dirac hop flips the chirality bit χ, temporarily violating constraint R2, while the link variable injects a charge-dependent phase. Two results follow automatically from the code's constraint structure:

- **Anomaly cancellation:** Σ Q = 0 (gravitational anomaly, required for quantum consistency)
- **Beta function coefficient:** Σ Q² = 3 × (16/3) = **16** (exact Standard Model value for 1-loop QED running)

We conjecture an **information saturation principle**: the total gauge entropy per lattice tick, summed over all charged species, saturates the Holevo capacity (1 bit) of a single lattice link:

> **S\_BE(α) × Σ Q² = 1 bit**

| Formula | α⁻¹ | Error |
|---|---|---|
| Shannon (tree-level): H(α) × 16 = 1 | 136.48 | 0.41% |
| Bose–Einstein (tree): S\_BE(α) × 16 = 1 | 136.68 | 0.26% |
| Bose + 1-loop: S\_BE(α) × 16 × (1+α/π) = 1 | 137.06 | **0.02%** |
| **Experimental** | **137.036** | — |

## Key Claim

No free parameters are introduced between the 8-bit code and the wave equations of quantum mechanics. The framework does not replace the Standard Model — it *derives* it, along with gravity (as information geometry), quantum mechanics (as the continuum limit of the CNOT walk), and dynamic dark energy (from the vacuum occupation fraction), from a single information-theoretic substrate. **The laws of physics are the error-correction rules of a computational universe.**

## Status & Testable Predictions

The paper contains 20 enumerated predictions including: the dark energy equation of state crossing w = −1 at z ≈ 0.41 (testable by DESI 5-year data), three sterile neutrinos with gravitational interaction only, the neutrino mass scale m_ν ~ 10⁻³ eV from the vacuum floor, massless bare leptons, and the information saturation prediction for α. The 23-page paper with full derivations and computational verifications is available on Github https://github.com/dgedge/itfrombit/blob/main/it-from-bit-final.pdf
