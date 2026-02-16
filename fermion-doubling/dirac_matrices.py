"""
Circlette Dirac algebra in the χ ⊗ I₃ coin space.

The 8-bit circlette has two kinematic bits (χ = chirality, I₃ = isospin)
whose tensor product gives a 4-component space supporting the 2+1D
Dirac equation. The continuum limit of the CNOT quantum walk reproduces
the exact 3+1D Dirac equation (see main paper, Part IV).

Matrices:
    α₁ = σ_x^χ ⊗ σ_x^{I₃}   — spatial hop in x-direction
    α₂ = σ_x^χ ⊗ σ_y^{I₃}   — spatial hop in y-direction
    β  = σ_z^χ ⊗ I^{I₃}      — mass matrix (rest energy)

These satisfy the Clifford algebra Cl(2,1):
    {αᵢ, αⱼ} = 2δᵢⱼ I₄
    {αᵢ, β}  = 0
    αᵢ² = β² = I₄
"""

import numpy as np

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

# Circlette Dirac matrices (4×4)
alpha1 = np.kron(sigma_x, sigma_x)   # σ_x^χ ⊗ σ_x^{I₃}
alpha2 = np.kron(sigma_x, sigma_y)   # σ_x^χ ⊗ σ_y^{I₃}
beta   = np.kron(sigma_z, I2)        # σ_z^χ ⊗ I^{I₃}

# Derived quantities
alpha12 = alpha1 @ alpha2  # = i(I^χ ⊗ σ_z^{I₃})
I4 = np.eye(4, dtype=complex)


def verify_clifford():
    """Verify the full Clifford algebra relations."""
    def anticomm(A, B):
        return A @ B + B @ A

    checks = {
        "{α₁, α₂} = 0": np.max(np.abs(anticomm(alpha1, alpha2))),
        "{α₁, β} = 0":   np.max(np.abs(anticomm(alpha1, beta))),
        "{α₂, β} = 0":   np.max(np.abs(anticomm(alpha2, beta))),
        "α₁² = I":        np.max(np.abs(alpha1 @ alpha1 - I4)),
        "α₂² = I":        np.max(np.abs(alpha2 @ alpha2 - I4)),
        "β² = I":          np.max(np.abs(beta @ beta - I4)),
        "α₁α₂ = i(I⊗σz)": np.max(np.abs(alpha12 - 1j * np.kron(I2, sigma_z))),
    }

    print("Clifford Algebra Verification")
    print("=" * 50)
    all_pass = True
    for name, err in checks.items():
        status = "✓" if err < 1e-14 else "✗"
        if err >= 1e-14:
            all_pass = False
        print(f"  {status} {name:<25} (error: {err:.1e})")
    print()
    return all_pass


if __name__ == "__main__":
    verify_clifford()

    print("Matrix representations:")
    print(f"\nα₁ = σ_x^χ ⊗ σ_x^{{I₃}} =\n{np.real(alpha1).astype(int)}")
    print(f"\nα₂ = σ_x^χ ⊗ σ_y^{{I₃}} =\n{alpha2}")
    print(f"\nβ = σ_z^χ ⊗ I^{{I₃}} =\n{np.real(beta).astype(int)}")
    print(f"\nα₁α₂ = i(I^χ ⊗ σ_z^{{I₃}}) =\n{alpha12}")
    print()
    print("Key observation:")
    print("  α₁α₂ acts on the ISOSPIN bit (σ_z^{I₃})")
    print("  β acts on the CHIRALITY bit (σ_z^χ)")
    print("  These are DIFFERENT operators on DIFFERENT bits.")
    print("  This is why the NNN two-hop product cannot generate a Wilson mass.")
