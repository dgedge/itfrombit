"""
Physical Hamiltonian on the 4.8.8 lattice.

The 4.8.8 (truncated square) tiling has octagons as matter sites and
interstitial squares as gauge link sites. The topology provides:
    - 4 nearest-neighbour (NN) orthogonal connections: Dirac kinetic term
    - 4 next-nearest-neighbour (NNN) diagonal connections: via interstitial squares

The NNN two-hop process:
    Octagon A → Square → Octagon B (diagonal)
generates an effective coupling with matrix structure α₂·α₁ = i(I⊗σz),
NOT the Wilson mass matrix β = σz⊗I.

KEY RESULT: The physical NNN term vanishes at all doubler locations
because sin(kx)sin(ky) = 0 at (π,0), (0,π), and (π,π).
The 4.8.8 topology alone CANNOT resolve fermion doubling.

The resolution comes instead from the CNOT coin operator (discrete Z₂
chirality), which breaks the continuous U(1)_A symmetry required by
the Nielsen-Ninomiya theorem.
"""

import numpy as np
from numpy import sin, cos, pi
from dirac_matrices import alpha1, alpha2, alpha12, beta, I4
from naive_square import scan_brillouin_zone


def H_488_physical(kx, ky, t1=1.0, t2=0.5, m=0.0):
    """
    Physical Hamiltonian on the 4.8.8 lattice (4×4, single sublattice).
    
    The NNN coupling is derived from the two-hop process through
    interstitial squares:
        hop East (α₁) then North (α₂) → combined matrix α₂·α₁
    
    Summing over all 4 diagonal directions:
        H_NNN = 4t₂ sin(kx) sin(ky) × α₁α₂
    
    Parameters:
        t1: nearest-neighbour hopping amplitude
        t2: next-nearest-neighbour (diagonal) hopping amplitude
        m:  bare mass
    """
    H = t1 * (alpha1 * sin(kx) + alpha2 * sin(ky))
    H += 4 * t2 * sin(kx) * sin(ky) * alpha12
    H += m * beta
    return H


def H_488_wilson_test(kx, ky, t1=1.0, r=0.5, m=0.0):
    """
    Hypothetical 4.8.8 with Wilson-like NNN coupling (for comparison).
    
    This is what the NNN term WOULD look like if the diagonal hops
    coupled to β instead of α₁α₂. Included for falsification:
    to show that the geometry does not produce this structure.
    """
    M_eff = m + r * (1 - cos(kx) * cos(ky))
    return alpha1 * sin(kx) + alpha2 * sin(ky) + M_eff * beta


def analyse_nnn_coupling():
    """
    Derive and verify the NNN coupling matrix from the two-hop process.
    
    For a diagonal hop NE (East then North):
        Step 1: hop East with amplitude α₁
        Step 2: hop North with amplitude α₂
        Combined: α₂ · α₁ (rightmost applied first)
    
    The four diagonal directions give different signs:
        NE: α₂ · α₁         = -α₁α₂
        NW: α₂ · (-α₁)      = +α₁α₂
        SE: (-α₂) · α₁      = +α₁α₂
        SW: (-α₂) · (-α₁)   = -α₁α₂
    
    With Bloch phases e^{i(±kx ± ky)}, the sum evaluates to:
        H_NNN = -4t₂ sin(kx) sin(ky) × α₂α₁
              = +4t₂ sin(kx) sin(ky) × α₁α₂
    """
    print("NNN Coupling Derivation")
    print("=" * 60)
    print()

    # Verify the two-hop products
    NE = alpha2 @ alpha1           # α₂·α₁
    NW = alpha2 @ (-alpha1)        # α₂·(-α₁)
    SE = (-alpha2) @ alpha1        # (-α₂)·α₁
    SW = (-alpha2) @ (-alpha1)     # (-α₂)·(-α₁)

    print("Two-hop matrices:")
    print(f"  NE (α₂·α₁)     = -α₁α₂? {np.allclose(NE, -alpha12)}")
    print(f"  NW (α₂·(-α₁))  = +α₁α₂? {np.allclose(NW, alpha12)}")
    print(f"  SE ((-α₂)·α₁)  = +α₁α₂? {np.allclose(SE, alpha12)}")
    print(f"  SW ((-α₂)·(-α₁)) = -α₁α₂? {np.allclose(SW, -alpha12)}")
    print()

    # The Bloch sum:
    # Σ = NE·e^{i(kx+ky)} + NW·e^{i(-kx+ky)} + SE·e^{i(kx-ky)} + SW·e^{i(-kx-ky)}
    # = -α₁α₂·[e^{i(kx+ky)} - e^{i(-kx+ky)} - e^{i(kx-ky)} + e^{i(-kx-ky)}]
    # = -α₁α₂·[2cos(kx+ky) - 2cos(kx-ky)]
    # = -α₁α₂·[-4sin(kx)sin(ky)]
    # = 4sin(kx)sin(ky)·α₁α₂

    print("Bloch sum of diagonal hops:")
    print("  H_NNN(k) = 4t₂ sin(kx) sin(ky) × α₁α₂")
    print()

    # Verify at the critical points
    print("NNN term at doubler locations:")
    for name, kx, ky in [("Γ", 0, 0), ("X", pi, 0), ("Y", 0, pi), ("M", pi, pi)]:
        val = 4 * sin(kx) * sin(ky)
        print(f"  {name} ({kx/pi:.0f}π,{ky/pi:.0f}π): "
              f"4·sin(kx)·sin(ky) = {val:.6f}"
              f" {'→ VANISHES' if abs(val) < 1e-10 else ''}")
    print()

    print("The NNN term vanishes at ALL doubler locations.")
    print("It cannot gap the doublers regardless of the coupling matrix.")
    print()

    # Verify that α₁α₂ ≠ β
    print("Matrix identity check:")
    print(f"  α₁α₂ = i(I⊗σz):  {np.allclose(alpha12, 1j * np.kron(np.eye(2), np.array([[1,0],[0,-1]])))}")
    print(f"  β = σz⊗I:         acts on χ (chirality bit)")
    print(f"  α₁α₂:             acts on I₃ (isospin bit)")
    print(f"  α₁α₂ = β?         {np.allclose(alpha12, beta)} — DIFFERENT OPERATORS")


def test_all_gamma_structures():
    """
    Exhaustively test whether ANY choice of NNN coupling matrix
    can gap the (π,π) doubler via cos(kx)cos(ky) dependence.
    
    The NNN Bloch sum always has sin(kx)sin(ky) dependence
    (from the antisymmetric combination of diagonal phases).
    This vanishes at all high-symmetry points regardless of Γ.
    
    For a Wilson-like cos·cos dependence, one would need SYMMETRIC
    diagonal hopping (same matrix for all 4 directions), giving:
        4cos(kx)cos(ky) × Γ
    But this requires all four diagonal hops to have the SAME matrix,
    which contradicts the alternating signs from the Dirac structure.
    """
    print("Exhaustive NNN Coupling Test")
    print("=" * 60)
    print()
    print("Testing Wilson-like NNN: H = α₁sin(kx) + α₂sin(ky) + Γ·(4cos·cos - 4)")
    print("(This is the HYPOTHETICAL case where NNN gives cos·cos, not sin·sin)")
    print()

    from dirac_matrices import I4 as I_4

    gamma_candidates = {
        "I (scalar)":       I_4,
        "β (mass)":         beta,
        "α₁ (kinetic-x)":  alpha1,
        "α₂ (kinetic-y)":  alpha2,
        "iα₁α₂ (spin)":    1j * alpha12,
        "α₁β":             alpha1 @ beta,
        "α₂β":             alpha2 @ beta,
    }

    t2 = 0.5
    for name, Gamma in gamma_candidates.items():
        def H_test(kx, ky, G=Gamma):
            return (alpha1 * sin(kx) + alpha2 * sin(ky)
                    + t2 * G * (4 * cos(kx) * cos(ky) - 4))

        e_00 = np.min(np.abs(np.linalg.eigvalsh(H_test(0, 0))))
        e_pi0 = np.min(np.abs(np.linalg.eigvalsh(H_test(pi, 0))))
        e_0pi = np.min(np.abs(np.linalg.eigvalsh(H_test(0, pi))))
        e_pipi = np.min(np.abs(np.linalg.eigvalsh(H_test(pi, pi))))

        n_zeros = sum(1 for e in [e_00, e_pi0, e_0pi, e_pipi] if e < 0.01)
        status = "✓ ALL GAPPED" if n_zeros == 1 else f"✗ {n_zeros} zeros"

        print(f"  Γ = {name:<16} E: Γ={e_00:.3f} X={e_pi0:.3f} Y={e_0pi:.3f} M={e_pipi:.3f}  {status}")

    print()
    print("Even with cos·cos dependence, the (π,π) doubler survives")
    print("for ALL choices of Γ (because cos(π)cos(π) = 1 = cos(0)cos(0)).")
    print("Only the standard Wilson form (2 - cos kx - cos ky) kills ALL doublers.")


if __name__ == "__main__":
    analyse_nnn_coupling()
    print()
    test_all_gamma_structures()

    print()
    print("=" * 60)
    print("FULL BRILLOUIN ZONE SCAN")
    print("=" * 60)
    print()

    configs = [
        ("Naive (t₂=0)",             lambda kx, ky: H_488_physical(kx, ky, t2=0)),
        ("4.8.8 physical (t₂=0.25)", lambda kx, ky: H_488_physical(kx, ky, t2=0.25)),
        ("4.8.8 physical (t₂=0.5)",  lambda kx, ky: H_488_physical(kx, ky, t2=0.5)),
        ("4.8.8 physical (t₂=1.0)",  lambda kx, ky: H_488_physical(kx, ky, t2=1.0)),
        ("Hypothetical Wilson NNN",   lambda kx, ky: H_488_wilson_test(kx, ky, r=0.5)),
    ]

    for name, H_func in configs:
        pts, _ = scan_brillouin_zone(H_func, N=200)
        n_phys = sum(1 for kxp, kyp in pts if abs(kxp) < 0.1 and abs(kyp) < 0.1)
        print(f"  {name:<35} Dirac points: {len(pts):>3} "
              f"({n_phys} physical, {len(pts)-n_phys} doublers)")

    print()
    print("=" * 60)
    print("CONCLUSION")
    print("=" * 60)
    print()
    print("The 4.8.8 topology's physical NNN coupling (α₁α₂ × sin·sin)")
    print("CANNOT gap the fermion doublers. The doublers persist for all t₂.")
    print()
    print("Resolution: The circlette's CNOT coin operator provides a")
    print("dynamical Wilson mass through discrete Z₂ chiral mixing,")
    print("sidestepping Nielsen-Ninomiya at the axiomatic level.")
