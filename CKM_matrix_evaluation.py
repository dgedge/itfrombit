#!/usr/bin/env python3
"""
CKM Matrix from the Circlette 8-Bit Code
==========================================

The CKM matrix is the mismatch between mass eigenstates and weak
eigenstates in the quark sector. In the Standard Model, it's a 3Ã3
unitary matrix parameterised by 3 angles and 1 CP-violating phase.

Experimental CKM matrix (magnitudes):
    |V_CKM| â [[0.974, 0.225, 0.004],
                [0.225, 0.973, 0.041],
                [0.009, 0.040, 0.999]]

Key features: nearly diagonal, small off-diagonal elements,
hierarchical structure |V_us| >> |V_cb| >> |V_ub|.

Strategy:
=========
1. The CNOT gate acts on bits (LQ, Iâ) â flips Iâ when LQ=1.
   This is the weak interaction: it converts uâd within a generation.
   In the WEAK basis, generations don't mix.

2. The MASS eigenstates are the eigenstates of the full Hamiltonian,
   which includes not just the CNOT but also the generation-dependent
   mass term. The mass term must come from how the generation bits
   (Gâ, Gâ) affect the CNOT frequency.

3. If the generation sector has Zâ symmetry (as claimed via Koide),
   then the mass matrix in generation space should reflect this
   symmetry. The most general Zâ-symmetric mass matrix is the
   "democratic matrix" plus a diagonal perturbation.

4. The CKM matrix is V = U_uâ  Â· U_d, where U_u and U_d diagonalise
   the up-type and down-type mass matrices respectively.

Author: D.G. Elliman / Neuro-Symbolic Ltd
Date: February 2026
"""

import numpy as np
from fractions import Fraction
from itertools import product
import scipy.linalg as la

# ============================================================
# 1. The generation sector and its symmetries
# ============================================================

def generation_analysis():
    """Analyse the Zâ structure of the generation bits."""
    
    print("â" * 70)
    print("GENERATION SECTOR ANALYSIS")
    print("â" * 70)
    print()
    
    # Generation encoding
    gens = [(0, 0, 1), (1, 0, 2), (0, 1, 3)]  # (Gâ, Gâ, gen_number)
    # (1,1) is forbidden by R1
    
    print("  Generation encoding on the ring:")
    print("    (Gâ, Gâ) = (0,0) â Gen 1  (e.g. u, d, e, Îœâ)")
    print("    (Gâ, Gâ) = (1,0) â Gen 2  (e.g. c, s, ÎŒ, ÎœÎŒ)")
    print("    (Gâ, Gâ) = (0,1) â Gen 3  (e.g. t, b, Ï, ÎœÏ)")
    print("    (Gâ, Gâ) = (1,1) â FORBIDDEN by R1")
    print()
    
    # The Zâ symmetry: binary counting mod 3
    # (0,0) â (1,0) â (0,1) â (1,1)=forbidden â back to (0,0)
    # This is NOT a clean cyclic group in binary.
    # The increment operation in binary is:
    #   (0,0) + 1 = (1,0)    [flip Gâ]
    #   (1,0) + 1 = (0,1)    [flip Gâ, carry to Gâ]
    #   (0,1) + 1 = (1,1)    [forbidden â wraps to (0,0)]
    
    print("  Binary increment (Zâ generator):")
    print("    Gen 1 â Gen 2:  flip Gâ           (0,0)â(1,0)")
    print("    Gen 2 â Gen 3:  flip Gâ, carry Gâ  (1,0)â(0,1)")
    print("    Gen 3 â Gen 1:  wrap (mod R1)       (0,1)â(0,0)")
    print()
    
    # The Zâ generator as a matrix on the 3D generation space
    # |gen 1â© â |gen 2â© â |gen 3â© â |gen 1â©
    Z3 = np.array([[0, 0, 1],
                    [1, 0, 0],
                    [0, 1, 0]], dtype=complex)
    
    eigenvalues = np.linalg.eigvals(Z3)
    print(f"  Zâ generator eigenvalues: {eigenvalues}")
    print(f"  = 1, Ï, ÏÂ²  where Ï = exp(2Ïi/3)")
    print()
    
    return Z3


# ============================================================
# 2. The CNOT as weak interaction in generation space
# ============================================================

def weak_basis_construction():
    """
    In the WEAK basis, the W boson couples u_L â d_L within each
    generation without mixing generations. The weak eigenstates are:
    
      (uâ, uâ, uâ) and (dâ, dâ, dâ)
    
    where the subscript is the generation index.
    The weak interaction vertex is diagonal: W couples u_n to d_n only.
    """
    
    print("â" * 70)
    print("WEAK BASIS (CNOT INTERACTION BASIS)")
    print("â" * 70)
    print()
    print("  The CNOT gate flips Iâ when LQ=1, INDEPENDENTLY of (Gâ,Gâ).")
    print("  In the weak basis, generations are uncoupled:")
    print("    Wâº: d_n â u_n  (for each n = 1,2,3)")
    print("    Wâ»: u_n â d_n  (for each n = 1,2,3)")
    print()
    print("  Weak coupling matrix (3Ã3 in generation space):")
    print("    V_weak = Iâ = diag(1, 1, 1)")
    print()
    print("  If the mass eigenstates were identical to the weak eigenstates,")
    print("  the CKM matrix would be the identity (no mixing).")
    print("  CKM mixing requires the mass matrix to be non-diagonal in")
    print("  generation space.")
    print()


# ============================================================
# 3. Mass matrix from Zâ symmetry
# ============================================================

def democratic_mass_matrix():
    """
    The Koide relation implies the mass matrix has Zâ symmetry.
    The most general Zâ-symmetric matrix is:
    
        M = aÂ·I + bÂ·D + cÂ·Zâ + c*Â·Zââ 
    
    where D is the democratic matrix D_ij = 1 for all i,j,
    and Zâ is the cyclic permutation matrix.
    
    In the Koide parametrisation:
        âm_n = A + BÂ·cos(Îž + 2Ïn/3)
    
    The mass-squared matrix Mâ M has eigenvectors that are the
    Fourier modes of the Zâ group.
    """
    
    print("â" * 70)
    print("MASS MATRIX FROM Zâ SYMMETRY")
    print("â" * 70)
    print()
    
    omega = np.exp(2j * np.pi / 3)
    
    # The Zâ Fourier basis (generation mass eigenstates)
    # These are the eigenvectors of ANY Zâ-symmetric matrix
    F = np.array([
        [1, 1, 1],             # Zâ eigenvalue 1
        [1, omega, omega**2],   # Zâ eigenvalue Ï
        [1, omega**2, omega],   # Zâ eigenvalue ÏÂ²
    ], dtype=complex) / np.sqrt(3)
    
    print("  Zâ Fourier basis (normalised):")
    print("    |fââ© = (1, 1, 1)/â3           â 'democratic' mode")
    print("    |fââ© = (1, Ï, ÏÂ²)/â3          â 1st harmonic")
    print("    |fââ© = (1, ÏÂ², Ï)/â3          â 2nd harmonic")
    print(f"    where Ï = exp(2Ïi/3) = {omega:.4f}")
    print()
    
    # If the mass matrix is Zâ-symmetric, its eigenvectors ARE the
    # Fourier modes, and the eigenvalues are the Fourier coefficients.
    # The diagonalising matrix is U = Fâ .
    
    # For the Koide parametrisation âm_n = A + BÂ·cos(Îž + 2Ïn/3):
    # This is a Zâ-symmetric function. The mass matrix in the
    # generation basis has the form:
    #   M_ij = Î£_n âm_n Â· (Fâ )_ni Â· (F)_nj
    
    print("  Key insight: if both up-type and down-type mass matrices")
    print("  are Zâ-symmetric, they share the SAME eigenvectors (the")
    print("  Fourier modes). Then:")
    print("    U_u = Fâ   and  U_d = Fâ ")
    print("    V_CKM = U_uâ  Â· U_d = F Â· Fâ  = I")
    print()
    print("  â EXACT Zâ SYMMETRY GIVES V_CKM = IDENTITY (no mixing)")
    print()
    print("  This means: CKM mixing requires Zâ BREAKING.")
    print("  The question is: what breaks Zâ differently for up-type")
    print("  and down-type quarks?")
    print()
    
    return F


# ============================================================
# 4. Zâ breaking from the CNOT structure
# ============================================================

def z3_breaking_analysis():
    """
    Investigate sources of Zâ breaking in the circlette code.
    
    The generation bits (Gâ, Gâ) are adjacent on the ring:
        ... - Gâ - Gâ - Câ - Câ - LQ - Iâ - Ï - W - (back to Gâ)
    
    The CNOT acts on (LQ, Iâ): flips Iâ when LQ=1.
    This is generation-independent â no direct Zâ breaking.
    
    But the MASS is the CNOT frequency, and the frequency might
    depend on the local constraint environment of each generation.
    
    Possible Zâ breaking mechanisms:
    1. The (1,1) exclusion (R1) breaks the Zâ of 2-bit space to Zâ.
       This introduces an asymmetry: gen 1 = (0,0) is the "origin",
       while gens 2 and 3 are related by GââGâ swap.
    
    2. The ring topology creates nearest-neighbour interactions.
       Gâ is adjacent to W (closing the ring) and to Gâ.
       Gâ is adjacent to Gâ and to Câ.
       These are DIFFERENT environments: Gâ touches the electroweak
       sector, Gâ touches the colour sector.
    
    3. Higher-order effects: the CNOT frequency for a given generation
       depends on the quantum walk's overlap with other generations,
       which depends on the generation-dependent phases.
    """
    
    print("â" * 70)
    print("Zâ BREAKING MECHANISMS")
    print("â" * 70)
    print()
    
    print("  The ring topology:")
    print("  ... W â Gâ â Gâ â Câ â Câ â LQ â Iâ â Ï â W ...")
    print()
    print("  Gâ neighbours: W (electroweak) and Gâ (generation)")
    print("  Gâ neighbours: Gâ (generation) and Câ (colour)")
    print()
    print("  This asymmetry means Gâ and Gâ are NOT equivalent,")
    print("  even though the Zâ cycles them symmetrically.")
    print()
    
    # The key Zâ-breaking effect comes from the ring adjacency.
    # The constraint R1 excludes (1,1). This creates a potential
    # energy landscape in generation space:
    #   (0,0): neither bit is "on" â lowest energy (?) 
    #   (1,0): Gâ on, adjacent to W sector â intermediate
    #   (0,1): Gâ on, adjacent to C sector â highest energy (?)
    
    # This is speculative but physically motivated: the generation
    # bits couple to different sectors through ring adjacency,
    # so their "cost" (contribution to the mass) differs.
    
    print("  Hypothesis: ring adjacency creates a generation-dependent")
    print("  mass correction:")
    print("    Gen 1 (0,0): no generation bits active â lightest")
    print("    Gen 2 (1,0): Gâ active, couples to W â intermediate")
    print("    Gen 3 (0,1): Gâ active, couples to C â heaviest")
    print()
    print("  This gives a DIAGONAL but non-degenerate mass matrix,")
    print("  which still produces V_CKM = I (no mixing).")
    print()
    print("  For off-diagonal mixing, we need TRANSITIONS between")
    print("  generations, i.e. processes that flip Gâ or Gâ.")
    print()


# ============================================================
# 5. Generation transitions from ring topology
# ============================================================

def generation_transition_matrix():
    """
    On the 8-bit ring, the generation bits Gâ and Gâ are adjacent.
    Information propagating around the ring creates effective 
    couplings between the generation sector and other sectors.
    
    The CNOT gate acts on (LQ, Iâ), but during the quantum walk,
    the wavefunction spreads across all bits. Virtual processes
    (off-shell intermediate states) can flip generation bits.
    
    The effective generation transition operator comes from the
    ring topology: 
    - Flipping Gâ alone: gen 1 â gen 2  (and gen 3 â forbidden)
    - Flipping Gâ alone: gen 1 â gen 3  (and gen 2 â forbidden)
    - Flipping both: gen 2 â gen 3       (and gen 1 â forbidden)
    
    R1 forbids (1,1), so transitions TO (1,1) are suppressed.
    This creates an asymmetric transition matrix.
    """
    
    print("â" * 70)
    print("GENERATION TRANSITION OPERATORS")
    print("â" * 70)
    print()
    
    # Label: gen 1 = (0,0), gen 2 = (1,0), gen 3 = (0,1)
    # Index: 0, 1, 2
    
    # Ï_x on Gâ (flip Gâ):
    # (0,0) â (1,0): gen1 â gen2  â
    # (1,0) â (0,0): gen2 â gen1  â
    # (0,1) â (1,1): gen3 â FORBIDDEN
    # (1,1) â (0,1): FORBIDDEN â gen3
    
    # Projected onto the valid 3-generation subspace:
    # Ï_x^Gâ flips gen 1 â gen 2 but ANNIHILATES gen 3
    
    sigma_G0 = np.array([
        [0, 1, 0],  # gen1 â gen2
        [1, 0, 0],  # gen2 â gen1
        [0, 0, 0],  # gen3 â forbidden (projected out)
    ], dtype=complex)
    
    # Ï_x on Gâ (flip Gâ):
    # (0,0) â (0,1): gen1 â gen3  â
    # (0,1) â (0,0): gen3 â gen1  â
    # (1,0) â (1,1): gen2 â FORBIDDEN
    # (1,1) â (1,0): FORBIDDEN â gen2
    
    sigma_G1 = np.array([
        [0, 0, 1],  # gen1 â gen3
        [0, 0, 0],  # gen2 â forbidden (projected out)
        [1, 0, 0],  # gen3 â gen1
    ], dtype=complex)
    
    # Combined flip (both Gâ and Gâ):
    # (0,0) â (1,1): gen1 â FORBIDDEN
    # (1,0) â (0,1): gen2 â gen3  â
    # (0,1) â (1,0): gen3 â gen2  â
    # (1,1) â (0,0): FORBIDDEN â gen1
    
    sigma_G0G1 = np.array([
        [0, 0, 0],  # gen1 â forbidden
        [0, 0, 1],  # gen2 â gen3
        [0, 1, 0],  # gen3 â gen2
    ], dtype=complex)
    
    print("  Ï_Gâ (flip Gâ, projected onto valid subspace):")
    print(f"    {sigma_G0}")
    print(f"    Connects: gen 1 â gen 2 only (gen 3 isolated)")
    print()
    print("  Ï_Gâ (flip Gâ, projected):")
    print(f"    {sigma_G1}")
    print(f"    Connects: gen 1 â gen 3 only (gen 2 isolated)")
    print()
    print("  Ï_GâGâ (flip both, projected):")
    print(f"    {sigma_G0G1}")
    print(f"    Connects: gen 2 â gen 3 only (gen 1 isolated)")
    print()
    
    # KEY OBSERVATION: each transition operator connects exactly
    # one pair of generations, and the THIRD generation is isolated.
    # The R1 constraint (forbidding (1,1)) is what creates this
    # structure â it prevents certain transitions.
    
    print("  â KEY: R1 creates PARTIAL transitions â each operator")
    print("    connects only 2 of 3 generations. This is precisely")
    print("    the structure needed for hierarchical mixing:")
    print("    â¢ Ï_Gâ controls 1â2 mixing (Cabibbo angle)")
    print("    â¢ Ï_Gâ controls 1â3 mixing (V_ub)")
    print("    â¢ Ï_GâGâ controls 2â3 mixing (V_cb)")
    print()
    
    return sigma_G0, sigma_G1, sigma_G0G1


# ============================================================
# 6. Effective mass matrix construction
# ============================================================

def build_mass_matrix(epsilon_01, epsilon_10, epsilon_11, 
                       m1, m2, m3, phase_01=0, phase_10=0, phase_11=0):
    """
    Build the most general mass matrix consistent with the 
    circlette ring topology.
    
    M = diag(m1, m2, m3) + ÎµââÂ·e^{iÏââ}Â·Ï_Gâ 
                          + ÎµââÂ·e^{iÏââ}Â·Ï_Gâ 
                          + ÎµââÂ·e^{iÏââ}Â·Ï_GâGâ
    
    The diagonal part comes from the base CNOT frequency.
    The off-diagonal part comes from virtual generation transitions.
    """
    
    sigma_G0 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=complex)
    sigma_G1 = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]], dtype=complex)
    sigma_G0G1 = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]], dtype=complex)
    
    M = np.diag([m1, m2, m3]).astype(complex)
    M += epsilon_01 * np.exp(1j * phase_01) * sigma_G0
    M += epsilon_10 * np.exp(1j * phase_10) * sigma_G1
    M += epsilon_11 * np.exp(1j * phase_11) * sigma_G0G1
    
    return M


def ckm_from_mass_matrices(M_u, M_d):
    """
    Compute the CKM matrix from up-type and down-type mass matrices.
    
    V_CKM = U_uâ  Â· U_d
    
    where M_f = U_f Â· diag(m_f1, m_f2, m_f3) Â· U_fâ 
    """
    
    # Diagonalise Mâ M (mass-squared matrix, always Hermitian)
    MuMu = M_u @ M_u.conj().T
    MdMd = M_d @ M_d.conj().T
    
    eigvals_u, U_u = np.linalg.eigh(MuMu)
    eigvals_d, U_d = np.linalg.eigh(MdMd)
    
    # Sort by eigenvalue (ascending = lightest first)
    idx_u = np.argsort(eigvals_u)
    idx_d = np.argsort(eigvals_d)
    U_u = U_u[:, idx_u]
    U_d = U_d[:, idx_d]
    
    # CKM matrix
    V_ckm = U_u.conj().T @ U_d
    
    # Fix phase convention: make diagonal elements real and positive
    for i in range(3):
        phase = np.angle(V_ckm[i, i])
        V_ckm[i, :] *= np.exp(-1j * phase)
    
    masses_u = np.sqrt(np.abs(eigvals_u[idx_u]))
    masses_d = np.sqrt(np.abs(eigvals_d[idx_d]))
    
    return V_ckm, masses_u, masses_d


# ============================================================
# 7. Physical parameter exploration
# ============================================================

def ring_coupling_model():
    """
    Model where the off-diagonal couplings are determined by ring
    adjacency distances.
    
    On the ring: Gâ - Gâ - Câ - Câ - LQ - Iâ - Ï - W - (Gâ)
    
    Gâ is 1 hop from Gâ, and 1 hop from W (closing the ring).
    Gâ is 1 hop from Gâ, and 1 hop from Câ.
    
    The effective coupling for flipping a generation bit depends on
    how strongly that bit couples to the dynamic sector (Iâ, Ï).
    
    Path from Gâ to Iâ: Gâ â W â Ï â Iâ  (3 hops via ring closure)
                     or: Gâ â Gâ â Câ â Câ â LQ â Iâ  (5 hops)
    
    Path from Gâ to Iâ: Gâ â Câ â Câ â LQ â Iâ  (4 hops)
                     or: Gâ â Gâ â W â Ï â Iâ    (4 hops)
    
    If coupling strength falls exponentially with distance:
      Îµ ~ e^{-d/ÎŸ}
    
    then Gâ (3 hops to Iâ) couples more strongly than Gâ (4 hops).
    """
    
    print("â" * 70)
    print("RING ADJACENCY COUPLING MODEL")
    print("â" * 70)
    print()
    
    print("  Ring: Gâ â Gâ â Câ â Câ â LQ â Iâ â Ï â W â (Gâ)")
    print()
    print("  Shortest path distances to Iâ (the CNOT target):")
    
    # Ring distances (both directions)
    bits = ['Gâ', 'Gâ', 'Câ', 'Câ', 'LQ', 'Iâ', 'Ï', 'W']
    I3_idx = 5
    for i, name in enumerate(bits):
        d_cw = (I3_idx - i) % 8   # clockwise
        d_ccw = (i - I3_idx) % 8   # counterclockwise
        d = min(d_cw, d_ccw)
        print(f"    {name:3s} â Iâ: {d} hops (min of {d_cw} CW, {d_ccw} CCW)")
    
    print()
    
    # Gâ is 3 hops from Iâ (via W, Ï)
    # Gâ is 3 hops from Iâ (via Câ, Câ, LQ or via Gâ, W, Ï â 4 hops)
    # Wait: let me recount
    # Gâ(0) â Iâ(5): CW = 5, CCW = 3 â min = 3
    # Gâ(1) â Iâ(5): CW = 4, CCW = 4 â min = 4
    
    d_G0 = 3  # Gâ to Iâ
    d_G1 = 4  # Gâ to Iâ
    
    print(f"  Gâ to Iâ: {d_G0} hops")
    print(f"  Gâ to Iâ: {d_G1} hops")
    print()
    
    # Generation transition couplings:
    # Ï_Gâ (flip Gâ): requires Gâ to couple to the dynamic sector
    #   â coupling â e^{-d_Gâ/ÎŸ}
    # Ï_Gâ (flip Gâ): requires Gâ to couple
    #   â coupling â e^{-d_Gâ/ÎŸ}
    # Ï_GâGâ (flip both): requires BOTH to couple
    #   â coupling â e^{-(d_Gâ+d_Gâ)/ÎŸ}
    
    print("  Coupling strengths (exponential decay with ring distance):")
    print(f"    Îµââ (Ï_Gâ):    ~ exp(-{d_G0}/ÎŸ)")
    print(f"    Îµââ (Ï_Gâ):    ~ exp(-{d_G1}/ÎŸ)")
    print(f"    Îµââ (Ï_GâGâ):  ~ exp(-{d_G0 + d_G1}/ÎŸ)")
    print()
    print("  This gives a HIERARCHY: Îµââ > Îµââ > Îµââ")
    print("  Which matches the CKM hierarchy: |V_us| > |V_ub| > ... ")
    print("  Wait â experimentally: |V_us| > |V_cb| > |V_ub|")
    print("  Our model gives: Îµââ > Îµââ > Îµââ")
    print("  Mapping: Îµââ ~ V_us, Îµââ ~ V_ub, Îµââ ~ V_cb")
    print()
    print("  Experimental hierarchy: |V_us| â 0.225 >> |V_cb| â 0.041 >> |V_ub| â 0.004")
    print("  Our hierarchy:          Îµââ > Îµââ > Îµââ")
    print("  Mismatch: we predict Îµââ > Îµââ but experiment has |V_cb| > |V_ub|")
    print("  i.e. 2â3 mixing > 1â3 mixing, but we predict the opposite.")
    print()
    print("  HOWEVER: V_ub and V_cb are EFFECTIVE â they include")
    print("  indirect paths (1â2â3 contributes to effective 1â3).")
    print("  The CKM matrix elements are NOT the bare couplings.")
    print()
    
    return d_G0, d_G1


def numerical_exploration():
    """
    Numerically explore the parameter space to find what coupling
    values reproduce the experimental CKM.
    """
    
    print("â" * 70)
    print("NUMERICAL CKM EXPLORATION")
    print("â" * 70)
    print()
    
    # Experimental CKM magnitudes (PDG 2024)
    V_exp = np.array([
        [0.97435, 0.22500, 0.00369],
        [0.22486, 0.97349, 0.04182],
        [0.00857, 0.04110, 0.99912]
    ])
    
    print("  Experimental |V_CKM| (PDG 2024):")
    for i in range(3):
        print(f"    [{V_exp[i,0]:.5f}  {V_exp[i,1]:.5f}  {V_exp[i,2]:.5f}]")
    print()
    
    # Experimental quark masses (MS-bar at 2 GeV, in MeV)
    m_u_phys = np.array([2.16, 1270, 172760])   # u, c, t (MeV)
    m_d_phys = np.array([4.67, 93.4, 4180])     # d, s, b (MeV)
    
    print("  Experimental quark masses (MeV):")
    print(f"    Up-type:   u = {m_u_phys[0]},  c = {m_u_phys[1]},  t = {m_u_phys[2]}")
    print(f"    Down-type: d = {m_d_phys[0]},  s = {m_d_phys[1]},  b = {m_d_phys[2]}")
    print()
    
    # ââ Model A: Pure ring-distance coupling ââ
    print("â" * 70)
    print("  MODEL A: Ring-distance exponential coupling")
    print("â" * 70)
    print()
    
    # Use physical masses as diagonal, explore off-diagonal
    # The off-diagonal couplings are small perturbations
    
    best_fit = None
    best_error = 1e10
    
    for xi in np.linspace(0.5, 5.0, 50):
        eps_12 = np.exp(-3 / xi)   # Gâ flip, d=3
        eps_13 = np.exp(-4 / xi)   # Gâ flip, d=4
        eps_23 = np.exp(-7 / xi)   # both flip, d=7
        
        for scale in np.logspace(-2, 2, 50):
            for delta_phase in np.linspace(0, 2*np.pi, 20):
                e12 = scale * eps_12
                e13 = scale * eps_13
                e23 = scale * eps_23
                
                # Build mass matrices (use sqrt of physical masses for Yukawa)
                M_u = build_mass_matrix(e12, e13, e23, 
                                         m_u_phys[0], m_u_phys[1], m_u_phys[2],
                                         phase_01=0, phase_10=delta_phase, phase_11=0)
                M_d = build_mass_matrix(e12, e13, e23,
                                         m_d_phys[0], m_d_phys[1], m_d_phys[2],
                                         phase_01=0, phase_10=delta_phase, phase_11=0)
                
                V, mu, md = ckm_from_mass_matrices(M_u, M_d)
                V_abs = np.abs(V)
                
                error = np.sum((V_abs - V_exp)**2)
                
                if error < best_error:
                    best_error = error
                    best_fit = {
                        'xi': xi, 'scale': scale, 'phase': delta_phase,
                        'V': V, 'V_abs': V_abs, 'mu': mu, 'md': md,
                        'e12': e12, 'e13': e13, 'e23': e23
                    }
    
    print(f"  Best fit (ÎŸ = {best_fit['xi']:.2f}, scale = {best_fit['scale']:.4f}, "
          f"ÎŽ = {best_fit['phase']:.3f}):")
    print(f"  Couplings: Îµââ = {best_fit['e12']:.6f}, Îµââ = {best_fit['e13']:.6f}, "
          f"Îµââ = {best_fit['e23']:.6f}")
    print()
    print(f"  |V_CKM| (model):")
    for i in range(3):
        print(f"    [{best_fit['V_abs'][i,0]:.5f}  {best_fit['V_abs'][i,1]:.5f}  "
              f"{best_fit['V_abs'][i,2]:.5f}]")
    print()
    print(f"  |V_CKM| (experiment):")
    for i in range(3):
        print(f"    [{V_exp[i,0]:.5f}  {V_exp[i,1]:.5f}  {V_exp[i,2]:.5f}]")
    print()
    print(f"  RMS error: {np.sqrt(best_error/9):.6f}")
    print()
    
    # ââ Model B: Direct fit with ring-motivated structure ââ
    print("â" * 70)
    print("  MODEL B: Direct fit with ring-motivated hierarchy")
    print("â" * 70)
    print()
    
    # Use the STRUCTURE (which generations couple to which) from the
    # ring topology, but fit the coupling strengths freely.
    # The ring topology predicts:
    #   - Ï_Gâ connects gen 1â2 only
    #   - Ï_Gâ connects gen 1â3 only
    #   - Ï_GâGâ connects gen 2â3 only
    # This is the SAME structure as the Wolfenstein parametrisation!
    
    from scipy.optimize import minimize
    
    def ckm_error(params):
        e12, e13, e23, ph12, ph13, ph23 = params
        M_u = build_mass_matrix(e12, e13, e23,
                                 m_u_phys[0], m_u_phys[1], m_u_phys[2],
                                 phase_01=ph12, phase_10=ph13, phase_11=ph23)
        M_d = build_mass_matrix(e12, e13, e23,
                                 m_d_phys[0], m_d_phys[1], m_d_phys[2],
                                 phase_01=ph12, phase_10=ph13, phase_11=ph23)
        V, _, _ = ckm_from_mass_matrices(M_u, M_d)
        return np.sum((np.abs(V) - V_exp)**2)
    
    # Multi-start optimisation
    best_result = None
    best_err = 1e10
    
    np.random.seed(42)
    for trial in range(500):
        x0 = np.random.uniform([-100, -100, -100, 0, 0, 0],
                                 [100, 100, 100, 2*np.pi, 2*np.pi, 2*np.pi])
        try:
            result = minimize(ckm_error, x0, method='Nelder-Mead',
                            options={'maxiter': 5000, 'xatol': 1e-10, 'fatol': 1e-12})
            if result.fun < best_err:
                best_err = result.fun
                best_result = result
        except:
            continue
    
    if best_result is not None:
        e12, e13, e23, ph12, ph13, ph23 = best_result.x
        M_u = build_mass_matrix(e12, e13, e23,
                                 m_u_phys[0], m_u_phys[1], m_u_phys[2],
                                 phase_01=ph12, phase_10=ph13, phase_11=ph23)
        M_d = build_mass_matrix(e12, e13, e23,
                                 m_d_phys[0], m_d_phys[1], m_d_phys[2],
                                 phase_01=ph12, phase_10=ph13, phase_11=ph23)
        V, mu, md = ckm_from_mass_matrices(M_u, M_d)
        
        print(f"  Best fit parameters:")
        print(f"    Îµââ = {e12:.4f} MeV  (Gâ flip: gen 1â2)")
        print(f"    Îµââ = {e13:.4f} MeV  (Gâ flip: gen 1â3)")
        print(f"    Îµââ = {e23:.4f} MeV  (GâGâ flip: gen 2â3)")
        print(f"    Ïââ = {ph12 % (2*np.pi):.4f} rad")
        print(f"    Ïââ = {ph13 % (2*np.pi):.4f} rad")
        print(f"    Ïââ = {ph23 % (2*np.pi):.4f} rad")
        print()
        print(f"  |V_CKM| (model):")
        for i in range(3):
            print(f"    [{np.abs(V[i,0]):.5f}  {np.abs(V[i,1]):.5f}  {np.abs(V[i,2]):.5f}]")
        print()
        print(f"  |V_CKM| (experiment):")
        for i in range(3):
            print(f"    [{V_exp[i,0]:.5f}  {V_exp[i,1]:.5f}  {V_exp[i,2]:.5f}]")
        print()
        print(f"  RMS error: {np.sqrt(best_err/9):.6f}")
        print()
        
        # Check coupling hierarchy
        print(f"  Coupling hierarchy: |Îµââ| = {abs(e12):.4f}, |Îµââ| = {abs(e13):.4f}, |Îµââ| = {abs(e23):.4f}")
        couplings = sorted([(abs(e12), '1â2'), (abs(e13), '1â3'), (abs(e23), '2â3')], reverse=True)
        print(f"  Order: {couplings[0][1]} > {couplings[1][1]} > {couplings[2][1]}")
        
        # Compare with ring distances
        print()
        print(f"  Ring distance prediction: Îµââ (d=3) > Îµââ (d=4) > Îµââ (d=7)")
        print(f"  Fit result:               |Îµââ| = {abs(e12):.2f}, |Îµââ| = {abs(e13):.2f}, |Îµââ| = {abs(e23):.2f}")
        
        if abs(e12) > abs(e13) > abs(e23):
            print(f"  â HIERARCHY MATCHES RING DISTANCE PREDICTION!")
        else:
            print(f"  â Hierarchy does not match ring distance prediction")
        
        # Extract Wolfenstein parameters
        lambda_w = np.abs(V[0, 1])
        A_w = np.abs(V[1, 2]) / lambda_w**2
        
        print()
        print(f"  Wolfenstein parameters:")
        print(f"    Î» = |V_us| = {lambda_w:.5f}  (exp: 0.22500)")
        print(f"    A = |V_cb|/Î»Â² = {A_w:.5f}  (exp: 0.826)")
        
        # CP-violating phase
        J = np.imag(V[0,0] * V[1,1] * np.conj(V[0,1]) * np.conj(V[1,0]))
        print(f"    Jarlskog invariant J = {J:.6e}  (exp: 3.18Ã10â»âµ)")
        print()


# ============================================================
# 8. Separate up/down coupling strengths
# ============================================================

def separate_sector_model():
    """
    Model where up-type and down-type quarks have DIFFERENT
    off-diagonal couplings (because they have different Iâ values
    and thus different CNOT dynamics).
    """
    
    print("â" * 70)
    print("MODEL C: SEPARATE UP/DOWN COUPLING STRENGTHS")
    print("â" * 70)
    print()
    print("  Physical motivation: up-type (Iâ=0) and down-type (Iâ=1)")
    print("  quarks experience the CNOT differently. The CNOT flips Iâ")
    print("  when LQ=1, so uâd conversion occurs at each tick.")
    print("  The generation coupling strength may depend on Iâ.")
    print()
    
    from scipy.optimize import minimize
    
    m_u_phys = np.array([2.16, 1270, 172760])
    m_d_phys = np.array([4.67, 93.4, 4180])
    
    V_exp = np.array([
        [0.97435, 0.22500, 0.00369],
        [0.22486, 0.97349, 0.04182],
        [0.00857, 0.04110, 0.99912]
    ])
    
    def ckm_error_sep(params):
        eu12, eu13, eu23, ed12, ed13, ed23, ph_u, ph_d = params
        M_u = build_mass_matrix(eu12, eu13, eu23,
                                 m_u_phys[0], m_u_phys[1], m_u_phys[2],
                                 phase_10=ph_u)
        M_d = build_mass_matrix(ed12, ed13, ed23,
                                 m_d_phys[0], m_d_phys[1], m_d_phys[2],
                                 phase_10=ph_d)
        V, _, _ = ckm_from_mass_matrices(M_u, M_d)
        return np.sum((np.abs(V) - V_exp)**2)
    
    best_result = None
    best_err = 1e10
    
    np.random.seed(123)
    for trial in range(1000):
        x0 = np.random.uniform(
            [-500, -500, -500, -500, -500, -500, 0, 0],
            [500, 500, 500, 500, 500, 500, 2*np.pi, 2*np.pi])
        try:
            result = minimize(ckm_error_sep, x0, method='Nelder-Mead',
                            options={'maxiter': 10000, 'xatol': 1e-12, 'fatol': 1e-14})
            if result.fun < best_err:
                best_err = result.fun
                best_result = result
        except:
            continue
    
    if best_result is not None:
        eu12, eu13, eu23, ed12, ed13, ed23, ph_u, ph_d = best_result.x
        M_u = build_mass_matrix(eu12, eu13, eu23,
                                 m_u_phys[0], m_u_phys[1], m_u_phys[2],
                                 phase_10=ph_u)
        M_d = build_mass_matrix(ed12, ed13, ed23,
                                 m_d_phys[0], m_d_phys[1], m_d_phys[2],
                                 phase_10=ph_d)
        V, mu, md = ckm_from_mass_matrices(M_u, M_d)
        
        print(f"  Best fit RMS error: {np.sqrt(best_err/9):.8f}")
        print()
        print(f"  Up-type couplings:   Îµââ={eu12:.2f}, Îµââ={eu13:.2f}, Îµââ={eu23:.2f} MeV")
        print(f"  Down-type couplings: Îµââ={ed12:.2f}, Îµââ={ed13:.2f}, Îµââ={ed23:.2f} MeV")
        print(f"  Phases: Ï_u={ph_u%(2*np.pi):.4f}, Ï_d={ph_d%(2*np.pi):.4f} rad")
        print()
        print(f"  |V_CKM| (model):")
        for i in range(3):
            print(f"    [{np.abs(V[i,0]):.5f}  {np.abs(V[i,1]):.5f}  {np.abs(V[i,2]):.5f}]")
        print()
        print(f"  |V_CKM| (experiment):")
        for i in range(3):
            print(f"    [{V_exp[i,0]:.5f}  {V_exp[i,1]:.5f}  {V_exp[i,2]:.5f}]")
        print()
        
        # Jarlskog invariant (CP violation measure)
        J = np.imag(V[0,0] * V[1,1] * np.conj(V[0,1]) * np.conj(V[1,0]))
        print(f"  Jarlskog invariant J = {J:.6e}  (exp: 3.18Ã10â»âµ)")
        print()
        
        # Unitarity check
        VVdag = V @ V.conj().T
        print(f"  Unitarity check |VÂ·Vâ  - I|:")
        unitarity_err = np.max(np.abs(VVdag - np.eye(3)))
        print(f"    Max deviation: {unitarity_err:.2e}")
        print()
        
        # Check mass eigenvalues
        print(f"  Mass eigenvalues (MeV):")
        print(f"    Up-type:   {mu[0]:.2f}, {mu[1]:.1f}, {mu[2]:.0f}  "
              f"(exp: 2.16, 1270, 172760)")
        print(f"    Down-type: {md[0]:.2f}, {md[1]:.1f}, {md[2]:.0f}  "
              f"(exp: 4.67, 93.4, 4180)")


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    
    print()
    print("âââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ")
    print("â  CKM MATRIX FROM THE CIRCLETTE 8-BIT CODE                      â")
    print("â  D.G. Elliman â Neuro-Symbolic Ltd â February 2026              â")
    print("âââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââââ")
    print()
    
    # Analysis steps
    Z3 = generation_analysis()
    weak_basis_construction()
    F = democratic_mass_matrix()
    z3_breaking_analysis()
    sigmas = generation_transition_matrix()
    ring_coupling_model()
    numerical_exploration()
    separate_sector_model()
    
    # Final summary
    print("â" * 70)
    print("SUMMARY AND CONCLUSIONS")
    print("â" * 70)
    print()
    print("  1. The circlette ring topology provides a SPECIFIC structure")
    print("     for generation mixing: three partial transition operators")
    print("     Ï_Gâ, Ï_Gâ, Ï_GâGâ, each connecting exactly 2 of 3")
    print("     generations. The R1 constraint (no (1,1)) is what creates")
    print("     this partial connectivity â it is structurally essential.")
    print()
    print("  2. If both up-type and down-type mass matrices have exact")
    print("     Zâ symmetry, V_CKM = I (no mixing). CKM mixing requires")
    print("     Zâ BREAKING, which the ring topology naturally provides")
    print("     through asymmetric adjacency (Gâ neighbours W, Gâ")
    print("     neighbours Câ).")
    print()
    print("  3. The ring-distance model predicts a coupling hierarchy")
    print("     Îµââ > Îµââ > Îµââ, based on the hop distances d=3,4,7.")
    print("     The experimental CKM hierarchy is |V_us| > |V_cb| > |V_ub|.")
    print()
    print("  4. The transition operator structure (Ï_Gâ connects 1â2,")
    print("     Ï_Gâ connects 1â3, Ï_GâGâ connects 2â3) maps naturally")
    print("     onto the Wolfenstein parametrisation of the CKM matrix.")
    print()
    print("  5. CP violation requires complex phases in the off-diagonal")
    print("     couplings. The ring topology provides a natural source:")
    print("     the two paths from Gâ to Iâ (clockwise vs counterclockwise)")
    print("     have different lengths, creating a path-dependent phase.")
    print()
    print("  OPEN QUESTIONS:")
    print("  â¢ Can the coupling strengths be derived from first principles")
    print("    (not just fitted)?")
    print("  â¢ Does the ring-distance decay constant ÎŸ have a physical")
    print("    interpretation?")
    print("  â¢ Can the CP-violating phase ÎŽ be predicted from the ring")
    print("    geometry (path asymmetry)?")
    print()
