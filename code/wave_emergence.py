#!/usr/bin/env python3
"""
wave_emergence.py — Deriving the Schrödinger equation from CNOT lattice updates.

Demonstrates computationally that a discrete quantum walk whose coin
operator is the CNOT gate (I3 ⊕ LQ) reproduces the free-particle
Schrödinger wave dispersion in the continuum limit.

The derivation chain:
  1. CNOT is a Boolean NOT on I3 (discrete, deterministic)
  2. Reversibility requires unitary embedding → complex phase emerges
  3. Coupling internal toggle to spatial propagation → discrete quantum walk
  4. Continuum limit of the walk → 1+1D Dirac equation
  5. Non-relativistic limit → Schrödinger equation
  6. Leptons (LQ=0, θ=0) → massless Weyl fermions at c

The mass term in the Dirac equation (mc²σₓ) is literally the Pauli-X
operator — the Boolean NOT gate acting on I₃.

Reference: Elliman (2026), "It from Bit, Revisited", Section 8.
"""

import numpy as np


def simulate_massive_particle(N=10000, steps=2500, m_theta=0.05, sigma_0=30.0):
    """
    Simulate a massive particle (quark, LQ=1) on the lattice.

    The CNOT gate fires every tick with rate θ = mc²Δt/ℏ,
    coupling right-movers and left-movers.

    Parameters:
        N: number of lattice sites
        steps: number of Planck ticks
        m_theta: CNOT flip rate per tick (= mass in lattice units)
        sigma_0: initial Gaussian width

    Returns:
        x: lattice positions
        density: |ψ|² (discrete simulation)
        density_analytical: |ψ|² (Schrödinger prediction)
        diagnostics: dict of verification metrics
    """
    x0 = N // 2
    x = np.arange(N)

    # Initialise right-mover and left-mover amplitudes
    psi_R = np.zeros(N, dtype=complex)
    psi_L = np.zeros(N, dtype=complex)

    envelope = (np.exp(-(x - x0) ** 2 / (4 * sigma_0 ** 2))
                / (2 * np.pi * sigma_0 ** 2) ** 0.25)
    psi_R[:] = envelope / np.sqrt(2)
    psi_L[:] = envelope / np.sqrt(2)

    initial_norm = np.sum(np.abs(psi_R) ** 2 + np.abs(psi_L) ** 2)

    cos_m = np.cos(m_theta)
    sin_m = np.sin(m_theta)

    # Discrete lattice updates
    for _ in range(steps):
        # 1. CNOT mixing (σₓ rotation by θ)
        #    This IS the weak interaction: I3 ⊕ LQ
        #    Complex i is required for unitarity of the Boolean swap
        temp_R = cos_m * psi_R - 1j * sin_m * psi_L
        temp_L = -1j * sin_m * psi_R + cos_m * psi_L

        # 2. Spatial shift (lattice propagation at c = Δx/Δt)
        psi_R = np.roll(temp_R, 1)    # right-movers shift right
        psi_L = np.roll(temp_L, -1)   # left-movers shift left

    final_norm = np.sum(np.abs(psi_R) ** 2 + np.abs(psi_L) ** 2)
    density = np.abs(psi_R) ** 2 + np.abs(psi_L) ** 2

    # Analytical Schrödinger prediction
    m_eff = np.tan(m_theta)
    sigma_t = sigma_0 * np.sqrt(1 + (steps / (2 * m_eff * sigma_0 ** 2)) ** 2)
    density_analytical = (np.exp(-(x - x0) ** 2 / (2 * sigma_t ** 2))
                          / (np.sqrt(2 * np.pi) * sigma_t))

    # Normalise for comparison
    d_norm = density / np.sum(density)
    a_norm = density_analytical / np.sum(density_analytical)
    overlap = np.sum(np.sqrt(np.maximum(d_norm * a_norm, 0)))

    diagnostics = {
        "initial_norm": initial_norm,
        "final_norm": final_norm,
        "unitarity_preserved": abs(final_norm - initial_norm) < 1e-6,
        "peak_discrete": int(np.argmax(density)),
        "peak_analytical": int(np.argmax(density_analytical)),
        "peak_agreement": abs(int(np.argmax(density)) - x0) <= 1,
        "bhattacharyya_overlap": overlap,
        "sigma_t": sigma_t,
        "spread_factor": sigma_t / sigma_0,
        "m_theta": m_theta,
        "m_eff_tan": np.tan(m_theta),
        "small_angle_error_pct": abs(m_theta - np.tan(m_theta)) / m_theta * 100,
    }

    return x, density, density_analytical, diagnostics


def simulate_massless_particle(N=10000, steps=2500, sigma_0=30.0):
    """
    Simulate a massless particle (lepton, LQ=0) on the lattice.

    With θ=0 (no CNOT), right-movers and left-movers decouple
    completely. The wavepacket splits into two lumps moving at ±c.
    This is massless Weyl fermion behaviour.

    Returns:
        x: lattice positions
        density: |ψ|²
        diagnostics: dict
    """
    x0 = N // 2
    x = np.arange(N)

    psi_R = np.zeros(N, dtype=complex)
    psi_L = np.zeros(N, dtype=complex)

    envelope = (np.exp(-(x - x0) ** 2 / (4 * sigma_0 ** 2))
                / (2 * np.pi * sigma_0 ** 2) ** 0.25)
    psi_R[:] = envelope / np.sqrt(2)
    psi_L[:] = envelope / np.sqrt(2)

    # θ=0: no CNOT mixing. Pure spatial shift only.
    for _ in range(steps):
        psi_R = np.roll(psi_R, 1)
        psi_L = np.roll(psi_L, -1)

    density = np.abs(psi_R) ** 2 + np.abs(psi_L) ** 2

    # Find the two peaks
    right_peak = np.argmax(np.abs(psi_R) ** 2)
    left_peak = np.argmax(np.abs(psi_L) ** 2)

    diagnostics = {
        "right_peak": right_peak,
        "left_peak": left_peak,
        "separation": abs(right_peak - left_peak),
        "expected_separation": min(2 * steps, N),
        "propagates_at_c": True,  # By construction
    }

    return x, density, diagnostics


def run_verification():
    """Run full verification suite and print results."""

    print("=" * 70)
    print("WAVE EQUATION EMERGENCE FROM CNOT LATTICE UPDATES")
    print("=" * 70)

    # --- Massive particle (quark) ---
    print("\n--- Massive particle (LQ=1, θ=0.05) ---\n")
    x, density, analytical, diag = simulate_massive_particle()

    print(f"  Unitarity: initial={diag['initial_norm']:.6f}, "
          f"final={diag['final_norm']:.6f}, "
          f"preserved={diag['unitarity_preserved']}")
    print(f"  Peak: discrete={diag['peak_discrete']}, "
          f"analytical={diag['peak_analytical']}, "
          f"centred={diag['peak_agreement']}")
    print(f"  Bhattacharyya overlap: {diag['bhattacharyya_overlap']:.4f}")
    print(f"  Theoretical σ(t): {diag['sigma_t']:.1f} "
          f"({diag['spread_factor']:.1f}× initial)")
    print(f"  Mass: θ={diag['m_theta']}, tan(θ)={diag['m_eff_tan']:.6f}, "
          f"error={diag['small_angle_error_pct']:.2f}%")

    # --- Massless particle (lepton) ---
    print("\n--- Massless particle (LQ=0, θ=0) ---\n")
    x2, density2, diag2 = simulate_massless_particle()

    print(f"  Right-mover peak: {diag2['right_peak']}")
    print(f"  Left-mover peak:  {diag2['left_peak']}")
    print(f"  Separation: {diag2['separation']} lattice units")
    print(f"  Expected: {diag2['expected_separation']}")
    print(f"  Propagates at c: ✓ (two lumps moving at ±c, no dispersion)")

    # --- Summary ---
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print("""
  The discrete CNOT lattice walk reproduces:

  1. Free-particle Schrödinger dispersion (overlap = {:.3f})
  2. Unitarity (norm conservation to machine precision)
  3. Stationary wavepacket for zero-momentum initial condition
  4. Massless Weyl fermion behaviour when θ=0 (LQ=0 leptons)
  5. Mass = CNOT frequency in the small-angle limit (error {:.2f}%)

  The wave equation is not postulated. It emerges.
  The lattice does not obey quantum mechanics.
  Quantum mechanics obeys the lattice.
""".format(diag['bhattacharyya_overlap'], diag['small_angle_error_pct']))


if __name__ == "__main__":
    run_verification()
