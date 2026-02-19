"""
Comparative band structure plots for the fermion doubling analysis.

Generates:
    1. Band structures along Γ→X→M→Γ for all configurations
    2. Brillouin zone dispersion surface maps
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dirac_matrices import alpha1, alpha2, alpha12, beta
from naive_square import H_naive
from wilson_fermions import H_Wilson
from lattice_488 import H_488_physical


def compute_bands(H_func, n_points=400):
    """Compute band structure along Γ→X→M→Γ path."""
    segments = [
        ((0, 0), (pi, 0)),      # Γ → X
        ((pi, 0), (pi, pi)),    # X → M
        ((pi, pi), (0, 0)),     # M → Γ
    ]

    all_k, all_E = [], []
    k_ticks = [0]
    k_pos = 0.0

    for (k1x, k1y), (k2x, k2y) in segments:
        for t in np.linspace(0, 1, n_points, endpoint=False):
            kx = k1x + t * (k2x - k1x)
            ky = k1y + t * (k2y - k1y)
            evals = np.sort(np.linalg.eigvalsh(H_func(kx, ky)))
            all_k.append(k_pos)
            all_E.append(evals)
            dk = sqrt((k2x - k1x)**2 + (k2y - k1y)**2) / n_points
            k_pos += dk
        k_ticks.append(k_pos)

    return np.array(all_k), np.array(all_E), k_ticks


def plot_band_comparison(output_path="bands_comparison.png"):
    """Plot band structures for all four configurations."""
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

    configs = [
        ("Naive square lattice\n(1 physical + 3 doublers)",
         lambda kx, ky: H_naive(kx, ky)),
        ("Wilson fermions (r=1)\n(doublers gapped by hand)",
         lambda kx, ky: H_Wilson(kx, ky, r=1.0)),
        ("4.8.8 lattice, t₂ = 0.25\n(physical NNN: α₁α₂ × sin·sin)",
         lambda kx, ky: H_488_physical(kx, ky, t2=0.25)),
        ("4.8.8 lattice, t₂ = 0.5\n(physical NNN: α₁α₂ × sin·sin)",
         lambda kx, ky: H_488_physical(kx, ky, t2=0.5)),
    ]

    for ax, (title, H_func) in zip(axes.flat, configs):
        k_arr, E_arr, k_ticks = compute_bands(H_func)

        for b in range(E_arr.shape[1]):
            touches_zero = np.min(np.abs(E_arr[:, b])) < 0.15
            color = '#E8553D' if touches_zero else '#2B6B8E'
            lw = 2.5 if touches_zero else 1.5
            ax.plot(k_arr, E_arr[:, b], color=color, linewidth=lw, alpha=0.8)

        ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
        for kt in k_ticks:
            ax.axvline(kt, color='gray', linewidth=0.5, linestyle='--')

        ax.set_xticks(k_ticks)
        ax.set_xticklabels(['Γ', 'X', 'M', 'Γ'], fontsize=16)
        ax.set_ylabel('E', fontsize=16)
        ax.set_title(title, fontsize=15, fontweight='bold')
        ax.set_ylim(-4, 4)
        ax.tick_params(labelsize=13)

    plt.suptitle("Fermion Doubling: Band Structure Comparison",
                 fontsize=20, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"Band comparison saved to {output_path}")


def plot_dispersion_surfaces(output_path="dispersion_surface.png"):
    """Plot |E|_min across the Brillouin zone for three configurations."""
    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))

    configs = [
        ("Naive (doublers at X, Y, M)", lambda kx, ky: H_naive(kx, ky)),
        ("4.8.8 t₂ = 0.25",            lambda kx, ky: H_488_physical(kx, ky, t2=0.25)),
        ("Wilson r = 1 (fixed)",         lambda kx, ky: H_Wilson(kx, ky, r=1.0)),
    ]

    N = 150
    kx_arr = np.linspace(-pi, pi, N)
    ky_arr = np.linspace(-pi, pi, N)
    KX, KY = np.meshgrid(kx_arr, ky_arr)

    for ax, (title, H_func) in zip(axes, configs):
        E_min = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                evals = np.linalg.eigvalsh(H_func(kx_arr[j], ky_arr[i]))
                E_min[i, j] = np.min(np.abs(evals))

        im = ax.pcolormesh(KX / pi, KY / pi, E_min, cmap='hot_r',
                           vmin=0, vmax=2, shading='auto')
        ax.set_xlabel('kx / π', fontsize=14)
        ax.set_ylabel('ky / π', fontsize=14)
        ax.set_title(title, fontsize=15, fontweight='bold')
        ax.set_aspect('equal')
        ax.tick_params(labelsize=12)
        plt.colorbar(im, ax=ax, label='|E|_min', shrink=0.8)

        # Mark high-symmetry points
        for (kxp, kyp) in [(0, 0), (1, 0), (0, 1), (1, 1),
                           (-1, 0), (0, -1), (-1, -1), (-1, 1), (1, -1)]:
            ax.plot(kxp, kyp, 'wo', markersize=5,
                    markeredgecolor='black', markeredgewidth=1)

    plt.suptitle("Minimum Energy |E| in Brillouin Zone (dark = Dirac point)",
                 fontsize=17, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"Dispersion surfaces saved to {output_path}")


if __name__ == "__main__":
    plot_band_comparison()
    plot_dispersion_surfaces()
    print("\nAll plots generated.")
