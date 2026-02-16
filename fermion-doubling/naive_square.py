"""
Fermion doubling on the naive square lattice.

The standard Dirac Hamiltonian on a square lattice:
    H(k) = α₁ sin(kx) + α₂ sin(ky) + m β

has zeros (massless Dirac cones) at:
    (0, 0)  — the physical fermion
    (π, 0)  — doubler
    (0, π)  — doubler
    (π, π)  — doubler

This is the Nielsen-Ninomiya fermion doubling problem.
"""

import numpy as np
from numpy import sin, cos, pi
from scipy.ndimage import label
from dirac_matrices import alpha1, alpha2, beta


def H_naive(kx, ky, m=0.0):
    """Naive Dirac Hamiltonian on the square lattice."""
    return alpha1 * sin(kx) + alpha2 * sin(ky) + m * beta


def scan_brillouin_zone(H_func, N=300, threshold=0.05):
    """
    Scan the Brillouin zone for Dirac points (zero-energy modes).
    
    Returns:
        dirac_points: list of (kx/π, ky/π) for each Dirac point
        E_min: 2D array of minimum |eigenvalue| at each k-point
    """
    kx_range = np.linspace(-pi, pi, N)
    ky_range = np.linspace(-pi, pi, N)

    E_min = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            evals = np.linalg.eigvalsh(H_func(kx_range[j], ky_range[i]))
            E_min[i, j] = np.min(np.abs(evals))

    zero_mask = E_min < threshold
    labeled_arr, n_clusters = label(zero_mask)

    dirac_points = []
    for c in range(1, n_clusters + 1):
        pts = np.argwhere(labeled_arr == c)
        centre = pts.mean(axis=0)
        kx_c = kx_range[int(centre[1])]
        ky_c = ky_range[int(centre[0])]
        dirac_points.append((kx_c / pi, ky_c / pi))

    return dirac_points, E_min


if __name__ == "__main__":
    print("=" * 60)
    print("NAIVE SQUARE LATTICE — FERMION DOUBLING")
    print("=" * 60)
    print()
    print("H(k) = α₁ sin(kx) + α₂ sin(ky)")
    print()

    # Analytical check at key points
    key_points = [
        ("Γ (0,0)",   0,  0),
        ("X (π,0)",  pi,  0),
        ("Y (0,π)",   0, pi),
        ("M (π,π)",  pi, pi),
    ]

    print("Spectrum at high-symmetry points (m=0):")
    for name, kx, ky in key_points:
        evals = np.sort(np.linalg.eigvalsh(H_naive(kx, ky)))
        n_zero = np.sum(np.abs(evals) < 1e-10)
        print(f"  {name}: E = {np.round(evals, 6)}, zero modes: {n_zero}")
    print()

    # Full BZ scan
    print("Full Brillouin zone scan:")
    dirac_pts, E_min = scan_brillouin_zone(H_naive)
    print(f"  Dirac points found: {len(dirac_pts)}")
    for i, (kxp, kyp) in enumerate(dirac_pts):
        label_str = "PHYSICAL" if abs(kxp) < 0.1 and abs(kyp) < 0.1 else "DOUBLER"
        print(f"    ({kxp:+.3f}π, {kyp:+.3f}π) — {label_str}")
    print()

    n_doublers = len(dirac_pts) - 1
    print(f"Result: {len(dirac_pts)} Dirac points = 1 physical + {n_doublers} doublers")
    print("The naive lattice Dirac equation has fermion doubling.")
