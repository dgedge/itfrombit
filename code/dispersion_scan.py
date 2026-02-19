"""
Full Brillouin zone scanning for Dirac point detection.

Provides systematic scanning of the BZ with configurable resolution
and threshold, producing heatmaps and Dirac point catalogues.
"""

import numpy as np
from numpy import pi
from scipy.ndimage import label
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from dirac_matrices import alpha1, alpha2, alpha12, beta
from naive_square import H_naive
from wilson_fermions import H_Wilson
from lattice_488 import H_488_physical


def full_bz_scan(H_func, N=300, threshold=0.05):
    """
    Scan the full Brillouin zone [-π,π]² and catalogue all Dirac points.

    Returns:
        results: dict with keys 'E_min', 'dirac_points', 'n_clusters'
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
        e_at_centre = E_min[int(centre[0]), int(centre[1])]
        dirac_points.append({
            'kx': kx_c, 'ky': ky_c,
            'kx_pi': kx_c / pi, 'ky_pi': ky_c / pi,
            'E_min': e_at_centre,
            'n_pixels': len(pts),
        })

    return {
        'E_min': E_min,
        'dirac_points': dirac_points,
        'n_clusters': n_clusters,
        'kx_range': kx_range,
        'ky_range': ky_range,
    }


def classify_point(kxp, kyp, tol=0.15):
    """Classify a Dirac point as physical or doubler."""
    if abs(kxp) < tol and abs(kyp) < tol:
        return "Γ (physical)"
    elif abs(abs(kxp) - 1) < tol and abs(kyp) < tol:
        return "X (doubler)"
    elif abs(kxp) < tol and abs(abs(kyp) - 1) < tol:
        return "Y (doubler)"
    elif abs(abs(kxp) - 1) < tol and abs(abs(kyp) - 1) < tol:
        return "M (doubler)"
    else:
        return "other"


if __name__ == "__main__":
    print("=" * 70)
    print("SYSTEMATIC BRILLOUIN ZONE SCAN")
    print("=" * 70)
    print()

    configs = [
        ("Naive square lattice",      lambda kx, ky: H_naive(kx, ky)),
        ("Wilson r=1 (by hand)",       lambda kx, ky: H_Wilson(kx, ky, r=1.0)),
        ("4.8.8 t₂=0 (NN only)",      lambda kx, ky: H_488_physical(kx, ky, t2=0)),
        ("4.8.8 t₂=0.25",             lambda kx, ky: H_488_physical(kx, ky, t2=0.25)),
        ("4.8.8 t₂=0.5",              lambda kx, ky: H_488_physical(kx, ky, t2=0.5)),
        ("4.8.8 t₂=1.0",              lambda kx, ky: H_488_physical(kx, ky, t2=1.0)),
    ]

    print(f"{'Configuration':<30} {'Points':>7} {'Physical':>9} {'Doublers':>9}")
    print("-" * 60)

    for name, H_func in configs:
        res = full_bz_scan(H_func, N=200)
        n_phys = sum(1 for p in res['dirac_points']
                     if "physical" in classify_point(p['kx_pi'], p['ky_pi']))
        n_doub = sum(1 for p in res['dirac_points']
                     if "doubler" in classify_point(p['kx_pi'], p['ky_pi']))
        print(f"  {name:<28} {res['n_clusters']:>7} {n_phys:>9} {n_doub:>9}")

    print()
    print("Observation: The physical NNN coupling from the 4.8.8 topology")
    print("does not reduce the doubler count for any value of t₂.")
    print("Only the Wilson mass (added by hand) removes doublers.")
