"""
Wilson fermions: the standard fix for fermion doubling.

The Wilson Hamiltonian adds a momentum-dependent mass term:
    H(k) = α₁ sin(kx) + α₂ sin(ky) + [m + r(2 - cos kx - cos ky)] β

The Wilson mass M(k) = r(2 - cos kx - cos ky):
    M(0,0) = 0      — physical fermion remains massless
    M(π,0) = 2r     — doubler gapped
    M(0,π) = 2r     — doubler gapped
    M(π,π) = 4r     — doubler gapped

This works but the Wilson term is added BY HAND — it has no
fundamental origin in standard lattice gauge theory.
"""

import numpy as np
from numpy import sin, cos, pi
from dirac_matrices import alpha1, alpha2, beta
from naive_square import scan_brillouin_zone


def H_Wilson(kx, ky, m=0.0, r=1.0):
    """Wilson Hamiltonian with momentum-dependent mass."""
    M_eff = m + r * (2 - cos(kx) - cos(ky))
    return alpha1 * sin(kx) + alpha2 * sin(ky) + M_eff * beta


if __name__ == "__main__":
    print("=" * 60)
    print("WILSON FERMIONS — STANDARD DOUBLER FIX")
    print("=" * 60)
    print()
    print("H(k) = α₁ sin(kx) + α₂ sin(ky) + M(k) β")
    print("M(k) = r(2 - cos kx - cos ky)")
    print()

    # Analytical Wilson masses
    print("Wilson mass at high-symmetry points (r=1):")
    for name, kx, ky in [("Γ", 0, 0), ("X", pi, 0), ("Y", 0, pi), ("M", pi, pi)]:
        M = 1.0 * (2 - cos(kx) - cos(ky))
        evals = np.sort(np.linalg.eigvalsh(H_Wilson(kx, ky, r=1.0)))
        E_min = np.min(np.abs(evals))
        status = "massless ✓" if M < 0.01 else f"GAPPED (M={M:.0f}r) ✓"
        print(f"  {name} ({kx/pi:.0f}π,{ky/pi:.0f}π): M = {M:.1f}r, |E|_min = {E_min:.4f} — {status}")
    print()

    # Full BZ scan
    print("Full Brillouin zone scan (r=1):")
    dirac_pts, _ = scan_brillouin_zone(lambda kx, ky: H_Wilson(kx, ky, r=1.0))
    print(f"  Dirac points: {len(dirac_pts)}")
    if len(dirac_pts) == 0:
        print("  All doublers removed ✓")
    print()

    # Scan for different r values
    print("Doubler count vs Wilson parameter r:")
    for r in [0.0, 0.1, 0.3, 0.5, 1.0, 2.0]:
        pts, _ = scan_brillouin_zone(lambda kx, ky, r=r: H_Wilson(kx, ky, r=r), N=200)
        # The r=0 case is the naive lattice
        n_physical = sum(1 for kxp, kyp in pts if abs(kxp) < 0.1 and abs(kyp) < 0.1)
        n_doublers = len(pts) - n_physical
        print(f"  r = {r:.1f}: {len(pts)} Dirac points ({n_physical} physical, {n_doublers} doublers)")

    print()
    print("Conclusion: Wilson mass removes all doublers for any r > 0.")
    print("But the Wilson term is an ad-hoc addition with no fundamental origin.")
