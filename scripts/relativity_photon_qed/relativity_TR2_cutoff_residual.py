#!/usr/bin/env python3
"""T-R2 photon cutoff residual.

The object-identity theorem establishes the physical photon as the
Gauss-projected Wilson/Maxwell cohomology, not the raw K6 T1u/Eg branch.
This script checks the separate question that remains: can the Wilson
cochain spacing be literally the QCD matter-cell spacing a0 ~= 1/Lambda?

It cannot, if the same lattice is meant as the photon's UV support.  The
Wilson/Maxwell anisotropy is harmless at optical wavelengths, but a femtometre
lattice has a Brillouin support scale of O(1 GeV), not TeV.  Therefore the
closed object identity leaves a named UV bridge:

    either a_photon << a_QCD is derived as a refined/renormalised gauge web,
    or the Wilson lattice is only an IR effective cohomology, not the literal
    microscopic photon cutoff.
"""

from __future__ import annotations

import math


HBARC_GEV_FM = 0.1973269804
LAMBDA_QCD_GEV = 0.332
A_QCD_FM = HBARC_GEV_FM / LAMBDA_QCD_GEV


def wilson_group_velocity(k_fm_inv: float, direction: tuple[float, float, float], a_fm: float) -> float:
    """d omega / d |p| for SC Wilson/Maxwell dispersion in a given direction.

    Dimensionless lattice momentum is K = a |p|.  For direction n,
    omega(K) = 2 sqrt(sum_i sin(K n_i / 2)^2).  This finite-difference
    derivative is enough for the small-k coefficient gate below.
    """

    norm = math.sqrt(sum(x * x for x in direction))
    n = tuple(x / norm for x in direction)
    K = a_fm * k_fm_inv

    def omega(K_local: float) -> float:
        return 2.0 * math.sqrt(sum(math.sin(K_local * ni / 2.0) ** 2 for ni in n))

    h = max(1e-7, 1e-5 * max(1.0, abs(K)))
    return (omega(K + h) - omega(K - h)) / (2.0 * h)


def anisotropy_group(K: float) -> float:
    """Relative [100]-[111] group-velocity anisotropy at dimensionless |K|."""

    a = 1.0
    k = K
    v100 = wilson_group_velocity(k, (1.0, 0.0, 0.0), a)
    v111 = wilson_group_velocity(k, (1.0, 1.0, 1.0), a)
    return abs(v111 - v100) / (0.5 * (v111 + v100))


def a_bound_from_anisotropy(E_GeV: float, bound: float, coeff: float = 1.0 / 12.0) -> float:
    """Maximum lattice spacing in fm from coeff * (a E/hbarc)^2 <= bound."""

    return math.sqrt(bound / coeff) * HBARC_GEV_FM / E_GeV


def a_bound_from_nyquist(E_GeV: float) -> float:
    """Maximum lattice spacing in fm from E < pi hbar c / a."""

    return math.pi * HBARC_GEV_FM / E_GeV


def main() -> None:
    print("[1] SC Wilson/Maxwell anisotropy coefficient")
    for K in (2e-2, 1e-2, 5e-3):
        measured = anisotropy_group(K)
        expected = K * K / 12.0
        print(
            f"    K={K:.3e}: Delta v/v={measured:.6e}, "
            f"K^2/12={expected:.6e}, ratio={measured / expected:.6f}"
        )
        assert abs(measured / expected - 1.0) < 5e-4

    print("\n[2] Optical SME-scale check for a_QCD")
    # 500 nm photon: E = hbar c * 2 pi / lambda.
    lambda_visible_m = 500e-9
    E_visible_GeV = (2.0 * math.pi * 0.1973269804e-15) / lambda_visible_m
    K_visible = A_QCD_FM * E_visible_GeV / HBARC_GEV_FM
    dv_visible = K_visible * K_visible / 12.0
    a_sme_fm = a_bound_from_anisotropy(E_visible_GeV, 1e-17)
    print(f"    a_QCD = hbar c / Lambda_QCD = {A_QCD_FM:.6f} fm")
    print(f"    500 nm: E={E_visible_GeV:.3e} GeV, a_QCD*k={K_visible:.3e}")
    print(f"    Wilson artifact = (a k)^2/12 = {dv_visible:.3e}")
    print(f"    optical 1e-17 anisotropy bound would allow a <= {a_sme_fm:.3f} fm")
    assert dv_visible < 1e-17
    assert A_QCD_FM < a_sme_fm

    print("\n[3] UV support check")
    for E_GeV, label in ((1.0, "1 GeV"), (1e3, "1 TeV"), (1e5, "100 TeV")):
        K_qcd = A_QCD_FM * E_GeV / HBARC_GEV_FM
        a_nyq = a_bound_from_nyquist(E_GeV)
        ratio = A_QCD_FM / a_nyq
        print(
            f"    {label:>7}: a_QCD*k={K_qcd:.3e}, "
            f"a_Nyquist<={a_nyq:.3e} fm, a_QCD/a_Nyquist={ratio:.3e}"
        )
    assert A_QCD_FM * 1e3 / HBARC_GEV_FM > math.pi

    print("\n[4] THEOREM STATUS")
    print(
        "    Object identity: CLOSED under Wilson/Gauss premises.\n"
        "    Cutoff identity: OPEN.  If a_photon = a_QCD, the photon Brillouin\n"
        "    scale is only about pi*Lambda_QCD = 1.04 GeV, so TeV photons are not\n"
        "    supportable as lattice modes.  Therefore T-R2 now requires a separate\n"
        "    UV bridge: derive a_photon << a_QCD, or state that the Wilson lattice\n"
        "    is an IR cohomology/effective action whose cutoff is not the matter-cell\n"
        "    spacing.  Without that bridge, the femtometre-cutoff photon reading is\n"
        "    falsified by high-energy photon support."
    )
    print("exit 0")


if __name__ == "__main__":
    main()
