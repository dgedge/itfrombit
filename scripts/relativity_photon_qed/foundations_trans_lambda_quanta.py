#!/usr/bin/env python3
"""Foundational trans-Lambda quanta fork.

This script formalises the consequence of the T-R2 photon cutoff closure.

Assumptions being tested:
  S1. The substrate has one spatial scale a0 = hbar c / Lambda_QCD.
  S2. The cosmological-constant route needs no independent UV oscillator density
      above Lambda_QCD.
  S3. The SC Wilson/Maxwell photon is the IR cohomology/envelope, not a UV lattice.

Question:
  How can the framework represent Lorentz-invariant quanta with omega >> Lambda_QCD?

Result:
  An elementary single-particle trans-Lambda massless mode is impossible under S1+S2.
  A finer photon UV lattice supports the mode but reopens the CC.  Brillouin aliasing
  destroys monotonic E=|p|.  The only CC-compatible representation left is a
  collinear multi-occupancy/null-bundle state: total four-momentum remains exactly
  null while each constituent lies below the single-scale support.  That is not yet
  a photon theorem; it needs an irreducible "bundle-as-one event" interaction map.

exit 0 = the no-go/fork arithmetic is internally consistent and the surviving route
         is explicitly named.
"""

from __future__ import annotations

import math


HBARC_GEV_FM = 0.1973269804
LAMBDA_GEV = 0.332
A0_FM = HBARC_GEV_FM / LAMBDA_GEV
E_BZ_GEV = math.pi * LAMBDA_GEV


def uv_scale_for_support(E_GeV: float) -> tuple[float, float, float]:
    """Return (required spacing fm, UV scale GeV, CC inflation factor).

    Nyquist support is the weakest possible condition: E <= pi hbar c/a.
    The corresponding oscillator cutoff scale is hbar c/a = E/pi.
    """

    a_req = math.pi * HBARC_GEV_FM / E_GeV
    uv = HBARC_GEV_FM / a_req
    cc = (uv / LAMBDA_GEV) ** 4
    return a_req, uv, cc


def null_bundle_mass_fraction(angles: list[float], weights: list[float]) -> float:
    """Invariant mass / energy for massless constituents in the x-y plane."""

    px = sum(w * math.cos(th) for w, th in zip(weights, angles))
    py = sum(w * math.sin(th) for w, th in zip(weights, angles))
    e = sum(weights)
    m2 = max(0.0, e * e - px * px - py * py)
    return math.sqrt(m2) / e


print("[0] Canon single-scale input")
print(f"    Lambda_QCD = {LAMBDA_GEV:.3f} GeV")
print(f"    a0 = hbar c / Lambda_QCD = {A0_FM:.6f} fm")
print(f"    one-scale Brillouin support ceiling = pi Lambda_QCD = {E_BZ_GEV:.3f} GeV")

print("\n[1] Elementary trans-Lambda massless mode: NO-GO under S1")
for E, label in ((1.0, "1 GeV"), (10.0, "10 GeV"), (1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
    K = E / LAMBDA_GEV
    supported = E <= E_BZ_GEV
    print(f"    {label:>7}: E/Lambda={K:9.1f}, E/(pi Lambda)={E / E_BZ_GEV:9.1f}, supported={supported}")
assert 1.0e3 > E_BZ_GEV

print("\n[2] Finer elementary photon lattice: supports the mode but costs the CC")
for E, label in ((10.0, "10 GeV"), (1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
    a_req, uv, cc = uv_scale_for_support(E)
    print(
        f"    {label:>7}: a <= {a_req:.3e} fm, "
        f"Lambda_UV >= {uv:.3e} GeV = {uv / LAMBDA_GEV:.3e} Lambda_QCD, "
        f"vacuum-density factor >= {cc:.3e}"
    )
assert uv_scale_for_support(1.0e3)[2] > 1.0e11

print("\n[3] Brillouin aliasing cannot be the answer")
for E in (10.0, 1.0e3):
    K = E / LAMBDA_GEV
    folded = ((K + math.pi) % (2.0 * math.pi)) - math.pi
    print(f"    E={E:g} GeV: physical K={K:.3e}, folded K={folded:.3e}")
print("    Folding keeps momenta in the first Brillouin zone, so it cannot preserve monotonic E=|p|.")
assert abs(((1.0e3 / LAMBDA_GEV + math.pi) % (2.0 * math.pi)) - math.pi) < math.pi

print("\n[4] CC-compatible survivor: collinear null bundle")
for E, label in ((10.0, "10 GeV"), (1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
    n_lambda = math.ceil(E / LAMBDA_GEV)
    n_bz = math.ceil(E / E_BZ_GEV)
    print(
        f"    {label:>7}: needs at least N={n_bz} BZ-edge quanta, "
        f"or N={n_lambda} Lambda-scale quanta"
    )
assert math.ceil(1.0e3 / LAMBDA_GEV) == 3013

same_direction_mass = null_bundle_mass_fraction([0.0, 0.0, 0.0], [1.0, 2.0, 3.0])
spread_1mrad_mass = null_bundle_mass_fraction([-5e-4, 0.0, 5e-4], [1.0, 1.0, 1.0])
spread_1deg_mass = null_bundle_mass_fraction(
    [-math.radians(0.5), 0.0, math.radians(0.5)], [1.0, 1.0, 1.0]
)
print(f"    exactly collinear bundle: M/E = {same_direction_mass:.3e}")
print(f"    1 mrad full angular spread: M/E = {spread_1mrad_mass:.3e}")
print(f"    1 degree full angular spread: M/E = {spread_1deg_mass:.3e}")
assert same_direction_mass == 0.0
assert spread_1mrad_mass > 0.0

print("\n[5] Verdict")
print(
    "    Under the single-scale + CC premise, there is no elementary trans-Lambda\n"
    "    massless Bloch photon.  A finer elementary UV lattice is mutually exclusive\n"
    "    with the CC cutoff premise, and Brillouin folding cannot preserve Lorentz\n"
    "    momentum.  The only remaining CC-compatible representation is a locked\n"
    "    collinear null bundle: many sub-cutoff lightlike records with total\n"
    "    P^mu=(E,E n).  This preserves the Lorentz invariant P^2=0 without adding\n"
    "    UV oscillator density.\n\n"
    "    OPEN THEOREM: prove a bundle-as-one event map.  The bundle must be a\n"
    "    gauge-invariant irreducible Wilson/service-current object whose absorption,\n"
    "    emission, and pair-production amplitudes depend on the total P^mu, not on\n"
    "    N independent soft photons.  Without that theorem, the single-scale\n"
    "    framework has no satisfactory account of observed omega >> Lambda_QCD\n"
    "    photons."
)
print("exit 0")
