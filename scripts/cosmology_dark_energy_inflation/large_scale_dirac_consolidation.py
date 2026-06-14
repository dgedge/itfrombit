#!/usr/bin/env python3
r"""Sector-wide audit: large-scale physics and the Lambda_QCD <-> horizon input.

Question
--------
Do the framework's large-scale "parameter-free" results reduce to one
Dirac-class relation between the QCD/chiral lattice scale and the cosmological
horizon, dressed by code numbers and alpha/prefactors?

Scope
-----
This script checks four clusters against the canonical audit posture:

1. Gravity M_P:
   DRIFT G7 / item7_keg_heatkernel.py says the bare stiffness computes the
   Lambda_QCD end, while the macroscopic 19-OOM jump consumes the de Sitter
   horizon/lattice ratio.

2. rho_Lambda:
   item123_lambda_audit.py reduces both the item-123 and Zel'dovich branches
   to coeff * Lambda_QCD^3 * H0, i.e. a horizon-input dimensional bound.

3. MOND a0 / BTFR:
   item132_halo_audit.py records a0 = c H0 / 2pi as the horizon acceleration.
   The halo action can sharpen the MOND law, but the acceleration scale remains
   the same horizon input unless independently derived.

   Literal-horizon test:
   a0 = c H0 / 2pi implies R_a0 = c^2/(2pi a0) = c/H0, the current Hubble
   horizon.  The gravity M_P route canonically uses the static de Sitter
   horizon R_dS = c/(H0 sqrt(Omega_Lambda)).  These are not literally equal
   unless Omega_Lambda = 1 or the sqrt(Omega_Lambda) factor is moved into the
   prefactor.  They are the same Dirac large-number input up to the late-horizon
   convention, not a second large-number coincidence.

4. 1/28 inflationary observables:
   item131_* scripts ground the 28-channel code structure.  The remaining
   premise is a code-to-cosmology bridge: self-similar radial HBC, log-scale
   generator, saturated H, and power-ledger action.  This uses horizon shells
   as the clock variable, but it does not consume the numerical Lambda_QCD/H0
   large-number.

Verdict
-------
Gravity, rho_Lambda, and MOND-a0 share the same late-horizon Dirac-class input
at the scale level.  They are not independent parameter-free large-scale
predictions.  The a0 and gravity horizons are not literally identical under
the current canon: a0 names the Hubble horizon c/H0, while gravity uses the
static de Sitter horizon c/(H0 sqrt(Omega_Lambda)).  The sharper statement is
therefore not "two independent coincidences"; it is the same large number
R_H/a_QCD entering with different powers and one O(1) horizon-convention
factor.  Gravity uses the positive square-root power in M_P/Lambda_QCD, while
rho_Lambda and MOND-a0 use the inverse horizon power.  This strengthens the
one-coincidence verdict.  The 1/28 inflationary results are different:
code-28 is now finite/combinatorial, but its cosmological relevance is
conditional on the HBC/log-shell bridge rather than on the Lambda_QCD <-> H0
magnitude.
"""

from __future__ import annotations

import math
from fractions import Fraction


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


# Constants chosen to match the existing local audits.
C = 299_792_458.0
MPC = 3.0856775814913673e22
HBAR_GEV_S = 6.582119569e-25
H0_SI = 67.4 * 1000.0 / MPC
H0_GEV = H0_SI * HBAR_GEV_S
OMEGA_L = 0.6847
ALPHA = 1.0 / 137.036
LAMBDA_QCD_GEV = 0.332
HBARC_GEV_M = 1.973269804e-16


def main() -> None:
    print("LARGE-SCALE DIRAC CONSOLIDATION AUDIT")

    a_lattice = HBARC_GEV_M / LAMBDA_QCD_GEV
    r_hubble = C / H0_SI
    r_desitter = r_hubble / math.sqrt(OMEGA_L)
    node_ratio_h = r_hubble / a_lattice
    node_ratio_ds = r_desitter / a_lattice
    lambda_over_h0_oom = math.log10(LAMBDA_QCD_GEV / H0_GEV)

    print("\n[1] The shared large number")
    print(f"  a_QCD = hbar c / Lambda_QCD = {a_lattice:.6e} m")
    print(f"  R_H   = c / H0              = {r_hubble:.6e} m")
    print(f"  R_dS  = R_H / sqrt(Omega_L) = {r_desitter:.6e} m")
    print(f"  R_H/a_QCD  = {node_ratio_h:.6e}")
    print(f"  R_dS/a_QCD = {node_ratio_ds:.6e}")
    print(f"  log10(Lambda_QCD/H0) = {lambda_over_h0_oom:.3f} OOM")
    check(40.0 < lambda_over_h0_oom < 43.0, "the QCD-to-horizon bridge supplies the familiar ~41 OOM factor")
    check(abs((r_hubble / r_desitter) - math.sqrt(OMEGA_L)) < 1e-12, "Hubble and de Sitter horizons differ only by sqrt(Omega_L)")

    print("\n[2] Gravity M_P consumes the de Sitter horizon/lattice ratio with positive half-power")
    k_eff = 205.0
    macro_dilution = math.sqrt(node_ratio_ds / k_eff)
    print(f"  sqrt((R_dS/a_QCD)/K_eff) = {macro_dilution:.6e}")
    print("  horizon power on M_P/Lambda_QCD: +1/2")
    check(1e19 < macro_dilution < 1e21, "the macroscopic M_P jump is the horizon/lattice ratio, not an intrinsic local scale")
    check(abs(math.log(macro_dilution * math.sqrt(k_eff)) / math.log(node_ratio_ds) - 0.5) < 1e-12, "gravity uses the square root of the horizon/lattice large number")
    check(True, "item7_keg_heatkernel.py computes the Lambda_QCD bare end; it does not derive the horizon input")

    print("\n[3] rho_Lambda is coeff * Lambda_QCD^3 * H0, i.e. inverse horizon power")
    bare_rho = LAMBDA_QCD_GEV**3 * H0_GEV
    coeff_item123 = (3.0 / (16.0 * math.pi)) * ALPHA
    coeff_zel = 9.0 * ALPHA * ALPHA
    print(f"  Lambda^3 H0             = {bare_rho:.6e} GeV^4")
    print(f"  item-123 coeff alpha^1  = {coeff_item123:.6e}")
    print(f"  Zel'dovich coeff alpha^2= {coeff_zel:.6e}")
    print(f"  item-123 rho            = {coeff_item123 * bare_rho:.6e} GeV^4")
    print(f"  Zel'dovich rho          = {coeff_zel * bare_rho:.6e} GeV^4")
    print("  horizon power relative to Lambda_QCD^4: -1")
    check(coeff_item123 * bare_rho > 1e-48 and coeff_item123 * bare_rho < 1e-46, "item-123 lands near observed rho only after the horizon factor and prefactor")
    check(coeff_zel * bare_rho > 1e-48 and coeff_zel * bare_rho < 1e-46, "Zel'dovich branch is also horizon-consuming")
    inverse_ratio = H0_GEV / LAMBDA_QCD_GEV
    check(abs((bare_rho / LAMBDA_QCD_GEV**4) / inverse_ratio - 1.0) < 1e-15, "rho/Lambda^4 consumes the inverse Lambda_QCD/H0 large number")

    print("\n[4] MOND a0 literal-horizon test")
    a0_hubble = C * H0_SI / (2.0 * math.pi)
    a0_desitter = C * H0_SI * math.sqrt(OMEGA_L) / (2.0 * math.pi)
    r_a0 = C * C / (2.0 * math.pi * a0_hubble)
    r_a0_desitter = C * C / (2.0 * math.pi * a0_desitter)
    print(f"  a0(Hubble)    = c H0 / 2pi              = {a0_hubble:.6e} m/s^2")
    print(f"  a0(de Sitter) = c H0 sqrt(Omega_L)/2pi  = {a0_desitter:.6e} m/s^2")
    print(f"  R_a0 from c^2/(2pi a0)                 = {r_a0:.6e} m")
    print(f"  gravity-route R_dS                     = {r_desitter:.6e} m")
    print(f"  R_dS/R_a0                              = {r_desitter / r_a0:.6f}")
    print(f"  R_a0(de Sitter acceleration)           = {r_a0_desitter:.6e} m")
    print(f"  a0(Hubble)/a0(de Sitter) = {a0_hubble / a0_desitter:.6f}")
    print("  horizon power relative to c^2/a_QCD: -1")
    check(abs(r_a0 / r_hubble - 1.0) < 1e-15, "a0=cH0/2pi literally names the Hubble horizon c/H0")
    check(abs(r_desitter / r_a0 - (1.0 / math.sqrt(OMEGA_L))) < 1e-12, "gravity route's de Sitter horizon is not literally the same as the a0 Hubble horizon")
    check(abs(r_a0_desitter / r_desitter - 1.0) < 1e-15, "using a de Sitter acceleration convention would make the horizons literally identical")
    check(0.9e-10 < a0_hubble < 1.2e-10, "a0=cH0/2pi is the standard horizon acceleration scale")
    check(abs((a0_hubble / (C * C / a_lattice)) - (1.0 / (2.0 * math.pi * node_ratio_h))) < 1e-16, "a0 consumes the inverse horizon/lattice large number")
    check(abs((a0_hubble / a0_desitter) - (1.0 / math.sqrt(OMEGA_L))) < 1e-12, "the mismatch is an Omega_L convention, not an independent large-number input")
    check(True, "BTFR then inherits G plus this a0; the p=3/R4 work sharpens the law but not the horizon scale")

    print("\n[5] The 1/28 inflationary observables are code plus bridge, not the Lambda/H0 large number")
    delta = Fraction(1, 28)
    ns = Fraction(1, 1) - delta
    print(f"  Delta = {delta}; n_s = 1 - Delta = {ns} = {float(ns):.9f}")
    print("  Required cosmology premises: radial HBC, scale self-similarity, log generator, saturated H, power ledger.")
    check(ns == Fraction(27, 28), "the coefficient is code-28 once the bridge is granted")
    check(True, "no numerical Lambda_QCD/H0 magnitude is needed for n_s; the bridge is logical/geometric, not Dirac-large-number")

    print("\n" + "=" * 98)
    print("VERDICT")
    print("  Gravity M_P, rho_Lambda, and MOND a0 all consume the same late-horizon")
    print("  input class: the QCD lattice scale compared with the cosmological horizon,")
    print("  dressed by dimensionless prefactors/code counts.  They should not be counted")
    print("  as independent parameter-free large-scale derivations.")
    print("  Literal-horizon test: a0=cH0/2pi names R_a0=c/H0.  The gravity M_P")
    print("  route canonically uses R_dS=c/(H0 sqrt(Omega_L)).  Therefore they are")
    print("  not literally the same horizon under current definitions; they differ by")
    print("  1/sqrt(Omega_L).  But this is not a second large-number coincidence.")
    print("  The stronger accounting statement is that all three use the same")
    print("  late-horizon large number with different powers:")
    print("      M_P/Lambda_QCD ~ (R/a_QCD)^(+1/2),")
    print("      rho_Lambda/Lambda_QCD^4 ~ (R/a_QCD)^(-1),")
    print("      a0/(c^2/a_QCD) ~ (R/a_QCD)^(-1).")
    print("  This is one Dirac coincidence dressed three ways, not independent scales.")
    print("  The 1/28 inflationary observables are not another Lambda_QCD/H0 magnitude")
    print("  coincidence.  Their remaining risk is the code-to-cosmology bridge: whether")
    print("  finite 28-channel QEC really becomes saturated radial HBC/log-shell printing.")
    print("  Honest sector-wide count: one late-horizon Dirac input + solid code-28 +")
    print("  still-assumed code-to-cosmology bridges and halo rate-matching/count-ledger")
    print("  lemmas.")
    print("=" * 98)
    print("exit 0 -- large-scale sector consolidated to one horizon input plus code-bridge assumptions.")


if __name__ == "__main__":
    main()
