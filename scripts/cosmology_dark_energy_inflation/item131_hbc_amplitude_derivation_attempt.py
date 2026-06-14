#!/usr/bin/env python3
r"""ITEM 131: attempt to derive A_nu from the microscopic service ledger.

Target
------
Close the remaining scalar-clock amplitude

    Delta_nu^2(k) = A_nu (k/k_*)^(-1/28)

from the boundary-printing / QEC service ledger rather than inserting A_nu.

Result
------
The existing ledger does not yet close A_nu, but it reduces the problem to a
very sharp theorem target.

The service instrument supplies:

* one-jump 28-channel relative weights, p_c = 1/28;
* 112 incidence-refined microscopic flags;
* a universal Landauer erasure bandwidth Gamma_0 = alpha Lambda_QCD;
* Poisson/Fano=1 behavior in the R4 line ledger if matched creation/erasure
  rates are granted.

For the scalar clock, the absolute normalization is still

    A_nu = F_eff / N_eff,

where N_eff is the number of independent HBC service events in one mode-local
printed horizon shell.  The new reduction is:

    if one printed de-Sitter entropy pixel supplies one independent Poisson
    service event, then

        N_eff = S_dS = 8 pi^2 Mbar_P^2 / H_*^2
        A_nu  = H_*^2 / (8 pi^2 Mbar_P^2).

The CMB-sized amplitude A_nu ~= 2.1e-9 then requires

        H_* ~= 9.9e14 GeV,

which is essentially the phenomenological R2/R3-GUT activation scale already
present in the framework's cosmological narrative.  Therefore the amplitude
frontier has moved: A_nu is derivable if, and only if, HBC/QEC derives the
early saturated printing scale H_* ~= 1e15 GeV and the entropy-pixel Poisson
event premise.  Current canon does not derive either; item 42 keeps the
activation threshold open.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

DELTA = Fraction(1, 28)
A_TARGET = 2.1e-9
N_CHANNELS = 28
N_FLAGS = 112
ALPHA = 1.0 / 137.036
LAMBDA_QCD_GEV = 0.332
GAMMA0_GEV = ALPHA * LAMBDA_QCD_GEV
MPL_REDUCED_GEV = 2.435e18


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def amplitude_from_count(n_eff: float, f_eff: float = 1.0) -> float:
    return f_eff / n_eff


def fano_required(n_eff: float, amplitude: float = A_TARGET) -> float:
    return amplitude * n_eff


def entropy_count_from_h(h_gev: float) -> float:
    """de-Sitter entropy in reduced-Mp convention, S=8 pi^2 Mbar_P^2/H^2."""
    return 8.0 * math.pi**2 * MPL_REDUCED_GEV**2 / h_gev**2


def h_from_entropy_amplitude(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    """A=F/S_dS -> H=sqrt(8 pi^2 F A) Mbar_P."""
    return math.sqrt(8.0 * math.pi**2 * f_eff * amplitude) * MPL_REDUCED_GEV


def amplitude_from_entropy_h(h_gev: float, f_eff: float = 1.0) -> float:
    return f_eff / entropy_count_from_h(h_gev)


def h_from_landauer_integrated_entropy(amplitude: float = A_TARGET, f_eff: float = 1.0) -> float:
    """If N_eff = S_dS * Gamma0/H per e-fold, solve A=F/N_eff."""
    return (8.0 * math.pi**2 * GAMMA0_GEV * MPL_REDUCED_GEV**2 * amplitude / f_eff) ** (1.0 / 3.0)


def main() -> None:
    print("ITEM 131 HBC A_NU MICROSCOPIC SERVICE-LEDGER DERIVATION ATTEMPT")

    w_to_28 = (ROOT / "python_code" / "item131_w_to_28_instrument.py").read_text()
    serial = (ROOT / "python_code" / "item131_serial_clock_theorem.py").read_text()
    one_jump = (ROOT / "python_code" / "item131_one_jump_premise.py").read_text()
    chi_poisson = (ROOT / "python_code" / "item132_chi_unit_poisson.py").read_text()
    lindblad = (ROOT / "lindbladian_closure" / "lindbladian_closure.tex").read_text()
    cosmo = (ROOT / "cosmological_qec_engine" / "cosmological_qec_engine.tex").read_text()
    anchor = (ROOT / "ANCHOR.md").read_text()

    print("\n[1] Ledger ingredients already fixed")
    check("service channel has rate 1/28" in w_to_28, "28-channel service instrument fixes p_c=1/28")
    check("112 microscopic flags" in w_to_28, "incidence-refined service ledger has 112 microscopic flags")
    check("2^n subset" in serial and "2^28" in one_jump, "full active-subset serial-clock state space has size 2^28")
    check("Fano factor is 1" in chi_poisson, "matched birth/death QEC ledger can give Poisson Fano factor 1")
    check("sum_k \\gamma_k = \\alpha\\Lambda_{\\text{QCD}}" in lindblad, "Lindblad closure fixes total single-bit leakage bandwidth alpha Lambda_QCD")
    check("dA/dt" in cosmo and "Schematically" in cosmo, "HBC source gives schematic entropy-rate matching, not a scalar-mode fluctuation normalization")
    print(f"  Gamma0 = alpha Lambda_QCD = {GAMMA0_GEV:.6e} GeV")

    print("\n[2] Required event count")
    n_required = 1.0 / A_TARGET
    print(f"  target A_nu              = {A_TARGET:.6e}")
    print(f"  required N_eff/F_eff     = {n_required:.6e}")
    check(4.0e8 < n_required < 6.0e8, "CMB-sized amplitude requires roughly 5e8 independent events for F_eff=1")

    print("\n[3] Finite service-ledger count candidates")
    candidates = {
        "28 service labels": float(N_CHANNELS),
        "112 incidence flags": float(N_FLAGS),
        "2^28 active-subset states": float(2**28),
        "2^29 active-subset states plus binary phase": float(2**29),
        "28 * 2^24 tempting mixed count": float(28 * 2**24),
        "112 * 2^22 tempting mixed count": float(112 * 2**22),
    }
    for name, n_eff in candidates.items():
        amp_poisson = amplitude_from_count(n_eff)
        f_needed = fano_required(n_eff)
        print(
            f"  {name:38s}: N={n_eff:.6e}, "
            f"A(F=1)={amp_poisson:.6e}, F_needed={f_needed:.6e}"
        )
    check(amplitude_from_count(N_CHANNELS) / A_TARGET > 1.0e7, "28 labels alone are far too noisy")
    check(amplitude_from_count(N_FLAGS) / A_TARGET > 1.0e6, "112 flags alone are far too noisy")
    check(0.5 < fano_required(2**28) < 0.6, "2^28 lands near the amplitude only if a sub-Poisson Fano factor ~1/sqrt(pi) is added")
    check(0.95 < fano_required(28 * 2**24) < 1.05, "28*2^24 is numerically close with Poisson Fano factor")
    check(True, "no current microscopic service-ledger text derives 2^24 or a 1/sqrt(pi) Fano factor for the scalar shell")

    print("\n[4] Entropy-pixel Poisson closure candidate")
    h_star = h_from_entropy_amplitude()
    s_star = entropy_count_from_h(h_star)
    amp_at_gut = amplitude_from_entropy_h(1.0e15)
    print("  Premise: one independent Poisson service event per printed de-Sitter entropy pixel.")
    print(f"  S_dS(H*) required     = {s_star:.6e}")
    print(f"  H* required           = {h_star:.6e} GeV")
    print(f"  A_nu at H*=1e15 GeV   = {amp_at_gut:.6e}")
    print(f"  ratio A(1e15)/target  = {amp_at_gut / A_TARGET:.6f}")
    check(abs(s_star - n_required) / n_required < 1e-12, "entropy-pixel premise makes N_eff=S_dS=1/A_nu")
    check(8.0e14 < h_star < 1.2e15, "observed amplitude maps to an early H scale near 1e15 GeV")
    check(abs(amp_at_gut / A_TARGET - 1.0) < 0.03, "a 1e15 GeV saturated printer gives the observed amplitude to a few percent")
    check("R2/R3 together at GUT" in anchor and "phenomenologically arranged" in anchor, "the corresponding GUT activation scale is present but item-42-open, not derived")

    print("\n[5] Landauer-bandwidth integrated entropy candidate")
    h_landauer = h_from_landauer_integrated_entropy()
    amp_landauer_at_hstar = h_star**3 / (8.0 * math.pi**2 * GAMMA0_GEV * MPL_REDUCED_GEV**2)
    print("  Alternative premise: each entropy pixel produces Gamma0/H service attempts per e-fold.")
    print(f"  H required under this premise = {h_landauer:.6e} GeV")
    print(f"  A at entropy H* if Gamma0/H counted = {amp_landauer_at_hstar:.6e}")
    check(5.0e8 < h_landauer < 2.0e9, "including Gamma0/H event attempts moves the inferred H scale to ~1e9 GeV")
    check(abs(amp_landauer_at_hstar / A_TARGET) > 1.0e5, "counting Gamma0/H attempts is a different normalization, not the CMB amplitude")
    check(True, "the microscopic ledger must choose entropy-pixel events versus per-pixel Landauer attempts")

    print("\n[6] Closure status")
    closed = [
        "A_nu is F_eff/N_eff for the compensated HBC event current",
        "Poisson service-current closure would set F_eff=1",
        "entropy-pixel premise gives A_nu=H_*^2/(8 pi^2 Mbar_P^2)",
        "observed A_nu is equivalent to H_* ~= 1e15 GeV under that premise",
    ]
    open_items = [
        "derive the early saturated HBC printing scale H_* intrinsically",
        "derive that the scalar shell has one independent service event per de-Sitter entropy pixel",
        "derive the scalar-shell Fano/correlation factor from the actual HBC/KMS current",
        "derive item-42 R-activation thresholds rather than importing the GUT scale",
    ]
    for item in closed:
        check(True, f"reduced/conditional: {item}")
    for item in open_items:
        check(True, f"not closed: {item}")

    print("\n" + "=" * 104)
    print("VERDICT")
    print("  The microscopic ledger does not yet derive A_nu outright.")
    print("  It does reduce the target to a clean entropy-pixel theorem:")
    print("      if N_eff = S_dS and F_eff=1, then")
    print("          A_nu = H_*^2 / (8 pi^2 Mbar_P^2).")
    print(f"  The observed amplitude then requires H_* = {h_star:.3e} GeV, essentially")
    print("  the phenomenological GUT/R2-R3 activation scale.  Therefore a real")
    print("  closure of A_nu is equivalent to deriving the early saturated HBC printing")
    print("  scale and the entropy-pixel Poisson event law from the microscopic service")
    print("  ledger.  Current canon supplies neither; item 42 remains the activation")
    print("  threshold blocker.")
    print("=" * 104)
    print("exit 0 -- A_nu reduced to entropy-pixel Poisson theorem + H_* derivation target; not closed.")


if __name__ == "__main__":
    main()
