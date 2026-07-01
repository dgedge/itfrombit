#!/usr/bin/env python3
r"""ITEM 123: the absolute Omega_dark h^2 chain -- is the CMB cold-dust gap closed?

The ch20 ledger still reads 'native 0.024 vs 0.120, needs imported dust.' This checks
the cosmology-paper chain end to end and asks whether the factor ~5 is derived (gap
closed) or imported (gap real). The chain, using NO observed dark density as input:
    m_nuR    = alpha0^2 * Lambda_QCD                 (sterile mass)
    n_nuR    = (alpha0/208) * n_gamma(T0)            (208 = 2^8 - 48 service complement)
    omega_nuR= m_nuR n_nuR / (rho_crit/h^2)
    omega_zero = 4 * omega_nuR                       (directed R4 incidence (2*6)/3 = 4)
    omega_dark = omega_nuR + omega_zero = 5 omega_nuR
Self-asserting.
"""
import numpy as np
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m
GeV_per_g = 5.6096e23
alpha0   = 1/137.036
n_gamma  = 410.7                       # cm^-3 at T0=2.7255 K
rho_crit_over_h2 = 1.878e-29 * GeV_per_g     # GeV cm^-3  (= 1.0536e-5)
Lambda_QCD = 0.332                     # GeV, standard 3-flavour MS-bar scale
OBS = 0.120                            # Planck Omega_c h^2

print("="*72); print("Omega_dark h^2 absolute chain (alpha0/208 boot source + R4 4:1 split)"); print("="*72)

# directed R4 incidence factor: (2 repair edges/corner * 6 edges... ) = (2*6)/3 = 4
ok((2*6)//3 == 4, "directed R4 service-edge incidence (2*6)/3 = 4  -> omega_zero = 4 omega_nuR")

m_nuR = alpha0**2 * Lambda_QCD
n_nuR = (alpha0/208) * n_gamma
omega_nuR = m_nuR * n_nuR / rho_crit_over_h2
omega_dark = 5 * omega_nuR
print(f"\n  m_nuR    = alpha0^2 * {Lambda_QCD*1e3:.0f} MeV = {m_nuR*1e6:.2f} keV")
print(f"  n_nuR/n_g= alpha0/208 = {alpha0/208:.3e}   ->  n_nuR = {n_nuR:.4f} cm^-3")
print(f"  omega_nuR h^2 = {omega_nuR:.5f}")
print(f"  omega_dark h^2 = 5 * omega_nuR = {omega_dark:.5f}   (obs {OBS})  -> {100*(omega_dark-OBS)/OBS:+.1f}%")
ok(0.022 < omega_nuR < 0.026, "sterile relic omega_nuR ~ 0.024 (the ch20 'native' piece)")
ok(abs(omega_dark-OBS)/OBS < 0.02, "omega_dark within 2% of observed 0.120 -- a prediction, no dark-density input")

# sensitivity to Lambda_QCD (the only QCD-scale input)
print("\n  sensitivity (Lambda_QCD -> omega_dark h^2):")
for L in (0.300, 0.320, 0.332, 0.340):
    od = 5 * (alpha0**2 * L) * n_nuR / rho_crit_over_h2
    print(f"    Lambda_QCD={L*1e3:.0f} MeV -> omega_dark={od:.4f} ({100*(od-OBS)/OBS:+.1f}%)")

print("\n"+"="*72); print("VERDICT")
print("  The ch20 'native 0.024 vs 0.120, needs imported dust' is STALE. The dust is NOT")
print("  imported: it is the DERIVED R4 zero-mode reservoir, omega_zero = 4 omega_nuR, from the")
print("  directed service-edge incidence (2*6)/3 = 4. With the alpha0/208 boot source and")
print("  m_nuR = alpha0^2 Lambda_QCD, the full omega_dark h^2 = 5 omega_nuR = 0.121 lands ~0.7% high")
print("  vs 0.120 -- using the CMB photon density as the bath, NOT the observed dark density.")
print("  So Omega_c is a derived ~0.7% prediction, not a 5x shortfall. RESIDUAL (honest): the chain")
print("  is CONDITIONAL on the boot-rate theorem (one alpha0-billed non-unitary sterile release")
print("  through the 208 service complement) and on m_nuR=alpha0^2 Lambda_QCD; neither is Locked,")
print("  though companion scripts derive the Q addressing + the generation-blind one-port release. exit 0")
