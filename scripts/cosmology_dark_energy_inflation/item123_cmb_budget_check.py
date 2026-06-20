#!/usr/bin/env python3
"""
item123_cmb_budget_check.py

Cosmology / CMB third-peak flag (DRIFT M15): the framework's MOND+nu_R+R4-DE
cosmology fails the third acoustic peak unless a cold pressureless component of
omega ~ 0.096 exists.  M15's live mechanism: an R4 "zero-mode reservoir" supplies
it, with absolute density from the boot source law n_nuR/n_gamma = alpha0/208,
m_nuR = alpha0^2 * Lambda_QCD, and a directed-incidence 4:1 split
(omega_zero = 4 omega_nuR), giving omega_dark = 5 omega_nuR.

This script verifies the load-bearing ARITHMETIC of that mechanism from the
framework's own derived inputs, and checks whether it (i) restores matter-
radiation equality z_eq (the third-peak epoch), (ii) gives a LCDM-like baryon
loading, (iii) closes flat.  Honest about what remains (a Boltzmann run + the
dust-charge derivation).

Inputs that are FRAMEWORK-DERIVED (not fitted):
  eta = (3/14) alpha0^4                         (baryogenesis, item 126)
  n_nuR/n_gamma = alpha0/208, m_nuR=alpha0^2 Lambda   (dark source law, item 123)
  omega_dark = 5 omega_nuR                       (nu_R 20% + zero-mode 80%, M15)
  Omega_Lambda = 12 pi / 55, H0 = 67.43          (Bekenstein C=55/8, item 22)
Standard physical constants (CODATA / Planck) are cited inline.

Self-asserting; numpy only.
"""
import sys
import numpy as np

a0 = 1 / 137.036                 # alpha0 (fine-structure; bare 1/137 within 0.03%)
Lam_keV = 331.7e3               # Lambda_QCD = Lambda_p = m_p/(2 sqrt2) = 0.3317 GeV
n_gamma = 410.7                  # CMB photons /cm^3 (T=2.7255 K)
rho_c100 = 1.0537e4             # critical density /h^2 in eV/cm^3
omega_r = 4.15e-5               # radiation density omega_r = omega_gamma+omega_nu (3 active nu)
_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# ---- baryons from baryogenesis eta ----------------------------------------
eta = (3 / 14) * a0**4
m_p_eV = 938.272e6
omega_b = eta * n_gamma * m_p_eV / rho_c100
print(f"eta = (3/14) alpha0^4 = {eta:.3e}  ->  omega_b = {omega_b:.4f}  (std ~0.0224)")
check("omega_b from baryogenesis eta matches standard (within 5%)",
      abs(omega_b - 0.0224) / 0.0224 < 0.05)

# ---- sterile nu_R relic from the boot source law --------------------------
n_nuR = (a0 / 208) * n_gamma
m_nuR_keV = a0**2 * Lam_keV
omega_nuR = (n_nuR * m_nuR_keV * 1e3) / rho_c100      # keV->eV
print(f"m_nuR = alpha0^2 Lambda = {m_nuR_keV:.2f} keV;  "
      f"n_nuR/n_gamma = alpha0/208 = {a0/208:.3e}")
print(f"omega_nuR = {omega_nuR:.5f}  (canon 0.02418)")
check("omega_nuR from (alpha0/208, alpha0^2 Lambda) matches canon 0.0242",
      abs(omega_nuR - 0.02418) / 0.02418 < 0.02)

# ---- dark budget: nu_R (20%) + R4 zero-mode reservoir (80%) ----------------
omega_dark = 5 * omega_nuR                            # omega_zero = 4 omega_nuR
omega_m = omega_b + omega_dark
print(f"omega_dark = 5 omega_nuR = {omega_dark:.4f}  (nu_R 0.024 + zero-mode 0.097)")
print(f"omega_m = omega_b + omega_dark = {omega_m:.4f}")
check("omega_dark ~ 0.121 (the pressureless 0.096 component + nu_R)",
      abs(omega_dark - 0.121) < 0.004)

# ---- the kill-shot: matter-radiation equality z_eq ------------------------
z_eq_full = omega_m / omega_r - 1
z_eq_nuRonly = (omega_b + omega_nuR) / omega_r - 1
print(f"\nz_eq with full dark budget = {z_eq_full:.0f}  (LCDM ~3400; third-peak epoch OK)")
print(f"z_eq with nu_R only (the M15 kill-shot) = {z_eq_nuRonly:.0f}  "
      f"(~z_rec 1100 -> third peak FAILS)")
check("full budget RESTORES z_eq into the LCDM band (3000-3600)",
      3000 < z_eq_full < 3600)
check("nu_R-only collapses z_eq to ~recombination (the failure being fixed)",
      z_eq_nuRonly < 1300)

# ---- baryon loading (third-peak height driver) ----------------------------
ratio_bc = omega_b / omega_dark
print(f"\nbaryon loading omega_b/omega_c = {ratio_bc:.3f}  (LCDM ~0.187)")
check("baryon loading is LCDM-like (third-peak heights ~ LCDM)",
      abs(ratio_bc - 0.187) / 0.187 < 0.10)

# ---- flatness cross-check (secondary; see caveat) -------------------------
H0 = 67.43
h = H0 / 100
Omega_m = omega_m / h**2
Omega_L = 12 * np.pi / 55
print(f"\nflatness check (H0={H0}, Omega_Lambda=12pi/55={Omega_L:.4f}):")
print(f"  Omega_m = {Omega_m:.4f};  Omega_m + Omega_Lambda = {Omega_m+Omega_L:.4f}")
check("budget closes FLAT: Omega_m + Omega_Lambda = 1 (within 1%)",
      abs(Omega_m + Omega_L - 1.0) < 0.01)
print("  [caveat: H0 via the Part-20 M_P relation may itself assume flatness, so")
print("   this closure can be partly built-in; z_eq + baryon-loading above are NOT.]")

# ---- verdict --------------------------------------------------------------
print("\n--- VERDICT ---")
print("The M15 third-peak FAILURE is resolved IN PRINCIPLE by the R4 zero-mode")
print("reservoir: the framework's OWN source law (alpha0/208, alpha0^2 Lambda, x5)")
print(f"gives omega_dark={omega_dark:.3f}, restoring z_eq {z_eq_nuRonly:.0f}->{z_eq_full:.0f}")
print("and a LCDM-like baryon loading -> the matter budget is now third-peak-viable.")
print("The flag moves from 'fails, no native escape' to CONDITIONAL/LIVE.")
print()
print("STILL OPEN (necessary, not yet done):")
print(" - a full Boltzmann (CAMB/CLASS) run: z_eq + loading are necessary, NOT")
print("   sufficient for the actual peak heights.")
print(" - the zero-mode as GENUINE conserved pressureless dust: the Brown-Kuchar")
print("   Hamiltonian form is derived in canon, but the absolute source law")
print("   alpha0/208 is a candidate (denominators 208-211 all land <1%), and the")
print("   R4-DE law currently contradicts a homogeneous dust charge.")
print(" - halo double-counting (zero-mode CDM vs active MOND) -- the phenomenology gate.")
print("So: NOT closed, but the budget is now numerically self-consistent and LCDM-")
print("like; the residual is a derivation (dust charge) + a Boltzmann run, not a")
print("factor-5 deficit.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
