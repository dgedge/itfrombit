#!/usr/bin/env python3
"""
item87_neutrino_mass_spectrum.py

PMNS / neutrino sector update.  The PMNS ANGLE running is open and the delta_nu
CP phase is a geometric-primitive no-go (item87_MR_derivation_attempt.py).  This
script checks the OTHER half: is the geometric primitive delta_nu=1/3 at least
ECONOMICAL -- does it (with R_nu=1) predict the neutrino MASS observables?

Neutrino Koide circulant (canon ch10 / item 83): masses
    m_n proportional to (1 + R_nu cos theta_n)^2,  theta_n = delta_nu + 2 pi n/3,
with the framework's geometric inputs R_nu = 1 and delta_nu = 1/3 (radians).

Predictions are PARAMETER-FREE given those two inputs (one overall scale mu is
fixed by Delta m^2_31).  Compared to NuFIT/Planck.  Self-asserting; numpy only.
"""
import sys
import numpy as np

_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# ---- the two geometric primitives -----------------------------------------
R_nu = 1.0
delta_nu = 1.0 / 3.0                       # radians (the item-86 defect-fraction primitive)
n = np.array([0, 1, 2])
theta = delta_nu + 2 * np.pi * n / 3
m_raw = (1 + R_nu * np.cos(theta))**2      # mass values up to overall scale
m_raw = np.sort(m_raw)                      # ascending -> (m1,m2,m3), normal ordering
print(f"theta_n = {np.round(theta,4)};  m_raw (sorted) = {np.round(m_raw,4)}")

# ---- Koide ratio Q_nu (scale-free) ----------------------------------------
Q_nu = m_raw.sum() / (np.sqrt(m_raw).sum())**2
print(f"Koide Q_nu = {Q_nu:.6f}  (framework claim 1/2)")
check("Q_nu = 1/2 exactly (analytic for R_nu=1, any delta)", abs(Q_nu - 0.5) < 1e-9)

# ---- mass-squared-splitting ratio (parameter-free; depends only on delta_nu)
ratio_pred = (m_raw[1]**2 - m_raw[0]**2) / (m_raw[2]**2 - m_raw[0]**2)
dm21_meas, dm31_meas = 7.42e-5, 2.515e-3   # eV^2, NuFIT 5.x normal ordering
ratio_meas = dm21_meas / dm31_meas
print(f"\nDelta m^2_21 / Delta m^2_31 = {ratio_pred:.4f}  vs measured "
      f"{ratio_meas:.4f}  ({100*abs(ratio_pred-ratio_meas)/ratio_meas:.1f}%)")
check("Delta m^2 ratio predicted from delta_nu=1/3 matches data within 8%",
      abs(ratio_pred - ratio_meas) / ratio_meas < 0.08)

# ---- absolute spectrum (fix the one scale mu by Delta m^2_31) --------------
mu2 = dm31_meas / (m_raw[2]**2 - m_raw[0]**2)
mu = np.sqrt(mu2)
m = mu * m_raw                              # eV
print(f"absolute masses (meV): m1={m[0]*1e3:.2f}, m2={m[1]*1e3:.2f}, "
      f"m3={m[2]*1e3:.2f}  (ch10: 0.79, 8.71, 50.12)")
check("absolute masses match canon ch10 (0.79, 8.71, 50.1 meV) within 3%",
      abs(m[0]*1e3 - 0.79) < 0.05 and abs(m[1]*1e3 - 8.71) < 0.3
      and abs(m[2]*1e3 - 50.1) < 1.5)

dm21_pred = m[1]**2 - m[0]**2
print(f"Delta m^2_21 = {dm21_pred:.3e} eV^2 vs measured {dm21_meas:.3e} "
      f"({100*abs(dm21_pred-dm21_meas)/dm21_meas:.1f}%)")
check("Delta m^2_21 absolute matches data within 8%",
      abs(dm21_pred - dm21_meas) / dm21_meas < 0.08)

# ---- sum of masses (cosmology) + ordering ---------------------------------
Sigma = m.sum()
print(f"Sigma m_nu = {Sigma*1e3:.1f} meV  (Planck bound < 120 meV; DESI/Euclid testable)")
check("Sigma m_nu below the cosmological bound (0.12 eV)", Sigma < 0.12)
check("normal ordering (m1<m2<m3)", m[0] < m[1] < m[2])

# ---- m_betabeta (0nubb): mixing-DEPENDENT -> flags the open part -----------
# with TBM mixing |U_e1|^2=2/3, |U_e2|^2=1/3, |U_e3|^2=0 and no Majorana phases:
mbb_tbm = abs(2 / 3 * m[0] + 1 / 3 * m[1])
print(f"\nm_betabeta (TBM, no Majorana phases) = {mbb_tbm*1e3:.2f} meV  "
      f"(ch10 ~3.7; LEGEND-1000/nEXO target) -- MIXING/phase-dependent (open)")

# ===========================================================================
print("\n--- VERDICT ---")
print("UPDATE to the PMNS/neutrino item: the sector splits cleanly.")
print(" * MASSES (this script): the single geometric primitive delta_nu=1/3 (with")
print(f"   R_nu=1) PREDICTS the whole spectrum parameter-free -- Q_nu=1/2 (exact),")
print(f"   Dm^2 ratio {ratio_pred:.4f} vs {ratio_meas:.4f} ({100*abs(ratio_pred-ratio_meas)/ratio_meas:.1f}%),")
print(f"   masses (0.79, 8.72, 50.2) meV, Sigma={Sigma*1e3:.0f} meV, normal ordering.")
print("   So delta_nu, though NOT walk-derivable (the item-87 no-go), is ECONOMICAL:")
print("   one un-derived geometric input -> the full neutrino mass observables.")
print(" * MIXING angles: still OPEN -- both substrate routes overshoot theta_13")
print("   (cross-talk 29 deg, seesaw ~30 deg, vs measured 8.6 deg).")
print(" * CP phase: root-caused -- walk J=0 (I3/chi-gated); a complex Majorana M_R")
print("   supplies CP with sign tracking sign(delta_nu); its derivation is the no-go.")
print("So the neutrino MASS sector is a clean parameter-free win; the MIXING + the")
print("CP-phase derivation remain the open frontier.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
