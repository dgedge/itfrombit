#!/usr/bin/env python3
r"""CMB DARK-RESERVOIR — completion status + end-to-end arithmetic verification (M15 / item 123).

The recorded next target. Reading M15 shows it is NOT a bare open binary but a LIVE CONDITIONAL
MECHANISM (2026-06-15): the R4 zero-mode reservoir exists (conserved Brown-Kuchar dust, w=c_s^2=0),
the 4:1 relative split is derived, the dust Hamiltonian / uniform-Q / singlet-source / continuum
lift are closed, the absolute density is a conditional alpha_0/208 theorem,
and the external CAMB TT check has now been run.  This script keeps the
end-to-end arithmetic visible and classifies the remaining open gates.

VERIFIED canon (M15, DRIFT l.1893-1913; owner scripts all exit 0):
  third peak needs Omega_c h^2 ~ 0.120; nu_R relic supplies only 0.024 (the kill-shot:
  baryons+nu_R alone give z_eq ~ z_rec). The R4 zero-mode reservoir supplies the rest IF the
  alpha_0/208 source law + the 4:1 directed-incidence split hold.  The
  companion CAMB scripts then show that the pressureless slot restores the
  relative third peak and that the full framework TT spectrum matches Planck
  peak heights/lensing at ~1% RMS; the residual is the acoustic-scale joint
  fit, not the third-peak height.
"""
ALPHA0   = 1.0 / 137.0
LAMBDA   = 0.332            # GeV, Lambda_QCD (canon anchor)
N_GAMMA  = 411.0           # cm^-3, CMB photon number density today
RHO_CRIT = 1.0537e4        # eV/cm^3, rho_crit / h^2
OMEGA_B  = 0.0224          # baryon omega_b h^2
# radiation: omega_r = omega_gamma (1 + 0.2271 N_eff)
OMEGA_R  = 2.473e-5 * (1.0 + 0.2271 * 3.046)
OMEGA_C_LCDM = 0.120       # the CDM the third peak requires

def z_eq(omega_m):
    return omega_m / OMEGA_R - 1.0

# ---------------------------------------------------------------------------
# (1) absolute nu_R density from the source law  n_nuR/n_gamma = alpha_0/208,  m_nuR = alpha_0^2 Lambda
# ---------------------------------------------------------------------------
m_nuR    = ALPHA0**2 * LAMBDA * 1e9          # eV  (alpha_0^2 Lambda_QCD)
ratio    = ALPHA0 / 208.0                     # n_nuR / n_gamma
n_nuR    = ratio * N_GAMMA                     # cm^-3
omega_nuR = (n_nuR * m_nuR) / RHO_CRIT        # Omega_nuR h^2
# ---------------------------------------------------------------------------
# (2) the 4:1 directed-incidence split: omega_zero = 4 * omega_nuR  (2x6 directed edges / 3 corners)
# ---------------------------------------------------------------------------
omega_zero = 4.0 * omega_nuR
omega_dark = omega_nuR + omega_zero

print("=" * 92)
print("CMB DARK-RESERVOIR — live conditional route, calculated end to end (M15 / item 123)")
print("=" * 92)
print(f"  source law:  m_nuR = alpha_0^2 Lambda   = {m_nuR/1e3:.2f} keV   (a WDM-scale sterile)")
print(f"               n_nuR/n_gamma = alpha_0/208 = {ratio:.3e}")
print(f"  => omega_nuR        = {omega_nuR:.5f}   (canon 0.02418)")
print(f"     omega_zero = 4x  = {omega_zero:.5f}   (4:1 directed-incidence split)")
print(f"     omega_dark       = {omega_dark:.5f}   (canon 0.12089; vs LCDM omega_c={OMEGA_C_LCDM})")

# ---------------------------------------------------------------------------
# (3) the z_eq kill-shot: does the dust close the third-peak epoch?
# ---------------------------------------------------------------------------
z_nodust = z_eq(OMEGA_B + omega_nuR)                       # baryons + nu_R only
z_dust   = z_eq(OMEGA_B + omega_dark)                       # baryons + nu_R + zero-mode dust
z_lcdm   = z_eq(OMEGA_B + OMEGA_C_LCDM)                     # LCDM reference
print(f"\n  z_eq (no dust: b+nu_R, omega_m={OMEGA_B+omega_nuR:.4f}) = {z_nodust:.0f}   <- KILL-SHOT (~z_rec 1100)")
print(f"  z_eq (dust:  b+dark, omega_m={OMEGA_B+omega_dark:.4f}) = {z_dust:.0f}   <- closes (~LCDM {z_lcdm:.0f})")

# ---------------------------------- assertions ----------------------------------
assert abs(omega_nuR - 0.0242) < 0.001,        "nu_R absolute density from alpha_0/208 ~ 0.024 (the 20% share)"
assert abs(omega_zero - 4 * omega_nuR) < 1e-9, "zero-mode = 4 x nu_R (directed-incidence 4:1)"
assert abs(omega_dark - 0.121) < 0.003,        "total dark = 0.121 (vs LCDM 0.120, +0.74%)"
assert 1000 < z_nodust < 1250,                 "no-dust z_eq collapses to ~recombination (the kill-shot)"
assert 3200 < z_dust < 3600,                   "dust route restores z_eq ~ 3400 (third-peak epoch closes)"
assert z_dust / z_nodust > 2.8,                "the dust shifts z_eq by ~3x -- the decisive third-peak fix"

# ---------------------------------------------------------------------------
# (4) completion roadmap: what is derived vs conditional vs open
# ---------------------------------------------------------------------------
status = {
    "zero-mode reservoir EXISTS (conserved, Brown-Kuchar dust, w=c_s^2=0)": "derived (l.1905-1913)",
    "relative split omega_zero = 4 omega_nuR (directed incidence)":          "derived (l.1905)",
    "dust Hamiltonian form (rest-count, w=c_s^2=0)":                         "derived (l.1907)",
    "uniform Q-addressing (|Q|=208, Evans-Frigerio)":                        "derived (l.1913)",
    "generation-blind singlet source (S_3 -> bright singlet)":               "derived (l.1913)",
    "continuum lift (Brown-Kuchar geodesic dust)":                           "derived (l.1913)",
    "absolute density (alpha_0/208 source law -> omega_dark=0.121)":         "conditional theorem (§5.9/item-79 tier; 0.74% high)",
    "    -> one alpha_0 non-unitary billing per sterile release":            "conditional theorem (sterile_release_billing.py)",
    "CAMB third-peak / TT peak-height run":                                  "closed (peak heights+lensing; theta_* selector gate remains)",
    "theta_* repair surface":                                                "quantified (H0=66.318, or Omega_k=0.00212, or phantom-flipped w)",
    "63/64 acoustic pre-latch projection":                                   "theorem target (numerically closes theta_*; not promoted)",
    "galaxy halo non-double-counting":                                       "branch choice: zero-mode CDM halo OR derive >95.9% depletion for active MOND",
}
print(f"\n{'-'*92}\n  COMPLETION ROADMAP:")
for k, v in status.items():
    v_lower = v.lower()
    tag = "OK  " if v_lower.startswith(("derived", "closed")) else ("cond" if v_lower.startswith("conditional") else "OPEN")
    print(f"    [{tag}] {k:62s} {v}")

print(f"""
{"=" * 92}
VERDICT (status + arithmetic, exit 0):  LIVE CONDITIONAL MECHANISM — the dust route closes z_eq.

  Calculated end to end: the alpha_0/208 source law + the 4:1 split give omega_dark = {omega_dark:.4f}
  (vs LCDM 0.120), which SHIFTS z_eq from {z_nodust:.0f} (the kill-shot, ~recombination, third peak
  wrong) to {z_dust:.0f} (~LCDM, third-peak epoch restored). So the live route is numerically on
  target -- the third-peak no-go is closed by the zero-mode dust, conditional on the source law.

  REMAINING TO COMPLETE (2 sharpened gates):
    (a) the selector/acoustic-scale gate: exact theta_* closure is numerically
        achieved by H_CMB=(63/64)H_selector, but no current service-ledger
        theorem proves the acoustic ruler sees one pre-latch depth-6 slot less;
    (b) halo branch choice: zero-mode as CDM-like halo is non-double-counted, but retaining
        active MOND/RAR needs >95.9% galaxy depletion/screening of the zero-mode.
  Plus the standing honesty flag (M15): this is the framework's MOST AT-RISK sector -- Planck has
  already measured the acoustic spectrum, so a failed joint-fit/halo gate is real falsification,
  not a soft gap.
{"=" * 92}""")
print(f"exit 0 -- CMB dark-reservoir LIVE CONDITIONAL: omega_dark={omega_dark:.4f} shifts z_eq {z_nodust:.0f}->{z_dust:.0f} "
      f"(kill-shot closed); CAMB peak-height gate closed; residual: selector/theta_* gate + halo branch choice.")
