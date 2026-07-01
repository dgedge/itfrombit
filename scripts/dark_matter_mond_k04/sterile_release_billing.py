#!/usr/bin/env python3
r"""BILLING THEOREM (pre-registered): the alpha_0-power of the sterile nu_R release (CMB last gate).

The CMB dark-reservoir absolute density rides on the source law  n_nuR/n_gamma = alpha_0^p / 208.
M15 derived two of the three factors; the OPEN piece (l.1909) is the POWER p — "the release carries
one alpha_0 non-unitary billing through the Q-complement". This pre-registers and tests p.

================================ FROZEN PRE-REGISTRATION ================================
ALLOWED LEDGER (only these):
  - §5.9 billing rule (ADOPTED canon): one alpha_0 per non-unitary syndrome firing/erasure
    (w=alpha_0 Lambda = irreversibility-fraction x energy/tick; item 79 / mass sector / item 119).
  - one S_3-singlet release PORT (DERIVED, M15 l.1913: the 3 sterile corners -> bright singlet).
  - uniform Q-addressing 1/|Q|=1/208 (DERIVED, M15 l.1913: Evans-Frigerio on the connected Q-graph).
  - the sterile mass m_nuR = alpha_0^2 Lambda (separate canon claim; needed for the ABSOLUTE density).
FORBIDDEN: shopping p to fit; cross-sector imports; tuning m_nuR or the 208.
ALTERNATIVES tested: p in {0, 1, 2}  (0 = unitary/free release; 1 = one §5.9 firing; 2 = two firings).
§5.9 PREDICTION: p = 1 (one singlet port = one non-unitary firing).
KILL CRITERIA:
  PASS (conditional theorem)  iff the §5.9-predicted p=1 is also the data-required power AND the
        alternatives p=0,2 are excluded; residual = the adopted §5.9 rule.
  FAIL (§5.9 refuted here)     iff the data requires p != 1.
  UNDERSPECIFIED               iff the release port / billing is not canonically fixed (it is: l.1913).
========================================================================================
"""
ALPHA0   = 1.0 / 137.0
LAMBDA   = 0.332             # GeV
N_GAMMA  = 411.0            # cm^-3
RHO_CRIT = 1.0537e4         # eV/cm^3 (rho_crit/h^2)
M_NUR    = ALPHA0**2 * LAMBDA * 1e9     # eV  (= alpha_0^2 Lambda)
REQUIRED = 0.0242           # omega_nuR h^2 required (the 20% dark share; from z_eq closure 0.024+0.096=0.120)

def omega_nuR(p):
    n_nuR = (ALPHA0**p / 208.0) * N_GAMMA          # cm^-3
    return n_nuR * M_NUR / RHO_CRIT                 # Omega_nuR h^2

w = {p: omega_nuR(p) for p in (0, 1, 2)}
P_PRED = 1                                          # §5.9 prediction

# data-required power = the p whose omega matches REQUIRED
best_p = min((0, 1, 2), key=lambda p: abs(w[p] / REQUIRED - 1.0))

# --- the power selection: §5.9 predicts 1; data must agree; alternatives excluded by 1/alpha_0=137x ---
assert best_p == P_PRED, "the data-required alpha_0-power must match the §5.9 prediction p=1"
assert abs(w[1] / REQUIRED - 1.0) < 0.05, "p=1 matches the required omega_nuR ~ 0.0242"
assert w[0] / REQUIRED > 100, "p=0 (free/unitary release) OVERPRODUCES by ~1/alpha_0 = 137x -> excluded"
assert w[2] / REQUIRED < 0.02, "p=2 (two firings) UNDERPRODUCES by ~alpha_0 = 1/137x -> excluded"
# the separation between adjacent powers is exactly 1/alpha_0 (robust to m_nuR):
assert abs((w[0] / w[1]) - 1.0 / ALPHA0) < 1e-6 and abs((w[1] / w[2]) - 1.0 / ALPHA0) < 1e-6, \
    "adjacent powers separated by 1/alpha_0 = 137x -> the power is robustly selected, independent of m_nuR"

# --- downstream: the absolute dark budget + z_eq (conditional on p=1 AND m_nuR=alpha_0^2 Lambda) ---
omega_zero = 4.0 * w[1]                             # 4:1 directed-incidence split (derived)
omega_dark = w[1] + omega_zero
OMEGA_R = 2.473e-5 * (1.0 + 0.2271 * 3.046)
z_eq_dust = (0.0224 + omega_dark) / OMEGA_R - 1.0
assert abs(omega_dark - 0.121) < 0.004 and 3300 < z_eq_dust < 3600, "dark budget 0.121, z_eq ~ 3400 (third peak)"

bar = "=" * 92
print(bar)
print("STERILE nu_R RELEASE — alpha_0-power BILLING THEOREM (pre-registered; CMB last native gate)")
print(bar)
print(f"  m_nuR = alpha_0^2 Lambda = {M_NUR/1e3:.2f} keV ;  required omega_nuR h^2 = {REQUIRED}")
print(f"\n  p   n_nuR/n_gamma = alpha_0^p/208      omega_nuR h^2      vs required")
print(f"  {'-'*70}")
for p in (0, 1, 2):
    tag = "  <- §5.9 PREDICTION + data-selected" if p == 1 else ("  (overproduces 137x)" if p == 0 else "  (underproduces 137x)")
    print(f"  {p}   {ALPHA0**p/208.0:.4e}                 {w[p]:.5f}        x{w[p]/REQUIRED:7.3f}{tag}")
print(f"\n  power selection: §5.9 predicts p=1 (one singlet port = one non-unitary firing); data requires p={best_p};")
print(f"  adjacent powers separated by 1/alpha_0 = 137x -> p=1 robustly selected (independent of m_nuR).")
print(f"  downstream (conditional on p=1 + m_nuR): omega_dark = {omega_dark:.4f}, z_eq = {z_eq_dust:.0f} (third peak closes)")
print(f"""
{bar}
VERDICT (billing theorem, exit 0):  CONDITIONAL THEOREM — alpha_0^1 is §5.9-PREDICTED and data-selected.

  The sterile-release billing POWER is alpha_0^1, derived from the adopted §5.9 rule (one alpha_0
  per non-unitary firing) applied to the one S_3-singlet port (derived) -- NOT shopped: §5.9
  PREDICTS p=1 before the comparison, and the data then requires exactly p=1, with the alternatives
  p=0 (free release) and p=2 (two firings) each excluded by 1/alpha_0 = 137x. So the source law
  n_nuR/n_gamma = alpha_0/208 has all three factors accounted: the power (§5.9-predicted+selected),
  the /208 (Evans-Frigerio, derived), the one port (S_3, derived).

  This CONVERTS the CMB absolute-density gate from "conditional candidate" (M15 l.1909) to a
  conditional THEOREM: the alpha_0-power is no longer an open premise but a §5.9 prediction the data
  confirms. RESIDUALS (honest): (a) the adopted §5.9 billing rule itself (item-79 tier, conditional);
  (b) the ABSOLUTE density (omega_nuR=0.024) additionally needs m_nuR=alpha_0^2 Lambda -- the power
  selection is robust to m_nuR (137x steps), but the absolute landing rides on it; (c) the standing
  Boltzmann/halo gates (M15) are unaffected. The CMB sector's last NATIVE billing gate is closed
  conditional on §5.9.
{bar}""")
print(f"exit 0 -- sterile-release billing: alpha_0^1 §5.9-predicted + data-selected (p=0,2 excluded 137x); "
      f"CMB absolute-density gate -> conditional theorem; residual = adopted §5.9 rule + m_nuR.")
