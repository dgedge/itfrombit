#!/usr/bin/env python3
r"""THE BOUNDARY-PRINTING LEDGER — expansion as a printer, and its one decisive exposure.

Consolidates the canon claim that space is not stretched but PRINTED at the causal boundary
(ANCHOR item 123 CLOSURE II, Holographic Boundary Crystallization: "the universe does not
stretch; it prints new boundary nodes to accommodate its own computational exhaust"). If so,
inflation, dark energy, and horizon entropy are one process — "printer dynamics" — at
different clock rates (item 128: inflation = the saturated/high-bandwidth phase; dark energy =
the steady state). This script ties the printer's TWO control parameters,

    * the EVENT UNIT     u_event = ln 2      (one Landauer bit per printed boundary cell)
    * the STOP RULE      N_shell * alpha0^4 = C_F = 4/3   (saturated-printer/current law,
                         the lambda_shell = C_F saturation identity, item131_saturation_closure.py)

to FOUR observables, and then states the single falsifier the whole picture rests on.

WHAT THIS SCRIPT IS: a SYNTHESIS + one new computation (the tensor exposure), not a fresh
derivation. The scalar amplitude A_nu=(3/4)alpha0^4, the tilt n_s=27/28, and the saturation
stop rule are ALREADY canon (cosmology.tex Amplitude section; items 131/128;
item131_{saturation_closure,primordial_tilt_logscale,scalar_mode_projector}.py, all exit 0).
The genuinely new content here is to QUANTIFY the tensor-to-scalar exposure that canon only
flags qualitatively, and to assemble the unified ledger with exact tiers.

exit 0 = the printer chain (event unit + stop rule -> A_s, n_s, H*, e-folds, w0) reproduces the
         measured scalars within their stated pulls; the steady-state readouts are cited; and the
         tensor exposure is correctly COMPUTED and FLAGGED as the open falsifier (exit 0 does NOT
         claim the tensor test passes -- it asserts the exposure is real).
"""
import math

# ---- inputs (compute-don't-recall) ----
ALPHA0 = 1 / 137.0
MPBAR = 2.435e18           # reduced Planck mass, GeV
C_F = 4 / 3               # SU(3) fundamental Casimir (the saturation budget per shell)
U_EVENT = math.log(2)     # one Landauer bit per printed boundary cell
# measured references
A_S_OBS, A_S_ERR = 2.10e-9, 0.03e-9          # Planck 2018 scalar amplitude
NS_OBS, NS_ERR = 0.9649, 0.0042              # Planck 2018 spectral index
R_BOUND = 0.036                              # BICEP/Keck 2021 tensor-to-scalar upper bound

print("[1] THE PRINTER'S TWO CONTROL PARAMETERS:")
N_shell = C_F * ALPHA0 ** -4                  # stop rule: lambda_shell = N_shell*alpha0^4 = C_F
print(f"    event unit  u_event = ln2            = {U_EVENT:.5f}  (one Landauer bit / printed cell)")
print(f"    stop rule   N_shell*alpha0^4 = C_F   = {C_F:.5f}  -> N_shell = {N_shell:.4e}")
assert abs(N_shell * ALPHA0 ** 4 - C_F) < 1e-12

print("\n[2] FOUR SCALAR READOUTS OF THE SATURATED (INFLATION) PHASE -- vs Planck:")
A_nu = ALPHA0 ** 4 / C_F                       # = F_eff*S_j/N_shell with F_eff=S_j=1 -> (3/4)alpha0^4
ns = 27 / 28                                   # 28-clock log-scale generator (item 131)
H_star = math.sqrt(6 * math.pi ** 2 / U_EVENT) * ALPHA0 ** 2 * MPBAR   # from N_shell=S_dS/u_event
n_efolds = 2 * 28                              # Markov relaxation (items 128/131)
pull_A = (A_nu - A_S_OBS) / A_S_ERR
pull_ns = (ns - NS_OBS) / NS_ERR
print(f"    A_s   = (3/4) alpha0^4              = {A_nu:.4e}   ({pull_A:+.2f} sigma vs Planck)")
print(f"    n_s   = 27/28                      = {ns:.5f}     ({pull_ns:+.2f} sigma vs Planck)")
print(f"    N_e   = 2 x 28                     = {n_efolds}          (in the 50-60 window)")
print(f"    H*    = sqrt(6 pi^2/ln2) alpha0^2 Mpbar = {H_star:.3e} GeV  (from N_shell = S_dS/u_event)")
assert abs(pull_A) < 1.5 and abs(pull_ns) < 0.5            # both land within ~1 sigma
assert 50 <= n_efolds <= 60
assert 1.0e15 < H_star < 1.4e15

print("\n[3] STEADY-STATE READOUTS OF THE SAME PRINTER (dark-energy phase, cited):")
w0 = -1 + 1 / 28                               # w(a) = -1 + a/28 (item 131); rho_Lambda item 118
print(f"    w0    = -1 + 1/28                  = {w0:.4f}    (item 131; the printer at steady state)")
print(f"    rho_Lambda = alpha Lambda^4 a0/(4pi L_H) ~ 3.2e-47 GeV^4  (item 118; Landauer exhaust)")
print("    -> inflation (saturated) and dark energy (steady state) are ONE printer (item 128).")
assert abs(w0 + 27 / 28) < 1e-12

print("\n[4] THE DECISIVE EXPOSURE -- the tensor ledger (NEW computation; the open falsifier):")
# standard de Sitter tensor vacuum: P_t = (2/pi^2)(H/Mpbar)^2 ; r = P_t / A_s
P_t = (2 / math.pi ** 2) * (H_star / MPBAR) ** 2
r_naive = P_t / A_S_OBS
suppression_needed = r_naive / R_BOUND
print(f"    de Sitter graviton vacuum at H* gives P_t = (2/pi^2)(H*/Mpbar)^2 = {P_t:.3e}")
print(f"    r_naive = P_t / A_s = {r_naive:.1f}   vs observational bound r < {R_BOUND}")
print(f"    -> the printer's H* OVERPRODUCES tensors by ~{suppression_needed:.0f}x (power) at face value.")
assert r_naive > R_BOUND                       # the exposure is REAL -- this does NOT pass naively
print("    THIS IS THE OPEN FALSIFIER (not asserted to pass): EITHER")
print("      (a) the printer suppresses tensor production by ~%.0fx relative to the dS vacuum" % suppression_needed)
print("          (plausible-but-unproven: the printer is a SCALAR boundary process and gravity is")
print("           emergent, so there may be no spin-2 printing channel), OR")
print("      (b) the event-unit identification u_event = ln2 (hence H*) is wrong.")

print(f"""
[verdict] BOUNDARY PRINTING UNIFIES THE SCALAR SECTOR -- with one sharp exposure.
  ONE printer, two control parameters (event unit u_event=ln2; stop rule N_shell*alpha0^4=C_F),
  reads out FOUR scalar numbers that all land at <=1 sigma: the amplitude A_s=(3/4)alpha0^4
  (+1.0s), the tilt n_s=27/28 (-0.1s), N_e=56 e-folds, and -- as the SAME printer at steady
  state -- w0=-27/28 and rho_Lambda. Inflation and dark energy are one process at two clock
  rates (item 128). That is a genuine tightening of n_s + amplitude + rho_Lambda + the
  expansion ledger onto two parameters.
  TIERS (honest): the amplitude form (F_eff=1, Pi_k fixed) and the tilt are CLOSED; the stop
  rule is CONDITIONALLY closed by the saturation identity lambda_shell=C_F
  (item131_saturation_closure.py) but its C_F=Delta*(8/12) and no-horizon-covariance premises
  remain; the event unit u_event=ln2 is well-motivated (Landauer) but not Locked.
  THE ONE FALSIFIER: the tensor ledger. H*~1.2e15 GeV gives r_naive~{r_naive:.0f}, ~{suppression_needed:.0f}x
  above r<{R_BOUND}. The whole picture stands or falls on whether the (emergent-gravity, scalar)
  printer suppresses tensors -- the decisive open theorem, sharper than any scalar success.
exit 0""")
print("ALL ASSERTIONS PASSED -- scalar chain reproduced; steady-state cited; tensor exposure computed and flagged (not passed).")
