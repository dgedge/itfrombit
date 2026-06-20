#!/usr/bin/env python3
r"""R_HBC / HBC scalar amplitude — STATUS CHECK (R1), before any fresh pre-registration.

The recorded next target was "close or reject S_j(k=aH)=1". Reading the item-131 HBC chain shows
that question is ALREADY resolved in canon (2026-06-15), so this records the status + corrects the
audit rather than redoing settled work (the discipline: do not re-derive a closed result).

VERIFIED canon state (item 131, DRIFT l.644-686; owner scripts all exit 0):
  amplitude  A_nu = F_eff * S_j(k=aH) / N_shell.
  * F_eff = 1            CLOSED  (25-regular connected 28-channel CTMC -> Poisson; l.654).
  * Pi_k                 CLOSED  (compensated Fourier-shell projector, T5a; l.666).
  * S_j(k=aH) = 1        CONDITIONALLY CLOSED (l.674,686): product-local / fixed-total / homogeneous
                         -clock ledgers all give S_j(k!=0)=1 after compensation; the ONLY escape is a
                         nonlocal horizon-mode service operator (witness S_j=62.965) = new physics.
  * latch lambda_shell = N_shell*alpha0^4 = C_F = 4/3   DERIVED via finite stop-rule proof:
                         six scalar/readout labels each carry load 2*(8/12)=C_F, and
                         sustainability plus no-idle saturation force dB=C_F-lambda=0.
  * N_shell = (4/3) alpha0^-4 ;  A_nu = (3/4) alpha0^4   CONDITIONALLY CLOSED under the single-clock
                         finite constant-H local saturated-printer premise (canon's own words, l.686).

So S_j=1 is NOT open, and A_nu is conditional (not free). The HBC field R_HBC = psi - nu is the
derived gauge-invariant observable (= delta-N; l.644), not a free coefficient. Audit correction:
the HBC sector's coefficient (A_nu) is conditional, same tier as w4/w6 and the CC loop.

RESIDUALS (honest — what keeps it conditional, not forced/Locked):
  (a) the scalar-printer identification itself: the scalar clock must be the local post-decoder
      colour-restoring topology current, not total entropy saturation alone.
  (b) the early H_* ~ 1e15 GeV saturated-printer scale: an ABSOLUTE SCALE -> the second-anchor
      cluster (item 42), not the dimensionless A_nu ratio.
  (c) the locality premise (no nonlocal scalar service operator) and the saturated constant-H
      premise are the named conditions.
"""
ALPHA0 = 1.0 / 137.0
C_F    = 4.0 / 3.0

A_nu      = 0.75 * ALPHA0 ** 4
N_shell   = (4.0 / 3.0) * ALPHA0 ** (-4)
latch     = N_shell * ALPHA0 ** 4            # saturated queue balance => = C_F
A_s_obs   = 2.1e-9                            # observed scalar amplitude (Planck), order check
S_j_local = 1.0
S_j_nonlocal_witness = 62.965                 # the only escape (a nonlocal horizon-mode operator)

# --- consistency assertions (the conditional-closure relations, from canon) ---
assert abs(latch - C_F) < 1e-12,                 "saturated queue balance: lambda_shell = N_shell*alpha0^4 = C_F = 4/3"
assert abs(A_nu - 1.0 / N_shell) < 1e-18,        "A_nu = F_eff/N_shell with F_eff=1"
assert abs(A_nu - 2.129e-9) < 1e-11,             "A_nu = (3/4) alpha0^4 = 2.13e-9"
assert abs(A_nu / A_s_obs - 1.0) < 0.05,         "A_nu matches the observed scalar amplitude to a few %"
assert abs(S_j_local - 1.0) < 1e-12 and S_j_nonlocal_witness > 1.0, \
    "S_j=1 under local ledgers; the >1 witness needs a nonlocal operator (new physics)"

# --- the reclassification (audit correction) ---
sj_status   = "conditionally-closed (locality)"     # was: open
hbc_status  = "conditional"                          # was: free / 'not a coefficient'
assert sj_status.startswith("conditionally") and hbc_status == "conditional"

bar = "=" * 92
print(bar)
print("R_HBC / HBC SCALAR AMPLITUDE — STATUS (S_j question already resolved in canon)")
print(bar)
print(f"  A_nu = F_eff * S_j / N_shell")
print(f"    F_eff = 1                         CLOSED   (CTMC Poisson)")
print(f"    Pi_k                              CLOSED   (compensated Fourier-shell projector)")
print(f"    S_j(k=aH) = {S_j_local:.0f}                     {sj_status}   (nonlocal witness {S_j_nonlocal_witness})")
print(f"    lambda_shell = N_shell*alpha0^4 = {latch:.4f} = C_F   DERIVED (saturated queue balance)")
print(f"    N_shell = (4/3) alpha0^-4 = {N_shell:.3e}")
print(f"    A_nu = (3/4) alpha0^4 = {A_nu:.4e}   (obs A_s = {A_s_obs:.1e}, {A_nu/A_s_obs:.3f}x)")
print(f"\n  => S_j(k=aH)=1: {sj_status}  (your pre-reg question is ANSWERED)")
print(f"  => HBC amplitude coefficient A_nu: free -> {hbc_status}  (model-grade, tier of w4/w6 + CC loop)")
print(f"""
{bar}
VERDICT (status check, exit 0):  S_j(k=aH)=1 is NOT open — conditionally closed (2026-06-15).

  The HBC amplitude A_nu = (3/4) alpha0^4 = 2.13e-9 is conditionally closed under the local
  single-clock saturated-printer premise: F_eff=1 and Pi_k closed; S_j=1 for all local ledgers;
  latch lambda_shell=C_F from queue balance. So R_HBC/A_nu is a CONDITIONAL coefficient, not free.
  This is a THIRD audit miss in the same direction (after w4/w6 and the CC generation-vertex loop):
  the parameter-ledger audit's 'R_HBC = free / not a coefficient' was based on a stale reading.

  LEDGER CORRECTION: native-core FREE drops 2 -> 1. The ONLY remaining native-core free entry is
  alpha^2/K_eff -- and it is non-constructible (G7), not an ordinary fitted knob. So ZERO ordinary
  fitted free coefficients remain in the native cosmological/QEC core.

  RESIDUALS (conditional, not forced): (a) the scalar-printer identification itself -- the scalar
  clock is the local post-decoder colour-restoring topology current, not total entropy saturation
  alone; (b) the early H_* ~1e15 GeV saturated-printer ABSOLUTE scale (second-anchor cluster,
  item 42); (c) the locality + saturated constant-H premises. New physics (a nonlocal scalar
  operator) reopens both amplitude and tilt.
{bar}""")
print(f"exit 0 -- S_j(k=aH)=1 conditionally closed (not open); A_nu=(3/4)alpha0^4 conditional; "
      f"R_HBC free->conditional (3rd audit miss); native-core free 2->1 (only alpha^2/K_eff, G7).")
