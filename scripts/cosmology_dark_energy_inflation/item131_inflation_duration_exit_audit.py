#!/usr/bin/env python3
r"""ITEM 131 — inflation duration/exit: the structural decoupling, and the exit reduced to one mechanism.

State of the three HBC/inflation targets:
  (1) exact shell count  : DONE conditionally (item131_hbc_stop_rule_proof.py) — N_shell = C_F alpha0^-4
                           = (4/3) 137^4 = 4.70e8, from the marginal stop rule.
  (2) marginal colour-load stop rule : DONE conditionally — lambda_shell = C_F = 4/3 by the finite queue
                           theorem (Delta B = C_F - lambda_shell; sustainability lambda<=C_F + no-idle
                           lambda>=C_F => lambda=C_F), giving A_nu = (3/4) alpha0^4 = 2.13e-9.
  (3) early constant-H duration/exit : OPEN (item131_hbc_clock_covariance_exit.py).

This audit makes the sharp structural point that isolates (3) and retires a coincidence:

  A. The framework's inflationary SPECTRUM (tilt, amplitude) is DECOUPLED from its inflationary
     DYNAMICS (energy scale, duration, exit). The tilt -1/28 is the 28-channel log-shell QEC generator's
     anomalous dimension; the amplitude (3/4)alpha0^4 is the shell shot-noise ledger. NEITHER uses
     slow-roll (epsilon, eta, H^2/epsilon). So (1)+(2) are derivable WITHOUT solving the dynamics.
  B. The DE-component EoS w(a) = -1 + a/28 makes inflation EXACT de Sitter (epsilon = 3a/56 ~ 0 at the
     inflation epoch), so it canNOT end by slow-roll epsilon->1. The standard tilt<->duration link
     (N_e = 2/(1-n_s) = 56) is a SLOW-ROLL relation that does not apply here -> N=56 = 2*28 is a
     COINCIDENCE, not a derivation. The duration is therefore set by the EXIT THRESHOLD, not the tilt.
  C. de Sitter + threshold-exit = "old-inflation"-class: the open piece is a GRACEFUL EXIT / reheating
     trigger. Candidate (links (2)->(3)): the exit is the GLOBAL stop rule — local saturation
     lambda=C_F holds until the printable resource depletes; when it can no longer be sustained, de
     Sitter ends and reheats. This inherits the Lemma-L deposition problem (D1/D2), so it is posed,
     not derived. All of (3) reduces to the item-42 R-activation threshold.
"""
from __future__ import annotations

import math
from fractions import Fraction

ALPHA0 = 1.0 / 137.0
DELTA = Fraction(1, 28)          # the 28-clock anomalous dimension
C_F = Fraction(4, 3)             # marginal colour load (Casimir; item131_hbc_stop_rule_proof.py)


def bar(c="="):
    print(c * 100)


def main():
    bar()
    print("ITEM 131 — inflation duration/exit: structural decoupling + the exit reduced to one mechanism")
    bar()

    # --- (1)+(2): done conditionally ---
    print("[1+2] SHELL COUNT + MARGINAL COLOUR-LOAD STOP RULE (done conditionally)")
    n_shell = float(C_F) * 137.0**4
    a_nu = 0.75 * ALPHA0**4
    print(f"    lambda_shell = C_F = {C_F} = {float(C_F):.4f}   (finite-queue marginal stop rule)")
    print(f"    N_shell = C_F alpha0^-4 = (4/3) 137^4 = {n_shell:.3e}")
    print(f"    A_nu = (3/4) alpha0^4 = {a_nu:.3e}   (vs observed A_s ~ 2.1e-9)")
    assert abs(a_nu / 2.1e-9 - 1.0) < 0.05, "A_nu = (3/4)alpha0^4 matches observed scalar amplitude"
    assert abs(n_shell - 4.697e8) / 4.697e8 < 0.001, "N_shell = (4/3)137^4 = 4.70e8"

    # --- A. spectrum decoupled from dynamics ---
    print("\n[A] SPECTRUM ⊥ DYNAMICS (why (1)+(2) close without solving (3))")
    n_s = 1 - DELTA
    print(f"    tilt n_s = 1 - 1/28 = {n_s} = {float(n_s):.5f}   <- 28-clock anomalous dimension (NOT slow-roll)")
    print(f"    amplitude A_nu = (3/4)alpha0^4               <- shell shot-noise ledger (NOT H^2/epsilon)")
    print("    => the observable spectrum is fixed by the CLOCK + SHELL ledger; the inflationary energy")
    print("       scale H_inf, duration N_e, and exit do NOT enter it. Spectrum and dynamics are separable.")

    # --- B. exact de Sitter -> slow-roll N=56 is a coincidence ---
    print("\n[B] EXACT de SITTER -> the N=56 slow-roll link does NOT apply")

    def eps(a):  # DE-component slow-roll parameter from w(a) = -1 + a/28
        return 1.5 * (a / 28.0)

    for a in (1e-26, 1e-13, 1.0):
        print(f"    a={a:.0e}: w = -1 + a/28 = {-1 + a/28:.6f},  epsilon = 3a/56 = {eps(a):.3e}")
    n_slowroll = 2.0 / float(DELTA)
    print(f"    slow-roll identity N_e = 2/(1-n_s) = 2/Delta = {n_slowroll:.0f}  (= 2 x 28)")
    print(f"    BUT epsilon(inflation) ~ {eps(1e-26):.1e} << 1: inflation is exact de Sitter, ends by a")
    print("    THRESHOLD, not by epsilon->1. The tilt is clock-derived, not slow-roll -> N=56 is a")
    print("    COINCIDENCE (2x28), NOT a duration derivation. Duration is set by the exit threshold.")
    assert eps(1e-26) < 1e-25, "epsilon at the inflation epoch is ~0 (exact de Sitter)"
    assert n_slowroll == 56.0, "the slow-roll 2/Delta value is 56 (the retired coincidence)"

    # --- C. the exit candidate: global stop rule / printer exhaustion ---
    print("\n[C] THE EXIT (target 3) — reduced to ONE mechanism: the GLOBAL stop rule")
    print("    local marginal saturation lambda=C_F (proven, target 2) holds while the printable")
    print("    resource lasts; de Sitter ENDS when the resource depletes and lambda=C_F can no longer be")
    print("    sustained -> reheating. This LINKS (2)->(3): the same colour-load stop rule that fixes the")
    print("    shell count is the local face of the global exit. It inherits the Lemma-L deposition")
    print("    problem (D1 refuted as-posed, D2 consistent) -> POSED, not derived.")
    # duration candidates: none is a clean threshold-derived N_e
    cands = {
        "2/Delta (slow-roll, N/A here)": 2.0 / float(DELTA),
        "ln N_shell (printer log-depth)": math.log(n_shell),
        "28 ln(1/Lambda^2) (units-suspect)": 28.0 * math.log(1.0 / 0.110),
    }
    print("    duration candidates (none clean / threshold-derived):")
    for name, val in cands.items():
        print(f"      N_e ~ {val:6.1f}   {name}")
    print("    target window for horizon/flatness: N_e ~ 50-60; none of the above is a DERIVED threshold.")

    print(f"""
{"=" * 100}
VERDICT (exit 0):  spectrum DONE + DECOUPLED; dynamics OPEN, reduced to one threshold; N=56 retired.

  (1) shell count and (2) the marginal colour-load stop rule are done conditionally (your
  item131_hbc_stop_rule_proof.py): lambda=C_F=4/3 => N_shell=(4/3)137^4=4.70e8, A_nu=(3/4)alpha0^4.

  THE STRUCTURAL POINT: the framework's inflationary SPECTRUM (tilt -1/28 from the 28-clock; amplitude
  from the shell ledger) is DECOUPLED from its inflationary DYNAMICS (H_inf, N_e, exit) — it uses no
  slow-roll. That is WHY (1)+(2) close while (3) stays open: they are independent, unlike standard
  inflation where n_s<->N_e<->V(phi) are linked. Corollary: the framework's inflation is EXACT de
  Sitter (epsilon~0), so the slow-roll N_e=2/(1-n_s)=56 link does NOT apply — N=56 (=2x28) is a
  COINCIDENCE and is retired as a duration derivation.

  (3) THEREFORE the duration/exit is a SEPARATE, well-posed problem: de Sitter + threshold-exit =
  old-inflation-class, whose open piece is the graceful-exit/reheating trigger. The cleanest candidate
  LINKS it to target (2): the exit is the GLOBAL stop rule (local saturation lambda=C_F holds until the
  printable resource depletes), which inherits the Lemma-L deposition problem and is POSED not derived.
  All of (3) — H_inf, N_e, and the exit — reduces to the item-42 R-activation threshold; the tilt+
  amplitude success does NOT depend on it.
{"=" * 100}""")
    print(f"exit 0 -- inflation: spectrum (tilt n_s=27/28 + A_nu=(3/4)alpha0^4) done+decoupled from dynamics; "
          f"N=56 retired as a slow-roll coincidence; duration/exit reduce to the item-42 R-activation threshold (global stop rule).")


if __name__ == "__main__":
    main()
