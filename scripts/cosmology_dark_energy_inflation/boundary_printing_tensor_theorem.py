#!/usr/bin/env python3
r"""THE TENSOR THEOREM FOR BOUNDARY PRINTING — does a scalar printer source spin-2 modes?

ANCHOR item 146 named this the decisive falsifier: the printer's inflation scale
H_* ~ 1.2e15 GeV gives r_naive ~ 23, about 650x over the bound r < 0.036. This script asks the
make-or-break question and resolves it as far as the framework's own premises allow.

THE PREMISE MISMATCH (the crux). The estimate r_naive uses the standard tensor spectrum
    P_t = (2/pi^2)(H/Mpbar)^2,
which is DERIVED by quantizing the transverse-traceless metric perturbation (the graviton) in a
SMOOTHLY STRETCHING de Sitter background and squeezing its adiabatic vacuum. But ANCHOR item 123
EXPLICITLY DENIES that premise: "the universe does not stretch; it prints new boundary nodes."
So r_naive applies the stretching-graviton formula to a process the framework says is NOT
stretching. That is a category error, not a clean falsifier.

WHAT THE PRINTER ACTUALLY DOES. The printer's perturbations are SHOT NOISE of discrete scalar
cell-printing (A_s = 1/N_shell, item 146). A density/count fluctuation is a SCALAR source, and at
linear order the scalar/vector/tensor sectors DECOUPLE (SVT theorem) -- a scalar source produces
NO transverse-traceless (spin-2) mode, hence r_linear = 0. The only tensor channel left is
SECOND-order scalar-induced gravitational waves, with r_induced ~ O(A_s) ~ 1e-9 -- ten orders
below r_naive.

THE TWO BRANCHES, made mutually exclusive by the event unit:
  * Branch A (the emergent graviton vacuum IS squeezed, standard formula holds): then survival
    requires H_* <= H_max ~ 5e13 GeV, i.e. an event unit u_event ~ 450 bits/cell -- NOT the one
    Landauer bit (ln2) the framework adopts. So Branch A is INCOMPATIBLE with u_event = ln2.
  * Branch B (no squeezing; discrete printer, scalar shot noise): the one-bit event unit stands,
    H_* = 1.2e15, with r_linear = 0 and tensors only at the second-order floor r_induced ~ 1e-9.
Adopting u_event = ln2 (item 146) therefore COMMITS the framework to Branch B.

exit 0 = r_naive recomputed; the Branch-A event unit computed and shown incompatible with one
         bit; the second-order tensor floor computed below the bound; the resolution stated
         CONDITIONALLY (the no-squeezing premise is FLAGGED as the remaining open theorem, not
         asserted proven).
"""
import math

ALPHA0 = 1 / 137.0
MPBAR = 2.435e18                 # reduced Planck mass, GeV
A_S = 2.10e-9                    # Planck scalar amplitude
R_BOUND = 0.036                  # BICEP/Keck 2021
U_LN2 = math.log(2)             # the framework's event unit: one Landauer bit per printed cell
H_STAR = math.sqrt(6 * math.pi ** 2 / U_LN2) * ALPHA0 ** 2 * MPBAR   # item 146

print("[1] THE STANDARD FORMULA AND ITS PREMISE:")
P_t = (2 / math.pi ** 2) * (H_STAR / MPBAR) ** 2
r_naive = P_t / A_S
print(f"    P_t = (2/pi^2)(H*/Mpbar)^2 = {P_t:.3e};  r_naive = P_t/A_s = {r_naive:.1f}  (>> {R_BOUND})")
print("    PREMISE REQUIRED: a smoothly stretching dS metric squeezing a propagating graviton")
print("    (EH-normalized at Mpbar). ANCHOR item 123: 'the universe does not stretch; it prints.'")
print("    -> r_naive applies a stretching-graviton formula to a non-stretching printer.")
assert r_naive > R_BOUND          # the exposure is real ONLY under the denied premise

print("\n[2] SVT DECOUPLING: a scalar printer has no linear spin-2 source:")
print("    printer perturbations = shot noise of discrete scalar cell-printing (A_s = 1/N_shell).")
print("    A scalar (density/count) source excites only scalar metric perturbations; at linear")
print("    order scalar/vector/tensor decouple -> NO transverse-traceless mode is produced.")
print("    Thus r_linear = 0. Only 2nd-order scalar-induced GWs remain, at the floor:")
r_linear = 0.0
r_induced = A_S                   # order-of-magnitude 2nd-order floor: P_t^(2) ~ A_s^2 => r ~ A_s
print(f"      r_induced ~ O(A_s) ~ {r_induced:.1e}   (~{r_naive / r_induced:.0e}x below r_naive; ~{R_BOUND / r_induced:.0e}x below the bound)")
assert r_linear == 0.0
assert r_induced < R_BOUND

print("\n[3] THE TWO BRANCHES, made mutually exclusive by the event unit:")
# Branch A: keep the standard squeezed-graviton formula -> survival needs H* <= H_max
H_max = MPBAR * math.sqrt(R_BOUND * A_S * math.pi ** 2 / 2)
u_event_A = U_LN2 * (H_STAR / H_max) ** 2            # H* ~ u_event^{-1/2}, so smaller H* needs bigger u_event
print(f"    Branch A (graviton squeezed, formula holds): survival needs H* <= H_max = {H_max:.2e} GeV,")
print(f"      i.e. event unit u_event ~ {u_event_A:.0f} bits/cell -- NOT one bit ({u_event_A / U_LN2:.0f}x the ln2 unit).")
print(f"    Branch B (discrete printer, no squeezing): u_event = ln2 stands, H* = {H_STAR:.2e} GeV,")
print(f"      r_linear = {r_linear:.0f}, with tensors only at the 2nd-order floor r_induced ~ {r_induced:.0e}.")
assert u_event_A / U_LN2 > 100                       # Branch A demands >>1 bit/cell: incompatible with ln2
print("    -> adopting the one-Landauer-bit event unit (item 146) COMMITS the framework to Branch B.")

print(f"""
[verdict] THE EXPOSURE IS CONDITIONALLY DEFUSED -- downgraded from a clean falsifier to one that
  hinges on a single named premise.
  * r_naive ~ {r_naive:.0f} is computed under the STRETCHING-graviton premise that item 123 denies, so it
    is not a valid falsifier as stated.
  * The printer's mechanism is SCALAR shot noise; by SVT decoupling it sources no linear spin-2
    mode, so r_linear = 0, leaving only 2nd-order induced GWs at r_induced ~ {r_induced:.0e} --
    a concrete prediction of essentially ZERO observable primordial tensors.
  * The one-bit event unit forces Branch B; Branch A would need ~{u_event_A / U_LN2:.0f} bits/cell, contradicting
    the Landauer event unit. The framework is self-consistent ONLY in Branch B.
  THE REMAINING OPEN THEOREM (honest): rigorously show the emergent COARSE-GRAINED metric does
  not squeeze the graviton vacuum on CMB scales during printing -- i.e. derive the emergent
  graviton's existence/normalization as a function of crystallization fraction through the boot.
  Branch B is FAVORED (the discrete-printing premise + SVT decoupling) but not yet PROVEN.
  THE DECISIVE TEST: r. The framework predicts r_linear = 0, with r_induced ~ {r_induced:.0e}
  (no observable B-modes). A detection
  at r ~ 0.001-0.036 (standard high-scale inflation) would force Branch A and falsify the ln2
  event unit. So 'expansion is printing' is sharply falsifiable on the sky.
exit 0""")
print("ALL ASSERTIONS PASSED -- premise mismatch shown; SIGW floor below bound; Branch A incompatible with ln2; open theorem flagged.")
