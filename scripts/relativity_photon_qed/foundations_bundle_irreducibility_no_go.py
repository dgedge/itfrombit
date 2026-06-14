#!/usr/bin/env python3
"""No-go for making a trans-Lambda null bundle into one QED photon by Wilson mechanics alone.

Previous result:
  foundations_trans_lambda_quanta.py showed that, under the single-scale
  a0 = hbar c / Lambda_QCD plus CC premise, the only kinematically viable
  representation of omega >> Lambda_QCD quanta is a collinear null bundle:
  many sub-cutoff lightlike records with total P^2 = 0.

Question here:
  Can an ordinary gauge-invariant Wilson/service-current construction make that
  bundle detector-irreducible, so QED amplitudes see only total P^mu and not N
  independent soft photons?

Answer:
  No, not under the existing Wilson/QED rules.  Gauge invariance can make a
  composite operator, but the S-matrix still sees N external photon legs.  The
  standard j.A vertex is linear in A, so one local detector vertex annihilates
  one photon.  Absorbing N sub-cutoff records as one observed high-energy event
  requires an N-photon operator with N independent Ward identities and N powers
  of the gauge coupling.  That is multi-photon QED, not a single photon.

Consequence:
  The null-bundle route is kinematically possible but dynamically open.  A
  positive closure would require a new canon-derived service-current primitive
  B_N whose LSZ/detector reduction is linear in B_N while B_N carries total
  P^mu and no independent UV oscillator density.  That object is not present
  in the current Wilson/Gauss photon construction.
"""

from __future__ import annotations

import math


ALPHA0 = 1.0 / 137.0
E_CHARGE = math.sqrt(4.0 * math.pi * ALPHA0)
LAMBDA_GEV = 0.332
E_BZ_GEV = math.pi * LAMBDA_GEV


def n_bz(E_GeV: float) -> int:
    return math.ceil(E_GeV / E_BZ_GEV)


print("[0] Setup")
print(f"    alpha0=1/137, e=sqrt(4 pi alpha0)={E_CHARGE:.6f}")
print(f"    Lambda_QCD={LAMBDA_GEV:.3f} GeV, BZ-edge energy=pi Lambda={E_BZ_GEV:.3f} GeV")

print("\n[1] External-leg counting")
for E, label in ((10.0, "10 GeV"), (1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
    n = n_bz(E)
    print(f"    {label:>7}: minimum BZ-edge records N={n}")
assert n_bz(1.0e3) == 959

print("\n[2] Linear detector vertex cannot see the bundle as one photon")
print("    Standard QED/Wilson coupling: H_int = integral j_mu A^mu.")
print("    Matrix element <detector excited; N-1 photons | H_int | detector ground; N photons>")
print("    is nonzero: one photon is absorbed.")
print("    Matrix element <detector excited; 0 photons | H_int | detector ground; N photons>")
print("    is zero for N>1: one local A field cannot annihilate N external photons.")
for n in (1, 2, 3, n_bz(10.0), n_bz(1.0e3)):
    first_order_all_absorb = n == 1
    print(f"      N={n:4d}: first-order all-bundle absorption possible? {first_order_all_absorb}")
assert not (n_bz(1.0e3) == 1)

print("\n[3] N-photon absorption is a different amplitude, not total-P single-photon QED")
for E, label in ((10.0, "10 GeV"), (1.0e3, "1 TeV")):
    n = n_bz(E)
    log_coupling_penalty = (n - 1) * math.log10(E_CHARGE)
    print(
        f"    {label:>7}: N={n}, coupling-order relative to one photon = e^(N-1), "
        f"log10 = {log_coupling_penalty:.1f}"
    )
print("    The scaling depends explicitly on N.  A single high-energy photon amplitude is O(e)")
print("    with one transverse polarisation vector at total P, not an N-leg tensor.")

print("\n[4] Ward identities also see the legs")
print("    A true N-photon external state has N independent gauge redundancies:")
print("      epsilon_i -> epsilon_i + c_i k_i,  i=1..N.")
print("    A one-photon state at total P has only one Ward identity.")
print("    Collapsing N Ward identities to one loses N-1 constraints unless a new")
print("    gauge-invariant composite field B_N is introduced and proven to be the detector variable.")
for n in (2, 10, n_bz(1.0e3)):
    print(f"      N={n:4d}: extra Ward constraints relative to one photon = {n - 1}")
assert n_bz(1.0e3) - 1 > 0

print("\n[5] Wilson operators do not evade this")
print("    Closed Wilson loops and service-current exponentials are gauge-invariant, but")
print("    exp(i integral A) expands into 0,1,2,... photon components.  Selecting the")
print("    N-photon component gives an N-leg operator; leaving the exponential intact gives")
print("    a coherent/Poisson object, not a fixed-energy one-photon external state.")

print("\n[6] Verdict")
print(
    "    The requested irreducible bundle theorem cannot be proved from the current\n"
    "    Wilson/Gauss/QED mechanics.  The opposite is derivable: a collinear bundle\n"
    "    has null total kinematics, but ordinary detectors and QED amplitudes see\n"
    "    N external photon legs, not one total-P photon.  Therefore the single-scale\n"
    "    framework still lacks a trans-Lambda photon representation.\n\n"
    "    A positive closure would need a new canon-derived object B_N with three\n"
    "    properties: (i) gauge-invariant composite of sub-cutoff records, (ii) one\n"
    "    detector/LSZ leg carrying total P^mu, and (iii) no independent UV oscillator\n"
    "    density that reopens the CC.  That is a new service-current primitive, not\n"
    "    a consequence of ordinary Wilson loops."
)
print("exit 0")
