#!/usr/bin/env python3
"""B_N service-packet derivability gate.

Question:
  Can the current QEC service ledger derive a non-oscillator null packet B_N
  that represents omega >> Lambda_QCD quanta without a finer UV lattice?

Required gates:
  G1 momentum: an unaliased generator P_label with spectrum N Lambda n^mu,
      not a compact lattice phase and not a Bloch momentum.
  G2 detector: matter couples linearly to B_N as one detector/LSZ leg, not as
      N insertions of A.
  G3 vacuum: the B_N tower is event-conditioned/non-oscillator, so it adds no
      vacuum zero-point density above Lambda_QCD.

Result:
  Current canon can support G3 and the scalar integer count part of G1, but it
  cannot derive G1 or G2.  The service ledger counts monitored register touches
  and event addresses; it does not define a noncompact spacetime-translation
  generator, and it does not supply a detector coupling replacing j.A.
"""

from __future__ import annotations

import math


LAMBDA_GEV = 0.332
HBARC_GEV_FM = 0.1973269804
A0_FM = HBARC_GEV_FM / LAMBDA_GEV
E_BZ_GEV = math.pi * LAMBDA_GEV


def support_need(E_GeV: float) -> tuple[int, int]:
    return math.ceil(E_GeV / E_BZ_GEV), math.ceil(E_GeV / LAMBDA_GEV)


def gate(name: str, passed: bool, reason: str) -> None:
    status = "PASS" if passed else "FAIL"
    print(f"    {name:<12} {status}: {reason}")


print("[0] Target object")
print("    B_N(n): event-conditioned null service packet")
print("    desired P^mu = N Lambda_QCD (1,n) + residual soft cohomology")
print("    desired detector coupling: one B_N leg, not N Maxwell legs")
print("    desired vacuum status: no oscillator zero-point tower")

print("\n[1] Existing service-ledger ingredients in canon")
ingredients = {
    "integer event count": True,
    "active-address demux/FIFO": True,
    "monitored billing, event-conditioned": True,
    "Gauss-projected IR photon cohomology": True,
    "noncompact translation generator": False,
    "linear detector coupling to B_N": False,
}
for k, v in ingredients.items():
    print(f"    {k:<38} {'present' if v else 'absent'}")

print("\n[2] Three hard gates")
gate(
    "Momentum",
    False,
    "ledger supplies scalar counts/addresses, but no noncompact P_label; compact spatial phases still alias",
)
gate(
    "Detector",
    False,
    "canon detector/QED coupling remains j.A; no one-leg matter vertex for B_N is defined",
)
gate(
    "Vacuum",
    True,
    "event-conditioned Kraus/service counts are not harmonic oscillator modes, so no zero-point tower if B_N exists",
)

print("\n[3] Why scalar service count is insufficient")
for E, label in ((10.0, "10 GeV"), (1.0e3, "1 TeV"), (1.0e5, "100 TeV")):
    n_bz, n_lambda = support_need(E)
    print(
        f"    {label:>7}: service count can name N~{n_lambda}, "
        f"but without P_label this is an energy bill, not a Lorentz momentum eigenstate"
    )

print("\n[4] Consequence")
print(
    "    Canon cannot currently derive B_N from the QEC service ledger.  It can\n"
    "    provide the right kind of non-oscillator bookkeeping (G3), and it can\n"
    "    count events, but it lacks the two load-bearing structures: an unaliased\n"
    "    spacetime momentum label and a linear detector/LSZ coupling to that label.\n"
    "    Therefore B_N remains a proposed new primitive, not a derived object.\n\n"
    "    Minimal positive theorem still needed:\n"
    "      (i) construct P_label on the monitored service history algebra,\n"
    "     (ii) prove [P_label, translations] gives total null four-momentum,\n"
    "    (iii) construct a matter jump operator J_N with <f|J_N|i> depending on\n"
    "          total P^mu and not on N Maxwell insertions,\n"
    "     (iv) prove these histories are source-conditioned so their vacuum\n"
    "          contribution is absent."
)
print("exit 0")
