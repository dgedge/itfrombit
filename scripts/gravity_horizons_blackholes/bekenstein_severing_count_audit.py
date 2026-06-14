#!/usr/bin/env python3
r"""Bekenstein severing-channel count audit.

Item 22 was reformulated as a microscopic rate statement:

    records per horizon node per tick = H0 M_P^2/(16 Lambda^3)
                                      = C alpha_0^2,

where alpha_0^2 has a canon-native home in pair severing: two partners, each
alpha-resolved.  The remaining open number is C ~= 6.87.

This script tests the sharpest code-native closure candidate:

    C_code = 7

from 8-bit entanglement monogamy.  A boundary area element supplies one outward
normal ledger leg; severing that leg from an 8-bit register has exactly seven
non-self monogamy partners.  The monitored-channel theorem then forces a
uniform measure over those seven channels if the partner graph is connected.

The audit also checks the one sharper quotient candidate now suggested by the
strain-decoder/CPT theorem:

    C_code = (8*7 - 1)/8 = 55/8 = 6.875.

Here the numerator is the 56 directed monogamy incidences in an 8-bit register,
with one global-complement/CPT blind slot removed.  That lands to 0.08% with the
bare-alpha constants, but the quotient map is not yet derived.

Exit 0 means: the seven-channel theorem is a real partial closure and the
55/8 quotient is the best sharp target, but exact C=6.87 is not Locked.
"""

from __future__ import annotations

import math

import numpy as np


LAMBDA_QCD_GEV = 0.332
ALPHA0 = 1.0 / 137.0
ALPHA_PHYS = 1.0 / 137.035999084
M_P_GEV = 1.220890e19
H0_KM_S_MPC = 67.4
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25
OMEGA_L = 0.685


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def observed_c(alpha: float = ALPHA0, lam: float = LAMBDA_QCD_GEV) -> tuple[float, float, float]:
    h0_gev = H0_KM_S_MPC / MPC_KM * HBAR_GEV_S
    s_node = M_P_GEV * M_P_GEV / (16.0 * lam * lam)
    n_ticks = lam / h0_gev
    rate = s_node / n_ticks
    c = rate / (alpha * alpha)
    return c, rate, s_node


def monitored_k_graph_stationary(n: int = 7) -> float:
    """Return max stationary error for a monitored complete graph on n channels."""

    a = np.ones((n, n), dtype=float) - np.eye(n)
    gen = a - np.diag(a.sum(axis=1))
    evals = np.sort(np.linalg.eigvalsh(gen))
    check(evals[-2] < -1e-12, f"K_{n} partner graph is connected with a unique fixed point")
    uniform = np.ones(n) / n
    check(np.linalg.norm(gen @ uniform) < 1e-12, f"monitored K_{n} fixed point is uniform")
    # Exact spectral relaxation of a corner start at t=4 is already negligible.
    start = np.zeros(n)
    start[0] = 1.0
    eigvals, eigvecs = np.linalg.eigh(gen)
    prop = eigvecs @ np.diag(np.exp(eigvals * 4.0)) @ eigvecs.T
    return float(np.max(np.abs(prop @ start - uniform)))


def main() -> None:
    print("BEKENSTEIN SEVERING-CHANNEL COUNT AUDIT")

    c_obs, rate_obs, s_node = observed_c()
    print("\n[1] Reproduce the item-22 inversion")
    print(f"  entropy per horizon node       = {s_node:.6e} nats")
    print(f"  required records/node/tick     = {rate_obs:.12e}")
    print(f"  C_obs in units of alpha0^2     = {c_obs:.9f}")
    print(f"  Part-20 identity control        = (3 pi/2)/Omega_L = {1.5 * math.pi / OMEGA_L:.9f}")
    check(6.7 < c_obs < 7.1, "C_obs is the registered 6.87 O(1) factor")

    print("\n[2] Code-native seven-channel theorem candidate")
    n_bits = 8
    unordered_pairs = math.comb(n_bits, 2)
    partner_channels = n_bits - 1
    incidence_per_bit = 2 * unordered_pairs / n_bits
    print(f"  register bits                  = {n_bits}")
    print(f"  unordered monogamy pairs        = C(8,2) = {unordered_pairs}")
    print(f"  pair incidences per bit         = 2*C(8,2)/8 = {incidence_per_bit:.0f}")
    print(f"  one boundary-normal leg sees    = {partner_channels} non-self partners")
    check(partner_channels == 7 and incidence_per_bit == 7, "8-bit monogamy gives seven partner channels per boundary leg")
    err = monitored_k_graph_stationary(7)
    print(f"  monitored K7 relaxation control = max|p - 1/7| at t=4 = {err:.3e}")
    print("  Thus, under the boundary-leg reading, the canon-native count is C_code = 7.")

    print("\n[3] Does C_code=7 close the Bekenstein rate?")
    c_code = 7.0
    rate_code = c_code * ALPHA0 * ALPHA0
    duty_req = c_obs / c_code
    print(f"  predicted rate 7 alpha0^2      = {rate_code:.12e}")
    print(f"  observed/required rate          = {rate_obs:.12e}")
    print(f"  C_code/C_obs - 1                = {(c_code / c_obs - 1.0) * 100:+.3f}%")
    print(f"  residual duty needed            = C_obs/7 = {duty_req:.9f}")
    check(0.97 < duty_req < 0.99, "seven-channel count leaves a real 1-2% residual")

    print("\n[4] Stronger quotient candidate: remove the one global CPT/complement blind slot")
    directed_incidences = n_bits * (n_bits - 1)
    c_quotient = (directed_incidences - 1) / n_bits
    c_obs_physical_alpha, _, _ = observed_c(alpha=ALPHA_PHYS)
    print(f"  directed partner incidences      = 8*7 = {directed_incidences}")
    print("  strain-decoder blind slot        = 1 global complement / CPT class")
    print(f"  quotient candidate               = (56-1)/8 = {c_quotient:.9f}")
    print(f"  miss vs C_obs, bare alpha0        = {(c_quotient / c_obs - 1.0) * 100:+.4f}%")
    print(f"  miss vs C_obs, physical alpha     = {(c_quotient / c_obs_physical_alpha - 1.0) * 100:+.4f}%")
    check(abs(c_quotient / c_obs - 1.0) < 0.001, "55/8 is a sub-per-mille landing against the inversion")
    print("  Status: promising theorem target only.  Canon has the global-complement")
    print("  blind spot, but not yet the map proving it removes exactly one directed")
    print("  Bekenstein severing incidence from the 56-channel ledger.")

    print("\n[5] Alternative normalisations and why they do not close")
    controls = [
        ("one boundary leg: 8-1", 7.0),
        ("directed minus CPT: 55/8", c_quotient),
        ("all unordered bit pairs", float(unordered_pairs)),
        ("Q3 edge checks", 12.0),
        ("six axis-pole modes", 6.0),
        ("CSS X checks", 4.0),
        ("Eg survivors / 8 bits", 2.0 / 8.0),
        ("2 pi", 2.0 * math.pi),
    ]
    print(f"  {'candidate':<28s} {'C':>11s} {'miss vs C_obs':>15s}")
    best = None
    for label, val in controls:
        miss = val / c_obs - 1.0
        print(f"  {label:<28s} {val:11.6f} {miss:14.2%}")
        if best is None or abs(miss) < abs(best[2]):
            best = (label, val, miss)
    assert best is not None
    check(best[0] == "directed minus CPT: 55/8", "the 55/8 quotient is the best code-native count target")

    print("\n[6] What would be required to make seven exact?")
    lam_req = LAMBDA_QCD_GEV * (c_obs / c_code) ** (1.0 / 3.0)
    h_req_ratio = c_code / c_obs
    omega_req = 1.5 * math.pi / c_code
    c_phys_alpha, _, _ = observed_c(alpha=ALPHA_PHYS)
    print(f"  Lambda needed for C=7 at fixed H0,Mp = {lam_req:.6f} GeV ({(lam_req / LAMBDA_QCD_GEV - 1.0) * 100:+.3f}%)")
    print(f"  H0 needed for C=7 at fixed Lambda,Mp  = {h_req_ratio:.6f} x current H0 ({(h_req_ratio - 1.0) * 100:+.3f}%)")
    print(f"  Omega_L needed in Part-20 form        = 3pi/14 = {omega_req:.9f} ({(omega_req / OMEGA_L - 1.0) * 100:+.3f}% vs current)")
    print(f"  using dressed physical alpha instead = C = {c_phys_alpha:.9f} (does not solve it)")
    check(abs(lam_req / LAMBDA_QCD_GEV - 1.0) < 0.01, "C=7 is within a percent-scale chiral-anchor shift")
    check(abs(h_req_ratio - 1.0) > 0.01, "exact closure still needs a new finite-cosmology/chiral correction")

    print(
        """
VERDICT
  Progress: the best channel-count candidate is now a real theorem target, not a
  loose near-integer.  If a horizon area element severs one boundary-normal
  ledger leg, 8-bit monogamy gives exactly seven partner channels, and the
  monitored-channel theorem gives the uniform measure over those channels.

  But this does not fully close item 22.  With canon's pinned constants,

      C_obs  = 6.870...
      C_code = 7

  so the seven-channel theorem overshoots by 1.9%.  The remaining factor
  0.98135 could come from a finite-horizon duty, a slightly shifted chiral
  anchor, or a more detailed severing operator, but none is derived here.

  The sharper candidate is:

      C_code = (8*7 - 1)/8 = 55/8 = 6.875,

  i.e. seven directed partner channels per bit, quotiented by the single global
  complement/CPT blind slot already known in the strain decoder.  This lands at
  sub-per-mille level, but it is not promoted because the quotient map has not
  been derived for the Bekenstein severing ledger.

  Honest status:
      alpha0^2 power      derived/home supplied by pair severing
      O(1) channel count  seven derived under boundary-leg reading
      best exact target   55/8 if the CPT-blind quotient is proven
      exact C=6.87        still an inversion / residual, not Locked

  Promotion criterion:
      derive the missing 0.98135 factor, or derive a different fixed severing
      operator whose channel trace is 6.870... without consuming the observed
      horizon relation.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
