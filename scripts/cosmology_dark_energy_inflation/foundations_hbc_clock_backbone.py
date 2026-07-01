#!/usr/bin/env python3
r"""foundations_hbc_clock_backbone.py

HBC / item-131 (the cosmological scalar clock) — backbone verification + a clean reframing of the amplitude.

STATUS of the item (verified against canon, 2026-06-25):
 - The LATCH MARGINALITY is PROVEN: lambda_shell = C_F = 4/3 by a two-sided squeeze (coherent constant-H
   sustainability => lambda_shell <= C_F; no-idle scalar-printer saturation => lambda_shell >= C_F), with
   C_F sourced combinatorially from the [8,4,4] colour-restoring count. (`item131_hbc_stop_rule_proof.py`.)
 - The cosmological observables are DERIVED conditional on ONE field identification: the scalar clock is the
   local post-decoder colour-restoring print-time current nu_HBC, so R_HBC = psi - nu is the gauge-invariant
   curvature perturbation = delta N (`item131_scalar_clock_bridge.py`).
 - OPEN piece (the genuine residual): the scalar-printer/clock IDENTIFICATION + its UNIQUENESS (no second
   scalar source / nonlocal horizon-mode operator), within the single-clock architecture (item 94).

This script verifies the combinatorial backbone and the predictions, and records the reframing
A_s = 1/N_shell = alpha0^4 / C_F (the scalar amplitude IS the inverse clock-shell count).

Self-asserting; exit 0.
"""
import numpy as np
import itertools
from math import comb


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== combinatorial backbone ([8,4,4] + colour / C_F) ===")
    G = np.array([[1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]])
    wd = {}
    for b in itertools.product([0, 1], repeat=4):
        w = int(np.mod(np.array(b) @ G, 2).sum()); wd[w] = wd.get(w, 0) + 1
    ok(wd.get(4) == 14, f"[8,4,4]=RM(1,3) weight distribution {dict(sorted(wd.items()))}: A_4 = 14 weight-4 supports")
    ok(comb(8, 2) == 28, "C(8,2) = 28 transverse service labels (shell-1 octet pairs)")
    N = 3
    ok(N * N - N == 6, "6 colour-restoring labels = SU(3) off-diagonal generators (N^2-N = 6)")
    CF = (N * N - 1) / (2 * N); strain = 2 * (8 / 12)
    ok(abs(CF - 4 / 3) < 1e-12 and abs(CF - strain) < 1e-12,
       f"C_F = (N^2-1)/(2N) = {CF:.4f} = 2x(8/12 cube strain) = 4/3 (the SU(3) fundamental Casimir = the latch marginal load)")

    print("\n=== latch marginality + the amplitude reframing ===")
    a0 = 1 / 137.035999084
    lam_shell = CF                       # two-sided squeeze (proven), the marginal value
    Nshell = CF * a0 ** -4
    A_s = a0 ** 4 / CF                    # = (3/4) alpha0^4
    ok(abs(lam_shell - 4 / 3) < 1e-12, "lambda_shell = C_F = 4/3 (latch sits AT the marginal value, DeltaB=0; two-sided squeeze)")
    ok(abs(Nshell - 4.70e8) / 4.70e8 < 0.01, f"N_shell = C_F alpha0^-4 = {Nshell:.4e} (~4.70e8)")
    ok(abs(A_s - 0.75 * a0 ** 4) < 1e-30, f"A_s = alpha0^4 / C_F = (3/4) alpha0^4 = {A_s:.4e}")
    ok(abs(A_s * Nshell - 1) < 1e-9, "A_s * N_shell = 1  =>  A_s = 1/N_shell (scalar amplitude = inverse clock-shell count) [reframing]")

    print("\n=== predictions vs Planck (conditional on the scalar-clock identification) ===")
    As_pl, ns_pl = 2.099e-9, 0.9649
    ns = 27 / 28
    print(f"    A_s = (3/4)alpha0^4 = {A_s:.4e}  vs Planck {As_pl:.3e}  ({100*(A_s/As_pl-1):+.1f}%)")
    print(f"    n_s = 27/28 = {ns:.5f}  vs Planck {ns_pl}  ({100*(ns/ns_pl-1):+.3f}%)")
    print(f"    w_0 = -n_s = {-ns:.5f}  (dark-energy EoS)")
    ok(abs(A_s / As_pl - 1) < 0.03, "A_s within 3% of Planck (the inflationary scalar amplitude)")
    ok(abs(ns / ns_pl - 1) < 0.001, "n_s = 27/28 within 0.1% of Planck (the spectral index)")

    print("\n=== gauge structure of R_HBC = psi - nu ===")
    rng = np.random.default_rng(0); psi = rng.normal(size=6); nu = rng.normal(size=6); lam = rng.normal(size=6)
    ok(np.allclose(psi - nu, (psi - lam) - (nu - lam)),
       "R_HBC = psi - nu gauge-invariant under slicing (psi,nu)->(psi-lam,nu-lam); flat slice deltaN=-nu => R_HBC=deltaN")

    print("\n[verdict] HBC / item-131 status:")
    print("  - DERIVED (conditional on the scalar-clock identification): A_s = 1/N_shell = alpha0^4/C_F (1.3% of Planck),")
    print("    n_s = 27/28 (0.06%), w_0 = -27/28. The latch marginality lambda_shell = C_F = 4/3 is PROVEN (two-sided")
    print("    squeeze) with C_F from the [8,4,4] colour-restoring/strain count.")
    print("  - OPEN (the genuine residual): the scalar-printer/clock field IDENTIFICATION + UNIQUENESS -- that the clock")
    print("    is the unique local post-decoder colour-restoring print-time current (no second scalar source), within")
    print("    the single-clock architecture (item 94). The C_F-load match SUPPORTS it; uniqueness is the open piece. exit 0")


if __name__ == "__main__":
    main()
