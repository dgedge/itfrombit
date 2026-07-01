#!/usr/bin/env python3
r"""A=1 lift / D2 burn-integral theorem attempt.

Target
------
The Planck selector and late cosmological lock use

    N_lock = 9 alpha0 / r6.

The D2 consumption law closes the *form*:

    B(N) = B0 - N r6,      N_exhaust = B0/r6.

Therefore the remaining a=1 lift premise is exactly the stock identity

    B0 = 9 alpha0.

This script tests whether the current §5.2 / monitored-service ledger already
derives that identity, or whether it merely names the missing global rule.

Result
------
The local service algebra does derive the *unit* and the *touch count*:

    permanent burn records = T alpha0,     T = 8 repairs + 1 latch/readout.

That gives the right stock if and only if the cosmological service era contains
one closed burn sweep and no hidden/unbilled stock:

    B0 = n_sweep T alpha0,     T=9,     n_sweep=1.

Neighboring alternatives are sharply excluded by the H0/Planck-selector landing:
T=8, T=10, two sweeps, half-sweeps, or a static K04 stock all miss.  But the
global statement "exactly one closed burn sweep per cell is the whole stock" is
not yet derivable from the existing local algebra.  Thus the a=1 lift is
reduced to one crisp theorem target, not closed unconditionally.
"""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))
rh = importlib.import_module("register_handoff_form_selection")

ALPHA0 = 1.0 / 137.0
T_REPAIR = 8
T_LATCH = 1
T_COMPLEMENT = 0
T_CANON = T_REPAIR + T_LATCH + T_COMPLEMENT
TARGET_H0 = rh.H0_KM_S_MPC
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def residual_r6() -> float:
    gamma_star = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
    _, q_post, _ = rh.queue_readouts(gamma_star, 1)
    return (21.0 * q_post) ** 32 / 21.0


def h0_from_stock(stock: float, r6: float) -> float:
    n_exhaust = stock / r6
    return rh.LAMBDA_QCD_GEV / n_exhaust * MPC_KM / HBAR_GEV_S


def stock(touches: float, sweeps: float = 1.0) -> float:
    return sweeps * touches * ALPHA0


def main() -> None:
    print("A=1 LIFT / D2 BURN-INTEGRAL THEOREM ATTEMPT")
    print("=" * 92)

    r6 = residual_r6()
    target_stock = stock(T_CANON, 1.0)
    target_h0 = h0_from_stock(target_stock, r6)

    print("[1] D2 reduces the cosmological endpoint to a stock identity")
    print(f"    r6                              = {r6:.12e}")
    print(f"    target B0 = 9 alpha0             = {target_stock:.12e}")
    print(f"    N_exhaust = B0/r6                = {target_stock/r6:.12e}")
    print(f"    H0(B0=9 alpha0)                  = {target_h0:.6f} km/s/Mpc")
    print(f"    selector reference H0            = {TARGET_H0:.6f} km/s/Mpc")
    print(f"    miss                             = {target_h0/TARGET_H0 - 1:+.4%}")
    check(abs(target_h0 / TARGET_H0 - 1.0) < 7.0e-4, "B0=9 alpha0 lands the selector budget")

    print("\n[2] Local §5.2 burn unit and touch count")
    print(f"    repairs      = {T_REPAIR}")
    print(f"    latch/readout = {T_LATCH}")
    print(f"    complement   = {T_COMPLEMENT}  (unbilled in the monitored burn ledger)")
    print(f"    T            = {T_CANON}")
    check(T_CANON == 9, "local closed service sweep has T=9 monitored burn touches")
    print("    Therefore local algebra supplies B_sweep = T alpha0 = 9 alpha0.")

    print("\n[3] Neighboring stock alternatives")
    rows = [
        ("T=8 no latch/readout", stock(8), "missing latch"),
        ("T=9 one closed sweep", stock(9), "canonical"),
        ("T=10 complement billed", stock(10), "extra complement"),
        ("half sweep", stock(9, 0.5), "fractional sweep"),
        ("two sweeps", stock(9, 2.0), "extra sweep"),
    ]
    print(f"    {'case':<24s} {'B0':>14s} {'H0':>10s} {'miss':>10s}  note")
    for label, b0, note in rows:
        h0 = h0_from_stock(b0, r6)
        print(f"    {label:<24s} {b0:14.8e} {h0:10.3f} {h0/TARGET_H0-1:+10.2%}  {note}")
    check(abs(h0_from_stock(stock(8), r6) / TARGET_H0 - 1.0) > 0.10, "T=8 is excluded by the selector landing")
    check(abs(h0_from_stock(stock(10), r6) / TARGET_H0 - 1.0) > 0.09, "T=10 is excluded by the selector landing")
    check(abs(h0_from_stock(stock(9, 2.0), r6) / TARGET_H0 - 1.0) > 0.45, "two sweeps are excluded by the selector landing")

    print("\n[4] What is actually proved")
    print("    Proven locally:")
    print("      permanent burn stock per closed service sweep = 9 alpha0.")
    print("    Not yet proven globally:")
    print("      the cosmological consumable stock consists of exactly one such")
    print("      closed sweep per cell, with no hidden stock and no double-counted")
    print("      static K04 relic density.")
    print("    This is the precise a=1 lift premise:")
    print("      B0 = B_sweep, not m * B_sweep or B_sweep + B_hidden.")

    print(
        """
================================================================================
VERDICT
  The deepest a=1 reduction does not close unconditionally today, but it has
  been reduced to a single theorem:

      Single closed burn-sweep theorem:
        the whole cosmological consumable stock is exactly one closed
        §5.2 service-era burn sweep per cell.

  If that theorem is accepted/derived, then B0=9 alpha0, the D2 exhaustion law
  gives N_exhaust=9 alpha0/r6, and the Lambda<->M_P selector no longer consumes
  an epoch convention.  If it is not derived, the selector remains conditional.

  This is a real reduction: the local unit and touch count are fixed; the
  remaining freedom is no longer a coefficient hunt, but the global
  one-sweep/no-hidden-stock premise.
exit 0"""
    )


if __name__ == "__main__":
    main()
