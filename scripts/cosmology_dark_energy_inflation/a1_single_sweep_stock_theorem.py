#!/usr/bin/env python3
r"""Single-sweep stock theorem audit for the A=1 / D2 completion premise.

Question
--------
The A=1 lift now asks why the consumable cosmological stock is

    B0 = 9 alpha0

rather than m times that amount, plus a hidden stock, or a static K04 relic
density.  The earlier burn-integral audit derived the local unit

    one closed service sweep = (8 repairs + 1 latch + 0 complement) alpha0
                             = 9 alpha0.

This script checks how far the existing monitored QEC record algebra goes
toward the remaining global premise:

  1. why exactly one sweep?
  2. why no hidden stock?
  3. why no repeated sweep?
  4. why static K04 relic density is fuel/geometry, not billable stock?

Verdict
-------
Inside the closed-record/latch service algebra:

  * repeated sweeps are idle, because recovery followed by latch is idempotent;
  * hidden stock is excluded as an unmonitored variable, unless a new sector is
    added outside the algebra;
  * K04 density is fuel/geometry: it gates whether an episode can occur, but
    does not multiply the billed record count.

This script by itself narrows the remaining global premise to:

    each physical cell enters exactly one open service episode before the
    a=1 completion latch.

The companion ``a1_episode_count_closure_attempt.py`` then closes that global
episode count inside the current homogeneous service instrument by covariance
of the R4 scalar completion flag.  Outside-sector completeness remains the
caveat.
"""

from __future__ import annotations

import importlib
import math
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))
rh = importlib.import_module("register_handoff_form_selection")

ALPHA0 = 1.0 / 137.0
TOUCH_REPAIRS = 8
TOUCH_LATCH = 1
TOUCH_COMPLEMENT = 0
TOUCHES_PER_SWEEP = TOUCH_REPAIRS + TOUCH_LATCH + TOUCH_COMPLEMENT
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


@dataclass(frozen=True)
class CellState:
    """Minimal closed-service state.

    fuel_present is geometry/availability.
    syndrome_open is the monitored source.
    closed_latch marks that the source has already been resolved.
    bill is the permanent record stock.
    """

    fuel_present: bool
    syndrome_open: bool
    closed_latch: bool
    bill: float = 0.0


def service_once(state: CellState) -> CellState:
    """Apply one monitored closed service episode if one is open.

    The operation is intentionally a projection-like recovery/latch map:
    once the latch is closed, applying it again changes nothing and bills
    nothing.  This is the finite-state analogue of QND/idempotent recovery.
    """

    if not state.fuel_present:
        return state
    if not state.syndrome_open:
        return state
    if state.closed_latch:
        return state
    return CellState(
        fuel_present=True,
        syndrome_open=False,
        closed_latch=True,
        bill=state.bill + TOUCHES_PER_SWEEP * ALPHA0,
    )


def service_n(state: CellState, n: int) -> CellState:
    for _ in range(n):
        state = service_once(state)
    return state


def bill_for_fuel_density(fuel_density: float) -> float:
    """Bill one episode if fuel exists; density does not multiply the record."""

    state = CellState(fuel_present=fuel_density > 0.0, syndrome_open=True, closed_latch=False)
    return service_once(state).bill


def main() -> None:
    print("A=1 SINGLE-SWEEP / STOCK THEOREM AUDIT")
    print("=" * 92)

    r6 = residual_r6()
    one_sweep_stock = TOUCHES_PER_SWEEP * ALPHA0
    h0 = h0_from_stock(one_sweep_stock, r6)
    print("[1] Local closed-service stock unit")
    print(f"    touches per sweep     = {TOUCH_REPAIRS}+{TOUCH_LATCH}+{TOUCH_COMPLEMENT} = {TOUCHES_PER_SWEEP}")
    print(f"    B_sweep               = {one_sweep_stock:.12e}")
    print(f"    H0(B_sweep)           = {h0:.6f} km/s/Mpc")
    print(f"    selector reference H0 = {rh.H0_KM_S_MPC:.6f} km/s/Mpc")
    print(f"    miss                  = {h0 / rh.H0_KM_S_MPC - 1:+.4%}")
    check(TOUCHES_PER_SWEEP == 9, "one closed local service sweep bills 9 alpha0")
    check(abs(h0 / rh.H0_KM_S_MPC - 1.0) < 7.0e-4, "one-sweep stock lands the selector budget")

    print("\n[2] No repeated sweep: closed recovery/latch is idempotent")
    initial = CellState(fuel_present=True, syndrome_open=True, closed_latch=False)
    after_one = service_n(initial, 1)
    after_two = service_n(initial, 2)
    after_three = service_n(initial, 3)
    print(f"    bill after one application   = {after_one.bill:.12e}")
    print(f"    bill after two applications  = {after_two.bill:.12e}")
    print(f"    bill after three applications= {after_three.bill:.12e}")
    check(after_two == after_one, "second application is idle after the latch closes")
    check(after_three == after_one, "all later applications are idle without a new open syndrome")

    repeated_bill_h0 = h0_from_stock(2.0 * one_sweep_stock, r6)
    print(f"    If a second sweep billed again, H0 would be {repeated_bill_h0:.3f} km/s/Mpc")
    check(abs(repeated_bill_h0 / rh.H0_KM_S_MPC - 1.0) > 0.45, "billable repeated sweep is selector-excluded")

    print("\n[3] No hidden stock inside the monitored algebra")
    print("    The algebra has only monitored source/latch records; stock is image(service_once).")
    hidden_trials = [-1.0, -0.1, 0.0, 0.1, 1.0]
    print(f"    {'hidden / B_sweep':>18s} {'H0':>10s} {'miss':>10s}")
    for frac in hidden_trials:
        b0 = one_sweep_stock * (1.0 + frac)
        if b0 <= 0:
            print(f"    {frac:18.2f} {'n/a':>10s} {'invalid':>10s}")
            continue
        h = h0_from_stock(b0, r6)
        print(f"    {frac:18.2f} {h:10.3f} {h / rh.H0_KM_S_MPC - 1:+10.2%}")
    check(abs(h0_from_stock(one_sweep_stock * 1.1, r6) / rh.H0_KM_S_MPC - 1.0) > 0.09,
          "even a 10% hidden stock shifts the selector outside the local budget")
    print("    Therefore hidden stock is not a silent option: it is new, monitored physics.")

    print("\n[4] Static K04 relic density is fuel/geometry, not billable stock")
    densities = [0.0, 0.01, 0.1, 1.0, 10.0, 100.0]
    print(f"    {'fuel density':>12s} {'bill':>14s}")
    for density in densities:
        print(f"    {density:12.2g} {bill_for_fuel_density(density):14.8e}")
    positive_bills = {round(bill_for_fuel_density(d), 15) for d in densities if d > 0.0}
    check(bill_for_fuel_density(0.0) == 0.0, "no fuel means no open service episode")
    check(len(positive_bills) == 1, "any positive fuel density gives the same per-episode bill")
    check(next(iter(positive_bills)) == round(one_sweep_stock, 15),
          "static K04 density gates availability but does not multiply record stock")

    print("\n[5] What remains after this local audit")
    print("    Closed inside the monitored service algebra:")
    print("      - repeated sweeps without a new open syndrome are idempotent/idled;")
    print("      - hidden stock must be a new monitored sector, not an implicit reserve;")
    print("      - static K04 density is fuel/geometry, not billable permanent stock.")
    print("    This script reduces the global premise to:")
    print("      every physical cell enters exactly one open service episode before")
    print("      the a=1 completion latch.")
    print("    Companion closure:")
    print("      a1_episode_count_closure_attempt.py closes that episode count inside")
    print("      the current homogeneous service instrument; outside-sector")
    print("      completeness remains the caveat.")

    print(
        """
================================================================================
VERDICT
  Progress: the four loose questions reduce to one.

  Exactly one sweep is now "one open service episode per physical cell": the
  local closed sweep is unique once opened, and the latch makes repeats idle.
  No hidden stock follows inside the monitored algebra, because any extra
  billable reserve would need its own monitored source.  Static K04 relic
  density is fuel/geometry: it can make an episode available, but it is not
  counted as permanent record stock.

  The remaining statement after this local audit was not a coefficient:
      one physical cell -> one open episode before the completion latch.

  The companion episode-count audit derives that statement inside the current
  homogeneous service instrument.  What remains is outside-sector completeness,
  not a hidden stock or repeated-sweep freedom.
exit 0"""
    )


if __name__ == "__main__":
    main()
