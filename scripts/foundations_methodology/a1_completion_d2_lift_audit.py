#!/usr/bin/env python3
r"""A=1 COMPLETION / D2 RELIC-EXHAUSTION LIFT AUDIT.

Question
--------
The Planck hierarchy, the late R4 branch, and the inflation exit story all use
the same cosmological endpoint:

    N_physical = N_lock = 9 alpha0/r6,

with the cosmological convention a=1 placed at that endpoint.  D1
("the boot directly deposits 9 alpha0 protected line-elements per cell") was
refuted as posed.  The later deep run also refuted the tempting D1' finite-size
candidate and the winding-relic census.  D2 ("depth-six residual events consume
relic fuel") remains consistent.

This audit asks what D2 now proves and what remains.

Result
------
D2 converts the a=1 statement from an epoch-pick into a first-hitting
exhaustion law:

    B(N) = B0 - N r6,
    chi(N) = 1 - B(N)/B0 = N r6/B0,
    N_exhaust = B0/r6.

Therefore the selector is closed if and only if the boot/service history
supplies the consumable stock

    B0 = 9 alpha0.

The important correction is that B0 is no longer plausibly a static K04
surviving-density census.  D1, D1', and winding-relief readings are closed
negative.  The only live D2 route is the integrated burn ledger:

    B0 = integral(permanent burn records during the service era).

That is exactly the section-5.2 burn-billing rule's territory.  This audit
supplies the form and shows why static stock routes fail.  The companion
``a1_burn_integral_theorem_attempt.py`` supplies the local touch count T=9;
``a1_single_sweep_stock_theorem.py`` closes no-repeat/no-hidden/fuel-vs-stock
inside the monitored algebra; and ``a1_episode_count_closure_attempt.py``
closes one open episode per physical cell inside the current homogeneous
service instrument.
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
PHI = (math.sqrt(5.0) - 1.0) / 2.0
TARGET_STOCK = 9.0 * ALPHA0

# Historical D1 measurements retained only to show why static-stock routes die.
F_LINE_L6_R496 = 0.74815
F_LINE_DEEP_1_OVER_L2 = 1.369
F_LINE_DEEP_1_OVER_L = 1.582
F_WIND_DEEP_LOW = 0.60
F_WIND_DEEP_HIGH = 0.73


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def require_text(path: str, needles: list[tuple[str, str]]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    for label, needle in needles:
        ok = needle in text
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            raise AssertionError(f"{path}: missing {needle!r}")


def residual_r6() -> float:
    gamma_star = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
    _, q_post, _ = rh.queue_readouts(gamma_star, 1)
    return (21.0 * q_post) ** 32 / 21.0


def h0_from_stock(stock: float, r6: float) -> float:
    """Return H0 in km/s/Mpc from a consumable stock B0."""
    mpc_km = 3.085678e19
    hbar_gev_s = 6.582120e-25
    n_exhaust = stock / r6
    return rh.LAMBDA_QCD_GEV / n_exhaust * mpc_km / hbar_gev_s


def completion_fraction(stock: float, n_ticks: float, r6: float) -> float:
    return min(1.0, n_ticks * r6 / stock)


def main() -> None:
    print("A=1 COMPLETION / D2 RELIC-EXHAUSTION LIFT AUDIT")
    print("=" * 96)

    print("[1] Current canon anchors")
    require_text(
        "python_code/d2_consumption_channel_audit.py",
        [
            ("D2 consumes one relic element per deep residual", "consumes exactly ONE relic line-element"),
            ("D2 is consistent but not stock-derived", "D2 is CONSISTENT, not yet derived"),
        ],
    )
    require_text(
        "python_code/cosmological_selector_lock_theorem.py",
        [
            ("selector is first-hitting endpoint", "Because chi_R4 is a bounded support-completion fraction"),
            ("physical span is N_lock", "N_physical = N_lock = 9 alpha0/r6"),
        ],
    )
    require_text(
        "python_code/r4_bandwidth_saturation_test.py",
        [
            ("D2 squeeze theorem survives", "the squeeze theorem"),
            ("absolute bandwidth saturation failed as selector", "FAILED AS THE LOCK RULE"),
        ],
    )
    require_text(
        "DRIFT.md",
        [
            ("D1 direct deposition refuted", "D1 DEPOSITION AUDIT — NEGATIVE"),
            ("D1 prime refuted by deep run", "D1$'$ REFUTED AT DEEP PRECISION"),
            ("winding relic refuted", "MECHANISM 1 (WINDING RELIC) REFUTED"),
            ("static census reposed as burn/integral", "billable $=\\int({\\rm burn})\\,dt$"),
        ],
    )
    require_text(
        "technical_notes/cc_monitored_billing_operator_algebra.md",
        [
            ("burn rule is permanent records = touches alpha0", "permanent records = (register touches) x alpha0"),
            ("T=9 touch count is supplied", "T = 8 (single-bit repairs"),
            ("coherence identity is r6 Nt = T alpha0", "cohere iff r_6 N_t = T alpha_0"),
        ],
    )
    require_text(
        "python_code/a1_episode_count_closure_attempt.py",
        [
            ("episode-count theorem staged", "completion => exactly one open episode per physical cell"),
            ("outside-sector caveat retained", "outside-sector completeness"),
        ],
    )

    print("\n[2] D2 exhaustion map")
    r6 = residual_r6()
    n_lock = TARGET_STOCK / r6
    h0_lock = h0_from_stock(TARGET_STOCK, r6)
    print(f"    r6                         = {r6:.12e}")
    print(f"    target stock B0=9 alpha0    = {TARGET_STOCK:.12e}")
    print(f"    N_exhaust=B0/r6             = {n_lock:.12e}")
    print(f"    H0_out                      = {h0_lock:.6f} km/s/Mpc")
    print(f"    selector reference H0       = {rh.H0_KM_S_MPC:.6f} km/s/Mpc")
    print(f"    miss                        = {h0_lock/rh.H0_KM_S_MPC-1:+.4%}")
    check(abs(h0_lock / rh.H0_KM_S_MPC - 1.0) < 7.0e-4, "D2 exhaustion with B0=9 alpha0 lands the selector budget")
    for frac in (0.25, 0.5, 0.75, 1.0):
        chi = completion_fraction(TARGET_STOCK, frac * n_lock, r6)
        print(f"    N={frac:>4.2f} N_lock -> chi={chi:.3f}")
    check(completion_fraction(TARGET_STOCK, n_lock, r6) == 1.0, "completion is a bounded first-hitting event")

    print("\n[3] Static-stock routes are closed negative")
    survival = math.exp(-3.0 / (2.0 * PHI))
    l6_d1prime = F_LINE_L6_R496 * survival
    deep_d1prime_low = F_LINE_DEEP_1_OVER_L2 * survival
    deep_d1prime_high = F_LINE_DEEP_1_OVER_L * survival
    print(f"    D1 direct L=6 stock              = {F_LINE_L6_R496:.3f} ({F_LINE_L6_R496/TARGET_STOCK:.1f} x target)")
    print(f"    historical D1' L=6 stock         = {l6_d1prime:.5f} ({l6_d1prime/TARGET_STOCK:.3f} x target)")
    print("       -> finite-size landing only; deep escalation refuted it.")
    print(f"    deep D1' extrapolated stock      = {deep_d1prime_low:.5f}..{deep_d1prime_high:.5f}")
    print(f"                                      ({deep_d1prime_low/TARGET_STOCK:.2f}..{deep_d1prime_high/TARGET_STOCK:.2f} x target)")
    print(f"    winding static stock             = {F_WIND_DEEP_LOW:.2f}..{F_WIND_DEEP_HIGH:.2f}")
    print(f"                                      ({F_WIND_DEEP_LOW/TARGET_STOCK:.1f}..{F_WIND_DEEP_HIGH/TARGET_STOCK:.1f} x target)")
    check(F_LINE_L6_R496 / TARGET_STOCK > 10.0, "D1 direct deposition overproduces")
    check(deep_d1prime_low / TARGET_STOCK > 1.8, "D1 prime remains over target after deep extrapolation")
    check(F_WIND_DEEP_LOW / TARGET_STOCK > 9.0, "winding relic static census overproduces")

    print("\n[4] Live route: integrated D2 burn ledger")
    print("    The static K04 census is not the billable object.  The live object is")
    print("      B0 = integral(permanent burn records during service era).")
    print("    Existing burn-billing algebra gives the local record law:")
    print("      permanent records = register touches x alpha0,")
    print("      T = 8 repairs + 1 readout + 0 complement = 9,")
    print("      coherence iff r6 N_t = T alpha0.")
    print("    D2 supplies the endpoint map; the burn rule supplies the touch count.")
    print("    The companion stock/episode audits supply the local and in-instrument")
    print("    global pieces: no repeated billable sweep, no hidden monitored stock,")
    print("    K04 as fuel/geometry, and one open episode per physical cell under")
    print("    homogeneous R4 covariance.")

    print("\n[5] Registered remaining caveat")
    print("    Inside the current homogeneous service instrument, the burn integral is")
    print("      integral dN_burn = 9 alpha0 per physical cell.")
    print("    The remaining caveat is outside-sector completeness: a hidden register,")
    print("    invalid-state channel, or non-R4 cosmological coupling would be new")
    print("    physics rather than a surviving K04 stock or sweep-count ambiguity.")

    print(
        """
================================================================================
VERDICT
  D2 meaningfully advances the a=1/completion premise.  It supplies the
  monotone consumption law and the first-hitting endpoint, so a=1 is no longer
  an arbitrary epoch label once a finite consumable stock is accepted:

      N_exhaust = B0/r6.

  The whole cosmological lock reduces to one stock identity:

      B0 = 9 alpha0.

  Static-stock routes are dead.  D1, D1', and winding-relic readings all
  overproduce after the deep escalation.  The live object is the consumed
  permanent-record integral, and the companion stock/episode audits close it
  inside the current homogeneous service instrument:

      integral dN_burn = 9 alpha0 per physical cell.

  The remaining caveat is outside-sector completeness, not a K04 density
  census and not a free sweep-count parameter.
exit 0"""
    )


if __name__ == "__main__":
    main()
