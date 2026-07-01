#!/usr/bin/env python3
r"""Dressed-alpha positive-closure gate.

Question
--------
Can the remaining inverse-alpha residual

    137 -> 137.035999...

be closed positively from the current substrate machinery?

This script tests the only remaining non-fit closure classes after the
Ward/Kubo and charge-weighted no-go audits:

  1. The physical EM observable is the Ward/Kubo Peierls-current Hessian.
     Can a finite multiplicative normalisation of that Hessian hit the
     measured residual, with the multiplier selected by existing canon?

  2. Failing that, can an additive local F^2 counterterm be pinned by an
     existing service/geometry number rather than by the observed alpha?

  3. Are the attractive numerical factors sparse enough to carry evidence,
     or are they drawn from a dense small-rational search space?

Result
------
The gate does not close.  The required multiplicative Ward/Kubo factor is

    Z_fin = 2.537609,     q_eff = sqrt(Z_fin) = 1.592988,

which is not the unit Peierls charge and is not selected by the present
service algebra.  The additive alternative is a finite counterterm equal to
0.605928 of the observed residual, again unselected.  Small rationals near
both numbers are dense enough that promoting one would repeat the N1=31
failure mode.

So the positive derivation still needs genuinely new physics: either an EM
sector-billing theorem or a finite normal-ordering rule that fixes the
Ward/Kubo term before comparison with alpha.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import dressed_alpha_bridge_web_open_system as bw
import dressed_alpha_ward_kubo_peierls_observable_audit as wk


ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
DELTA_TARGET = ALPHA_PHYS_INV - ALPHA0_INV


@dataclass(frozen=True)
class Hit:
    err: float
    name: str
    value: float


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def rational_hits(x: float, max_den: int, max_num: int, tol: float) -> list[Hit]:
    hits: list[Hit] = []
    seen: set[tuple[int, int]] = set()
    for q in range(1, max_den + 1):
        for p in range(1, max_num + 1):
            g = math.gcd(p, q)
            pr, qr = p // g, q // g
            if (pr, qr) in seen:
                continue
            seen.add((pr, qr))
            value = pr / qr
            err = abs(value / x - 1.0)
            if err <= tol:
                hits.append(Hit(err, f"{pr}/{qr}", value))
    return sorted(hits, key=lambda h: h.err)


def named_factor_hits(x: float) -> list[Hit]:
    candidates = {
        "sqrt(2*pi)": math.sqrt(2.0 * math.pi),
        "8/pi": 8.0 / math.pi,
        "pi^2/4": math.pi * math.pi / 4.0,
        "5/2": 2.5,
        "81/32": 81.0 / 32.0,
        "137/54": 137.0 / 54.0,
        "33/13": 33.0 / 13.0,
        "sqrt(2*pi)*137/136": math.sqrt(2.0 * math.pi) * 137.0 / 136.0,
    }
    return sorted(
        (Hit(abs(v / x - 1.0), name, v) for name, v in candidates.items()),
        key=lambda h: h.err,
    )


def main() -> None:
    print("DRESSED-ALPHA POSITIVE-CLOSURE GATE")
    print("=" * 96)

    required_kernel = DELTA_TARGET * 2.0 * math.pi / ALPHA0
    ward = wk.ward_kubo_dispersion()
    ward_ratio = float(ward["ratio_web"])
    ward_delta = ward_ratio * DELTA_TARGET
    missing_delta = DELTA_TARGET - ward_delta
    z_fin = 1.0 / ward_ratio
    q_eff = math.sqrt(z_fin)
    counterterm_fraction = missing_delta / DELTA_TARGET

    print("[1] Target residual and physical Ward/Kubo slot")
    print(f"  alpha0^-1                         = {ALPHA0_INV:.0f}")
    print(f"  alpha_phys^-1                     = {ALPHA_PHYS_INV:.9f}")
    print(f"  target residual delta(alpha^-1)   = {DELTA_TARGET:.9f}")
    print(f"  one-loop-style required kernel    = {required_kernel:.6f}")
    print(f"  Ward/Kubo web-dressed ratio       = {ward_ratio:.6f}")
    print(f"  Ward/Kubo residual contribution   = {ward_delta:.9f}")
    print(f"  missing residual after Ward slot  = {missing_delta:.9f}")
    check(30.8 < required_kernel < 31.1, "observed residual is the historical N1~31 kernel")
    check(0.35 < ward_ratio < 0.45, "physical Ward/Kubo slot undershoots the observed residual")

    print("\n[2] Multiplicative finite-normalisation gate")
    print(f"  Z_fin needed                      = {z_fin:.6f}")
    print(f"  equivalent Peierls charge q_eff   = {q_eff:.6f}")
    print("  Existing Peierls current fixes the unit-charge derivative.  Therefore")
    print("  q_eff != 1 is not a derived renormalisation unless a new sector rule")
    print("  selects it independently of alpha.")
    check(abs(q_eff - 1.0) > 0.5, "required q_eff is a new charge normalisation, not unit Peierls charge")

    print("\n  Named factors near Z_fin:")
    for hit in named_factor_hits(z_fin):
        print(f"    {hit.name:<22} value={hit.value:.9f}  rel.err={hit.err:.6e}")

    rational = rational_hits(z_fin, max_den=64, max_num=180, tol=5.0e-3)
    print(f"\n  Small-rational hits within 0.5% of Z_fin (denominator <= 64): {len(rational)}")
    for hit in rational[:10]:
        print(f"    {hit.name:<7} value={hit.value:.9f}  rel.err={hit.err:.6e}")
    check(len(rational) > 20, "nearby small rationals are dense; numerical closeness is not evidence")

    print("\n[3] Additive F^2 counterterm gate")
    print(f"  counterterm needed                = {missing_delta:.9f}")
    print(f"  counterterm / observed residual   = {counterterm_fraction:.6f}")
    add_hits = rational_hits(counterterm_fraction, max_den=64, max_num=80, tol=5.0e-3)
    print(f"  Small-rational hits within 0.5% of counterterm fraction: {len(add_hits)}")
    for hit in add_hits[:10]:
        print(f"    {hit.name:<7} value={hit.value:.9f}  rel.err={hit.err:.6e}")
    print("  Gauge invariance permits a local F^2 finite term, but the current canon")
    print("  has no rule choosing this coefficient.  Selecting 20/33, 17/28, etc.")
    print("  would be a fit unless derived from a new normal-ordering theorem.")
    check(0.55 < counterterm_fraction < 0.65, "additive rescue is a large finite term")
    check(len(add_hits) > 5, "nearby additive rational hits are also dense")

    print("\n[4] Closure decision")
    print(
        """
  The positive derivation does not close under the current axioms.

  What is closed:
    * alpha0^-1 = 137 as the monitored record-pair service weight.
    * the old N1=31 route as a structural no-go for physical dressing.
    * the Ward/Kubo observable identity: physical EM dressing bills the
      Peierls current-current / photon-self-energy slot.

  What remains unclosed:
    * a theorem selecting Z_fin = 2.537609 in that Ward/Kubo slot; or
    * a theorem selecting an additive finite F^2 counterterm equal to
      0.605928 of the observed residual; or
    * a sector-billing theorem proving that QED alpha measures the monitored
      service occupation instead of the Ward/Kubo self-energy.

  Since none of those selectors is present, the measured residual remains a
  bounded electromagnetic residual.  The next real route is not another number
  search; it is a new EM sector-billing or normal-ordering theorem.
exit 0 -- positive closure rejected under current canon.
"""
    )


if __name__ == "__main__":
    main()
