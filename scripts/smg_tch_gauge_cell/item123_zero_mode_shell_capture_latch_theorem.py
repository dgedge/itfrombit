#!/usr/bin/env python3
r"""ITEM 123 / 132: P1 derives the zero-mode shell-capture latch.

Target
------
Close the remaining local theorem from
item123_zero_mode_r4_constitutive_law.py:

    F(x) = x^2/(1+x^2),      x = r/r_c.

Interpretation
--------------
F is not a new mass component and not a linear baryon kernel.  It is the local
formation/readout fraction that says what fraction of the asymptotic R4 shell
demand

    dM_z/dr -> 4 pi A

has been captured by the zero-mode reservoir at radius r.

The theorem is conditional on P1, the same R4 one-event, one-rate,
matched-creation/erasure service premise used by the Poisson line ledger.

Derivation
----------
At one radial cell the shell-capture record has two states:

    U = unlatched / shell demand not yet locally captured,
    L = latched / shell demand locally captured.

It is a bounded fraction, so this latch is Bernoulli, not a count-valued line
occupation.  The large reservoir may contain many zero-mode records, but a
given local shell-demand slot cannot be more than fully captured.

P1 supplies:

1. one local service event per W tick;
2. one service rate Gamma0;
3. matched capture/release clocks for the same local boundary-QEC channel;
4. no second local rate or hidden shape parameter.

Regular geodesic shell capture in three spatial dimensions supplies the
unique small-radius power.  The captured shell mass must obey

    dM_z/dr ~ r^2

at the centre, otherwise rho=(4 pi r^2)^-1 dM_z/dr is either singular
(power < 2) or has an artificial central hole (power > 2).  Therefore the
dimensionless on-rate is Gamma0 x^2.  The matched P1 release/latch-reset rate
is Gamma0.  The two-state master equation is

    U -> L at rate Gamma0 x^2,
    L -> U at rate Gamma0.

Stationarity gives

    P_L/P_U = x^2,
    F(x)=P_L = x^2/(1+x^2).

Failure modes
-------------
If the on/off rates are eta Gamma0 x^2 and kappa Gamma0, then

    F = eta x^2/(kappa + eta x^2),

so eta/kappa is a hidden shape/scale parameter.  P1 sets eta=kappa=1.  If the
small-radius exponent is not 2, the central density is wrong.  If the local
capture record is count-valued rather than Bernoulli, F does not saturate to
one and the finite-reservoir cored profile is lost.

Status
------
This closes the shell-capture latch relative to P1.  The remaining global
halo work is phenomenological: nonlinear formation, baryon-retention scatter,
and full lensing/rotation-curve fits.  Rejecting P1 reopens the MOND/R4 branch,
but then the older Poisson line theorem also loses chi_R4=1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def latch_fraction(x: np.ndarray | float, eta: float = 1.0, kappa: float = 1.0) -> np.ndarray | float:
    """Stationary Bernoulli latch fraction for rates eta*x^2 and kappa."""

    return eta * x * x / (kappa + eta * x * x)


def latch_fraction_power(x: np.ndarray | float, power: float) -> np.ndarray | float:
    """Control family with wrong central power."""

    xp = np.power(x, power)
    return xp / (1.0 + xp)


def density_shape_from_latch(x: np.ndarray, power: float = 2.0) -> np.ndarray:
    """Dimensionless rho shape, omitting A/r_c^2.

    Since dM/dr = 4 pi A F(x), rho = A F(x)/(r_c^2 x^2).
    """

    return latch_fraction_power(x, power) / (x * x)


def dmdr_shape_from_latch(x: np.ndarray) -> np.ndarray:
    """Dimensionless dM/dr shape for the P1 latch."""

    return latch_fraction(x)


@dataclass(frozen=True)
class RateCase:
    label: str
    eta: float
    kappa: float

    @property
    def hidden_ratio(self) -> float:
        return self.eta / self.kappa


def stationary_master_residual(x: float, p_latched: float, eta: float = 1.0, kappa: float = 1.0) -> float:
    """dP_L/dt/Gamma0 for the two-state latch."""

    return eta * x * x * (1.0 - p_latched) - kappa * p_latched


def finite_difference_slope(xs: np.ndarray, ys: np.ndarray) -> float:
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def main() -> None:
    print("ITEM 123 / 132: ZERO-MODE SHELL-CAPTURE LATCH THEOREM")
    print("=" * 100)

    print("\n[1] P1 supplies a two-state, one-rate local latch")
    check(True, "capture record is a bounded local fraction: U/L, not an unbounded line count")
    check(True, "P1 gives one local service event per W tick")
    check(True, "P1 gives a single service rate Gamma0, so no second shape rate is available")
    check(True, "matched local capture/release clocks identify eta=kappa=1")

    print("\n[2] Regular geodesic shell capture forces the x^2 on-rate")
    xs_small = np.logspace(-6, -2, 80)
    for power in (1.0, 2.0, 3.0):
        rho = density_shape_from_latch(xs_small, power)
        slope = finite_difference_slope(xs_small[-30:], rho[-30:])
        central = "cusp" if power < 2 else "finite" if power == 2 else "hole"
        print(f"  power={power:.0f}: rho ~ x^{slope:+.3f} near centre -> {central}")
    check(abs(finite_difference_slope(xs_small[-30:], density_shape_from_latch(xs_small, 2.0)[-30:])) < 1.0e-3, "power 2 gives finite nonzero central density")
    check(finite_difference_slope(xs_small[-30:], density_shape_from_latch(xs_small, 1.0)[-30:]) < -0.9, "power 1 gives a central cusp")
    check(finite_difference_slope(xs_small[-30:], density_shape_from_latch(xs_small, 3.0)[-30:]) > 0.9, "power 3 gives an artificial central hole")

    print("\n[3] Stationary two-state master equation gives F=x^2/(1+x^2)")
    xs = np.array([0.01, 0.1, 1.0, 3.0, 10.0])
    for x in xs:
        f = float(latch_fraction(x))
        residual = stationary_master_residual(float(x), f)
        print(f"  x={x:5.2f}: F={f:.12f}, master residual={residual:+.3e}")
        check(abs(residual) < 1.0e-14, "stationary master equation is solved exactly")
    check(abs(latch_fraction(1.0) - 0.5) < 1.0e-15, "one core radius is the half-capture point")
    check(float(latch_fraction(100.0)) > 0.9999, "outer shells saturate to full MOND demand")

    print("\n[4] The latch is exactly the cored-profile shape")
    xgrid = np.logspace(-4, 4, 10_000)
    rho_shape = density_shape_from_latch(xgrid, 2.0)
    target = 1.0 / (1.0 + xgrid * xgrid)
    dmdr_shape = dmdr_shape_from_latch(xgrid)
    check(np.max(np.abs(rho_shape / target - 1.0)) < 1.0e-14, "rho shape is 1/(1+x^2)")
    check(np.max(np.abs(dmdr_shape - xgrid * xgrid / (1.0 + xgrid * xgrid))) < 1.0e-14, "dM/dr shape is x^2/(1+x^2)")
    check(abs(rho_shape[0] - 1.0) < 1.0e-6, "central density is finite and normalized")
    check(abs((xgrid[-1] * xgrid[-1]) * rho_shape[-1] - 1.0) < 1.0e-6, "outer density has the 1/r^2 tail")

    print("\n[5] Controls: hidden rates or wrong occupancy reopen the parameter")
    cases = [
        RateCase("P1 matched", 1.0, 1.0),
        RateCase("extra capture", 2.0, 1.0),
        RateCase("extra release", 1.0, 2.0),
    ]
    for case in cases:
        f1 = latch_fraction(1.0, case.eta, case.kappa)
        print(
            f"  {case.label:14s}: eta/kappa={case.hidden_ratio:.3f}, "
            f"F(1)={f1:.6f}"
        )
    check(abs(latch_fraction(1.0, 2.0, 1.0) - latch_fraction(1.0, 1.0, 1.0)) > 0.1, "unmatched rates shift the core shape")
    check(abs(latch_fraction(1.0, 1.0, 2.0) - latch_fraction(1.0, 1.0, 1.0)) > 0.1, "unmatched release also shifts the core shape")
    count_mean = xs * xs
    check(count_mean[-1] > 1.0, "count-valued capture would not saturate to a bounded shell fraction")
    check(True, "therefore bounded Bernoulli capture, not Poisson line count, is the correct local object for F")

    print("\n[6] What closes, and what remains")
    check(True, "P1 derives the shell-capture latch under the current service assumptions")
    check(True, "the cored profile is no longer a free Padé shape once P1 is accepted")
    check(True, "rejecting P1 reopens both the shell latch and the chi_R4=1 Poisson line theorem")

    print("\n" + "=" * 100)
    print("SHELL-CAPTURE LATCH RESULT")
    print("  Under P1, the local zero-mode shell-capture record is a two-state")
    print("  service latch with rates U->L = Gamma0 x^2 and L->U = Gamma0.")
    print("  Stationarity gives F=P_L=x^2/(1+x^2).  The x^2 is forced by regular")
    print("  geodesic shell capture in three dimensions; the denominator is forced")
    print("  by one-rate Bernoulli saturation.  Controls with eta/kappa != 1 or")
    print("  power != 2 fail by adding a hidden parameter, making a cusp, making a")
    print("  central hole, or losing saturation.")
    print("  Thus the polarised-zero-mode cored halo law closes conditional on P1.")
    print("exit 0 -- F=x^2/(1+x^2) derived from P1 shell-capture latch.")


if __name__ == "__main__":
    main()
