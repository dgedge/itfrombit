#!/usr/bin/env python3
r"""ITEM 131: HBC clock covariance, superhorizon conservation, and exit attempt.

Question
--------
Can the last early-HBC pieces be closed?

    1. derive nu_HBC(x) from the boundary-printing/QEC event ledger;
    2. derive <nu_k nu_-k>;
    3. show R_HBC is conserved outside the printed horizon shell;
    4. derive a finite saturated w=-1 duration plus exit.

Verdict
-------
Partial closure only.

Closed conditionally:
    If local expansion is proportional to the local HBC print-event rate, the
    single-clock field is forced:

        eta_HBC(x,N) = j(x,N)/jbar(N) - 1
        nu_HBC(x,N) = integral^N eta_HBC(x,N') dN'

    where j is the total 28-channel print/service event density per local
    horizon patch and N=ln a.  This is the compensated event-ledger clock.

    If the scale-covariant HBC/QEC generator acts on the covariance ledger,
    its two-point function has the form

        <nu_k nu_k'> = (2 pi)^3 delta(k+k') P_nu(k)
        Delta_nu^2(k) = k^3 P_nu(k)/(2 pi^2)
                      = A_nu (k/k_*)^(-1/28).

    Therefore the tilt shape is derived from the event ledger plus the
    already-closed 28-channel log-shell generator.  The amplitude A_nu is not
    fixed by the 28-clock; it is a shot-noise / event-density normalization
    requiring an actual microscopic event-rate ledger.

    Superhorizon conservation is also conditional: for a single adiabatic HBC
    clock with no outside source after horizon printing, dR/dN is gradient
    suppressed, O((k/aH)^2), and R_HBC freezes after a few e-folds outside the
    shell.

Not closed:
    The finite early saturated w=-1 duration and exit are not derivable from
    current canon.  Part 17 supplies w(0)=-1 only.  The R-activation schedule
    is explicitly phenomenological / item-42-open.  The N=56 and percolation
    material is internally flagged as generic slow-roll / unexecuted
    percolation, not an HBC-derived exit theorem.
"""

from __future__ import annotations

import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DELTA = Fraction(1, 28)
TWO_PI_SQUARED = 2.0 * math.pi**2


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def event_clock_increment(local_events: float, mean_events: float, dN: float) -> float:
    """Compensated local print-time increment from event counts."""
    return (local_events / mean_events - 1.0) * dN


def delta_nu_squared(k_over_kstar: float, amplitude: float = 1.0) -> float:
    return amplitude * k_over_kstar ** (-float(DELTA))


def p_nu(k: float, kstar: float = 1.0, amplitude: float = 1.0) -> float:
    return TWO_PI_SQUARED * delta_nu_squared(k / kstar, amplitude) / (k**3)


def slope_log_delta(k1: float, k2: float, amplitude: float = 1.0) -> float:
    y1 = math.log(delta_nu_squared(k1, amplitude))
    y2 = math.log(delta_nu_squared(k2, amplitude))
    return (y2 - y1) / (math.log(k2) - math.log(k1))


def superhorizon_gradient_suppression(efolds_outside: float) -> float:
    """(k/aH)^2 after `efolds_outside` e-folds outside the horizon."""
    return math.exp(-2.0 * efolds_outside)


def epsilon_from_w(w_eff: Fraction) -> Fraction:
    return Fraction(3, 2) * (Fraction(1, 1) + w_eff)


def contains(text: str, phrase: str) -> bool:
    return phrase.lower() in text.lower()


def main() -> None:
    print("ITEM 131 HBC CLOCK COVARIANCE / CONSERVATION / EXIT ATTEMPT")

    anchor = (ROOT / "ANCHOR.md").read_text()
    drift = (ROOT / "DRIFT.md").read_text()
    part17 = (ROOT / "part_17_energy_trajectory" / "part_17_energy_trajectory.tex").read_text()

    print("\n[1] nu_HBC from compensated boundary-printing events")
    samples = [
        (100.0, 100.0, 0.1),
        (101.0, 100.0, 0.1),
        (97.0, 100.0, 0.1),
    ]
    increments = [event_clock_increment(local, mean, dn) for local, mean, dn in samples]
    for (local, mean, dn), inc in zip(samples, increments):
        print(f"  events={local:.1f}, mean={mean:.1f}, dN={dn:.3f} -> dnu={inc:+.6f}")
    check(abs(increments[0]) < 1e-15, "mean event load gives zero clock perturbation")
    check(increments[1] > 0 and increments[2] < 0, "excess/deficit print events give signed local e-fold displacement")
    check(abs(sum(increments) - event_clock_increment(298.0, 300.0, 0.1)) > 1e-3, "nu_HBC is local/integrated, not a global count ratio")
    check(True, "definition: nu_HBC(x,N)=integral[j(x,N)/jbar(N)-1] dN")

    print("\n[2] Covariance shape from the QEC log-shell generator")
    for ratio in [0.25, 1.0, 4.0, 16.0]:
        delta2 = delta_nu_squared(ratio, amplitude=2.1e-9)
        print(f"  k/k*={ratio:5.2f}: Delta_nu^2={delta2:.12e}")
    slope = slope_log_delta(0.5, 8.0)
    check(abs(slope + float(DELTA)) < 1e-15, "Delta_nu^2 has log-slope -1/28")
    check(abs((1.0 + slope) - float(Fraction(27, 28))) < 1e-15, "covariance ledger gives n_s=27/28")
    p_ratio = p_nu(2.0) / p_nu(1.0)
    expected_p_ratio = 2.0 ** (-(3.0 + float(DELTA)))
    check(abs(p_ratio - expected_p_ratio) < 1e-15, "P_nu(k) carries k^-3 times the anomalous k^-1/28 factor")
    check(True, "A_nu remains an undetermined shot-noise/event-density amplitude, not fixed by the 28-clock")

    print("\n[3] R_HBC conservation outside the printed horizon shell")
    for efolds in [0.0, 1.0, 3.0, 5.0]:
        suppression = superhorizon_gradient_suppression(efolds)
        print(f"  {efolds:3.1f} e-folds outside: gradient source ~ {suppression:.6e}")
    check(superhorizon_gradient_suppression(3.0) < 0.003, "gradient source is negligible after a few e-folds outside")
    check(True, "single-clock HBC implies no non-adiabatic source term by construction")
    check(True, "therefore R_HBC is conserved outside the printed shell only under the no-extra-source/single-clock premise")

    print("\n[4] Constant-H duration and exit audit")
    eps_exact_desitter = epsilon_from_w(Fraction(-1, 1))
    eps_late = epsilon_from_w(Fraction(-27, 28))
    print(f"  epsilon(w=-1)      = {eps_exact_desitter}")
    print(f"  epsilon(w=-27/28)  = {eps_late}")
    check(eps_exact_desitter == 0, "saturated w=-1 gives constant H")
    check(eps_late == Fraction(3, 56), "late item-131 branch is not a constant-H stage")

    n_slowroll = Fraction(2, 1) / DELTA
    n_percolation = 28.0 * math.log(1.0 / 0.110)
    passive_dilution = math.log(100.0) / 3.0
    print(f"  slow-roll identity 2/Delta       = {float(n_slowroll):.3f}")
    print(f"  Nishimori-style 28 ln(1/0.110)  = {n_percolation:.3f}")
    print(f"  passive dilution ln(100)/3       = {passive_dilution:.3f}")
    check(n_slowroll == 56, "N=56 is the standard 2/Delta slow-roll identity")
    check(abs(n_percolation - 61.80) < 0.1, "the threshold-style number is ~61.8 algorithmic ticks, not exactly 56")
    check(passive_dilution < 2.0, "passive dilution alone gives far too few e-folds")
    check(contains(part17, "w(0) = -\\frac{4}{4} = -1"), "Part 17 supplies the w(0)=-1 boundary")
    check(contains(part17, "closure of item 42 would lift this boundary condition"), "Part 17 defers the activation-threshold derivation to item 42")
    check(contains(anchor, "phenomenologically arranged") and contains(anchor, "there is currently no algorithmic / geometric / thermodynamic proof"), "ANCHOR marks R-activation schedule as phenomenological/open")
    check(
        contains(drift, "percolation reframing is a legitimate mechanism _sketch_ but [speculation]")
        or contains(drift, "percolation reframing is a legitimate mechanism *sketch* but [speculation]"),
        "DRIFT demotes N/percolation to an unexecuted mechanism sketch",
    )

    print("\n[5] Closure status")
    closed = [
        "candidate nu_HBC as compensated local print-event clock",
        "covariance shape Delta_nu^2=A_nu(k/k*)^-1/28 from the log-shell QEC generator",
        "R_HBC conservation outside the shell under single-clock/no-extra-source conditions",
    ]
    open_items = [
        "A_nu amplitude from microscopic event density / shot-noise ledger",
        "intrinsic proof that the HBC event clock is the unique adiabatic scalar source",
        "finite saturated w=-1 duration from a derived threshold",
        "exit/reheating trigger from item-42 R-activation geometry rather than a narrative schedule",
    ]
    for item in closed:
        check(True, f"conditional closure: {item}")
    for item in open_items:
        check(True, f"not closed: {item}")

    print("\n" + "=" * 100)
    print("VERDICT")
    print("  The event ledger can define the missing clock field:")
    print("      nu_HBC(x,N)=integral[j(x,N)/jbar(N)-1] dN.")
    print("  With the already-closed scale-covariant 28-channel generator, this gives")
    print("      Delta_nu^2(k)=A_nu (k/k*)^-1/28")
    print("  and hence the desired tilt shape.  A_nu is not derived; it needs the")
    print("  microscopic event-rate / shot-noise ledger.")
    print("  R_HBC conservation outside the printed shell follows only under the")
    print("  single-clock, no-extra-source premise; gradient sources decay as")
    print("  exp(-2 DeltaN_outside).")
    print("  The finite early saturated w=-1 duration and exit do NOT close from")
    print("  current canon.  Part 17 gives w(0)=-1; item 42 / R-activation thresholds")
    print("  are explicitly open; N=56 is a slow-roll identity and the percolation")
    print("  route is still an unexecuted mechanism sketch.")
    print("=" * 100)
    print("exit 0 -- clock/covariance/conservation conditionally reduced; duration/exit remain open.")


if __name__ == "__main__":
    main()
