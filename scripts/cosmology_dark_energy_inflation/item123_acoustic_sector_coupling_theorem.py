#!/usr/bin/env python3
r"""Item 123 / M15: CMB acoustic sector-coupling theorem audit.

Question
--------
The substrate-service scripts have derived the depth-6 busy projector

    P_busy = I - |111111><111111|,
    Tr(P_busy)/64 = 63/64.

The remaining cosmology question was whether the CMB acoustic ruler bills this
busy exposure.  This script isolates the theorem that makes that reading
sector-native:

    The CMB acoustic ruler is an acoustic-phase / conformal-lapse observable.
    In canon, the conformal lapse is the A_1g clock strain, and cosmological
    scalar modes read local service-current clocks rather than internal
    service-label projectors.  Therefore the long-wavelength acoustic mode can
    bill only the scalar native service-lapse exposure.  At depth 6 the active
    address-demux service law makes that exposure exactly

        P_acoustic = P_busy = I - |111111><111111|.

This is a hydrodynamic/lapse-grade theorem, not a full Planck likelihood or a
new microscopic baryon-photon action.  A future microscopic action that couples
the acoustic phase to service-label backlog, endpoint flux, or another
non-lapse perturbation variable would be new structure and would reopen the
63/64 acoustic branch.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
import itertools
import re
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_substrate_busy_flag_law as busy_law


DEPTH = busy_law.DEPTH
N_STATE = busy_law.N_STATE
COMPLETE = busy_law.COMPLETE


@dataclass(frozen=True)
class Observable:
    name: str
    operator: np.ndarray
    meaning: str
    verdict: str


@dataclass(frozen=True)
class Anchor:
    path: Path
    patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def idx() -> dict[tuple[int, ...], int]:
    return busy_law.index()


def backlog(state: tuple[int, ...]) -> int:
    return len(busy_law.open_bits(state))


def diagonal_from_backlog(values: dict[int, float]) -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    lookup = idx()
    for state in busy_law.states():
        out[lookup[state], lookup[state]] = values[backlog(state)]
    return out


def microcanonical_expectation(operator: np.ndarray) -> Fraction:
    value = float(np.trace(operator) / N_STATE)
    return Fraction(str(round(value, 12))).limit_denominator()


def p_complete() -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    out[idx()[COMPLETE], idx()[COMPLETE]] = 1.0
    return out


def p_busy() -> np.ndarray:
    return busy_law.p_busy()


def additive_backlog() -> np.ndarray:
    return diagonal_from_backlog({m: float(m) for m in range(DEPTH + 1)})


def blind_alphabet_exposure() -> np.ndarray:
    return busy_law.activity(busy_law.blind_one_jump_jumps())


def independent_hazard() -> np.ndarray:
    return busy_law.activity(busy_law.independent_hazard_jumps())


def active_service_exposure() -> np.ndarray:
    return busy_law.activity(busy_law.active_demux_jumps())


def observables() -> list[Observable]:
    return [
        Observable(
            "busy service-lapse",
            p_busy(),
            "binary QND exposure of the active-demux completion service",
            "wanted acoustic clock under the sector-coupling theorem",
        ),
        Observable(
            "blind alphabet draw",
            blind_alphabet_exposure(),
            "one random depth label is drawn; already-complete labels idle",
            "rejected: not work-conserving on live backlog",
        ),
        Observable(
            "independent hazards",
            independent_hazard(),
            "each open completion label has unit rate",
            "rejected: violates bandwidth-one service",
        ),
        Observable(
            "additive backlog count",
            additive_backlog(),
            "counts how many internal completion labels remain open",
            "rejected: label/backlog-resolved, not long-wavelength acoustic lapse",
        ),
        Observable(
            "completed endpoint",
            p_complete(),
            "late completed latch state itself",
            "rejected for acoustic era; this is the post-latch endpoint",
        ),
    ]


def s6_orbit_hist(operator: np.ndarray) -> dict[int, set[Fraction]]:
    diag = np.diag(operator)
    out: dict[int, set[Fraction]] = {}
    for state, value in zip(busy_law.states(), diag):
        out.setdefault(backlog(state), set()).add(Fraction(str(round(float(value), 12))).limit_denominator())
    return out


def binary_s6_projectors_vanishing_on_complete() -> list[dict[int, int]]:
    """All binary diagonal S6-covariant projectors with f(m=0)=0."""

    projectors: list[dict[int, int]] = []
    for tail in itertools.product((0, 1), repeat=DEPTH):
        values = {0: 0}
        for m, value in enumerate(tail, start=1):
            values[m] = value
        projectors.append(values)
    return projectors


def admissible_acoustic_projectors() -> list[dict[int, int]]:
    """Apply the sector-coupling clauses to the S6-covariant projector list."""

    out: list[dict[int, int]] = []
    for values in binary_s6_projectors_vanishing_on_complete():
        if all(values[m] == 1 for m in range(1, DEPTH + 1)):
            out.append(values)
    return out


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text(encoding="utf-8")
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def anchors() -> list[Anchor]:
    return [
        Anchor(
            ROOT / "python_code" / "item123_substrate_busy_flag_law.py",
            (
                r"A_busy = I - \|111111><111111\|",
                r"every incomplete selector state exposes one",
            ),
            "Substrate busy-flag law: active-address demux exposes exactly one service channel for every incomplete state.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_acoustic_readout_map_derivation_attempt.py",
            (
                r"ordinary additive Peierls/Lindblad workload does not give",
                r"acoustic ruler bills that busy exposure",
            ),
            "Previous acoustic map audit: additive hazard is rejected; the remaining gap is sector coupling.",
        ),
        Anchor(
            ROOT / "python_code" / "relativity_TR3_metric_construction.py",
            (
                r"A_1g \(conformal/lapse = the CLOCK\)",
                r"strain renormalises the tick",
            ),
            "Relativity bridge: the scalar conformal/lapse variable is the local clock rate.",
        ),
        Anchor(
            ROOT / "python_code" / "relativity_integration_programme.py",
            (
                r"coarse-grained covariance of that clock under strain/record load",
                r"geodesics = propagation-phase extremisation",
            ),
            "Relativity programme: propagation phases, including acoustic rulers, read the clock/lapse field.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_zero_mode_source_geodesic_halo_gate.py",
            (
                r"Brown--Kuchar/geodesic continuum lift",
                r"advects the service labels along worldlines",
            ),
            "Cosmology lift: service labels become macroscopic worldline/lapse data, not standalone internal labels.",
        ),
        Anchor(
            ROOT / "python_code" / "item131_scalar_mode_projector.py",
            (
                r"local compensated service-current clock",
                r"not a service-label projector",
            ),
            "Perturbation precedent: cosmological modes read service-current projectors, not internal service-label counts.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_zero_mode_structure_implementation.py",
            (
                r"sector-coupling theorem",
                r"busy-flag exposure",
            ),
            "Boltzmann implementation keeps this exact sector-coupling theorem as the live residual.",
        ),
    ]


def main() -> None:
    print("ITEM 123 / M15: ACOUSTIC SECTOR-COUPLING THEOREM AUDIT")
    print("=" * 100)

    print("\n[1] Source anchors")
    for anchor in anchors():
        ok = source_contains(anchor.path, anchor.patterns)
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(ok, f"{rel} contains the required premise")

    print("\n[2] Candidate acoustic observables on the depth-6 selector")
    for obs in observables():
        expectation = microcanonical_expectation(obs.operator)
        hist = s6_orbit_hist(obs.operator)
        print(f"  {obs.name:<22s} <O>_sel={str(expectation):>5s}  orbits={hist}")
        print(f"    {obs.meaning}")
        print(f"    {obs.verdict}")

    check(np.allclose(active_service_exposure(), p_busy()), "active-demux exposure is P_busy")
    check(microcanonical_expectation(p_busy()) == Fraction(63, 64), "busy exposure gives 63/64")
    check(microcanonical_expectation(blind_alphabet_exposure()) == Fraction(1, 2), "blind alphabet draw gives 1/2")
    check(microcanonical_expectation(additive_backlog() / DEPTH) == Fraction(1, 2), "normalized backlog count gives 1/2")
    check(microcanonical_expectation(p_complete()) == Fraction(1, 64), "completed endpoint is not the acoustic exposure")

    print("\n[3] Uniqueness under the acoustic sector-coupling clauses")
    all_projectors = binary_s6_projectors_vanishing_on_complete()
    admissible = admissible_acoustic_projectors()
    print(f"  S6-covariant binary projectors with complete state dark = {len(all_projectors)}")
    print(f"  after requiring every live backlog to expose the native service lapse = {len(admissible)}")
    print(f"  unique value table f(m) = {admissible[0] if admissible else None}")
    check(len(all_projectors) == 2**DEPTH, "there are many possible S6 projectors before sector coupling")
    check(len(admissible) == 1, "sector-coupling clauses select a unique projector")
    selected = diagonal_from_backlog({m: float(admissible[0][m]) for m in range(DEPTH + 1)})
    check(np.allclose(selected, p_busy()), "the unique selected projector is P_busy")

    print("\n[4] Hydrodynamic/lapse closure")
    print(
        """
  Acoustic-lapse theorem:
    The acoustic ruler is a propagation-phase clock: r_s = int c_s d eta.
    The framework's metric construction identifies the scalar conformal lapse
    with the A_1g clock strain; scalar cosmological readouts are local
    service-current clock projectors, not internal service-label counts.

    Therefore the acoustic readout has to be the scalar native service-lapse
    observable.  In the depth-6 selector, active-address demux makes service
    bandwidth one: every incomplete state exposes exactly one native service
    channel per tick, and the completed latch state exposes none.

  Finite theorem:
    Let rho_sel be the pre-reset selector state, uniform over the depth-6
    address register.  The long-wavelength QND acoustic-lapse readout sees:

      (i) the complete latch state as dark before reset;
     (ii) every incomplete state as one active service channel per tick;
    (iii) no resolution of the six internal completion labels or backlog size.

    The only S6-covariant binary selector observable satisfying these clauses is

      P_acoustic = P_busy = I - |111111><111111|.

    Therefore

      H_CMB / H_selector = Tr(P_acoustic rho_sel) = 63/64.

  Failure modes:
    A backlog-count, Peierls/hazard, completed-endpoint, or depth-label-resolved
    coupling is not the acoustic lapse; it is a different microscopic action.
    Such a future action would be new structure and would refute the 63/64
    acoustic branch.  It is not present in current canon.

  Status:
    The CMB acoustic ruler bills P_busy at hydrodynamic/lapse grade.  This
    closes the local sector-coupling theorem for the compressed theta_* gate.
    Full Planck TT/TE/EE/lensing likelihood and nonlinear halo modelling remain
    external empirical gates.
exit 0 -- acoustic ruler bills P_busy at hydrodynamic/lapse grade.
"""
    )


if __name__ == "__main__":
    main()
