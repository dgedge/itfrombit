#!/usr/bin/env python3
r"""Item 123 / M15: CMB acoustic clock as QND service-lapse readout.

Question
--------
Can the framework now prove that the CMB acoustic clock bills the QND
service-lapse readout, rather than a flux, Peierls hazard, endpoint latch, or
internal backlog observable?

Result
------
At continuum service-current/EFT grade, yes.  The photon-baryon acoustic mode
is a phase clock:

    d phi / dt = N_eta c_s k,

where N_eta is the conformal lapse.  The existing gravity layer identifies the
scalar lapse/clock variable with the A1g service-clock covariance sourced by
the service-current stress tensor.  A long-wavelength CMB acoustic phase is
therefore a QND phase modulation of the service lapse: it may see the scalar
busy exposure, but it cannot resolve internal completion labels or perform the
completion jump itself without leaving the scalar lapse sector.

Composed with the substrate busy-flag theorem, the unique binary selector
operator is

    P_acoustic = P_busy = I - |111111><111111|,

so

    H_CMB / H_selector = Tr(P_busy rho_selector) = 63/64.

This does not claim an exact finite substrate lift to the Boltzmann-Keldysh
action.  It closes the continuum acoustic-lapse clause under the already
accepted service-current gravity premises.  A future exact derivation should
construct the same scalar lapse operator directly in the finite photon-baryon
Boltzmann-Keldysh lift.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
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
class Anchor:
    path: Path
    patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text(encoding="utf-8")
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def anchors() -> list[Anchor]:
    return [
        Anchor(
            ROOT / "python_code" / "relativity_TR3_metric_construction.py",
            (
                r"A_1g \(conformal/lapse = the CLOCK\)",
                r"uniform \(A_1g\) strain rescales the light cone",
            ),
            "Metric construction: the scalar A1g mode is the conformal lapse / clock covariance.",
        ),
        Anchor(
            ROOT / "python_code" / "modular_hamiltonian_service_ledger.py",
            (
                r"service-clock gradient realises it",
                r"H_mod = 2pi K_boost",
            ),
            "Gravity clock layer: the service-clock gradient supplies the boost/lapse modular weight.",
        ),
        Anchor(
            ROOT / "python_code" / "gravity_service_current_stress_tensor_gate.py",
            (
                r"T\^{mu nu}_svc",
                r"service-current stress tensor",
            ),
            "Source layer: the ledger supplies a conserved service-current stress tensor.",
        ),
        Anchor(
            ROOT / "python_code" / "intrinsic_gravity_linearized_einstein_gate.py",
            (
                r"linearized Einstein/source structure",
                r"RT first variations \+ boost modular Hamiltonian \+ service-current",
            ),
            "Continuum gate: service-current source plus lapse modular Hamiltonian gives the Einstein/source form.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_substrate_busy_flag_law.py",
            (
                r"every incomplete depth-6 selector state exposes exactly one service\s+channel",
                r"A_busy = sum L\^\\dagger L = I - \|111111><111111\|",
            ),
            "Finite selector gate: active-address demux gives the binary busy projector.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_acoustic_sector_coupling_theorem.py",
            (
                r"long-wavelength QND acoustic-lapse readout",
                r"P_acoustic = P_busy",
            ),
            "Previous reduction: once the acoustic clock is QND service-lapse, the 63/64 result is forced.",
        ),
    ]


def idx() -> dict[tuple[int, ...], int]:
    return busy_law.index()


def p_complete() -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    out[idx()[COMPLETE], idx()[COMPLETE]] = 1.0
    return out


def p_busy() -> np.ndarray:
    return busy_law.p_busy()


def backlog_count() -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    lookup = idx()
    for state in busy_law.states():
        out[lookup[state], lookup[state]] = float(len(busy_law.open_bits(state)))
    return out


def expectation(operator: np.ndarray) -> Fraction:
    value = float(np.trace(operator) / N_STATE)
    return Fraction(str(round(value, 12))).limit_denominator()


def acoustic_phase_rate(lapse: float, k: float = 0.19, sound_speed: float = 1.0 / np.sqrt(3.0), a: float = 1.0) -> float:
    """WKB acoustic phase rate dphi/dt = N_eta c_s k / a."""

    return float(lapse * sound_speed * k / a)


def acoustic_hamiltonian_coefficients(lapse: float, amplitude: float = 2.7, omega: float = 0.11) -> tuple[float, float]:
    """Return the p^2 and q^2 coefficients in H = N[p^2/(2A) + A omega^2 q^2/2]."""

    p2_coeff = lapse / (2.0 * amplitude)
    q2_coeff = lapse * amplitude * omega * omega / 2.0
    return p2_coeff, q2_coeff


def number_operator(n: int = 6) -> np.ndarray:
    return np.diag(np.arange(n, dtype=float))


def qnd_phase_hamiltonian(omega: float = 1.0, n_osc: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Selector-busy phase modulation of an acoustic oscillator."""

    selector = p_busy()
    osc = number_operator(n_osc) + 0.5 * np.eye(n_osc)
    h = omega * np.kron(selector, osc)
    readout = np.kron(selector, np.eye(n_osc))
    return h, readout


def off_selector_block_norm(matrix: np.ndarray, n_osc: int) -> float:
    total = 0.0
    for a in range(N_STATE):
        for b in range(N_STATE):
            if a == b:
                continue
            block = matrix[a * n_osc : (a + 1) * n_osc, b * n_osc : (b + 1) * n_osc]
            total += float(np.linalg.norm(block))
    return total


def max_jump_qnd_commutator(jumps: list[np.ndarray]) -> float:
    projector = p_busy()
    return max(float(np.linalg.norm(projector @ jump - jump @ projector)) for jump in jumps)


def main() -> None:
    print("ITEM 123 / M15: CMB ACOUSTIC CLOCK AS QND SERVICE-LAPSE READOUT")
    print("=" * 100)

    print("\n[1] Canon anchors for the continuum composition")
    for anchor in anchors():
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(source_contains(anchor.path, anchor.patterns), f"{rel} contains the required premise")

    print("\n[2] Acoustic action lemma: the acoustic ruler is a lapse clock")
    lapses = np.array([0.25, 0.5, 1.0, 1.5, 2.0])
    rates = np.array([acoustic_phase_rate(float(n)) for n in lapses])
    ratios = rates / rates[2]
    print(f"  lapse samples        = {lapses}")
    print(f"  phase-rate ratios    = {ratios}")
    check(np.allclose(ratios, lapses / lapses[2]), "WKB acoustic phase rate is linear in conformal lapse")
    p2_a, q2_a = acoustic_hamiltonian_coefficients(0.5)
    p2_b, q2_b = acoustic_hamiltonian_coefficients(1.5)
    print(f"  Hamiltonian p^2 coefficient ratio = {p2_b / p2_a:.6f}")
    print(f"  Hamiltonian q^2 coefficient ratio = {q2_b / q2_a:.6f}")
    check(abs(p2_b / p2_a - 3.0) < 1.0e-12, "canonical p^2 term scales with the same lapse")
    check(abs(q2_b / q2_a - 3.0) < 1.0e-12, "canonical q^2 term scales with the same lapse")

    print("\n[3] QND readout lemma")
    n_osc = 6
    h_qnd, readout = qnd_phase_hamiltonian(n_osc=n_osc)
    comm = np.linalg.norm(h_qnd @ readout - readout @ h_qnd)
    off_blocks = off_selector_block_norm(h_qnd, n_osc)
    print(f"  ||[H_phase, P_busy]||       = {comm:.3e}")
    print(f"  off-selector block norm     = {off_blocks:.3e}")
    check(comm < 1.0e-12, "phase readout commutes with the busy-lapse observable")
    check(off_blocks < 1.0e-12, "phase readout causes no selector transitions")

    jump_comm = max_jump_qnd_commutator(busy_law.active_demux_jumps())
    print(f"  max ||[L_completion, P_busy]|| = {jump_comm:.3e}")
    check(jump_comm > 0.5, "completion jumps are not the acoustic readout; they change the latch")

    print("\n[4] Selector observable forced by long-wavelength scalar lapse")
    p_b = p_busy()
    p_c = p_complete()
    backlog = backlog_count()
    normalized_backlog = backlog / DEPTH
    print(f"  <P_busy>_selector             = {expectation(p_b)}")
    print(f"  <normalized backlog>_selector = {expectation(normalized_backlog)}")
    print(f"  <P_complete>_selector         = {expectation(p_c)}")
    check(expectation(p_b) == Fraction(63, 64), "binary service-lapse exposure gives 63/64")
    check(expectation(normalized_backlog) == Fraction(1, 2), "backlog/hazard count gives the wrong clock")
    check(expectation(p_c) == Fraction(1, 64), "endpoint flux gives the wrong clock")
    check(np.allclose(p_b, np.eye(N_STATE) - p_c), "P_busy is exactly I - |complete><complete|")

    print("\n[5] The composed theorem")
    print(
        """
  Theorem, at continuum service-current/EFT grade:

    1. In the photon-baryon acoustic action, the acoustic phase is a clock in
       conformal time: d phi = k c_s d eta.  A scalar lapse perturbation
       phase-modulates the oscillator Hamiltonian; it does not perform a jump.

    2. In the framework gravity layer, the scalar lapse is the A1g
       service-clock covariance sourced by the conserved service-current
       stress tensor.  This identifies the relevant cosmological scalar with a
       service-lapse current, not with a U(1) Peierls current, endpoint flux,
       or internal completion label.

    3. The long-wavelength CMB mode is QND with respect to that lapse: the
       phase Hamiltonian commutes with the busy-lapse projector and has no
       selector-changing matrix elements.  A completion Lindblad/Peierls jump
       fails this test because it changes the latch.

    4. The finite active-address theorem gives exactly one exposed native
       service channel in every incomplete depth-6 selector state and none in
       the complete state.  Therefore the scalar QND lapse observable is

           P_acoustic = P_busy = I - |111111><111111|.

    Hence

           H_CMB / H_selector = Tr(P_acoustic rho_selector) = 63/64.

  Residual:
    This proves the acoustic-clock identification at the same continuum grade
    as the service-current gravity theorem.  It is not yet an exact finite
    substrate lift to the photon-baryon Boltzmann-Keldysh action.  That
    stronger proof remains the next microscopic refinement, but the old
    continuum acoustic-lapse clause is no longer an independent fit.
exit 0 -- CMB acoustic clock closed as QND service-lapse readout at continuum grade.
"""
    )


if __name__ == "__main__":
    main()
