#!/usr/bin/env python3
r"""Item 123 / M15: substrate busy-flag law for the acoustic pre-latch.

Goal
----
Prove the service law needed by the acoustic readout map:

    every incomplete depth-6 selector state exposes exactly one service
    channel per substrate tick, while the complete state exposes none.

This is stronger than the old blind one-jump scheduler.  A blind scheduler
chooses from the whole alphabet and idles when it hits an inactive address,
giving exposure m/6 for m open completion bits.  The required law is the
active-address demux law: the decoder tracks the active set and the finite
bandwidth service acts on one true active address per tick.  S6 covariance then
forces uniform weight 1/m over the m open bits.

The resulting Lindblad/Kraus activity operator is

    A_busy = sum L^\dagger L = I - |111111><111111|,

so the normalized activity trace is 63/64.

Status
------
This script proves the busy-flag law from the active-address demux service
primitive already selected in the register-handoff audits.  It does not, by
itself, prove the final sector-coupling statement that the CMB acoustic ruler
must bill this service exposure rather than another perturbation variable.
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


DEPTH = 6
N_STATE = 2**DEPTH
COMPLETE = (1,) * DEPTH


@dataclass(frozen=True)
class Anchor:
    path: Path
    patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def states() -> list[tuple[int, ...]]:
    return list(itertools.product((0, 1), repeat=DEPTH))


def index() -> dict[tuple[int, ...], int]:
    return {state: i for i, state in enumerate(states())}


def open_bits(state: tuple[int, ...]) -> list[int]:
    return [i for i, bit in enumerate(state) if bit == 0]


def complete_one(state: tuple[int, ...], bit: int) -> tuple[int, ...]:
    out = list(state)
    out[bit] = 1
    return tuple(out)


def active_demux_jumps() -> list[np.ndarray]:
    """One true active completion address per tick.

    From a state with m open bits, the active demux selects one of those m
    active addresses with probability/rate 1/m.  The jump amplitude is therefore
    sqrt(1/m), making the total outgoing activity exactly one.
    """

    idx = index()
    jumps: list[np.ndarray] = []
    for state in states():
        active = open_bits(state)
        if not active:
            continue
        amp = 1.0 / np.sqrt(len(active))
        for bit in active:
            op = np.zeros((N_STATE, N_STATE), dtype=float)
            op[idx[complete_one(state, bit)], idx[state]] = amp
            jumps.append(op)
    return jumps


def blind_one_jump_jumps() -> list[np.ndarray]:
    """Blind alphabet draw over all six completion labels.

    If the selected label is already complete, the service action is idle.  The
    non-idle activity from a state with m open bits is m/6.
    """

    idx = index()
    jumps: list[np.ndarray] = []
    amp = 1.0 / np.sqrt(DEPTH)
    for state in states():
        for bit in open_bits(state):
            op = np.zeros((N_STATE, N_STATE), dtype=float)
            op[idx[complete_one(state, bit)], idx[state]] = amp
            jumps.append(op)
    return jumps


def independent_hazard_jumps() -> list[np.ndarray]:
    """Independent local hazards, one unit rate per open completion bit."""

    idx = index()
    jumps: list[np.ndarray] = []
    for state in states():
        for bit in open_bits(state):
            op = np.zeros((N_STATE, N_STATE), dtype=float)
            op[idx[complete_one(state, bit)], idx[state]] = 1.0
            jumps.append(op)
    return jumps


def activity(jumps: list[np.ndarray]) -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    for jump in jumps:
        out += jump.T @ jump
    return out


def p_busy() -> np.ndarray:
    idx = index()
    out = np.eye(N_STATE)
    out[idx[COMPLETE], idx[COMPLETE]] = 0.0
    return out


def exposure_by_backlog(operator: np.ndarray) -> dict[int, set[Fraction]]:
    diag = np.diag(operator)
    out: dict[int, set[Fraction]] = {}
    for state, value in zip(states(), diag):
        m = len(open_bits(state))
        out.setdefault(m, set()).add(Fraction(str(round(float(value), 12))))
    return out


def transition_matrix_active_demux() -> np.ndarray:
    """Discrete one-tick stochastic matrix for the active-demux service."""

    idx = index()
    mat = np.zeros((N_STATE, N_STATE), dtype=float)
    for state in states():
        active = open_bits(state)
        col = idx[state]
        if not active:
            mat[col, col] = 1.0
            continue
        share = 1.0 / len(active)
        for bit in active:
            mat[idx[complete_one(state, bit)], col] += share
    return mat


def permute_state(state: tuple[int, ...], perm: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(state[i] for i in perm)


def is_s6_covariant_exposure(operator: np.ndarray) -> bool:
    diag = np.diag(operator)
    idx = index()
    for perm in itertools.permutations(range(DEPTH)):
        for state in states():
            if abs(diag[idx[state]] - diag[idx[permute_state(state, perm)]]) > 1e-12:
                return False
    return True


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text(encoding="utf-8")
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def anchors() -> list[Anchor]:
    return [
        Anchor(
            ROOT / "python_code" / "item131_one_jump_premise.py",
            (
                r"one channel serviced per\s+tick",
                r"conv\{single-channel resets\}",
                r"finite-bandwidth",
            ),
            "Finite bandwidth gives at most one service action per tick.",
        ),
        Anchor(
            ROOT / "python_code" / "record_content_from_syndrome.py",
            (
                r"QND ledger's increment IS a vertex label",
                r"increment alphabet is 8",
            ),
            "QND increments create tracked active addresses rather than raw snapshots only.",
        ),
        Anchor(
            ROOT / "python_code" / "register_handoff_demux_decision_audit.py",
            (
                r"active-address demultiplexer is operationally selected",
                r"one true active address\s+is serviced",
                r"not on a blind channel",
            ),
            "The service form selected by the register handoff is active-address demux.",
        ),
        Anchor(
            ROOT / "python_code" / "alpha0_count_rate_theorem.py",
            (
                r"one label per interrogation",
                r"one elementary monitor interrogation per cell tick",
            ),
            "The native substrate tick supplies one monitored service interrogation.",
        ),
        Anchor(
            ROOT / "python_code" / "item123_acoustic_readout_map_derivation_attempt.py",
            (
                r"support of the completion dissipator is exactly I-\|complete><complete\|",
                r"required per-link rate by backlog m",
            ),
            "The preceding acoustic audit identified the missing 1/m backlog law.",
        ),
    ]


def main() -> None:
    print("ITEM 123 / M15: SUBSTRATE BUSY-FLAG LAW")
    print("=" * 96)

    print("\n[1] Canon service premises")
    for anchor in anchors():
        ok = source_contains(anchor.path, anchor.patterns)
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(ok, f"{rel} contains the required premise")

    print("\n[2] Three possible service laws on the depth-6 completion register")
    laws = [
        ("active-address demux", active_demux_jumps()),
        ("blind one-jump alphabet", blind_one_jump_jumps()),
        ("independent local hazards", independent_hazard_jumps()),
    ]
    for name, jumps in laws:
        op = activity(jumps)
        exposure = exposure_by_backlog(op)
        trace = Fraction(str(round(float(np.trace(op)), 12)))
        print(f"  {name:<27} jumps={len(jumps):3d} trace={trace} exposure={exposure}")

    busy = activity(active_demux_jumps())
    blind = activity(blind_one_jump_jumps())
    independent = activity(independent_hazard_jumps())
    target = p_busy()

    check(np.allclose(busy, target), "active-address demux activity is exactly the busy projector")
    check(not np.allclose(blind, target), "blind one-jump idles at partial backlog and fails the busy law")
    check(not np.allclose(independent, target), "independent hazards count open bits and fail the busy law")

    print("\n[3] Busy law theorem")
    mat = transition_matrix_active_demux()
    col_sums = mat.sum(axis=0)
    diag_busy = np.diag(busy)
    print(f"  transition matrix column-sum range = {col_sums.min():.6f} .. {col_sums.max():.6f}")
    print(f"  busy exposure values               = {sorted(set(Fraction(str(round(float(x), 12))) for x in diag_busy))}")
    print(f"  Tr(P_busy)                         = {int(np.trace(busy))}")
    print(f"  Tr(P_busy)/64                      = {Fraction(int(np.trace(busy)), N_STATE)}")
    check(np.allclose(col_sums, 1.0), "active-demux service is a CPTP/stochastic one-tick channel")
    check(Fraction(int(np.trace(busy)), N_STATE) == Fraction(63, 64), "busy projector has normalized trace 63/64")
    check(is_s6_covariant_exposure(busy), "busy exposure is S6-covariant")
    check(diag_busy[index()[COMPLETE]] == 0.0, "complete state exposes no service channel")
    check(
        all(abs(diag_busy[index()[state]] - 1.0) < 1e-12 for state in states() if state != COMPLETE),
        "every incomplete state exposes exactly one service channel",
    )

    print("\n[4] Why the 1/m local rate is forced")
    print("  Let B(x) be the open-bit set in state x and m=|B(x)|.")
    print("  Finite bandwidth: total non-idle service per tick is at most one.")
    print("  Active-address demux/no-idle: if m>0, the service acts on a true active address.")
    print("  S6 covariance: all active completion addresses in B(x) have equal weight.")
    print("  Therefore each open bit has rate 1/m, and the total exposure is m*(1/m)=1.")
    for m in range(1, DEPTH + 1):
        print(f"    m={m}: per-open-bit rate={Fraction(1, m)}, total={m}*{Fraction(1, m)}=1")

    print("\n[5] Status")
    print(
        """
  Proved at the substrate-service level:
    The active-address demux primitive gives a work-conserving, bandwidth-one
    completion service.  On a depth-6 selector register, its activity operator
    is exactly

        A_busy = I - |111111><111111|.

    Hence every incomplete selector state exposes one, and only one, native
    service channel per substrate tick; the complete state exposes none.

  Rejected alternatives:
    A blind one-jump alphabet gives exposure m/6.  Independent local hazards
    give exposure m.  Both are the wrong observable for the 63/64 acoustic map.

  Remaining boundary:
    This proves the busy-flag law as a substrate service statement.  The final
    cosmology-sector sentence still has to say why the CMB acoustic ruler bills
    this busy service exposure.  If that sector-coupling is accepted, the
    acoustic readout map follows without an extra scalar fit.
exit 0 -- substrate busy-flag law proved from active-address demux; acoustic sector coupling remains explicit.
"""
    )


if __name__ == "__main__":
    main()
