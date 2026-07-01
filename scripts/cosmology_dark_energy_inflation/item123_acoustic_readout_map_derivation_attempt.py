#!/usr/bin/env python3
r"""Item 123 / M15: attempt to derive the acoustic readout map.

Target map
----------
Derive, rather than assign,

    P_acoustic = I - |complete><complete|

on the 64-state depth-6 selector register.

Result
------
There is a clean derivation from the active-address substrate completion
service supplied by item123_substrate_busy_flag_law.py.  If the CMB acoustic
ruler couples to the *binary pre-latch support* of that completion dissipator,
then the support projector is exactly

    supp(sum L^\dagger L) = I - |111111><111111| .

But the ordinary additive Peierls/Lindblad workload does not give that
operator.  It gives a hazard proportional to the number of unlatched bits,

    A_hazard = sum_i (1 - n_i),

whose normalized trace is 1/2, not 63/64.  To obtain the projector itself,
one must add a saturated busy-flag/backlog primitive: each incomplete register
state counts as one acoustically open state regardless of how many completion
links remain.

The remaining boundary is now narrower: the busy-flag law is proved at the
substrate-service level, while the CMB-sector coupling still has to say why
the acoustic ruler bills that busy exposure.
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

import item123_cmb_lock_attempt as lock


DEPTH = 6
N_STATE = 2**DEPTH
COMPLETE = (1,) * DEPTH
PRE_LATCH = Fraction(N_STATE - 1, N_STATE)


@dataclass(frozen=True)
class SourceAnchor:
    path: Path
    patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def basis(depth: int = DEPTH) -> list[tuple[int, ...]]:
    return list(itertools.product((0, 1), repeat=depth))


def index_map(states: list[tuple[int, ...]]) -> dict[tuple[int, ...], int]:
    return {state: i for i, state in enumerate(states)}


def open_bits(state: tuple[int, ...]) -> int:
    return sum(1 - bit for bit in state)


def completion_jumps(*, saturated: bool) -> list[np.ndarray]:
    """Completion-service jumps on the Boolean lattice.

    Unsaturated jumps have unit rate on every unlatched bit.  Saturated jumps
    divide the per-link amplitude by sqrt(number of open bits), so each
    incomplete state has total workload one.  That is the extra busy-flag
    primitive under test.
    """

    states = basis()
    idx = index_map(states)
    jumps: list[np.ndarray] = []
    for state in states:
        m = open_bits(state)
        if m == 0:
            continue
        amplitude = 1.0 / np.sqrt(m) if saturated else 1.0
        for bit in range(DEPTH):
            if state[bit] == 1:
                continue
            target = list(state)
            target[bit] = 1
            op = np.zeros((N_STATE, N_STATE), dtype=float)
            op[idx[tuple(target)], idx[state]] = amplitude
            jumps.append(op)
    return jumps


def activity_operator(jumps: list[np.ndarray]) -> np.ndarray:
    out = np.zeros((N_STATE, N_STATE), dtype=float)
    for jump in jumps:
        out += jump.T @ jump
    return out


def support_projector(operator: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    evals, evecs = np.linalg.eigh(operator)
    keep = evals > tol
    return (evecs[:, keep] @ evecs[:, keep].T).round(12)


def ideal_p_acoustic() -> np.ndarray:
    states = basis()
    idx = index_map(states)
    p = np.eye(N_STATE)
    p[idx[COMPLETE], idx[COMPLETE]] = 0.0
    return p


def acoustic_numerics() -> tuple[float, float, float, float, float]:
    omega_b, _omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    selector_theta = lock.theta100_framework(lock.H0_SELECTOR, omega_b, omega_dark)
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_pre = lock.H0_SELECTOR * float(PRE_LATCH)
    theta_pre = lock.theta100_framework(h0_pre, omega_b, omega_dark)
    return target, selector_theta, h0_root, h0_pre, theta_pre


def eigenvalue_hist(diagonal: np.ndarray) -> dict[float, int]:
    hist: dict[float, int] = {}
    for value in diagonal:
        key = float(round(value, 12))
        hist[key] = hist.get(key, 0) + 1
    return dict(sorted(hist.items()))


def local_rate_solution() -> list[tuple[int, Fraction]]:
    """Per-link rate needed to make total workload one at each backlog m."""

    return [(m, Fraction(1, m)) for m in range(1, DEPTH + 1)]


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text(encoding="utf-8")
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def anchors() -> list[SourceAnchor]:
    return [
        SourceAnchor(
            ROOT / "python_code" / "item123_acoustic_prelatch_operator_audit.py",
            (
                r"rank-63 projector is not S64-invariant",
                r"adds new physical structure",
            ),
            "The prior audit already separates the mathematical projector from a licensed selector repair.",
        ),
        SourceAnchor(
            ROOT / "python_code" / "item123_cmb_lock_attempt.py",
            (
                r"THEOREM TARGET ONLY",
                r"no operator proving",
            ),
            "The CMB lock attempt keeps 63/64 as a target rather than a theorem.",
        ),
        SourceAnchor(
            ROOT / "python_code" / "dressed_alpha_bridge_latch_audit.py",
            (
                r"Peierls bridge readout has an analogous local\s+one-bit latch",
                r"Peierls bridge readout does not inherit",
            ),
            "The bridge-local Peierls latch route is already refuted for the alpha bridge.",
        ),
        SourceAnchor(
            ROOT / "python_code" / "dressed_alpha_monitor_web_cptp_audit.py",
            (
                r"full Peierls/Wilson web width operator",
                r"promising, not closed",
            ),
            "The best CPTP bridge/web precedent is promising but still not a licensing theorem.",
        ),
        SourceAnchor(
            ROOT / "python_code" / "item123_substrate_busy_flag_law.py",
            (
                r"substrate busy-flag law proved from active-address demux",
                r"every incomplete selector state exposes one",
            ),
            "The substrate-service busy-flag law has now been supplied by active-address demux.",
        ),
    ]


def main() -> None:
    print("ITEM 123 / M15: ACOUSTIC READOUT MAP DERIVATION ATTEMPT")
    print("=" * 100)
    print(f"CAMB version: {getattr(lock.camb, '__version__', 'unknown')}")

    print("\n[1] Acoustic target still wants trace 63/64")
    target, selector_theta, h0_root, h0_pre, theta_pre = acoustic_numerics()
    print(f"  Planck/LCDM reference 100 theta_* = {target:.6f}")
    print(f"  selector H0={lock.H0_SELECTOR:.3f}: theta = {selector_theta:.6f} ({(selector_theta / target - 1.0) * 100:+.3f}%)")
    print(f"  CAMB acoustic root H0             = {h0_root:.6f}")
    print(f"  trace-factor H0=(63/64)selector   = {h0_pre:.6f}")
    print(f"  theta at trace-factor H0          = {theta_pre:.6f} ({(theta_pre / target - 1.0) * 100:+.4f}%)")
    check(abs(theta_pre / target - 1.0) < 2.0e-5, "63/64 remains the acoustic-scale target")
    check(abs(h0_root / h0_pre - 1.0) < 1.0e-4, "CAMB root is essentially the trace-factor candidate")

    print("\n[2] Completion-service Lindblad support")
    p_target = ideal_p_acoustic()
    jumps = completion_jumps(saturated=False)
    hazard = activity_operator(jumps)
    support = support_projector(hazard)
    diag_hazard = np.diag(hazard)
    print(f"  completion register states          = {N_STATE}")
    print(f"  local completion jumps              = {len(jumps)}")
    print(f"  unsaturated hazard eigenvalue hist  = {eigenvalue_hist(diag_hazard)}")
    print(f"  Tr support                          = {np.trace(support):.0f}")
    print(f"  Tr support / 64                     = {Fraction(int(np.trace(support)), N_STATE)}")
    print(f"  Tr hazard / (6*64)                  = {np.trace(hazard) / (DEPTH * N_STATE):.6f}")
    check(np.allclose(support, p_target), "support of the completion dissipator is exactly I-|complete><complete|")
    check(Fraction(int(np.trace(support)), N_STATE) == PRE_LATCH, "support projector has normalized trace 63/64")
    check(abs(np.trace(hazard) / (DEPTH * N_STATE) - 0.5) < 1e-12, "additive hazard does not give 63/64")

    print("\n[3] Saturated busy-flag primitive")
    saturated_hazard = activity_operator(completion_jumps(saturated=True))
    diag_sat = np.diag(saturated_hazard)
    print(f"  saturated hazard eigenvalue hist    = {eigenvalue_hist(diag_sat)}")
    print(f"  Tr saturated hazard                 = {np.trace(saturated_hazard):.0f}")
    print(f"  Tr saturated hazard / 64            = {Fraction(int(np.trace(saturated_hazard)), N_STATE)}")
    check(np.allclose(saturated_hazard, p_target), "backlog-normalized busy flag gives P_acoustic exactly")
    print("  This is the desired map, but the normalization is global:")
    print(f"  required per-link rate by backlog m = {local_rate_solution()}")
    check(len(set(rate for _m, rate in local_rate_solution())) == DEPTH, "exact projector needs backlog-dependent per-link rates 1/m")

    print("\n[4] Bridge/web and substrate licensing")
    for anchor in anchors():
        ok = source_contains(anchor.path, anchor.patterns)
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(ok, f"{rel} still records the required guardrail")

    print("\n[5] Derivation status")
    print(
        """
  What is derived:
    The monotone completion Lindbladian has a unique dark/absorbing state,
    |complete>.  Therefore the support of its activity operator is exactly

        supp(A_completion) = I - |complete><complete|.

    If the acoustic ruler is a QND busy-flag readout of this support, then
    P_acoustic is derived and Tr(P_acoustic)/64 = 63/64.

  What is not derived from the existing bridge/web primitive:
    The Peierls-current or ordinary local Lindblad workload reads the hazard
    A_completion = sum_i (1-n_i), not its support.  That hazard weights states
    by the number of unlatched bits and has normalized trace 1/2.  It is a
    different observable.

  What is now supplied:
    item123_substrate_busy_flag_law.py proves the 1/m backlog normalization
    from the active-address demux service primitive:

        every incomplete depth-6 selector state exposes one substrate service
        channel per clock tick, independent of remaining backlog m;
        the all-complete state exposes none.

    Therefore the completion-service activity operator is P_acoustic.

  Verdict:
    P_acoustic is derived as a substrate-service readout from active-address
    demux, while the additive Peierls/hazard reading remains rejected.  The
    remaining explicit boundary is the cosmology-sector coupling: why the CMB
    acoustic ruler bills this busy service exposure.
exit 0 -- acoustic map derived at substrate-service level; acoustic sector coupling remains explicit.
"""
    )


if __name__ == "__main__":
    main()
