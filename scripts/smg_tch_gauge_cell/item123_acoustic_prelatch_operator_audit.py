#!/usr/bin/env python3
r"""Item 123 / M15: acoustic pre-latch operator audit.

Question
--------
Can the CMB selector repair

    H_CMB = (63/64) H_selector

be upgraded from a numerical pre-latch target into an operator statement?

This audit tests the strongest local construction available:

    P_pre = I - |111111><111111|

on a depth-6 binary completion register.  It has trace 63 on the 64 basis
states, is a genuine projector, and is invariant under permutations of the six
completion coordinates.  If this were an already-licensed acoustic readout
register, it would implement the desired normalized trace factor exactly.

Verdict
-------
The operator can be constructed mathematically only after adding a labelled
six-bit completion register and declaring the all-complete latch state
acoustically invisible.  The current Canon does not yet license that register
as the CMB acoustic ruler.  Under the existing unlabelled S_64 depth-slot
ledger, a rank-63 invariant projector is impossible.  Therefore the
selector-lock repair is rejected as a current theorem and kept as a sharp
operator target.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import itertools
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))

import item123_cmb_lock_attempt as lock


DEPTH = 6
N_SLOT = 2**DEPTH
PRE_LATCH_FACTOR = Fraction(N_SLOT - 1, N_SLOT)
COMPLETE_STATE = (1,) * DEPTH


@dataclass(frozen=True)
class CanonAnchor:
    path: Path
    required_patterns: tuple[str, ...]
    meaning: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def depth_register_basis(depth: int = DEPTH) -> list[tuple[int, ...]]:
    return list(itertools.product((0, 1), repeat=depth))


def pre_latch_mask(state: tuple[int, ...]) -> int:
    return 0 if state == COMPLETE_STATE else 1


def pre_latch_diagonal() -> list[int]:
    return [pre_latch_mask(state) for state in depth_register_basis()]


def permute_state(state: tuple[int, ...], perm: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(state[i] for i in perm)


def is_s6_invariant() -> bool:
    for perm in itertools.permutations(range(DEPTH)):
        for state in depth_register_basis():
            if pre_latch_mask(permute_state(state, perm)) != pre_latch_mask(state):
                return False
    return True


def s64_invariant_ranks() -> tuple[Fraction, Fraction]:
    """Ranks of diagonal projectors invariant under the unlabelled S_64 action.

    A full permutation of the 64 slot labels is transitive.  Any diagonal
    projector invariant under it must assign the same eigenvalue to every slot,
    so only the zero and identity projectors survive.
    """

    return Fraction(0, 1), Fraction(1, 1)


def acoustic_numerics() -> tuple[float, float, float, float, float]:
    omega_b, _omega_nur, omega_dark = lock.framework_matter_budget()
    target = lock.theta100_planck()
    selector_theta = lock.theta100_framework(lock.H0_SELECTOR, omega_b, omega_dark)
    h0_root = lock.find_h0_root(omega_b, omega_dark, target)
    h0_pre = lock.H0_SELECTOR * float(PRE_LATCH_FACTOR)
    theta_pre = lock.theta100_framework(h0_pre, omega_b, omega_dark)
    return target, selector_theta, h0_root, h0_pre, theta_pre


def source_contains(path: Path, patterns: tuple[str, ...]) -> bool:
    text = path.read_text()
    return all(re.search(pattern, text, flags=re.IGNORECASE | re.DOTALL) for pattern in patterns)


def canon_anchors() -> list[CanonAnchor]:
    return [
        CanonAnchor(
            ROOT / "python_code" / "item123_cmb_lock_attempt.py",
            (
                r"THEOREM TARGET ONLY",
                r"no operator proving",
                r"acoustic propagation omits exactly one depth-6 slot",
            ),
            "CMB lock attempt records 63/64 as a target, not a lock.",
        ),
        CanonAnchor(
            ROOT / "python_code" / "cmb_dark_reservoir_status.py",
            (
                r"63/64 acoustic pre-latch projection",
                r"no current service-ledger\s+theorem proves",
                r"one pre-latch depth-6 slot less",
            ),
            "Dark-reservoir status keeps the acoustic pre-latch as a residual gate.",
        ),
        CanonAnchor(
            ROOT / "python_code" / "item123_cmb_selector_halo_accounting_audit.py",
            (
                r"symmetric depth-6 ledger",
                r"current symmetric ledger cannot derive a rank-63 projector",
            ),
            "The consolidated selector/halo audit rejects derivation from symmetric slots.",
        ),
        CanonAnchor(
            ROOT / "python_code" / "register_handoff_demux_decision_audit.py",
            (
                r"edge-sweep/latch\s+dwell candidate is dead",
                r"keep tick 13 demoted",
            ),
            "The older latch-dwell route is not available as a canon shortcut.",
        ),
        CanonAnchor(
            ROOT / "python_code" / "rho_lambda_129_boundary_derivation_attempt.py",
            (
                r"No current term .*licenses the extra '\+1' leg",
                r"129/128 mean shift is NOT derived",
            ),
            "Analogous finite-depth boundary numerology is rejected without a licensed leg.",
        ),
    ]


def main() -> None:
    print("ITEM 123 / M15: ACOUSTIC PRE-LATCH OPERATOR AUDIT")
    print("=" * 96)
    print(f"CAMB version: {getattr(lock.camb, '__version__', 'unknown')}")

    print("\n[1] Acoustic scale target")
    target, selector_theta, h0_root, h0_pre, theta_pre = acoustic_numerics()
    print(f"  Planck/LCDM reference 100 theta_* = {target:.6f}")
    print(
        f"  selector H0={lock.H0_SELECTOR:.3f}: theta = "
        f"{selector_theta:.6f} ({(selector_theta / target - 1.0) * 100:+.3f}%)"
    )
    print(f"  CAMB acoustic root H0             = {h0_root:.6f}")
    print(f"  trace-factor H0=(63/64)selector   = {h0_pre:.6f}")
    print(
        f"  theta at trace-factor H0          = {theta_pre:.6f} "
        f"({(theta_pre / target - 1.0) * 100:+.4f}%)"
    )
    print(f"  root / trace-factor H0            = {h0_root / h0_pre:.8f}")
    check(abs(selector_theta / target - 1.0) > 0.002, "selector-H0 acoustic offset is real")
    check(abs(theta_pre / target - 1.0) < 2.0e-5, "63/64 trace factor numerically closes theta_*")
    check(abs(h0_root / h0_pre - 1.0) < 1.0e-4, "CAMB root is essentially the trace-factor candidate")

    print("\n[2] Strongest constructible pre-latch projector")
    diagonal = pre_latch_diagonal()
    trace = sum(diagonal)
    rank_fraction = Fraction(trace, len(diagonal))
    excluded = [state for state in depth_register_basis() if pre_latch_mask(state) == 0]
    print("  construction: P_pre = I - |111111><111111| on (C^2)^{tensor 6}")
    print(f"  basis size                         = {len(diagonal)}")
    print(f"  excluded state                     = {excluded[0]}")
    print(f"  trace(P_pre)                       = {trace}")
    print(f"  normalized trace                   = {rank_fraction}")
    check(len(diagonal) == N_SLOT, "depth-6 binary register has 64 basis states")
    check(set(diagonal) <= {0, 1}, "P_pre is diagonal with projector eigenvalues")
    check(all(x * x == x for x in diagonal), "P_pre is idempotent")
    check(trace == 63, "P_pre has trace 63")
    check(rank_fraction == PRE_LATCH_FACTOR, "normalized trace is exactly 63/64")
    check(is_s6_invariant(), "P_pre is invariant under the S6 coordinate symmetry")

    print("\n[3] Why this is not yet the current Canon")
    zero_rank, full_rank = s64_invariant_ranks()
    print(f"  unlabelled S_64 invariant ranks    = {zero_rank}, {full_rank}")
    print(f"  candidate rank                     = {rank_fraction}")
    check(rank_fraction not in s64_invariant_ranks(), "rank 63 breaks unlabelled S_64 slot symmetry")

    assumptions = [
        (
            "six labelled completion bits",
            "the 64 selector slots must be a physical completion register, not just a symmetric count",
        ),
        (
            "pre-latch acoustic coupling",
            "the Peierls/acoustic ruler must couple to every incomplete state and omit only |111111>",
        ),
        (
            "late selector completion",
            "late H0 must bill the completed 64/64 span while the acoustic ruler bills 63/64",
        ),
        (
            "no boundary-factor shortcut",
            "the missing slot must be a derived readout primitive, not an added '+1' or '-1' count",
        ),
    ]
    for name, requirement in assumptions:
        print(f"  needed primitive: {name:24s} {requirement}")

    print("\n[4] Current-canon anchors")
    for anchor in canon_anchors():
        ok = source_contains(anchor.path, anchor.required_patterns)
        rel = anchor.path.relative_to(ROOT)
        print(f"  {rel}: {anchor.meaning}")
        check(ok, f"{rel} still records the required guardrail")

    print("\n[5] Verdict")
    print(
        """
  Constructed:
    If the depth-6 selector is promoted to a labelled six-bit completion
    register, the acoustic pre-latch projector

        P_pre = I - |111111><111111|

    is a clean trace-63 operator.  It is idempotent, S6-invariant, and its
    normalized trace gives exactly 63/64.

  Rejected as a selector-lock repair under the current Canon:
    The existing service ledger has only symmetric 64 depth slots for this
    purpose.  A rank-63 projector is not S64-invariant.  The construction above
    therefore adds new physical structure: a labelled completion register plus
    an acoustic rule that makes the all-complete latch state invisible to the
    CMB ruler but visible to late H0.

  Live theorem target:
    Derive the acoustic readout map from a bridge/web or substrate service
    primitive.  The required statement is not "there are 64 slots and we count
    63"; it is:

        theta_* bills Tr(P_acoustic rho_selector) with
        P_acoustic = I - |complete><complete|,
        while the late selector bills Tr(I rho_selector).

    Until that map is derived, the 63/64 selector repair remains numerically
    compelling but canonically unlicensed.
exit 0 -- pre-latch projector constructed conditionally; selector-lock repair rejected as current theorem.
"""
    )


if __name__ == "__main__":
    main()
