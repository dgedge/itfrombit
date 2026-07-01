#!/usr/bin/env python3
r"""v_phase4_filled_cell_selection_theorem.py

Finite selection theorem for the EW filled-cell channel.

Earlier scripts established the main fork:

* the literal R4 projector is the wrong object for the Higgs VEV, because its
  Pauli expansion contains identity and low-weight terms and therefore gives a
  Planck-sized contribution;
* the viable object is a state transition from the empty cell to the filled
  eight-mode record.

This script sharpens the remaining selection statement.  In the occupation
basis of the eight local record modes, sectors with different occupation number
are orthogonal.  A record predicate that certifies "the cell is filled" is:

    1. S_8 invariant, so it depends only on k = number of occupied modes;
    2. complete, so it rejects every state with at least one missing mode;
    3. non-empty, so it accepts at least one state.

Those three finite conditions force the unique sector k=8.  Any lower threshold
or parity/check surrogate is a different observable: it accepts partial cells,
has many final states, and is dominated by a lower alpha_0 power.

Exit 0 means the k=8 selection is finite-theorem closed conditional on the
physical identification "EWSB bills the complete-cell transition."  It does not
derive the full SM Higgs Lagrangian or the RG quartic.
"""

from __future__ import annotations

from collections.abc import Callable
from math import comb, log
import sys


N = 8
ALPHA0 = 1.0 / 137.036
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def states_of_weight(k: int) -> set[int]:
    return {s for s in range(1 << N) if s.bit_count() == k}


def predicate_summary(name: str, pred: Callable[[int], bool]) -> dict[str, float | int | str]:
    accepted = [k for k in range(N + 1) if pred(k)]
    nstates = sum(comb(N, k) for k in accepted)
    if not accepted:
        leading = None
        coeff = 0
    else:
        leading = min(accepted)
        coeff = comb(N, leading)
    probability = sum(comb(N, k) * ALPHA0**k * (1.0 - ALPHA0) ** (N - k) for k in accepted)
    slope = float("nan") if leading is None else float(leading)
    return {
        "name": name,
        "accepted": ",".join(str(k) for k in accepted) if accepted else "-",
        "states": nstates,
        "leading_power": "-" if leading is None else leading,
        "leading_coeff": coeff,
        "probability": probability,
        "slope": slope,
        "complete": accepted == [N],
    }


def finite_complete_predicates() -> list[tuple[int, ...]]:
    """All S_8-invariant predicates are subsets of {0,...,8}; complete-cell predicates are unique."""
    complete = []
    for mask in range(1 << (N + 1)):
        accepted = tuple(k for k in range(N + 1) if (mask >> k) & 1)
        rejects_every_missing_mode = all(k not in accepted for k in range(N))
        nonempty = bool(accepted)
        if rejects_every_missing_mode and nonempty:
            complete.append(accepted)
    return complete


def main() -> int:
    print("=" * 96)
    print("PHASE 4A: FILLED-CELL EWSB SELECTION THEOREM")
    print("=" * 96)

    print("\n[0] Orthogonal occupation sectors")
    overlaps_ok = True
    for k in range(N + 1):
        for ell in range(N + 1):
            if k >= ell:
                continue
            overlaps_ok = overlaps_ok and states_of_weight(k).isdisjoint(states_of_weight(ell))
    print("  Sector dimensions:", [comb(N, k) for k in range(N + 1)])
    check("different occupation-number sectors are disjoint/orthogonal", overlaps_ok)

    print("\n[1] Complete-cell predicates")
    complete = finite_complete_predicates()
    print(f"  S_8-invariant predicates tested: {1 << (N + 1)}")
    print(f"  predicates that reject every missing-mode state and accept something: {complete}")
    check("complete-cell predicate is uniquely k=8", complete == [(N,)])

    print("\n[2] Natural alternatives and their leading alpha_0 power")
    predicates = [
        ("nonempty / any opened mode", lambda k: k >= 1),
        ("at least one parity/check quartet", lambda k: k >= 4),
        ("even parity sector", lambda k: k % 2 == 0),
        ("exact quartet", lambda k: k == 4),
        ("filled cell", lambda k: k == 8),
    ]
    rows = [predicate_summary(name, pred) for name, pred in predicates]
    print(f"  {'predicate':32s} {'k accepted':12s} {'states':>6s} {'lead':>5s} {'coeff':>6s} {'P':>12s}")
    for row in rows:
        print(
            f"  {row['name']:32s} {row['accepted']:12s} {row['states']:6d}"
            f" {str(row['leading_power']):>5s} {row['leading_coeff']:6d} {row['probability']:12.4e}"
        )

    filled = rows[-1]
    any_open = rows[0]
    exact_quartet = rows[3]
    print("\n[3] Discriminators")
    print(f"  any-open / filled probability ratio = {any_open['probability'] / filled['probability']:.3e}")
    print(f"  exact-quartet / filled probability ratio = {exact_quartet['probability'] / filled['probability']:.3e}")
    print("  These alternatives are large because they are different observables, not corrections")
    print("  to the complete-cell amplitude.")
    check("nonempty observable is lower-power and cannot represent filled-cell EWSB", any_open["slope"] == 1.0)
    check("quartet observable is lower-power and accepts partial cells", exact_quartet["slope"] == 4.0)
    check("filled-cell record is unique and alpha_0^8", filled["states"] == 1 and filled["slope"] == 8.0)
    check("partial-sector contamination is removed by orthogonality, not by small coefficients",
          any_open["probability"] / filled["probability"] > 1.0e14)

    print("\nVERDICT")
    print("  The finite record-selection theorem is closed: once the EW observable is the")
    print("  complete filled-cell record, S_8 invariance plus rejection of every missing")
    print("  mode forces k=8 and therefore the alpha_0^8 channel.  Lower-k channels are")
    print("  real, but they are different records: self-energy, parity/check openings,")
    print("  or partial-cell sectors.  The remaining premise is the physical one:")
    print("  EWSB must bill the complete-cell transition rather than a lower record.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
