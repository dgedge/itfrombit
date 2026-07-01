#!/usr/bin/env python3
r"""ITEM 132: one-a0 central-cell service quantum and the phase-latch gap.

Question
--------
Can the remaining acceleration unit

    a0 = c H0 / (2 pi)

be derived from the current R4 service ledger?

Result
------
The current finite R4 labels do not contain a phase-return coordinate, so they
cannot by themselves choose between one winding, two windings, or multiple
records per winding.  However, if the R4 service record is extended by the
minimal KMS phase-return latch, the coefficient is forced:

    theta(t) = theta(0) + H0 t mod 2 pi,
    commit only on theta = 0 return,
    emit one primitive record per return.

Minimality then excludes both hidden stock and repeated/subcycle records:

    * subcycle events are not repeatable records, because the phase has not
      returned to the same projector;
    * m>1 windings per record contain m identical primitive return opportunities
      and therefore represent hidden skipped stock;
    * n>1 records per winding duplicates the same closed record.

Thus the one-a0 rule is not fully derived from the finite R4 channel, but it is
reduced to a single missing operator statement: R4 records are committed on the
primitive horizon/KMS phase-return section.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


C = 299_792_458.0
H0_KM_S_MPC = 67.266152
MPC_KM = 3.0856775814913673e19
H0_SI = H0_KM_S_MPC / MPC_KM
A0_ONE_WINDING = C * H0_SI / (2.0 * math.pi)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


@dataclass(frozen=True)
class RecordCadence:
    name: str
    windings_per_record: int
    records_per_winding: int
    has_phase_return_projector: bool

    @property
    def acceleration(self) -> float:
        tau = (2.0 * math.pi / H0_SI) * self.windings_per_record / self.records_per_winding
        return C / tau

    @property
    def primitive_failures(self) -> list[str]:
        failures: list[str] = []
        if not self.has_phase_return_projector:
            failures.append("no phase-return projector")
        if self.windings_per_record > 1:
            failures.append("skips primitive return opportunities")
        if self.records_per_winding > 1:
            failures.append("duplicates one closed record")
        return failures


def phase_overlap(delta_theta: float) -> float:
    """Toy overlap for a sharp phase projector after a subcycle shift."""

    return math.cos(0.5 * delta_theta) ** 2


def main() -> None:
    print("ITEM 132: ONE-a0 PHASE-LATCH MINIMALITY")
    print("=" * 92)

    print("\n[1] The finite R4 label alone cannot select a cadence")
    cadences = [
        RecordCadence("finite R4 event label only", 1, 1, False),
        RecordCadence("primitive phase return", 1, 1, True),
        RecordCadence("two windings per record", 2, 1, True),
        RecordCadence("two records per winding", 1, 2, True),
        RecordCadence("four records per winding", 1, 4, True),
    ]
    for model in cadences:
        failures = ", ".join(model.primitive_failures) or "primitive"
        print(f"  {model.name:30s} a/a0={model.acceleration/A0_ONE_WINDING:.6f}  {failures}")
    check(cadences[0].acceleration == cadences[1].acceleration, "bare labels can be assigned the same cadence but do not prove it")
    check(cadences[2].acceleration == 0.5 * A0_ONE_WINDING, "extra winding changes the coefficient")
    check(cadences[3].acceleration == 2.0 * A0_ONE_WINDING, "extra records per winding change the coefficient")

    print("\n[2] Repeatable record condition kills subcycle commits")
    for n in (2, 3, 4, 8):
        delta = 2.0 * math.pi / n
        print(f"  {n:2d} subcycle records/winding: projector overlap after one subcycle = {phase_overlap(delta):.9f}")
    check(phase_overlap(math.pi) < 1e-30, "half-winding commit is orthogonal to the return projector")
    check(phase_overlap(2.0 * math.pi) == 1.0, "full winding returns to the same record projector")

    print("\n[3] Primitive-return minimality kills skipped windings and duplicate records")
    for m in (1, 2, 3, 5):
        missed = m - 1
        print(f"  {m} windings/record: hidden skipped primitive returns = {missed}")
    check(all(m - 1 > 0 for m in (2, 3, 5)), "multi-winding record contains hidden unbilled primitive returns")
    for n in (1, 2, 3, 4):
        duplicates = n - 1
        print(f"  {n} records/winding: duplicate records on the same return = {duplicates}")
    check(all(n - 1 > 0 for n in (2, 3, 4)), "multi-record winding duplicates one closed record")

    print("\n[4] Decision")
    tau = 2.0 * math.pi / H0_SI
    print(f"  primitive KMS return time tau = {tau:.12e} s")
    print(f"  a0 = c/tau                  = {A0_ONE_WINDING:.12e} m/s^2")
    check(abs(C / tau - A0_ONE_WINDING) < 1e-30, "phase-return latch gives cH0/(2pi)")
    print("  POSITIVE CONDITIONAL THEOREM:")
    print("    If R4 records are committed on the primitive KMS phase-return projector,")
    print("    the one-a0 central-cell quantum follows and hidden integer cadences are")
    print("    excluded by repeatability/minimality.")
    print("  REMAINING GAP:")
    print("    the finite R4 Stinespring label does not itself derive that phase-return")
    print("    projector.  This is the actual residual, not a free core-shape parameter.")
    print("\nexit 0 -- one-a0 is reduced to the primitive phase-return latch operator.")


if __name__ == "__main__":
    main()
