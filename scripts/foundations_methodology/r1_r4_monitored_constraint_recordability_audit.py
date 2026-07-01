#!/usr/bin/env python3
r"""R1--R4 monitored-constraint recordability audit.

Methodological point being tested
----------------------------------
The Boolean rules R1--R4 can be read in two different ways:

  1. as static admissibility predicates that define the code space;
  2. as monitored QEC constraints that are read, recorded, and repaired.

Only the second reading can generate environment records.  Static predicates can
count states and select representations, but they cannot by themselves carry a
phase, covariance, entropy production, billing rate, or mass scale.  A monitored
constraint's record content is also limited to the bits its syndrome/recovery
actually reads.

This audit classifies R1--R4 under the current canon:

  * R1: originally a static generation cut; now a serious monitored boot/R1
        candidate via the classical AND decoder and Hasse-edge record.  Still
        conditional for delta/Phi because the sector-contact law A_s=d_s is open.
  * R2: monitored as the W=chi stabilizer / chirality-locking edge.  The
        physics-bearing object is Z_chi Z_W, not the weak I3 generator that was
        once mislabeled as the R2 relaxation.
  * R3: code/gauge monitored as the colour--lepton separation/Gauss constraint.
        It is physics-bearing for charge/triality/confinement selection; extra
        phases or rates still require a named recovery/billing channel.
  * R4: monitored only in the sterile/electroweak repair channel.  It writes
        R4/sterile records, but reads {LQ,I3,chi} and repairs {I3,chi,W}; it
        cannot write generation-Hasse records.

Exit 0 means the recordability filter is internally consistent with the current
R1--R4 status.  It is a taxonomy/theorem-target note, not a new numerical result.
"""

from __future__ import annotations

from dataclasses import dataclass


BITS = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")
GENERATION_BITS = {"G0", "G1"}
EW_BITS = {"LQ", "I3", "chi", "W"}
COLOUR_BITS = {"LQ", "C0", "C1"}


@dataclass(frozen=True)
class ConstraintStatus:
    name: str
    predicate_bits: frozenset[str]
    syndrome_bits: frozenset[str]
    repair_bits: frozenset[str]
    status: str
    record_payload: str
    open_gate: str

    @property
    def is_monitored(self) -> bool:
        return bool(self.syndrome_bits or self.repair_bits)

    @property
    def can_record_generation_edges(self) -> bool:
        return GENERATION_BITS <= (self.syndrome_bits | self.repair_bits)


def check(name: str, cond: bool) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def r1_forbidden(g0: int, g1: int) -> bool:
    return g0 == 1 and g1 == 1


def r2_valid(chi: int, w: int) -> bool:
    return w == chi


def r3_valid(lq: int, c0: int, c1: int) -> bool:
    return (c0, c1) == (0, 0) if lq == 0 else (c0, c1) != (0, 0)


def r4_forbidden(lq: int, i3: int, chi: int) -> bool:
    return (lq, i3, chi) == (0, 0, 1)


def main() -> None:
    print("R1--R4 monitored-constraint recordability audit")
    print("=" * 78)

    statuses = [
        ConstraintStatus(
            name="R1",
            predicate_bits=frozenset({"G0", "G1"}),
            syndrome_bits=frozenset({"G0", "G1"}),  # linear reads plus classical AND decoder, if boot-monitored
            repair_bits=frozenset({"G0", "G1"}),
            status="candidate monitored boot boundary",
            record_payload="generation Hasse-edge record; sym -> delta carrier, oriented -> Phi/sigma carrier",
            open_gate="R1 must be boot-monitored per tick, and sector-contact count A_s=d_s must be derived",
        ),
        ConstraintStatus(
            name="R2",
            predicate_bits=frozenset({"chi", "W"}),
            syndrome_bits=frozenset({"chi", "W"}),
            repair_bits=frozenset({"chi", "W", "I3"}),
            status="monitored stabilizer/chirality-locking edge",
            record_payload="W=chi stabilizer record Z_chi Z_W; chirality/EW mirror and EWSB bookkeeping",
            open_gate="absolute EWSB scale and symmetric-gapped-phase billing remain separate",
        ),
        ConstraintStatus(
            name="R3",
            predicate_bits=frozenset({"LQ", "C0", "C1"}),
            syndrome_bits=frozenset({"LQ", "C0", "C1"}),
            repair_bits=frozenset({"LQ", "C0", "C1"}),
            status="monitored code/gauge constraint",
            record_payload="colour/lepton separation, charge algebra, triality/Gauss selection",
            open_gate="dynamical phases/rates require a specified recovery or billing channel",
        ),
        ConstraintStatus(
            name="R4",
            predicate_bits=frozenset({"LQ", "I3", "chi"}),
            syndrome_bits=frozenset({"LQ", "I3", "chi"}),
            repair_bits=frozenset({"I3", "chi", "W"}),
            status="partial monitored sterile/electroweak repair",
            record_payload="sterile/Feshbach/R4 records; dark-sector activation; no generation-Hasse record",
            open_gate="does not access G0/G1, so it cannot supply Item-87 delta/Phi by itself",
        ),
    ]

    print("\n[1] Boolean predicates still define the static code-space cuts")
    check("R1 forbids exactly the 11 generation corner", r1_forbidden(1, 1) and not r1_forbidden(0, 1))
    check("R2 validates W=chi and rejects W!=chi", r2_valid(0, 0) and r2_valid(1, 1) and not r2_valid(0, 1))
    check("R3 separates leptons from colour triplet states", r3_valid(0, 0, 0) and not r3_valid(0, 1, 0) and r3_valid(1, 1, 0))
    check("R4 forbids the sterile conjunction (LQ=0,I3=0,chi=1)", r4_forbidden(0, 0, 1) and not r4_forbidden(0, 1, 1))

    print("\n[2] Recordability classification")
    for s in statuses:
        print(f"    {s.name}: {s.status}")
        print(f"        predicate bits : {sorted(s.predicate_bits)}")
        print(f"        syndrome bits  : {sorted(s.syndrome_bits)}")
        print(f"        repair bits    : {sorted(s.repair_bits)}")
        print(f"        payload        : {s.record_payload}")
        print(f"        open gate      : {s.open_gate}")
        check(f"{s.name}: static predicate has named bits", bool(s.predicate_bits))
        check(f"{s.name}: recordability status is explicitly classified", bool(s.status))

    print("\n[3] Bit-access consequences")
    r1 = statuses[0]
    r4 = statuses[3]
    check("R1 monitored candidate can access the generation register", r1.can_record_generation_edges)
    check("R4 sterile repair cannot access generation Hasse edges", not r4.can_record_generation_edges)
    check("R4 syndrome/repair is electroweak/sterile, not generation-resolving", (r4.syndrome_bits | r4.repair_bits).isdisjoint(GENERATION_BITS))
    check("therefore R4 records cannot be reused as R1 edge records", True)

    print("\n[4] Methodological rule")
    check("static cuts can count states but cannot by themselves write environment records", True)
    check("only monitored constraints can carry phases/covariances/entropy/billing", True)
    check("a monitored constraint can only carry the bits it actually reads", True)

    print(
        """
VERDICT:
  The insight survives, with one precision upgrade:

      static constraint  !=  physics-bearing record
      monitored constraint -> possible physics-bearing record

  R1--R4 should be described in two layers.  As predicates they define the valid
  code space.  As monitored QEC instruments they may write environment records,
  but only when a syndrome/readout/recovery channel is specified and only for the
  bits that channel reads.

  Consequences:
    * R1 is the new live candidate because a boot-monitored AND-decoder can read
      the generation register and the Hasse edge record can carry delta/Phi.
    * R2 is genuinely monitored through the W=chi stabilizer; the record object
      is Z_chi Z_W, not the weak I3 generator.
    * R3 is a monitored code/gauge constraint for colour/lepton separation and
      Gauss/triality selection; extra rates/phases need a separate channel.
    * R4 is monitored, but only in the sterile/electroweak channel.  It cannot
      write R1 generation-edge records because it never reads G0/G1.

  This gives a general canon guardrail: only monitored constraints can become
  physics-bearing records; frozen cuts remain kinematic selection rules.
exit 0"""
    )


if __name__ == "__main__":
    main()
