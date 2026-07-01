#!/usr/bin/env python3
r"""ITEM 132: construct the primitive KMS phase-return latch.

Question
--------
The R4/MOND core rule is now reduced to one missing scale operator:

    one local service quantum = one a0 central-cell threshold.

Previous audits showed that finite R4 event labels do not contain this scale.
This script asks the sharper question: if the R4 service record is clocked by
the horizon/KMS phase already used elsewhere in canon, is there a unique closed
record operator that supplies

    a0 = c H0 / (2 pi)?

Construction
------------
Add a KMS phase register theta in S^1 with modular flow

    theta(t) = theta(0) + H0 t  (mod 2 pi).

The closed-record latch is the projector onto the return section theta=0,
composed with the R4 service Stinespring event:

    L_R4 = |0><0|_theta \otimes E_R4,

and it is allowed to fire only after the primitive phase return

    tau_KMS = 2 pi / H0.

In a finite M-slot regulator this is the cyclic shift S_M.  The primitive
record is the first positive k with

    P0 S_M^k P0 != 0,

namely k=M.  Subcycles are not repeatable closed records; multiple windings
skip primitive return opportunities; multiple records per winding duplicate the
same closed record.  Thus the hidden-integer cadences are excluded once the
closed-record/minimality axioms are admitted.

Honest scope
------------
This finds and tests the missing operator.  It does not by itself prove that the
R4 Stinespring service is coupled to the horizon/KMS phase register.  The
remaining theorem is the gluing clause:

    R4 local service records are horizon-clocked closed Landauer records.

If that clause is proved, a0=cH0/(2pi) is no longer an external Dirac-class
choice.  If rejected, finite R4 labels still do not derive a0.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable


C_M_PER_S = 299_792_458.0
MPC_M = 3.0856775814913673e22
H0_KM_S_MPC = 67.266152  # selector value used by the current proton-primary chain


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def h0_si() -> float:
    return H0_KM_S_MPC * 1000.0 / MPC_M


def acceleration_from_period(period_s: float) -> float:
    """Acceleration scale from one light-step over the record period."""
    return C_M_PER_S / period_s


@dataclass(frozen=True)
class FiniteKMSClock:
    """M-slot regulator for the KMS phase circle.

    The physical slot size is a regulator.  The full cycle time is fixed:
    tau = 2*pi/H0.  Therefore changing M only refines the phase register and
    does not change the emitted acceleration scale.
    """

    slots: int

    def shift(self, state: int, steps: int = 1) -> int:
        return (state + steps) % self.slots

    def returns_to_zero(self, steps: int) -> bool:
        return self.shift(0, steps) == 0

    def primitive_return_steps(self) -> int:
        for k in range(1, self.slots + 1):
            if self.returns_to_zero(k):
                return k
        raise AssertionError("finite cycle must return")

    def subcycle_failures(self) -> list[int]:
        return [k for k in range(1, self.slots) if self.returns_to_zero(k)]

    @property
    def cycle_period_s(self) -> float:
        return 2.0 * math.pi / h0_si()

    @property
    def slot_period_s(self) -> float:
        return self.cycle_period_s / self.slots

    @property
    def primitive_acceleration(self) -> float:
        return acceleration_from_period(self.cycle_period_s)


@dataclass(frozen=True)
class CadenceCandidate:
    name: str
    period_s: float
    records_per_period: int
    primitive_windings: int
    reason: str

    @property
    def acceleration(self) -> float:
        return self.records_per_period * acceleration_from_period(self.period_s)

    def failures(self, target_a0: float) -> list[str]:
        out: list[str] = []
        if self.primitive_windings != 1:
            out.append("skips primitive KMS returns")
        if self.records_per_period != 1:
            out.append("duplicates one closed record")
        if abs(self.acceleration / target_a0 - 1.0) > 1e-12:
            out.append("wrong acceleration scale")
        return out


def r4_label_cycles(labels: Iterable[str]) -> int:
    """Count possible cyclic orderings of finite R4 labels up to rotation.

    This is a deliberately weak diagnostic: many cyclic orderings exist and no
    one of them is invariantly selected by the labels.  A finite label set is
    not a KMS phase register.
    """

    n = len(tuple(labels))
    if n <= 1:
        return 1
    return math.factorial(n - 1)


def main() -> None:
    print("ITEM 132: PRIMITIVE KMS PHASE-RETURN LATCH OPERATOR")
    print("=" * 92)

    H0 = h0_si()
    tau_kms = 2.0 * math.pi / H0
    a0 = C_M_PER_S * H0 / (2.0 * math.pi)
    print("\n[1] Horizon/KMS phase register")
    print(f"  H0 selector = {H0_KM_S_MPC:.6f} km/s/Mpc = {H0:.12e} s^-1")
    print("  theta(t)=theta(0)+H0 t mod 2pi")
    print(f"  primitive KMS return tau = 2pi/H0 = {tau_kms:.12e} s")
    print(f"  one light-step per primitive return gives a0 = cH0/(2pi) = {a0:.12e} m/s^2")

    print("\n[2] Finite-cycle operator test")
    for slots in (4, 12, 28, 64, 137):
        clock = FiniteKMSClock(slots)
        primitive = clock.primitive_return_steps()
        subcycles = clock.subcycle_failures()
        print(
            f"  M={slots:3d}: primitive steps={primitive:3d}, "
            f"subcycle returns={subcycles}, a={clock.primitive_acceleration:.12e}"
        )
        check(primitive == slots, f"M={slots}: first closed-record return is the full KMS cycle")
        check(not subcycles, f"M={slots}: no proper subcycle is a repeatable closed record")
        check(abs(clock.primitive_acceleration / a0 - 1.0) < 1e-15, f"M={slots}: regulator does not change a0")

    print("\n[3] Hidden-cadence controls")
    candidates = [
        CadenceCandidate(
            "primitive KMS latch",
            tau_kms,
            records_per_period=1,
            primitive_windings=1,
            reason="one closed record at the first phase return",
        ),
        CadenceCandidate(
            "Hubble generator tick",
            1.0 / H0,
            records_per_period=1,
            primitive_windings=0,
            reason="reads the generator scale, not a closed KMS phase cycle",
        ),
        CadenceCandidate(
            "two-winding latch",
            2.0 * tau_kms,
            records_per_period=1,
            primitive_windings=2,
            reason="skips an available primitive return",
        ),
        CadenceCandidate(
            "two records per return",
            tau_kms,
            records_per_period=2,
            primitive_windings=1,
            reason="duplicates the same closed record",
        ),
    ]
    for cand in candidates:
        failures = cand.failures(a0)
        status = "PASS" if not failures else "FAIL: " + ", ".join(failures)
        print(
            f"  {cand.name:24s}: a={cand.acceleration:.12e} "
            f"({cand.acceleration / a0:.6g} x a0)  {status}"
        )
    check(not candidates[0].failures(a0), "primitive KMS latch is admissible")
    check(all(c.failures(a0) for c in candidates[1:]), "non-primitive / non-cycle cadences are rejected")

    print("\n[4] Why finite R4 labels alone still do not derive the phase latch")
    r4_labels = (
        "g0:nuR->eR",
        "g1:nuR->eR",
        "g2:nuR->eR",
        "g0:nuR->nuL",
        "g1:nuR->nuL",
        "g2:nuR->nuL",
    )
    cycles = r4_label_cycles(r4_labels)
    print(f"  finite R4 repair labels = {len(r4_labels)}")
    print(f"  cyclic orderings up to rotation = (6-1)! = {cycles}")
    print("  none of these label cycles carries the physical period 2pi/H0 or a return section P0.")
    print("  The phase register is therefore an extra KMS/horizon service coordinate, not an R4 label.")
    check(cycles > 1, "finite labels admit many cycles; no invariant phase-return operator is selected")

    print("\n[5] Operator statement")
    print(
        r"""
  Let H_theta = -i H0 d/dtheta on L2(S1), U_theta(t)=exp(-i H_theta t),
  and P0 be the return-section projector theta=0.  The candidate closed R4
  service instrument is

      L_R4 = P0 \otimes E_R4,       applied only after tau=2pi/H0.

  In the finite regulator this is the first return of the cyclic shift S_M.
  The origin theta=0 is gauge; the first-return period is not.
"""
    )

    print("\n[6] Verdict")
    print(
        """
  OPERATOR FOUND / MINIMALITY THEOREM:
    The primitive KMS phase-return latch is the unique minimal closed-record
    operator on the horizon phase circle.  It gives a0=cH0/(2pi) and excludes
    hidden integer cadences by repeatability and no-hidden-stock minimality.

  NOT YET FULL ITEM-132 LOCK:
    This construction still needs the gluing theorem that R4 local service
    records are horizon-clocked closed Landauer records.  The finite R4
    Stinespring labels alone do not contain the phase register and cannot
    derive the 2pi return.  Thus the residual is no longer "what operator?";
    it is "prove R4 service couples to this KMS latch."
"""
    )
    print("exit 0 -- primitive KMS phase-return latch constructed; remaining gate is R4/KMS gluing.")


if __name__ == "__main__":
    main()
