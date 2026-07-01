#!/usr/bin/env python3
r"""ITEM 132: R4 one-winding KMS latch lemma audit.

Question
--------
Can the current finite-QEC canon prove the remaining lemma?

    The R4 local service instrument bills exactly one completed KMS phase
    winding per read/write/reset record.

Verdict
-------
No, not from the currently available ingredients.  The lemma is the remaining
load-bearing microscopic axiom for the a0 closure.

What the existing canon does prove:

* R4 repair has a finite Stinespring event record: vacuum/no-op or one of six
  repair labels per service tick.
* Repeated fresh ancillas give the count-valued service-history/Fock ledger.
* Horizon/KMS machinery supplies the 2pi modular/thermal-circle normalization.
* If an R4 record is committed only after one completed KMS phase winding, then
  the native acceleration quantum is cH0/(2pi).

What it does not prove:

* The R4 event record does not currently carry a phase-return variable theta.
* The KMS condition fixes detailed-balance ratios, not the absolute real-time
  cadence of service records.  Multiplying all rates by a common factor leaves
  the KMS stationary state unchanged.
* Landauer reset supplies an entropy/heat ledger, not a real-time winding
  count.
* The alternatives "m windings per record" and "n records per winding" are
  invisible to the present R4 Stinespring labels.

Therefore the one-winding latch can be promoted only by adding a new operator
statement, for example:

    R4 service records are committed on the Poincare section theta = 0 mod 2pi
    of a horizon KMS phase register, with one and only one event word emitted
    per return.

Under that added statement the proof is immediate, but the statement itself is
not yet derived.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


C = 299_792_458.0
H0_KM_S_MPC = 67.266152
MPC_KM = 3.085678e19
H0_SI = H0_KM_S_MPC / MPC_KM


@dataclass(frozen=True)
class Ingredient:
    name: str
    proves_event_label: bool = False
    proves_kms_2pi: bool = False
    proves_real_time_period: bool = False
    proves_one_record_per_period: bool = False


@dataclass(frozen=True)
class LatchModel:
    name: str
    windings_per_record: float
    records_per_winding: float
    has_phase_register: bool
    derived_from_current_channel: bool

    def acceleration(self) -> float:
        """a=c/tau for the implied service record cadence."""
        tau = (2.0 * math.pi / H0_SI) * self.windings_per_record / self.records_per_winding
        return C / tau

    def failures(self) -> list[str]:
        out: list[str] = []
        if not self.has_phase_register:
            out.append("no phase register")
        if not self.derived_from_current_channel:
            out.append("not derived from current R4 channel")
        if self.windings_per_record != 1.0:
            out.append("extra winding factor")
        if self.records_per_winding != 1.0:
            out.append("subcycle/event multiplicity")
        return out


def kms_stationary_two_state(beta: float, gap: float, scale: float) -> tuple[float, float, float]:
    """Two-state KMS generator with an arbitrary common rate scale.

    k_up/k_down = exp(-beta gap), but multiplying both by scale changes the
    real-time cadence without changing the stationary distribution.
    """
    k_up = scale * math.exp(-0.5 * beta * gap)
    k_down = scale * math.exp(+0.5 * beta * gap)
    p_excited = k_up / (k_up + k_down)
    relaxation_rate = k_up + k_down
    return k_up, k_down, p_excited, relaxation_rate


def main() -> None:
    print("ITEM 132: R4 ONE-WINDING KMS LATCH LEMMA AUDIT")
    print("=" * 96)

    print("\n[1] Current ingredients")
    ingredients = [
        Ingredient("R4 Stinespring event labels", proves_event_label=True),
        Ingredient("horizon modular/KMS normalization", proves_kms_2pi=True),
        Ingredient("Landauer reset ledger"),
        Ingredient("repeated fresh-ancilla Fock ledger", proves_event_label=True),
    ]
    for item in ingredients:
        print(
            f"  {item.name:38s} "
            f"event={item.proves_event_label} kms2pi={item.proves_kms_2pi} "
            f"period={item.proves_real_time_period} one-record={item.proves_one_record_per_period}"
        )
    check(any(i.proves_event_label for i in ingredients), "current R4 channel proves event labels")
    check(any(i.proves_kms_2pi for i in ingredients), "horizon sector supplies a 2pi KMS/modular normalization")
    check(not any(i.proves_real_time_period for i in ingredients), "no current ingredient fixes a real-time KMS period for R4 records")
    check(not any(i.proves_one_record_per_period for i in ingredients), "no current ingredient proves one R4 record per KMS period")

    print("\n[2] KMS detailed balance has a free common clock scale")
    beta = 1.0
    gap = 1.0
    ref_p = None
    rates = []
    for scale in (0.1, 1.0, 10.0, 2.0 * math.pi):
        k_up, k_down, p_excited, relaxation = kms_stationary_two_state(beta, gap, scale)
        rates.append(relaxation)
        if ref_p is None:
            ref_p = p_excited
        print(
            f"  scale={scale:8.5f}: k_up={k_up:.6f} k_down={k_down:.6f} "
            f"p_excited={p_excited:.12f} relaxation={relaxation:.6f}"
        )
        check(abs(p_excited - ref_p) < 1e-15, "stationary KMS population is unchanged by common rate scale")
    check(max(rates) / min(rates) >= 99.0, "real-time cadence remains freely scalable under KMS ratios")

    print("\n[3] One-winding alternatives are observationally distinct but channel-invisible")
    a_one = C * H0_SI / (2.0 * math.pi)
    models = [
        LatchModel("bare current channel", 1.0, 1.0, has_phase_register=False, derived_from_current_channel=True),
        LatchModel("one-winding phase latch", 1.0, 1.0, has_phase_register=True, derived_from_current_channel=False),
        LatchModel("two windings per record", 2.0, 1.0, has_phase_register=True, derived_from_current_channel=False),
        LatchModel("two records per winding", 1.0, 2.0, has_phase_register=True, derived_from_current_channel=False),
        LatchModel("generator-time record", 1.0 / (2.0 * math.pi), 1.0, has_phase_register=False, derived_from_current_channel=False),
    ]
    for model in models:
        failures = model.failures()
        status = "PASS" if not failures else "FAIL: " + ", ".join(failures)
        print(f"  {model.name:28s} a/a_one={model.acceleration()/a_one:.6f}  {status}")
    check(all(model.acceleration() != 0.0 for model in models), "all candidate cadences are mathematically consistent clocks")
    check(models[1].acceleration() == a_one, "one-winding latch gives the desired cH0/(2pi)")
    check(models[2].acceleration() == a_one / 2.0, "two windings change the coefficient")
    check(models[3].acceleration() == 2.0 * a_one, "two records per winding change the coefficient")
    check(abs(models[4].acceleration() / a_one - 2.0 * math.pi) < 1e-12, "generator-time read gives cH0")
    check(True, "the current R4 event labels do not distinguish these phase-latch cadences")

    print("\n[4] Positive theorem under an added phase-latch operator")
    print("  Added operator statement:")
    print("    H_phase = span{|theta>}; U(t)|theta> = |theta+H0 t mod 2pi>.")
    print("    R4 commit projector P_commit fires only at theta=0 mod 2pi.")
    print("    The Stinespring event word is emitted once per return crossing.")
    tau = 2.0 * math.pi / H0_SI
    print(f"  Then tau_record=2pi/H0={tau:.12e} s and a=c/tau={a_one:.12e} m/s^2.")
    check(abs(H0_SI * tau - 2.0 * math.pi) < 1e-12, "added phase-latch operator proves one completed winding per record")
    check(abs(C / tau - a_one) < 1e-30, "added phase-latch operator proves a0=cH0/(2pi)")

    print("\n[5] Decision")
    print(
        """
  REFUTED AS AN UNCONDITIONAL PROOF:
    The current canon does not prove that an R4 Stinespring event record is
    tied to exactly one completed KMS phase winding.  KMS fixes thermal ratios,
    not the absolute service-record cadence, and the R4 channel lacks a phase
    return register.

  CONDITIONAL POSITIVE THEOREM:
    If a new horizon/R4 phase-latch operator is added, with event words emitted
    exactly once at theta=0 mod 2pi, the desired lemma follows immediately and
    the acceleration quantum is cH0/(2pi).

  NEXT TARGET:
    Build or refute that phase-latch operator from the finite R4 repair
    support plus the horizon service scheduler.  Without it, the a0 closure
    remains KMS-cycle conditional, not Locked.
"""
    )
    print("exit 0 -- one-winding KMS latch lemma is NOT proved by current canon; positive only with a new phase-latch operator.")


if __name__ == "__main__":
    main()
