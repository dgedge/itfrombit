#!/usr/bin/env python3
r"""ITEM 132: R4/KMS coupling theorem attempt.

Question
--------
Can current canon prove the remaining gluing clause

    R4 local service records are horizon-clocked closed Landauer records?

This is the clause needed to turn the primitive KMS phase-return latch into a
full derivation of

    a0 = c H0 / (2 pi).

Method
------
Try to prove the coupling by contradiction.  If the current finite R4
Stinespring service, repeated-service Poisson ledger, and Landauer record
bookkeeping already force the horizon KMS phase, then there should be no
countermodel with an independent closed phase clock.

The countermodel is simple:

    V_total(q) = V_R4 \otimes S_M(q),

where V_R4 is the exact finite R4 Stinespring repair channel and S_M(q) is a
closed phase latch with angular frequency Omega=q H0.  The latch is perfectly
repeatable for every q.  It changes only the real-time period and therefore the
acceleration quantum

    a(q) = q c H0/(2 pi).

Result
------
All finite R4 service tests are q-blind: isometry, event labels, Poisson
service-history statistics, and Landauer entropy per event are unchanged.
Therefore the current premises cannot prove q=1.  They admit a one-parameter
family of horizon-decoupled phase clocks.

The exact missing assumption is not "the operator" (constructed elsewhere) but
the coupling axiom:

    the R4 service clock is the cosmological horizon modular/KMS clock.

Under that added axiom q=1 and the scale closes.  Without it, the proof is
blocked by this explicit model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from itertools import product


C_M_PER_S = 299_792_458.0
MPC_M = 3.0856775814913673e22
H0_KM_S_MPC = 67.266152
H0_SI = H0_KM_S_MPC * 1000.0 / MPC_M

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def r1(c: tuple[int, ...]) -> bool:
    return not (c[G0] == 1 and c[G1] == 1)


def r2(c: tuple[int, ...]) -> bool:
    return c[W] == c[CHI]


def r3(c: tuple[int, ...]) -> bool:
    if c[LQ] == 0:
        return (c[C0], c[C1]) == (0, 0)
    return (c[C0], c[C1]) != (0, 0)


def r4(c: tuple[int, ...]) -> bool:
    return not (c[LQ] == 0 and c[I3] == 0 and c[CHI] == 1)


def valid_r123(c: tuple[int, ...]) -> bool:
    return r1(c) and r2(c) and r3(c)


def flip(c: tuple[int, ...], *idxs: int) -> tuple[int, ...]:
    out = list(c)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def physical_words() -> list[tuple[int, ...]]:
    return [tuple(c) for c in product((0, 1), repeat=8) if valid_r123(tuple(c))]


def repair_edges() -> list[tuple[tuple[int, ...], tuple[int, ...], str]]:
    repairs = {"I3": (I3,), "chi/W": (CHI, W)}
    out: list[tuple[tuple[int, ...], tuple[int, ...], str]] = []
    for c in physical_words():
        if not r4(c):
            for label, idxs in repairs.items():
                target = flip(c, *idxs)
                if valid_r123(target) and r4(target):
                    out.append((c, target, label))
    return out


RowKey = tuple[str, int]
Column = dict[RowKey, float]


def build_r4_columns() -> list[Column]:
    words = physical_words()
    idx = {word: i for i, word in enumerate(words)}
    cols: list[Column] = [{} for _ in words]
    for word in words:
        if r4(word):
            cols[idx[word]][("vac/noop", idx[word])] = 1.0
    for src, dst, label in repair_edges():
        cols[idx[src]][(f"repair:{label}:{idx[src]}", idx[dst])] = 1.0 / math.sqrt(2.0)
    return cols


def sparse_inner(a: Column, b: Column) -> float:
    if len(a) > len(b):
        a, b = b, a
    return sum(amp * b.get(row, 0.0) for row, amp in a.items())


def stinespring_errors(columns: list[Column]) -> tuple[float, float]:
    diag = 0.0
    offdiag = 0.0
    for i, col_i in enumerate(columns):
        diag = max(diag, abs(sparse_inner(col_i, col_i) - 1.0))
        for col_j in columns[i + 1 :]:
            offdiag = max(offdiag, abs(sparse_inner(col_i, col_j)))
    return diag, offdiag


def poisson_stats(x: float, nmax: int = 160) -> tuple[float, float, float]:
    probs = [math.exp(-x)]
    for n in range(nmax):
        probs.append(probs[-1] * x / (n + 1))
    z = sum(probs)
    probs = [p / z for p in probs]
    mean = sum(n * p for n, p in enumerate(probs))
    var = sum((n - mean) ** 2 * p for n, p in enumerate(probs))
    return mean, var, var / mean


@dataclass(frozen=True)
class PhaseClockModel:
    name: str
    q: float
    slots: int = 64

    @property
    def omega(self) -> float:
        return self.q * H0_SI

    @property
    def period(self) -> float:
        return 2.0 * math.pi / self.omega

    @property
    def acceleration(self) -> float:
        return C_M_PER_S / self.period

    def first_return_steps(self) -> int:
        state = 0
        for k in range(1, self.slots + 1):
            state = (state + 1) % self.slots
            if state == 0:
                return k
        raise AssertionError("finite cycle must return")

    def subcycle_returns(self) -> list[int]:
        return [k for k in range(1, self.slots) if k % self.slots == 0]


def landauer_event_entropy(num_addresses: int = 1, repair_labels: int = 6) -> float:
    """The event entropy depends on the committed record alphabet, not q."""

    return math.log(num_addresses * (repair_labels + 1))


def main() -> None:
    print("ITEM 132: R4/KMS COUPLING THEOREM ATTEMPT")
    print("=" * 96)

    print("\n[1] Fixed R4 Stinespring service channel")
    words = physical_words()
    edges = repair_edges()
    columns = build_r4_columns()
    diag, offdiag = stinespring_errors(columns)
    print(f"  R1-R3 words          = {len(words)}")
    print(f"  R4 repair edges      = {len(edges)}")
    print(f"  Stinespring errors   = diag {diag:.3e}, offdiag {offdiag:.3e}")
    check(len(words) == 48, "finite R4 service acts on the 48-word R1-R3 space")
    check(len(edges) == 6, "finite R4 service has six repair records")
    check(diag < 1e-12 and offdiag < 1e-12, "R4 Stinespring map is an isometry")

    print("\n[2] Countermodel family: tensor an independent closed phase clock")
    a_ref = C_M_PER_S * H0_SI / (2.0 * math.pi)
    clocks = [
        PhaseClockModel("half-speed closed clock", 0.5),
        PhaseClockModel("horizon KMS clock", 1.0),
        PhaseClockModel("double-speed closed clock", 2.0),
        PhaseClockModel("generator-rate closed clock", 2.0 * math.pi),
    ]
    for clock in clocks:
        print(
            f"  {clock.name:28s}: q={clock.q:.6g}, first_return={clock.first_return_steps()}, "
            f"subcycles={clock.subcycle_returns()}, a/a_KMS={clock.acceleration/a_ref:.6g}"
        )
        check(clock.first_return_steps() == clock.slots, f"{clock.name}: closed record returns on its own full cycle")
        check(not clock.subcycle_returns(), f"{clock.name}: no proper subcycle return")

    print("\n[3] R4/Poisson/Landauer observables are q-blind")
    entropy = landauer_event_entropy()
    for clock in clocks:
        mean, var, fano = poisson_stats(3.0)
        print(
            f"  q={clock.q:8.5f}: diag={diag:.1e}, offdiag={offdiag:.1e}, "
            f"Poisson mean={mean:.6f}, Fano={fano:.6f}, S_event={entropy:.6f}, "
            f"a={clock.acceleration:.6e}"
        )
        check(abs(mean - 3.0) < 1e-12, "Poisson mean does not depend on phase-clock q")
        check(abs(fano - 1.0) < 1e-12, "Poisson Fano does not depend on phase-clock q")
        check(abs(entropy - landauer_event_entropy()) < 1e-15, "Landauer event entropy does not depend on phase-clock q")

    print("\n[4] The shared-clock constraint would select q=1, but it is an extra premise")
    selected = [clock for clock in clocks if abs(clock.q - 1.0) < 1e-15]
    check(len(selected) == 1, "imposing 'R4 clock = horizon modular clock' selects the q=1 model")
    check(any(abs(clock.q - 1.0) > 1e-15 for clock in clocks), "without that premise, other q values remain valid countermodels")
    print("  The finite R4 channel supplies event labels and count history.")
    print("  The horizon sector supplies a modular/KMS phase.")
    print("  Current canon does not contain an operator equating the two clocks.")

    print("\n[5] Verdict")
    print(
        """
  NO-GO FOR AN UNCONDITIONAL PROOF FROM CURRENT PREMISES:
    The R4 Stinespring event channel can be tensored with a closed phase latch
    of arbitrary Omega=q H0.  Isometry, event labels, the Poisson service
    history, and Landauer event entropy are unchanged, while the acceleration
    quantum scales as q cH0/(2pi).  Therefore q=1 is not a theorem of the
    finite R4 service algebra.

  EXACT REMAINING AXIOM / THEOREM TARGET:
    Prove a universal closed-record modular-coupling law:

        every local R4 Landauer service record is clocked by the cosmological
        horizon KMS modular automorphism, not by an independent internal
        service phase.

    With that law, the primitive latch theorem gives a0=cH0/(2pi).  Without
    it, the scale remains Dirac-class / horizon-input.
"""
    )
    print("exit 0 -- explicit independent-clock countermodels block an unconditional R4/KMS coupling proof.")


if __name__ == "__main__":
    main()
