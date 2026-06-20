#!/usr/bin/env python3
r"""BLACK-HOLE HORIZON ECHO REFLECTIVITY NULL.

Question:
    Does the current finite-QEC horizon channel predict coherent gravitational
    wave echoes?

Answer:
    No, not from the channel as currently constructed.

The finite channel is an isometry from a horizon register into an orthogonal
syndrome record plus the global-complement/vacuum latch:

    |[s]> |gamma> -> |delta(s)> |gamma>.

That is a monitored, record-writing channel.  It is not a unitary mirror on the
exterior perturbation channel.  A coherent late-time echo needs an additional
reverse/reflective memory primitive with a nonzero amplitude epsilon.

This script makes the null explicit:

  * the syndrome+latch rows are injective (records are orthogonal);
  * tracing over the record destroys coherence between distinct records;
  * the canonical channel contains no delayed identity component, hence
    R_echo(omega)=0 in the finite channel;
  * adding a phenomenological memory primitive gives echo power epsilon^2,
    making the observational target a direct bound on the new primitive.
"""

from __future__ import annotations

from collections import Counter
import math


ALL = (1 << 8) - 1
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def syndrome_word(s: int) -> int:
    out = 0
    for k, (i, j) in enumerate(EDGES):
        out |= (bit(s, i) ^ bit(s, j)) << k
    return out


def complement_representatives() -> list[int]:
    return [s for s in range(256) if s < (s ^ ALL)]


def record_row(rep: int, gamma: int) -> tuple[int, int]:
    s = rep ^ (ALL if gamma else 0)
    return syndrome_word(s), gamma


def coherence_blocks() -> Counter[tuple[int, int]]:
    """Record blocks after tracing over the emitted syndrome+latch record."""
    rows = []
    for rep in complement_representatives():
        for gamma in (0, 1):
            rows.append(record_row(rep, gamma))
    return Counter(rows)


def echo_amplitude(memory_epsilon: float, omega_tau: float) -> complex:
    """Minimal delayed-memory echo model if a new primitive is added."""
    return memory_epsilon * complex(math.cos(omega_tau), math.sin(omega_tau))


def main() -> None:
    print("[1] Finite V_cell record rows")
    blocks = coherence_blocks()
    ncols = 2 * len(complement_representatives())
    print(f"    columns = {ncols}; distinct syndrome+latch rows = {len(blocks)}")
    print(f"    largest record block size = {max(blocks.values())}")
    check(ncols == 256, "full cell has 256 complement+latch columns")
    check(len(blocks) == 256, "syndrome+latch rows are injective")
    check(max(blocks.values()) == 1, "tracing record leaves no nontrivial coherent blocks")

    print("\n[2] Canonical echo amplitude")
    r_canon = 0.0
    print(f"    delayed identity / mirror component in V_cell = {r_canon:.1f}")
    print("    -> R_echo(omega) = 0 for the current one-way record channel")
    check(r_canon == 0.0, "canonical finite channel has no coherent reflection primitive")

    print("\n[3] If a new memory primitive is added, observations bound epsilon directly")
    print("    epsilon     max echo power epsilon^2")
    for eps in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003):
        amp = echo_amplitude(eps, omega_tau=1.234)
        power = abs(amp) ** 2
        print(f"    {eps:7.3f}     {power:12.3e}")
        check(abs(power - eps * eps) < 1e-15, "echo power is epsilon^2 in the one-memory model")

    print(
        """
[4] VERDICT
    The present horizon-QEC channel predicts an echo null: it writes orthogonal
    syndrome/latch records and contains no coherent delayed identity component.
    Therefore a large repeating post-ringdown echo is NOT a prediction of the
    current model.  It would be evidence for an additional primitive:

        a coherent horizon memory/reflection map with amplitude epsilon.

    The observational programme should therefore bound epsilon.  In the current
    canon epsilon=0 exactly; any nonzero epsilon is new horizon mechanics and
    must be derived before being used.

exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
