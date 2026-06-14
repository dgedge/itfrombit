#!/usr/bin/env python3
"""
Open-link center Gauss orbits for the matter-coupled CSS/TCH operator.

Motivation:
  css_matter_coupled_gauss_proxy.py found a clean Z2 story for U(1)/SU(2), but
  a negative SU(3) diagnostic: Pauli CSS X-sector signs are Z2 data, while the
  SU(3) center is Z3. This script isolates that point before any Wilson
  character truncation.

Setup:
  Four oriented links around one plaquette carry center exponents

      a_l in Z_N,   U_l = exp(2 pi i a_l / N).

  Use the orientation compatible with U0 U1 U2^dag U3^dag:

      e0: v0 -> v1
      e1: v1 -> v2
      e2: v3 -> v2
      e3: v0 -> v3

  A vertex Gauss transformation h_v in Z_N acts by

      a_l -> a_l + h_tail(l) - h_head(l).

  The plaquette center flux

      F = a0 + a1 - a2 - a3  mod N

  is invariant. The physical one-plaquette sectors are Gauss orbits labeled by F.

What this says about candidate TCH operators:
  * For N=2, this reproduces the Pauli CSS sign story: two orbits, exactly the
    P_X = prod_l X_l = +/-1 closed-flux sectors. The candidate projector is the
    familiar P_B=(1+B)/2.

  * For N=3, the analogous operator is not Pauli. It is a Z3 clock/Bianchi
    constraint P_B=(1+B+B^dag)/3 with B carrying the oriented center flux.
    A binary Pauli subset is not closed under the Z3 Gauss action.

This is still center-only. It does not build the full SU(3) Wilson/colour
Hilbert space, but it identifies a necessary algebraic feature of the missing
matter-coupled TCH operator.

pure Python; self-asserting.
"""

import itertools


EDGES = [
    (0, 1),  # U0
    (1, 2),  # U1
    (3, 2),  # U2, so the plaquette uses U2^dag
    (0, 3),  # U3, so the plaquette uses U3^dag
]


def hd(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def all_link_states(order):
    return list(itertools.product(range(order), repeat=4))


def flux(state, order):
    a0, a1, a2, a3 = state
    return (a0 + a1 - a2 - a3) % order


def transform(state, vertex_shift, order):
    out = []
    for value, (tail, head) in zip(state, EDGES):
        out.append((value + vertex_shift[tail] - vertex_shift[head]) % order)
    return tuple(out)


def orbit(seed, order):
    seen = set()
    stack = [seed]
    generators = []
    for vertex in range(4):
        shift = [0, 0, 0, 0]
        shift[vertex] = 1
        generators.append(tuple(shift))

    while stack:
        state = stack.pop()
        if state in seen:
            continue
        seen.add(state)
        for generator in generators:
            target = transform(state, generator, order)
            if target not in seen:
                stack.append(target)
    return seen


def orbit_decomposition(order):
    remaining = set(all_link_states(order))
    orbits = []
    while remaining:
        seed = next(iter(remaining))
        current = orbit(seed, order)
        orbits.append(current)
        remaining -= current
    return orbits


def binary_subset_closed(order, binary_value):
    subset = set(itertools.product([0, binary_value], repeat=4))
    leaked_examples = []
    for state in subset:
        for vertex in range(4):
            shift = [0, 0, 0, 0]
            shift[vertex] = 1
            target = transform(state, tuple(shift), order)
            if target not in subset:
                leaked_examples.append((state, vertex, target))
                if len(leaked_examples) >= 3:
                    return False, leaked_examples
    return True, leaked_examples


def report_order(order):
    hd(f"Z_{order} Center Gauss Orbits")
    orbits = orbit_decomposition(order)
    sizes = sorted(len(current) for current in orbits)
    fluxes = sorted({next(iter({flux(state, order) for state in current})) for current in orbits})
    print(f"  number of link-source sectors : {order ** 4}")
    print(f"  number of Gauss orbits        : {len(orbits)}")
    print(f"  orbit sizes                   : {sizes}")
    print(f"  orbit flux labels             : {fluxes}")

    for current in orbits:
        values = {flux(state, order) for state in current}
        assert len(values) == 1
    assert len(orbits) == order
    assert sizes == [order**3] * order
    assert fluxes == list(range(order))
    return orbits


def main():
    hd("Scope")
    print(
        "This is the open-link center/Gauss skeleton behind the candidate "
        "matter-coupled TCH operator. It checks closure of link-source sectors "
        "before compressing to one plaquette character variables."
    )

    report_order(2)
    report_order(3)

    hd("CSS Pauli Signs Versus SU(3) Center")
    closed_01, leaks_01 = binary_subset_closed(order=3, binary_value=1)
    closed_02, leaks_02 = binary_subset_closed(order=3, binary_value=2)
    print(f"  binary embedding {{0,1}}^4 closed under Z3 Gauss action: {closed_01}")
    for state, vertex, target in leaks_01:
        print(f"    leak: state={state}, vertex={vertex} -> {target}")
    print(f"  binary embedding {{0,2}}^4 closed under Z3 Gauss action: {closed_02}")
    for state, vertex, target in leaks_02:
        print(f"    leak: state={state}, vertex={vertex} -> {target}")

    assert not closed_01
    assert not closed_02

    hd("Candidate Operator Consequence")
    print(
        "  Z2/U(1)/SU(2) center shadow: the Pauli CSS closed-flux operator "
        "P_B=(1+B)/2 is algebraically consistent because the link-source sectors "
        "are closed under Z2 Gauss transformations."
    )
    print(
        "  SU(3): the analogous center constraint is Z3-valued, "
        "P_B=(1+B+B^dag)/3 with B carrying oriented plaquette center flux. "
        "A Pauli-sign CSS subset is not closed under the required Gauss action, "
        "so the actual item-13 operator must add a qutrit/clock center layer or "
        "an equivalent colour-orbit compensator."
    )
    print("\nALL ASSERTS PASSED.")


if __name__ == "__main__":
    main()
