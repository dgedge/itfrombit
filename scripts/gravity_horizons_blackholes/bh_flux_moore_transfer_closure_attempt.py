#!/usr/bin/env python3
r"""Black-hole flux: can the Moore-alphabet transfer be made unavoidable?

Question
--------
The conditional 10/27 route says

    item-120 Moore shell + V_cell latch + outward horizon face
      -> (9+1)/(26+1) = 10/27.

This script tests the remaining epistemic gap: does horizon-local structure
force the full 26-slot Moore shell, or does it still require transferring the
item-120 Landauer-Moore service alphabet?

Result
------
Local cubic symmetry alone does not force the full Moore shell.  The
Chebyshev-nearest shell splits into three O_h orbits:

    F = 6 face contacts, E = 12 edge contacts, C = 8 corner contacts.

The monitored-channel connectedness criterion rules out some subsets, but not
all.  F+E and E+C are connected local alphabets with different flux
coefficients.  Therefore the full F+E+C Moore alphabet is not forced by
locality + cubic symmetry + monitored connectedness alone.

The exact closure target is now sharper:

    prove horizon severing is a closed-cell Landauer erasure across every
    non-empty nearest contact of the cubic cell.

That all-contact rule locally and uniquely selects F+E+C, after which the
existing 10/27 theorem follows.  Without it, the 10/27 route remains a strong
transfer theorem, not a locked black-hole-only derivation.
"""

from __future__ import annotations

from collections import deque
from itertools import combinations, product
import math


PHI = (math.sqrt(5.0) - 1.0) / 2.0
EPS_F = 1.0 / (2.0 * PHI)
ALPHA0 = 1.0 / 137.0
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def moore_shell() -> list[tuple[int, int, int]]:
    return [v for v in product((-1, 0, 1), repeat=3) if v != (0, 0, 0)]


def orbit_name(v: tuple[int, int, int]) -> str:
    nonzero = sum(1 for x in v if x != 0)
    return {1: "F", 2: "E", 3: "C"}[nonzero]


def orbit_union(names: set[str]) -> list[tuple[int, int, int]]:
    return [v for v in moore_shell() if orbit_name(v) in names]


def connected(nodes: list[tuple[int, int, int]]) -> bool:
    if not nodes:
        return False
    node_set = set(nodes)
    seen = {nodes[0]}
    queue: deque[tuple[int, int, int]] = deque([nodes[0]])
    while queue:
        a = queue.popleft()
        for b in node_set:
            if b not in seen and sum(abs(a[i] - b[i]) for i in range(3)) == 1:
                seen.add(b)
                queue.append(b)
    return len(seen) == len(nodes)


def mean_f(beta_eff: float) -> float:
    raw = {f: g * math.exp(-beta_eff * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return sum(f * value / z for f, value in raw.items())


def required_gamma_for_stefan(beta_eff: float = 1.0) -> float:
    rho_lambda = beta_eff / (2.0 * math.pi * EPS_F)
    target_c = 1.0 / (15360.0 * math.pi)
    denom = (27.0 * math.pi / 32.0) * EPS_F * mean_f(beta_eff) * rho_lambda**4
    return target_c / denom


def model_over_target(prefactor: float) -> float:
    return prefactor * ALPHA0 / required_gamma_for_stefan()


def emitted_prefactor(nodes: list[tuple[int, int, int]], latch: bool = True) -> tuple[int, int, float]:
    outward = [v for v in nodes if v[2] == 1]
    numerator = len(outward) + (1 if latch else 0)
    denominator = len(nodes) + (1 if latch else 0)
    return numerator, denominator, numerator / denominator


def all_orbit_subsets() -> list[set[str]]:
    names = ("F", "E", "C")
    out: list[set[str]] = []
    for r in range(1, len(names) + 1):
        out.extend(set(c) for c in combinations(names, r))
    return out


def main() -> None:
    print("BLACK-HOLE FLUX: MOORE-TRANSFER CLOSURE ATTEMPT")
    print("=" * 92)

    shell = moore_shell()
    orbit_counts = {name: sum(1 for v in shell if orbit_name(v) == name) for name in ("F", "E", "C")}
    print("[1] Local cubic shell orbit inventory")
    for name, count in orbit_counts.items():
        print(f"    {name}: {count} contacts")
    check(orbit_counts == {"F": 6, "E": 12, "C": 8}, "Moore shell splits as 6 face + 12 edge + 8 corner contacts")

    print("\n[2] All O_h-invariant nearest-contact alphabets")
    rows = []
    for names in all_orbit_subsets():
        nodes = orbit_union(names)
        num, den, pref = emitted_prefactor(nodes, latch=True)
        rows.append((tuple(sorted(names)), len(nodes), connected(nodes), num, den, pref, model_over_target(pref)))
    for names, size, is_conn, num, den, pref, ratio in rows:
        print(
            f"    {''.join(names):<3s} size={size:2d} connected={str(is_conn):5s} "
            f"emit={num:2d}/{den:2d} pref={pref:.9f} P/P_SB={ratio:.6f}"
        )

    connected_rows = [row for row in rows if row[2]]
    connected_names = [row[0] for row in connected_rows]
    check(len(rows) == 7, "local cubic symmetry permits seven non-empty orbit-union alphabets")
    check(("E", "F") in connected_names, "F+E is a connected monitored local alphabet")
    check(("C", "E") in connected_names, "E+C is a connected monitored local alphabet")
    check(("C", "E", "F") in connected_names, "full Moore F+E+C is connected")
    check(len(connected_rows) > 1, "connectedness does not uniquely select full Moore")

    print("\n[3] Full Moore is selected only by an all-contact severing rule")
    full = orbit_union({"F", "E", "C"})
    num, den, pref = emitted_prefactor(full, latch=True)
    print(f"    all-contact alphabet: emit={num}/{den}, pref={pref:.12f}, P/P_SB={model_over_target(pref):.9f}")
    check(len(full) == 26, "all non-empty nearest contacts give the 26-slot Moore shell")
    check((num, den) == (10, 27), "outward face plus latch gives 10/27 only for the full shell")
    check(abs(pref - 10.0 / 27.0) < 1.0e-15, "full Moore transfer reproduces exactly 10/27")
    check(abs(model_over_target(pref) - 1.0) < 4.0e-3, "10/27 alpha is within 0.4 percent of the Stefan-Hawking coefficient")

    print(
        """
[4] VERDICT
    Not closed from existing local data alone.

    Locality + cubic symmetry + monitored connectedness still leave multiple
    admissible local alphabets.  The full 26-slot Moore shell is forced only if
    horizon severing is proven to be closed-cell Landauer erasure through every
    non-empty nearest contact of a cubic cell.

    Therefore the 10/27 coefficient is not a theorem of V_cell/V_Sch alone and
    not yet a black-hole-only locked coefficient.  The route is, however, now
    maximally sharp:

        all-contact local Landauer severing + V_cell latch
          -> 26+1 slots, 9+1 emitted slots
          -> Gamma0 = (10/27) alpha0 Lambda_QCD.

    The remaining proof obligation is exactly the all-contact severing rule.
ALL ASSERTIONS PASSED -- transfer premise sharpened; not eliminated."""
    )


if __name__ == "__main__":
    main()
