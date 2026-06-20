#!/usr/bin/env python3
r"""Black-hole horizon severing: all-contact theorem audit.

Question
--------
Can current canon prove

    horizon severing is closed-cell Landauer erasure through every non-empty
    nearest contact?

This is the exact remaining premise behind the black-hole 10/27 flux route.
The earlier Moore-transfer audit showed that

    all non-empty nearest contacts + latch -> (9+1)/(26+1) = 10/27.

This script tests whether the *existing* black-hole objects force that
all-contact rule, or whether they remain compatible with narrower surface /
radial-contact alphabets.

Verdict
-------
The current canon proves a strong horizon record channel:

* V_cell is a finite isometry over the 12-edge Q3 strain syndrome plus latch.
* V_Sch(M) is a direct sum over radial shell cells/bonds crossing the frozen
  horizon shell.
* The severing ledger is partner-asymmetric and radial/surface-local.

But those objects do not contain a spatial contact alphabet.  Tensoring V_cell
with any orthogonal spatial contact label preserves isometry.  Local cubic
symmetry and connectedness still admit multiple nearest-contact alphabets, not
only the full Moore shell.

Therefore the all-contact statement is not proved from current canon.  The
precise theorem still needed is a service-class identification:

    horizon severing is a closed-cell scalar Landauer event over every
    non-empty nearest contact, rather than a codimension-one radial pair
    severing event.

If that theorem is added or derived, the 10/27 coefficient follows immediately.
Without it, 10/27 remains a strong Moore-transfer premise, not a locked
black-hole-only coefficient.
"""

from __future__ import annotations

from collections import Counter, deque
from itertools import combinations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANCHOR = (ROOT / "ANCHOR.md").read_text(encoding="utf-8")
DRIFT = (ROOT / "DRIFT.md").read_text(encoding="utf-8")


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def text_has(*needles: str) -> bool:
    corpus = ANCHOR + "\n" + DRIFT
    return all(needle in corpus for needle in needles)


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


def all_orbit_subsets() -> list[set[str]]:
    names = ("F", "E", "C")
    out: list[set[str]] = []
    for r in range(1, len(names) + 1):
        out.extend(set(c) for c in combinations(names, r))
    return out


def emitted_prefactor(names: set[str], latch: bool = True) -> tuple[int, int]:
    nodes = orbit_union(names)
    outward = [v for v in nodes if v[2] == 1]
    return len(outward) + int(latch), len(nodes) + int(latch)


def main() -> None:
    print("BLACK-HOLE ALL-CONTACT SEVERING THEOREM AUDIT")
    print("=" * 88)

    print("[1] Canon anchors currently present")
    anchors = {
        "V_cell finite isometry": text_has("V_cell"),
        "V_Sch radial shell isometry": (
            text_has("V_{\\rm Sch}(M)", "direct-summed over frozen-coin shell cells")
            or text_has("V_{\\rm Sch}", "direct sum over frozen-coin horizon-shell cells")
            or text_has("V_Sch(M) is the direct-sum isometry")
        ),
        "12-edge strain syndrome": text_has("12-edge strain syndrome"),
        "partner asymmetry": text_has("partner-asymmetric"),
        "all-contact residual named": text_has("all-contact"),
    }
    for name, ok in anchors.items():
        check(ok, name)

    print("\n[2] Contact alphabets allowed by cubic nearest-neighbour locality")
    counts = Counter(orbit_name(v) for v in moore_shell())
    check(counts == {"F": 6, "E": 12, "C": 8}, "Moore shell splits into face/edge/corner orbits 6/12/8")

    rows: list[tuple[str, int, bool, int, int]] = []
    for names in all_orbit_subsets():
        label = "".join(sorted(names))
        num, den = emitted_prefactor(names)
        rows.append((label, len(orbit_union(names)), connected(orbit_union(names)), num, den))
    for label, size, is_connected, num, den in rows:
        print(f"    {label:<3s} size={size:2d} connected={str(is_connected):5s} emit={num:2d}/{den:2d}")

    connected_labels = {label for label, _, is_connected, _, _ in rows if is_connected}
    check("EF" in connected_labels, "face+edge is a connected local alphabet")
    check("CE" in connected_labels, "edge+corner is a connected local alphabet")
    check("CEF" in connected_labels, "full Moore face+edge+corner is connected")
    check(len(connected_labels) > 1, "locality + connectedness do not uniquely force full Moore")

    print("\n[3] Why V_cell/V_Sch do not by themselves select full Moore")
    print("    V_cell record labels: syndrome class + latch; no spatial contact coordinate.")
    print("    V_Sch shell labels: radial shell cell/bond membership; no face/edge/corner orbit selector.")
    print("    Any orthogonal contact alphabet can be tensored onto V_cell without spoiling isometry.")
    print("    Thus the isometry proof is compatible with EF, CE, CEF, etc.; it cannot select CEF.")

    full_num, full_den = emitted_prefactor({"C", "E", "F"})
    check((full_num, full_den) == (10, 27), "closed-cell all-contact rule gives exactly 10/27 with latch")

    print(
        """
[4] VERDICT
    No: the all-contact rule is not proved by current canon.

    What current canon proves is a radial/shell horizon record channel:
    V_cell plus the 12-edge strain syndrome, lifted by V_Sch(M) across the
    frozen horizon shell.  That is enough for a finite KMS ladder and a
    standard greybody transfer, but it does not force a closed-cell Moore
    contact alphabet.

    The exact missing theorem is now:

        horizon severing is a scalar closed-cell Landauer event over every
        non-empty nearest contact of the shell cell.

    If that theorem is proved, the 10/27 flux coefficient follows.  If horizon
    severing remains only a codimension-one radial pair severing, the 10/27
    route remains conditional.

ALL ASSERTIONS PASSED -- proof obstruction explicit; closure criterion exact.
"""
    )


if __name__ == "__main__":
    main()
