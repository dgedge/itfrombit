#!/usr/bin/env python3
r"""Foundation disclosed-input count K.

Purpose
-------
Restate "parameter-free" honestly as:

    no continuous fitted parameters; K discrete structural postulates; D scale anchors.

This is an accounting ledger, not a new physics derivation.  It follows the
E2/item-99 input-vs-derived audits:

* item 99: 4.8.8 vertex figure is derived from five disclosed inputs plus one
  presentation convention, not taken as a root geometry.
* E2: charge weights and the colour one-hot normalisation are anomaly-forced for
  charge bookkeeping, not hand-set root inputs.
* item 133: generation count is input through R1, not derived from |O_h|/16.

The strict public count is the atomic count below.  The bundled count is supplied
only to compare with older prose that treats R1-R4 as one "constraint set" and
S*C as one "walk operator".
"""

from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANCHOR = (ROOT / "ANCHOR.md").read_text(encoding="utf-8")
DRIFT = (ROOT / "DRIFT.md").read_text(encoding="utf-8")


@dataclass(frozen=True)
class Entry:
    key: str
    statement: str
    source_checks: tuple[str, ...]
    bundle: str


POSTULATES = [
    Entry(
        "K01",
        "Binary/qubit substrate: the microscopic state alphabet is a finite two-state register.",
        ("3-regular qubit network", "8-qubit register"),
        "substrate",
    ),
    Entry(
        "K02",
        "Graph locality / spatial extent: the register lives on a network with adjacency and extended tiling.",
        ("spatial extent", "tile Euclidean space", "inter-cluster gauge bridges"),
        "substrate",
    ),
    Entry(
        "K03",
        "Trivalence / degree-3 budget of the pre-geometric network.",
        ("3-regular qubit network", "Trivalence"),
        "substrate",
    ),
    Entry(
        "K04",
        "Topological crystallisation rule: favour 4-cycles and 6-cycles while penalising degree defects.",
        ("4-cycles", "6-cycles", "deg(v) - 3"),
        "substrate",
    ),
    Entry(
        "K05",
        "Eight-face/octet matter register: one 8-bit code cell with Q3/octal addressing.",
        ("8 triangular faces", "000 → G₀", "111 → W"),
        "register",
    ),
    Entry(
        "K06",
        "Fixed semantic bit ordering / 4+4 partition of the 8-bit register.",
        ("000 → G₀", "100 → C₁", "4+4 bit partition"),
        "register",
    ),
    Entry(
        "K07",
        "R1 generation constraint: G0 * G1 != 1.",
        ("R1", "G", "limits to three generations"),
        "constraints",
    ),
    Entry(
        "K08",
        "R2 weak/chirality lock: W = chi.",
        ("R2", "W = χ", "chirality locks"),
        "constraints",
    ),
    Entry(
        "K09",
        "R3 lepton/colour constraint: colourless iff LQ = 0.",
        ("R3", "LQ = 0", "colour separates quarks from leptons"),
        "constraints",
    ),
    Entry(
        "K10",
        "R4 right-handed-neutrino exclusion constraint.",
        ("R4", "right-handed-neutrino exclusion", "forbidding the three right-handed-neutrino"),
        "constraints",
    ),
    Entry(
        "K11",
        "Weak coin: zero-controlled CNOT with control chi and target I3.",
        ("zero-controlled CNOT", "flip I", "χ = 0"),
        "walk",
    ),
    Entry(
        "K12",
        "Shift/bridge rule: eight body-diagonal face-complement channels in one O_h orbit.",
        ("eight body-diagonal directions", "face", "joins face"),
        "walk",
    ),
    Entry(
        "K13",
        "Bitwise CPT / global matter-antimatter doublet tag.",
        ("bitwise inversion", "Z_p", "matter/antimatter"),
        "cpt",
    ),
]


CONVENTIONS = [
    Entry(
        "C01",
        "R4-edge axis anchoring convention used for the horizontal/vertical presentation of item 99.",
        ("R4-edge reading", "axis-anchoring", "presentation-level"),
        "presentation",
    ),
]


BUNDLED_POSTULATES = [
    "B01 binary/qubit substrate",
    "B02 graph locality / spatial extent",
    "B03 trivalence / degree-3 budget",
    "B04 topological crystallisation rule",
    "B05 eight-face/octet register",
    "B06 fixed semantic bit ordering / 4+4 partition",
    "B07 R1-R4 constraint set",
    "B08 walk operator S*C",
    "B09 bitwise CPT / matter-antimatter doublet",
]


NOT_ROOT_INPUTS = [
    "4.8.8 vertex figure",
    "electric-charge weights {1/2,1/3,1/2}",
    "colour one-hot normalisation for charge bookkeeping",
    "three generations from |O_h|/16",
    "bare alpha_0^-1 = 137",
]


def check_sources() -> None:
    corpus = ANCHOR + "\n" + DRIFT
    missing = []
    for entry in POSTULATES + CONVENTIONS:
        for needle in entry.source_checks:
            if needle not in corpus:
                missing.append((entry.key, needle))
    if missing:
        detail = "\n".join(f"  {key}: {needle!r}" for key, needle in missing)
        raise AssertionError(f"source check failed:\n{detail}")


def main() -> None:
    check_sources()
    atomic_k = len(POSTULATES)
    bundled_k = len(BUNDLED_POSTULATES)
    conventions = len(CONVENTIONS)
    scale_anchors = 1  # ANCHOR 16.2 Single Dimensionful Anchor Principle.

    print("FOUNDATION DISCLOSED-INPUT COUNT")
    print("=" * 42)
    print(f"K_atomic  = {atomic_k}")
    print(f"K_bundle  = {bundled_k}")
    print(f"C         = {conventions} presentation convention")
    print(f"D         = {scale_anchors} dimensionful scale anchor (outside K)")
    print()
    for entry in POSTULATES:
        print(f"{entry.key} [{entry.bundle:11s}] {entry.statement}")
    print()
    print("Bundled postulates:")
    for statement in BUNDLED_POSTULATES:
        print(f"  - {statement}")
    print()
    print("Not counted as root inputs:")
    for statement in NOT_ROOT_INPUTS:
        print(f"  - {statement}")
    print()
    print(
        "Canonical phrasing: no continuous fitted parameters; "
        "K=13 disclosed discrete structural postulates (or K=9 bundled modules) "
        "plus D=1 disclosed dimensionful anchor. The G3 annealing coefficients are "
        "not used as fitted real parameters in this count; if a later result needs "
        "specific w4:w6:lambda values, those become additional continuous inputs."
    )


if __name__ == "__main__":
    main()
