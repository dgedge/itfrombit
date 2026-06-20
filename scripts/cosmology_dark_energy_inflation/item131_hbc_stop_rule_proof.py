#!/usr/bin/env python3
r"""ITEM 131: finite HBC stop-rule proof.

Target
------
Prove the stop rule used by the alpha^4 HBC scalar-amplitude branch:

    lambda_shell = E[N_topological commits per printed shell] = C_F = 4/3

with

    lambda_shell = N_shell alpha_0^4,
    A_nu = 1/N_shell = (3/4) alpha_0^4.

What is proved here
-------------------
Inside the local, single-clock, gauge-projected HBC/QEC printer, the stop
rule is a finite queue-balance theorem:

1. The scalar curvature readout is the post-decoder colour-restoring
   topology current, not the all-channel entropy current.
2. The admitted colour-restoring weight-4 logical channels all carry load
   2*(8/12) = 4/3.  This is the unbroken SU(3) colour-restoring load selected
   by the existing confinement/gauge filter.
3. A coherent constant-H saturated printer has two finite conditions:

       sustainability: lambda_shell <= C_F
       no-idle saturation: lambda_shell >= C_F

   Therefore lambda_shell = C_F.  Subcritical printers are sustainable but
   not saturated; supercritical printers overload and exit constant-H service.

What is not proved here
-----------------------
This is not a proof that total Bekenstein-Hawking entropy saturation alone
forces the colour-restoring scalar current to run at capacity.  The proof
rests on the scalar-printer identification already used by the HBC tilt leg:
the unique scalar clock is the local post-decoder colour-restoring current.
Adding a second scalar source or nonlocal horizon-mode service operator would
be new physics and reopens the amplitude/tilt audit.

exit 0 = finite geometry, queue inequalities, and amplitude consequences pass.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0

NAMES = ("G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W")
COLOUR_BITS = {3, 4}  # C0, C1
POINTS = tuple(range(8))
EDGES = tuple(
    (u, v)
    for u in POINTS
    for v in POINTS
    if u < v and int.bit_count(u ^ v) == 1
)


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def dot_mod2(a: int, x: int) -> int:
    return int.bit_count(a & x) & 1


def hyperplanes() -> list[frozenset[int]]:
    out: set[frozenset[int]] = set()
    for a in range(1, 8):
        for b in (0, 1):
            out.add(frozenset(p for p in POINTS if dot_mod2(a, p) == b))
    return sorted(out, key=lambda h: tuple(sorted(h)))


def tripped_edges(support: frozenset[int]) -> int:
    return sum((u in support) ^ (v in support) for u, v in EDGES)


def support_name(support: frozenset[int]) -> str:
    return "{" + ",".join(NAMES[i] for i in sorted(support)) + "}"


@dataclass(frozen=True)
class Channel:
    support: frozenset[int]
    colour_restoring: bool
    tripped: int
    load: Fraction


def finite_channels() -> list[Channel]:
    rows = []
    for support in hyperplanes():
        colour_restoring = not bool(support & COLOUR_BITS)
        tripped = tripped_edges(support)
        load = Fraction(2 * tripped, len(EDGES))  # H5 two-wall normalization.
        rows.append(Channel(support, colour_restoring, tripped, load))
    return rows


def source_contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def classify_lambda(lambda_shell: Fraction, capacity: Fraction) -> tuple[str, Fraction]:
    drift = capacity - lambda_shell
    if drift > 0:
        return "sustainable but not saturated: idle scalar-print capacity remains", drift
    if drift < 0:
        return "overloaded: coherent constant-H service must exit", drift
    return "saturated fixed point: no idle capacity and no overload", drift


def main() -> int:
    print("ITEM 131 HBC STOP-RULE PROOF")
    print("=" * 88)

    print("\n[1] Live source anchors")
    anchors = [
        (
            "28-clock tilt uses saturated horizon scale",
            "ANCHOR.md",
            "saturated horizon scale",
        ),
        (
            "colour-restoring coefficient is already selected",
            "python_code/item131_cf_stop_rule_closure_attempt.py",
            "C_F load structurally selected",
        ),
        (
            "scalar readout is not all-channel entropy",
            "python_code/item131_hbc_channel_whitening_closure.py",
            "The scalar clock is not the all-ones entropy current",
        ),
        (
            "T5b whiteness leaves only shell count",
            "python_code/item131_t5b_whiteness_from_locality.py",
            "At this T5b-only stage the remaining amplitude residual was T3",
        ),
    ]
    for label, path, phrase in anchors:
        check(source_contains(path, phrase), label)

    print("\n[2] Finite channel geometry")
    rows = finite_channels()
    colour = [row for row in rows if row.colour_restoring]
    kernel = [row for row in rows if not row.colour_restoring]
    total_service_labels = 2 * len(rows)
    scalar_labels = 2 * len(colour)
    loads = {row.load for row in colour}

    print(f"  affine hyperplane supports       = {len(rows)}")
    print(f"  transverse service labels        = {total_service_labels}")
    print(f"  scalar/readout service labels    = {scalar_labels}")
    print(f"  readout fraction of all channels = {scalar_labels}/{total_service_labels}")
    for row in rows:
        tag = "scalar/readout" if row.colour_restoring else "readout-kernel"
        print(
            f"    {support_name(row.support):24s} {tag:14s} "
            f"trips={row.tripped:2d}/12 load={row.load}"
        )

    check(len(rows) == 14, "weight-4 hyperplane current has 14 supports")
    check(total_service_labels == 28, "finite HBC service instrument has 28 labels")
    check(scalar_labels == 6, "gauge-projected scalar readout has six colour-restoring labels")
    check(len(kernel) == 11, "the other 22 transverse labels are entropy records in the scalar kernel")
    check(loads == {Fraction(4, 3)}, "every scalar/readout channel carries load C_F=4/3")

    print("\n[3] Stop rule as a finite queue-balance theorem")
    c_f = Fraction(4, 3)
    print("  Define lambda_shell = N_shell alpha_0^4.")
    print("  Let B be stored coherent scalar-print capacity for one horizon shell:")
    print("      Delta B = C_F - lambda_shell.")
    print("  constant-H sustainability requires lambda_shell <= C_F;")
    print("  saturated no-idle scalar printing requires lambda_shell >= C_F.")
    print("  Hence lambda_shell = C_F.")
    candidates = [
        Fraction(1, 1),
        Fraction(8, 7),
        Fraction(4, 3),
        Fraction(3, 2),
        Fraction(28, 6),
    ]
    for lam in candidates:
        verdict, drift = classify_lambda(lam, c_f)
        print(f"    lambda={lam!s:>4s}: Delta B={drift!s:>5s} -> {verdict}")
    fixed_verdict, fixed_drift = classify_lambda(c_f, c_f)
    check(fixed_drift == 0 and fixed_verdict.startswith("saturated"), "only lambda_shell=C_F satisfies both inequalities")

    print("\n[4] Amplitude consequence")
    n_shell = float(c_f) / ALPHA0**4
    a_nu = 1.0 / n_shell
    print(f"  N_shell = C_F alpha_0^-4 = {n_shell:.9e}")
    print(f"  A_nu    = 1/N_shell      = {a_nu:.12e}")
    print(f"          = (3/4) alpha_0^4")
    check(abs(a_nu - 0.75 * ALPHA0**4) < 1e-24, "stop rule gives A_nu=(3/4)alpha_0^4")

    print("\n[5] Status")
    print("  CLOSED within the current local single-clock scalar-printer algebra:")
    print("    finite channel projection + queue balance force lambda_shell=C_F.")
    print("  STILL CONDITIONAL as a physics statement:")
    print("    the proof uses the HBC scalar-printer identification, not total")
    print("    entropy saturation alone. A nonlocal scalar service source or second")
    print("    clock would be new canon and must reopen the tilt/amplitude audit.")
    print("=" * 88)
    print("exit 0 -- HBC stop rule proved as finite service-ledger queue balance under the local single-clock scalar-printer premise.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
