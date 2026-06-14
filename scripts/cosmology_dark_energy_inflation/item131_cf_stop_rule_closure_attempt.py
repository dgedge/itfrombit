#!/usr/bin/env python3
"""ITEM 131 / 42: coefficient and stop-rule closure attempt.

Target
------
After the gauge-filtered R2/R3 route, the remaining threshold equation is

    N_print(H_*) alpha_0^4 = C_F,      C_F = 4/3.

This audit attacks the two remaining clauses:

  1. coefficient: why the threshold couples to the unbroken SU(3)
     colour-restoring load C_F=4/3, not the broken Pati-Salam loads;
  2. stopping rule: why the expected number of topology-changing commits per
     printed shell equals that load.

Result
------
The coefficient is substantially improved by the actual weight-4
Kraus/current geometry:

  * after the confinement/gauge filter, the admitted colour-silent logical
    channels all trip 8 of the 12 strain checks;
  * the canon string-tension load uses the two-wall factor Delta=2 times the
    violated-check fraction, so Delta*(8/12)=4/3;
  * this equals the SU(3) fundamental Casimir computed independently.

That is a real structural selection of C_F for the colour-restoring load.
It also explains why the broken Pati-Salam loads are the wrong object: they
belong to coset-changing transitions, while the post-filtered threshold is a
confinement/restoring load on the unbroken colour triplet.

The stop rule does NOT close from the finite Kraus geometry alone.  The
finite table supplies the event current and the load per admitted event, but
it does not supply a dynamical marginality theorem saying that the saturated
printer exits exactly when the Poisson mean topological action equals the
restoring load.  That remains one sharply named lemma:

    Marginal colour-load balance:
        lambda_shell = E[N_topological] = C_F.

If that lemma is adopted/proved, the threshold equation closes.  Without it,
nearby internally meaningful stop rules give different constants.
"""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)
COLOUR_BITS = {C0, C1}
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


def bits3(x: int) -> tuple[int, int, int]:
    return ((x >> 2) & 1, (x >> 1) & 1, x & 1)


def dot(a: int, x: int) -> int:
    return int.bit_count(a & x) & 1


def hyperplanes() -> list[frozenset[int]]:
    out = set()
    for a in range(1, 8):
        for b in (0, 1):
            out.add(frozenset(p for p in POINTS if dot(a, p) == b))
    return sorted(out, key=lambda h: tuple(sorted(h)))


def support_name(support: frozenset[int]) -> str:
    return "{" + ",".join(NAMES[i] for i in sorted(support)) + "}"


def tripped_edges(support: frozenset[int]) -> int:
    """Number of Q3 strain edges crossing the support boundary."""
    return sum((u in support) ^ (v in support) for u, v in EDGES)


def su_n_fundamental_casimir(n: int) -> F:
    return F(n * n - 1, 2 * n)


def source_contains(path: str, needle: str) -> bool:
    return needle in (ROOT / path).read_text(encoding="utf-8")


def ps_broken_loads() -> dict[str, F]:
    c_su2 = su_n_fundamental_casimir(2)
    c_su3 = su_n_fundamental_casimir(3)
    c_su4 = su_n_fundamental_casimir(4)
    u1_quark = F(1, 24)  # T15^2 for quark entries in SU(4)->SU(3)xU(1)
    u1_lepton = F(3, 8)
    su4_broken_quark = c_su4 - c_su3 - u1_quark
    su4_broken_lepton = c_su4 - u1_lepton
    su4_broken_average = (3 * su4_broken_quark + su4_broken_lepton) / 4
    su2r_broken_right = c_su2 - F(1, 4)
    su2r_broken_spinor_average = su2r_broken_right / 2
    return {
        "SU(3) fundamental C_F": c_su3,
        "SU(4)/SU(3)xU1 broken, quark": su4_broken_quark,
        "SU(4)/SU(3)xU1 broken, lepton": su4_broken_lepton,
        "SU(4) broken average over 4": su4_broken_average,
        "SU(2)_R broken on right doublet": su2r_broken_right,
        "SU(2)_R broken spinor average": su2r_broken_spinor_average,
        "full PS broken average over 16": su4_broken_average + su2r_broken_spinor_average,
    }


def main() -> None:
    print("ITEM 131 / 42 C_F COEFFICIENT AND STOP-RULE CLOSURE ATTEMPT")

    print("\n[1] Source checks")
    check(source_contains("python_code/item126_weight4_event_algebra.py", "the 3 colourless channels all"),
          "weight-4 event algebra identifies the confinement-unvetoed colour-silent channels")
    check(source_contains("python_code/item126_weight4_event_algebra.py", "trips exactly 8"),
          "event algebra records that the colourless channels trip exactly 8 strain checks")
    check(source_contains("DRIFT.md", "sqrtσ=(4/3)Λ") or source_contains("DRIFT.md", "sqrt\\sigma"),
          "canon carries the string/Casimir load entry")
    check(source_contains("DRIFT.md", "2·(2/3) = 4/3"),
          "H5 records the canonical two-wall times violation-fraction load")
    check(source_contains("python_code/item131_r2r3_gauge_filtered_threshold_route.py", "post-decoder topology-changing current"),
          "previous audit narrowed alpha^4 to the post-decoder topology-changing current")

    print("\n[2] Actual weight-4 strain-load geometry")
    planes = hyperplanes()
    check(len(planes) == 14 and all(len(h) == 4 for h in planes), "there are 14 weight-4 hyperplane logical supports")
    loads = []
    for h in planes:
        trips = tripped_edges(h)
        colour_silent = not bool(h & COLOUR_BITS)
        load = F(2 * trips, len(EDGES))  # Delta=2 times violated-check fraction.
        loads.append((h, colour_silent, trips, load))
    count_by_trips = Counter(trips for _, _, trips, _ in loads)
    check(count_by_trips == Counter({4: 6, 8: 6, 12: 2}), "all-channel strain trips split as 6x4, 6x8, 2x12")
    colour_silent = [(h, trips, load) for h, silent, trips, load in loads if silent]
    check(len(colour_silent) == 3, "three colour-silent logical channels survive the colour firewall")
    check(all(trips == 8 for _, trips, _ in colour_silent), "each colour-silent channel trips 8 of 12 strain checks")
    check(all(load == F(4, 3) for _, _, load in colour_silent), "two-wall strain-load normalization gives 2*(8/12)=4/3")
    for h, silent, trips, load in loads:
        tag = "colour-silent" if silent else "coloured"
        print(f"  {support_name(h):24s} {tag:14s} trips={trips:2d}/12  Delta*p={load}={float(load):.9f}")
    all_channel_mean_load = sum(load for *_, load in loads) / len(loads)
    print(f"  all-channel mean load             = {all_channel_mean_load} = {float(all_channel_mean_load):.9f}")
    print("  confinement-filtered load per admitted channel = 4/3")
    check(all_channel_mean_load == F(8, 7), "averaging all 14 channels would give 8/7, not 4/3")

    print("\n[3] Representation cross-check: SU(3) Casimir versus broken PS loads")
    loads_ps = ps_broken_loads()
    for label, value in loads_ps.items():
        print(f"  {label:42s} = {value} = {float(value):.9f}")
    check(loads_ps["SU(3) fundamental C_F"] == F(4, 3), "SU(3) fundamental Casimir equals 4/3")
    check(loads_ps["SU(4) broken average over 4"] == F(3, 4), "SU(4) broken-coset average is 3/4")
    check(loads_ps["SU(2)_R broken on right doublet"] == F(1, 2), "SU(2)_R broken load is 1/2")
    check(loads_ps["full PS broken average over 16"] == F(1, 1), "full broken PS average is 1")
    print("  Interpretation: the actual confinement-filtered strain geometry selects")
    print("  the same object as the unbroken SU(3) fundamental Casimir.  Broken-coset")
    print("  loads describe transitions through the Pati-Salam-invalid sector, not")
    print("  the post-filter colour-restoring threshold load.")

    print("\n[4] Stop-rule alternatives")
    c_f = F(4, 3)
    alternatives = [
        ("raw one logical event", F(1, 1), "N alpha^4 = 1"),
        ("colour-load marginality", c_f, "N alpha^4 = C_F"),
        ("unit-normalized colour load", F(1, 1) / c_f, "C_F N alpha^4 = 1"),
        ("all-channel mean strain load", all_channel_mean_load, "N alpha^4 = <Delta p>_14"),
        ("colour-silent branch count if alpha^4 means all 14 commits", F(14, 3) * c_f, "(3/14) N alpha^4 = C_F"),
    ]
    for label, constant, equation in alternatives:
        amp_coeff = F(1, 1) / constant
        print(
            f"  {label:54s} C={constant}={float(constant):.9f}  "
            f"A coeff=1/C={amp_coeff}={float(amp_coeff):.9f}  [{equation}]"
        )
    check(alternatives[1][1] == F(4, 3), "the desired threshold is the colour-load marginality rule")
    check(len({constant for _, constant, _ in alternatives}) == len(alternatives), "nearby meaningful stop rules are distinct")
    print("  The finite Kraus/current geometry does not choose among these stop rules.")
    print("  It supplies the load; it does not by itself supply the dynamical")
    print("  marginality law saying that the Poisson mean shell action equals the load.")

    print("\n[5] Gate verdict")
    gates = [
        (
            "alpha^4 current",
            "CONDITIONAL-CLOSED",
            "post-decoder topology-changing current from previous audit",
        ),
        (
            "C_F coefficient",
            "SUBSTANTIALLY CLOSED",
            "confinement-filtered colour-silent logicals give 2*(8/12)=4/3 and match SU(3) C_F",
        ),
        (
            "broken-coset exclusion",
            "CLOSED UNDER FILTER",
            "broken PS loads belong to coset-changing transitions, not the colour-restoring post-filter load",
        ),
        (
            "stop rule",
            "OPEN",
            "requires marginal colour-load balance: E[N_topological]=C_F",
        ),
    ]
    for name, status, note in gates:
        print(f"  {name:24s} {status:22s} {note}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  The coefficient clause can now be promoted one notch: it is no longer")
    print("  just 'R3 has an SU(3) triplet'.  In the actual weight-4 Kraus/current")
    print("  geometry, the confinement-unvetoed colour-silent logicals trip 8/12")
    print("  strain checks; the canonical two-wall string/load normalization gives")
    print("  2*(8/12)=4/3, exactly the SU(3) fundamental Casimir.")
    print("")
    print("  The stop rule still does not honestly close from finite geometry alone.")
    print("  To finish the threshold one must prove the marginal colour-load balance")
    print("  lemma:")
    print("      lambda_shell = E[N_topological commits per printed shell] = C_F.")
    print("  With that lemma, N_print alpha_0^4 = 4/3 follows.  Without it, several")
    print("  internally meaningful thresholds give different constants.  So the sticky")
    print("  problem is reduced to one dynamical marginality theorem, not solved away.")
    print("=" * 108)
    print("exit 0 -- C_F load structurally selected; stop-rule marginality remains open.")


if __name__ == "__main__":
    main()
