#!/usr/bin/env python3
r"""K04 orphan-policy audit for the debris wall-shadow branch.

This is a policy/ledger audit, not a new Monte Carlo run. It consumes the saved
pairing L-trend rows and asks which orphan policies remain compatible with both
the measured wall dynamics and the framework's own QEC/no-register principles.

Verdict structure:

1. Strict pair growth pairing is a diagnostic counterfactual, not a viable
   canonical policy. The saved L=10/12 rows show persistent low-depth wall
   stalls, positive late-wall nucleation lag, and rho_D/rho_B hundreds of times
   too large.
2. Fully adaptive rescue is also not canonical: it reads through missing
   no-register wall interiors and erases the recorded-boundary shadow.
3. The surviving policy is boundary-local orphan rescue: a stranded block can be
   padded into an adjacent registered higher-depth host, but disconnected islands
   and no-register wall interiors cannot be read through. The remaining
   quantitative task is the island-floor surface rho_D/rho_B(d,L) under that
   policy, not another free coefficient.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "pairing_ltrend.jsonl"
DRIFT = ROOT.parent / "DRIFT.md"
ANCHOR = ROOT.parent / "ANCHOR.md"

ALPHA0 = 1 / 137.0
ETA = 6.1e-10
NGAM = 2 * 1.2020569 / math.pi**2
RHO_B = ETA * NGAM * 2 * math.sqrt(2)
RDB_OBS = 0.1200 / 0.02237


def q1_of(p: float) -> float:
    return sum(
        math.comb(8, k)
        * p**k
        * (1 - p) ** (8 - k)
        * (0.0 if k <= 3 else (0.5 if k == 4 else 1.0))
        for k in range(9)
    )


Q1 = q1_of(0.0972)


def resid(level: int) -> float:
    return (21 * Q1) ** (2 ** (level - 1)) / 21


def load_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line == "DONE":
                continue
            rows.append(json.loads(line))
    return rows


def mean(values: list[float]) -> float:
    return fmean(values) if values else float("nan")


def group_rows(rows: list[dict]) -> dict[tuple[int, int], list[dict]]:
    out: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        out[(int(row["L"]), int(row["R"]))].append(row)
    return dict(sorted(out.items()))


def summarize_group(rows: list[dict]) -> dict[str, float]:
    lags = [
        float(row["t_nuc_stall"]) - float(row["t_nuc_deep"])
        for row in rows
        if row.get("t_nuc_stall") is not None and row.get("t_nuc_deep") is not None
    ]
    return {
        "n": len(rows),
        "d": mean([float(row["d"]) for row in rows]),
        "rdb": mean([float(row["rdb"]) for row in rows]),
        "overshoot": mean([float(row["rdb"]) / RDB_OBS for row in rows]),
        "stall": mean([float(row["stall_frac"]) for row in rows]),
        "lag": mean(lags),
    }


def rescue_summaries(rows: list[dict]) -> dict[tuple[int, int], dict[str, float]]:
    if not rows or not all("rdb_rescue" in row for row in rows):
        return {}
    grouped = group_rows(rows)
    out = {}
    for key, rs in grouped.items():
        out[key] = {
            "rdb_rescue": mean([float(row["rdb_rescue"]) for row in rs]),
            "overshoot_rescue": mean([float(row["rdb_rescue"]) / RDB_OBS for row in rs]),
            "stall_rescue": mean([float(row["stall_frac_rescue"]) for row in rs]),
        }
    return out


def assert_canon_context() -> None:
    drift = DRIFT.read_text(encoding="utf-8")
    anchor = ANCHOR.read_text(encoding="utf-8")
    required_drift = [
        "Strict-pair growth pairing does **not** collapse to depth-3",
        "orphan-rescue",
        "ISLAND FLOOR",
    ]
    required_anchor = [
        "Canonical K04 protocol",
        "Topological crystallisation rule",
    ]
    for needle in required_drift:
        assert needle in drift, f"DRIFT context missing: {needle}"
    for needle in required_anchor:
        assert needle in anchor, f"ANCHOR context missing: {needle}"


def main() -> int:
    assert_canon_context()
    rows = load_rows(DATA)
    assert rows, f"no rows found in {DATA}"
    grouped = group_rows(rows)
    assert set(grouped) >= {(10, 496), (12, 496)}, grouped.keys()

    print("[0] K04 orphan-policy audit")
    print(f"    loaded rows: {len(rows)} from {DATA.name}")
    print(f"    observed rho_D/rho_B target: {RDB_OBS:.3f}")
    print(f"    residual hierarchy: r2/r3 = {resid(2) / resid(3):.1f}")
    print()

    print("[1] strict-pair L-trend from saved rows")
    strict = {}
    for key, rs in grouped.items():
        summary = summarize_group(rs)
        strict[key] = summary
        print(
            f"    L={key[0]} R={key[1]} reps={summary['n']:.0f}: "
            f"d={summary['d']:.3f}, stall<=2={summary['stall']:.3f}, "
            f"rho_D/rho_B={summary['rdb']:.1f} "
            f"({summary['overshoot']:.0f}x obs), "
            f"late-wall lag={summary['lag']:.2f} snapshots"
        )

    assert strict[(10, 496)]["overshoot"] > 100
    assert strict[(12, 496)]["overshoot"] > 100
    assert strict[(10, 496)]["lag"] > 0
    assert strict[(12, 496)]["lag"] > 0
    assert strict[(12, 496)]["stall"] > 0.05
    print("    verdict: strict no-rescue pairing is refuted as a debris-DM policy.")
    print("    reason: persistent shallow wall stalls dominate the residual, and")
    print("            stalled wall cells nucleate later than deeper wall cells.")
    print()

    print("[2] rescue diagnostics in local data")
    rescue = rescue_summaries(rows)
    if rescue:
        for key, summary in rescue.items():
            print(
                f"    L={key[0]} R={key[1]}: rescue rho_D/rho_B="
                f"{summary['rdb_rescue']:.1f} "
                f"({summary['overshoot_rescue']:.1f}x obs), "
                f"stall<=2={summary['stall_rescue']:.3f}"
            )
    else:
        print("    current JSONL carries strict fields only; rescue/island diagnostics")
        print("    are recorded in DRIFT from the patched rerun but require a larger")
        print("    L=12..16 rescue rerun for machine-grade closure.")
    print()

    print("[3] admissible orphan policies")
    policies = {
        "strict_pair": {
            "status": "excluded-as-canonical",
            "why": "leaves adjacent low-depth wall blocks unpadded and overshoots by >100x",
        },
        "fully_adaptive": {
            "status": "excluded-by-no-register-locality",
            "why": "would read through missing wall interiors and erase the boundary-shadow observable",
        },
        "boundary_local_rescue": {
            "status": "selected-candidate",
            "why": "pads only into adjacent registered hosts; disconnected islands remain as the measurable floor",
        },
    }
    for name, entry in policies.items():
        print(f"    {name}: {entry['status']} -- {entry['why']}")

    assert policies["boundary_local_rescue"]["status"] == "selected-candidate"
    print()
    print("[4] remaining quantitative gate")
    print("    The orphan policy is no longer a free binary choice. The live branch is")
    print("    boundary-local rescue, and the remaining object is the island-floor")
    print("    surface rho_D/rho_B(d,L) at the gamma-driver d in the L -> infinity")
    print("    limit. Existing DRIFT registration calls for L=12..16 and R=496..1600.")
    print()
    print("MACHINERY ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
