#!/usr/bin/env python3
r"""K04 cooling-driver and energy-fraction bookkeeping audit.

This script consumes the embedded KZ sweep, then evaluates two cooling-rate
branches without re-running the Monte Carlo:

1. Canon/paper geometric anneal:
   crystallisation.tex pins gamma = 0.995 per sweep, so the KZ worker's
   T: 6 -> 0.5 ramp corresponds to
       R = ln(0.5 / 6) / ln(0.995).

2. Conditional expansion branch:
   if the boot substrate temperature redshifts as T proportional to a^-1 and
   the HBC/QEC radial print clock is identified with 28 ticks per d ln a, then
       R = 28 ln(6 / 0.5).
   If the bipartite two-step physical cycle is charged, use 56 ticks per d ln a.

The second branch is a clock-map estimate, not a canon-locked derivation.  The
KZ exponents in the underlying data are finite-size/toy-lattice slopes, so all
interpolated d_trapped values are lattice-unit estimates.

The final section sets up, but does not close, the vertex-fraction to energy-
fraction map shared by the debris-DM and R4/nu_R budget branches.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import fmean, pstdev


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "k04_kz_results.jsonl"

T_START = 6.0
T_END = 0.5
GAMMA_CANON = 0.995
DELTA_N = math.log(T_START / T_END)
REFERENCE_MATTER_DARK_FRACTION = 0.84


def load_rows(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def grouped_means(rows: list[dict]) -> dict[int, dict[int, dict[str, float]]]:
    grouped: dict[tuple[int, int], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(int(row["L"]), int(row["sweeps"]))].append(row)

    out: dict[int, dict[int, dict[str, float]]] = defaultdict(dict)
    for (L, R), rs in grouped.items():
        ds = [float(r["d_final"]) for r in rs]
        masses = [
            float(r["e_per_debris_v"])
            for r in rs
            if float(r.get("d_final", 0.0)) > 0.0
            and r.get("e_per_debris_v") is not None
        ]
        out[L][R] = {
            "d": fmean(ds),
            "d_se": (pstdev(ds) / math.sqrt(len(ds))) if len(ds) > 1 else 0.0,
            "n": len(ds),
            "eD": fmean(masses) if masses else float("nan"),
        }
    return {L: dict(sorted(v.items())) for L, v in sorted(out.items())}


def loglog_interp(series: dict[int, dict[str, float]], R: float) -> tuple[float, str]:
    """Interpolate d(R) in log-log space between neighbouring positive points."""
    pts = sorted((float(r), float(v["d"])) for r, v in series.items() if v["d"] > 0)
    if len(pts) < 2:
        raise ValueError("need at least two positive KZ points")

    if R <= pts[0][0]:
        (r0, d0), (r1, d1) = pts[0], pts[1]
        flag = "extrapolated-low"
    elif R >= pts[-1][0]:
        (r0, d0), (r1, d1) = pts[-2], pts[-1]
        flag = "extrapolated-high"
    else:
        flag = "interpolated"
        for left, right in zip(pts, pts[1:]):
            if left[0] <= R <= right[0]:
                (r0, d0), (r1, d1) = left, right
                break
        else:  # pragma: no cover
            raise AssertionError("bracket not found")

    slope = math.log(d1 / d0) / math.log(r1 / r0)
    d = d0 * (R / r0) ** slope
    return max(0.0, min(1.0, d)), flag


def logR_interp_field(series: dict[int, dict[str, float]], R: float, field: str) -> float:
    """Interpolate a positive diagnostic linearly in log R."""
    pts = sorted(
        (float(r), float(v[field]))
        for r, v in series.items()
        if field in v and math.isfinite(float(v[field]))
    )
    if not pts:
        return float("nan")
    if len(pts) == 1:
        return pts[0][1]
    if R <= pts[0][0]:
        (r0, y0), (r1, y1) = pts[0], pts[1]
    elif R >= pts[-1][0]:
        (r0, y0), (r1, y1) = pts[-2], pts[-1]
    else:
        for left, right in zip(pts, pts[1:]):
            if left[0] <= R <= right[0]:
                (r0, y0), (r1, y1) = left, right
                break
        else:  # pragma: no cover
            raise AssertionError("bracket not found")
    t = math.log(R / r0) / math.log(r1 / r0)
    return y0 + t * (y1 - y0)


def percolation_status(d: float) -> str:
    crystal = 1.0 - d
    if crystal < 0.18:
        return "FAIL: crystal backbone below measured span floor"
    if crystal < 0.40:
        return "MARGINAL: measured percolation crossover"
    return "PASS: spanning crystal backbone in L=8 audit"


def live_window_slope(series: dict[int, dict[str, float]]) -> float | None:
    pts = [
        (math.log(float(R)), math.log(v["d"]))
        for R, v in series.items()
        if 0.02 < v["d"] < 0.6
    ]
    if len(pts) < 3:
        return None
    sx = sum(x for x, _ in pts)
    sy = sum(y for _, y in pts)
    sxx = sum(x * x for x, _ in pts)
    sxy = sum(x * y for x, y in pts)
    n = len(pts)
    return (n * sxy - sx * sy) / (n * sxx - sx * sx)


def energy_fraction(d: float, chi_visible: float) -> float:
    """f_D for chi_visible = e_visible_per_good_vertex / e_debris_per_bad_vertex."""
    return d / (d + (1.0 - d) * chi_visible)


def required_chi_for_target_fraction(d: float, target_fraction: float) -> float:
    """Solve target f_D = d / (d + (1-d) chi) for chi."""
    if not 0.0 < d < 1.0:
        return float("nan")
    f = target_fraction
    return d * (1.0 - f) / (f * (1.0 - d))


def main() -> int:
    rows = load_rows(DATA)
    means = grouped_means(rows)

    all_masses = [
        float(r["e_per_debris_v"])
        for r in rows
        if float(r.get("d_final", 0.0)) > 0.0
        and r.get("e_per_debris_v") is not None
    ]
    eD_mean = fmean(all_masses)

    branches = [
        (
            "canon-geometric-gamma-0.995",
            math.log(T_END / T_START) / math.log(GAMMA_CANON),
            "pinned by crystallisation.tex: T <- gamma T per sweep",
        ),
        (
            "expansion-radial-28tick",
            28.0 * DELTA_N,
            "conditional: T proportional a^-1 and one 28-clock tick block per e-fold",
        ),
        (
            "expansion-bipartite-56tick",
            56.0 * DELTA_N,
            "conditional: same, charging the two-step bipartite physical cycle",
        ),
    ]

    print("[0] K04 cooling-driver audit")
    print(f"    rows loaded: {len(rows)} from {DATA.name}")
    print(f"    ramp span: T={T_START:g} -> {T_END:g}; Delta ln a = ln(12) = {DELTA_N:.6f}")
    print(f"    debris mass ledger: raw <e_D> = {eD_mean:.3f} w6/debris-vertex")
    print()

    print("[1] finite-size KZ slopes already imply a caveat")
    for L, series in means.items():
        slope = live_window_slope(series)
        if slope is None:
            print(f"    L={L}: insufficient live-window points for a slope")
        else:
            print(f"    L={L}: d_trapped ~ R^({slope:.3f}) in 0.02 < d < 0.6")
    print("    caveat: these are toy finite-size slopes, not a converged KZ exponent.")
    print()

    print("[2] cooling branches, interpolated on the measured KZ curves")
    branch_results: list[tuple[str, float, list[tuple[float, float, float]]]] = []
    for name, R, note in branches:
        print(f"    {name}: R = {R:.1f} sweeps")
        print(f"      note: {note}")
        vals = []
        for L, series in means.items():
            d, flag = loglog_interp(series, R)
            eD = logR_interp_field(series, R, "eD")
            rhoD = d * eD
            vals.append((d, eD, rhoD))
            crystal = 1.0 - d
            print(
                f"      L={L}: d={d:.3f}, crystal={crystal:.3f}, "
                f"e_D={eD:.3f} w6, d*e_D={rhoD:.3f} w6/site, "
                f"{percolation_status(d)} ({flag})"
            )
        ds = [d for d, _, _ in vals]
        rhos = [rhoD for _, _, rhoD in vals]
        print(
            f"      finite-size envelope: d in [{min(ds):.3f}, {max(ds):.3f}], "
            f"d*e_D in [{min(rhos):.3f}, {max(rhos):.3f}] w6/site"
        )
        branch_results.append((name, R, vals))
    print()

    print("[3] energy-fraction bookkeeping shared by debris and later budget reconciliation")
    print("    define e_D = debris excess energy per debris vertex in w6 units")
    print("    define e_V = visible/register energy per crystallized vertex in the same units")
    print("    then rho_D = n_v d e_D w6 and rho_V = n_v (1-d) e_V w6")
    print("    so f_D = d / [d + (1-d) chi], chi = e_V/e_D")
    print(
        f"    reference target: f_D={REFERENCE_MATTER_DARK_FRACTION:.2f} "
        "for the observed matter dark fraction"
    )
    print("    open bridges: chi from the register ledger, n_v/V, and w6 -> physical scale.")
    for name, _, vals in branch_results:
        d_mid = fmean(d for d, _, _ in vals)
        rho_mid = fmean(rhoD for _, _, rhoD in vals)
        eD_mid = fmean(eD for _, eD, _ in vals)
        f_chi_1 = energy_fraction(d_mid, 1.0)
        f_chi_3 = energy_fraction(d_mid, 3.0)
        chi_for_half = required_chi_for_target_fraction(d_mid, 0.5)
        chi_for_target = required_chi_for_target_fraction(
            d_mid, REFERENCE_MATTER_DARK_FRACTION
        )
        eV_for_target = chi_for_target * eD_mid
        print(
            f"    {name}: mean d={d_mid:.3f}, mean e_D={eD_mid:.3f} w6, "
            f"mean d*e_D={rho_mid:.3f} w6/site; "
            f"f_D(chi=1)={f_chi_1:.3f}, f_D(chi=3)={f_chi_3:.3f}, "
            f"chi needed for f_D=1/2 is {chi_for_half:.3f}, "
            f"chi needed for f_D=0.84 is {chi_for_target:.3f} "
            f"(e_V={eV_for_target:.3f} w6/good-vertex)"
        )
    print()

    assert 490.0 < branches[0][1] < 500.0
    assert 69.0 < branches[1][1] < 70.0
    assert 139.0 < branches[2][1] < 140.0
    assert 2.0 < eD_mean < 2.3
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
