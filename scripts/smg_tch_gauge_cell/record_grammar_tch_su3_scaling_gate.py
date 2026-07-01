#!/usr/bin/env python3
r"""Scaling gate for the full-SU(3) bond-bipyramid scale-setting runs.

This is a post-processor for ``record_grammar_tch_su3_bulk_scale_setting.py``.
It asks the one question that separates a strong-coupling structural result
from a continuum QCD-scale result:

    does T_c / sqrt(sigma) stabilise as Nt and spatial volume increase?

The script deliberately does not repair or reinterpret the Monte Carlo output.
It reports three separate diagnostics:

1. Nt stability of the largest-volume ``T_c/sqrt(sigma)`` summary.
2. finite-volume drift inside each Nt row.
3. sigma-estimator quality, especially whether the chosen value is still only
   a 1x1 Creutz estimator and whether it is sampled near the FSS beta_c.

Verdict vocabulary
------------------
``SCALE-SET-CANDIDATE``
    The measured Nt spread is below the requested tolerance and no hard warning
    is present.  This is still only a candidate until production statistics and
    continuum extrapolation are documented.

``NOT-SCALE-SET``
    The measured Nt spread is too large, or the sigma estimates/volume drift are
    too unstable to promote.

This script is intentionally conservative.  A positive result should survive
this gate; a negative result should not be softened by prose.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class NtPoint:
    source: str
    profile: str
    nt: int
    nvol: int
    largest_size: tuple[int, int, int]
    beta_c_fss: float
    beta_sigma: float
    a_sqrt_sigma: float | None
    tc_over_sqrt_sigma: float | None
    slope_1_over_nv: float | None
    volume_tc_values: tuple[float, ...]
    volume_beta_c_values: tuple[float, ...]
    sigma_labels_at_beta: tuple[str, ...]


def finite(value: object) -> bool:
    return isinstance(value, (int, float)) and math.isfinite(float(value))


def rel_spread(values: Iterable[float]) -> float | None:
    vals = [float(v) for v in values if math.isfinite(float(v))]
    if len(vals) < 2:
        return None
    lo, hi = min(vals), max(vals)
    mid = 0.5 * (lo + hi)
    if mid == 0:
        return None
    return (hi - lo) / mid


def adjacent_spreads(points: list[NtPoint]) -> list[tuple[int, int, float]]:
    out: list[tuple[int, int, float]] = []
    ordered = [p for p in sorted(points, key=lambda p: p.nt) if p.tc_over_sqrt_sigma]
    for left, right in zip(ordered, ordered[1:]):
        a = float(left.tc_over_sqrt_sigma)
        b = float(right.tc_over_sqrt_sigma)
        mid = 0.5 * (a + b)
        if mid:
            out.append((left.nt, right.nt, abs(a - b) / mid))
    return out


def linear_extrapolate(points: list[NtPoint]) -> tuple[float, float] | None:
    """Fit y = intercept + slope / Nt^2 using unweighted least squares."""

    xs: list[float] = []
    ys: list[float] = []
    for point in points:
        if point.tc_over_sqrt_sigma and point.nt > 0:
            xs.append(1.0 / (point.nt * point.nt))
            ys.append(float(point.tc_over_sqrt_sigma))
    if len(xs) < 3:
        return None
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if denom == 0:
        return None
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return intercept, slope


def load_points(paths: list[str]) -> list[NtPoint]:
    expanded: list[Path] = []
    for item in paths:
        matches = [Path(p) for p in glob.glob(item)]
        expanded.extend(matches if matches else [Path(item)])
    points: list[NtPoint] = []
    seen: set[tuple[str, int]] = set()
    for path in sorted(expanded):
        if not path.exists() or path.suffix != ".json":
            continue
        data = json.loads(path.read_text())
        profile = str(data.get("profile", "unknown"))
        runs_by_nt = data.get("runs", {})
        for summary in data.get("summaries", []):
            nt = int(summary["nt"])
            key = (str(path), nt)
            if key in seen:
                continue
            seen.add(key)
            runs = runs_by_nt.get(str(nt), [])
            volume_tc_values = tuple(
                float(run["tc_over_sqrt_sigma_near"])
                for run in runs
                if finite(run.get("tc_over_sqrt_sigma_near"))
            )
            volume_beta_c_values = tuple(
                float(run["beta_c_chi"])
                for run in runs
                if finite(run.get("beta_c_chi"))
            )
            labels: list[str] = []
            beta_sigma = float(summary["beta_sigma"]) if finite(summary.get("beta_sigma")) else float("nan")
            for run in runs:
                for row in run.get("rows", []):
                    if finite(row.get("beta")) and abs(float(row["beta"]) - beta_sigma) < 1e-9:
                        label = row.get("creutz_label")
                        if label:
                            labels.append(str(label))
            points.append(
                NtPoint(
                    source=str(path),
                    profile=profile,
                    nt=nt,
                    nvol=len(runs),
                    largest_size=tuple(int(x) for x in summary["largest_size"]),
                    beta_c_fss=float(summary["beta_c_fss"]),
                    beta_sigma=beta_sigma,
                    a_sqrt_sigma=(
                        float(summary["a_sqrt_sigma"])
                        if finite(summary.get("a_sqrt_sigma"))
                        else None
                    ),
                    tc_over_sqrt_sigma=(
                        float(summary["tc_over_sqrt_sigma"])
                        if finite(summary.get("tc_over_sqrt_sigma"))
                        else None
                    ),
                    slope_1_over_nv=(
                        float(summary["slope_1_over_nv"])
                        if finite(summary.get("slope_1_over_nv"))
                        else None
                    ),
                    volume_tc_values=volume_tc_values,
                    volume_beta_c_values=volume_beta_c_values,
                    sigma_labels_at_beta=tuple(labels),
                )
            )
    return points


def profile_groups(points: list[NtPoint]) -> dict[str, list[NtPoint]]:
    groups: dict[str, list[NtPoint]] = {}
    for point in points:
        groups.setdefault(point.profile, []).append(point)
    return {key: sorted(value, key=lambda item: item.nt) for key, value in groups.items()}


def fmt(value: float | None, digits: int = 4) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{digits}f}"


def verdict_for(points: list[NtPoint], nt_tol: float, volume_tol: float, beta_tol: float) -> tuple[str, list[str]]:
    warnings: list[str] = []
    tc_values = [p.tc_over_sqrt_sigma for p in points if p.tc_over_sqrt_sigma]
    nt_spread = rel_spread(float(v) for v in tc_values if v is not None)
    if len(tc_values) < 3:
        warnings.append("fewer than three Nt values with sigma")
    elif nt_spread is None or nt_spread > nt_tol:
        warnings.append(f"Nt spread {fmt(nt_spread, 3)} exceeds tolerance {nt_tol:.3f}")
    for point in points:
        v_spread = rel_spread(point.volume_tc_values)
        if v_spread is not None and v_spread > volume_tol:
            warnings.append(f"Nt={point.nt} volume Tc/sqrt spread {v_spread:.3f} exceeds {volume_tol:.3f}")
        if point.a_sqrt_sigma is None:
            warnings.append(f"Nt={point.nt} has no usable sigma estimate")
        if point.sigma_labels_at_beta and all(label == "1x1" for label in point.sigma_labels_at_beta):
            warnings.append(f"Nt={point.nt} sigma is still only 1x1 Creutz at beta_sigma")
        if math.isfinite(point.beta_c_fss) and math.isfinite(point.beta_sigma):
            mismatch = abs(point.beta_sigma - point.beta_c_fss)
            if mismatch > beta_tol:
                warnings.append(f"Nt={point.nt} |beta_sigma-beta_c_FSS|={mismatch:.3f} exceeds {beta_tol:.3f}")
    return ("NOT-SCALE-SET" if warnings else "SCALE-SET-CANDIDATE"), warnings


def print_group(profile: str, points: list[NtPoint], nt_tol: float, volume_tol: float, beta_tol: float) -> dict:
    print(f"\n[{profile}]")
    print("  Nt  nvol largest    beta_c_FSS beta_sigma  a√sigma  Tc/√sigma  vol_spread  labels")
    for point in points:
        labels = ",".join(point.sigma_labels_at_beta) or "-"
        v_spread = rel_spread(point.volume_tc_values)
        print(
            f"  {point.nt:2d}  {point.nvol:4d} "
            f"{'x'.join(str(x) for x in point.largest_size):>9s} "
            f"{fmt(point.beta_c_fss):>10s} {fmt(point.beta_sigma):>10s} "
            f"{fmt(point.a_sqrt_sigma):>8s} {fmt(point.tc_over_sqrt_sigma):>10s} "
            f"{fmt(v_spread, 3):>10s}  {labels}"
        )

    tc_values = [p.tc_over_sqrt_sigma for p in points if p.tc_over_sqrt_sigma]
    nt_spread = rel_spread(float(v) for v in tc_values if v is not None)
    print(f"  Nt spread in T_c/sqrt(sigma): {fmt(nt_spread, 3)}")
    adj = adjacent_spreads(points)
    if adj:
        print("  adjacent Nt jumps:", ", ".join(f"{a}->{b}: {jump:.3f}" for a, b, jump in adj))
    extrap = linear_extrapolate(points)
    if extrap:
        intercept, slope = extrap
        print(f"  naive y=c+s/Nt^2 fit: c={intercept:.4f}, s={slope:.4f}")
    verdict, warnings = verdict_for(points, nt_tol, volume_tol, beta_tol)
    print(f"  verdict: {verdict}")
    for warning in warnings:
        print(f"    WARN: {warning}")
    return {
        "profile": profile,
        "nt_spread": nt_spread,
        "adjacent_spreads": adj,
        "linear_1_over_nt2": extrap,
        "verdict": verdict,
        "warnings": warnings,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="JSON result paths or glob patterns")
    parser.add_argument("--nt-tol", type=float, default=0.15, help="relative Nt spread tolerance")
    parser.add_argument("--volume-tol", type=float, default=0.25, help="relative finite-volume spread tolerance")
    parser.add_argument("--beta-tol", type=float, default=0.50, help="allowed |beta_sigma-beta_c_FSS|")
    parser.add_argument("--json", default=None, help="optional machine-readable summary path")
    args = parser.parse_args()

    points = load_points(args.paths)
    if not points:
        raise SystemExit("no completed JSON scale-setting result files found")

    print("Full-SU(3) bond-bipyramid scaling gate")
    print("=" * 72)
    print(f"loaded {len(points)} Nt point(s) from {len(set(p.source for p in points))} file(s)")
    reports = [
        print_group(profile, group, args.nt_tol, args.volume_tol, args.beta_tol)
        for profile, group in profile_groups(points).items()
    ]

    if args.json:
        payload = {
            "points": [point.__dict__ for point in points],
            "reports": reports,
            "nt_tol": args.nt_tol,
            "volume_tol": args.volume_tol,
            "beta_tol": args.beta_tol,
        }
        Path(args.json).write_text(json.dumps(payload, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
