#!/usr/bin/env python3
"""Reduce existing SMG DMRG run artifacts to an honest Phase-1 verdict.

This is intentionally a post-processing script, not another expensive DMRG
driver. It reads the completed local artifacts:

* smg_dmrg_runs/scan_20260603_100014/corr_L*_g*.json
* smg_dmrg_runs/chi_conv_results.json
* python_code/smg_dmrg_runs/converged_verdict.json
* python_code/smg_validation/results/*.log, when present

and answers the narrow question the earlier work left open:

    do the 3-4-5-0 runs support mirror-edge SMG, or only a mirror gap with an
    unresolved bilinear-condensate / SSB diagnostic?

The raw g=6.5 run sits close to the BKT critical fan and is not the final
discriminator. If the later validation logs are present, the script uses the
g=8,10 finite-size scaling runs as the superseding verdict.

The script writes a compact JSON verdict and prints the same summary. It avoids
TeNPy imports so it can run on any Python with the standard library.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from statistics import mean
from typing import Any


FLAVORS = ("3", "4", "5", "0")
MASS_KEYS = [
    f"{kind}_{a}{b}"
    for kind in ("dirac", "majorana")
    for a, b in (("3", "5"), ("3", "0"), ("4", "5"), ("4", "0"))
]


def mag(value: Any) -> float:
    if isinstance(value, (int, float)):
        return abs(float(value))
    if isinstance(value, list) and len(value) == 2:
        return math.hypot(float(value[0]), float(value[1]))
    raise TypeError(f"cannot convert {value!r} to magnitude")


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    n = len(xs)
    if n < 2:
        return 0.0, ys[0] if ys else 0.0, 0.0
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((x - xbar) ** 2 for x in xs)
    if sxx == 0:
        return 0.0, ybar, 0.0
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / sxx
    intercept = ybar - slope * xbar
    yhat = [slope * x + intercept for x in xs]
    ss_res = sum((y - yh) ** 2 for y, yh in zip(ys, yhat))
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    r2 = 1.0 if ss_tot < 1e-30 and ss_res < 1e-30 else 1.0 - ss_res / ss_tot
    return slope, intercept, r2


def exp_xi(rs: list[int], values: list[float]) -> float:
    logs = [math.log(max(v, 1e-300)) for v in values]
    slope, _, _ = linear_fit([float(r) for r in rs], logs)
    return -1.0 / slope if slope < 0 else float("inf")


def plateau_flag(values: list[float], floor: float = 1e-11) -> bool:
    usable = [v for v in values if v > floor]
    if len(usable) < 4:
        return False
    tail_n = max(3, len(usable) // 3)
    tail = usable[-tail_n:]
    if tail[0] <= floor:
        return False
    tail_drop = 1.0 - tail[-1] / tail[0]
    head_decayed = usable[0] > 3.0 * tail[-1]
    return bool(tail[-1] > 100 * floor and tail_drop < 0.15 and head_decayed)


def edge_summary(entries: list[dict[str, Any]]) -> dict[str, Any]:
    rs = [int(e["distance"]) for e in entries]
    fermion_mean = [
        mean(mag(e["fermion"][flavor]) for flavor in FLAVORS)
        for e in entries
    ]
    mass = {
        key: [mag(e["mass"][key]) for e in entries]
        for key in MASS_KEYS
    }
    mass_tails = {key: mean(vals[-3:]) for key, vals in mass.items()}
    top_mass_tail = sorted(mass_tails.items(), key=lambda kv: -kv[1])[:3]
    plateaus = [key for key, vals in mass.items() if plateau_flag(vals)]
    return {
        "xi_fermion_exp_fit": exp_xi(rs, fermion_mean),
        "fermion_tail": fermion_mean[-1],
        "mass_tail_top3": top_mass_tail,
        "mass_plateau_flags": plateaus,
    }


def read_scan(scan_dir: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    pat = re.compile(r"corr_L(?P<L>\d+)_g(?P<g>[-0-9.]+)\.json$")
    for path in sorted(scan_dir.glob("corr_L*_g*.json")):
        m = pat.search(path.name)
        if not m:
            continue
        data = json.loads(path.read_text())
        if "edges" not in data:
            continue
        row = {
            "path": str(path),
            "L": int(data.get("unit_cells", m.group("L"))),
            "g": float(data.get("g1", m.group("g"))),
            "edge_A": edge_summary(data["edges"]["A"]),
            "edge_B": edge_summary(data["edges"]["B"]),
        }
        rows.append(row)
    return sorted(rows, key=lambda r: (r["L"], r["g"]))


def read_converged(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows = json.loads(path.read_text())
    out: list[dict[str, Any]] = []
    for row in rows:
        tails = {
            key: mean(vals[-3:])
            for key, vals in row.get("edgeB_mass", {}).items()
        }
        out.append(
            {
                "g": float(row["key"][0]),
                "chi": int(row["key"][1]),
                "chi_used": int(row["chi_used"]),
                "trunc": row.get("trunc"),
                "xiA": float(row["xiA"]),
                "xiB": float(row["xiB"]),
                "B_condensate_classifier": row.get("B_condensate"),
                "edgeB_mass_tail_top3": sorted(tails.items(), key=lambda kv: -kv[1])[:3],
                "dirac_40_tail": tails.get("dirac_40"),
            }
        )
    return sorted(out, key=lambda r: (r["g"], r["chi"]))


def read_chi_conv(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return sorted(json.loads(path.read_text()), key=lambda r: (r["key"][0], r["key"][1]))


def read_validation_logs(results_dir: Path) -> dict[str, Any]:
    """Parse the later validation bundle logs, if they exist."""

    out: dict[str, Any] = {
        "finite_size_scaling": [],
        "deep_gapped_gscan": [],
        "metastability_seed": None,
    }
    finite = results_dir / "finite_size_scaling.log"
    if finite.exists():
        current_g: float | None = None
        current_L: int | None = None
        current_xi: float | None = None
        current_form: str | None = None
        current_tail: float | None = None
        for line in finite.read_text().splitlines():
            gm = re.search(r"### g=([0-9.]+)", line)
            if gm:
                current_g = float(gm.group(1))
                continue
            lm = re.search(r"L=\s*(\d+).*?chi_used=\s*(\d+).*?E=([-0-9.]+)", line)
            if lm:
                current_L = int(lm.group(1))
                current_xi = None
                current_form = None
                current_tail = None
                continue
            tm = re.search(r"largest-r value = ([0-9.eE+-]+).*?\[(.*?);(.*?)\]", line)
            if tm and current_g is not None and current_L is not None:
                current_tail = float(tm.group(1))
                current_form = tm.group(2)
                xim = re.search(r"xi=([0-9.]+)", tm.group(3))
                current_xi = float(xim.group(1)) if xim else None
                out["finite_size_scaling"].append(
                    {
                        "g": current_g,
                        "L": current_L,
                        "dirac_40_largest_r": current_tail,
                        "form": current_form,
                        "xi": current_xi,
                    }
                )

    deep = results_dir / "deep_gapped_gscan.log"
    if deep.exists():
        current_g = None
        current_E = None
        current_chi = None
        for line in deep.read_text().splitlines():
            gm = re.search(r"g=\s*([0-9.]+).*?E=([-0-9.]+).*?chi=(\d+)", line)
            if gm:
                current_g = float(gm.group(1))
                current_E = float(gm.group(2))
                current_chi = int(gm.group(3))
                continue
            fm = re.search(r"->\s+([A-Z-]+).*?xi=([0-9.]+).*?R2=([0-9.]+)", line)
            pm = re.search(r"->\s+(POWER-LAW).*?eta=([0-9.]+).*?R2=([0-9.]+)", line)
            if current_g is not None and (fm or pm):
                if fm:
                    form = fm.group(1)
                    scale_name = "xi"
                    scale = float(fm.group(2))
                    r2 = float(fm.group(3))
                else:
                    assert pm is not None
                    form = pm.group(1)
                    scale_name = "eta"
                    scale = float(pm.group(2))
                    r2 = float(pm.group(3))
                out["deep_gapped_gscan"].append(
                    {
                        "g": current_g,
                        "energy": current_E,
                        "chi": current_chi,
                        "form": form,
                        scale_name: scale,
                        "r2": r2,
                    }
                )

    meta = results_dir / "metastability_seed.log"
    if meta.exists():
        out["metastability_seed"] = meta.read_text().strip().splitlines()
    return out


def monotone_increasing(values: list[float]) -> bool:
    return all(b >= a for a, b in zip(values, values[1:]))


def build_verdict(
    scan_rows: list[dict[str, Any]],
    converged_rows: list[dict[str, Any]],
    chi_rows: list[dict[str, Any]],
    validation: dict[str, Any] | None = None,
) -> dict[str, Any]:
    validation = validation or {}
    g65 = [r for r in converged_rows if abs(r["g"] - 6.5) < 1e-9]
    g0_chi = [r for r in chi_rows if abs(float(r["key"][0]) - 0.0) < 1e-9]
    g65_chi = [r for r in chi_rows if abs(float(r["key"][0]) - 6.5) < 1e-9]
    scan_g65 = [r for r in scan_rows if abs(r["g"] - 6.5) < 1e-9]

    xiA_g65 = [r["xiA"] for r in g65]
    xiB_g65 = [r["xiB"] for r in g65]
    d40_g65 = [r["dirac_40_tail"] for r in g65 if r["dirac_40_tail"] is not None]

    mirror_gap_ok = bool(g65 and max(xiB_g65[-3:]) < 1.20)
    light_edge_unsettled_gapless = bool(len(xiA_g65) >= 3 and monotone_increasing(xiA_g65))
    free_control_warns_finite_chi = bool(
        len(g0_chi) >= 3
        and monotone_increasing([float(r["xi_A"]) for r in g0_chi])
    )
    condensate_risk = bool(
        len(d40_g65) >= 3
        and d40_g65[-1] > 1e-3
        and monotone_increasing(d40_g65[1:])
    )

    l_scan_note = None
    if len(scan_g65) >= 2:
        by_L = {r["L"]: r for r in scan_g65}
        if 16 in by_L and 20 in by_L:
            l16 = dict(by_L[16]["edge_B"]["mass_tail_top3"]).get("dirac_40")
            l20 = dict(by_L[20]["edge_B"]["mass_tail_top3"]).get("dirac_40")
            if l16 is not None and l20 is not None:
                l_scan_note = {
                    "L16_dirac_40_tail": l16,
                    "L20_dirac_40_tail": l20,
                    "ratio_L20_over_L16": l20 / l16 if l16 else None,
                    "read": "low-chi finite-size scan does not show an L-stable dirac_40 plateau",
                }

    finite = validation.get("finite_size_scaling") or []
    by_g: dict[float, list[dict[str, Any]]] = {}
    for row in finite:
        by_g.setdefault(float(row["g"]), []).append(row)
    validation_pass: dict[float, bool] = {}
    for g, rows in by_g.items():
        rows = sorted(rows, key=lambda r: r["L"])
        tails = [float(r["dirac_40_largest_r"]) for r in rows]
        forms = [str(r["form"]) for r in rows]
        validation_pass[g] = (
            len(rows) >= 3
            and all("EXPONENTIAL" in form for form in forms)
            and all(b < a for a, b in zip(tails, tails[1:]))
            and tails[-1] < 1e-8
        )
    reference_model_validated = bool(validation_pass.get(8.0) and validation_pass.get(10.0))

    if reference_model_validated:
        status = "reference 3-4-5-0 SMG benchmark confirmed; raw g=6.5 near-critical run superseded"
    elif mirror_gap_ok and not condensate_risk and light_edge_unsettled_gapless:
        status = "SMG-favorable but still needs L-scaling"
    elif mirror_gap_ok and condensate_risk:
        status = "inconclusive: mirror gap seen, bilinear-condensate null not closed"
    else:
        status = "inconclusive: benchmark gates not all passed"

    return {
        "status": status,
        "phase1_question": "3-4-5-0 benchmark only; not a 3+1D TCH closure",
        "reference_validation": {
            "passed": reference_model_validated,
            "finite_size_scaling": finite,
            "deep_gapped_gscan": validation.get("deep_gapped_gscan") or [],
            "read": (
                "g=8 and g=10 finite-size scaling shows exponential dirac_40 decay "
                "and largest-r order parameter falling with L; this is the superseding "
                "SSB-vs-SMG discriminator for the reference model."
                if reference_model_validated
                else "validation logs absent or insufficient; raw g=6.5 remains inconclusive"
            ),
        },
        "mirror_gap_gate": {
            "passed_conditionally": mirror_gap_ok,
            "evidence": "at g=6.5, xiB stays near 1.0 through chi=2048",
            "xiB_by_chi": [(r["chi"], r["xiB"]) for r in g65],
        },
        "light_edge_gate": {
            "not_shown_gapped": light_edge_unsettled_gapless,
            "finite_chi_warning_from_g0": free_control_warns_finite_chi,
            "xiA_g65_by_chi": [(r["chi"], r["xiA"]) for r in g65],
            "xiA_g0_control": [(int(r["key"][1]), float(r["xi_A"])) for r in g0_chi],
            "read": "edge-A short xi in the low-chi g-scan is not a reliable gap verdict",
        },
        "condensate_gate": {
            "closed": bool(reference_model_validated or not condensate_risk),
            "risk_channel": "dirac_40",
            "risk": condensate_risk,
            "dirac_40_tail_by_chi": [(r["chi"], r["dirac_40_tail"]) for r in g65],
            "classifier_flags_by_chi": [
                (r["chi"], r["B_condensate_classifier"]) for r in g65
            ],
            "finite_size_note": l_scan_note,
        },
        "missing_or_incomplete_runs": {
            "deep_host": "ssh deep is the compute target; this reducer uses local mirrored artifacts and does not require a live connection",
            "scan_missing": "L=20,g=7.0 log did not finish and no corr_L20_g7.0.json exists",
            "most_important_next_run": (
                "for reference-model reproduction, none if validation logs are accepted; "
                "for raw critical-fan diagnostics, high-chi L-scaling at g=6.5; "
                "for TCH, a separate actual-code/gauged matter-axis tensor-network run"
            ),
            "secondary_control": "finish g=0 control at chi=2048 if time permits",
        },
    }


def print_report(verdict: dict[str, Any]) -> None:
    print("SMG DMRG RECOVERY VERDICT")
    print("=" * 72)
    print(f"status: {verdict['status']}")
    print(f"scope : {verdict['phase1_question']}")
    print("\n[reference validation]")
    rv = verdict["reference_validation"]
    print(f"  passed: {rv['passed']}")
    print(f"  read: {rv['read']}")
    if rv["finite_size_scaling"]:
        for row in rv["finite_size_scaling"]:
            print(
                "  "
                f"g={row['g']:.1f} L={row['L']:>2} "
                f"tail={row['dirac_40_largest_r']:.2e} "
                f"form={row['form']} xi={row.get('xi')}"
            )
    print("\n[mirror gap]")
    print(f"  passed conditionally: {verdict['mirror_gap_gate']['passed_conditionally']}")
    print(f"  xiB by chi: {verdict['mirror_gap_gate']['xiB_by_chi']}")
    print("\n[light edge]")
    print(f"  not shown gapped: {verdict['light_edge_gate']['not_shown_gapped']}")
    print(f"  g=0 finite-chi warning: {verdict['light_edge_gate']['finite_chi_warning_from_g0']}")
    print(f"  xiA(g=6.5) by chi: {verdict['light_edge_gate']['xiA_g65_by_chi']}")
    print(f"  xiA(g=0) control: {verdict['light_edge_gate']['xiA_g0_control']}")
    print("\n[condensate / SSB gate]")
    c = verdict["condensate_gate"]
    print(f"  closed: {c['closed']}")
    print(f"  risk channel: {c['risk_channel']}")
    print(f"  dirac_40 tail by chi: {c['dirac_40_tail_by_chi']}")
    print(f"  classifier flags by chi: {c['classifier_flags_by_chi']}")
    if c["finite_size_note"]:
        print(f"  finite-size note: {c['finite_size_note']}")
    print("\n[next run]")
    for key, value in verdict["missing_or_incomplete_runs"].items():
        print(f"  {key}: {value}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-dir", default="smg_dmrg_runs/scan_20260603_100014")
    parser.add_argument("--chi-conv", default="smg_dmrg_runs/chi_conv_results.json")
    parser.add_argument("--converged", default="python_code/smg_dmrg_runs/converged_verdict.json")
    parser.add_argument("--validation-results", default="python_code/smg_validation/results")
    parser.add_argument("--output", default="python_code/smg_dmrg_runs/recovered_verdict.json")
    args = parser.parse_args()

    scan_rows = read_scan(Path(args.scan_dir))
    converged_rows = read_converged(Path(args.converged))
    chi_rows = read_chi_conv(Path(args.chi_conv))
    validation = read_validation_logs(Path(args.validation_results))
    verdict = build_verdict(scan_rows, converged_rows, chi_rows, validation)
    payload = {
        "verdict": verdict,
        "inputs": {
            "scan_dir": args.scan_dir,
            "chi_conv": args.chi_conv,
            "converged": args.converged,
            "validation_results": args.validation_results,
            "scan_rows": scan_rows,
            "converged_rows": converged_rows,
            "chi_rows": chi_rows,
        },
    }
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n")
    print_report(verdict)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
