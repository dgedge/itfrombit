#!/usr/bin/env python3
r"""Aggregate the Phase-0 validation + Phase-1 ladder per-rung JSONs into the single
deliverable smg_escalation_phase01_results.json (rows + metadata block).

Pure aggregation/IO. Reads the JSONs the runners already wrote; computes no physics.
"""
from __future__ import annotations
import json
from pathlib import Path

RUN_DIR = Path("/home/dave/octahedrons/python_code/smg_dmrg_runs")
OUT = Path("/home/dave/octahedrons/python_code/smg_escalation_phase01_results.json")
C3 = 4.0 / 3.0
DSMG = 2.0

# Phase-1 per-rung files written by smg_phase1_ladder.py (the unified Delta_raw + full_gap(t) rows).
PHASE1_FILES = [
    "phase1_1plaq_minimal.json",
    "phase1_1plaq_extended.json",
    "phase1_2plaq_minimal.json",
    "phase1_2plaq_extended_anchor.json",   # N=2 ext, the anchor (may be absent if still running)
    "phase1_3plaq_minimal.json",
    "phase1_3plaq_extended.json",          # only if it was run within budget
]

# Phase-0 validation anchors (the runner's own charged-gap JSONs, t=0 only).
PHASE0_VALIDATION = {
    "ext_2plaq_106460": "smg_region2_2plaq_extended_magnetic_gap.json",
    "min_3plaq_141435": "smg_region2_3plaq_minimal_magnetic_gap.json",
}


def softening(row):
    fg = row.get("full_gap_by_t", {})
    g0 = fg.get("0.0", {}).get("full_gap")
    out = {}
    for tk, d in fg.items():
        if tk != "0.0" and g0:
            out[tk] = (g0 - d["full_gap"]) / g0
    return out


def main():
    rows = []
    present = []
    for fn in PHASE1_FILES:
        p = RUN_DIR / fn
        if not p.exists():
            continue
        present.append(fn)
        blob = json.loads(p.read_text())
        for r in blob["rows"]:
            fg = r.get("full_gap_by_t", {})
            rows.append({
                "N": r["n_plaq"],
                "rep_set": r["rep_set"],
                "beta": r["beta"],
                "t_values": sorted(float(k) for k in fg.keys()),
                "dim": r["dim"],
                "delta_raw_t0": r["delta_raw"],          # conserved-sector meson gap at t=0
                "delta_mir_t0": r["delta_mir"],          # delta_raw - C3/beta
                "C3_over_beta": r["C3_over_beta"],
                "pred_2Dsmg_plus_C3_over_beta": r["pred_2Dsmg_plus_C3"],
                "full_gap_t0": fg.get("0.0", {}).get("full_gap"),
                "full_gap_by_t": {k: v["full_gap"] for k, v in fg.items()},
                "full_gap_softening_vs_t0": softening(r),
            })

    # sort rows for readability: rep_set, beta, N
    rows.sort(key=lambda x: (x["rep_set"], x["beta"], x["N"]))

    validation = {}
    for key, fn in PHASE0_VALIDATION.items():
        p = RUN_DIR / fn
        if p.exists():
            blob = json.loads(p.read_text())
            cell = blob[0] if isinstance(blob, list) else blob
            b1 = next((rr for rr in cell["rows"] if abs(rr["beta"] - 1.0) < 1e-6), None)
            validation[key] = {
                "dim": cell["dim"],
                "source_json": fn,
                "beta1_charged_gap_delta_raw": (b1["charged_gap"] if b1 else None),
                "pred_2Dsmg_plus_C3_at_beta1": 2 * DSMG + C3 / 1.0,
                "delta_mir_at_beta1": (b1["charged_gap"] - C3 / 1.0) if b1 else None,
            }

    meta = {
        "what": "SMG spin-network escalation, Phase 0 (validation gate) + Phase 1 (volume ladder).",
        "observable": {
            "delta_raw": "E0(n_charged>=1) - E0(n_charged=0) at t=0 (conserved charged sector; the meson gap).",
            "delta_mir": "delta_raw - C3/beta, C3=4/3 (electric-string-subtracted mirror offset).",
            "full_gap_softening": ("With t>0 the hop breaks n_charged conservation, so the charged-vs-vacuum "
                                   "sector gap is not clean. The hopping leg therefore reports the FULL spectral "
                                   "gap e1-e0 on the SAME basis (identical quantity to tch_2plaq_extended_hop.py) "
                                   "and softening=(full_gap(t=0)-full_gap(t))/full_gap(t=0). NB full_gap(t=0) != "
                                   "delta_raw in general: the full gap's first excited state is the cheapest "
                                   "excitation of ANY sector (here a vacuum-sector magnetic mode ~5.03 at beta=1), "
                                   "not the charged-sector gap (~5.40 at beta=1)."),
        },
        "constants": {"DSMG": DSMG, "C3": C3},
        "code_reused": {
            "build_geometry/build_basis/make_sparse_wilson/diag_arrays/run_cell/make_robust_vbasis/load_machinery":
                "smg_region2_nplaq_magnetic_gap.py",
            "sparse_hop/hop_tail/hop_head (gauge-covariant neutral fermion hop, validated |dHOP|<1e-10 vs oracle)":
                "tch_2plaq_extended_hop.py (defs only; STEP1-3 driver not executed)",
            "K_corner/F_emb/vbasis/REPS_MIN/REPS_EXT/DSMG/geometry":
                "tch_finish_port_and_swap.py (loaded via load_machinery, GATE A/B run on import)",
            "C3 + Delta_mir decomposition": "smg_electric_subtracted_gap_observable_theorem.py",
        },
        "validation_gate": validation,
        "phase1_files_present": present,
        "commands": [
            "# Phase 0 validation:",
            "/home/dave/tenpy-env/bin/python smg_region2_nplaq_magnetic_gap.py --n 2 --reps extended --betas 0.25 0.5 1.0 2.0",
            "/home/dave/tenpy-env/bin/python smg_region2_nplaq_magnetic_gap.py --n 3 --reps minimal  --betas 0.5 1.0 2.0",
            "# Phase 1 ladder (each rung):",
            "/home/dave/tenpy-env/bin/python smg_phase1_ladder.py --n {1,2,3} --reps {minimal,extended} --betas 0.5 1.0 --ts 0.0 1.0",
        ],
        "notes": [
            "t>0 fixed hopping value t=1.0 taken from tch_2plaq_extended_hop.py's own t-scan (t in {0.0,0.2,1.0}); t=1.0 is its largest.",
            "Phase-0 ext-2plaq beta=1.0 delta_raw reproduces 2*DSMG+C3/beta=5.3333 to 0.76% (got 5.37402); minimal-3plaq reproduces stored JSON to all digits.",
            "Tiny N=1 charge sectors use dense eigvalsh (the runner's sector_ground eigsh crashes for <=7-state sectors: ARPACK k>=N-1).",
        ],
    }

    OUT.write_text(json.dumps({"metadata": meta, "rows": rows}, indent=2) + "\n")
    print(f"wrote {OUT} with {len(rows)} rows from {len(present)} rung files")
    print("rungs present:", present)
    print("\nDelta_mir(N) summary (t=0):")
    for rep in ("minimal", "extended"):
        for beta in (0.5, 1.0):
            sub = [r for r in rows if r["rep_set"] == rep and abs(r["beta"] - beta) < 1e-9]
            if sub:
                s = "  ".join(f"N={r['N']}:{r['delta_mir_t0']:.4f}(d{r['dim']})" for r in sub)
                print(f"  {rep:9s} beta={beta}: {s}")
    print("\nfull-gap softening at t=1.0 (volume trend):")
    for rep in ("minimal", "extended"):
        for beta in (0.5, 1.0):
            sub = [r for r in rows if r["rep_set"] == rep and abs(r["beta"] - beta) < 1e-9]
            if sub:
                s = "  ".join(f"N={r['N']}:{100*r['full_gap_softening_vs_t0'].get('1.0',0):.2f}%" for r in sub)
                print(f"  {rep:9s} beta={beta}: {s}")


if __name__ == "__main__":
    main()
