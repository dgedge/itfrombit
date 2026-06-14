#!/usr/bin/env python3
r"""EMBEDDED-BRACKET ANALYSIS — the five reading gates, in registered order.

GATE 1: locate the ordering window per (L, ratio): T_melt (cold branch leaves the
        crystal), T_freeze (hot branch reaches order), hysteresis interval; T_c-
        candidates = equilibrated-and-merged ordered points (same branch-merge
        instrument as the toy: |d - d_late| < 0.08 per branch, |d_cold - d_hot| < 0.10).
GATE 2: d_eq at T_c- (only where Gate 1 produces a merged ordered point).
GATE 3: mass ledger — require E_excess per debris vertex > 0 there.
GATE 4: durability ledger — heal_spec dominated by extensive trails, CHECKED FOR
        L-SCALING (max trail / bonds at comparable window position, L=4 vs 6).
        Barrier numbers from walked trails are UPPER BOUNDS (labelled; the minimal
        ridge search is a registered follow-on).
GATE 5: abundance discussion is EMITTED ONLY IF Gates 3 and 4 both pass; otherwise
        the script prints the failed/blocked gates and stops. Machine-readable gates,
        no interpretive slack. exit 0 = data integrity (not gate success)."""
import json, math
import numpy as np
from collections import defaultdict

PATH = __file__.replace("k04_embedded_analysis.py", "k04_embedded_results.jsonl")
rows = [json.loads(l) for l in open(PATH)]
print(f"[0] {len(rows)} embedded results loaded "
      f"({sum(1 for r in rows if 'heal_spec' in r)} carry healing spectra)")
G = defaultdict(list)
for r in rows:
    G[(r["L"], r["w4"], round(r["T"], 4), r["start"])].append(r)

def branch(L, w4, T, st):
    rs = G.get((L, w4, T, st), [])
    if not rs:
        return None
    d = float(np.mean([r["d"] for r in rs]))
    dl = float(np.mean([r["d_late"] for r in rs]))
    return d, dl, rs

verdicts = {}
for (L, w4) in sorted({(l, w) for (l, w, *_ ) in G}):
    Ts = sorted({t for (l, w, t, s) in G if l == L and w == w4})
    print(f"\n=== L = {L}, ratio = {w4} ===")
    print(f"    {'T':>5} | {'d_cold':>7} {'d_hot':>7} | eq? merged?")
    T_melt = T_freeze = None
    merged_ordered = []
    for T in Ts:
        bc, bh = branch(L, w4, T, "cold"), branch(L, w4, T, "hot")
        if bc is None or bh is None:
            continue
        eq = abs(bc[0] - bc[1]) < 0.08 and abs(bh[0] - bh[1]) < 0.08
        mg = abs(bc[0] - bh[0]) < 0.10
        print(f"    {T:5.2f} | {bc[0]:7.3f} {bh[0]:7.3f} | {'Y' if eq else 'n'}   {'Y' if mg else 'n'}")
        if T_melt is None and bc[0] > 0.10:
            T_melt = T
        if bh[0] < 0.50:
            T_freeze = T
        if eq and mg and 0.0 <= (bc[0] + bh[0]) / 2 < 0.5:
            merged_ordered.append((T, (bc[0] + bh[0]) / 2, bc[2] + bh[2]))
    # GATE 1
    print(f"    GATE 1: cold leaves crystal at T ~ {T_melt}; hot orders up to T ~ {T_freeze}; "
          f"hysteresis {'[' + str(T_freeze) + ', ' + str(T_melt) + ']' if T_melt and T_freeze and T_freeze < T_melt else 'narrow/none'}")
    if not merged_ordered:
        print(f"    GATE 1: no equilibrated merged ordered point in this coarse grid -> fine grid target")
        verdicts[(L, w4)] = dict(gate1=None)
        continue
    Tc, d_eq, pool = max(merged_ordered)
    # GATE 2 + 3
    epd = [r.get("e_per_debris_v", 0.0) for r in pool if r.get("d_final", 0) > 0]
    e_pos = bool(epd) and float(np.mean(epd)) > 0
    print(f"    GATE 2: T_c- = {Tc}, d_eq = {d_eq:.4f}  (n = {len(pool)} runs)")
    if epd:
        print(f"    GATE 3: E_excess/debris vertex = {np.mean(epd):.3f} w6 -> {'PASS' if e_pos else 'FAIL (negative)'}")
    else:
        print(f"    GATE 3: N/A — zero surviving debris at termination in every pooled run")
        print(f"            (the crystal healed completely: itself the finding, not a failure)")
    # GATE 4 (needs spectra at/near Tc)
    specs = [r["heal_spec"] for r in pool if r.get("heal_spec")]
    if specs:
        fracs = [s[0] / (1.5 * L ** 3) for s in specs]
        print(f"    GATE 4: max-trail / bonds = {np.mean(fracs):.3f} (n = {len(specs)} spectra)")
    else:
        print(f"    GATE 4: no spectra at this point (pre-patch rows) -> fine runs supply")
    verdicts[(L, w4)] = dict(gate1=Tc, d_eq=d_eq, e_pos=e_pos,
                             frac=float(np.mean(fracs)) if specs else None)

# GATE 4 L-scaling + GATE 5
print("\n[GATE 4 L-SCALING + GATE 5]")
for w4 in sorted({w for (_, w) in verdicts}):
    vs = {L: v for (L, w), v in verdicts.items() if w == w4 and v.get("gate1") is not None}
    if len(vs) < 2:
        print(f"  ratio {w4}: insufficient merged points across L -> abundance NOT discussed (gates blocked)")
        continue
    fr = {L: v["frac"] for L, v in vs.items() if v["frac"] is not None}
    extensive = len(fr) >= 2 and all(f > 0.10 for f in fr.values())
    mass_ok = all(v["e_pos"] for v in vs.values())
    if mass_ok and extensive:
        print(f"  ratio {w4}: GATES 3+4 PASS -> dark-defect abundance readout licensed:")
        for L, v in sorted(vs.items()):
            print(f"      L={L}: d_eq = {v['d_eq']:.4f}")
    else:
        why = [] if mass_ok else ["mass ledger"]
        if not extensive:
            why.append("durability ledger (spectra missing or non-extensive)")
        print(f"  ratio {w4}: abundance NOT discussed — blocked on: {', '.join(why)}")

print("\nexit 0 — gates evaluated; interpretation lives in the session record.")
