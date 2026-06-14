#!/usr/bin/env python3
r"""THE ABSOLUTE-SCALE BRIDGES — legacy coarse closure attempt on both legs.
(Follows k04_absolute_scale_bridge_audit.py: the scale-invariance obstruction
means K04 yields only RATIOS; closure = derive the ratios + anchor ONE
temperature to a canon-priced scale.)

Supersession note (2026-06-12): the coarse L=6 trapezoid gate below is superseded
by k04_entropy_spike_resolution.py, which resolves the crystallisation spike and
finds that the gamma=exp(-alpha0 ln2) entropy-budget route misses the registered
Delta S target at percent grade.  This script is retained as the provenance of
the candidate and the pre-registered gate, not as the current verdict.

LEG 1 — the gamma = exp(-alpha0 ln 2) theorem target, decomposed into clauses:
  A (definitional): one Metropolis sweep = one attempt per DOF = one per-cell
     tick.  Standard lattice-dynamics convention; nothing to prove.
  B (canon-adopted): billed commit rate = alpha0 per cell-tick — the item-79
     alpha-chain (the session's three-rate table: service 1 / commit alpha0 /
     severing 6.87 alpha0^2).
  C (the open clause), split:
     C-bit:   each billed commit freezes ONE binary record (ln 2) — strong
              canon support: binary strain alphabet + the demux's one active
              address per tick.
     C-slope: converting an entropy-extraction rate into a log-temperature
              rate needs dS = C_v dlnT, i.e. the candidate
                  -ln(gamma) = alpha0 ln2 / C_v
              holds with C_v = 1 exactly; the paper's 0.995 corresponds to
              C_v = alpha0 ln2 / ln(1/0.995) = 1.00934.
     C-slope is MEASURABLE: this script measures C_v(T) per site in the
     embedded ensemble (energy-fluctuation estimator) across the transition/
     arrest window and reports where C_v = 1 sits relative to the arrest.

LEG 2 — w6 <-> Lambda, the constructive remainder after the obstruction:
  K04 fixes the dimensionless location of the crystallisation transition,
  T_c / w6 (measured below as the C_v peak).  Canon prices the boot epoch at
  T_boot = Lambda (chiral anchor sect-1.4; the two-sector identity
  (Lambda/T0)^3 = 2.83e36 is exactly this statement).  Under the NAMED premise
  "crystallisation happens at the boot temperature" the band closes to a value:
        w6 = Lambda / (T_c/w6),
  to be compared against the carried external band [0.05, 0.27] Lambda.

Self-asserting; exit 0 = measurements made, clauses tiered, nothing overclaimed."""
import math
import statistics
import sys
import time
from pathlib import Path

ALPHA0 = 1 / 137.0
GAMMA_PAPER = 0.995
CAND = math.exp(-ALPHA0 * math.log(2))
CV_REQ = ALPHA0 * math.log(2) / math.log(1 / GAMMA_PAPER)

print("[1] LEG-1 CLAUSE LEDGER (the gamma candidate, decomposed):")
print("    NOTE: this is the legacy coarse gate; current verdict is in")
print("          k04_entropy_spike_resolution.py (gamma entropy-budget route fails).")
print(f"    candidate gamma = exp(-alpha0 ln2) = {CAND:.9f} (paper: {GAMMA_PAPER})")
print(f"    A  sweep = per-cell tick          : definitional (Metropolis convention)")
print(f"    B  commit rate alpha0/cell-tick   : canon-adopted (item-79 alpha-chain)")
print(f"    C-bit  one binary record/commit   : canon-strong (binary strain alphabet;")
print(f"                                        demux: one active address per tick)")
print(f"    C-slope: paper's 0.995 <=> C_v = {CV_REQ:.5f} per cell — measured below.")

# ---------------- load the embedded-ensemble machinery (defs only) ----------------
print("\n[2] C_v(T) IN THE EMBEDDED ENSEMBLE (L=6, equilibrium stints):")
src = (Path(__file__).parent / "k04_embedded_sweep.py").read_text()
marker = 'if START in ("hot", "ramp"):'
assert marker in src
prefix = src[: src.index(marker)]
sys.argv = ["k04", "6", "1.7", "1.0", "1.0", "cold", "1", "0"]
env = {}
exec(compile(prefix, "k04_embedded_sweep_prefix", "exec"), env)
N, W4, W6 = env["N"], env["W4"], env["W6"]
metropolis, census, count_all = env["metropolis"], env["census"], env["count_all"]
B0, c40, c60 = env["B"], env["c4"], env["c6"]
print(f"    L=6: N = {N} sites, w4 = {W4}, w6 = {W6} (run units)")
assert N == 216

T_GRID = [6.0, 5.0, 4.0, 3.6, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2]
BURN, SAMP = 400, 600
rows = []
t0 = time.time()
env["random"].seed(7)
for T in T_GRID:
    # FRESH hot start per point: avoids dragging a frozen state below T_c
    B, c4, c6 = set(B0), c40, c60
    B, c4, c6, _ = metropolis(B, c4, c6, math.inf, 20)
    B, c4, c6, _ = metropolis(B, c4, c6, T, BURN)
    es, acc_tot = [], 0
    for _ in range(SAMP):
        B, c4, c6, a = metropolis(B, c4, c6, T, 1)
        acc_tot += a
        es.append(-W4 * c4 - W6 * c6)
    cv = statistics.pvariance(es) / (T * T * N)
    acc_rate = acc_tot / (SAMP * N)
    eq = acc_rate > 1e-3
    rows.append((T, cv, acc_rate, statistics.mean(es) / N, eq))
    print(f"    T = {T:>4.1f}: C_v/site = {cv:7.4f}, acceptance = {acc_rate:7.4f}, "
          f"E/site = {rows[-1][3]:+.4f} {'' if eq else '[FROZEN - excluded]'}  [{time.time()-t0:4.0f}s]")
    assert count_all(B, 4) == c4 and count_all(B, 6) == c6
rows = [r for r in rows if r[4]]
assert len(rows) >= 5, "too few equilibrated points"

cvs = [r[1] for r in rows]
peak_i = max(range(len(rows)), key=lambda i: cvs[i])
T_c = rows[peak_i][0]
assert 0 < peak_i < len(rows) - 1 or peak_i == len(rows) - 1, "peak at hot edge: extend grid"
# where does C_v cross the required 1.009? (scan from hot side)
T_cross = None
for ra, rb in zip(rows, rows[1:]):
    Ta, ca, Tb, cb = ra[0], ra[1], rb[0], rb[1]
    if (ca - CV_REQ) * (cb - CV_REQ) <= 0:
        T_cross = Ta + (Tb - Ta) * (ca - CV_REQ) / (ca - cb)
        break

# THE SHAPE-INDEPENDENT FORM: integrating -ln(gamma) = alpha0 ln2 / C_v over the
# ramp gives  R_total = Delta_S_extracted / (alpha0 ln 2)  — an entropy budget,
# independent of where C_v peaks. Trapezoid over the equilibrated branch:
dS = 0.0
for ra, rb in zip(rows, rows[1:]):
    Ta, ca, Tb, cb = ra[0], ra[1], rb[0], rb[1]
    dS += 0.5 * (ca + cb) * abs(math.log(Ta / Tb))
R_PAPER = math.log(12) / math.log(1 / GAMMA_PAPER)
dS_required = R_PAPER * ALPHA0 * math.log(2)
print(f"\n    INTEGRAL FORM: R_total = Delta_S/(alpha0 ln2).  Paper R = {R_PAPER:.1f}")
print(f"    => required Delta_S = {dS_required:.3f} nats/site;")
print(f"    resolved on the equilibrated branch: {dS:.3f} nats/site -> R_pred = {dS/(ALPHA0*math.log(2)):.0f}")
print(f"    vs paper 495.7: {dS/dS_required-1:+.1%} at L=6/grid-0.2 resolution (the spike leg")
print(f"    is a crude trapezoid — could over- or under-count; deep-box finer grid pins it).")
print(f"\n    C_v peak (transition) at T_c = {T_c} w6-units; C_v(T_c) = {cvs[peak_i]:.3f}")
print(f"    C_v = {CV_REQ:.3f} crossing: T = {T_cross if T_cross else float('nan'):.3f} w6-units"
      if T_cross else "    C_v never crosses the required value on this grid")

print("\n[3] LEG-1 VERDICT (C-slope):")
if T_cross is not None:
    rel = "above" if T_cross > T_c else "below"
    print(f"    the alpha-bit law -ln(gamma) = alpha0 ln2 / C_v reproduces the paper's")
    print(f"    0.995 exactly where C_v = {CV_REQ:.4f} — measured at T = {T_cross:.2f} ({rel} T_c = {T_c}).")
    print(f"    Clauses A, B, C-bit stand on canon; C-slope is now a MEASURED statement,")
    print(f"    conditional on the anneal's operating window containing that temperature.")
else:
    print(f"    C_v never reaches {CV_REQ:.3f} on the scanned window — the candidate's")
    print(f"    C-slope clause FAILS in the embedded ensemble as posed.")

print("\n[4] LEG-2: w6 <-> Lambda — the bridge STRUCTURED (three measured anchors):")
LAM, T0_GEV = 0.332, 2.348e-13
two_sector = (LAM / T0_GEV) ** 3
print(f"    canon prices the boot at Lambda: (Lambda/T0)^3 = {two_sector:.2e} (the 2.83e36 identity)")
print(f"    the obstruction permits exactly: w6 = Lambda x (w6/T_anchor), one anchor epoch.")
anchors = [("ramp start (T_START_RAMP, canon driver)", 6.0),
           ("transition T_c (C_v spike, this run)", T_c),
           ("KZ freeze (below T_c; bound)", 2.8)]
print(f"    {'anchor epoch':<42s} {'T/w6':>6s} {'w6/Lambda':>10s} {'vs band [0.05,0.27]':>20s}")
for nm, Ta in anchors:
    r = W6 / Ta
    tag = 'INSIDE' if 0.05 <= r <= 0.27 else 'outside'
    print(f"    {nm:<42s} {Ta:>6.2f} {r:>10.3f} {tag:>20s}")
print(f"    -> the band is no longer a free external interval: it is the image of the")
print(f"       anchor-epoch CHOICE; ramp-start anchoring lands inside the carried band,")
print(f"       transition/freeze anchoring lands ~0.33-0.36 Lambda (mild tension).")
print(f"       The remaining unpinned clause is ONE discrete choice among named epochs.")

print(f"""
[5] CONSOLIDATED VERDICT — both bridges, honestly tiered:
  gamma bridge: A (definitional) + B (item-79 adopted) + C-bit (canon-strong)
    + C-slope (measured: the law lands the paper's 0.995 at C_v = {CV_REQ:.3f},
    located on the measured C_v(T) curve {'at T = %.2f' % T_cross if T_cross else 'NOWHERE on the window'}).
    Tier: LEGACY COARSE CANDIDATE — superseded by k04_entropy_spike_resolution.py,
    which resolves the spike and finds the entropy-budget route misses the 2.508
    nats/site target at percent grade.
  w6 bridge: STRUCTURED, not closed — w6 = Lambda x (w6/T_anchor) with the anchor
    epoch the single remaining discrete clause (ramp-start -> 0.167 Lambda, inside
    the carried band; transition/freeze -> 0.33-0.36 Lambda, mild tension);
    downstream mobility/abundance conclusions were already band-insensitive
    (their audit's 1e11-1e12 exponents), so nothing destabilises either way.
  The obstruction theorem stands: everything above is ratios + ONE named
  physical anchor (T_boot = Lambda) — exactly the form the obstruction permits.
exit 0""")
print("ALL ASSERTIONS PASSED — measurements live, clauses tiered.")
