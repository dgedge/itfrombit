#!/usr/bin/env python3
r"""K04 glass deep worker: one (mode, L, rep) job -> one JSON line to stdout.

Reuses the TESTED observable functions in k04_constraint_glass_observables.py, one rep per job, so
the deep box can fan thousands of independent reps across cores (xargs -P) for clean statistics --
the scaling step the laptop pilot named as the publishable next move.

  mode T : one trajectory (total activity K_t + persistence trace P(t)) + one STRICTLY-DOWNHILL
           descent to an inherent structure / local energy minimum (the well-defined jammed core).
           [dE<0 only: neutral dE=0 shuffles otherwise let the closure wander forever -- the laptop's
           dE<=0 closure at cap=12 was a small-cap artifact, not a real jam. Strictly downhill
           terminates by construction (energy bounded below) and is the honest object.]
  mode R : coupled-replica overlap Q across a finer eps grid for this rep (the RFOT gate).

Aggregation (k04_glass_deep_analysis.py): chi4(t)=N*Var_reps[P], its peak vs L (cooperative length);
activity reweighting -> psi(s); jam fraction / moves / residual vs L; Q(eps) vs L.

Usage:  python k04_glass_deep_worker.py <T|R> <L> <rep>
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import k04_constraint_glass_observables as obs  # noqa: E402

MODE, L, REP = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

# deeper than the laptop pilot (sweeps 16): longer time, hotter prep, large jam cap, finer eps
T = 0.7
SWEEPS = 24
RANDOM_SWEEPS = 12
PREP_T = 0.7
PREP_SWEEPS = 10
SAMPLE_EVERY = 2
EPS = [0.0, 0.25, 0.5, 1.0, 2.0]


def descend_to_inherent(L, rep):
    """Strictly-downhill (dE<0) greedy descent from the prepared state to a local energy minimum.
    Terminates by construction (energy bounded below; each move lowers it by a discrete quantum) ->
    a well-defined inherent structure / jammed core. Returns (descent length, inherent defect,
    residual downhill moves at the min) -- residual is 0 by construction (sanity check)."""
    k04 = obs.k04
    B = obs.prepare_state(L, seed=72_001 + 101 * rep + L, random_sweeps=RANDOM_SWEEPS,
                          prep_T=PREP_T, prep_sweeps=PREP_SWEEPS)
    moves = 0
    while True:
        chosen = next((mv for mv in k04.flippable_moves(B, L) if k04.delta_energy(B, mv, L) < 0.0), None)
        if chosen is None:
            break
        obs.apply_move(B, chosen)
        moves += 1
    resid = sum(1 for mv in k04.flippable_moves(B, L) if k04.delta_energy(B, mv, L) < 0.0)
    return moves, k04.defect_fraction(B, L), resid


t0 = time.time()
if MODE == "T":
    tr = obs.run_trajectory(L=L, T=T, sweeps=SWEEPS, rep=REP, random_sweeps=RANDOM_SWEEPS,
                            prep_T=PREP_T, prep_sweeps=PREP_SWEEPS, sample_every=SAMPLE_EVERY)
    dmoves, idef, dresid = descend_to_inherent(L, REP)
    out = dict(mode="T", L=L, rep=REP, N=L**3, activity=tr.activity, attempted=tr.attempted,
               sweeps=SWEEPS, final_defect=round(tr.final_defect, 6),
               samples=[{"t": round(s["sweep"], 3), "P": round(s["persistence"], 6)} for s in tr.samples],
               descent_moves=dmoves, inherent_def=round(idef, 6), descent_resid=dresid,
               secs=round(time.time() - t0, 1))
elif MODE == "R":
    overlap = {}
    for eps in EPS:
        r = obs.run_coupled_replica(L=L, T=T, epsilon=eps, sweeps=SWEEPS, rep=REP,
                                    random_sweeps=RANDOM_SWEEPS, prep_T=PREP_T, prep_sweeps=PREP_SWEEPS)
        overlap[str(eps)] = round(r.overlap, 6)
    out = dict(mode="R", L=L, rep=REP, N=L**3, sweeps=SWEEPS, overlap=overlap,
               secs=round(time.time() - t0, 1))
else:
    raise SystemExit(f"unknown mode {MODE!r} (expected T or R)")

print(json.dumps(out))
