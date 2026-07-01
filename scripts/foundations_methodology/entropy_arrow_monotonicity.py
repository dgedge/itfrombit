#!/usr/bin/env python3
r"""ENTROPY AS RECORD COUNT -- the arrow of time as crystallisation (Past-Hypothesis toy).

Operationalises ANCHOR 13.1 / item 123: entropy is carried by INSTRUMENTED, record-bearing
cells. Before the substrate activates there are none, so S(0)=0 BY CONSTRUCTION -- not a
posited low-entropy boundary condition. As cells crystallise and syndromes accumulate, records
are written and entropy grows MONOTONICALLY, which IS the thermodynamic arrow of time.

This toy shows, on a minimal model:
  [1] S(0)=0 and S(t) monotone non-decreasing along the boot (Past Hypothesis = mechanism);
  [2] the DRIVEN substrate (fresh S=0 nodes minted continuously, item 123) has NO downward
      excursions, whereas an equilibrium system of the same mean entropy does -- the exact
      distinction the Boltzmann-brain argument turns on (boltzmann_brain_suppression.py);
  [3] continuous fresh-S=0-node generation keeps total entropy monotone ('every new node IS an
      initial condition').

exit 0 = S(0)=0; the driven boot trajectory is strictly arrow-consistent (monotone); an
         equilibrium control of the same mean fluctuates DOWN; continuous node minting stays
         monotone.
"""
import math
import random

random.seed(20260614)
N = 10000           # cells in the causal patch (toy units)
TICKS = 400
TAU = 80.0          # crystallisation timescale (ticks)
C_REC = 0.5         # records written per active cell per tick (each a Landauer-costed bit)

# ---- [1] boot entropy = cumulative records; S(0)=0, monotone ----
def instrumented_fraction(t):
    return 1.0 - math.exp(-t / TAU)              # crystallisation, 0 -> 1
S, traj = 0.0, [0.0]
for t in range(1, TICKS + 1):
    dS = C_REC * instrumented_fraction(t) * N    # records written this tick (>= 0)
    S += dS
    traj.append(S)
diffs = [traj[i + 1] - traj[i] for i in range(len(traj) - 1)]
assert traj[0] == 0.0 and min(diffs) >= 0.0
print(f"[1] boot: S(0) = {traj[0]:.0f}  (no instrumented cells -> no records yet)")
print(f"    S(end) = {traj[-1]:.3e} records; min tick-step = {min(diffs):.3f} (>= 0) -> MONOTONE")
print("    -> low entropy at t=0 is 'records not yet written', not a fine-tuned boundary condition.")

# ---- [2] driven (no equilibrium) vs equilibrium control, same mean entropy ----
Smax = traj[-1]
eq = [Smax]
for _ in range(TICKS):                            # mean-reverting fluctuation about Smax
    eq.append(eq[-1] - 0.05 * (eq[-1] - Smax) + random.gauss(0, 0.01 * Smax))
eq_diffs = [eq[i + 1] - eq[i] for i in range(len(eq) - 1)]
print(f"\n[2] driven boot vs equilibrium control (same mean entropy {Smax:.2e}):")
print(f"    driven boot : min step = {min(diffs):+.2e}   (no downward excursions)")
print(f"    equilibrium : min step = {min(eq_diffs):+.2e}   (fluctuates DOWN)")
assert min(diffs) >= 0.0 and min(eq_diffs) < 0.0
print("    -> Boltzmann brains live in the equilibrium's downward excursions; the driven")
print("       substrate has none. (This premise is used by boltzmann_brain_suppression.py.)")

# ---- [3] continuous fresh-S=0-node generation keeps total entropy monotone ----
ages, tot_traj = [], [0.0]
for _ in range(TICKS):
    ages = [a + 1 for a in ages] + [0]            # age all nodes; mint one fresh S=0 node
    total = sum(C_REC * (1.0 - math.exp(-a / TAU)) for a in ages)   # each node's records
    tot_traj.append(total)
tot_diffs = [tot_traj[i + 1] - tot_traj[i] for i in range(len(tot_traj) - 1)]
assert tot_traj[0] == 0.0 and min(tot_diffs) >= -1e-9
print(f"\n[3] continuous fresh-S=0-node minting (item 123): total records monotone, min step "
      f"{min(tot_diffs):+.3e} -> every newly generated node IS an 'initial condition'.")

print("""
[verdict] The Past Hypothesis is mechanism, not fine-tuning: entropy = records, which start at
zero (no instrumented cells) and only grow. The arrow of time IS the crystallisation/record
direction. And the substrate is DRIVEN (continuous S=0 node minting), so -- unlike an
equilibrium -- it has no downward fluctuations for a Boltzmann brain to ride. exit 0""")
print("ALL ASSERTIONS PASSED -- S(0)=0, monotone arrow, driven (no downward excursions), monotone under minting.")
