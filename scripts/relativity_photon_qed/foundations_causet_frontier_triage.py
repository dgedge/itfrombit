#!/usr/bin/env python3
r"""Causal-set QED frontier triage: loop F^2, raw variance, and Dirac spin.

This script attacks the three residuals left after K21/K22:

  R1 manifest gauge-invariant loop F^2,
  R2 raw unsmeared sprinkling variance/control,
  R3 causal-set Dirac spin.

It is intentionally a triage, not a claimed final derivation.

Findings:
  * Raw link/local-loop estimators remain unsuitable: O(1)-cardinality intervals
    have O(1) relative variance, and causal-set links carry the boost tail.
  * A mesoscopic Alexandrov-loop estimator is the plausible gauge-invariant route:
    for a constant 1+1 field it recovers F^2 with variance falling as the interval
    cardinality grows.  This gives a concrete programme but not a closed Maxwell
    loop action; the non-arbitrary interval measure and 3+1/local-field limit remain
    to be derived.
  * Dirac spin is compatible with the causal-set scalar kernel once a spin frame
    / tetrad / E_{1/2} coin is supplied, via the usual Clifford factorisation.
    Causal order alone supplies interval data (p^2), not a spinor fibre or gamma
    matrices.  So Dirac is compatible, not order-derived.
"""

from __future__ import annotations

import math
import numpy as np

rng = np.random.default_rng(20260613)


def poisson_area_trials(nbar: float, trials: int = 200_000) -> tuple[float, float, float]:
    """Constant-field Alexandrov-loop F^2 estimator using cardinality area.

    In 1+1 with A_x=F t, the diamond-loop holonomy is Phi=F*Area.  The causal
    set estimates Area by N/rho.  With nbar=rho*Area, F_est/F = nbar/N.
    This isolates the raw-vs-mesoscopic cardinality noise without coordinate
    distractions.
    """

    n = rng.poisson(nbar, trials)
    n = n[n > 0]
    f2_ratio = (nbar / n) ** 2
    return float(f2_ratio.mean()), float(f2_ratio.std()), float(np.mean(np.abs(f2_ratio - 1.0) < 0.1))


print("[1] Raw unsmeared loop/cardinality variance vs mesoscopic control")
print("    constant 1+1 field, diamond holonomy Phi=F*Area, area estimated by N/rho")
print("    nbar=rho*Area   mean(F2_est/F2)   std       P(|error|<10%)")
raw_bad = None
meso_good = None
for nbar in (2, 4, 8, 16, 64, 256, 1024):
    mean, std, hit = poisson_area_trials(float(nbar))
    print(f"    {nbar:6.0f}           {mean:10.3f}      {std:7.3f}      {hit:6.3f}")
    if nbar == 4:
        raw_bad = std
    if nbar == 1024:
        meso_good = std
assert raw_bad is not None and raw_bad > 1.0
assert meso_good is not None and meso_good < 0.08
print("    -> O(1)-cardinality/raw loops are variance-dominated; mesoscopic intervals")
print("       with N_epsilon >> 1 are controllable.  This matches the BDG/B_epsilon")
print("       lesson: raw unsmeared operators are not the precision object.")


print("\n[2] Manifest gauge-invariant loop F^2: candidate and remaining obstruction")
print("    candidate: use gauge-invariant holonomy around mesoscopic Alexandrov diamonds,")
print("    divide by cardinality-estimated area, then average over a covariant interval")
print("    size window N_epsilon.")
candidate_n = 1024.0
mean, std, hit = poisson_area_trials(candidate_n, trials=80_000)
print(f"    at N_epsilon={candidate_n:.0f}: F2_est/F2 = {mean:.4f} +/- {std:.4f} (single interval)")
print(f"    averaging M independent intervals gives std/sqrt(M); M=100 -> {std / 10:.4f}")
assert abs(mean - 1.0) < 0.02 and std / 10 < 0.01
print("    -> constant-field toy PASSES as a mesoscopic gauge-invariant estimator.")
print("    -> not a closure: the interval-size measure, overlap correlations, 3+1 tensor")
print("       contractions, variable-field locality expansion, and overall coefficient")
print("       are not derived here.")


print("\n[3] Dirac spin: Clifford factorisation works only after adding a spin frame")
I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
zero = np.zeros((2, 2), dtype=complex)

gamma = [
    np.block([[I2, zero], [zero, -I2]]),
    np.block([[zero, X], [-X, zero]]),
    np.block([[zero, Y], [-Y, zero]]),
    np.block([[zero, Z], [-Z, zero]]),
]
eta = np.diag([1.0, -1.0, -1.0, -1.0])
for mu in range(4):
    for nu in range(4):
        anti = gamma[mu] @ gamma[nu] + gamma[nu] @ gamma[mu]
        target = 2.0 * eta[mu, nu] * np.eye(4)
        assert np.linalg.norm(anti - target) < 1e-12

p = np.array([2.3, 0.4, -0.7, 1.1])
m = 0.9
slash = sum(gamma[mu] * p[mu] for mu in range(4))
factor = (slash - m * np.eye(4)) @ (slash + m * np.eye(4))
p2_minus_m2 = p @ eta @ p - m * m
err = np.linalg.norm(factor - p2_minus_m2 * np.eye(4))
print(f"    Clifford anticommutators verified; factorisation error = {err:.3e}")
assert err < 1e-12
print("    -> with a tetrad/spin frame, the scalar kernel factors into a Dirac operator.")
print("    -> causal order alone gives p^2 / intervals, not the gamma matrices, spinor")
print("       fibre, or spin connection.  The framework's E_{1/2} coin is a candidate")
print("       supplier of that missing spin frame, but that is added structure.")


print(
    r"""
[4] VERDICT
  R1 loop F^2:
    * raw link-loop action remains obstructed by link nonlocality and raw variance;
    * mesoscopic Alexandrov-loop holonomies are the viable manifestly gauge-invariant
      candidate, and the constant-field toy has the right convergence;
    * still open: non-arbitrary interval measure, 3+1 tensor contraction, variable-field
      locality expansion, overlap correlations, and normalisation.

  R2 raw variance/control:
    * raw unsmeared O(1)-cardinality loops are excluded as precision estimators;
    * the control condition is N_epsilon -> infinity while the physical interval size
      remains mesoscopic, exactly the B_epsilon-type smoothing logic.

  R3 Dirac spin:
    * compatible once a spin frame/tetrad/E_{1/2} coin is supplied;
    * not derived from causal order alone.  A positive theorem must construct a spinor
      fibre and spin connection from the framework's extra structure, not from order-only
      causal-set data.

  Net progress: the frontiers are sharpened, not solved.  The remaining hard theorem is
  a mesoscopic, gauge-invariant Alexandrov-loop Maxwell action plus a spin-frame lift.
exit 0"""
)
print("ALL ASSERTIONS PASSED -- frontier triage completed; no overclaim of loop-F2/Dirac closure.")
