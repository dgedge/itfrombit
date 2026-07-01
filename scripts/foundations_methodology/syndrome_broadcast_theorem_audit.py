#!/usr/bin/env python3
r"""Quantum Darwinism / syndrome broadcast theorem audit.

This script checks the finite inequalities used by the noisy syndrome-broadcast
paper.  It is not a simulation of a microscopic bath.  It is a theorem witness
for the abstract stabilizer-syndrome extraction model:

  * a classical syndrome variable S with distribution p_s;
  * N conditionally independent record fragments;
  * fragment i reports the correct syndrome with probability 1-eps_i and a
    uniformly wrong label otherwise;
  * reset/Landauer costs are computed from the fragment record distributions.

Exit 0 means:

  1. the noisy state is within the claimed total-variation / trace-distance
     bound from the perfect spectrum-broadcast state;
  2. each fragment's mutual information obeys the Fano lower bound;
  3. conditional independence I(F_A:F_B|S)=0 holds;
  4. the reset ledger separates local reset cost from globally correlated reset
     cost, and both are bounded by Shannon entropies.
"""

from __future__ import annotations

from itertools import product
import math

import numpy as np


TOL = 1.0e-12


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def h(probs: np.ndarray) -> float:
    probs = np.asarray(probs, dtype=float)
    probs = probs[probs > 1.0e-15]
    return float(-np.sum(probs * np.log(probs)))


def h2(e: float) -> float:
    if e <= 0.0 or e >= 1.0:
        return 0.0
    return -e * math.log(e) - (1.0 - e) * math.log(1.0 - e)


def fano_deficit(e: float, m: int) -> float:
    if m <= 1:
        return 0.0
    return h2(e) + e * math.log(m - 1)


def channel(m: int, eps: float) -> np.ndarray:
    """Rows are syndrome s; columns are record y."""

    w = np.full((m, m), eps / (m - 1), dtype=float)
    np.fill_diagonal(w, 1.0 - eps)
    return w


def mutual_info_joint(joint: np.ndarray) -> float:
    ps = np.sum(joint, axis=1)
    py = np.sum(joint, axis=0)
    out = 0.0
    for s in range(joint.shape[0]):
        for y in range(joint.shape[1]):
            if joint[s, y] > 0.0:
                out += joint[s, y] * math.log(joint[s, y] / (ps[s] * py[y]))
    return float(out)


def joint_distribution(p: np.ndarray, channels: list[np.ndarray]) -> np.ndarray:
    """Return P(s,y_1,...,y_N)."""

    m = len(p)
    shape = (m,) + (m,) * len(channels)
    out = np.zeros(shape, dtype=float)
    for s in range(m):
        for ys in product(range(m), repeat=len(channels)):
            prob = p[s]
            for i, y in enumerate(ys):
                prob *= channels[i][s, y]
            out[(s,) + ys] = prob
    return out


def perfect_distribution(p: np.ndarray, n: int) -> np.ndarray:
    m = len(p)
    shape = (m,) + (m,) * n
    out = np.zeros(shape, dtype=float)
    for s in range(m):
        out[(s,) + (s,) * n] = p[s]
    return out


def total_variation(a: np.ndarray, b: np.ndarray) -> float:
    return 0.5 * float(np.sum(np.abs(a - b)))


def marginal(dist: np.ndarray, keep: tuple[int, ...]) -> np.ndarray:
    axes = tuple(i for i in range(dist.ndim) if i not in keep)
    return np.sum(dist, axis=axes)


def conditional_mutual_info_fragments(dist: np.ndarray, a_axis: int, b_axis: int) -> float:
    """Compute I(F_a:F_b|S) for a classical joint distribution."""

    s_axis = 0
    ps = marginal(dist, (s_axis,))
    out = 0.0
    m = dist.shape[0]
    for s in range(m):
        if ps[s] <= 0.0:
            continue
        # axes are S,F1,...,FN.  Slice S=s, then select fragment axes shifted down by 1.
        cond = np.take(dist, s, axis=s_axis) / ps[s]
        fa = a_axis - 1
        fb = b_axis - 1
        pab = marginal(cond, (fa, fb))
        pa = marginal(cond, (fa,))
        pb = marginal(cond, (fb,))
        for a in range(pab.shape[0]):
            for b in range(pab.shape[1]):
                if pab[a, b] > 0.0:
                    out += ps[s] * pab[a, b] * math.log(pab[a, b] / (pa[a] * pb[b]))
    return float(out)


def main() -> int:
    print("SYNDROME BROADCAST THEOREM AUDIT")
    print("=" * 88)

    p = np.array([0.55, 0.25, 0.15, 0.05], dtype=float)
    eps = np.array([0.01, 0.03, 0.05, 0.06, 0.12, 0.20], dtype=float)
    m = len(p)
    n = len(eps)
    channels = [channel(m, float(e)) for e in eps]
    Hs = h(p)

    print(f"\n[0] syndrome model: m={m}, N={n}, H(S)={Hs:.9f} nats")
    print(f"    eps_i = {', '.join(f'{x:.2f}' for x in eps)}")

    print("\n[1] local Fano bounds")
    witness_count = 0
    delta = 0.25
    for i, w in enumerate(channels, 1):
        joint = p[:, None] * w
        info = mutual_info_joint(joint)
        e = 1.0 - float(np.sum(p * np.diag(w)))
        deficit = fano_deficit(e, m)
        lower = Hs - deficit
        is_witness = info >= (1.0 - delta) * Hs
        witness_count += int(is_witness)
        print(
            f"  F{i}: err={e:.6f}, I(S:F_i)={info:.9f}, "
            f"Fano lower={lower:.9f}, witness_delta={is_witness}"
        )
        check(info + TOL >= lower, f"fragment {i} obeys Fano lower bound")
    check(witness_count >= 4, "at least four single fragments meet the 25 percent redundancy threshold")

    print("\n[2] approximate spectrum-broadcast trace-distance bounds")
    noisy = joint_distribution(p, channels)
    perfect = perfect_distribution(p, n)
    tv = total_variation(noisy, perfect)
    exact_any_error = 1.0 - float(np.prod(1.0 - eps))
    union_bound = float(np.sum(eps))
    print(f"  D(noisy, perfect SBS) = {tv:.9f}")
    print(f"  exact any-error probability = {exact_any_error:.9f}")
    print(f"  union bound sum eps_i      = {union_bound:.9f}")
    check(abs(tv - exact_any_error) < TOL, "trace distance equals probability of at least one wrong copy")
    check(tv <= union_bound + TOL, "trace distance obeys union bound")
    for i, e in enumerate(eps, 1):
        local_noisy = marginal(noisy, (0, i))
        local_perfect = marginal(perfect, (0, i))
        d_local = total_variation(local_noisy, local_perfect)
        print(f"  local D(SF{i}, ideal) = {d_local:.9f}")
        check(abs(d_local - e) < TOL, f"local trace distance equals eps_{i}")

    print("\n[3] conditional independence / strong-Darwinism structure")
    cmi_12 = conditional_mutual_info_fragments(noisy, 1, 2)
    cmi_25 = conditional_mutual_info_fragments(noisy, 2, 5)
    print(f"  I(F1:F2 | S) = {cmi_12:.3e} nats")
    print(f"  I(F2:F5 | S) = {cmi_25:.3e} nats")
    check(abs(cmi_12) < TOL, "fragments 1 and 2 are conditionally independent given syndrome")
    check(abs(cmi_25) < TOL, "fragments 2 and 5 are conditionally independent given syndrome")

    print("\n[4] reset / Landauer ledger")
    local_reset = 0.0
    for i, w in enumerate(channels, 1):
        q = p @ w
        hi = h(q)
        local_reset += hi
        print(f"  F{i}: H(record_i)={hi:.9f} nats")
    joint_records = marginal(noisy, tuple(range(1, n + 1)))
    global_reset = h(joint_records.reshape(-1))
    perfect_local_reset = n * Hs
    perfect_global_reset = Hs
    print(f"  local reset ledger   sum_i H(F_i) = {local_reset:.9f} nats")
    print(f"  global reset ledger  H(F_1...F_N) = {global_reset:.9f} nats")
    print(f"  perfect-copy local/global = {perfect_local_reset:.9f} / {perfect_global_reset:.9f} nats")
    check(local_reset + TOL >= global_reset, "local reset cost is at least globally compressed reset cost")
    check(perfect_local_reset > perfect_global_reset, "perfect redundant records cost more if reset locally")
    check(global_reset >= Hs - TOL, "joint record reset at least stores the syndrome entropy")

    print("\n[5] Decision")
    print(
        """
  The finite noisy-copy model obeys the theorem ledger:
    * single fragments carry syndrome information bounded below by Fano;
    * the global state is close to perfect spectrum broadcast by the any-error
      probability, and locally close by eps_i;
    * fragments are conditionally independent given the syndrome;
    * reset costs are ledger-dependent: local resets cost sum_i H(F_i), while
      a globally compressed reset costs H(F_1...F_N).
"""
    )
    print("ALL ASSERTIONS PASSED -- noisy syndrome-broadcast theorem inequalities verified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
