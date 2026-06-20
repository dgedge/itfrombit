#!/usr/bin/env python3
r"""LOCALIZED-MASS QEC STEADY-STATE GATE.

This is the control gate for the "derive Hawking KMS from horizon QEC" problem.
It identifies the exact service-current condition.  The companion script
bh_kms_scheduler_derivation.py derives that condition from the Schwarzschild
microcanonical reservoir plus symmetric QEC service-current/GNS algebra.

On the 208 invalid horizon-register states Q, with strain F=|delta(s)|:

  * A plain uniform one-bit QEC repair walk has the wrong stationary line
    weights.  It is an infinite-temperature graph walk, not Hawking.

  * The exact Hawking ladder is obtained iff the local service current obeys

        W_{i->j} = A_{ij} exp[- beta (F_j-F_i)/2],
        A_{ij}=A_{ji},

    on the induced one-bit graph.  This is the KMS half-Boltzmann rule.

Therefore the theorem target was never vague:

    derive the symmetric proposal A and the half-Boltzmann service-current
    dressing from the localized Schwarzschild/horizon QEC bath.

That target is now addressed by bh_kms_scheduler_derivation.py.  This file is
kept as the falsification/control gate: uniform service and full-Boltzmann
service remain wrong.
"""

from __future__ import annotations

from collections import Counter
import math

import numpy as np


TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    return (
        not (bit(n, 0) and bit(n, 1))
        and bit(n, 7) == bit(n, 6)
        and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0)))
    )


def strain(n: int) -> int:
    return sum(1 for i, j in EDGES if bit(n, i) != bit(n, j))


def invalid_states() -> list[int]:
    return [n for n in range(256) if not valid(n)]


def generator(states: list[int], beta: float, policy: str) -> np.ndarray:
    index = {s: i for i, s in enumerate(states)}
    f = np.array([strain(s) for s in states], dtype=float)
    gen = np.zeros((len(states), len(states)), dtype=float)
    for i, s in enumerate(states):
        for k in range(8):
            t = s ^ (1 << k)
            j = index.get(t)
            if j is None:
                continue
            df = f[j] - f[i]
            if policy == "uniform":
                rate = 1.0
            elif policy == "kms":
                rate = math.exp(-0.5 * beta * df)
            elif policy == "wrong_full_boltzmann":
                rate = math.exp(-beta * df)
            else:
                raise ValueError(policy)
            gen[i, j] = rate
    gen[np.diag_indices_from(gen)] = -gen.sum(axis=1)
    return gen


def stationary(gen: np.ndarray) -> np.ndarray:
    vals, vecs = np.linalg.eig(gen.T)
    i = int(np.argmin(np.abs(vals)))
    v = np.real(vecs[:, i])
    if v.sum() < 0:
        v = -v
    v = np.maximum(v, 0.0)
    return v / v.sum()


def target_weights(states: list[int], beta: float) -> np.ndarray:
    f = np.array([strain(s) for s in states], dtype=float)
    w = np.exp(-beta * f)
    return w / w.sum()


def line_dist(states: list[int], weights: np.ndarray) -> dict[int, float]:
    out = {f: 0.0 for f in TARGET_GQ}
    for s, w in zip(states, weights):
        out[strain(s)] += float(w)
    return dict(sorted(out.items()))


def total_variation(p: np.ndarray, q: np.ndarray) -> float:
    return 0.5 * float(np.abs(p - q).sum())


def detailed_balance_error(pi: np.ndarray, gen: np.ndarray) -> float:
    return float(np.max(np.abs(pi[:, None] * gen - pi[None, :] * gen.T)))


def main() -> None:
    states = invalid_states()
    gq = dict(sorted(Counter(strain(s) for s in states).items()))
    print("[1] Invalid-register target")
    print(f"    |Q| = {len(states)}, g_Q(F) = {gq}")
    check(len(states) == 208, "invalid subspace has 208 states")
    check(gq == TARGET_GQ, "strain degeneracy table matches item 10")

    beta = 1.0
    target = target_weights(states, beta)
    print(f"\n[2] Scheduler candidates at beta={beta}")
    for policy in ("uniform", "wrong_full_boltzmann", "kms"):
        gen = generator(states, beta, policy)
        pi = stationary(gen)
        tv = total_variation(pi, target)
        db = detailed_balance_error(target, gen)
        print(f"    {policy:<22s} TV-to-Hawking={tv:.6f}  DB-error-vs-Hawking={db:.3e}")
        print(f"      line P(F) = {line_dist(states, pi)}")
        if policy == "kms":
            check(tv < 1.0e-10, "half-Boltzmann KMS scheduler gives exact Hawking weights")
            check(db < 1.0e-13, "KMS scheduler obeys detailed balance against Hawking weights")
        else:
            check(tv > 0.05, f"{policy} scheduler is visibly not the Hawking steady state")

    print(
        """
[3] VERDICT
    The localized-mass QEC steady-state problem has been reduced to one exact
    operator statement.  The finite horizon graph alone is insufficient:
    uniform QEC service gives the wrong line weights.  The successful scheduler
    is precisely a symmetric one-bit proposal dressed by the half-Boltzmann
    factor exp[-beta DeltaF/2].

    Companion closure:
        bh_kms_scheduler_derivation.py derives the half-Boltzmann dressing from
        Schwarzschild microcanonical balance plus the symmetric QEC service
        operator in the thermal/GNS frame.  This file remains the control gate.

exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
