#!/usr/bin/env python3
r"""Hawking-ladder wavepacket dynamical verification for item 10.

The exact spectral target was already fixed by ``bh_qec_observables.py``:
invalid horizon-register configurations form a 208-state subspace Q, with
integer Q3 strain F and degeneracies

    g_Q(F) = {0:1, 3:11, 4:22, 5:38, 6:54, 7:41, 8:25, 9:14, 12:2}.

This script checks the remaining "wavepacket" bracket without adding new
physics: a monitored local horizon wavepacket is dephased into populations on Q,
and the local single-bit KMS jump process on that same Q relaxes to the exact
Boltzmann ladder P(F) proportional to g_Q(F) exp(-beta F).

Scope:
  * closes the finite-state dynamical consistency check for the exact ladder;
  * does not derive the KMS assumption itself, nor a continuum Planck spectrum.
"""
from __future__ import annotations

from collections import Counter
import math

import numpy as np
from scipy.linalg import eigvals, expm


TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    """Canonical item-10 valid-register predicate from bh_qec_observables.py."""
    return (
        not (bit(n, 0) and bit(n, 1))
        and bit(n, 7) == bit(n, 6)
        and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0)))
    )


def strain(n: int) -> int:
    return sum(1 for i, j in EDGES if bit(n, i) != bit(n, j))


def invalid_subspace() -> list[int]:
    return [n for n in range(256) if not valid(n)]


def line_distribution(states: list[int], weights: np.ndarray) -> dict[int, float]:
    out: dict[int, float] = {f: 0.0 for f in TARGET_GQ}
    for n, w in zip(states, weights):
        out[strain(n)] += float(w)
    return dict(sorted(out.items()))


def local_kms_generator(states: list[int], beta: float) -> np.ndarray:
    """Row-convention CTMC generator for local one-bit horizon updates.

    For neighboring configurations i,j with one bit flipped, the rate
    W_ij = exp[-beta(F_j-F_i)/2] obeys detailed balance against
    pi_i proportional to exp[-beta F_i].  The graph is the induced 8-cube graph
    on Q, so this verifies the spectral ladder under local register moves.
    """
    index = {n: i for i, n in enumerate(states)}
    fvals = np.array([strain(n) for n in states], dtype=float)
    gen = np.zeros((len(states), len(states)), dtype=float)
    for i, n in enumerate(states):
        for bit_index in range(8):
            m = n ^ (1 << bit_index)
            j = index.get(m)
            if j is None:
                continue
            gen[i, j] = math.exp(-0.5 * beta * (fvals[j] - fvals[i]))
    gen[np.diag_indices_from(gen)] = -gen.sum(axis=1)
    return gen


def connected_one_bit_graph(states: list[int]) -> tuple[bool, int, int]:
    index = set(states)
    seen = {states[0]}
    stack = [states[0]]
    degrees = []
    while stack:
        n = stack.pop()
        neigh = [n ^ (1 << i) for i in range(8) if (n ^ (1 << i)) in index]
        degrees.append(len(neigh))
        for m in neigh:
            if m not in seen:
                seen.add(m)
                stack.append(m)
    return len(seen) == len(states), min(degrees), max(degrees)


def stationary_weights(states: list[int], beta: float) -> np.ndarray:
    fvals = np.array([strain(n) for n in states], dtype=float)
    raw = np.exp(-beta * fvals)
    return raw / raw.sum()


def deterministic_wavepacket(states: list[int], seed: int) -> np.ndarray:
    """A reproducible complex wavepacket; monitoring reads its probabilities."""
    rng = np.random.default_rng(seed)
    amp = rng.normal(size=len(states)) + 1j * rng.normal(size=len(states))
    amp /= np.linalg.norm(amp)
    probs = np.abs(amp) ** 2
    return probs / probs.sum()


def evolve_to_ladder(states: list[int], beta: float, seed: int) -> tuple[float, dict[int, float]]:
    gen = local_kms_generator(states, beta)
    pi = stationary_weights(states, beta)

    # Detailed balance gate.
    db = np.max(np.abs(pi[:, None] * gen - pi[None, :] * gen.T))
    assert db < 2.0e-14, f"detailed balance failed: {db}"

    # Row-generator spectrum: one zero eigenvalue, all other real parts negative.
    ev = np.sort(np.real(eigvals(gen)))
    gap = -ev[-2]
    assert gap > 1.0e-10, f"not mixing: spectral gap {gap}"
    assert abs(ev[-1]) < 1.0e-9, f"stationary eigenvalue not zero: {ev[-1]}"

    p0 = deterministic_wavepacket(states, seed)
    pt = p0 @ expm(gen * (45.0 / gap))
    pt = np.maximum(pt, 0.0)
    pt /= pt.sum()

    err = float(np.max(np.abs(pt - pi)))
    return err, line_distribution(states, pt)


def expected_line_distribution(beta: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return {f: v / z for f, v in sorted(raw.items())}


def main() -> None:
    states = invalid_subspace()
    gq = dict(sorted(Counter(strain(n) for n in states).items()))
    gv = dict(sorted(Counter(strain(n) for n in range(256) if valid(n)).items()))

    print("[1] Exact item-10 spectral target")
    print(f"    invalid states = {len(states)}; valid states = {256 - len(states)}")
    print(f"    g_invalid(F) = {gq}")
    print(f"    g_valid(F)   = {gv}")
    assert len(states) == 208
    assert gq == TARGET_GQ
    assert min(f for f in gq if f > 0) == 3
    assert 1 not in gq and 2 not in gq

    connected, min_deg, max_deg = connected_one_bit_graph(states)
    print("\n[2] Local wavepacket move graph")
    print(f"    induced one-bit graph connected = {connected}")
    print(f"    local degree range on Q = {min_deg}..{max_deg}")
    assert connected
    assert min_deg > 0

    print("\n[3] KMS wavepacket relaxation to the exact ladder")
    for beta, seed in ((0.25, 11), (0.70, 12), (1.00, 13), (1.60, 14)):
        err, got = evolve_to_ladder(states, beta, seed)
        want = expected_line_distribution(beta)
        line_err = max(abs(got[f] - want[f]) for f in want)
        ratio_43 = got[4] / got[3]
        expected_ratio_43 = (TARGET_GQ[4] / TARGET_GQ[3]) * math.exp(-beta)
        print(f"    beta={beta:>4.2f}: max state error {err:.2e}; max line error {line_err:.2e}")
        print(f"      P_F = {got}")
        print(
            f"      I4/I3 = {ratio_43:.8f}; "
            f"target (22/11) exp(-beta) = {expected_ratio_43:.8f}"
        )
        assert err < 5.0e-11
        assert line_err < 5.0e-11
        assert abs(ratio_43 / expected_ratio_43 - 1.0) < 1.0e-10

    print(
        "\nVERDICT: item 10's finite Hawking ladder has a local wavepacket "
        "dynamics: monitored wavepackets on the 208-state Q subspace, evolved "
        "by the induced one-bit KMS jump process, relax to the exact "
        "g_Q(F) exp(-beta F) spectrum. The emitting spectral gap is F=3; "
        "F=1 and F=2 lines are dynamically absent because no such strain "
        "states exist in Q. This verifies the dynamical bracket, conditional "
        "only on the same KMS horizon assumption already used by item 10."
    )


if __name__ == "__main__":
    main()
