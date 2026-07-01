#!/usr/bin/env python3
r"""Two-plaquette correlation toy: when loop records stay compact.

This script addresses the first "correlation triage" question for the record
grammar.

One plaquette is easy: the closed flux W = exp(i Phi) is a readable record.
Two plaquettes introduce a new question:

    when can their loop records be treated independently,
    and when must they be promoted to a collective record?

Toy ensemble
------------
Use two compact U(1) fluxes Phi_1, Phi_2 with action

    S = beta[(1-cos Phi_1) + (1-cos Phi_2)]
        + kappa[1-cos(Phi_1-Phi_2)].

The beta terms are local plaquette actions.  The kappa term is a local coupling
between neighbouring loop records.

Record-grammar reading
----------------------

    kappa = 0:    independent records; connected covariance is zero.
    small kappa: bounded weak correlation; keep an error budget.
    large kappa: Phi_1 ~= Phi_2; promote to one collective record.

This is not a lattice gauge simulation.  It is a tiny exact-grid calculation
showing how "neglecting correlations" must be replaced by factorisation,
explicit bounds, or promotion to a joint record.
"""

from __future__ import annotations

import math

import numpy as np


def flux_grid(n: int) -> np.ndarray:
    """Midpoint compact grid on [-pi, pi)."""

    return np.array([-math.pi + (i + 0.5) * 2.0 * math.pi / n for i in range(n)], dtype=float)


def action(phi1: np.ndarray, phi2: np.ndarray, beta: float, kappa: float) -> np.ndarray:
    return beta * ((1.0 - np.cos(phi1)) + (1.0 - np.cos(phi2))) + kappa * (1.0 - np.cos(phi1 - phi2))


def ensemble(beta: float, kappa: float, n: int = 320) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return Phi1 mesh, Phi2 mesh, and normalized probability grid."""

    phis = flux_grid(n)
    p1, p2 = np.meshgrid(phis, phis, indexing="ij")
    s = action(p1, p2, beta, kappa)
    weights = np.exp(-(s - float(np.min(s))))
    probs = weights / float(np.sum(weights))
    return p1, p2, probs


def mean(probs: np.ndarray, value: np.ndarray) -> float:
    return float(np.sum(probs * value))


def moments(beta: float, kappa: float, n: int = 320) -> dict[str, float]:
    p1, p2, probs = ensemble(beta, kappa, n=n)
    c1 = np.cos(p1)
    c2 = np.cos(p2)
    s1 = np.sin(p1)
    s2 = np.sin(p2)
    c12 = np.cos(p1 - p2)
    m1 = mean(probs, c1)
    m2 = mean(probs, c2)
    return {
        "cos1": m1,
        "cos2": m2,
        "sin1": mean(probs, s1),
        "sin2": mean(probs, s2),
        "cos12": mean(probs, c12),
        "cov_cos": mean(probs, c1 * c2) - m1 * m2,
        "cov_sin": mean(probs, s1 * s2) - mean(probs, s1) * mean(probs, s2),
        "mean_abs_delta": mean(probs, np.abs(np.angle(np.exp(1j * (p1 - p2))))),
    }


def one_plaquette_mean_cos(beta: float, n: int = 320) -> float:
    phis = flux_grid(n)
    weights = np.exp(beta * np.cos(phis))
    probs = weights / float(np.sum(weights))
    return float(np.sum(probs * np.cos(phis)))


def assert_close(name: str, value: float, target: float, tol: float = 1e-10) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def main() -> None:
    print("Two-plaquette correlation toy: compact records, bounded tails")
    print("=" * 92)

    beta = 2.0
    print("\n[1] kappa=0 factorises exactly on the grid")
    m0 = moments(beta=beta, kappa=0.0)
    one = one_plaquette_mean_cos(beta)
    assert_close("<cos Phi1>", m0["cos1"], one)
    assert_close("<cos Phi2>", m0["cos2"], one)
    assert_close("Cov(cos Phi1, cos Phi2)", m0["cov_cos"], 0.0, tol=1e-14)
    assert_close("Cov(sin Phi1, sin Phi2)", m0["cov_sin"], 0.0, tol=1e-14)
    print("  -> with no coupling, two loop records are independent; no hidden covariance is being ignored.")

    print("\n[2] Weak coupling produces a small, measured connected correlation")
    weak = moments(beta=beta, kappa=0.1)
    print(
        f"  kappa=0.1: <cos12>={weak['cos12']:.9f}, "
        f"cov_cos={weak['cov_cos']:.9f}, cov_sin={weak['cov_sin']:.9f}"
    )
    assert_less("|cov_cos| at weak coupling", abs(weak["cov_cos"]), 0.025)
    assert_less("|cov_sin| at weak coupling", abs(weak["cov_sin"]), 0.025)
    print("  -> weak off-grammar correlation is not discarded; it is bounded.")

    print("\n[3] Correlations grow monotonically with local coupling")
    kappas = [0.0, 0.1, 0.5, 2.0, 8.0]
    previous_cos12 = -1.0
    rows: list[tuple[float, dict[str, float]]] = []
    for kappa in kappas:
        m = moments(beta=beta, kappa=kappa)
        rows.append((kappa, m))
        print(
            f"  kappa={kappa:<4.1f}  <cos(Phi1-Phi2)>={m['cos12']:.9f}  "
            f"cov_cos={m['cov_cos']:.9f}  cov_sin={m['cov_sin']:.9f}  "
            f"<|delta|>={m['mean_abs_delta']:.9f}"
        )
        if m["cos12"] + 1e-12 < previous_cos12:
            raise AssertionError("<cos(Phi1-Phi2)> decreased as kappa increased")
        previous_cos12 = m["cos12"]
    print("  -> the local coupling turns two independent records into a correlated pair continuously.")

    print("\n[4] Strong coupling promotes the pair to a collective loop record")
    strong = rows[-1][1]
    assert_less("mean wrapped phase separation at kappa=8", strong["mean_abs_delta"], 0.35)
    assert_close("Phi1 and Phi2 one-point equality", strong["cos1"], strong["cos2"], tol=1e-12)
    print("  -> for large kappa, the compact description is one collective record Phi1 ~= Phi2.")

    print("\n[5] Connected correlation triage")
    independent = rows[0][1]
    weak = rows[1][1]
    strong = rows[-1][1]
    assert_close("independent connected covariance", independent["cov_cos"], 0.0, tol=1e-14)
    assert_less("weak connected covariance kept as tail", abs(weak["cov_cos"]), 0.025)
    if abs(strong["cov_sin"]) < 0.1:
        raise AssertionError("strong-coupling sine covariance too small to require promotion")
    print("  independent: exact factorisation")
    print("  weak:        bounded connected tail")
    print("  strong:      promote to joint/collective loop record")

    print(
        """
Verdict
-------
The two-plaquette toy gives the first executable correlation triage:

  kappa = 0    -> exact factorisation, connected covariance zero;
  small kappa  -> explicit bounded connected tail;
  large kappa  -> collective record Phi_1 ~= Phi_2.

So the record grammar does not license simply ignoring correlations.  It
requires one of three outcomes: prove factorisation, bound the residual, or
promote the correlated variables to a joint record.

Boundary
--------
This is a compact U(1) two-variable ensemble, not a many-body gauge simulation.
Its purpose is to make the triage rule executable before applying it to larger
record networks.
"""
    )


if __name__ == "__main__":
    main()
