#!/usr/bin/env python3
r"""Multi-plaquette correlation-length toy.

The two-plaquette script established a local triage rule:

    factorize, bound the connected tail, or promote to a joint record.

This script asks the next question: what happens in a chain of loop records?
Do correlations explode combinatorially, or do they decay so the record
description remains compact?

Toy ensemble
------------
Use a periodic chain of compact U(1) plaquette fluxes Phi_i with action

    S = beta sum_i [1 - cos(Phi_i)]
        + kappa sum_i [1 - cos(Phi_i - Phi_{i+1})].

The beta term pins each plaquette toward small flux.  The kappa term couples
neighbouring loop records.

Transfer-matrix calculation
---------------------------
The symmetric transfer matrix is

    T(Phi,Phi') =
      exp[ beta(cos Phi + cos Phi')/2 + kappa cos(Phi-Phi') ],

up to constants that cancel in normalized expectations.  On a periodic chain,

    <f(Phi_0) g(Phi_r)>
      = Tr[ F T^r G T^(L-r) ] / Tr[T^L].

That gives connected correlations without Monte Carlo noise.

Record-grammar reading
----------------------

    kappa = 0:        exact factorisation across all distances;
    finite short xi:  off-grammar tails are bounded by the correlation length;
    xi comparable L:  local records are no longer enough; promote the aligned
                      chain mode to a collective record.

This is still a compact U(1) toy, not a many-body QFT.  Its purpose is to make
the correlation-length criterion executable.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np


ArrayFn = Callable[[np.ndarray], np.ndarray]


def flux_grid(n: int) -> np.ndarray:
    """Midpoint compact grid on [-pi, pi)."""

    return np.array([-math.pi + (i + 0.5) * 2.0 * math.pi / n for i in range(n)], dtype=float)


@dataclass
class PlaquetteChain:
    """Exact-grid periodic U(1) chain solved by a symmetric transfer matrix."""

    beta: float
    kappa: float
    length: int = 64
    n_grid: int = 96
    phis: np.ndarray = field(init=False)
    eigvals: np.ndarray = field(init=False)
    eigvecs: np.ndarray = field(init=False)
    _powers: dict[int, np.ndarray] = field(default_factory=dict, init=False)

    def __post_init__(self) -> None:
        self.phis = flux_grid(self.n_grid)
        p1, p2 = np.meshgrid(self.phis, self.phis, indexing="ij")
        transfer = np.exp(0.5 * self.beta * (np.cos(p1) + np.cos(p2)) + self.kappa * np.cos(p1 - p2))
        vals, vecs = np.linalg.eigh(transfer)
        order = np.argsort(vals)[::-1]
        vals = vals[order]
        vecs = vecs[:, order]
        if vals[0] <= 0.0:
            raise ValueError("transfer matrix must have a positive Perron root")
        self.eigvals = vals / vals[0]
        self.eigvecs = vecs

    def power(self, exponent: int) -> np.ndarray:
        if exponent < 0:
            raise ValueError("negative transfer-matrix power")
        if exponent not in self._powers:
            self._powers[exponent] = (self.eigvecs * (self.eigvals**exponent)) @ self.eigvecs.T
        return self._powers[exponent]

    def partition(self) -> float:
        return float(np.sum(self.eigvals**self.length))

    def expectation(self, f: ArrayFn, g: ArrayFn, r: int) -> float:
        """Return <f(Phi_0) g(Phi_r)> on the periodic chain."""

        if r < 0 or r > self.length:
            raise ValueError("separation outside periodic chain")
        f_diag = np.diag(f(self.phis))
        g_diag = np.diag(g(self.phis))
        numerator = np.trace(f_diag @ self.power(r) @ g_diag @ self.power(self.length - r))
        return float(np.real_if_close(numerator / self.partition()))

    def one_point(self, f: ArrayFn) -> float:
        return self.expectation(f, lambda x: np.ones_like(x), 0)

    def covariance(self, f: ArrayFn, g: ArrayFn, r: int) -> float:
        return self.expectation(f, g, r) - self.one_point(f) * self.one_point(g)

    def spectral_correlation_length(self) -> float:
        """Largest transfer-matrix correlation length in lattice spacings."""

        ratio = float(abs(self.eigvals[1]))
        if ratio <= 0.0:
            return 0.0
        if ratio >= 1.0:
            return math.inf
        return -1.0 / math.log(ratio)

    def relative_alignment(self, r: int) -> float:
        """Return <cos(Phi_0 - Phi_r)> = <cos cos> + <sin sin>."""

        return self.expectation(np.cos, np.cos, r) + self.expectation(np.sin, np.sin, r)


def one_plaquette_mean_cos(beta: float, n_grid: int = 96) -> float:
    phis = flux_grid(n_grid)
    weights = np.exp(beta * np.cos(phis))
    probs = weights / float(np.sum(weights))
    return float(np.sum(probs * np.cos(phis)))


def fit_exponential_length(corrs: list[float], start_r: int = 1, stop_r: int = 8) -> tuple[float, float]:
    """Fit |C(r)| ~= A exp(-r/xi) over a short positive window."""

    rs = np.arange(start_r, stop_r + 1, dtype=float)
    ys = np.log(np.maximum(np.abs(np.array(corrs[start_r - 1 : stop_r], dtype=float)), 1e-300))
    slope, intercept = np.polyfit(rs, ys, 1)
    if slope >= 0.0:
        return math.inf, float(math.exp(intercept))
    return float(-1.0 / slope), float(math.exp(intercept))


def assert_close(name: str, value: float, target: float, tol: float = 1e-12) -> None:
    err = abs(value - target)
    print(f"  {name:<76s} value={value:.12g} target={target:.12g} err={err:.3e}")
    if err > tol:
        raise AssertionError(name)


def assert_less(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value < bound:
        raise AssertionError(name)


def assert_greater(name: str, value: float, bound: float) -> None:
    print(f"  {name:<76s} value={value:.12g} bound={bound:.12g}")
    if not value > bound:
        raise AssertionError(name)


def print_correlations(chain: PlaquetteChain, distances: list[int]) -> tuple[list[float], list[float]]:
    ccos: list[float] = []
    csin: list[float] = []
    for r in range(1, max(distances) + 1):
        ccos.append(chain.covariance(np.cos, np.cos, r))
        csin.append(chain.covariance(np.sin, np.sin, r))
    print("  r     C_cos(r)          C_sin(r)          <cos(Phi_0-Phi_r)>")
    for r in distances:
        print(
            f"  {r:<3d}  {ccos[r-1]: .12e}  {csin[r-1]: .12e}  "
            f"{chain.relative_alignment(r): .12e}"
        )
    return ccos, csin


def main() -> None:
    print("Plaquette-chain correlation-length toy: compactness by decay or promotion")
    print("=" * 96)

    length = 64
    n_grid = 96

    print("\n[1] kappa=0 factorises at every distance")
    independent = PlaquetteChain(beta=2.0, kappa=0.0, length=length, n_grid=n_grid)
    one = one_plaquette_mean_cos(beta=2.0, n_grid=n_grid)
    assert_close("<cos Phi> equals one-plaquette value", independent.one_point(np.cos), one)
    for r in [1, 2, 4, 8, 16, 32]:
        assert_close(f"C_cos({r}) at kappa=0", independent.covariance(np.cos, np.cos, r), 0.0, tol=5e-14)
        assert_close(f"C_sin({r}) at kappa=0", independent.covariance(np.sin, np.sin, r), 0.0, tol=5e-14)
    print("  -> no coupling means no hidden correlation tail anywhere on the chain.")

    print("\n[2] With local pinning, finite coupling gives a short correlation length")
    local = PlaquetteChain(beta=2.0, kappa=0.5, length=length, n_grid=n_grid)
    local_xi = local.spectral_correlation_length()
    local_ccos, local_csin = print_correlations(local, [1, 2, 4, 8, 16, 32])
    fit_xi, fit_a = fit_exponential_length(local_csin, 1, 8)
    print(f"  spectral xi={local_xi:.9f}; fitted sine-tail xi={fit_xi:.9f}; fitted A={fit_a:.9f}")
    assert_less("local spectral correlation length", local_xi, 0.75)
    assert_less("|C_sin(8)| local tail", abs(local_csin[7]), 1e-6)
    assert_less("|C_sin(16)| local tail", abs(local_csin[15]), 1e-10)
    print("  -> off-grammar correlations are not assumed away; they are bounded by xi.")

    print("\n[3] Increasing neighbour coupling lengthens the tail monotonically")
    previous_xi = -1.0
    for kappa in [0.1, 0.5, 2.0, 8.0]:
        chain = PlaquetteChain(beta=2.0, kappa=kappa, length=length, n_grid=n_grid)
        xi = chain.spectral_correlation_length()
        c8 = chain.covariance(np.sin, np.sin, 8)
        print(f"  beta=2.0 kappa={kappa:<3.1f}  xi={xi:.9f}  C_sin(8)={c8:.12e}")
        if xi + 1e-12 < previous_xi:
            raise AssertionError("correlation length decreased as kappa increased")
        previous_xi = xi
    print("  -> the transfer gap supplies a quantitative compactness criterion.")

    print("\n[4] Strong relative coupling with no local pinning creates a collective chain record")
    collective = PlaquetteChain(beta=0.0, kappa=20.0, length=length, n_grid=n_grid)
    collective_xi = collective.spectral_correlation_length()
    assert_close("<cos Phi> unpinned global angle", collective.one_point(np.cos), 0.0, tol=5e-13)
    assert_close("<sin Phi> unpinned global angle", collective.one_point(np.sin), 0.0, tol=5e-13)
    collective_ccos, collective_csin = print_correlations(collective, [1, 2, 4, 8, 16, 32])
    print(f"  spectral xi={collective_xi:.9f} on a length-{length} ring")
    assert_greater("collective spectral xi exceeds half-chain scale", collective_xi, length / 2.0)
    assert_greater("nearest-neighbour alignment", collective.relative_alignment(1), 0.95)
    assert_greater("half-chain alignment", collective.relative_alignment(length // 2), 0.60)
    assert_greater("half-chain C_sin", collective_csin[length // 2 - 1], 0.25)
    print("  -> relative phases are locked across the ring; the absolute flux remains unpinned.")

    print("\n[5] Record-grammar triage on a chain")
    assert_close("factorised chain covariance", independent.covariance(np.cos, np.cos, 8), 0.0, tol=5e-14)
    assert_less("short-xi chain tail", abs(local_csin[7]), 1e-6)
    assert_greater("collective chain long-distance covariance", collective_csin[length // 2 - 1], 0.25)
    print("  independent: exact factorisation")
    print("  short-xi:     bounded connected tail")
    print("  collective:   promote aligned chain mode to a record")

    print(
        """
Verdict
-------
The plaquette-chain toy turns correlation compactness into an executable test.

  kappa = 0:
      all connected correlations vanish across the chain.

  beta > 0 with finite transfer gap:
      connected correlations decay with a measured correlation length xi, so
      distant off-grammar tails can be bounded rather than guessed away.

  beta = 0 with strong neighbour coupling:
      the relative phases align across the whole chain while the global angle
      remains unpinned.  The compact description is then not 64 independent
      plaquette records, but one collective aligned-loop record plus its
      residual fluctuations.

Record-grammar rule:

  zero by factorisation;
  small by correlation-length bound;
  long-range by promotion to a collective record.

Boundary
--------
This is a one-dimensional compact U(1) transfer-matrix toy.  It is not a QFT
limit, not a non-abelian simulation, and not a proof that real matter always
has short correlations.  Its role is to make the required compression test
concrete: correlations must be shown to vanish, decay, or become named
collective records.
"""
    )


if __name__ == "__main__":
    main()
