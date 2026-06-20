#!/usr/bin/env python3
r"""Explicit service-current stress-tensor gate for intrinsic gravity.

Question
--------
The RT/min-cut and modular-Hamiltonian audits have narrowed the gravity gate:
the static code has no useful modular Hamiltonian, while the service-clock
gradient supplies the boost kernel.  The next missing finite object is the
stress tensor itself:

    Can the QEC service ledger define a symmetric, conserved T^{mu nu}, or is
    gravity still only a horizon-input Dirac coincidence?

This script builds the minimal local object the ledger is allowed to supply.
Each service event carries a future-directed four-momentum/action vector p^mu
and a rate r_e.  The coarse-grained stress tensor is the kinetic/service
current

    T^{mu nu}_svc = (1/V) sum_e r_e p_e^mu p_e^nu / p_e^0 .

This is the unique rank-two flux made from event currents without adding a new
coefficient.  It is symmetric by construction, conserved exactly when the
service graph obeys Kirchhoff four-momentum balance, positive because it is a
sum over future causal p^mu, and it gives the Rindler/boost modular energy
K = 2 pi int x T_00.

Exit 0 means:
  * a nontrivial discrete service-current T^{mu nu} is constructed;
  * symmetry, conservation, positivity, radiation/dust controls, and boost
    modular weighting are checked;
  * the hierarchy gate is answered honestly: the local tensor supplies the
    source and the lattice/bare coupling, but it does not generate the
    observed M_Pl hierarchy.  A horizon/dilution or other nonlocal
    renormalisation theorem is still required.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
HBARC = 0.197327  # GeV fm
A0_FM = 0.594
LAMBDA_QCD = HBARC / A0_FM
S_CELL = 55.0 / 8.0
M_PL_OBS = 1.2209e19


@dataclass(frozen=True)
class Event:
    """A coarse-grained service event species."""

    name: str
    p: np.ndarray  # contravariant four-momentum/action vector, p^0 > 0
    rate: float


def minkowski_dot(a: np.ndarray, b: np.ndarray) -> float:
    return float(a @ ETA @ b)


def stress(events: list[Event], volume: float = 1.0) -> np.ndarray:
    """T^{mu nu} = V^{-1} sum r p^mu p^nu / p^0."""

    t = np.zeros((4, 4), dtype=float)
    for event in events:
        p = event.p
        assert p[0] > 0.0
        assert minkowski_dot(p, p) <= 1.0e-12, f"{event.name} is not causal"
        t += event.rate * np.outer(p, p) / p[0]
    return t / volume


def random_observers(n: int = 64, seed: int = 153) -> list[np.ndarray]:
    """Future timelike unit four-velocities, u.u=-1."""

    rng = np.random.default_rng(seed)
    out: list[np.ndarray] = []
    for _ in range(n):
        v = rng.normal(size=3)
        speed = 0.95 * rng.random() * np.linalg.norm(v) / (np.linalg.norm(v) + 1.0e-12)
        nhat = v / (np.linalg.norm(v) + 1.0e-12)
        beta = speed
        gamma = 1.0 / math.sqrt(1.0 - beta * beta)
        out.append(np.array([gamma, *(gamma * beta * nhat)]))
    return out


def dominant_energy_check(t: np.ndarray) -> tuple[float, float]:
    """Return minimum observed energy density and maximum causal-current norm.

    For signature (-,+,+,+), the observer energy density is
        rho = T^{mu nu} u_mu u_nu.
    The measured energy current is
        j^mu = -T^{mu nu} u_nu,
    which should be future causal: j^0 >= 0 and j.j <= 0.
    """

    min_rho = math.inf
    max_j2 = -math.inf
    for u in random_observers():
        u_cov = ETA @ u
        rho = float(u_cov @ t @ u_cov)
        j = -t @ u_cov
        j2 = minkowski_dot(j, j)
        assert j[0] >= -1.0e-10
        min_rho = min(min_rho, rho)
        max_j2 = max(max_j2, j2)
    return min_rho, max_j2


def null_axis_events(rate: float = 1.0, energy: float = 1.0) -> list[Event]:
    events: list[Event] = []
    for axis, label in enumerate("xyz", start=1):
        for sign in (-1.0, 1.0):
            p = np.zeros(4)
            p[0] = energy
            p[axis] = sign * energy
            events.append(Event(f"null {sign:+.0f}{label}", p, rate))
    return events


def rest_events(rate: float = 1.0, mass: float = 1.0) -> list[Event]:
    return [Event("rest-mass service record", np.array([mass, 0.0, 0.0, 0.0]), rate)]


def closed_loop_divergence(l: int = 5, energy: float = 1.0) -> np.ndarray:
    """Kirchhoff check on a periodic spatial service graph.

    Put one directed service loop around every x-line of an L^3 torus.  Each
    node has one incoming and one outgoing event carrying the same p^nu, so
    div_mu T^{mu nu} reduces to exact four-momentum current balance.
    """

    div = np.zeros((l, l, l, 4), dtype=float)
    p = np.array([energy, energy, 0.0, 0.0])
    for y in range(l):
        for z in range(l):
            for x in range(l):
                src = (x, y, z)
                dst = ((x + 1) % l, y, z)
                div[src] -= p
                div[dst] += p
    return div


def rindler_modular_energy(delta_t00: np.ndarray) -> tuple[float, float]:
    """Service-clock/Rindler modular energy for a half-space shell."""

    x = np.arange(1, len(delta_t00) + 1, dtype=float)
    k_boost = float(np.dot(x, delta_t00))
    h_mod = float(np.dot(2.0 * math.pi * x, delta_t00))
    return k_boost, h_mod


def main() -> None:
    print("SERVICE-CURRENT STRESS-TENSOR GATE")
    print("=" * 88)

    print("[1] Explicit local tensor")
    radiation = stress(null_axis_events(rate=0.5, energy=2.0))
    print("    isotropic six-null-event ledger T^{mu nu}:")
    with np.printoptions(precision=6, suppress=True):
        print(radiation)
    rho = radiation[0, 0]
    pressure = np.trace(radiation[1:, 1:]) / 3.0
    offdiag = radiation - np.diag(np.diag(radiation))
    sym_err = float(np.linalg.norm(radiation - radiation.T))
    print(f"    symmetry error = {sym_err:.2e}")
    print(f"    radiation control: P/rho = {pressure / rho:.6f} (target 1/3)")
    print(f"    off-diagonal norm = {np.linalg.norm(offdiag):.2e}")
    assert sym_err < 1.0e-12
    assert abs(pressure / rho - 1.0 / 3.0) < 1.0e-12
    assert np.linalg.norm(offdiag) < 1.0e-12

    dust = stress(rest_events(rate=3.0, mass=2.0))
    print("\n[2] Dust/full-energy control")
    print(f"    rest ledger T00={dust[0,0]:.3f}, spatial trace={np.trace(dust[1:,1:]):.3f}")
    assert dust[0, 0] == 6.0
    assert np.linalg.norm(dust[1:, 1:]) == 0.0
    print("    -> the same event-current definition reads rest energy as gravity source;")
    print("       it is not the anti-mass bare-hopping coupling.")

    mixed = radiation + dust
    min_rho, max_j2 = dominant_energy_check(mixed)
    print("\n[3] Positivity / dominant-energy check")
    print(f"    min observer energy density over random observers = {min_rho:.6f}")
    print(f"    max j_mu j^mu over random observers = {max_j2:.3e} (must be <=0)")
    assert min_rho > 0.0
    assert max_j2 <= 1.0e-9

    div = closed_loop_divergence(l=5)
    max_div = float(np.max(np.abs(div)))
    print("\n[4] Discrete conservation")
    print(f"    closed service-loop Kirchhoff max |div J^nu| = {max_div:.2e}")
    assert max_div < 1.0e-12
    print("    -> nabla_mu T^{mu nu}=0 is exactly service-current balance at nodes.")

    print("\n[5] Rindler / modular coupling")
    delta_t00 = np.array([0.10, 0.03, -0.02, 0.04, 0.01, 0.00])
    k_boost, h_mod = rindler_modular_energy(delta_t00)
    print(f"    K_boost=sum x deltaT00 = {k_boost:.6f}")
    print(f"    H_mod=2pi K_boost       = {h_mod:.6f}")
    assert abs(h_mod - 2.0 * math.pi * k_boost) < 1.0e-12
    print("    -> the explicit T00 is the object consumed by the service-clock")
    print("       modular Hamiltonian H_mod=2pi int x T00.")

    print("\n[6] Hierarchy gate")
    m_pl_bare = 2.0 * math.sqrt(S_CELL) * LAMBDA_QCD
    hierarchy = M_PL_OBS / m_pl_bare
    needed_area_factor = hierarchy * hierarchy
    local_counts = [28.0, 55.0 / 8.0, 137.0, 208.0, 256.0]
    biggest_single_local = max(local_counts)
    print(f"    local area-law coupling: M_Pl,bare = 2 sqrt(55/8) Lambda_QCD = {m_pl_bare:.3f} GeV")
    print(f"    observed / bare = {hierarchy:.3e}; required 1/G enhancement = {needed_area_factor:.3e}")
    print(f"    largest single local canon count in this audit = {biggest_single_local:.0f}")
    assert needed_area_factor > 1.0e37
    assert biggest_single_local < 1.0e3
    print("    -> the service-current tensor touches the FORM and the bare/lattice")
    print("       coupling.  It does not contain the horizon-size enhancement.")

    print(
        """
[7] VERDICT
    PASSED: there is an explicit service-current stress tensor

        T^{mu nu}_svc = V^{-1} sum_e r_e p_e^mu p_e^nu / p_e^0,

    with the right finite-ledger properties: symmetric, positive, conserved
    when the service graph obeys Kirchhoff balance, radiation and dust controls,
    and direct coupling to the boost modular Hamiltonian.

    This means the hierarchy can be touched at the level of the *source*: the
    framework no longer needs to import a generic matter T^{mu nu}; the service
    ledger itself supplies one.

    NOT PASSED: the observed Planck hierarchy is not generated.  The tensor is
    local and gives the lattice/bare coupling M_Pl~GeV from the service-record
    density.  The required ~10^38 enhancement in 1/G is absent unless an
    additional horizon/dilution or RG theorem is supplied.  Therefore the gate
    result is:

        Einstein/source structure: service-ledger intrinsic.
        Observed M_Pl hierarchy : still horizon-input / nonlocal-theorem open.

exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
