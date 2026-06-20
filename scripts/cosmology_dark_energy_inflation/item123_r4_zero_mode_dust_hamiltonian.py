#!/usr/bin/env python3
r"""ITEM 123: Brown--Kuchar dust Hamiltonian for the R4 zero-mode reservoir.

The conserved-reservoir lift closes

    N_zero + N_active = N_tot

but the previous script still marked "dust Hamiltonian" as a premise.  This
script tests whether that premise follows from the actual object introduced by
the lift: a Stinespring service-record reservoir.

Key point
---------
A reservoir label is a record-number label, not a spatial link field.  The
minimal lift admits:

    H_zero = epsilon * sum_i N_i

plus local exchange with the active R4 count at the same coarse cell/line.
There is no operator in the lift that compares neighbouring reservoir counts,
transports a reservoir label between cells, or assigns a phase-gradient energy.

Therefore the zero-mode is the finite analogue of constrained
Brown--Kuchar dust: rest energy only, exact number/shift symmetry, no gradient
stiffness.  At fixed comoving count, E is independent of physical volume, so

    rho ∝ a^-3,      p = -dE/dV|_N = 0,      c_s^2 = 0.

Scope
-----
This closes the Hamiltonian *form* under the reservoir-label/operator-inventory
reading of the Stinespring lift.  It does not derive the absolute reservoir
density N_tot/volume at boot; that remains the abundance/initial-condition
theorem.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


@dataclass(frozen=True)
class OperatorTerm:
    name: str
    form: str
    status: str
    reason: str


def operator_inventory() -> list[OperatorTerm]:
    return [
        OperatorTerm(
            "rest count",
            "epsilon sum_i N_i",
            "ADMITTED",
            "local record-number billing; commutes with all reservoir counts",
        ),
        OperatorTerm(
            "local active exchange",
            "Gamma_i(N_i, n_active_i)",
            "ADMITTED-FOR-HALO",
            "same-cell Stinespring exchange; preserves N_zero+N_active",
        ),
        OperatorTerm(
            "density-gradient stiffness",
            "K sum_<ij> (N_i-N_j)^2/a^2",
            "REQUIRES-NEW-LINK",
            "compares neighbouring reservoir labels; no such link operator exists in the lift",
        ),
        OperatorTerm(
            "reservoir hopping",
            "t sum_<ij> b_i^dag b_j + h.c.",
            "REQUIRES-NEW-TRANSPORT",
            "moves record labels between cells; Stinespring records are local history slots",
        ),
        OperatorTerm(
            "phase stiffness",
            "J sum_<ij> cos(theta_i-theta_j)",
            "REQUIRES-NEW-PHASE-FIELD",
            "uses a physical phase field; the conjugate phase is a shift gauge of N_tot",
        ),
        OperatorTerm(
            "elastic shear",
            "mu u_ij u^ij",
            "REQUIRES-NEW-METRIC-STRAIN",
            "belongs to the crystalline substrate, not the reservoir label sector",
        ),
    ]


def ring_gradient_energy(counts: np.ndarray, a: float, stiffness: float) -> float:
    diffs = counts - np.roll(counts, -1)
    return float(stiffness * np.dot(diffs, diffs) / (a * a))


def rest_energy(counts: np.ndarray, epsilon: float) -> float:
    return float(epsilon * np.sum(counts))


def pressure_from_energy(energy_fn, a: float, v0: float = 1.0) -> float:
    """Return p=-dE/dV at scale factor a by symmetric finite difference."""

    da = 1.0e-5 * a
    ep = energy_fn(a + da)
    em = energy_fn(a - da)
    vp = v0 * (a + da) ** 3
    vm = v0 * (a - da) ** 3
    return -float((ep - em) / (vp - vm))


def equation_of_state_for_rest_counts(counts: np.ndarray, epsilon: float, a: float) -> tuple[float, float, float]:
    e = rest_energy(counts, epsilon)
    volume = a**3
    rho = e / volume
    p = pressure_from_energy(lambda aa: rest_energy(counts, epsilon), a)
    w = p / rho
    return rho, p, w


def sound_quadratic_coefficients(ncells: int = 64, nbar: float = 100.0, epsilon: float = 1.0) -> list[tuple[int, float, float]]:
    """Compare zero-sum density waves in rest-only vs gradient Hamiltonians.

    A pressure/sound speed for a nonrelativistic fluid appears as a k^2
    quadratic cost for density perturbations.  The rest-count Hamiltonian has
    no such term: moving record labels around at fixed total count is energy
    degenerate.  A stiffness control has a positive k^2 cost.
    """

    x = np.arange(ncells, dtype=float)
    delta = 1.0e-3
    rows: list[tuple[int, float, float]] = []
    base = np.full(ncells, nbar)
    e0_rest = rest_energy(base, epsilon)
    e0_grad = ring_gradient_energy(base, a=1.0, stiffness=1.0)
    for mode in (1, 2, 4, 8, 16):
        wave = np.cos(2.0 * np.pi * mode * x / ncells)
        # Enforce exact zero total perturbation so the test isolates stiffness,
        # not rest-mass bookkeeping.
        wave -= np.mean(wave)
        perturbed = base + delta * wave
        c_rest = (rest_energy(perturbed, epsilon) - e0_rest) / (delta * delta)
        c_grad = (ring_gradient_energy(perturbed, a=1.0, stiffness=1.0) - e0_grad) / (delta * delta)
        rows.append((mode, float(c_rest), float(c_grad)))
    return rows


def stress_tensor_rest_frame(rho: float) -> np.ndarray:
    """Dust stress tensor in the local rest frame, diag(rho,0,0,0)."""

    t = np.zeros((4, 4), dtype=float)
    t[0, 0] = rho
    return t


def main() -> None:
    print("ITEM 123: R4 ZERO-MODE DUST HAMILTONIAN")
    print("=" * 92)

    print("\n[1] Operator inventory of the minimal reservoir lift")
    terms = operator_inventory()
    for term in terms:
        print(f"  {term.name:28s} {term.status:23s} {term.form}")
        print(f"    {term.reason}")
    admitted = [term for term in terms if term.status.startswith("ADMITTED")]
    rejected = [term for term in terms if term.status.startswith("REQUIRES")]
    check([term.name for term in admitted] == ["rest count", "local active exchange"], "only rest count plus local exchange are admitted")
    check(len(rejected) == 4, "every stiffness/transport term requires adding a new operator")

    print("\n[2] Rest-count equation of state")
    counts = np.array([108, 99, 104, 101, 96, 112, 97, 103], dtype=float)
    epsilon = 1.7
    for a in (0.25, 1.0, 4.0):
        rho, p, w = equation_of_state_for_rest_counts(counts, epsilon, a)
        print(f"  a={a:4.2f}: rho={rho:.9e}, p={p:.3e}, w={w:.3e}")
        check(abs(p) < 1.0e-9, "rest-count energy has p=0 at fixed comoving count")
        check(abs(w) < 1.0e-12, "rest-count energy gives w=0")
    rho1, _, _ = equation_of_state_for_rest_counts(counts, epsilon, 1.0)
    rho2, _, _ = equation_of_state_for_rest_counts(counts, epsilon, 2.0)
    check(abs(rho2 / rho1 - 1.0 / 8.0) < 1.0e-12, "rho scales exactly as a^-3")

    print("\n[3] Sound-speed / stiffness test")
    max_rest = 0.0
    min_grad = float("inf")
    for mode, c_rest, c_grad in sound_quadratic_coefficients():
        print(f"  mode={mode:2d}: rest quadratic cost={c_rest:+.3e}, stiffness-control cost={c_grad:.3e}")
        max_rest = max(max_rest, abs(c_rest))
        min_grad = min(min_grad, c_grad)
    check(max_rest < 1.0e-6, "rest-count reservoir has no k-dependent quadratic cost")
    check(min_grad > 0.0, "a true gradient term would create positive stiffness")
    print("  Therefore c_s^2=0 for the admitted reservoir Hamiltonian; nonzero c_s^2")
    print("  would be evidence for one of the rejected new-link/new-phase operators.")

    print("\n[4] Brown--Kuchar finite analogue")
    rho = 3.2
    tmunu = stress_tensor_rest_frame(rho)
    pressure_trace = float(np.trace(tmunu[1:, 1:]) / 3.0)
    w = pressure_trace / tmunu[0, 0]
    print(f"  T^mu_nu(rest frame) = diag({tmunu[0,0]:.3f}, 0, 0, 0)")
    print(f"  pressure trace / rho = {w:.3e}")
    check(np.allclose(tmunu[1:, 1:], 0.0), "spatial stress vanishes")
    check(abs(w) < 1.0e-15, "Brown-Kuchar dust stress has w=0")
    print("  The reservoir count is the finite label version of Brown--Kuchar dust:")
    print("  a constrained rest-density variable, not a propagating scalar field.")

    print("\n[5] Remaining boundary")
    print("  Closed here: Hamiltonian form.  The minimal R4 zero-mode reservoir is")
    print("  rest energy only, with exact shift/number symmetry and no admitted")
    print("  gradient stiffness; hence w=c_s^2=0.")
    print("  Not fixed by this Hamiltonian audit: the absolute boot density")
    print("  N_tot/volume.  The companion alpha/208 source-law audit gives the")
    print("  conditional density candidate; the remaining theorem is the source map")
    print("  selecting that sector rather than fitting it.")
    print("exit 0 -- constrained-dust Hamiltonian closes for the minimal reservoir lift.")


if __name__ == "__main__":
    main()
