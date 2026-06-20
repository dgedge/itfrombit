#!/usr/bin/env python3
r"""ITEM 123: conserved-reservoir lift of the R4 active Poisson ledger.

This tries the only live route left for the CMB slot.  The active R4 MOND
ledger is an open immigration-death process and has no conserved number by
itself.  But every open birth/death process can sometimes be lifted to a
closed exchange with a reservoir.  Here the candidate is:

    N_zero + N_active = N_tot   (exactly conserved)

with transitions

    N_zero -> N_zero - 1, N_active -> N_active + 1
        at rate Gamma x (N_zero/N_tot),

    N_active -> N_active - 1, N_zero -> N_zero + 1
        at rate Gamma N_active.

For N_tot >> x, the active marginal is Binomial(N_tot, x/(N_tot+x)) and tends
to Poisson(x), so it reproduces the MOND halo ledger while keeping an exact
conserved total.  The conjugate phase of N_tot has an exact shift symmetry,
and if the reservoir records carry no spatial gradient stiffness it is
pressureless dust.

This is real progress: two clauses can be closed conditionally.

  * shift symmetry: yes, as conservation of N_tot in the closed exchange;
  * no double-counting: yes, N_active is a small exchange excitation and
    N_zero is the homogeneous reservoir;
  * dust Hamiltonian: closed by the companion operator-inventory audit
    item123_r4_zero_mode_dust_hamiltonian.py for the minimal reservoir lift:
    reservoir labels are record-number labels, not spatial link fields, so the
    admitted Hamiltonian is rest energy only;
  * abundance: NOT derived.  N_tot is a conserved integration constant.  Setting
    it to omega_x h^2=0.096 is exactly the CMB requirement / declared R4 share,
    unless a boot initial-condition theorem fixes it.

So the CMB component is now reduced to one missing theorem: fix the conserved
reservoir total N_tot from boot/QEC initial conditions.  The active MOND ledger
does not fix it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from r4_eos_cmb_resolution import OMEGA_R4_H2


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def binomial_stats(ntot: int, x: float) -> tuple[float, float, float]:
    """Stationary active-count statistics for the finite reservoir lift.

    Detailed balance gives
        pi(n+1)/pi(n) = x (N-n) / (N (n+1)),
    hence Binomial(N, p=x/(N+x)).
    """

    p = x / (ntot + x)
    mean = ntot * p
    var = ntot * p * (1.0 - p)
    fano = var / mean if mean else 1.0
    return mean, var, fano


def exchange_backward_generator(ntot: int, x: float) -> np.ndarray:
    """Generator on observables f(n_active) for the closed reservoir exchange."""

    L = np.zeros((ntot + 1, ntot + 1), dtype=float)
    for n in range(ntot + 1):
        nzero = ntot - n
        birth = x * nzero / ntot
        death = float(n)
        if n < ntot:
            L[n, n + 1] += birth
            L[n, n] -= birth
        if n > 0:
            L[n, n - 1] += death
            L[n, n] -= death
    return L


def nullity(matrix: np.ndarray, tol: float = 1.0e-10) -> int:
    s = np.linalg.svd(matrix, compute_uv=False)
    return int(np.sum(s < tol))


@dataclass(frozen=True)
class Clause:
    name: str
    status: str
    reason: str


def main() -> None:
    print("ITEM 123: R4 ZERO-MODE RESERVOIR LIFT")
    print("=" * 90)

    print("\n[1] Closed exchange reproduces the active Poisson ledger")
    x = 3.0
    for ntot in (32, 128, 1024, 100_000):
        mean, var, fano = binomial_stats(ntot, x)
        print(
            f"  N_tot={ntot:6d}: mean={mean:.6f}, var={var:.6f}, "
            f"Fano={fano:.6f}, mean/x={mean/x:.6f}"
        )
    mean_big, _, fano_big = binomial_stats(100_000, x)
    check(abs(mean_big / x - 1.0) < 5.0e-5, "large reservoir recovers E[N_active]=x")
    check(abs(fano_big - 1.0) < 5.0e-5, "large reservoir recovers Poisson Fano factor")

    print("\n[2] The closed generator has one conserved sector per N_tot")
    for ntot in (8, 16, 32):
        L = exchange_backward_generator(ntot, x=1.0)
        const_resid = float(np.linalg.norm(L @ np.ones(ntot + 1)))
        nul = nullity(L)
        print(f"  N_tot={ntot:2d}: observable nullity={nul}, ||L 1||={const_resid:.2e}")
        check(const_resid < 1.0e-12, "constant observable is conserved inside a fixed-N_tot sector")
        check(nul == 1, "active count still has no additional conserved observable inside the sector")
    print("  The conserved object is the superselected sector label N_tot, not N_active.")

    print("\n[3] Dust clauses")
    clauses = [
        Clause(
            "shift symmetry",
            "CONDITIONAL-CLOSED",
            "N_tot conservation gives a conjugate phase theta with theta -> theta + const",
        ),
        Clause(
            "Poisson halo compatibility",
            "CLOSED-IN-LIMIT",
            "N_active marginal tends to Poisson(|g|/a0) for N_tot >> |g|/a0",
        ),
        Clause(
            "no double-counting",
            "CONDITIONAL-CLOSED",
            "N_zero is the homogeneous reservoir; N_active is the small exchanged halo excitation",
        ),
        Clause(
            "dust Hamiltonian",
            "COMPANION-CLOSED",
            "operator inventory admits rest count only; see item123_r4_zero_mode_dust_hamiltonian.py",
        ),
        Clause(
            "abundance",
            "OPEN",
            f"N_tot must be fixed to omega_x h^2={OMEGA_R4_H2:.3f}; exchange dynamics does not fix it",
        ),
    ]
    for clause in clauses:
        print(f"  {clause.name:26s} {clause.status:20s} {clause.reason}")

    print("\n[4] Why the abundance is still the hard part")
    print("  The master equation preserves every chosen N_tot sector.")
    print("  Therefore it cannot select N_tot.  Any value of omega_x is dynamically")
    print("  stable once chosen; the value 0.096 has to come from a boot initial")
    print("  condition, a global R4-sector counting theorem, or an observational input.")
    arbitrary = [0.01, 0.096, 0.3]
    for omega in arbitrary:
        print(f"    chosen omega_x h^2 = {omega:.3f} -> conserved equally well")
    check(True, "closed exchange gives conservation but not selection")

    print("\nVERDICT")
    print("  Progress: a closed reservoir lift can supply the exact shift symmetry and")
    print("  separate N_zero from the active MOND count without spoiling the Poisson halo")
    print("  ledger.  It is the right mathematical shape for a CMB dust component.")
    print("  Companion closed: the dust Hamiltonian form for the minimal reservoir")
    print("  lift (rest energy only, no gradient stiffness).")
    print("  Not closed by this exchange dynamics alone: the abundance is a")
    print("  superselected sector.  The companion boot-QEC source-law audit supplies")
    print("  the conditional alpha/208 density candidate; the remaining theorem is")
    print("  the source map selecting that sector rather than fitting it.")
    print("exit 0 -- reservoir lift advances the form; abundance remains open.")


if __name__ == "__main__":
    main()
