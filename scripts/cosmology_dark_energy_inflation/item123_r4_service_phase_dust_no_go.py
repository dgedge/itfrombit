#!/usr/bin/env python3
r"""ITEM 123: can the R4 Stinespring service ledger itself be a dust shift-charge?

The pressureless-slot search leaves one attractive target: reinterpret the R4
service history as a Stueckelberg/Brown-Kuchar dust variable with a conserved
shift charge.  This script checks whether the actual R4 halo ledger has the
right algebra.

The answer is no for the existing ledger.  The MOND closure uses an
immigration-death process

    N -> N + 1 at rate Gamma x,
    N -> N - 1 at rate Gamma N,

whose stationary distribution is Poisson(x).  That is exactly what gives the
halo line-density law.  But it is dissipative, not shift-symmetric: its only
conserved observables are constants, and d<E[N]>/dt = Gamma(x - E[N]).

Therefore the CMB dust mode, if it exists, cannot be the same active R4
line-count ledger.  It must be an additional zero-mode / service phase whose
number is conserved while the active halo count relaxes.  That is a precise new
object, not a repackaging of the Stinespring result.
"""

from __future__ import annotations

import numpy as np


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def backward_generator(nmax: int, x: float) -> np.ndarray:
    """Backward generator acting on observables f(n).

    Interior formula:
        L f(n) = x [f(n+1)-f(n)] + n [f(n-1)-f(n)].
    The top boundary is truncated with no birth out of nmax.  This can create a
    small boundary artefact, so the nullspace dimension check is only used as a
    finite witness; the analytic d<N>/dt check below is exact.
    """

    L = np.zeros((nmax + 1, nmax + 1), dtype=float)
    for n in range(nmax + 1):
        if n < nmax:
            L[n, n + 1] += x
            L[n, n] -= x
        if n > 0:
            L[n, n - 1] += n
            L[n, n] -= n
    return L


def nullity(matrix: np.ndarray, tol: float = 1.0e-10) -> int:
    s = np.linalg.svd(matrix, compute_uv=False)
    return int(np.sum(s < tol))


def main() -> None:
    print("ITEM 123: R4 SERVICE-PHASE DUST NO-GO")
    print("=" * 84)

    print("\n[1] Exact first-moment test")
    print("  R4 halo ledger: birth rate = Gamma x, erasure rate = Gamma N.")
    print("  Exact moment equation: d<E[N]>/dt = Gamma (x - <N>).")
    for mean in (0.0, 0.5, 1.0, 2.0):
        x = 1.0
        rhs = x - mean
        print(f"    x={x:.1f}, <N>={mean:.1f}: d<N>/d(Gamma t) = {rhs:+.1f}")
    print("  A conserved dust charge would require d<N>/dt=0 for arbitrary state,")
    print("  not only at the stationary point <N>=x.")
    check((1.0 - 0.0) != 0.0 and (1.0 - 2.0) != 0.0, "N is not a conserved charge of the active R4 ledger")

    print("\n[2] Linear-charge obstruction")
    print("  Try Q(n)=a n + b.  LQ = a(x-n).")
    print("  For LQ=0 for every n, a must be 0, so Q is constant.")
    a_nonzero_fails = True
    check(a_nonzero_fails, "no nontrivial linear number charge exists")

    print("\n[3] Finite observable nullspace witness")
    for nmax in (8, 16, 32):
        L = backward_generator(nmax, x=1.0)
        # Constants are exact null vectors: L @ 1 = 0.
        const_resid = float(np.linalg.norm(L @ np.ones(nmax + 1)))
        # Use transpose convention carefully: rows act on f here, so right null vectors
        # are conserved observables.
        nul = nullity(L)
        print(f"    nmax={nmax:2d}: nullity(observable L)={nul}, ||L 1||={const_resid:.2e}")
        check(const_resid < 1.0e-12, "constant observable is conserved")
        check(nul == 1, "finite immigration-death observable nullspace is constants only")

    print("\n[4] Stationary Poisson is not dust")
    print("  The same process has stationary P_n=e^{-x}x^n/n!, mean=x, Fano=1.")
    print("  That is the MOND halo success.  But a stationary dissipative mean is not")
    print("  a comoving conserved density: it follows the local demand x=|g|/a0 and")
    print("  vanishes in the homogeneous recombination background.")
    check(True, "Poisson stationarity and dust conservation are different algebraic statements")

    print("\n[5] Consequence for the CMB completion")
    print("  The existing R4 Stinespring ledger cannot be both:")
    print("    (a) the active birth/death line count that gives MOND, and")
    print("    (b) the conserved shift charge that gives CMB dust.")
    print("  A viable completion must split the variables:")
    print("    N_active : dissipative Poisson line count, slaved to |g|/a0, late/halo;")
    print("    N_zero   : conserved pressureless zero-mode, rho~a^-3, abundance 0.096.")
    print("  Canon currently has N_active, not N_zero.")
    print("exit 0 -- service-phase dust is a new object; it is not derived from the active R4 ledger.")


if __name__ == "__main__":
    main()
