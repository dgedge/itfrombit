#!/usr/bin/env python3
r"""v_phase7_broken_triad_service_scale.py

Phase 7A: derive the broken-triad normal-ordering scale from the EW service
algebra.

The Phase-6 CW audit found that the finite Coleman-Weinberg normal-ordering
scale

    kappa = 1/sqrt(3)

lands the EW VEV at the 0.1% level in the current one-loop convention.  This
script tests whether that scale is an arbitrary fit or the RMS service scale
of the electroweak broken generator algebra.

We build the Higgs-doublet service representation for

    SU(2)_L x U(1)_Y -> U(1)_em.

Using the Higgs vacuum direction |v> = (0,1), the four real Lie-generator
service vectors are

    T1|v>, T2|v>, T3|v>, Y|v>.

The realified service map has rank 3.  Its null vector is the unbroken
electromagnetic generator Q = T3 + Y; the active quotient has exactly three
broken service directions (W1, W2, and Z-like T3-Y).  The realification matters:
T1|v> and T2|v> differ by a complex phase in the doublet Hilbert space, but
they are distinct real Lie-algebra service directions.  If the scalar pointer
writes a single radial record distributed uniformly over the active broken
service subspace, the per-service-leg normal-ordering scale is the RMS
projection

    h_leg = h / sqrt(rank_broken) = h / sqrt(3).

Exit 0 means kappa=1/sqrt(3) is derived under a precise and testable service
normalisation theorem:

    finite CW normal ordering bills one equal RMS unit per broken EW service
    direction, with the unbroken photon service direction excluded.

This is still not the full multi-loop precision calculation; it supplies the
missing substrate scale used by that calculation.
"""

from __future__ import annotations

import math
import sys

import numpy as np


TOL = 1.0e-10
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


def main() -> int:
    print("=" * 96)
    print("PHASE 7A: BROKEN-TRIAD SERVICE SCALE")
    print("=" * 96)

    sigma1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    sigma2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
    sigma3 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
    ident = np.eye(2, dtype=complex)
    generators = {
        "T1": 0.5 * sigma1,
        "T2": 0.5 * sigma2,
        "T3": 0.5 * sigma3,
        "Y": 0.5 * ident,
    }
    names = list(generators)
    vac = np.array([0.0, 1.0], dtype=complex)
    service_vectors = np.column_stack([generators[name] @ vac for name in names])
    # Gauge service directions are real Lie-algebra directions.  Realify the
    # complex doublet action before taking the rank; otherwise T1 and T2 are
    # falsely collapsed as a single complex line.
    service_real = np.vstack([service_vectors.real, service_vectors.imag])
    gram = service_real.T @ service_real
    evals, evecs = np.linalg.eigh(gram)

    print("\n[0] Higgs-doublet service Gram matrix")
    print("  generator order:", names)
    print(np.array2string(gram, precision=6, suppress_small=True))
    print("  eigenvalues:", [round(float(x), 12) for x in evals])
    rank = int(np.sum(evals > TOL))
    null = evecs[:, int(np.argmin(evals))]
    # Fix sign convention for display.
    if null[names.index("T3")] < 0:
        null = -null
    print("  null vector coefficients:", {name: round(float(null[i]), 6) for i, name in enumerate(names)})
    check("EW service Gram matrix has rank 3", rank == 3)
    check("exactly one unbroken null direction exists", np.sum(evals < TOL) == 1)
    check("null direction is proportional to T3 + Y",
          abs(abs(null[names.index("T3")]) - 1.0 / math.sqrt(2.0)) < 1e-9
          and abs(abs(null[names.index("Y")]) - 1.0 / math.sqrt(2.0)) < 1e-9
          and abs(null[names.index("T1")]) < 1e-9
          and abs(null[names.index("T2")]) < 1e-9
          and null[names.index("T3")] * null[names.index("Y")] > 0)

    print("\n[1] Broken service directions")
    q = (generators["T3"] + generators["Y"]) @ vac
    z = (generators["T3"] - generators["Y"]) @ vac
    w1 = generators["T1"] @ vac
    w2 = generators["T2"] @ vac
    print(f"  ||Q|v>||^2        = {np.vdot(q, q).real:.12f}")
    print(f"  ||(T3-Y)|v>||^2   = {np.vdot(z, z).real:.12f}")
    print(f"  ||T1|v>||^2       = {np.vdot(w1, w1).real:.12f}")
    print(f"  ||T2|v>||^2       = {np.vdot(w2, w2).real:.12f}")
    check("Q annihilates the Higgs vacuum", np.vdot(q, q).real < TOL)
    check("three broken directions have nonzero service vectors",
          all(np.vdot(x, x).real > TOL for x in (z, w1, w2)))

    print("\n[2] RMS normal-ordering scale")
    kappa = 1.0 / math.sqrt(rank)
    print(f"  broken service rank      = {rank}")
    print(f"  kappa = 1/sqrt(rank)     = {kappa:.12f}")
    print("  Interpretation: a single scalar pointer amplitude h is distributed over")
    print("  the three broken service directions.  The per-leg normal-ordering scale")
    print("  is therefore h/sqrt(3), while the unbroken photon direction is excluded.")
    check("EW broken-triad service scale is 1/sqrt(3)", abs(kappa - 1.0 / math.sqrt(3.0)) < 1e-15)

    print("\nVERDICT")
    print("  The kappa=1/sqrt(3) scale is derived from the EW service algebra, under")
    print("  the explicit normalisation rule that finite CW normal ordering bills one")
    print("  equal RMS unit per broken generator direction.  The theorem is local and")
    print("  algebraic: SU(2)_L x U(1)_Y has four service directions, the Higgs vacuum")
    print("  leaves Q=T3+Y unbroken, and the active service subspace has rank 3.")
    print("  What remains for a fully locked VEV is the external precision calculation:")
    print("  fixed-scheme threshold matching and multi-loop SM running.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
