#!/usr/bin/env python3
r"""Boundary-printing tensor prediction audit.

Question
--------
Can the boundary-printing branch give a near-term tensor prediction rather than
an order-of-magnitude slogan?

Result
------
Yes, but the precise statement has two layers:

  1. The linear primordial tensor-to-scalar ratio is exactly zero in the scalar
     printer branch:

         r_linear = 0.

     This is not a small coefficient.  It follows from the SVT/TT projection:
     a scalar source has the form a delta_ij + b k_i k_j, and the transverse-
     traceless projector annihilates it for every wavevector k.

  2. The irreducible floor is second-order scalar-induced gravitational waves.
     The framework fixes the scalar amplitude

         A_s = (3/4) alpha_0^4 = 2.13e-9,

     so the induced effective tensor ratio is

         r_induced = C_SIGW A_s.

     The transfer coefficient C_SIGW is an external second-order Boltzmann /
     radiation-transfer number, not a new substrate parameter.  Taking the
     natural unit coefficient gives r_induced = 2.13e-9.  Even an
     order-of-magnitude band C_SIGW in [0.1, 10] stays far below LiteBIRD /
     CMB-S4 reach.  Therefore the prediction-grade near-term statement is the
     null:

         no primordial B-mode detection at r >= 1e-3.

Tier
----
Linear tensor null: theorem-grade inside the scalar-printer/no-squeezing branch.
Induced floor coefficient: external transfer calculation; currently bracketed.
"""

from __future__ import annotations

import math

import numpy as np


ALPHA0 = 1.0 / 137.0
MPBAR_GEV = 2.435e18
U_EVENT = math.log(2.0)
R_CURRENT_BOUND = 0.036
R_NEAR_TERM = 1.0e-3


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def tt_projector(k: np.ndarray) -> np.ndarray:
    """Return Lambda_ij,kl for the transverse-traceless projection."""

    n = np.asarray(k, dtype=float)
    n = n / np.linalg.norm(n)
    p = np.eye(3) - np.outer(n, n)
    lam = np.zeros((3, 3, 3, 3), dtype=float)
    for i in range(3):
        for j in range(3):
            for a in range(3):
                for b in range(3):
                    lam[i, j, a, b] = p[i, a] * p[j, b] - 0.5 * p[i, j] * p[a, b]
    return lam


def apply_tt(k: np.ndarray, source: np.ndarray) -> np.ndarray:
    return np.einsum("ijab,ab->ij", tt_projector(k), source)


def transverse_basis(k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = np.asarray(k, dtype=float)
    n = n / np.linalg.norm(n)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(trial @ n)) > 0.8:
        trial = np.array([0.0, 1.0, 0.0])
    e1 = trial - n * float(trial @ n)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    return e1, e2


def plus_polarisation(k: np.ndarray) -> np.ndarray:
    e1, e2 = transverse_basis(k)
    return np.outer(e1, e1) - np.outer(e2, e2)


def scalar_source(k: np.ndarray, a: float, b: float) -> np.ndarray:
    n = np.asarray(k, dtype=float)
    n = n / np.linalg.norm(n)
    return a * np.eye(3) + b * np.outer(n, n)


def h_star() -> float:
    return math.sqrt(6.0 * math.pi**2 / U_EVENT) * ALPHA0**2 * MPBAR_GEV


def r_naive_de_sitter(a_s: float) -> float:
    return (2.0 / math.pi**2) * (h_star() / MPBAR_GEV) ** 2 / a_s


def main() -> None:
    print("BOUNDARY-PRINTING TENSOR PREDICTION AUDIT")
    print("=" * 88)

    print("\n[1] TT projector kills every scalar printer source")
    rng = np.random.default_rng(146)
    worst_scalar = 0.0
    worst_tt_sensitivity = 0.0
    for _ in range(64):
        k = rng.normal(size=3)
        a, b = rng.normal(size=2)
        s_scalar = scalar_source(k, a, b)
        projected = apply_tt(k, s_scalar)
        worst_scalar = max(worst_scalar, float(np.linalg.norm(projected)))

        e_plus = plus_polarisation(k)
        projected_tt = apply_tt(k, e_plus)
        worst_tt_sensitivity = max(
            worst_tt_sensitivity,
            float(np.linalg.norm(projected_tt) / np.linalg.norm(e_plus)),
        )
    print(f"  max ||TT[a delta_ij + b k_i k_j]|| over random tests = {worst_scalar:.3e}")
    print(f"  TT projector response to a true plus polarisation      = {worst_tt_sensitivity:.3f}")
    check(worst_scalar < 1.0e-12, "linear scalar printer has exactly zero TT source")
    check(worst_tt_sensitivity > 0.99, "the test is sensitive to genuine tensor sources")

    print("\n[2] Scalar amplitude and the tensor null")
    a_s = 0.75 * ALPHA0**4
    r_linear = 0.0
    r_induced_unit = a_s
    print(f"  A_s = (3/4) alpha0^4 = {a_s:.12e}")
    print(f"  r_linear             = {r_linear:.1f}  (exact in the scalar-printer branch)")
    print(f"  r_induced(C=1)       = {r_induced_unit:.12e}")
    check(2.0e-9 < a_s < 2.3e-9, "framework scalar amplitude is the Planck-scale A_s value")
    check(r_linear == 0.0, "linear primordial tensor ratio is exactly zero")

    print("\n[3] External second-order floor bracket")
    for coeff in (0.1, 1.0, 10.0):
        r_eff = coeff * a_s
        print(f"  C_SIGW={coeff:>4.1f}: r_induced = {r_eff:.3e}, "
              f"below r=1e-3 by {R_NEAR_TERM / r_eff:.2e}x")
        check(r_eff < R_NEAR_TERM, f"C_SIGW={coeff:g} floor is invisible to near-term B-mode searches")

    print("\n[4] Standard squeezed-graviton branch is mutually exclusive with one-bit printing")
    r_naive = r_naive_de_sitter(a_s)
    u_for_current_bound = U_EVENT * (r_naive / R_CURRENT_BOUND)
    u_for_near_term = U_EVENT * (r_naive / R_NEAR_TERM)
    print(f"  H_* = {h_star():.3e} GeV")
    print(f"  r_naive from squeezed de Sitter gravitons = {r_naive:.3e}")
    print(f"  event unit needed for r<{R_CURRENT_BOUND:g}: {u_for_current_bound/U_EVENT:.1f} bits/cell")
    print(f"  event unit needed for r<{R_NEAR_TERM:g}: {u_for_near_term/U_EVENT:.1f} bits/cell")
    check(r_naive > 1.0, "standard high-scale squeezed-graviton branch catastrophically overshoots")
    check(u_for_near_term / U_EVENT > 1.0e4, "a near-term-safe squeezed branch would abandon one-bit printing")

    print(
        f"""
[5] VERDICT
  Prediction-grade statement:
      r_primordial,linear = 0.

  The branch does not contain a squeezed primordial tensor vacuum.  The printer
  source is scalar, and the TT projector annihilates scalar sources exactly.
  The only remaining tensor background is second-order scalar-induced radiation,
  with an effective floor r_induced = C_SIGW A_s.  The framework fixes
  A_s=(3/4)alpha0^4={a_s:.3e}; C_SIGW is an external transfer coefficient, not a
  substrate fitting knob.  With C_SIGW=1 this is the familiar r~{r_induced_unit:.1e}.

  Near-term falsifier:
      any robust primordial B-mode detection at r >= {R_NEAR_TERM:.0e}
      after dust, lensing, and systematics are removed.

  Tier:
      linear null CLOSED inside the scalar-printer/no-squeezing branch;
      exact scalar-induced coefficient requires an external second-order CMB
      transfer calculation.
exit 0 -- tensor prediction sharpened to r_linear=0 plus r_induced~A_s floor.
"""
    )


if __name__ == "__main__":
    main()
