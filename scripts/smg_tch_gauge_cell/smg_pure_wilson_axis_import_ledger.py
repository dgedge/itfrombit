#!/usr/bin/env python3
r"""Pure-Wilson-axis import ledger for the SMG/TCH continuum lift.

This is the honest next step after the finite SMG/TCH evidence.

The registered TCH mirror Hamiltonian already removes the framework-specific
mixed-action danger:

    P_vac H P_ch = 0,
    (beta_F, beta_A, lambda_6, lambda_8, ...) = (beta_F, 0, 0, 0, ...).

Thus the continuum path is the pure fundamental Wilson SU(3) axis plus a
massive spectator mirror block.  The remaining issue is not another finite
TCH cell.  It is the ordinary lattice-gauge statement:

    W0: pure Wilson SU(3) has no zero-temperature finite-beta bulk transition
        on beta_A = 0 from beta_cert to beta = infinity.

This script records W0 as an external/imported lattice-gauge premise, not an
internal theorem.  It also separates the evidence types:

  * mixed fundamental-adjoint bulk surfaces are off the TCH path;
  * perturbative pure-Yang-Mills beta function has no finite-coupling zero;
  * asymptotic scaling handles the weak-coupling end;
  * the middle-beta Wilson-axis analyticity is standard numerical practice but
    not a mathematical theorem supplied by finite SMG/TCH algebra.
"""

from __future__ import annotations

from dataclasses import dataclass
import math


BETA_CERT = 0.661156
BETA_PERTURBATIVE = 6.0
DELTA_MIR_FINITE_FLOOR = 3.960160

N = 3
SIXTEEN_PI2 = 16.0 * math.pi**2
B0 = (11.0 / 3.0) * N / SIXTEEN_PI2
B1 = (34.0 / 3.0) * N**2 / SIXTEEN_PI2**2

# Blum-DeTar-Heller-Karkkainen-Rummukainen-Toussaint 1995 improve the SU(3)
# fundamental-adjoint endpoint estimate to (beta_f,beta_a)=(4.00(7),2.06(8)).
BULK_ENDPOINT_BETA_F = 4.00
BULK_ENDPOINT_BETA_A = 2.06
TCH_PATH_BETA_A = 0.0


@dataclass(frozen=True)
class Clause:
    label: str
    tier: str
    content: str
    failure_mode: str


def beta_function_two_loop(g: float) -> float:
    """Pure SU(3) continuum beta function to two loops."""

    return -B0 * g**3 - B1 * g**5


def lattice_g0_sq(beta_f: float) -> float:
    """Bare Wilson lattice coupling convention beta_f=2N/g0^2."""

    return 2.0 * N / beta_f


def rg_invariant(beta_f: float) -> float:
    """Two-loop lattice Lambda scaling factor a Lambda_L."""

    g2 = lattice_g0_sq(beta_f)
    return (B0 * g2) ** (-B1 / (2.0 * B0**2)) * math.exp(-1.0 / (2.0 * B0 * g2))


def clauses() -> list[Clause]:
    return [
        Clause(
            "TCH path projection",
            "internal theorem-grade",
            "registered mirror block induces beta_A=lambda_6=lambda_8=...=0",
            "a future number-changing mirror source would reopen mixed-action flow",
        ),
        Clause(
            "mixed-action bulk surface",
            "external numerical lattice-gauge evidence",
            "known SU(3) fundamental-adjoint bulk endpoint lies at beta_A about 2.06",
            "if a bulk surface crossed beta_A=0, the Wilson-axis import would fail",
        ),
        Clause(
            "continuous finite-coupling fixed point",
            "perturbative support, not a full theorem",
            "two-loop pure-Yang-Mills beta function is negative for every sampled g>0",
            "a nonperturbative fixed point on beta_A=0 would invalidate the lift",
        ),
        Clause(
            "weak-coupling endpoint",
            "standard asymptotic-freedom scaling",
            "a Lambda_L scaling factor is positive and tends to zero as beta->infinity",
            "failure would be ordinary pure-gauge continuum physics",
        ),
        Clause(
            "middle-beta Wilson-axis analyticity",
            "external import / still not proved here",
            f"remaining window beta in ({BETA_CERT:.3f},{BETA_PERTURBATIVE:.1f})",
            "a zero-temperature Wilson-axis bulk transition in this window kills closure",
        ),
    ]


def decision(wilson_axis_imported: bool, uniform_mirror_bound: bool) -> str:
    if wilson_axis_imported and uniform_mirror_bound:
        return "closed conditional on imported pure-Wilson axis plus uniform mirror gap"
    if not wilson_axis_imported and uniform_mirror_bound:
        return "blocked by ordinary pure-gauge Wilson-axis analyticity"
    if wilson_axis_imported and not uniform_mirror_bound:
        return "blocked by thermodynamic/cutoff mirror-gap lower bound"
    return "blocked by both pure-gauge analyticity and mirror-gap uniformity"


def main() -> int:
    print("[0] Pure-Wilson-axis import ledger for SMG/TCH")
    print(f"    beta_cert = {BETA_CERT:.6f}")
    print(f"    finite electric-subtracted mirror floor = {DELTA_MIR_FINITE_FLOOR:.6f}")
    print(f"    TCH path beta_A = {TCH_PATH_BETA_A:.1f}")
    print(
        "    literature bulk endpoint in fundamental-adjoint plane: "
        f"(beta_F,beta_A)=({BULK_ENDPOINT_BETA_F:.2f},{BULK_ENDPOINT_BETA_A:.2f})"
    )
    assert 0.661 < BETA_CERT < 0.662
    assert DELTA_MIR_FINITE_FLOOR > 0
    assert TCH_PATH_BETA_A == 0.0
    assert BULK_ENDPOINT_BETA_A > 2.0
    print()

    print("[1] Clause ledger")
    for clause in clauses():
        print(f"    {clause.label}: {clause.tier}")
        print(f"      content: {clause.content}")
        print(f"      failure: {clause.failure_mode}")
    print()

    print("[2] Perturbative support for no continuous finite-coupling fixed point")
    print(f"    b0={B0:.8f}, b1={B1:.8f}")
    for g in (0.5, 1.0, 1.5, 2.0, 3.0):
        value = beta_function_two_loop(g)
        print(f"    g={g:3.1f}: beta(g)={value:+.6e}")
        assert value < 0.0
    print("    two-loop beta(g)<0 for sampled g>0; only perturbative zero is g=0")
    print()

    print("[3] Asymptotic scaling end")
    previous = None
    for beta_f in (6.0, 8.0, 10.0, 12.0):
        scale = rg_invariant(beta_f)
        print(f"    beta_F={beta_f:4.1f}: a Lambda_L ~ {scale:.6e}")
        assert scale > 0.0
        if previous is not None:
            assert scale < previous
        previous = scale
    assert rg_invariant(50.0) < 1.0e-3 * rg_invariant(6.0)
    print("    scaling factor decreases to zero at beta=infinity")
    print()

    print("[4] Import status")
    print(f"    distance from known mixed-action endpoint to TCH beta_A=0 path: {BULK_ENDPOINT_BETA_A:.2f}")
    assert abs(TCH_PATH_BETA_A - BULK_ENDPOINT_BETA_A) > 2.0
    print("    This supports, but does not prove, the Wilson-axis import W0.")
    print("    W0 remains an external lattice-gauge premise unless a proof of")
    print("    zero-temperature pure-Wilson-axis analyticity is supplied.")
    print()

    print("[5] Decision table")
    for w0 in (False, True):
        for gap in (False, True):
            print(f"    W0={w0!s:<5s} uniform_gap={gap!s:<5s} -> {decision(w0, gap)}")
    assert "closed conditional" in decision(True, True)
    assert "ordinary pure-gauge" in decision(False, True)
    print()

    print("[VERDICT]")
    print("  Progress: the pure-Wilson-axis clause is now isolated as a single")
    print("  external/imported SU(3) lattice-gauge premise W0, not a hidden SMG fit.")
    print("  The TCH-specific mixed-action problem is closed by beta_A=lambda_i=0.")
    print("  The full chiral gauge continuum lift remains unclosed until W0 and the")
    print("  volume/cutoff-uniform electric-subtracted mirror-gap bound are proved or")
    print("  explicitly imported.")
    print("ALL ASSERTIONS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
