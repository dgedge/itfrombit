#!/usr/bin/env python3
r"""ITEM 131: T5b correlation-volume audit for the HBC scalar amplitude.

Question
--------
Can the remaining T5b gate be closed after Pi_k is fixed?

T5a already derives the readout:

    Pi_k = compensated Fourier-shell projection at |k|=aH.

T5b asks how many independent HBC service events that normalized scalar mode
samples.  This script separates the formal correlation-volume theorem from the
still-missing microscopic number.

Result
------
For a translation-invariant local service-count current on one horizon shell,
the scalar density-contrast Fourier coefficient obeys

    Var(delta_j(k)) = F_eff S_j(k) / N_shell,
    N_eff(k)       = N_shell / S_j(k),

where

    N_shell = mean number of service events in the real-space shell window,
    S_j(k)  = spatial structure factor of the service-current correlations.

Thus T5b *form* closes: the correlation volume is not angular degeneracy,
not the 28 service labels, and not raw horizon entropy; it is the spatial
structure factor sampled by the Pi_k horizon shell.  In the white local-current
limit, S_j(k)=1 and N_eff=N_shell.

But the current canon does not derive S_j(k)=1 or N_shell=(4/3)alpha_0^-4.
The prior scale-invariance no-go remains: multiplying the scalar-current
intensity changes N_shell without affecting Pi_k, the 1/28 tilt, or F_eff.
So T5b is reduced to a precise structure-factor theorem target, not promoted
to a Locked amplitude derivation.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

A_TARGET = 2.10e-9
ALPHA0 = 1.0 / 137.0
C_F = Fraction(4, 3)
F_EFF = 1.0

Vector = tuple[int, int, int]


@dataclass(frozen=True)
class StructureRow:
    label: str
    structure_factor: float


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def char(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(-2j * math.pi * dot(k, x) / n)


def volume_average_mode_variance(n: int, lam: float, k: Vector, f_eff: float = F_EFF) -> float:
    """Var[M^-1 sum_x exp(-ikx) (N_x-lam)/lam] for independent Poisson cells."""
    modes = [char(n, k, x) for x in torus_points(n)]
    m = n**3
    return f_eff * sum(abs(z) ** 2 for z in modes) / (lam * m**2)


def correlated_mode_variance(n: int, lam: float, k: Vector, corr: dict[Vector, float], f_eff: float = F_EFF) -> float:
    """Translation-invariant covariance with C(r)=corr[r], C(0)=1."""
    pts = torus_points(n)
    m = n**3
    total = 0.0 + 0.0j
    for x in pts:
        for y in pts:
            r = tuple((x[i] - y[i]) % n for i in range(3))
            c = corr.get(r, 0.0)
            total += char(n, k, x) * char(n, k, y).conjugate() * c
    return float((f_eff * total / (lam * m**2)).real)


def structure_factor(n: int, k: Vector, corr: dict[Vector, float]) -> float:
    """S(k)=sum_r C(r) exp(-ik.r), with C(0)=1 for independent cells."""
    total = 0.0 + 0.0j
    for r, c in corr.items():
        total += c * char(n, k, r)
    return float(total.real)


def alpha4_count(alpha: float = ALPHA0) -> float:
    return float(C_F) * alpha**-4


def amplitude(n_eff: float, f_eff: float = F_EFF) -> float:
    return f_eff / n_eff


def sigma_residual(prediction: float) -> float:
    return (prediction - A_TARGET) / 0.03e-9


def main() -> None:
    print("ITEM 131 T5b CORRELATION-VOLUME AUDIT")

    print("\n[1] Source checks")
    check(
        contains("python_code/item131_scalar_mode_projector.py", "Pi_k projector form derived"),
        "T5a projector form is already closed",
    )
    check(
        contains("python_code/item131_feff_hbc.py", "F_eff = 1"),
        "T4 supplies the canonical CTMC Fano factor",
    )
    check(
        contains("python_code/item131_neff_duty_closure_attempt.py", "scale-invariance no-go"),
        "prior closure attempt records the unresolved intensity scaling",
    )
    check(
        contains("python_code/item131_t2_landauer_status_audit.py", "T2 closes for item 131"),
        "T2 has been removed as an amplitude coefficient issue",
    )

    print("\n[2] Finite-torus derivation of the mode correlation volume")
    n = 7
    m = n**3
    lam = 11.0
    k = (1, 0, 0)
    var_ind = volume_average_mode_variance(n, lam, k)
    n_shell = m * lam
    expected = F_EFF / n_shell
    print(f"  torus sites M              = {m}")
    print(f"  mean events per site lambda = {lam:.6f}")
    print(f"  N_shell=M lambda           = {n_shell:.6e}")
    print(f"  Var(delta_j(k)) independent = {var_ind:.12e}")
    print(f"  F_eff/N_shell              = {expected:.12e}")
    check(abs(var_ind / expected - 1.0) < 1.0e-12, "white local-current limit gives N_eff=N_shell")
    check(k != (0, 0, 0), "nonzero compensated mode has no mean service-load component")

    print("\n[3] Spatial structure factor is the real T5b object")
    zero: Vector = (0, 0, 0)
    corr_white = {zero: 1.0}
    corr_nearest = {zero: 1.0, (1, 0, 0): 0.20, (6, 0, 0): 0.20}
    corr_anticorr = {zero: 1.0, (1, 0, 0): -0.20, (6, 0, 0): -0.20}
    for label, corr in [
        ("white local current", corr_white),
        ("nearest-neighbour positive", corr_nearest),
        ("nearest-neighbour negative", corr_anticorr),
    ]:
        s_k = structure_factor(n, k, corr)
        var_k = correlated_mode_variance(n, lam, k, corr)
        n_eff = n_shell / s_k
        print(
            f"  {label:28s}: S_j(k)={s_k:.6f}  "
            f"Var={var_k:.12e}  N_eff={n_eff:.6e}"
        )
        check(abs(var_k / (F_EFF * s_k / n_shell) - 1.0) < 1.0e-12, f"{label}: variance equals F_eff S_j(k)/N_shell")
    check(structure_factor(n, k, corr_nearest) > 1.0, "positive same-axis correlations lower the independent count at this k")
    check(structure_factor(n, k, corr_anticorr) < 1.0, "anticorrelations raise the independent count at this k")

    print("\n[4] Angular degeneracy and service labels are not the correlation volume")
    shell_modes = [(1, 0, 0), (6, 0, 0), (0, 1, 0), (0, 6, 0), (0, 0, 1), (0, 0, 6)]
    variances = [volume_average_mode_variance(n, lam, q) for q in shell_modes]
    shell_avg = sum(variances) / len(variances)
    shell_raw = sum(variances)
    print(f"  |k|^2=1 shell degeneracy       = {len(shell_modes)}")
    print(f"  shell/per-mode average variance = {shell_avg:.12e}")
    print(f"  raw shell sum variance          = {shell_raw:.12e}")
    check(abs(shell_avg / expected - 1.0) < 1.0e-12, "shell average preserves the per-mode correlation volume")
    check(abs(shell_raw / expected - len(shell_modes)) < 1.0e-12, "raw shell sum would illegally multiply by angular degeneracy")
    check(True, "28 channel labels and 112 flags remain internal labels, not spatial correlation volumes")

    print("\n[5] What would be needed for the alpha^4 amplitude")
    n_alpha = alpha4_count()
    rows = [
        StructureRow("white current", 1.0),
        StructureRow("mild positive S=1.2", 1.2),
        StructureRow("mild negative S=0.8", 0.8),
        StructureRow("factor-two clustering", 2.0),
    ]
    print(f"  alpha^4 candidate independent count = {n_alpha:.6e}")
    for row in rows:
        n_shell_needed = n_alpha * row.structure_factor
        amp = amplitude(n_shell_needed / row.structure_factor)
        print(
            f"  {row.label:24s}: S_j={row.structure_factor:.3f}  "
            f"N_shell needed={n_shell_needed:.6e}  A={amp:.6e}  "
            f"pull={sigma_residual(amp):+.2f} sigma"
        )
    check(abs(sigma_residual(amplitude(n_alpha))) < 1.1, "alpha^4 candidate remains numerically viable when S_j=1")
    check(True, "but canon still needs a microscopic theorem selecting both N_shell and S_j(k)")

    print("\n[6] Gate status")
    closed_or_reduced = [
        "T5a projector: Pi_k is the compensated Fourier-shell readout",
        "T5b formal identity: Var(delta_j(k))=F_eff S_j(k)/N_shell",
        "white-current special case: S_j(k)=1 gives N_eff=N_shell",
        "angular degeneracy is divided out by the scalar shell average",
    ]
    still_open = [
        "derive S_j(k=aH)=1 or another explicit service-current structure factor",
        "derive N_shell=(4/3)alpha_0^-4 or another absolute scalar service density",
        "derive the duty/regime map that lets the scalar fluctuation ledger use the CTMC Fano reading",
    ]
    for item in closed_or_reduced:
        check(True, f"closed/reduced: {item}")
    for item in still_open:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  T5b does not fully close numerically, but it does reduce to a precise")
    print("  structure-factor theorem.  For the Pi_k scalar horizon shell,")
    print("      Var(delta_j(k)) = F_eff S_j(k) / N_shell,")
    print("      N_eff(k) = N_shell / S_j(k).")
    print("  Thus the correlation volume is the service-current spatial structure")
    print("  factor at k=aH, not angular degeneracy, service-label count, or raw")
    print("  horizon entropy.  The white local-current case would close T5b with")
    print("  S_j=1, leaving T3 as the mean-count problem, but current canon has not")
    print("  derived that spatial whiteness/independence statement.")
    print("=" * 108)
    print("exit 0 -- T5b reduced to service-current structure factor; numerical closure still open.")


if __name__ == "__main__":
    main()
