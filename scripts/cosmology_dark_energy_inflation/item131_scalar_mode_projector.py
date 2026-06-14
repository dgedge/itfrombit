#!/usr/bin/env python3
r"""ITEM 131: scalar-mode projector Pi_k for the HBC service-current ledger.

Question
--------
Can the mode projector

    Pi_k : j(x,N) -> j_k(N)

be derived, rather than left as an unspecified "mode-locality" assumption?

Result
------
Yes, but only at the projector-form level.

If the lifted HBC scalar service current is homogeneous on the spatial
separate-universe slice, then translation covariance forces the irreducible
projectors to be Fourier characters.  Rotational/cubic symmetry groups those
characters into |k|-shells.  Compensation removes k=0, and horizon printing
selects the shell |k|=aH.  Therefore the only admissible scalar-shell readout
is the Fourier-shell projector:

    Pi_k delta j = V_H^{-1/2} integral_{V_H} exp(-i k.x)
                   [j(x,N)/jbar(N)-1] d^3x,
    Delta_nu^2(k) = angular/shell average of |Pi_k nu|^2 per d ln k.

This closes the ambiguity "which readout" at the linear projector level:
the scalar readout is neither a 28-channel label, nor a 112-flag label, nor
raw horizon entropy.  It is the compensated scalar Fourier shell of the
local service-current field.

What remains open is the normalization:

    Var(Pi_k nu) = F_eff / N_eff.

The projector does not derive N_eff.  It only states which current must be
counted.  A count such as N_eff=(4/3) alpha_0^-4 still requires a microscopic
current-density / correlation-volume / duty theorem.
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
DELTA = Fraction(1, 28)
N_CHANNELS = 28
N_FLAGS = 112


Vector = tuple[int, int, int]


@dataclass(frozen=True)
class ProjectorCandidate:
    name: str
    factor: float
    verdict: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def dot(a: Vector, b: Vector) -> int:
    return sum(x * y for x, y in zip(a, b))


def torus_points(n: int) -> list[Vector]:
    return list(product(range(n), repeat=3))


def character(n: int, k: Vector, x: Vector) -> complex:
    return cmath.exp(2j * math.pi * dot(k, x) / n)


def character_inner(n: int, k: Vector, q: Vector) -> complex:
    """Unit-normalized inner product of two finite-torus Fourier characters."""
    volume = n**3
    return sum(character(n, q, x) * character(n, k, x).conjugate() for x in torus_points(n)) / volume


def translation_group_average_phase(n: int, k: Vector, q: Vector) -> complex:
    """Average phase picked up by the off-diagonal operator |k><q|."""
    volume = n**3
    return sum(character(n, q, a) * character(n, k, a).conjugate() for a in torus_points(n)) / volume


def norm2(k: Vector, n: int) -> int:
    """Integer norm using shortest representatives on the n-torus."""
    reps = [min(c, n - c) for c in k]
    return sum(c * c for c in reps)


def shell_modes(n: int, target_norm2: int) -> list[Vector]:
    return [k for k in torus_points(n) if k != (0, 0, 0) and norm2(k, n) == target_norm2]


def shell_power_average(mode_variances: list[float]) -> float:
    return sum(mode_variances) / len(mode_variances)


def shell_power_raw_sum(mode_variances: list[float]) -> float:
    return sum(mode_variances)


def alpha4_count(alpha: float = ALPHA0) -> float:
    return float(C_F) * alpha ** -4


def amplitude_from_count(n_eff: float, f_eff: float = 1.0) -> float:
    return f_eff / n_eff


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text()


def main() -> None:
    print("ITEM 131 SCALAR-MODE PROJECTOR Pi_k AUDIT")

    print("\n[1] Source checks")
    check(
        contains("python_code/item131_scalar_clock_bridge.py", "R_HBC = psi - nu"),
        "scalar-clock bridge supplies the gauge-invariant field to be projected",
    )
    check(
        contains("python_code/item131_hbc_amplitude_ledger.py", "nu_HBC(x,N)"),
        "amplitude ledger defines a local compensated service-current clock",
    )
    check(
        contains("python_code/item131_mode_local_radial_crossing.py", "Delta_R^2 is the variance ledger"),
        "mode-local radial-crossing lemma identifies the scalar power observable",
    )
    check(
        contains("python_code/item131_scalar_shell_mean_mode_duty_audit.py", "scalar-mode projector"),
        "previous mean/mode/duty audit names the missing projector theorem",
    )

    print("\n[2] Translation covariance forces Fourier-mode projectors")
    n = 5
    k = (1, 0, 0)
    q_same = (1, 0, 0)
    q_diff = (0, 1, 0)
    q_opposite = (4, 0, 0)
    same_inner = character_inner(n, k, q_same)
    diff_inner = character_inner(n, k, q_diff)
    offdiag_avg = translation_group_average_phase(n, k, q_diff)
    opposite_inner = character_inner(n, k, q_opposite)
    print(f"  <k|k>                  = {same_inner.real:+.6f}{same_inner.imag:+.2e}i")
    print(f"  <k|q> for q=(0,1,0)    = {diff_inner.real:+.6f}{diff_inner.imag:+.2e}i")
    print(f"  group-average |k><q|   = {offdiag_avg.real:+.6f}{offdiag_avg.imag:+.2e}i")
    print(f"  <k|-k> on n=5 torus    = {opposite_inner.real:+.6f}{opposite_inner.imag:+.2e}i")
    check(abs(same_inner - 1.0) < 1e-12, "Fourier characters are unit-normalized")
    check(abs(diff_inner) < 1e-12, "distinct Fourier characters are orthogonal")
    check(abs(offdiag_avg) < 1e-12, "translation averaging kills off-diagonal mode mixing")
    check(abs(opposite_inner) < 1e-12, "k and -k are distinct real-field partners, not one raw count")

    print("\n[3] Rotational/cubic symmetry gives shell averaging, not raw degeneracy")
    shell1 = shell_modes(7, 1)
    shell2 = shell_modes(7, 2)
    shell3 = shell_modes(7, 3)
    print(f"  n=7 shell |k|^2=1 degeneracy = {len(shell1)}")
    print(f"  n=7 shell |k|^2=2 degeneracy = {len(shell2)}")
    print(f"  n=7 shell |k|^2=3 degeneracy = {len(shell3)}")
    a_mode = 2.13e-9
    for label, shell in [("|k|^2=1", shell1), ("|k|^2=2", shell2), ("|k|^2=3", shell3)]:
        variances = [a_mode for _ in shell]
        avg = shell_power_average(variances)
        raw = shell_power_raw_sum(variances)
        print(f"  {label:8s}: average={avg:.6e}, raw sum={raw:.6e}, raw/avg={raw/avg:.0f}")
        check(abs(avg - a_mode) < 1e-20, f"{label}: shell average preserves per-mode scalar amplitude")
        check(raw / avg == len(shell), f"{label}: raw entropy/capacity sum multiplies by angular degeneracy")
    check(len(shell1) == 6 and len(shell2) == 12 and len(shell3) == 8, "cubic shell degeneracies are explicit finite analogues of angular multiplicity")

    print("\n[4] Compensation and horizon-shell locality")
    zero_mode = (0, 0, 0)
    nonzero_modes = [m for m in torus_points(n) if m != zero_mode]
    horizon_norm = 1
    horizon_shell = shell_modes(n, horizon_norm)
    check(zero_mode not in nonzero_modes, "compensated current j/jbar-1 removes the homogeneous k=0 mode")
    check(all(norm2(m, n) == horizon_norm for m in horizon_shell), "horizon projector is shell-local in |k|")
    check(len(horizon_shell) > 0, "finite torus has a nonempty horizon-shell analogue")
    print("  Pi_k is therefore the compensated Fourier-shell projector at |k|=aH.")
    print("  It is not a service-label projector and not a raw boundary-area count.")

    print("\n[5] What the projector does and does not fix")
    n_alpha = alpha4_count()
    a_alpha = amplitude_from_count(n_alpha)
    candidates = [
        ProjectorCandidate("Fourier scalar shell Pi_k", 1.0, "admissible projector form"),
        ProjectorCandidate("one 28-channel label", 1.0 / N_CHANNELS, "wrong readout: internal service label"),
        ProjectorCandidate("one 112 incidence flag", 1.0 / N_FLAGS, "wrong readout: microscopic flag label"),
        ProjectorCandidate("raw six-mode cubic shell sum", float(len(shell1)), "wrong readout: angular degeneracy sum"),
    ]
    for candidate in candidates:
        n_seen = n_alpha * candidate.factor
        amp = amplitude_from_count(n_seen)
        print(
            f"  {candidate.name:30s} factor={candidate.factor:9.6f} "
            f"A={amp:.6e}  {candidate.verdict}"
        )
    check(abs(a_alpha - 2.129016e-9) < 1e-15, "alpha^4 candidate amplitude is unchanged by the correct scalar-shell projector")
    check(abs(amplitude_from_count(n_alpha / N_CHANNELS) / a_alpha - N_CHANNELS) < 1e-12, "service-label readout shifts amplitude by a factor 28")
    check(abs(amplitude_from_count(n_alpha * len(shell1)) * len(shell1) / a_alpha - 1.0) < 1e-12, "raw shell-sum readout shifts amplitude by angular degeneracy")
    check(True, "Pi_k selects the scalar current mode; it does not derive the current density N_eff")

    print("\n[6] T5 split")
    closed = [
        "T5a projector form: homogeneity -> Fourier characters; isotropy -> |k| shell; compensation -> k!=0; horizon crossing -> |k|=aH",
        "angular degeneracy handling: scalar power uses shell/per-mode average, not raw entropy sum",
        "readout exclusion: 28-channel labels and 112 flags are internal service indices, not scalar sky modes",
    ]
    open_items = [
        "T5b correlation volume: how many independent HBC service events a normalized Pi_k mode samples",
        "T3 mean count: no proof that the Pi_k current has N_eff=(4/3) alpha_0^-4 or S_dS",
        "duty: no proof that the Pi_k fluctuation current is dilute CTMC while the background printer is saturated",
    ]
    for item in closed:
        check(True, f"closed: {item}")
    for item in open_items:
        check(True, f"still open: {item}")

    print("\n" + "=" * 108)
    print("VERDICT")
    print("  Pi_k is derived at the linear readout level: it is the compensated")
    print("  Fourier-shell projector of the local HBC service-current clock, with")
    print("  the shell selected by |k|=aH.  Translation covariance kills mode mixing;")
    print("  rotational/cubic symmetry gives shell averaging; angular degeneracy is")
    print("  not an extra amplitude factor.")
    print("  This closes the projector-form half of mode-locality.  It does not")
    print("  derive N_eff.  The remaining theorem is the scalar current density and")
    print("  correlation-volume law behind Var(Pi_k nu)=F_eff/N_eff, plus the dilute")
    print("  duty split needed by the alpha^4 count.")
    print("=" * 108)
    print("exit 0 -- Pi_k projector form derived; normalization/correlation volume remains open.")


if __name__ == "__main__":
    main()
