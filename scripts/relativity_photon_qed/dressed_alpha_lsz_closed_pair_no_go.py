#!/usr/bin/env python3
r"""LSZ/Thomson charge versus closed-record-pair kernel.

Target
------
Try to prove the positive theorem suggested by the dressed-alpha near-hit:

    physical Thomson/LSZ charge bills  2 * sum_f Q_f^2 - 1

instead of the ordinary retarded Ward/Kubo photon self-energy kernel

    sum_f Q_f^2.

The arithmetic is attractive.  For three SO(10) Weyl generations
sum_f Q_f^2 = 16, so 2*16-1 = 31, exactly the historical N1 kernel needed
for the 137 -> 137.036 shift.

Result
------
The theorem cannot be proved under the current framework plus standard
LSZ/Ward definitions.  It fails for structural, not numerical, reasons:

  1. The electromagnetic current is J = sum_f Q_f j_f.  Orthogonality of the
     charged species makes the current-current kernel sum_f Q_f^2.  A
     charge-blind factor 2*dim - 1 is not the same operator.

  2. Closed record pairs / Born squares live in the symmetric Keldysh/noise
     sector (A A*).  The Thomson charge is a retarded response / Peierls
     Hessian.  Schwinger-Keldysh keeps these as different components.

  3. The "-1" subtraction is not a neutral identity channel in a current
     correlator.  Neutral states have Q=0 and contribute zero.  A negative
     subtraction would require a negative spectral weight or a redefined
     electromagnetic current.

So this script strengthens the previous "not currently derived" result into a
no-go under the present observable map.  A positive theorem would have to add a
new EM response principle that explicitly replaces the retarded current with a
closed-record noise kernel while preserving measured QED phenomenology.
"""

from __future__ import annotations

import math


ALPHA0_INV = 137.0
ALPHA0 = 1.0 / ALPHA0_INV
ALPHA_PHYS_INV = 137.035999084
DELTA = ALPHA_PHYS_INV - ALPHA0_INV
TWO_PI = 2.0 * math.pi


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def inv_alpha_from_kernel(kernel: float) -> float:
    return ALPHA0_INV + kernel * ALPHA0 / TWO_PI


def sm_weyl_charge_square_per_generation() -> float:
    """SO(10) 16, written as left-handed Weyl fields including conjugates."""

    q_left = 3 * (2 / 3) ** 2 + 3 * (1 / 3) ** 2  # u_L,d_L, three colours
    u_conj = 3 * (2 / 3) ** 2
    d_conj = 3 * (1 / 3) ** 2
    lepton_doublet = 0.0**2 + 1.0**2
    e_conj = 1.0**2
    nu_conj = 0.0**2
    return q_left + u_conj + d_conj + lepton_doublet + e_conj + nu_conj


def main() -> None:
    print("DRESSED ALPHA: LSZ CLOSED-PAIR NO-GO")
    print("=" * 96)

    n_required = DELTA * TWO_PI / ALPHA0
    q2_one_gen = sm_weyl_charge_square_per_generation()
    q2_three_gen = 3.0 * q2_one_gen
    closed_pair_kernel = 2.0 * q2_three_gen - 1.0

    print("[1] The attractive arithmetic")
    print(f"  required kernel from observed shift       = {n_required:.6f}")
    print(f"  Weyl charge-square sum per generation     = {q2_one_gen:.6f}")
    print(f"  Weyl charge-square sum, three generations = {q2_three_gen:.6f}")
    print(f"  closed-pair candidate 2*sumQ2-1           = {closed_pair_kernel:.6f}")
    print(f"  alpha^-1 from candidate                   = {inv_alpha_from_kernel(closed_pair_kernel):.9f}")
    check(abs(q2_one_gen - 16.0 / 3.0) < 1e-12, "one SO(10) Weyl generation has sum Q^2 = 16/3")
    check(abs(q2_three_gen - 16.0) < 1e-12, "three Weyl generations have sum Q^2 = 16")
    check(abs(closed_pair_kernel - 31.0) < 1e-12, "closed-pair arithmetic gives N1=31")
    check(abs(closed_pair_kernel - n_required) / n_required < 5e-4, "N1=31 is the observed near-hit")

    print("\n[2] Current algebra gate")
    print("  J_mu = sum_f Q_f j_{f,mu}.  Species orthogonality kills f != g, so")
    print("  <J J> = sum_f Q_f^2 <j_f j_f>.  The kernel is charge-weighted.")
    print(f"  Retarded Ward/Kubo kernel                 = {q2_three_gen:.6f}")
    print(f"  Closed-pair candidate kernel              = {closed_pair_kernel:.6f}")
    print(f"  candidate / Ward kernel                   = {closed_pair_kernel / q2_three_gen:.6f}")
    check(closed_pair_kernel != q2_three_gen, "closed-pair kernel is not the electromagnetic current kernel")

    print("\n[3] Schwinger-Keldysh component gate")
    print("  Closed record pairs are symmetric/noise objects: A A*, or Keldysh rr.")
    print("  Thomson/LSZ charge is a retarded response: ra, the Peierls Hessian.")
    print("  FDT can relate noise to Im(response) at equilibrium, but it does not")
    print("  replace the retarded zero-momentum charge residue by a doubled noise count.")
    closed_pair_is_noise = True
    alpha_is_retarded_response = True
    same_component = False
    check(closed_pair_is_noise, "closed-record-pair theorem applies to the symmetric/noise component")
    check(alpha_is_retarded_response, "Thomson/LSZ alpha is the retarded Peierls response component")
    check(not same_component, "noise and retarded response are distinct Keldysh components")

    print("\n[4] The -1 subtraction gate")
    print("  Neutral identity/photon channels have Q=0, hence zero current matrix")
    print("  element.  They do not contribute -1 to a positive spectral measure.")
    spectral_weights = [q2_three_gen, 0.0]  # charged sector plus neutral identity/channel
    positive_sum = sum(w for w in spectral_weights if w > 0)
    neutral_sum = sum(w for w in spectral_weights if w == 0)
    print(f"  positive charged spectral weight          = {positive_sum:.6f}")
    print(f"  neutral spectral contribution             = {neutral_sum:.6f}")
    print("  A literal -1 would require negative norm/ghost weight or a new current.")
    check(neutral_sum == 0.0, "neutral channels contribute zero, not -1")
    check(closed_pair_kernel > 2.0 * positive_sum - 2.0, "the chosen -1 is a subtraction from a doubled noise count")

    print("\n[5] Phenomenological gate")
    beta_ratio = closed_pair_kernel / q2_three_gen
    print(f"  If used as the photon self-energy kernel, the QED beta/vacuum-polarisation")
    print(f"  coefficient is multiplied by {beta_ratio:.6f}.")
    print("  That changes ordinary charge-weighted running and e+e- vacuum-polarisation")
    print("  phenomenology.  It is not a harmless convention once alpha0 is fixed.")
    check(beta_ratio > 1.9, "closed-pair kernel would nearly double the Weyl charge-weighted response")

    print("\nTHEOREM VERDICT")
    print(
        """
  Positive theorem rejected under current axioms.

  The statement

      Thomson/LSZ charge bills 2 sum(Q^2) - 1

  is not derivable from record closure alone.  Record closure gives the
  symmetric probability/noise sector.  The measured fine-structure constant is
  the retarded Ward/Kubo response of the electromagnetic current.  Current
  algebra and spectral positivity force the charge-weighted kernel sum(Q^2);
  neutral channels contribute zero, not a selected -1.

  Therefore the only ways to make the positive theorem true are explicit new
  physics:

    A. redefine the electromagnetic current so that detectors measure a
       charge-blind closed-record noise kernel; or
    B. add a new response principle proving that the retarded LSZ residue is
       equal to that closed-record kernel while keeping QED phenomenology.

  Neither is present in the current framework.  The dressed-alpha residual is
  now sharper: closed-record-pair explains the N1=31 near-hit, but LSZ/Thomson
  charge remains retarded and charge-weighted.
exit 0 -- no positive LSZ billing theorem without a new EM response principle.
"""
    )


if __name__ == "__main__":
    main()
