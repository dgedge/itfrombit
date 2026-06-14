#!/usr/bin/env python3
r"""SMG/TCH no-bulk-transition theorem skeleton.

This is a theorem-design artifact, not a canon update.

Question:
    What theorem would carry the certified SMG/TCH mirror-gap domain from the
    current polymer endpoint beta ~= 0.661 to the weak-coupling SU(3)
    continuum limit?

Answer:
    The needed theorem is not another finite-row diagonalisation.  It is a
    Wilsonian stability statement:

      certified endpoint + exact charge superselection + gapped mirror block
      + no gauge-sector bulk transition
      => mirror sector remains decoupled along the RG path to beta = infinity.

The mirror-induced-coupling clause has since been sharpened by
smg_induced_lambda_bound.py: inside the registered Hamiltonian, P_vac H P_ch=0
exactly, so the induced lambdas vanish.  The remaining non-finite input is the
ordinary pure/fundamental Wilson SU(3) no-bulk-transition statement.
"""

from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class Gate:
    name: str
    status: str
    evidence: str
    failure_mode: str


# From smg_certified_domain_extension.py, FP/independent-neighbourhood criterion.
BETA_KP_UNIFORM = 0.519101
BETA_KP_TWO_TYPE = 0.563360
BETA_DOBRUSHIN = 0.630301
BETA_FP = 0.661156

# From smg_continuum_decoupling_argument.py, continuous finite-cell envelope.
T_PHYSICAL = 1.0
T_CELL_STRESS = 3.846019
BETA_CELL_STRESS = 1.913777


def theorem_statement() -> str:
    return r"""
Theorem candidate (conditional RG decoupling).

Let H_TCH(beta,t) be the gauge-invariant TCH SU(3) Hamiltonian with the
CSS/SMG mirror block and true charged hopping, restricted to the Gauss-law
sector.  Let beta_* = 0.661156 be the strongest current polymer-certified
endpoint and let t <= 1 be the physical hopping normalization.

Assume:
  A. Endpoint certificate:
     At beta = beta_* the mirror block is in a unique gapped, symmetry-preserving
     phase with exponential clustering in the charged sector.

  B. Exact sector bookkeeping:
     The charged/mirror number used in the finite cell is an exact superselection
     block for the local TCH operators, so the vacuum cannot mix with a charged
     mirror sector except through a genuine gap closure.

  C. Gauge-path no-bulk-transition:
     The TCH Wilson SU(3) gauge path obtained by increasing beta from beta_* to
     infinity stays in one confining/asymptotically-free continuum basin.

  D. Zero mirror backreaction:
     For the registered electric, matter, magnetic, and charged-hopping
     operators, P_vac H P_ch=0 exactly.  Integrating out the massive mirror
     block therefore induces no adjoint or irrelevant pure-gauge deformation:
     beta_A=lambda_6=lambda_8=...=0.

  E. Gapped-sector stability:
     Along that gauge path, the mirror correlation length remains finite and
     the quasi-adiabatic continuation/stability bound applies uniformly.

Then the mirror block remains massive and symmetry-preserving for all beta >=
beta_*, and the beta -> infinity continuum limit contains the desired gauge
sector with the mirror block decoupled.  Any counterexample must be one of:
  (i) a gauge-sector bulk transition,
  (ii) a charged mirror pole becoming massless,
  (iii) addition of a new number-changing operator absent from the registered
        Hamiltonian, or
  (iv) failure of the endpoint certificate in the thermodynamic limit.
"""


def gates() -> list[Gate]:
    return [
        Gate(
            "polymer endpoint",
            "PROVED within registered criterion",
            (
                f"uniform KP beta={BETA_KP_UNIFORM:.6f}; two-type KP beta="
                f"{BETA_KP_TWO_TYPE:.6f}; Dobrushin beta={BETA_DOBRUSHIN:.6f}; "
                f"FP beta_*={BETA_FP:.6f}"
            ),
            "cheap polymer sharpening stalls before beta=1",
        ),
        Gate(
            "finite charge superselection",
            "PROVED in finite TCH cell",
            "commutators with electric, matter, magnetic, and hopping blocks vanish",
            "if a larger cell adds a noncommuting charged operator, the theorem changes",
        ),
        Gate(
            "finite stress margin",
            "MEASURED finite-cell support",
            (
                f"physical t={T_PHYSICAL:.1f}; first finite-cell stress crossing "
                f"t_c={T_CELL_STRESS:.6f} at beta={BETA_CELL_STRESS:.6f}; "
                f"margin={T_CELL_STRESS / T_PHYSICAL:.3f}x"
            ),
            "finite cell is evidence only; it is not a thermodynamic theorem",
        ),
        Gate(
            "gauge-path no-bulk-transition",
            "EXTERNAL lattice-gauge input",
            (
                "must accept or prove that the pure/fundamental Wilson SU(3) path "
                "from beta_* to infinity has no finite-beta bulk transition"
            ),
            "a pure-gauge bulk transition on the Wilson axis kills adiabatic continuation",
        ),
        Gate(
            "mirror-induced lambdas",
            "CLOSED for registered Hamiltonian",
            (
                "smg_induced_lambda_bound.py gives ||P0 O Pch||=0 for electric, "
                "matter, magnetic, and hopping operators; lambda_induced=0"
            ),
            "a future pair-creation term would need its own eta^2/Delta bound",
        ),
        Gate(
            "uniform mirror stability",
            "OPEN theorem",
            (
                "needs a thermodynamic lower bound on the charged mirror gap, or a "
                "quasi-adiabatic continuation bound using finite correlation length"
            ),
            "a charged mirror pole becoming massless is the actual phase-transition signal",
        ),
    ]


def induced_coupling_bound_shape() -> str:
    return r"""
Resolved induced-coupling target:

  S_eff[U] = S_Wilson[U; beta_R] + sum_i lambda_i(beta) O_i[U],
  lambda_i = 0

for the registered number-conserving TCH/SMG Hamiltonian, because every local
operator commutes with N_ch and P_vac O P_ch=0.  Hence the coupling-space path is

  (beta_F, beta_A, lambda_6, lambda_8, ...) = (beta_F, 0, 0, 0, ...).

If a new number-changing perturbation V_pair is later introduced, then the
conservative second-order budget is |lambda_pair| <= eta^2/Delta_ch, with
Delta_ch >= 3.960160 in the current finite audit.  Without such a new operator,
there is no SMG-specific induced-coupling coordinate to tune.
"""


def main() -> int:
    print("[0] SMG/TCH no-bulk-transition theorem skeleton")
    print(f"    certified endpoint beta_* = {BETA_FP:.6f}")
    print(f"    finite-cell physical hopping t = {T_PHYSICAL:.1f}")
    print(f"    finite-cell stress threshold t_c = {T_CELL_STRESS:.6f}")
    print(f"    stress margin t_c/t = {T_CELL_STRESS / T_PHYSICAL:.3f}x")
    assert 0.661 < BETA_FP < 0.662
    assert T_CELL_STRESS / T_PHYSICAL > 3.0
    print()

    print("[1] theorem candidate")
    print(theorem_statement().strip())
    print()

    print("[2] clause ledger")
    for gate in gates():
        print(f"    {gate.name}: {gate.status}")
        print(f"      evidence: {gate.evidence}")
        print(f"      failure:  {gate.failure_mode}")
    print()

    print("[3] induced-coupling reduction")
    print(induced_coupling_bound_shape().strip())
    print()

    print("[VERDICT]")
    print("  The SMG-specific induced-coupling problem is closed for the registered")
    print("  Hamiltonian: lambda_i=0 exactly. The remaining theorem is the ordinary")
    print("  pure/fundamental Wilson SU(3) no-finite-beta-bulk-transition input plus")
    print("  thermodynamic endpoint stability. A new pair-creation operator would reopen")
    print("  the lambda budget under the explicit eta^2/Delta rule.")
    print("exit 0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
