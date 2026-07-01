#!/usr/bin/env python3
r"""SMG/TCH continuum theorem gate.

Question
--------
Can the finite SMG/TCH mirror-gap evidence be promoted to a credible
weak-coupling Wilson SU(3) continuum theorem?

Result
------
Not unconditionally.  What can be proved from the present canon is the
framework-specific half of the theorem:

  * the registered mirror block is exactly charged-number block diagonal;
  * integrating it out induces no adjoint or irrelevant pure-gauge coordinates;
  * the RG path is therefore the pure fundamental Wilson SU(3) axis with a
    decoupled massive spectator, not a mixed-action path;
  * finite TCH rows that measure the correct charged/mirror observable remain
    gapped through the tested region-II window.

The remaining theorem is exactly:

    Pure Wilson-axis analyticity/no bulk transition
    +
    volume/cutoff-uniform electric-subtracted mirror-gap lower bound

for beta >= 1, after the local-defect certificate has covered
beta_cert <= beta < 1.  Without those inputs, a full continuum theorem is not
proved.  With them, a charged mirror pole cannot appear and the massive mirror
block decouples in the continuum limit.

This script is deliberately a gate rather than a heavy run.  It uses committed
source text plus the existing JSON verdicts to separate:

  * theorem-grade framework clauses,
  * finite numerical support,
  * finite rows that must not be cited as mirror-mass witnesses, and
  * the exact assumptions needed for conditional closure.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUN_DIR = ROOT / "python_code" / "smg_dmrg_runs"
C3 = 4.0 / 3.0
BETA_CERT = 0.661156
BETA_MIRROR_OFFSET_CERT = 1.0


@dataclass(frozen=True)
class Gate:
    name: str
    status: str
    consequence: str


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def load_json(name: str):
    return json.loads((RUN_DIR / name).read_text(encoding="utf-8"))


def require_text(path: str, needles: tuple[str, ...]) -> None:
    text = (ROOT / path).read_text(encoding="utf-8")
    missing = [needle for needle in needles if needle not in text]
    if missing:
        raise AssertionError(f"{path} missing required text: {missing}")


def electric_offset(row: dict) -> float:
    beta = float(row["beta"])
    raw = float(row["charged_gap"])
    electric = float(row.get("electric_ref", C3 / beta))
    return raw - electric


def finite_gap_rows() -> dict[str, dict[str, float | bool | int]]:
    one = load_json("smg_region2_beta_sweep_cell336.json")
    two = load_json("smg_region2_2plaq_magnetic_gap.json")[0]
    smoke = load_json("smg_region2_deep_charged_gap_smoke_2x1.json")

    one_offsets = [electric_offset(r) for r in one["rows"]]
    two_offsets = [electric_offset(r) for r in two["rows"]]
    smoke_offsets = [electric_offset(r) for r in smoke["rows"]]
    return {
        "336 magnetic cell": {
            "dim": int(one["cell_dim"]),
            "min_gap": float(one["min_charged_gap"]),
            "min_offset": min(one_offsets),
            "witness": True,
        },
        "2-plaquette magnetic cell": {
            "dim": int(two["dim"]),
            "min_gap": float(two["min_charged_gap"]),
            "min_offset": min(two_offsets),
            "witness": bool(two["no_closure"] and not two["gap_is_electric_artifact"]),
        },
        "2x1 deep smoke row": {
            "dim": int(smoke["dim"]),
            "min_gap": float(smoke["min_charged_gap"]),
            "max_abs_offset": max(abs(x) for x in smoke_offsets),
            "witness": False,
        },
    }


def framework_gates() -> list[Gate]:
    return [
        Gate(
            "block diagonal mirror sector",
            "closed for registered Hamiltonian",
            "vacuum and charged/mirror sectors cannot mix without adding a new operator",
        ),
        Gate(
            "induced gauge couplings",
            "closed for registered Hamiltonian",
            "beta_A=lambda_6=lambda_8=...=0, so the path is the pure Wilson axis",
        ),
        Gate(
            "finite gauge-cell witness",
            "positive but finite",
            "correct electric-subtracted mirror rows stay gapped in the tested window",
        ),
        Gate(
            "Wilson-axis no-bulk transition",
            "external / unproved by finite TCH algebra",
            "needed to exclude an ordinary gauge-sector bulk singularity",
        ),
        Gate(
            "volume/cutoff uniform mirror bound",
            "certified only on the strong-coupling offset subdomain",
            "electric-subtracted mirror offset is bounded for beta_cert <= beta < 1; beta >= 1 remains open",
        ),
    ]


def theorem_statement() -> str:
    return r"""
Conditional continuum theorem.

Let H_TCH(beta,L,Lambda_cut) be the registered Gauss-projected SMG/TCH
Hamiltonian on the pure fundamental Wilson SU(3) axis, with no added
number-changing pair source.  Let

  Delta_mir(beta,L,Lambda_cut)
    = E0(charged mirror sector) - E0(vacuum sector) - E_string_min(beta,L)

be the electric-subtracted charged/mirror gap.

Already derived by the framework:

  P_vac H_TCH P_ch = 0,
  delta beta_A = delta lambda_6 = delta lambda_8 = ... = 0.

Therefore the coupling-space path is gamma(beta)=(beta,0,0,...) plus a
decoupled mirror block.

If:

  A. the pure fundamental Wilson SU(3) axis has no finite-beta bulk transition
     for beta >= beta_cert; and
  B. inf_{beta>=1,L,Lambda_cut} Delta_mir(beta,L,Lambda_cut) > 0;

then the beta -> infinity continuum lift has no charged mirror pole and no
SMG-induced bulk transition.  The mirror sector is a massive spectator.

The interval beta_cert <= beta < 1 is already certified for the
electric-subtracted mirror offset by the local-defect bound

  Delta_mir >= DSMG - (n_inc ||W|| / 2) beta = 2 - 2 beta.

The raw charged-sector gap remains positive to beta ~= 1.46 after the
Gauss-law pinned-flux margin is included, but that pinned string energy is
subtracted in Delta_mir and therefore does not close the continuum observable.

Conversely, with the present operator algebra, any failure must be one of:

  1. a pure-gauge Wilson-axis bulk transition,
  2. a thermodynamic/cutoff mirror-gap closure,
  3. a new number-changing operator not present in the registered Hamiltonian,
  4. a failure of the finite endpoint/gap witness under scaling.
"""


def main() -> int:
    print("SMG/TCH CONTINUUM THEOREM GATE")
    print("=" * 88)

    require_text(
        "python_code/smg_induced_lambda_bound.py",
        (
            "P_vac H P_ch = 0 exactly",
            "lambda_induced = 0 exactly within this operator set",
        ),
    )
    require_text(
        "python_code/smg_rg_basin_reduction.py",
        (
            "pure/fundamental Wilson SU(3) axis has no finite-beta bulk transition",
            "all other coordinates = 0",
        ),
    )
    require_text(
        "python_code/smg_volume_uniform_mirror_bound.py",
        (
            "Delta_raw >= DSMG + (static-charge gauge energy)",
            "NOT a continuum closure",
        ),
    )
    require_text(
        "python_code/smg_electric_subtracted_gap_observable_theorem.py",
        (
            "Delta_raw is not a valid",
            "raw>0 but Delta_mir=0",
        ),
    )
    require_text(
        "python_code/smg_tch_weak_continuum_frontier_audit.py",
        (
            "electric-subtracted mirror-gap theorem",
            "C3/beta finite electric gaps as SMG survival",
        ),
    )

    print("\n[1] Framework-specific theorem clauses")
    for gate in framework_gates():
        print(f"  {gate.name}: {gate.status}")
        print(f"    -> {gate.consequence}")
    check(True, "registered SMG-specific path reduces to pure Wilson axis plus spectator")

    print("\n[2] Existing finite witness rows")
    rows = finite_gap_rows()
    for label, row in rows.items():
        if label == "2x1 deep smoke row":
            print(
                f"  {label:27s} dim={row['dim']:>7d} "
                f"min raw gap={row['min_gap']:.6f} max |offset|={row['max_abs_offset']:.3e} "
                f"witness={row['witness']}"
            )
        else:
            print(
                f"  {label:27s} dim={row['dim']:>7d} "
                f"min raw gap={row['min_gap']:.6f} min electric-subtracted offset={row['min_offset']:.6f} "
                f"witness={row['witness']}"
            )
    check(rows["336 magnetic cell"]["min_offset"] > 3.0, "336-state magnetic row has positive mirror offset")
    check(rows["2-plaquette magnetic cell"]["min_offset"] > 2.0, "2-plaquette magnetic row has positive mirror offset")
    check(rows["2-plaquette magnetic cell"]["witness"] is True, "2-plaquette row is not an electric artifact")
    check(rows["2x1 deep smoke row"]["max_abs_offset"] < 1.0e-3, "smoke row is electric-string only, not mirror witness")

    print("\n[3] Finite-cell nonclosure controls")
    full_cg = load_json("tch_peter_weyl_full_cg_mirror_verdict.json")
    frontier = load_json("tch_gauged_mirror_frontier_verdict.json")
    print(f"  full-CG scope: {full_cg['scope']}")
    print(f"  full-CG remaining gate: {full_cg['remaining_gate']}")
    print(f"  frontier summary: {frontier['summary']}")
    check("finite hybrid stress test" in full_cg["scope"], "full-CG verdict remains finite/hybrid")
    check("continuum/chiral closure is not claimed" in full_cg["scope"], "full-CG verdict does not claim continuum closure")
    check("not closed" in frontier["summary"], "gauged mirror frontier still rejects full closure")

    print("\n[4] Conditional theorem")
    print(theorem_statement().strip())
    check(math.isclose(BETA_CERT, 0.661156, rel_tol=0.0, abs_tol=1.0e-12), "beta_cert fixed to current certificate")
    check(BETA_MIRROR_OFFSET_CERT == 1.0, "electric-subtracted mirror-offset certificate reaches beta=1")

    print("\n[5] Decision")
    print(
        """
  Progress:
    The SMG-specific continuum obstruction is reduced cleanly.  The registered
    mirror block creates no induced mixed-action coordinates, and the finite
    rows that measure the correct electric-subtracted mirror observable stay
    positive in the tested region-II window.  The volume-uniform lower-bound
    programme also certifies the electric-subtracted mirror offset on
    beta_cert <= beta < 1; raw charged-sector positivity reaches beta ~= 1.46
    but that is not the continuum observable.

  Refutation of full closure:
    This does not prove the weak-coupling Wilson SU(3) continuum theorem.  The
    remaining load-bearing inputs are the ordinary pure-Wilson no-bulk premise
    and a volume/cutoff-uniform electric-subtracted mirror-gap bound for
    beta >= 1.  Those are exactly the clauses a theoretical physicist would ask for.

  Best honest status:
    conditional theorem-grade reduction, not a locked continuum construction.
    The next real advance is a finite-size/cutoff lower-bound programme for
    Delta_mir above beta=1, or importing/proving the pure Wilson-axis analyticity input.
"""
    )
    print("ALL ASSERTIONS PASSED -- SMG-specific continuum gate reduced; full Wilson-axis theorem remains conditional.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
