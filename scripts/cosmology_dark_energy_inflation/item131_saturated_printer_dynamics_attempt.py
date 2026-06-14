#!/usr/bin/env python3
"""ITEM 131: saturated-printer dynamics derivation attempt.

Target
------
Derive the saturation identity needed by the Item-131 tilt/amplitude cleanup:

    lambda_shell = N_shell alpha0^4 = C_F,      C_F = 4/3.

What existing canon already supplies:

  * HBC: boundary printing supplies zero-entropy capacity and keeps the
    Bekenstein-Hawking entropy bound saturated against incoming Landauer flux.
  * QEC: the post-decoder topology-changing current is weight-4, so its
    event probability is alpha0^4.
  * C_F audit: the confinement-unvetoed colour-silent logicals carry the
    colour-restoring load C_F = 4/3.

Result
------
The current mechanics derive the *admissible capacity inequality*

    lambda_shell <= C_F,

and they derive the constrained-optimization corollary

    if the early printer maximizes coherent cell production subject to that
    QEC capacity, then lambda_shell = C_F.

But the maximizing/backlog-latch clause is not derived by the existing HBC/QEC
mechanics.  HBC saturates the entropy-area ledger; it does not by itself state
that the topology-changing colour-load ledger must run at its QEC capacity.
Subcritical colour-load printers are consistent with the currently stated
mechanics and simply predict a different scalar amplitude.

Therefore the requested saturation theorem is not Locked.  It is reduced to one
new dynamical axiom/theorem: maximal coherent boundary printing, or equivalently
a backlog feedback law that drives the topology-changing current to the
colour-restoring capacity during the finite early phase.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ALPHA0 = 1.0 / 137.0
C_F = 4.0 / 3.0


def check(condition: bool, message: str) -> None:
    print(f"  [{'PASS' if condition else 'FAIL'}] {message}")
    if not condition:
        raise AssertionError(message)


def contains(path: str, phrase: str) -> bool:
    return phrase in (ROOT / path).read_text(encoding="utf-8")


def amp_from_utilization(utilization: float) -> float:
    """A_nu if lambda_shell = utilization * C_F and F_eff*S_j=1."""
    n_shell = utilization * C_F / ALPHA0**4
    return 1.0 / n_shell


def hstar_from_utilization(utilization: float, u_event: float = math.log(2.0)) -> float:
    """H_* from N_shell = S_dS/u_event with reduced M_P in GeV."""
    mbar = 2.435e18
    n_shell = utilization * C_F / ALPHA0**4
    return math.sqrt(8.0 * math.pi**2 * mbar**2 / (u_event * n_shell))


def mm1_mean_backlog(rho: float) -> float:
    """Mean queue length for a stable M/M/1 queue."""
    if rho >= 1.0:
        return math.inf
    return rho / (1.0 - rho)


@dataclass(frozen=True)
class UtilizationCase:
    rho: float
    label: str
    status: str


def main() -> None:
    print("ITEM 131 SATURATED-PRINTER DYNAMICS DERIVATION ATTEMPT")

    print("\n[1] Source gates")
    check(
        contains("cosmological_qec_engine/cosmological_qec_engine.tex", "must furnish, at every instant, sufficient zero-entropy capacity"),
        "HBC requires sufficient zero-entropy capacity at every instant",
    )
    check(
        contains("cosmological_qec_engine/cosmological_qec_engine.tex", "Bekenstein--Hawking entropy bound saturation against the incoming Landauer flux"),
        "HBC identifies boundary printing with entropy-bound saturation",
    )
    check(
        contains("cosmological_qec_engine/cosmological_qec_engine.tex", "raising $S_{\\max}$ by exactly $\\ln 2$ per cell"),
        "HBC supplies the printed-cell entropy unit ln2",
    )
    check(
        contains("python_code/item131_r2r3_gauge_filtered_threshold_route.py", "The topology-changing sector is weight-4"),
        "gauge-filtered route identifies topology-changing current as weight-4",
    )
    check(
        contains("python_code/item131_cf_stop_rule_closure_attempt.py", "2*(8/12)=4/3"),
        "C_F audit selects the colour-restoring load 4/3",
    )
    check(
        contains("python_code/item131_saturation_residuals_audit.py", "derive saturated printer dynamics")
        and contains("python_code/item131_saturation_residuals_audit.py", "Lock both at once"),
        "cluster cleanup names saturated-printer dynamics as the remaining Lock target",
    )

    print("\n[2] What follows from QEC capacity")
    n_max = C_F / ALPHA0**4
    a_sat = amp_from_utilization(1.0)
    h_sat = hstar_from_utilization(1.0)
    print("  Define lambda_shell = N_shell alpha0^4.")
    print("  Coherent colour-restoring sustainability gives:")
    print(f"      lambda_shell <= C_F = {C_F:.12f}")
    print(f"      N_shell <= C_F alpha0^-4 = {n_max:.6e}")
    print("  If a separate maximal-throughput principle is supplied, the optimum is")
    print("  the boundary point:")
    print(f"      lambda_shell = C_F, A_nu = {a_sat:.6e}, H_* = {h_sat:.6e} GeV")
    check(abs(a_sat - 0.75 * ALPHA0**4) / a_sat < 1e-15, "boundary optimum gives A_nu=(3/4)alpha0^4")

    print("\n[3] Countermodels inside the current mechanics")
    cases = [
        UtilizationCase(0.25, "quarter-load printer", "sustainable but not amplitude target"),
        UtilizationCase(0.50, "half-load printer", "sustainable but different A_nu/H_*"),
        UtilizationCase(0.90, "near-capacity printer", "sustainable with finite queue backlog"),
        UtilizationCase(1.00, "critical capacity printer", "target equality, queue-critical"),
        UtilizationCase(1.10, "overloaded printer", "incoherent / exits or accumulates uncorrected load"),
    ]
    for case in cases:
        if case.rho <= 1.0:
            amp = amp_from_utilization(case.rho)
            h = hstar_from_utilization(case.rho)
            backlog = mm1_mean_backlog(case.rho)
            backlog_text = "infinite" if math.isinf(backlog) else f"{backlog:.3f}"
            print(
                f"  rho={case.rho:4.2f}  {case.label:25s} "
                f"A_nu={amp:.3e}  H_*={h:.3e} GeV  M/M/1 backlog={backlog_text}  [{case.status}]"
            )
        else:
            print(
                f"  rho={case.rho:4.2f}  {case.label:25s} "
                f"lambda_shell>C_F  [{case.status}]"
            )
    check(amp_from_utilization(0.5) != amp_from_utilization(1.0), "subcritical printers give different amplitudes")
    check(mm1_mean_backlog(0.9) < math.inf, "ordinary stable queue mechanics allows rho<1")
    check(math.isinf(mm1_mean_backlog(1.0)), "exact capacity equality is critical, not a stationary stable M/M/1 point")

    print("\n[4] Ledger separation")
    print("  HBC entropy saturation is:")
    print("      printed entropy capacity per shell matches incoming Landauer entropy flux.")
    print("  Item-131 amplitude saturation is:")
    print("      topology-changing colour load per scalar shell equals C_F.")
    print("  The current text does not derive an identity between these two ledgers.")
    print("  That separation is why lambda_shell<C_F remains a valid countermodel:")
    print("      the entropy/area ledger can be saturated while the colour topological")
    print("      commit ledger has unused capacity.")

    print("\n[5] Exact status under the T1-T9 thermodynamic protocol")
    gates = [
        ("T1 event algebra", "PASS", "weight-4 topology current and colour-silent C_F load are specified"),
        ("T3 mean current", "PARTIAL", "upper bound on N_shell is derived; absolute mean requires maximal-throughput law"),
        ("T4 covariance/Fano", "CONDITIONAL", "F_eff=1 and S_j=1 are separate conditional ledgers already audited"),
        ("T5 correlation volume", "CONDITIONAL", "mode-local shell and no-horizon-covariance assumptions remain"),
        ("T6 observable map", "CONDITIONAL", "scalar-clock/delta-N wrapper is formal, still HBC-field conditional"),
        ("T8 alternatives", "PASS", "subcritical, critical, and overloaded queue branches are enumerated"),
        ("T9 verdict", "OPEN", "saturation equality needs one new dynamical theorem"),
    ]
    for gate, status, note in gates:
        print(f"  {gate:22s} {status:12s} {note}")

    print("\n[6] What would close it")
    print("  A true saturated-printer theorem must add, and derive from HBC/QEC:")
    print("    A. a monotone boundary-growth objective: use every coherent print slot;")
    print("    B. a backlog/latch law: unused capacity drives cell precipitation until rho=1;")
    print("    C. an exit law: rho>1 destroys coherent printing or triggers handoff;")
    print("    D. a finite-phase law: exact criticality is allowed as a finite latch,")
    print("       not as an eternal stationary M/M/1 queue.")

    print("\nVERDICT")
    print("  No full derivation from current HBC/QEC mechanics.  What is derived is")
    print("  the capacity inequality lambda_shell <= C_F and the conditional theorem")
    print("  that a maximal coherent printer must sit at lambda_shell=C_F.  The missing")
    print("  piece is precisely the maximal-throughput/backlog-latch dynamics linking")
    print("  HBC entropy saturation to the topological colour-load ledger.")
    print("exit 0 -- saturated-printer dynamics reduced to one maximal-throughput/backlog-latch theorem; not Locked.")


if __name__ == "__main__":
    main()
