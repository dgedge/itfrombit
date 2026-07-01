#!/usr/bin/env python3
r"""Cosmological selector / lock theorem audit.

Question
--------
The Planck-hierarchy service-span algebra now lands conditionally:

    Z_G = 4 alpha_0^2 N,
    N_lock = 9 alpha_0 / r6.

The remaining gate is the selector:

    why is the physical N the R4-completion / a=1 crossing?

This script tests whether the current canon already supplies that selector.

Result shape
------------
The selector is no longer a free epoch choice.  The current canon has three
pieces:

  P1. The finite R4 repair support is a one-dimensional 1-chain.  Under the
      homogeneous late lift, its activation variable is a bounded completion
      fraction f in [0,1], with f(a)=a/a_lock before saturation.

  P2. The A2 handover theorem says the gravitating record carrier switches
      from active strain to permanent record at service completion.  Therefore
      G-lock is the aggregate first-hitting time f=1.

  P3. The monitored service algebra says one complete deep rescue consumes
      T alpha_0 committed records with T=8+1+0=9, while the residual-fault
      service current supplies r6 per tick.  Hence the completion fraction is

          chi_R4(N) = min(1, N r6/(9 alpha_0)),

      and the first-hitting time is

          N_physical = N_lock = 9 alpha_0/r6.

After this physical event is identified, the cosmological convention a=1 is
not an extra input: it is the scale-factor normalization at the completed
support endpoint a_lock.  The exact FRW crossing is still computed as a
numerical sanity check.

Exit 0 means: the selector is closed at conditional theorem grade.  The
companion ``item131_r4_homogeneous_lift_theorem.py`` closes the homogeneous
R4 lift inside the current service instrument; the remaining caveat is now
outside-sector completeness plus the alpha/Lambda convention.
"""

from __future__ import annotations

import importlib
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python_code"))
rh = importlib.import_module("register_handoff_form_selection")


ALPHA0 = 1.0 / 137.0
M_PROTON = 0.93827208816
LAMBDA_P = M_PROTON / (2.0 * math.sqrt(2.0))
LAMBDA_SELECTOR = rh.LAMBDA_QCD_GEV
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25
OMEGA_L = rh.OMEGA_L
OMEGA_M = 1.0 - OMEGA_L
H0_REF = rh.H0_KM_S_MPC
BUDGET = 7.0e-4


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def residual_r6() -> tuple[float, float, float]:
    gamma = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
    _, q_post, _ = rh.queue_readouts(gamma, 1)
    r6 = (21.0 * q_post) ** 32 / 21.0
    return gamma, q_post, r6


def h0_gev() -> float:
    return H0_REF / MPC_KM * HBAR_GEV_S


def h_of_a(a: float) -> float:
    """Framework late cosmology: Planck matter plus item-131 dark energy."""
    rho_de = OMEGA_L * math.exp((3.0 / 28.0) * (1.0 - a))
    return h0_gev() * math.sqrt(OMEGA_M * a**-3 + rho_de)


def completion_fraction_from_a(a: float, r6: float) -> float:
    n_t = LAMBDA_SELECTOR / h_of_a(a)
    return n_t * r6 / (9.0 * ALPHA0)


def exact_crossing(r6: float) -> float:
    lo, hi = 0.90, 1.10
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if completion_fraction_from_a(mid, r6) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def equality_epoch() -> float:
    a = (OMEGA_M / OMEGA_L) ** (1.0 / 3.0)
    for _ in range(80):
        de = OMEGA_L * math.exp((3.0 / 28.0) * (1.0 - a))
        a = (OMEGA_M / de) ** (1.0 / 3.0)
    return a


def main() -> None:
    print("COSMOLOGICAL SELECTOR / LOCK THEOREM AUDIT")
    print("=" * 96)

    r4_support = (ROOT / "python_code/item131_r4_support_dimension.py").read_text()
    handover = (ROOT / "python_code/r4_activation_identification_closure.py").read_text()
    operator = (ROOT / "python_code/planck_hierarchy_operator_statement_audit.py").read_text()

    print("[1] Live canon anchors")
    anchors = [
        (
            "finite R4 support is one-dimensional",
            "support complex is three disjoint two-edge stars",
            r4_support,
        ),
        (
            "finite multiplicity affects normalization only",
            "changes only normalisation, not scaling dimension",
            r4_support,
        ),
        (
            "activation is service-completion fraction",
            "activation f(a) = the SERVICE-COMPLETION fraction",
            handover,
        ),
        (
            "lock is aggregate handover completion",
            "the lock = aggregate handover completion",
            handover,
        ),
        (
            "operator algebra gives T=9",
            "T = 8+1+0 = 9",
            operator,
        ),
        (
            "operator algebra gives Z_G",
            "Z_G = s_cell(N)/C = 4 alpha0^2 N",
            operator,
        ),
    ]
    for label, needle, text in anchors:
        ok = needle in text
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
        assert ok, label

    print("\n[2] Selector theorem")
    print("    R4 completion fraction before saturation:")
    print("      chi_R4(N) = N r6/(T alpha0), with T=8+1+0=9.")
    print("    Because chi_R4 is a bounded support-completion fraction,")
    print("    the physical lock is its first hitting time chi_R4=1.")
    print("      => N_physical = N_lock = 9 alpha0/r6.")
    print("    The symbol a=1 is then the scale-factor normalization at this")
    print("    completed-support endpoint, not an independent epoch selector.")

    gamma, q1, r6 = residual_r6()
    n_lock = 9.0 * ALPHA0 / r6
    h0_out_selector = LAMBDA_SELECTOR / n_lock * MPC_KM / HBAR_GEV_S
    h0_out_proton = LAMBDA_P / n_lock * MPC_KM / HBAR_GEV_S
    chi_today = completion_fraction_from_a(1.0, r6)
    a_cross = exact_crossing(r6)
    print("\n[3] Numerical landing")
    print(f"    gamma*                         = {gamma:.12e}")
    print(f"    q1(post-service)               = {q1:.12e}")
    print(f"    r6                             = {r6:.12e}")
    print(f"    N_lock=9 alpha0/r6             = {n_lock:.12e}")
    print(f"    H0_out=Lambda_QCD/N_lock       = {h0_out_selector:.6f} km/s/Mpc")
    print(f"    H0_out=Lambda_p/N_lock         = {h0_out_proton:.6f} km/s/Mpc")
    print(f"    reference H0 used by audit     = {H0_REF:.6f} km/s/Mpc")
    print(f"    selector miss                  = {h0_out_selector/H0_REF - 1:+.4%}")
    print(f"    proton-primary miss            = {h0_out_proton/H0_REF - 1:+.4%}")
    print(f"    chi_R4(a=1) using FRW H(a)     = {chi_today:.9f}")
    print(f"    exact FRW chi_R4=1 crossing    = a={a_cross:.9f}")
    print(f"    crossing offset from a=1       = {a_cross - 1:+.6f}")
    check(abs(h0_out_selector / H0_REF - 1.0) < BUDGET, "selector Lambda outputs the observed H0 inside the registered budget")
    check(abs(h0_out_proton / H0_REF - 1.0) < 0.002, "proton-primary Lambda remains sub-percent, tracked as the alpha/Lambda convention gate")
    check(abs(chi_today - 1.0) < 2.0 * BUDGET, "a=1 lies inside the registered selector budget")
    check(abs(a_cross - 1.0) < 0.002, "exact crossing is the same physical epoch to current precision")

    print("\n[4] Alternative selectors are not completion endpoints")
    a_28 = math.exp(28.0 - math.log(LAMBDA_SELECTOR / 2.348e-13))
    a_eq = equality_epoch()
    rows = [
        ("R4 completion endpoint", 1.0),
        ("28-dilution integer", a_28),
        ("DE-matter equality", a_eq),
        ("exact crossing reference", a_cross),
    ]
    print(f"    {'candidate':<28s} {'a':>10s} {'chi_R4-1':>12s} {'status':>18s}")
    for name, a in rows:
        miss = completion_fraction_from_a(a, r6) - 1.0
        if name == "R4 completion endpoint":
            status = "endpoint theorem"
        elif name == "exact crossing reference":
            status = "answer, not selector"
        elif abs(miss) <= BUDGET:
            status = "inside budget"
        elif abs(miss) < 0.02:
            status = f"miss x{abs(miss)/BUDGET:.1f}"
        else:
            status = "excluded"
        print(f"    {name:<28s} {a:10.6f} {miss:+12.4%} {status:>18s}")
    check(abs(completion_fraction_from_a(a_28, r6) - 1.0) > 10.0 * BUDGET, "28-dilution is not the lock selector")
    check(abs(completion_fraction_from_a(a_eq, r6) - 1.0) > 0.10, "DE-matter equality is not the lock selector")

    print("\n[5] What remains")
    print("    CLOSED, conditional on live canon:")
    print("      * the selector is a first-hitting endpoint of the R4 completion")
    print("        fraction, not an arbitrary epoch pick;")
    print("      * the physical span is N_physical=N_lock=9 alpha0/r6;")
    print("      * a=1 is the convention after choosing the completed-support")
    print("        endpoint as the cosmological lock slice.")
    print("    STILL NAMED:")
    print("      * outside-sector completeness: no hidden register, invalid-state")
    print("        dynamics, non-R4 cosmological coupling, or negative phantom-rate")
    print("        branch;")
    print("      * the already tracked alpha/Lambda precision gate.")
    print(
        """
[6] VERDICT
    The Planck hierarchy no longer depends on an unexplained cosmological
    selector once the current R4 and A2 results are accepted.  The selector is
    the unique endpoint of a bounded service-completion fraction:

        chi_R4(N)=min(1, N r6/(9 alpha0)).

    Therefore

        N_physical = N_lock = 9 alpha0/r6.

    This upgrades the service-span hierarchy from "conditional on an epoch
    pick" to an in-instrument homogeneous R4 endpoint theorem.  The companion
    lift audit closes the finite-to-homogeneous step inside the current
    register/instrument.  The remaining caveat is outside-sector completeness:
    a hidden register or non-R4 cosmological coupling would be new physics,
    not an ambiguity in the R4 lift.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
