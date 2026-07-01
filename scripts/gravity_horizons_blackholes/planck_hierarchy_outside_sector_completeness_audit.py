#!/usr/bin/env python3
r"""Stress-test the remaining Planck-hierarchy outside-sector caveat.

Item 153's current route gives the observed Planck scale from

    M_Pl^2 = (4 C Lambda^2) * (4 alpha0^2 N_lock),
    N_lock = 9 alpha0/r6,

while the same span outputs

    H0 = Lambda/N_lock.

The live caveat is not the nine-touch rule or the service-span operator
inside the current instrument.  It is outside-sector completeness: could a
hidden register, invalid-state dynamics, or non-R4 cosmological coupling
quietly multiply the span or the gravitational susceptibility?

This audit makes that caveat quantitative.  Write

    eta_L      = chiral-anchor convention shift,
    eta_N      = hidden multiplicative change of the lock span,
    eta_Z      = hidden multiplicative change of gravitational susceptibility.

Then

    H0      -> eta_L eta_N^{-1} H0_base,
    M_Pl    -> eta_L sqrt(eta_N eta_Z) M_Pl_base.

Therefore a new outside sector is not a harmless knob: the joint M_Pl/H0
landing over-constrains it.  A hidden positive span requires a correlated
suppression of eta_Z roughly as eta_N^{-3} if H0 is kept fixed by a Lambda
shift.  That is not an in-instrument silent-stock freedom; it is new physics
with its own covariance/phantom/no-double-counting theorem.

Exit 0 means:
  * fixed-current-instrument alternatives (T=8, T=10, repeated sweep, 1%
    hidden stock) are precision-dead;
  * a continuous outside multiplier is allowed only as an explicitly
    correlated new sector, not as an unbilled hidden ledger inside item 153;
  * the hierarchy status is sharpened to "conditional on outside-sector
    no-extra-ledger/completeness", not "horizon-input open".
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
C_CELL = 55.0 / 8.0
M_PL_OBS = 1.220890e19
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25
H0_REF = rh.H0_KM_S_MPC


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def qec_residual() -> tuple[float, float]:
    gamma = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
    _, q_post, _ = rh.queue_readouts(gamma, 1)
    r6 = (21.0 * q_post) ** 32 / 21.0
    return q_post, r6


def observables(
    *,
    eta_lambda: float = 1.0,
    eta_span: float = 1.0,
    eta_z: float = 1.0,
) -> tuple[float, float]:
    """Return (M_Pl, H0) under multiplicative outside-sector deformations."""
    _, r6 = qec_residual()
    n_lock = 9.0 * ALPHA0 / r6
    m_bare_sq = 4.0 * C_CELL * LAMBDA_P**2
    z_g = 4.0 * ALPHA0**2 * n_lock
    m_pl = eta_lambda * math.sqrt(m_bare_sq * z_g * eta_span * eta_z)
    h0 = eta_lambda * LAMBDA_P / (eta_span * n_lock) * MPC_KM / HBAR_GEV_S
    return m_pl, h0


def miss_pair(m_pl: float, h0: float) -> tuple[float, float]:
    return m_pl / M_PL_OBS - 1.0, h0 / H0_REF - 1.0


def required_eta_z_if_h0_fixed_by_lambda(eta_span: float, m_base: float, h_base: float) -> tuple[float, float]:
    """For a chosen eta_span, fix H0 by eta_lambda and solve eta_z for M_Pl."""
    eta_lambda = (H0_REF / h_base) * eta_span
    eta_z = (M_PL_OBS / m_base) ** 2 / (eta_lambda**2 * eta_span)
    return eta_lambda, eta_z


def main() -> None:
    print("PLANCK HIERARCHY OUTSIDE-SECTOR COMPLETENESS AUDIT")
    print("=" * 96)

    print("[1] Owner-audit anchors")
    owner = (ROOT / "python_code/planck_hierarchy_operator_statement_audit.py").read_text()
    selector = (ROOT / "python_code/cosmological_selector_lock_theorem.py").read_text()
    lift = (ROOT / "python_code/item131_r4_homogeneous_lift_theorem.py").read_text()
    episode = (ROOT / "python_code/a1_episode_count_closure_attempt.py").read_text()
    anchors = [
        ("service-span susceptibility", "Z_G = s_cell(N)/C = 4 alpha0^2 N", owner),
        ("nine-touch lock span", "N_lock=T alpha0/r6", owner),
        ("selector lock endpoint", "N_physical = N_lock = 9 alpha0/r6", selector),
        ("homogeneous R4 lift", "chi_R4(a)=a/a_lock", lift),
        ("one episode per physical cell", "completion => exactly one open episode per physical cell", episode),
        ("outside-sector caveat named", "outside-sector completeness", owner + lift + episode),
    ]
    for label, needle, text in anchors:
        ok = needle in text
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}")
        if not ok:
            raise AssertionError(label)

    print("\n[2] Baseline joint landing")
    q_post, r6 = qec_residual()
    n_lock = 9.0 * ALPHA0 / r6
    m_base, h_base = observables()
    m_miss, h_miss = miss_pair(m_base, h_base)
    print(f"    q1(post-service)       = {q_post:.12e}")
    print(f"    r6                     = {r6:.12e}")
    print(f"    N_lock                 = {n_lock:.12e}")
    print(f"    M_Pl(base)             = {m_base:.12e} GeV ({m_miss:+.5%})")
    print(f"    H0(base)               = {h_base:.6f} km/s/Mpc ({h_miss:+.5%})")
    check(abs(m_miss) < 2.0e-3, "M_Pl lands at sub-percent grade")
    check(abs(h_miss) < 2.0e-3, "H0 lands at sub-percent grade from the same span")

    print("\n[3] Fixed-instrument alternatives are precision-dead")
    variants = [
        ("T=8 no post-service readout", 8.0 / 9.0, 1.0, 1.0),
        ("T=10 bills blind slot", 10.0 / 9.0, 1.0, 1.0),
        ("two repeated sweeps", 2.0, 1.0, 1.0),
        ("1% hidden positive stock", 1.01, 1.0, 1.0),
        ("1% hidden common stock+gravity", 1.01, 1.0, 1.01),
    ]
    print(f"    {'variant':<34s} {'eta_N':>9s} {'eta_Z':>9s} {'M miss':>11s} {'H0 miss':>11s}")
    for name, eta_n, eta_l, eta_z in variants:
        m, h = observables(eta_span=eta_n, eta_lambda=eta_l, eta_z=eta_z)
        dm, dh = miss_pair(m, h)
        print(f"    {name:<34s} {eta_n:9.5f} {eta_z:9.5f} {dm:+11.3%} {dh:+11.3%}")
    check(abs(miss_pair(*observables(eta_span=8.0 / 9.0))[1]) > 0.10, "T=8 fails the H0 span")
    check(abs(miss_pair(*observables(eta_span=10.0 / 9.0))[1]) > 0.09, "T=10 fails the H0 span")
    check(abs(miss_pair(*observables(eta_span=2.0))[1]) > 0.40, "a repeated sweep is catastrophic")
    check(abs(miss_pair(*observables(eta_span=1.01))[1]) > 0.01, "one-percent silent stock is already visible in H0")

    print("\n[4] Joint M_Pl/H0 constraints on an outside sector")
    eta_z_only = (M_PL_OBS / m_base) ** 2
    eta_n_h0_only = h_base / H0_REF
    print(f"    eta_Z required if only gravity susceptibility moves      = {eta_z_only:.12f}")
    print(f"    eta_N required if only the span moves and Lambda fixed   = {eta_n_h0_only:.12f}")
    check(eta_z_only < 1.0, "a positive hidden gravitational record worsens the current M_Pl landing")
    check(abs(eta_n_h0_only - 1.0) < 2.0e-3, "span-only correction is already pinned at the per-mille level")

    # If eta_Z is not allowed, the two observations solve eta_N and eta_Lambda.
    h_ratio = H0_REF / h_base
    m_ratio = M_PL_OBS / m_base
    eta_n_no_z = (m_ratio / h_ratio) ** (2.0 / 3.0)
    eta_l_no_z = h_ratio * eta_n_no_z
    print("\n    with eta_Z fixed to 1, the joint solve is:")
    print(f"      eta_N      = {eta_n_no_z:.12f}")
    print(f"      eta_Lambda = {eta_l_no_z:.12f}")
    check(abs(eta_n_no_z - 1.0) < 2.0e-3, "no-Z joint solve leaves no O(1) hidden span")
    check(abs(eta_l_no_z - 1.0) < 1.0e-3, "the remaining shift is an alpha/Lambda convention-sized correction")

    print("\n    if H0 is held fixed by moving Lambda, eta_Z must compensate cubically:")
    print(f"    {'eta_N':>8s} {'eta_Lambda':>13s} {'eta_Z required':>16s} {'reading':>28s}")
    for eta_n in [1.0, 1.001, 1.01, 1.10, 2.0]:
        eta_l, eta_z = required_eta_z_if_h0_fixed_by_lambda(eta_n, m_base, h_base)
        if eta_z < 0.0:
            reading = "phantom"
        elif eta_z < 0.9:
            reading = "large suppression"
        elif eta_z < 0.99:
            reading = "percent suppression"
        else:
            reading = "near current"
        print(f"    {eta_n:8.3f} {eta_l:13.9f} {eta_z:16.9f} {reading:>28s}")
    _, eta_z_two_sweeps = required_eta_z_if_h0_fixed_by_lambda(2.0, m_base, h_base)
    check(eta_z_two_sweeps < 0.2, "a two-sweep outside span would require an eightfold gravity suppression")

    print("\n[5] Canonical interpretation")
    print("    The current finite/instrumental theorems have removed the ordinary")
    print("    hidden-stock freedom: T=8, T=10, repeated sweeps, and percent-level")
    print("    silent reserves spoil M_Pl, H0, or both.  The only remaining caveat")
    print("    is a genuine outside sector with its own correlated contribution to")
    print("    eta_N, eta_Z, and/or eta_Lambda.")
    print("    Because H0 reads eta_Lambda/eta_N while M_Pl reads")
    print("    eta_Lambda*sqrt(eta_N eta_Z), such a sector must satisfy a joint")
    print("    covariance condition.  It cannot be introduced as an unbilled")
    print("    positive ledger without becoming visible.")

    print(
        """
================================================================================
VERDICT
  Inside the current R1-R4 service instrument, the observed Planck hierarchy is
  no longer horizon-input: the same monitored lock span gives H0 and supplies
  the susceptibility factor Z_G=4 alpha0^2 N_lock.

  The remaining residual is now sharply named:

      Outside-sector no-extra-ledger/completeness theorem.

  A hidden register, invalid-state dynamics, or non-R4 cosmological coupling is
  not a free coefficient in item 153.  It is new physics that must pass the
  joint M_Pl/H0 covariance equations above; ordinary positive hidden stock is
  already over-constrained.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
