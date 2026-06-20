#!/usr/bin/env python3
r"""Planck-hierarchy service-span gate.

The local service-current stress tensor gives the correct *source* and the
lattice/bare coupling,

    M_Pl,bare^2 = 4 C Lambda^2,   C = 55/8.

It does not by itself contain the observed 10^38 enhancement in 1/G.  This
script tests the only framework-native nonlocal completion currently strong
enough to carry that enhancement: service-span/RG accrual over the lock time

    Z_G = 4 alpha_0^2 N_lock,
    N_lock = 9 alpha_0 / r6,

where r6 is the depth-six post-service residual from the active-demux handoff
chain.  Then

    M_Pl^2 = M_Pl,bare^2 Z_G
           = (4 C Lambda^2)(4 alpha_0^2 N_lock)
           = 110 alpha_0^2 Lambda^2 N_lock
           = 990 alpha_0^3 Lambda^2 / r6.

This is algebraically the proton-primary route's bare-alpha Planck mass, but
now expressed as a service-current RG theorem target:

  local tensor source  x  nonlocal service-span susceptibility.

Exit 0 means the hierarchy can be resolved *conditionally* by the service-span
theorem without consuming a measured Hubble/horizon value; companion audits now
own the operator statement and homogeneous R4 lift.  This script remains the
numerical factorization gate.
"""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))
rh = importlib.import_module("register_handoff_form_selection")


ALPHA0 = 1.0 / 137.0
ALPHA_D = 1.0 / 137.035999084
M_PROTON = 0.93827208816
LAMBDA_P = M_PROTON / (2.0 * math.sqrt(2.0))
C_HORIZON = 55.0 / 8.0
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25
M_PL_OBS = 1.220890e19


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def depth_six_residual() -> tuple[float, float]:
    """Return q1 and r6 from the live active-demux stationary branch."""

    gamma = rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705)
    _, q_post, _ = rh.queue_readouts(gamma, 1)
    r6 = (21.0 * q_post) ** 32 / 21.0
    return q_post, r6


def main() -> None:
    print("PLANCK HIERARCHY SERVICE-SPAN GATE")
    print("=" * 92)

    q1, r6 = depth_six_residual()
    n_lock = 9.0 * ALPHA0 / r6
    h0_out = (LAMBDA_P / n_lock) * MPC_KM / HBAR_GEV_S

    print("[1] Local service-current coupling")
    m_bare = 2.0 * math.sqrt(C_HORIZON) * LAMBDA_P
    print(f"    Lambda_p = m_p/(2 sqrt2) = {LAMBDA_P:.9f} GeV")
    print(f"    C = 55/8 = {C_HORIZON:.6f}")
    print(f"    M_Pl,bare = 2 sqrt(C) Lambda_p = {m_bare:.6f} GeV")
    check(1.6 < m_bare < 1.9, "local service tensor gives only the GeV-scale bare Planck mass")

    print("\n[2] Nonlocal service-span candidate")
    print(f"    active-demux q1(post-service) = {q1:.12e}")
    print(f"    r6=(21 q1)^32/21             = {r6:.12e}")
    print(f"    N_lock=9 alpha0/r6           = {n_lock:.12e} ticks")
    print(f"    H0_out=Lambda_p/N_lock       = {h0_out:.6f} km/s/Mpc")
    check(2.0e41 < n_lock < 3.0e41, "N_lock is the required QCD-to-Hubble tick count")
    check(abs(h0_out - 67.27) < 0.2, "H0 is output on the Planck side")

    z_g = 4.0 * ALPHA0 * ALPHA0 * n_lock
    m_macro_from_rg = m_bare * math.sqrt(z_g)
    m_macro_direct = LAMBDA_P * math.sqrt(110.0 * ALPHA0 * ALPHA0 * n_lock)
    m_macro_r6 = LAMBDA_P * math.sqrt(990.0 * ALPHA0**3 / r6)
    print("\n[3] Hierarchy factorization")
    print(f"    Z_G = 4 alpha0^2 N_lock                  = {z_g:.12e}")
    print(f"    sqrt(Z_G)                                = {math.sqrt(z_g):.12e}")
    print(f"    M_bare sqrt(Z_G)                         = {m_macro_from_rg:.12e} GeV")
    print(f"    Lambda sqrt(110 alpha0^2 N_lock)         = {m_macro_direct:.12e} GeV")
    print(f"    Lambda sqrt(990 alpha0^3/r6)             = {m_macro_r6:.12e} GeV")
    print(f"    observed M_Pl                            = {M_PL_OBS:.12e} GeV")
    print(f"    relative miss vs observed                = {m_macro_from_rg / M_PL_OBS - 1:+.4%}")
    check(abs(m_macro_from_rg / m_macro_direct - 1.0) < 1.0e-14, "local x service-span factor equals the 110-form")
    check(abs(m_macro_direct / m_macro_r6 - 1.0) < 1.0e-14, "N_lock form equals r6 form")
    check(abs(m_macro_from_rg / M_PL_OBS - 1.0) < 0.002, "bare-alpha macro Planck mass lands at sub-percent grade")

    print("\n[4] Horizon-input audit")
    # A pure measured-horizon route would set N = Lambda/H0_obs.  Here the
    # script computes N from r6 first, and H0 is printed only after.
    n_from_planck_h0 = LAMBDA_P / (67.36 / MPC_KM * HBAR_GEV_S)
    print(f"    N_lock from QEC residual = {n_lock:.12e}")
    print(f"    Lambda/H0(Planck)        = {n_from_planck_h0:.12e}")
    print(f"    ratio                    = {n_lock / n_from_planck_h0:.6f}")
    check(abs(n_lock / n_from_planck_h0 - 1.0) < 0.003, "QEC-derived lock span matches the empirical horizon tick count")
    print("    -> the hierarchy is not resolved by a local tensor; it is resolved only")
    print("       if the residual-fault service span is accepted as the origin of N_lock.")

    print("\n[5] Load-bearing clauses")
    clauses = [
        "C=55/8 value-level horizon record count",
        "post-service active-demux q1 readout",
        "generation-vertex loop correction C_gen≈0.30356",
        "depth-six residual recurrence r6=(21q1)^32/21",
        "9-touch service-span rule N_lock=9 alpha0/r6",
        "service-span RG/accrual theorem Z_G=4 alpha0^2 N_lock",
        "bare-alpha convention, or a separate dressed-alpha settlement",
    ]
    for i, clause in enumerate(clauses, 1):
        print(f"    T{i}: {clause}")

    print(
        """
[6] VERDICT
    CONDITIONAL RESOLUTION FOUND.

    The Planck hierarchy has a clean service-ledger factorization:

        M_Pl^2 = (4 C Lambda^2) * (4 alpha0^2 N_lock).

    The first factor is the local service-current / area-law source.  The
    second factor is nonlocal service-span susceptibility.  With
    N_lock=9 alpha0/r6 from the depth-six residual ledger, the observed
    hierarchy is reproduced and H0 is an output, not an input.

    This is stronger than the old horizon-Dirac route because the large number
    is no longer inserted as R_H/a0; it is computed as 9 alpha0/r6.  The
    operator-statement audit now derives Z_G=4 alpha0^2 N_lock and the
    9-touch rule inside the monitored service algebra, and the homogeneous R4
    lift audit closes the finite-to-FRW endpoint inside the current
    instrument.

    The live residual is therefore not a missing selector coefficient.  It is
    outside-sector completeness plus the alpha/Lambda convention: a hidden
    register, invalid-state dynamics, or non-R4 cosmological coupling would be
    new physics and would reopen the fallback to a sharp proton-primary/horizon
    near-coincidence.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
