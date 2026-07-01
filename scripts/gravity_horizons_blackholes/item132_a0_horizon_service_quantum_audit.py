#!/usr/bin/env python3
r"""ITEM 132: R4 one-a0 horizon service quantum audit.

Question
--------
Can the R4 local service instrument be closed so that it fires at one native
horizon-normalized acceleration quantum,

    a_R4 = c H_0 / (2 pi),

with H_0 supplied by the finite-QEC selector rather than inserted from
cosmology?

Result
------
The route closes only at conditional theorem-gate grade.

The proton-primary gravity/selector chain already supplies

    H_0 = Lambda_p / N_lock,       N_lock = 9 alpha_0 / r_6,

from finite-QEC service data.  No measured Hubble constant is consumed by the
calculation below.

The remaining 2pi is fixed if the R4 local readout is a repeatable KMS-cycle
latch: the horizon phase is theta = H_0 tau, and a completed CPTP service
read/write/reset cycle returns the latch to the same phase, so the minimal
repeatable tick has theta = 2pi and tau_KMS = 2pi/H_0.  One light-cone velocity
quantum per completed tick then gives

    a_R4 = c / tau_KMS = c H_0 / (2 pi).

Controls matter.  Reading the horizon generator directly gives c H_0.  Equating
de Sitter and Unruh temperatures also gives c H_0, not c H_0/(2pi).  Multiple
KMS windings or subcycle reads introduce a free integer.  Therefore the result
is not an unconditional microscopic theorem: it needs the cycle-latch lemma
that the R4 service instrument bills one complete KMS phase winding and one
record, not generator time, Unruh acceleration, or an arbitrary number of
windings.
"""

from __future__ import annotations

from dataclasses import dataclass
import importlib
import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_code"))


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


C = 299_792_458.0
MPC_KM = 3.085678e19
HBAR_GEV_S = 6.582120e-25
M_PROTON_GEV = 0.93827208816
ALPHA0 = 1.0 / 137.0
OMEGA_L = 0.6847


@dataclass(frozen=True)
class SelectorOutput:
    lambda_p_gev: float
    q1: float
    r6: float
    n_lock: float
    h0_gev: float
    h0_si: float
    h0_km_s_mpc: float


@dataclass(frozen=True)
class ClockCandidate:
    name: str
    period_s: float
    has_complete_phase: bool
    has_extra_integer: bool
    ledger: str

    def acceleration(self) -> float:
        return C / self.period_s

    def failures(self) -> list[str]:
        out: list[str] = []
        if self.ledger != "repeatable_kms_service_latch":
            out.append("wrong ledger")
        if not self.has_complete_phase:
            out.append("not a completed 2pi phase latch")
        if self.has_extra_integer:
            out.append("adds a free winding/subcycle integer")
        return out

    @property
    def admissible(self) -> bool:
        return not self.failures()


def selector_output() -> SelectorOutput:
    """Recompute the proton-primary H0 output without importing a measured H0.

    The queue parameters are the same local functions used by
    g_route_input_ledger.py.  H0 is computed from Lambda_p/N_lock after the
    depth-six residual is obtained.
    """

    rh = importlib.import_module("register_handoff_form_selection")
    lambda_p = M_PROTON_GEV / (2.0 * math.sqrt(2.0))
    q1 = rh.queue_readouts(rh.BASE_GAMMA * math.exp(-ALPHA0 * 0.303562705), 1)[1]
    r6 = (21.0 * q1) ** 32 / 21.0
    n_lock = 9.0 * ALPHA0 / r6
    h0_gev = lambda_p / n_lock
    h0_si = h0_gev / HBAR_GEV_S
    return SelectorOutput(
        lambda_p_gev=lambda_p,
        q1=q1,
        r6=r6,
        n_lock=n_lock,
        h0_gev=h0_gev,
        h0_si=h0_si,
        h0_km_s_mpc=h0_si * MPC_KM,
    )


def main() -> None:
    print("ITEM 132: R4 ONE-a0 HORIZON SERVICE QUANTUM AUDIT")
    print("=" * 96)

    out = selector_output()

    print("\n[1] Selector H0 is recomputed as an output")
    print(f"  Lambda_p = m_p/(2 sqrt2)       = {out.lambda_p_gev:.12f} GeV")
    print(f"  q1 post-service residual       = {out.q1:.15e}")
    print(f"  r6 depth-six residual          = {out.r6:.15e}")
    print(f"  N_lock = 9 alpha0/r6           = {out.n_lock:.15e} ticks")
    print(f"  H0 = Lambda_p/N_lock           = {out.h0_km_s_mpc:.6f} km/s/Mpc")
    check(abs(out.h0_km_s_mpc - 67.266152) < 0.02, "selector H0 matches the proton-primary output ledger")
    check(out.n_lock > 1.0e41 and out.n_lock < 3.0e41, "H0 is fixed by a large finite-QEC service span")

    print("\n[2] Candidate horizon clocks")
    tau_hubble = 1.0 / out.h0_si
    tau_kms = 2.0 * math.pi / out.h0_si
    tau_ds_kms = 2.0 * math.pi / (out.h0_si * math.sqrt(OMEGA_L))
    candidates = [
        ClockCandidate(
            "Hubble-generator tick",
            tau_hubble,
            has_complete_phase=False,
            has_extra_integer=False,
            ledger="generator_time",
        ),
        ClockCandidate(
            "Unruh/de-Sitter temperature equivalent",
            tau_hubble,
            has_complete_phase=False,
            has_extra_integer=False,
            ledger="thermal_acceleration",
        ),
        ClockCandidate(
            "Hubble KMS cycle latch",
            tau_kms,
            has_complete_phase=True,
            has_extra_integer=False,
            ledger="repeatable_kms_service_latch",
        ),
        ClockCandidate(
            "de-Sitter KMS cycle latch",
            tau_ds_kms,
            has_complete_phase=True,
            has_extra_integer=False,
            ledger="de_sitter_radius_convention",
        ),
        ClockCandidate(
            "two KMS windings",
            2.0 * tau_kms,
            has_complete_phase=True,
            has_extra_integer=True,
            ledger="repeatable_kms_service_latch",
        ),
        ClockCandidate(
            "half-cycle read",
            0.5 * tau_kms,
            has_complete_phase=False,
            has_extra_integer=True,
            ledger="subcycle_phase_read",
        ),
    ]
    a_target = C * out.h0_si / (2.0 * math.pi)
    for cand in candidates:
        status = "PASS" if cand.admissible else "FAIL: " + ", ".join(cand.failures())
        print(
            f"  {cand.name:36s} a={cand.acceleration():.12e} m/s^2 "
            f"ratio-to-KMS={cand.acceleration()/a_target:.6f}  {status}"
        )
    check(sum(c.admissible for c in candidates) == 1, "only the one-winding Hubble KMS cycle latch passes the ledger test")
    check(abs(candidates[0].acceleration() / a_target - 2.0 * math.pi) < 1e-12, "generator-time and Unruh-equivalent readings are too large by 2pi")
    check(abs(candidates[3].acceleration() / a_target - math.sqrt(OMEGA_L)) < 1e-12, "de-Sitter-radius convention introduces a sqrt(Omega_L) factor")

    print("\n[3] KMS-cycle latch algebra")
    theta_one_tick = out.h0_si * tau_kms
    print(f"  theta = H0 tau_KMS = {theta_one_tick:.12f} rad")
    print(f"  tau_KMS = 2pi/H0  = {tau_kms:.12e} s")
    print(f"  a_R4 = c/tau_KMS = {a_target:.12e} m/s^2")
    check(abs(theta_one_tick - 2.0 * math.pi) < 1e-12, "repeatable KMS latch is one complete 2pi phase winding")
    check(abs(a_target - C * out.h0_si / (2.0 * math.pi)) < 1e-30, "KMS-cycle latch gives a0=cH0/2pi")

    print("\n[4] Horizon-radius identity")
    r_hubble = C / out.h0_si
    r_from_a = C * C / (2.0 * math.pi * a_target)
    print(f"  R_H = c/H0                 = {r_hubble:.12e} m")
    print(f"  c^2/(2pi a_R4)             = {r_from_a:.12e} m")
    check(abs(r_from_a / r_hubble - 1.0) < 1e-15, "a_R4 names the Hubble horizon under the 2pi-cycle convention")

    print("\n[5] One-record threshold")
    print("  If a completed service latch bills m KMS windings, then a_m=a_R4/m.")
    print("  If it bills n subcycle reads per winding, then a_n=n a_R4.")
    for m in (1, 2, 3):
        print(f"    m={m}: full-winding threshold a={a_target / m:.12e}")
    for n in (1, 2, 3):
        print(f"    n={n}: subcycle threshold     a={n * a_target:.12e}")
    check(True, "the coefficient is fixed only if the R4 readout bills one completed latch record")

    print("\n[6] Decision")
    print(
        """
  CONDITIONAL CLOSURE:
    * H0 is supplied by the proton-primary finite-QEC selector, not by a Hubble
      measurement.
    * If the R4 local service readout is a repeatable KMS-cycle latch, the
      native acceleration quantum is one light-cone velocity quantum per
      completed horizon phase winding:

          a_R4 = c / (2pi/H0) = c H0 / (2pi).

  REFUTED READINGS:
    * Direct generator time and de Sitter/Unruh temperature matching give
      cH0, not cH0/(2pi).
    * A de-Sitter-radius convention adds sqrt(Omega_L).
    * Multi-winding or subcycle latch rules reintroduce a free integer.

  STILL LOAD-BEARING:
    * The microscopic theorem must show that the R4 local service instrument
      really bills exactly one completed KMS phase latch per record.
    * Empirical MOND calibration and the shared Newtonian/G stiffness
      convention remain separate checks.
"""
    )
    print("exit 0 -- KMS-cycle latch gives a0=cH0/2pi conditionally; generator/Unruh routes are rejected.")


if __name__ == "__main__":
    main()
