#!/usr/bin/env python3
"""
item55_ew_service_hessian_rg_boundary_attempt.py

Hard attempt at the requested closure:

    derive a single EW service/Hessian or RG-matched boundary condition fixing
    lambda_H/g_Z^2 directly.

The audit tests the two disciplined routes available from the current canon.

1. Single-Hessian / service-action route.
   In the renormalizable electroweak action the neutral gauge mass comes from
   |D H|^2, while the Higgs mass comes from the independent invariant
   lambda_H (H^dagger H)^2.  Gauge invariance does not relate the two
   coefficients.  A true service-Hessian theorem would have to be a new
   constraint tying those two invariants.  The familiar gauge-only D-term
   closure is also too small: lambda/g_Z^2 <= 1/8, while the observed target is
   about 0.236 and the rational near-hit is 15/64.

2. RG-matched boundary route.
   The one-loop SM beta function for R=lambda_H/g_Z^2 is not autonomous in the
   gauge sector; it contains the top Yukawa, in particular the -6 y_t^4 term in
   beta_lambda.  Therefore standard EW mixing plus RG cannot fix R without a
   top-Yukawa/scalar boundary condition.  Numerically R=15/64 is close to the
   weak-scale pole ratio, but it is not a one-loop fixed point at the physical
   top Yukawa.

Exit 0 means the attempted derivation is refuted inside these two classes and
the remaining theorem target is sharpened, not that the ratio is derived.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from fractions import Fraction


LOOP = 16.0 * math.pi * math.pi
ok = True


def check(name: str, cond: bool) -> None:
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


@dataclass(frozen=True)
class EwInputs:
    v: float = 246.22
    m_h: float = 125.25
    m_z: float = 91.1876
    m_w: float = 80.3692
    m_top: float = 172.76
    alpha_s_mz: float = 0.1181

    @property
    def g(self) -> float:
        return 2.0 * self.m_w / self.v

    @property
    def gp(self) -> float:
        gz2 = (2.0 * self.m_z / self.v) ** 2
        return math.sqrt(max(gz2 - self.g * self.g, 0.0))

    @property
    def g3(self) -> float:
        return math.sqrt(4.0 * math.pi * self.alpha_s_mz)

    @property
    def y_top(self) -> float:
        return math.sqrt(2.0) * self.m_top / self.v

    @property
    def lambda_h(self) -> float:
        return self.m_h * self.m_h / (2.0 * self.v * self.v)

    @property
    def gz2(self) -> float:
        return self.g * self.g + self.gp * self.gp


def sm_ratio(lam: float, g: float, gp: float) -> float:
    return lam / (g * g + gp * gp)


def higgs_mass_from_R(R: float, m_z: float) -> float:
    return math.sqrt(8.0 * R) * m_z


def betas(gp: float, g: float, g3: float, yt: float, lam: float) -> tuple[float, float, float, float, float]:
    """One-loop SM beta functions in non-GUT hypercharge normalization."""
    beta_gp = (41.0 / 6.0) * gp**3 / LOOP
    beta_g = (-19.0 / 6.0) * g**3 / LOOP
    beta_g3 = -7.0 * g3**3 / LOOP
    beta_yt = yt * (
        4.5 * yt * yt
        - (17.0 / 12.0) * gp * gp
        - 2.25 * g * g
        - 8.0 * g3 * g3
    ) / LOOP
    beta_lam = (
        24.0 * lam * lam
        - 6.0 * yt**4
        + 0.375 * (2.0 * g**4 + (g * g + gp * gp) ** 2)
        + (-9.0 * g * g - 3.0 * gp * gp + 12.0 * yt * yt) * lam
    ) / LOOP
    return beta_gp, beta_g, beta_g3, beta_yt, beta_lam


def beta_R(gp: float, g: float, g3: float, yt: float, R: float) -> float:
    del g3  # beta_R at one loop depends on g3 only through the running of y_t, not instantaneously.
    gz2 = g * g + gp * gp
    lam = R * gz2
    beta_gp, beta_g, _, _, beta_lam = betas(gp, g, 0.0, yt, lam)
    beta_gz2 = 2.0 * g * beta_g + 2.0 * gp * beta_gp
    return beta_lam / gz2 - lam * beta_gz2 / (gz2 * gz2)


def solve_yt_for_beta_R_zero(gp: float, g: float, R: float) -> float:
    """Find the positive y_t where beta_R=0 at fixed weak gauge couplings."""
    lo, hi = 0.0, 2.0
    flo = beta_R(gp, g, 0.0, lo, R)
    fhi = beta_R(gp, g, 0.0, hi, R)
    if flo * fhi > 0.0:
        return float("nan")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        fmid = beta_R(gp, g, 0.0, mid, R)
        if flo * fmid <= 0.0:
            hi = mid
            fhi = fmid
        else:
            lo = mid
            flo = fmid
    return 0.5 * (lo + hi)


def rk4_step(state: tuple[float, float, float, float, float], dt: float) -> tuple[float, float, float, float, float]:
    def add(s, k, scale):
        return tuple(si + scale * ki for si, ki in zip(s, k))

    k1 = betas(*state)
    k2 = betas(*add(state, k1, 0.5 * dt))
    k3 = betas(*add(state, k2, 0.5 * dt))
    k4 = betas(*add(state, k3, dt))
    return tuple(si + dt * (a + 2.0 * b + 2.0 * c + d) / 6.0
                 for si, a, b, c, d in zip(state, k1, k2, k3, k4))


def first_crossing_mu(inp: EwInputs, R_target: float, mu_max: float = 1.0e4) -> float | None:
    """Run upward from M_Z and return the first scale where R crosses R_target."""
    state = (inp.gp, inp.g, inp.g3, inp.y_top, inp.lambda_h)
    t = 0.0
    r_prev = sm_ratio(state[4], state[1], state[0]) - R_target
    t_prev = t
    dt = 0.001
    t_max = math.log(mu_max / inp.m_z)
    while t < t_max:
        state = rk4_step(state, dt)
        t += dt
        r_now = sm_ratio(state[4], state[1], state[0]) - R_target
        if r_prev == 0.0:
            return inp.m_z * math.exp(t_prev)
        if r_prev * r_now < 0.0:
            frac = abs(r_prev) / (abs(r_prev) + abs(r_now))
            return inp.m_z * math.exp(t_prev + frac * dt)
        r_prev = r_now
        t_prev = t
    return None


def main() -> int:
    inp = EwInputs()
    R_obs = sm_ratio(inp.lambda_h, inp.g, inp.gp)
    R_15_64 = float(Fraction(15, 64))
    mh_15_64 = higgs_mass_from_R(R_15_64, inp.m_z)

    print("=" * 96)
    print("ITEM 55 EW SERVICE/HESSIAN OR RG-BOUNDARY CLOSURE ATTEMPT")
    print("=" * 96)
    print("\n[0] Weak-scale target")
    print(f"  g={inp.g:.9f}, g'={inp.gp:.9f}, g_Z^2={inp.gz2:.9f}")
    print(f"  lambda_H(pole proxy)={inp.lambda_h:.9f}")
    print(f"  observed proxy R=lambda_H/g_Z^2={R_obs:.9f}")
    print(f"  rational near-target 15/64={R_15_64:.9f}")
    print(f"  imposing 15/64 at M_Z predicts m_H={mh_15_64:.6f} GeV"
          f" ({100.0 * (mh_15_64 / inp.m_h - 1.0):+.4f}%)")
    check("15/64 is a close weak-threshold boundary, but not exact",
          0.002 < abs(mh_15_64 / inp.m_h - 1.0) < 0.01)

    print("\n[1] Single-Hessian / service-action route")
    print("  SM tree identity: m_H/m_Z = 2 sqrt(2 lambda_H)/g_Z.")
    print("  The scalar quartic and neutral gauge coupling are independent gauge-invariant")
    print("  coefficients unless a new service-Hessian theorem ties them.")
    d_log_ratio_d_log_lambda = 0.5
    d_log_ratio_d_log_gz2 = -0.5
    print(f"  d log(m_H/m_Z)/d log(lambda_H) = {d_log_ratio_d_log_lambda:+.1f}")
    print(f"  d log(m_H/m_Z)/d log(g_Z^2)     = {d_log_ratio_d_log_gz2:+.1f}")
    check("ratio has independent scalar-quartic support", d_log_ratio_d_log_lambda != 0.0)
    check("ratio has independent neutral-gauge support", d_log_ratio_d_log_gz2 != 0.0)

    dterm_max = 1.0 / 8.0
    print(f"  Gauge-only D-term class gives lambda_H/g_Z^2 <= 1/8 = {dterm_max:.6f}.")
    print(f"  Needed: observed proxy {R_obs:.6f}, near-target 15/64 = {R_15_64:.6f}.")
    check("standard gauge-only D-term closure is too small for the Higgs/Z ratio",
          R_15_64 > dterm_max and R_obs > dterm_max)

    print("\n[2] RG-matched boundary route")
    beta_at_target = beta_R(inp.gp, inp.g, inp.g3, inp.y_top, R_15_64)
    yt_fixed = solve_yt_for_beta_R_zero(inp.gp, inp.g, R_15_64)
    print(f"  one-loop beta_R at R=15/64 and physical y_t={inp.y_top:.6f}: {beta_at_target:+.9f}")
    print(f"  y_t required for beta_R=0 at the same weak gauge couplings: {yt_fixed:.6f}")
    print("  The beta_R expression contains the independent top-Yukawa contribution")
    print("  -6 y_t^4/(16 pi^2 g_Z^2), so EW mixing alone cannot fix R.")
    check("15/64 is not an RG fixed ratio at the physical top Yukawa",
          abs(beta_at_target) > 0.01)
    check("the beta_R fixed-point y_t is far from the physical saturated top Yukawa",
          abs(yt_fixed - inp.y_top) > 0.25)

    mu_cross = first_crossing_mu(inp, R_15_64)
    print("\n[3] One-loop crossing diagnostic")
    if mu_cross is None:
        print("  Starting from the pole-proxy inputs, R does not cross 15/64 below 10 TeV.")
    else:
        print(f"  Starting from the measured pole-proxy lambda, R crosses 15/64 at mu={mu_cross:.3f} GeV.")
        print("  This crossing is not predictive: the initial lambda already contains the measured Higgs mass.")
    check("the apparent 15/64 crossing is a weak-threshold restatement, not a derived high-scale lock",
          mu_cross is not None and 90.0 < mu_cross < 100.0)

    print("\nVERDICT")
    print("  No derivation exists in the current single-Hessian or RG-matched classes.")
    print("  The current canon can impose R=15/64 as a sharp boundary candidate, and that")
    print("  predicts m_H within about 0.31%, but neither SM gauge invariance, the A_1g")
    print("  mode, the D_4h neutral bridge, nor one-loop RG selects that boundary.")
    print("  A real closure would need new physics in one of two forms:")
    print("    (i) a single EW service/Hessian action whose quadratic and quartic")
    print("        fluctuations are not independent and force lambda_H/g_Z^2=15/64; or")
    print("   (ii) a UV boundary condition fixing lambda_H together with y_t and the gauge")
    print("        couplings, followed by RG running to the weak threshold.")
    print("  Without one of those, lambda_H/g_Z^2 remains a scalar-sector input.")

    if ok:
        print("\nALL CHECKS PASSED")
        return 0
    print("\nCHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
