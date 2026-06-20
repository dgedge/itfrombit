#!/usr/bin/env python3
r"""Black-hole absolute flux coefficient gate.

The previous freeze-surface script closes two structural pieces:

  * beta-one line scale: rho Lambda = phi/pi,
  * Hawking luminosity scaling: the Schwarzschild escape cone turns the local
    shell source into P proportional to M^-2.

This script asks whether the remaining absolute coefficient can be derived from
the QEC service ledger.

Result:

  * Matching the Stefan-Hawking coefficient

        P = M_P^4 / (15360 pi M^2)

    with the beta-one shell and outward-hemisphere escape cone requires

        Gamma0/Lambda = 0.002711306813
                      = 1.114347... * alpha0/3.

  * The nearest simple service-bandwidth form is

        Gamma0/Lambda = (10/27) alpha0,

    which gives 0.9971 of the Stefan-Hawking coefficient.  Equivalently, if
    this attempt rate is assumed and the shell is solved self-consistently, the
    near-beta-one root is beta=1.0014867.

  * That is not a derivation.  The actual finite Q graph exposes several
    nearby-looking but different counts (7/64, 14/13, weighted boundary current
    1.2219, etc.), none of which is 10/27.  Therefore the coefficient is
    reduced to a sharp attempt-rate theorem:

        derive Gamma0 = (10/27) alpha0 Lambda

    or replace it with the true spin/partial-wave greybody coefficient.
"""

from __future__ import annotations

from collections import Counter
import math


PHI = (math.sqrt(5.0) - 1.0) / 2.0
EPS_F = 1.0 / (2.0 * PHI)
ALPHA0 = 1.0 / 137.0
TARGET_GQ = {0: 1, 3: 11, 4: 22, 5: 38, 6: 54, 7: 41, 8: 25, 9: 14, 12: 2}
EDGES = [(i, j) for i in range(8) for j in range(i + 1, 8) if (i ^ j).bit_count() == 1]


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def bit(n: int, i: int) -> int:
    return (n >> i) & 1


def valid(n: int) -> bool:
    return (
        not (bit(n, 0) and bit(n, 1))
        and bit(n, 7) == bit(n, 6)
        and ((bit(n, 2) == 0) == ((bit(n, 3), bit(n, 4)) == (0, 0)))
    )


def strain(n: int) -> int:
    return sum(1 for i, j in EDGES if bit(n, i) != bit(n, j))


def invalid_states() -> list[int]:
    return [n for n in range(256) if not valid(n)]


def source_ladder(beta_eff: float) -> dict[int, float]:
    raw = {f: g * math.exp(-beta_eff * f) for f, g in TARGET_GQ.items()}
    z = sum(raw.values())
    return dict(sorted((f, v / z) for f, v in raw.items()))


def mean_f(beta_eff: float) -> float:
    src = source_ladder(beta_eff)
    return sum(f * p for f, p in src.items())


def beta_from_rho_lambda(rho_lambda: float) -> float:
    return 2.0 * math.pi * rho_lambda * EPS_F


def required_gamma_for_stefan(beta_eff: float, outward_hemisphere: bool = True) -> float:
    """Return Gamma0/Lambda needed for P=M_P^4/(15360 pi M^2).

    The near-horizon coefficient is

        C_model = (27 pi / 32) gamma eps_F <F> rho_lambda^4

    for the outward-hemisphere escape-cone convention.  If emission is isotropic
    over the full local sphere instead, only half the outward cone is populated,
    and the required gamma doubles.
    """

    rho_lambda = beta_eff / (2.0 * math.pi * EPS_F)
    target_c = 1.0 / (15360.0 * math.pi)
    denom = (27.0 * math.pi / 32.0) * EPS_F * mean_f(beta_eff) * rho_lambda**4
    gamma = target_c / denom
    return gamma if outward_hemisphere else 2.0 * gamma


def model_over_target(gamma: float, beta_eff: float, outward_hemisphere: bool = True) -> float:
    return gamma / required_gamma_for_stefan(beta_eff, outward_hemisphere)


def roots_for_gamma(gamma: float) -> list[float]:
    """Solve required_gamma_for_stefan(beta)=gamma over the useful shell range."""

    roots: list[float] = []
    lo_beta = 0.20
    hi_beta = 5.20
    n = 2000
    prev_x = lo_beta
    prev_y = required_gamma_for_stefan(prev_x) - gamma
    for i in range(1, n + 1):
        x = lo_beta + i * (hi_beta - lo_beta) / n
        y = required_gamma_for_stefan(x) - gamma
        if prev_y * y < 0.0:
            a, b = prev_x, x
            for _ in range(80):
                m = 0.5 * (a + b)
                if (required_gamma_for_stefan(a) - gamma) * (required_gamma_for_stefan(m) - gamma) <= 0:
                    b = m
                else:
                    a = m
            roots.append(0.5 * (a + b))
        prev_x, prev_y = x, y
    return roots


def finite_q_graph_controls(beta_eff: float = 1.0) -> dict[str, float]:
    q = invalid_states()
    qset = set(q)
    internal_edges = 0
    boundary_edges = 0
    weighted_boundary = 0.0
    z = sum(math.exp(-beta_eff * strain(s)) for s in q)
    for s in q:
        p = math.exp(-beta_eff * strain(s)) / z
        for k in range(8):
            t = s ^ (1 << k)
            if t in qset:
                internal_edges += 1
            else:
                boundary_edges += 1
                weighted_boundary += p * math.exp(-0.5 * beta_eff * (strain(t) - strain(s)))
    return {
        "directed_internal_per_Q": internal_edges / len(q),
        "directed_boundary_per_Q": boundary_edges / len(q),
        "boundary_fraction_all_cube_arrows": boundary_edges / (256.0 * 8.0),
        "weighted_boundary_current_raw": weighted_boundary,
    }


def main() -> None:
    print("BLACK-HOLE ABSOLUTE FLUX COEFFICIENT GATE")
    print("=" * 92)

    print("[1] Required attempt rate at the beta-one freeze shell")
    beta_eff = 1.0
    rho_lambda = beta_eff / (2.0 * math.pi * EPS_F)
    fbar = mean_f(beta_eff)
    gamma_req = required_gamma_for_stefan(beta_eff, outward_hemisphere=True)
    gamma_req_full = required_gamma_for_stefan(beta_eff, outward_hemisphere=False)
    print(f"    rho Lambda = {rho_lambda:.12f} = phi/pi")
    print(f"    <F>_beta=1 = {fbar:.12f}")
    print(f"    Gamma_req/Lambda, outward-hemisphere = {gamma_req:.12e}")
    print(f"    Gamma_req/Lambda, full-sphere        = {gamma_req_full:.12e}")
    print(f"    Gamma_req / (alpha0/3)               = {gamma_req / (ALPHA0 / 3.0):.12f}")
    check(abs(rho_lambda - PHI / math.pi) < 1.0e-15, "beta-one shell is phi/pi")

    print("\n[2] Candidate service bandwidths")
    candidates = {
        "alpha0/3": ALPHA0 / 3.0,
        "(10/27) alpha0": (10.0 / 27.0) * ALPHA0,
        "alpha0": ALPHA0,
        "Q-boundary 7/64 times alpha0": (7.0 / 64.0) * ALPHA0,
    }
    for name, gamma in candidates.items():
        print(
            f"    {name:<30s} gamma={gamma:.12e}  "
            f"P/P_SB={model_over_target(gamma, beta_eff):.9f}"
        )
    close = model_over_target((10.0 / 27.0) * ALPHA0, beta_eff)
    check(abs(close - 1.0) < 4.0e-3, "(10/27) alpha0 lands within 0.4 percent")
    check(abs(model_over_target(ALPHA0 / 3.0, beta_eff) - 1.0) > 0.05, "alpha0/3 alone is not the coefficient")

    print("\n[3] If Gamma0=(10/27) alpha0 Lambda, solve for the shell")
    gamma_1027 = (10.0 / 27.0) * ALPHA0
    roots = roots_for_gamma(gamma_1027)
    for beta in roots:
        rho = beta / (2.0 * math.pi * EPS_F)
        print(
            f"    beta={beta:.12f}, rho Lambda={rho:.12f}, "
            f"rho/(phi/pi)={rho/(PHI/math.pi):.12f}, <F>={mean_f(beta):.12f}"
        )
    check(any(abs(beta - 1.0) < 2.0e-3 for beta in roots), "10/27 alpha selects a near-beta-one root")

    print("\n[4] Actual finite-Q graph controls")
    controls = finite_q_graph_controls(beta_eff)
    for name, value in controls.items():
        print(f"    {name:<34s} {value:.12f}")
    print(f"    target attempt prefactor 10/27      {10.0/27.0:.12f}")
    check(abs(controls["boundary_fraction_all_cube_arrows"] - 7.0 / 64.0) < 1.0e-15, "raw boundary fraction is 7/64, not 10/27")
    check(all(abs(v - 10.0 / 27.0) > 1.0e-2 for v in controls.values()), "10/27 is not a raw finite-Q graph count")

    print(
        """
[5] VERDICT
    The absolute coefficient is reduced to a transfer-theorem gate.

    Progress:
      With the beta-one proper shell and escape-cone scaling, the standard
      Stefan-Hawking coefficient is equivalent to a local attempt rate

          Gamma0/Lambda = 0.002711306813.

      The service-bandwidth candidate (10/27) alpha0 Lambda lands at
      P/P_SB = 0.9971, and solving backward from that candidate selects a
      near-beta-one shell, beta=1.00149.

    What this script rules out:
      The actual finite Q graph does not supply 10/27 as a raw count.  The
      follow-up horizon service-stencil theorem derives
      Gamma0=(10/27) alpha0 Lambda only conditionally, by transferring the
      item-120 Landauer-Moore service alphabet to horizon severing.

    Remaining physics after that transfer:
      compute the true spin/partial-wave greybody coefficient and rerun this
      gate if the target is not the pure Stefan-Hawking scalar coefficient.

ALL ASSERTIONS PASSED"""
    )


if __name__ == "__main__":
    main()
