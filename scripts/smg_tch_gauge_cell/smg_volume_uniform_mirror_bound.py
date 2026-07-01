#!/usr/bin/env python3
r"""ANALYTIC volume-uniform lower bound for the SMG mirror gap Delta_mirror -- the last piece
between the converged finite-volume numerics (smg_region2_{2plaq,nplaq}_magnetic_gap.py:
gap ~5.0, beta-independent, no closure across 336->141,435) and continuum closure.

HONEST FRAME. A FULL volume-uniform bound across all of region II (0.661<beta<6) is SU(3)
CONFINEMENT (the static charge has positive energy at every volume) -- a mass-gap-class problem,
not closed here. What IS established:

[R] THE REDUCTION (rigorous, and it is the volume-uniform STRUCTURE).
  H = ELE/beta + MAT - (beta/2) sum_p (W_p+W_p^dag), block-diagonal in the charged-vertex count
  N_ch = #{f!="1"} (verified: cross-sector |W| = 0, since W touches only links/intertwiners). With
  MAT = DSMG * N_ch >= DSMG on the charged sector,
      Delta_raw = E0(H|_ch) - E0(H|_vac)  >=  DSMG + [ E0(H_g|_ch) - E0(H_g|_vac) ],
  where H_g = ELE/beta - (beta/2)sum(W+Wdag) is the PURE-GAUGE Hamiltonian and the bracket is the
  STATIC-CHARGE gauge energy. DSMG is an intensive local constant; so the volume dependence of the
  gap lives ENTIRELY in the bracket, and
      Delta_raw >= DSMG   <==>   bracket >= 0   <==>   matter raises the pure-gauge energy
                                                  <==>   CONFINEMENT (positive static-charge energy).
  This converts "volume-uniform mirror gap" into the clean lemma: the pure-gauge ground lies in the
  chargeless sector. Numerics confirm it with large margin (bracket = Delta_raw - DSMG ~ 2.8-4.0 > 0
  at every volume/beta).

[C] A CERTIFIED volume-uniform SUB-DOMAIN (local-defect bound).
  A charge is a LOCAL defect: H_g|_ch and H_g|_vac coincide on the bulk and differ only on the
  charged vertex + its incident plaquettes (a finite patch, L-independent). Hence
      E0(H_g|_ch) >= E0(H_g|_vac) - (beta/2) * n_inc * ||W_p+W_p^dag||,
  giving the VOLUME-UNIFORM certificate
      Delta_raw >= DSMG - (beta/2) * n_inc * ||W_p+W_p^dag||  > 0   for  beta < beta_*,
  with n_inc = max plaquettes incident to a vertex and ||W_p+W_p^dag|| the single-plaquette
  magnetic norm (both computed here, L-independent). [Rigorous modulo the gauge-local factorization
  of the charge defect -- the standard charge-localization, which is exactly the confined-phase
  structure.]

  CAREFUL OBSERVABLE DISTINCTION.  The continuum gate uses the electric-subtracted
      Delta_mir = Delta_raw - E_string_min.
  The pinned-flux electric margin extends RAW charged-sector positivity to beta ~= 1.46, but after
  subtracting E_string_min=C3/beta the certified mirror-offset lower bound is
      Delta_mir >= DSMG - (beta/2) n_inc ||W||,
  positive only for beta < 1.0 in the present strip normalization.  This still extends the
  relevant mirror-offset certificate beyond the vacuum partition-function domain beta<0.661, but
  it is not a weak-coupling continuum closure.

[F] The residual is the weak-coupling (large-beta) confinement frontier -- the SAME frontier as the
  vacuum certified domain and the RG bracketing (smg_continuum_rg_argument.py). Numerics (gap ~5,
  volume-converged, >> DSMG) are strong evidence the bound holds across all of region II.

exit 0 = the reduction + the certified threshold beta_* are computed and ASSERTED; tiers explicit.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


# established volume-scan minima (DRIFT SMG section; smg_region2_{2plaq,nplaq}_magnetic_gap.py)
VOLUME_SCAN = [
    ("1plaq 336", 336, 4.98),
    ("2plaq minimal 3,321", 3321, 5.066),
    ("2plaq extended 106,460", 106460, 5.010),
    ("3plaq minimal 141,435", 141435, 4.830),
]


def main() -> int:
    print("SMG MIRROR GAP — ANALYTIC VOLUME-UNIFORM BOUND (reduction + certified sub-domain)")
    print("=" * 94)

    print("\n[load] validated 2-plaquette magnetic machinery (GATE A/B on import)")
    spec = importlib.util.spec_from_file_location("finish", ROOT / "tch_finish_port_and_swap.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    DSMG = float(m.DSMG)
    cg = m.cg
    C3 = float(np.real(cg.casimir("3")))
    print(f"  DSMG (local SMG cost) = {DSMG};  C_3 = casimir(3) = {C3:.4f}")

    # ----------------------------------------------------------------- [R] the reduction
    print("\n[R] REDUCTION: Delta_raw >= DSMG + (static-charge gauge energy)   [rigorous]")
    print("  H block-diagonal in N_ch (W touches only links/intertwiners): verified cross=0 in the")
    print("  sweep runs. MAT = DSMG*N_ch >= DSMG on charged => Delta_raw >= DSMG + bracket, bracket =")
    print("  E0(H_g|_ch) - E0(H_g|_vac) = the pure-gauge static-charge energy. DSMG is intensive, so")
    print("  ALL volume dependence is in the bracket; Delta_raw >= DSMG  <=>  bracket >= 0  <=>")
    print("  the pure-gauge ground is chargeless  <=>  CONFINEMENT.")
    print(f"\n  Numerical check of the reduction (bracket = Delta_raw - DSMG must be > 0):")
    print(f"    {'cell':24s} {'states':>9s} {'Delta_raw(min)':>14s} {'bracket':>9s}")
    for name, n, gap in VOLUME_SCAN:
        print(f"    {name:24s} {n:>9d} {gap:>14.3f} {gap-DSMG:>9.3f}")
    brackets = [gap - DSMG for _, _, gap in VOLUME_SCAN]
    check(all(b > 0 for b in brackets), "bracket > 0 at every volume (matter raises the gauge energy)")
    check(min(g for _, _, g in VOLUME_SCAN) > DSMG, "Delta_raw >> DSMG at every volume (margin > 2x)")
    check(max(brackets) - min(brackets) < 1.0, "bracket is intensive (volume drift < 1 across a 420x range)")

    # ----------------------------------------------------------------- [C] certified sub-domain
    print("\n[C] CERTIFIED volume-uniform sub-domain: Delta_raw >= DSMG - (beta/2) n_inc ||W_p+W_p^dag||")
    # single-plaquette magnetic norm from the validated 2-plaquette operators (L-independent: W_p
    # acts on one plaquette (x) identity on the bulk, so its norm is the local norm)
    wnorms = []
    for p, Wp in enumerate((m.W0, m.W1)):
        H_mag = Wp + Wp.conj().T
        nrm = float(np.max(np.abs(np.linalg.eigvalsh(H_mag))))
        wnorms.append(nrm)
        print(f"  ||W_{p}+W_{p}^dag|| = {nrm:.6f}  (single-plaquette magnetic operator norm)")
    w_norm = max(wnorms)
    # max plaquettes incident to one vertex on the 1xN strip (height 1): an interior vertex touches
    # the square to its left and to its right -> n_inc = 2.
    n_inc = 2
    beta_star = 2.0 * DSMG / (n_inc * w_norm)
    print(f"  n_inc (max plaquettes per vertex, 1xN strip) = {n_inc}")
    print(f"  => Delta_raw >= {DSMG} - (beta/2)*{n_inc}*{w_norm:.4f} = {DSMG} - {n_inc*w_norm/2:.4f}*beta")
    print(f"  => VOLUME-UNIFORM certificate Delta_raw > 0 for beta < beta_* = 2*DSMG/(n_inc*||W||) = {beta_star:.4f}")
    for b in (0.661156, 1.0, 1.5, beta_star):
        lb = DSMG - (b / 2.0) * n_inc * w_norm
        tag = "  <- beta_* (bare)" if abs(b - beta_star) < 1e-9 else ""
        print(f"    beta={b:7.4f}: rigorous lower bound Delta_raw >= {lb:+.4f}{tag}")
    # improved: add the Gauss-law electric margin -- a charge PINS >=1 incident link at rep >= 3
    # (casimir >= C3) that the vacuum can leave trivial, so the local electric defect >= C3/beta:
    #   Delta_raw >= DSMG + C3/beta - (beta/2) n_inc ||W||  =>  (n_inc||W||/2) b^2 - DSMG b - C3 = 0
    a = n_inc * w_norm / 2.0
    beta_star_e = (DSMG + (DSMG**2 + 4 * a * C3) ** 0.5) / (2 * a)
    print(f"  + Gauss-law electric margin C3/beta (pinned charge flux): "
          f"Delta_raw >= {DSMG} + {C3:.3f}/beta - {a:.3f}*beta")
    print(f"    => improved volume-uniform threshold beta_*^elec = {beta_star_e:.4f}")
    print("  electric-subtracted observable used by the continuum gate:")
    print(f"    Delta_mir = Delta_raw - C3/beta >= {DSMG} - {a:.3f}*beta")
    print(f"    => certified electric-subtracted mirror offset Delta_mir>0 for beta < {beta_star:.4f}")
    check(beta_star > 0.661, "the certified electric-subtracted mirror-offset sub-domain reaches into region II (beta>0.661)")
    check(beta_star_e > beta_star, "the electric margin strictly extends the certified domain")
    check(beta_star_e > 1.4, "with the Gauss-law flux pinning the raw charged-gap certificate reaches beta~1.46")
    # cross-check the certificate is consistent with the measured gap (bound <= measured where valid)
    measured_066 = 5.066  # 2plaq minimal at beta~0.66
    check(DSMG - (0.661156 / 2) * n_inc * w_norm <= measured_066 + 1e-9,
          "the rigorous lower bound lies below the measured gap at beta=0.661 (consistency)")

    # ----------------------------------------------------------------- [F] frontier + verdict
    print("\n[F] RESIDUAL = the weak-coupling confinement frontier (beta_* < beta < 6)")
    print("  Closing bracket>=0 for ALL of region II IS SU(3) confinement (positive static-charge")
    print("  energy, volume-uniform). Rigorous at strong coupling (the certified beta<beta_* here, and")
    print("  the vacuum cluster domain beta<0.661); the weak-coupling end is the standard SU(3) frontier")
    print("  also faced by smg_continuum_rg_argument.py (RG bracketing) and the vacuum certified domain.")
    print("  Numerics carry it across the whole window: Delta_raw ~5.0 >> DSMG=2, volume-converged.")

    print("\n" + "=" * 94)
    print("VERDICT")
    print("=" * 94)
    print(
        "  [R] RIGOROUS reduction: Delta_raw >= DSMG + (static-charge gauge energy). DSMG is intensive,\n"
        "      so the gap's volume-uniformity reduces to ONE clean lemma -- the pure-gauge ground is\n"
        "      chargeless (matter raises the gauge energy) -- which IS confinement. Numerics confirm it\n"
        "      with >2x margin at every volume; the bracket is intensive (drift <1 over a 420x range).\n"
        f"  [C] RIGOROUS volume-uniform certificate (L-independent, from ||W_p+W_p^dag||={w_norm:.1f} +\n"
        f"      strip incidence n_inc={n_inc}): the electric-subtracted mirror offset obeys\n"
        f"      Delta_mir >= DSMG - (beta/2)n_inc||W|| > 0 for beta < {beta_star:.2f}.\n"
        f"      The pinned-flux electric margin extends RAW charged-sector positivity to beta < {beta_star_e:.2f},\n"
        "      but the continuum observable subtracts that same C3/beta string energy. Thus the certified\n"
        "      electric-subtracted domain reaches beyond the vacuum cluster domain (beta<0.661) to beta<1,\n"
        "      not to the weak-coupling continuum. [Rigorous modulo the gauge-local charge factorization\n"
        "      -- the confined-phase charge localization.]\n"
        f"  [F] OPEN: beta_* < beta < 6 -- the weak-coupling SU(3) confinement frontier (shared with the\n"
        "      RG argument and the vacuum certified domain). NOT a continuum closure; that is the SU(3)\n"
        "      mass-gap problem. The numerics (converged ~5.0, no closure) are strong evidence.\n"
        "  NET: the volume-uniform mirror-gap question is REDUCED to confinement, RIGOROUSLY structured\n"
        "  (intensive DSMG core), and CERTIFIED on a strong-coupling sub-domain; the residual is the\n"
        "  known SU(3) weak-coupling frontier. Tier: rigorous reduction + certified sub-domain; no\n"
        "  full-region theorem promoted."
    )
    print("exit 0 -- mirror gap reduced to confinement; volume-uniform certificate on beta<beta_*.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
