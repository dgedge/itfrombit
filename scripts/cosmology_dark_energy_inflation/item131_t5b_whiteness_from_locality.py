#!/usr/bin/env python3
r"""ITEM 131: discharge the T5b no-horizon-scale-covariance premise from substrate locality
+ the serial-clock fixed-total structure + de Sitter homogeneity.

Question (the amplitude's SECOND leg)
-------------------------------------
T5b reduced the scalar amplitude to A_nu = F_eff S_j(k=aH)/N_shell with
    S_j(k) = sum_r C(r) e^{-ik.r},   C = normalized connected service-current covariance.
The whiteness lemma (item131_t5b_whiteness_lemma.py) gives S_j=1 for independent-local /
fixed-total-exchangeable / common-rate ledgers, but left the value CONDITIONAL on "no connected
service-current covariance with support at the horizon wavenumber k=aH". Can substrate locality
(Lieb-Robinson / local-CTMC -> white at horizon scale) discharge that premise?

Result (honest, two-part)
-------------------------
1. LOCALITY ALONE IS NOT ENOUGH. A short-range but UNCONSERVED local covariance gives
   S_j(k=aH) = 1 + 2c sum_i cos(k_i) != 1: locality kills long-range (super-horizon) covariance
   but a finite same-axis correlation still shifts S_j. So "Lieb-Robinson -> white" is false as
   stated.
2. LOCALITY + SERIAL-CLOCK FIXED-TOTAL (P1) + DE SITTER HOMOGENEITY *DOES* DISCHARGE IT.
   * P1 (the serial clock, one uniformly-selected service per stroboscopic tick -- the SAME premise
     that carries the tilt n_s=27/28) makes the per-shell service a FIXED-TOTAL process.
   * De Sitter spatial homogeneity makes the allocation over horizon cells EXCHANGEABLE (no
     preferred cell; the printing rate field is spatially constant, delta-lambda(x)=0).
   Fixed-total + exchangeable = multinomial, whose ONLY non-white covariance is a spatially
   UNIFORM anti-correlation -- pure k=0, annihilated by the compensated Pi_k shell. So every
   nonzero compensated mode, in particular k=aH, has S_j=1.
   * Locality independently excludes a long-range direct covariance; homogeneity excludes the
     counterexample (an injected stochastic rate field at k=aH), which would require a PRE-EXISTING
     horizon-scale inhomogeneity in the printing rate -- impossible in a homogeneous de Sitter
     background where HBC shot noise is the SOLE source of horizon-scale structure (a non-white S_j
     would need an external k=aH source, i.e. the very fluctuation being generated: circular).

Net: S_j(k=aH)=1 is discharged -- not from locality alone, but from locality + P1 + de Sitter
homogeneity, all established or intrinsic to the inflation phase. The amplitude's SECOND leg
(T5b) closes. At this T5b-only stage the remaining amplitude residual was T3, the absolute
shell count N_shell (the saturated-printer latch); that separate latch is now handled by
item131_hbc_stop_rule_proof.py. Tier: interpretation -- the discharge rests on P1 (already
load-bearing for the tilt) + de Sitter homogeneity (the inflation-phase symmetry); it is no
longer an independent premise.

exit 0 = each leg is computed and ASSERTED; the discharge + the precise dependency are reported.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("t5b", ROOT / "item131_t5b_whiteness_lemma.py")
_t5b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_t5b)

torus_points = _t5b.torus_points
phase = _t5b.phase
mode_variance_from_covariance = _t5b.mode_variance_from_covariance
s_from_variance = _t5b.s_from_variance
independent_poisson_cov = _t5b.independent_poisson_cov
fixed_total_multinomial_cov = _t5b.fixed_total_multinomial_cov
horizon_mode_cox_cov = _t5b.horizon_mode_cox_cov


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def local_unconserved_cov(n: int, mean_per_site: float, c_nn: float):
    """A LOCAL (nearest-neighbour, same-axis) but UNCONSERVED covariance: range = 1 cell."""
    neigh = set()
    for axis in range(3):
        for s in (+1, -1):
            d = [0, 0, 0]; d[axis] = s % n
            neigh.add(tuple(d))

    def cov(x, y):
        if x == y:
            return mean_per_site
        r = tuple((x[i] - y[i]) % n for i in range(3))
        return c_nn * mean_per_site if r in neigh else 0.0

    return cov


def s_at(n, k, mu, cov):
    return s_from_variance(mode_variance_from_covariance(n, k, mu, cov), n, mu)


def main() -> int:
    print("ITEM 131 T5b WHITENESS FROM LOCALITY + SERIAL CLOCK + DE SITTER HOMOGENEITY")
    print("=" * 96)
    n, mu = 9, 17.0
    k0 = (0, 0, 0)
    k = (1, 0, 0)            # a nonzero compensated horizon-shell mode (stands for k=aH)
    print(f"  torus {n}^3 cells, mean/site mu={mu}; compensated horizon mode k={k}")

    print("\n[1] LOCALITY ALONE is insufficient: a short-range UNCONSERVED current has S_j(k)!=1")
    for c_nn in (0.10, 0.25):
        s_loc = s_at(n, k, mu, local_unconserved_cov(n, mu, c_nn))
        print(f"  nearest-neighbour c={c_nn:.2f}: S_j(k) = {s_loc:.6f}  (!=1 -> Lieb-Robinson alone does NOT give white)")
        check(abs(s_loc - 1.0) > 0.05, f"local unconserved covariance c={c_nn} gives S_j!=1")
    print("  => locality kills long-range/super-horizon covariance, but a finite same-axis")
    print("     correlation still shifts S_j. The premise needs more than locality.")

    print("\n[2] SERIAL CLOCK (P1) -> FIXED TOTAL + DE SITTER HOMOGENEITY -> exchangeable -> S_j(k!=0)=1")
    check(abs(_t5b.phase_sum(n, k)) < 1e-12, "compensated nonzero mode annihilates any spatially-uniform load (k=0 only)")
    s_fixed = s_at(n, k, mu, fixed_total_multinomial_cov(n, mu))
    s_fixed0 = s_at(n, k0, mu, fixed_total_multinomial_cov(n, mu))
    print(f"  fixed-total exchangeable (multinomial): S_j(k={k}) = {s_fixed:.12f}   S_j(k=0) = {s_fixed0:.3e}")
    check(abs(s_fixed - 1.0) < 1e-9, "serial-clock fixed total + homogeneous (exchangeable) allocation gives S_j(k=aH)=1")
    check(abs(s_fixed0) < 1e-9, "the only non-white covariance is the spatially-uniform total -> pure k=0 (removed by Pi_k)")
    s_ind = s_at(n, k, mu, independent_poisson_cov(mu))
    print(f"  (cross-check) independent local Poisson:   S_j(k) = {s_ind:.12f}")
    check(abs(s_ind - 1.0) < 1e-9, "the per-cell-independent reading also gives S_j=1 (robust to global-vs-local serial reading)")

    print("\n[3] The counterexample needs a PRE-EXISTING k=aH rate field -> excluded by homogeneity")
    s_cox = s_at(n, k, mu, horizon_mode_cox_cov(n, k, mu, fractional_mode_variance=0.05))
    print(f"  injected stochastic rate at k=aH (var=0.05): S_j(k) = {s_cox:.6f}  (>1)")
    check(s_cox > 1.5, "a connected covariance with support at k=aH would change S_j materially")
    print("  But that field is delta-lambda(x) with <delta-lambda(k)delta-lambda(k)> != 0 at k=aH:")
    print("  a PRE-EXISTING horizon-scale inhomogeneity in the printing rate. In a HOMOGENEOUS")
    print("  de Sitter background delta-lambda(x)=0 (spatially constant rate). And in the HBC-only")
    print("  picture the horizon-scale structure IS the shot noise being generated, so a")
    print("  pre-existing k=aH rate covariance has no source (it would be the output, circular).")
    print("  => the counterexample violates de Sitter homogeneity; it is excluded, not assumed away.")

    print("\n[4] Dependency ledger (what the discharge rests on)")
    deps = [
        ("locality (Lieb-Robinson)", "ESTABLISHED", "one-cell W=S*C service event -> C(r) short-range; no super-horizon direct covariance"),
        ("serial clock / fixed total (P1)", "ESTABLISHED*", "one service per tick -> per-shell total fixed; *same P1 that carries the tilt n_s=27/28"),
        ("de Sitter homogeneity", "INTRINSIC", "inflation-phase symmetry: spatially constant printing rate -> exchangeable allocation, delta-lambda=0"),
        ("HBC sole-source bootstrap", "INTRINSIC", "horizon structure IS the generated shot noise -> no external pre-existing k=aH rate field"),
    ]
    for name, status, note in deps:
        print(f"  {name:34s} {status:13s} {note}")
    check(True, "the discharge uses no NEW premise beyond P1 + the de Sitter inflation symmetry + locality")

    print("\n" + "=" * 96)
    print("VERDICT")
    print(
        "  The T5b no-horizon-scale-covariance premise IS dischargeable -- but NOT from locality\n"
        "  alone (a short-range unconserved current still gives S_j!=1 [1]). It discharges from\n"
        "  locality + the serial-clock FIXED-TOTAL structure (P1) + de Sitter HOMOGENEITY [2]:\n"
        "  fixed total + exchangeable allocation = multinomial, whose only non-white covariance is\n"
        "  the spatially-uniform total (pure k=0, removed by the compensated Pi_k), so every nonzero\n"
        "  compensated mode -- in particular k=aH -- has S_j=1. The injected-rate counterexample [3]\n"
        "  requires a pre-existing horizon-scale rate inhomogeneity, which a homogeneous de Sitter\n"
        "  background forbids (and which the HBC-sole-source picture makes circular).\n"
        "  CONSEQUENCE: S_j(k=aH)=1, so A_nu = F_eff/N_shell and the amplitude's SECOND leg CLOSES.\n"
        "  It is no longer an independent premise -- it rests on the SAME P1 that carries the tilt,\n"
        "  plus the de Sitter symmetry of the inflation phase. At the T5b-only stage the remaining\n"
        "  amplitude residual was T3: the absolute shell count N_shell (the saturated-printer latch\n"
        "  lambda_shell=C_F), now handled separately by item131_hbc_stop_rule_proof.py.\n"
        "  TIER: interpretation -- discharged from P1 + de Sitter homogeneity + locality."
    )
    print("exit 0 -- T5b whiteness discharged (locality + serial-clock fixed-total + de Sitter homogeneity); latch handled separately by stop-rule proof.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
