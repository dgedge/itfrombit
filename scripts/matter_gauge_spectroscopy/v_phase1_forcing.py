#!/usr/bin/env python3
r"""v_phase1_forcing.py -- Phase 1 of the EW-scale program: force the radiative reframing of v.

CLAIM (the forcing theorem, by elimination): in the dimensionless substrate the Higgs mass^2 mu^2 can
only be built from the available scales {Lambda, M_P} times dimensionless couplings. With a non-fine-
tuned quartic (lambda ~ O(0.1-1), as the Higgs has), the EW scale v = sqrt(-mu^2/lambda) then forces:

  (A) mu^2 ~ M_P^2 (O(1) coeff)  ->  v ~ M_P  -> v/M_P ~ 1, but observed 2e-17. MISS by ~17 orders.
  (B) mu^2 ~ Lambda^2 (O(1) coeff) -> v ~ Lambda -> v/Lambda ~ 1, but observed 741 -> needs lambda ~ 2e-6
      (absurd fine-tuning). EXCLUDED.
  (C) mu^2 RADIATIVELY/coupling-suppressed: mu^2 ~ alpha_0^n M_P^2  ->  v/M_P = alpha_0^(n/2)/sqrt(lambda).
      Viable: at lambda ~ lambda_Higgs, n ~ 16, i.e. v/M_P ~ alpha_0^8 -- the only option with no fine-tuning.

So the EW scale is FORCED to be a radiative (coupling-power-suppressed) transmutation scale, not an O(1)
tree input -- exactly Coleman-Weinberg / dimensional transmutation, like Lambda_QCD. This reframes
No-Go-1: v is not an arbitrary 2nd anchor; it is the radiative EW scale, and the open question is the
POWER n and the quartic lambda (Phases 2-3), not "why is there a 2nd scale".

HOPEFUL CONSISTENCY CHECK: the failed-alphabet hint was v/M_P ~ 2.50 * alpha_0^8 (exponent 7.81). Read
radiatively, v/M_P = alpha_0^8/sqrt(lambda) => the prefactor 2.50 = 1/sqrt(lambda) => lambda ~ 0.16,
the SAME ballpark as the measured Higgs quartic lambda_Higgs = m_H^2/2v^2 ~ 0.13. So the "2.50*alpha_0^8"
that looked like numerology is, radiatively, (mu^2 ~ alpha_0^16 M_P^2) with a Higgs-like quartic -- a
genuine sign the reframing is the right physics. The power n/2 = 8 is suggestively the [8,4,4] byte.

SCOPE: this is the forcing + the consistency reframing (Phase 1). It is SCHEMATIC (O(1) Coleman-Weinberg
factors uncomputed) and does NOT derive n or lambda -- those are Phases 2-3. The lambda ~ 0.16 vs 0.13
is a ballpark (~24%) match, not a precision result.
"""
import math

# constants (GeV)
v = 246.22                 # EW vev
LAMBDA = 0.332             # Lambda_QCD anchor (m_p / 2sqrt2)
M_P = 1.2209e19            # Planck mass (full)
ALPHA0 = 1.0 / 137.036
lam_Higgs = 125.25 ** 2 / (2 * v ** 2)   # measured Higgs quartic m_H^2/2v^2


def main():
    vMP = v / M_P
    vL = v / LAMBDA
    a8 = ALPHA0 ** 8
    print("=== Phase 1: forcing the radiative reframing of the EW scale v ===\n")
    print(f"  v = {v} GeV,  Lambda = {LAMBDA} GeV,  M_P = {M_P:.3e} GeV,  alpha_0 = 1/137.036")
    print(f"  v/M_P = {vMP:.3e}   (= alpha_0^{math.log(vMP)/math.log(ALPHA0):.2f});   v/Lambda = {vL:.1f}")
    print(f"  alpha_0^8 = {a8:.3e};   v/M_P / alpha_0^8 = {vMP/a8:.3f}  (the ew_nuclear 'prefactor 2.50')\n")

    print("  Elimination over mu^2 = (dimensionless coeff) x (scale)^2, with a non-fine-tuned lambda~O(0.1):")
    # (A) mu^2 ~ M_P^2
    v_from_MP = M_P * math.sqrt(1 / lam_Higgs)       # v if mu^2 = M_P^2
    print(f"    (A) mu^2 ~ M_P^2   -> v ~ {v_from_MP:.2e} GeV  (v/M_P~{v_from_MP/M_P:.2f}); "
          f"MISS by {abs(math.log10(v_from_MP/v)):.0f} orders -> EXCLUDED (no hierarchy).")
    # (B) mu^2 ~ Lambda^2 -> needs tiny lambda
    lam_needed_B = (LAMBDA / v) ** 2
    print(f"    (B) mu^2 ~ Lambda^2 -> needs lambda = (Lambda/v)^2 = {lam_needed_B:.2e}  "
          f"(vs Higgs {lam_Higgs:.3f}) -> EXCLUDED (absurd fine-tuning).")
    # (C) radiative mu^2 ~ alpha_0^n M_P^2
    n_at_higgs = 2 * math.log(vMP * math.sqrt(lam_Higgs)) / math.log(ALPHA0)
    lam_implied = (a8 / vMP) ** 2                    # from v/M_P = alpha_0^8/sqrt(lam)
    print(f"    (C) mu^2 ~ alpha_0^n M_P^2 -> v/M_P = alpha_0^(n/2)/sqrt(lambda):")
    print(f"          at lambda=lambda_Higgs={lam_Higgs:.3f}: n = {n_at_higgs:.1f}  (~16 = 2x8, the byte)")
    print(f"          OR back out lambda from the alpha_0^8 fit: lambda = {lam_implied:.3f}  "
          f"vs lambda_Higgs={lam_Higgs:.3f}  (ratio {lam_implied/lam_Higgs:.2f})")

    print("\n[verdict] PHASE 1 -- the radiative reframing is FORCED, and consistently so:")
    print("  - the O(1) tree options for mu^2 fail: ~M_P^2 misses by 17 orders, ~Lambda^2 needs a")
    print("    ~2e-6 quartic. With a non-fine-tuned lambda, ONLY a coupling-suppressed mu^2 ~ alpha_0^16")
    print("    M_P^2 works -> v is a RADIATIVE (Coleman-Weinberg / dimensional-transmutation) scale, like")
    print("    Lambda_QCD. So v is NOT an arbitrary 2nd anchor; it is the radiative EW scale.")
    print("  - the 'numerology' 2.50*alpha_0^8 is, radiatively, mu^2 ~ alpha_0^16 M_P^2 with the prefactor")
    print("    = 1/sqrt(lambda) => lambda ~ 0.16, the SAME ballpark as the measured Higgs quartic 0.13.")
    print("    The exponent n/2 = 8 is suggestively the [8,4,4] byte. Both the power AND the prefactor get")
    print("    a sensible radiative reading -- a real sign the reframing is the right physics.")
    print("  - NO-GO-1 reframed: the open question is no longer 'why a 2nd scale' but 'derive the power")
    print("    n=16 (the byte/loop structure) and the quartic lambda' -- Phases 2-3.")
    print("  SCOPE: schematic (O(1) Coleman-Weinberg factors uncomputed); does NOT derive n or lambda;")
    print("  lambda~0.16 vs 0.13 is a ~24% ballpark, not a precision match. Phase 1 = the forcing + the")
    print("  consistent reframing, the high-confidence piece. The number stays Phases 2-3 (input-limited).")

    # gates
    assert abs(math.log10(v_from_MP / v)) > 15, "mu^2~M_P^2 must miss by >15 orders (no hierarchy)"
    assert lam_needed_B < 1e-5, "mu^2~Lambda^2 must need an absurdly tiny quartic (excluded)"
    assert 15.0 < n_at_higgs < 17.0, "radiative power n must be ~16 (= 2x8 byte) at the Higgs quartic"
    assert 0.5 < lam_implied / lam_Higgs < 2.0, "lambda implied by the alpha_0^8 fit must be ~lambda_Higgs (ballpark)"
    print("\nGATES PASSED -- O(1) tree options excluded; radiative mu^2~alpha_0^16 M_P^2 forced; the alpha_0^8")
    print("prefactor maps to a Higgs-ballpark quartic (lambda~0.16 vs 0.13). Radiative reframing established;")
    print("n=16 and lambda remain to derive (Phases 2-3). exit 0")


if __name__ == "__main__":
    main()
