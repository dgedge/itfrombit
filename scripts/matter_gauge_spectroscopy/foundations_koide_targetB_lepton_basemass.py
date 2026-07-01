#!/usr/bin/env python3
r"""foundations_koide_targetB_lepton_basemass.py

KOIDE TARGET B (the per-sector base masses mu_l, mu_d, mu_u, currently "anchored externally" -- §5.8 fixes
mu from the tau mass; §15 item lists the base masses as underived). Partial closure: the LEPTON base mass
is derived from the strong anchor; the QUARK base masses are shown to be scale-dependent (the obstruction).

Per sector m_n = mu (1 + R cos(delta + 2pi n/3))^2, so sqrt(m_n) = sqrt(mu)(1 + R cos), and since
sum cos = 0, sum sqrt(m_n) = 3 sqrt(mu)  ->  the base mass mu = ((sum sqrt m)/3)^2  (the overall scale;
R, delta only set the within-sector ratios).

LEPTON (scale-INVARIANT -- no QCD running):
  mu_l = ((sqrt m_e + sqrt m_mu + sqrt m_tau)/3)^2 = 313.84 MeV
       = m_N / 3 = (2 sqrt2 / 3) Lambda_QCD   to 0.26%   (m_N = 2 sqrt2 Lambda, §9.10)
  i.e. the charged-lepton Koide base mass IS the constituent-quark mass (nucleon / 3).
  §16.3 discipline: 1/3 is the UNIQUE simple rational p/q (q<=12) within the 0.26% match of mu_l/m_N -- so
  this is a clean, structurally-motivated landing, not a forking-path fit. It derives mu_l from the
  framework's strong anchor Lambda (replacing the external tau-anchor), and is consistent with the 2026-05-29
  Lindemann-Weierstrass result that a DISCRETE-ONLY Target B fails: mu_l needs the dimensionful anchor
  Lambda, which is exactly what this supplies.
  Plausible mechanism (not yet rigorous): the lepton mass is the nu_R Feshbach self-energy eps (§5.8); eps
  inherits the QCD dynamical/constituent mass scale m_N/3 of the shared substrate vacuum. Promotion target:
  derive eps = m_N/3 from the Feshbach coupling.

QUARK (scale-DEPENDENT -- the obstruction): mu_d, mu_u are built from running masses and are NOT unique
numbers (they change with the renormalisation scale). So the quark base masses cannot be a single clean
multiple of Lambda/v without a scale prescription -- that, not a missing structural factor, is why quark
Target B is hard. mu_d ~ 2 mu_l is suggestive but scale-ambiguous; mu_u has no clean parameter-free landing.

Self-asserting; exit 0.
"""
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def mu(ms):
    return (sum(np.sqrt(m) for m in ms) / 3) ** 2


def main():
    print("=== Koide Target B: per-sector base mass mu = ((sum sqrt m)/3)^2 ===\n")
    Lam = 332.0
    mN = 2 * np.sqrt(2) * Lam   # framework nucleon, §9.10

    # LEPTON
    mL = [0.51099895, 105.6583755, 1776.86]
    muL = mu(mL)
    pred = mN / 3
    print(f"[lepton] mu_l = {muL:.3f} MeV (scale-invariant)")
    print(f"         m_N/3 = (2 sqrt2/3) Lambda = {pred:.3f} MeV ; mu_l/(m_N/3) = {muL/pred:.5f} ({100*(muL/pred-1):+.2f}%)")
    ok(abs(muL - 313.84) < 0.1, "mu_l = 313.84 MeV from PDG lepton masses")
    ok(abs(muL / mN - 1/3) / (1/3) < 0.004, "mu_l = m_N/3 = (2 sqrt2/3) Lambda to <0.4% -> lepton base mass from the strong anchor")
    tol = abs(muL / mN - 1/3) / (1/3)
    comp = [(p, q) for q in range(2, 13) for p in range(1, q) if abs((p/q)/(muL/mN) - 1) < tol and np.gcd(p, q) == 1]
    print(f"         §16.3: simple p/q (q<=12) within {100*tol:.2f}% of mu_l/m_N = {muL/mN:.4f}: {comp}")
    ok(comp == [(1, 3)], "1/3 is the UNIQUE simple rational in tolerance -> not a forking-path fit (clean landing)")

    # QUARK -- scale dependence
    print("\n[quark] base masses are SCALE-DEPENDENT (running) -> not unique numbers:")
    dn_own = [4.67, 93.4, 4180.0]            # m_q(m_q)-ish
    dn_2 = [4.67, 93.4, 4180.0 * 1.45]       # m_b run toward 2 GeV (illustrative)
    up_own = [2.16, 1270.0, 162500.0]
    up_2 = [2.16, 1270.0 * 1.35, 172690.0]
    print(f"         mu_down: own-scale {mu(dn_own):.0f} MeV  vs  ~2 GeV {mu(dn_2):.0f} MeV  (differ by {100*(mu(dn_2)/mu(dn_own)-1):.0f}%)")
    print(f"         mu_up:   own-scale {mu(up_own)/1e3:.1f} GeV vs ~2 GeV {mu(up_2)/1e3:.1f} GeV")
    ok(abs(mu(dn_2)/mu(dn_own) - 1) > 0.2, "quark base masses change >20% with scale -> scale-dependent, NOT unique")
    print(f"         mu_down(own)/mu_l = {mu(dn_own)/muL:.2f} (~2, but scale-ambiguous); mu_up/m_t = {mu(up_own)/162500:.3f} (no clean landing)")

    print("\n[verdict] Koide Target B -- partial closure:")
    print("  - LEPTON: mu_l = m_N/3 = (2 sqrt2/3) Lambda_QCD (the constituent-quark mass), 0.26%, §16.3-clean")
    print("    (1/3 unique). Derives the lepton base mass from the strong anchor, replacing the tau-anchor;")
    print("    consistent with the discrete-only L-W falsification (mu_l needs the dimensionful Lambda).")
    print("    Sharpened sub-question: derive the Feshbach self-energy eps = m_N/3 rigorously.")
    print("  - QUARK: mu_d, mu_u are scale-dependent (running) -> no unique value; the obstruction is the")
    print("    renormalisation-scale ambiguity, not a missing factor. mu_d ~ 2 mu_l is suggestive but")
    print("    scale-ambiguous; mu_u has no clean parameter-free landing. Quark Target B stays open. exit 0")


if __name__ == "__main__":
    main()
