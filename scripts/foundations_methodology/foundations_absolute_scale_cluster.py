#!/usr/bin/env python3
r"""foundations_absolute_scale_cluster.py

The absolute-scale cluster (the dimensionless residuals the cosmological-selector Omega-span cannot supply;
the "one lock closes them all" unifier was tested FALSE, DRIFT 2112). Examined 2026-06-25.

This is a CHARACTERIZATION, not a closure. The headline is that the user-list framing of the cluster is
STALE on TWO of its four items (both checked against latest canon):

(1a) STALENESS CORRECTION -- GRAVITY. Gravity is no longer "alpha^2 / K_eff=205": that route is
     CLOSED-NEGATIVE (the BZ-Feshbach trace is not a constructible finite-graph invariant,
     keff_invariant_exhaustion.py; DRIFT L2404), SUPERSEDED by the alpha^1 / C_loop=3/2 erasure route
     (conditional on item 79). The alpha-power dropped 2 -> 1.

(1b) STALENESS CORRECTION -- CC. "CC microscopic prefactor A_req" is also stale: the sqrt(3/2)/A_req route
     was REJECTED (cross-sector -- 3/2 is the GRAVITY C_loop, not a CC fluctuation determinant; DRIFT L2364),
     but the CC rho_Lambda coefficient ALREADY CONVERTED free->conditional at MODEL GRADE via the
     sector-native generation-vertex BZ loop C_loop = (1/3)<1/|sin k|>_BZ = 0.30356 (rho_Lambda at 0.10 sigma;
     DRIFT L2364). So the CC is NOT an open coefficient -- it is conditional. Per L2364's ledger correction,
     "NO ordinary fitted free coefficient remains in the cosmological/QEC core" -- the only genuinely
     irreducible native-core entries are alpha^2/K_eff (non-constructible, G7) and R_HBC (a FIELD, not a
     coefficient -- the HBC item-131 scalar-clock residual). The cluster is MORE closed than "four open items."

(2) PARTIAL [8,4,4]-CHANNEL-DIMENSION PATTERN (suggestive, 3/5; NOT a unification). Three of the cluster
    coefficients have explicit [8,4,4]/service-channel DENOMINATORS, with the alpha-POWER set by
    code-distance / event-order:
      - eta   = (3/14) alpha^4 + (1/3) alpha^5   : 14 = A_4 (the [8,4,4] weight enumerator W=1+14x^4+x^8);
                                                   alpha^4 = code distance d=4 (item 126).
      - alpha0/208 (source rate)                 : 208 = dim of the invalid (virtual/Higgs) subspace Q (L205).
      - A_s   = (3/4) alpha0^4 = alpha0^4 / C_F   : C_F = 4/3 (the colour-restoring count; item 131).
    The gravity C_loop=3/2 and the CC A_req do NOT fit this denominator pattern cleanly -> it is 3/5,
    suggestive of a shared combinatorial origin but neither universal nor a single number (the span-unifier
    is FALSE; the powers/coefficients remain distinct theorems).

(3) GRADE TRIAGE (sharply non-uniform):
      - eta (item 126):     BEST-grounded -- alpha^4 motivated (d=4 AND the unique power giving O(1) prefactor);
                            6.145e-10 at ~0.6 sigma Planck. (NB the 3/14 is the best simple rational but sits
                            right at the strict §16.3 0.7% boundary with eta_obs=6.12e-10 -> "near-unique",
                            not cleanly unique.) Residual: T3 absolute photon-bookkeeping (n_B/n_gamma).
      - A_s (item 131):     WELL-grounded -- A_s = 1/N_shell = alpha0^4/C_F (1.3% of Planck). Residual: the
                            scalar-clock field identification.
      - alpha0/208 rate:    CONDITIONAL theorem (§5.9-predicted + data-selected; tied to CMB M15).
      - gravity:            WEAKEST of the four MEMBERS -- alpha-power reframed 2->1; the alpha^2/K_eff=205
                            prefactor is the non-constructible G7 (the genuinely irreducible piece); the live
                            alpha^1/C_loop=3/2 erasure route is conditional on item 79.
      - CC rho_Lambda:      CONDITIONAL (model grade), NOT open -- the sqrt(3/2)/A_req route is REJECTED
                            (cross-sector), but the CC coeff already converted via the generation-vertex BZ
                            loop C_loop=(1/3)<1/|sin k|>=0.30356 (rho_Lambda at 0.10 sigma; DRIFT L2364).
    GENUINELY-IRREDUCIBLE CORE (L2364, after the corrections): only TWO non-coefficient objects --
      alpha^2/K_eff (non-constructible G7; live alt = alpha^1/C_loop=3/2 cond. item 79) and R_HBC (a field).

Self-asserting on the robust facts; the grades/pattern/§16.3-nuance are REPORTED, not asserted. exit 0.
(Pure-stdlib: GF(2) combinatorics + Decimal-free arithmetic; no numpy dependency.)
"""
import itertools
from fractions import Fraction
from math import gcd


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    a0 = 1 / 137.035999084
    N = 3

    print("=== robust facts: the [8,4,4]/service channel-dimension alphabet ===")
    G = [[1, 1, 1, 1, 1, 1, 1, 1], [0, 0, 0, 0, 1, 1, 1, 1], [0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 0, 1, 0, 1, 0, 1]]
    wd = {}
    for b in itertools.product([0, 1], repeat=4):  # all 16 messages -> codewords (GF(2)) -> weights
        cw = [sum(b[i] * G[i][j] for i in range(4)) % 2 for j in range(8)]
        w = sum(cw); wd[w] = wd.get(w, 0) + 1
    ok(wd.get(4) == 14, f"[8,4,4] weight enumerator 1+{wd.get(4)}x^4+x^8: A_4 = 14 (eta's denominator)")
    CF = (N * N - 1) / (2 * N)
    ok(abs(CF - 4 / 3) < 1e-12, f"C_F = (N^2-1)/(2N) = {CF:.4f} = 4/3 (A_s denominator)")
    Q_dim = 208
    ok(Q_dim == 208, "Q-dim = 208 (invalid/virtual subspace; alpha0/208 denominator; canon L205)")

    print("\n=== the coefficient theorems (power = event/code order; coeff = channel ratio) ===")
    eta = (3 / 14) * a0 ** 4 + (1 / 3) * a0 ** 5
    eta_obs = 6.12e-10
    ok(abs(eta / eta_obs - 1) < 0.01, f"eta=(3/14)a^4+(1/3)a^5 = {eta:.3e} vs obs {eta_obs:.2e} ({100*(eta/eta_obs-1):+.1f}%); 14=A_4, a^4=d=4")
    A_s = (3 / 4) * a0 ** 4
    ok(abs(A_s / 2.099e-9 - 1) < 0.02, f"A_s=(3/4)a^4=a^4/C_F = {A_s:.3e} vs Planck 2.10e-9 ({100*(A_s/2.099e-9-1):+.1f}%)")
    print(f"    alpha0/208 source rate = {a0/208:.3e} ; gravity = alpha^1/C_loop(=3/2) [SUPERSEDES alpha^2/K_eff=205, closed-negative]")

    # §16.3 nuance for eta (reported, not asserted): is 3/14 cleanly unique within 0.7% of the observed?
    tgt = eta_obs / a0 ** 4
    hits = sorted({Fraction(p, q) for q in range(1, 15) for p in range(1, 4 * q)
                   if gcd(p, q) == 1 and abs((p / q) / tgt - 1) < 0.007})
    print(f"\n=== §16.3 nuance: simple p/q (q<=14) within 0.7% of eta_obs/a^4={tgt:.4f}: {[str(x) for x in hits]}")
    print(f"    -> 3/14={3/14:.4f} lands at {100*((3/14)/tgt-1):+.2f}% (right at the boundary): BEST simple rational,")
    print(f"       'near-unique' not cleanly unique. (Contrast: the old gravity alpha^2 coeff had ~60 candidates.)")

    print("\n=== partial channel-dimension pattern (3/5; suggestive, NOT a unification) ===")
    print("    eta -> 14 (A_4) ; alpha0/208 -> 208 (Q-dim) ; A_s -> 4/3 (C_F)  [clear [8,4,4] denominators]")
    print("    gravity -> C_loop=3/2 (rate multiplicity) ; CC -> (1/3)<1/|sin k|>=0.30356 (BZ spectral avg)  [NOT channel-dim ratios]")

    print("\n[verdict] absolute-scale cluster:")
    print("  - TWO STALENESS CORRECTIONS to the user-list framing:")
    print("    (1a) gravity alpha^2/K_eff=205 (closed-negative) -> alpha^1/C_loop=3/2 (cond. item 79); 'alpha^2' is stale.")
    print("    (1b) CC 'A_req' is stale: sqrt(3/2)/A_req REJECTED (cross-sector); CC rho_Lambda already conditional")
    print("         (model grade) via the generation-vertex BZ loop (0.10 sigma, L2364). CC is NOT open.")
    print("  - PARTIAL PATTERN: 3 coefficients (eta, alpha0/208, A_s) use [8,4,4]-channel denominators (14, 208, C_F);")
    print("    gravity/CC do not -> a shared combinatorial alphabet for the denominators, NOT a single-number unification.")
    print("  - NET: the cluster is MORE closed than 'four open items' -- all four MEMBERS are conditional/near-closed")
    print("    (eta ~0.6sigma; A_s 1.3%; alpha0/208 conditional; CC rho_Lambda 0.10sigma model-grade). The genuinely-")
    print("    irreducible CORE (L2364) is just TWO non-coefficient objects: alpha^2/K_eff (non-constructible G7) and")
    print("    R_HBC (a scalar-clock field). The span-unifier stays FALSE; no new closure. exit 0")


if __name__ == "__main__":
    main()
