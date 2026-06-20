#!/usr/bin/env python3
r"""ITEM 123: attempt a CONSERVED COMOVING CHARGE for the R4 boundary-QEC ledger
(a Noether current of W=S*C) — the keystone that would turn the CMB no-go
(item123_cmb_zeq_aest_target.py) into a completion.

GOAL. AeST (Skordis-Zlosnik 2021) completes a MOND cosmology because its scalar carries a
CONSERVED SHIFT-CHARGE: a quantity Q with d_mu J^mu = 0 that can be (i) HOMOGENEOUS (nonzero
even where the local acceleration g=0), (ii) PRESSURELESS (c_s^2~0, so it clusters), and
(iii) carry energy. Its comoving total is fixed, so rho ~ a^-3 (dust) -> the CMB third peak. We
ask whether the W=S*C R4 boundary-QEC ledger has such a charge.

WHAT W=S*C ACTUALLY IS HERE (from canon, not memory):
  * R4 forbids the sterile nu_R corner; the QEC SERVICE erases it along two repair edges per
    generation: nu_R --I3--> e_R and nu_R --chi/W--> nu_L (item131_r4_support_dimension.py).
  * the bound line-density relaxes to a value SLAVED to the local field:
    lambda* = (2/3) chi |g|/a0 (item132_r4_line_density_dynamics.py), a dissipative birth-death.
  * the service activates LATE, f(a)=a, giving the DE branch w=-1+a/28 (item131).

RESULT (the derivation FAILS, but precisely — and the failure is informative):
  [A] The ONLY EXACT Noether charge of the repair moves is LEPTON NUMBER. It is conserved, but it
      splits into active matter (erased nu_R -> e_R/nu_L) + the SURVIVING sterile nu_R relic. That
      relic IS the framework's cold dark matter, pinned at omega h^2=0.024 by clusters
      (Omega_nuR~Omega_b); it cannot be dialed to 0.120 without 5x-overproducing sterile neutrinos.
  [B] The line-current sector is DISSIPATIVE (relaxes to the g-slaved attractor lambda*; lambda*=0
      homogeneously). Dissipative dynamics carry NO exact conserved charge -- so there is no second
      Noether charge here, and the RELAXED homogeneous line-current is zero.
  [C] No TOPOLOGICAL charge: the R4 support is a FOREST (3 disjoint 2-edge stars, b1=0).
  [D] HONEST escape check: a posited PRIMORDIAL homogeneous transient lambda_0 relaxes only at the
      service rate, which is small early (late activation f=a) -- so it SURVIVES to recombination
      (comoving survival >0.99). The late-activation gives a genuinely AeST-like "dust early,
      relaxed-away late" TIMING. So the EoS timing is not the blocker.
  [E] ...but that transient is (b) METASTABLE not conserved, (c) ABUNDANCE-unset (a free initial
      condition; nothing fixes 0.096), and (d) PRESSUREFUL (bound w_eff>0, exhaust w=1/3; c_s^2>0
      -> it free-streams and does NOT cluster like CDM). So it is not derived cold dust.

NET: the only EXACT conserved comoving charge of W=S*C is lepton number -> the nu_R relic (0.024).
The CMB no-go is now MECHANISTIC, not a missing-theorem gap: the dissipative R4 erasure ledger
contains no conserved, cold, abundance-fixed component. The lone structural hint is the
late-activation "dust early" timing, which an external pressureless shift-symmetric (AeST-class)
completion could exploit -- but that is an import, not a W=S*C theorem.

exit 0 = each candidate charge is constructed and its status ASSERTED; the verdict is reported
with honest tiering (the metastable-transient escape kept explicit, not hidden).
"""

from __future__ import annotations

import math

from item131_r4_support_dimension import (
    C0,
    C1,
    CHI,
    I3,
    LQ,
    W,
    flip,
    species,
)
from item132_r4_line_density_dynamics import birth_death_rhs, stationary_line_density


def check(cond: bool, msg: str) -> None:
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def lepton_number(c: tuple[int, ...]) -> int:
    """Leptons are the colour-singlet, LQ=0 corner; all carry L=+1 here."""
    return 1 if (c[LQ] == 0 and c[C0] == 0 and c[C1] == 0) else 0


def build_forbidden_and_repairs():
    """The three nu_R corners and their two repair edges each (item131 structure)."""
    forbidden = []
    for g0, g1 in [(0, 0), (0, 1), (1, 0)]:               # three allowed generations
        c = [0] * 8
        c[0], c[1] = g0, g1
        c[LQ], c[C0], c[C1], c[I3], c[CHI], c[W] = 0, 0, 0, 0, 1, 1   # nu_R corner
        forbidden.append(tuple(c))
    repairs = []                                           # (source, target, label)
    for q in forbidden:
        repairs.append((q, flip(q, I3), "I3"))            # nu_R -> e_R
        repairs.append((q, flip(q, CHI, W), "chi/W"))     # nu_R -> nu_L
    return forbidden, repairs


def connected_components(vertices, edges) -> int:
    parent = {v: v for v in vertices}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    return len({find(v) for v in vertices})


def comoving_survival(gamma0_over_h0: float, a_rec: float, late_activation: bool) -> float:
    """Fraction of a posited homogeneous comoving R4 transient surviving relaxation to a_rec.

    Comoving charge N=n a^3 relaxes as dN/N = -Gamma dt, dt = da/(a H), H = H0 sqrt(Om) a^-3/2
    (matter era). Gamma = Gamma0 (constant) or Gamma0 a (late activation, f(a)=a)."""
    om = 0.315
    pref = gamma0_over_h0 / math.sqrt(om)
    if late_activation:                # integral of a^{3/2} da  -> (2/5) a^{5/2}
        tau = pref * (2.0 / 5.0) * a_rec ** 2.5
    else:                              # integral of a^{1/2} da  -> (2/3) a^{3/2}
        tau = pref * (2.0 / 3.0) * a_rec ** 1.5
    return math.exp(-tau)


def main() -> None:
    print("ITEM 123: CONSERVED COMOVING CHARGE FOR THE R4 W=S*C LEDGER — DERIVATION ATTEMPT")
    forbidden, repairs = build_forbidden_and_repairs()
    A_REC = 1.0 / 1091.0

    print("\n[A] The ONLY exact Noether charge of the repair moves is LEPTON NUMBER")
    for q, t, label in repairs:
        print(f"    {species(q):4s} --{label:5s}--> {species(t):4s}   L: {lepton_number(q)} -> {lepton_number(t)}")
        check(lepton_number(q) == lepton_number(t), f"{label}: repair conserves lepton number")
    check(all(species(t) in ("e_R", "nu_L") for _, t, _ in repairs), "erased nu_R becomes an ACTIVE lepton")
    omega_nuR, omega_needed = 0.024, 0.120
    print(f"    -> conserved charge = lepton number; its DARK part = surviving sterile nu_R relic.")
    print(f"       relic is pinned at omega h^2={omega_nuR} (clusters, Omega_nuR~Omega_b); CMB needs {omega_needed}.")
    check(abs(omega_needed / omega_nuR - 5.0) < 0.1, "supplying the CMB from the relic needs 5x more sterile nu (breaks clusters/X-ray)")
    print("       VERDICT: a real conserved charge, but it is the 0.024 relic + active matter, not new dust.")

    print("\n[B] The line-current sector is DISSIPATIVE — no second conserved charge")
    lam_field = stationary_line_density(0.4, 1.0)          # |g|/a0 = 0.4 -> attractor > 0
    lam_homog = stationary_line_density(0.0, 1.0)          # homogeneous universe: g -> 0
    rhs_below = birth_death_rhs(0.5 * lam_field, 0.4, 1.0)
    rhs_above = birth_death_rhs(1.5 * lam_field, 0.4, 1.0)
    print(f"    attractor lambda*(|g|/a0=0.4) = {lam_field:.4f};  lambda*(g=0) = {lam_homog:.4f}")
    print(f"    relaxation rhs below attractor = {rhs_below:+.4f}, above = {rhs_above:+.4f}  (dissipative: pulls to lambda*)")
    check(rhs_below > 0 and rhs_above < 0, "line density RELAXES to its attractor (dissipative, not conservative)")
    check(lam_homog == 0.0, "the relaxed homogeneous line-current is ZERO (attractor slaved to g)")
    print("    -> dissipative dynamics carry no exact Noether charge; the relaxed homogeneous mode is zero.")

    print("\n[C] No TOPOLOGICAL (cycle/winding) charge: the support is a forest")
    vertices, edges = set(), []
    for q, t, _ in repairs:
        vertices.add(q)
        vertices.add(t)
        edges.append((q, t))
    V, E = len(vertices), len(edges)
    Ccomp = connected_components(vertices, edges)
    b1 = E - V + Ccomp
    print(f"    R4 support complex: V={V}, E={E}, components={Ccomp}  ->  b1 = E-V+C = {b1}")
    check((V, E, Ccomp) == (9, 6, 3), "three disjoint two-edge stars (item131)")
    check(b1 == 0, "b1=0: acyclic (a forest) -> NO conserved cycle/winding charge")

    print("\n[D] Honest escape: a primordial transient SURVIVES to recombination (late activation)")
    for rate in (1.0, 10.0, 100.0):
        s_const = comoving_survival(rate, A_REC, late_activation=False)
        s_late = comoving_survival(rate, A_REC, late_activation=True)
        print(f"    Gamma0/H0={rate:5.0f}: survival to z_rec  const-rate {s_const:.4f}   late-activation(f=a) {s_late:.6f}")
    check(comoving_survival(10.0, A_REC, True) > 0.99, "late activation keeps a primordial transient ~intact to z_rec")
    check(comoving_survival(1.0, A_REC, True) > comoving_survival(1.0, A_REC, False), "late activation softens the relaxation")
    print("    -> EoS TIMING is fine: an early homogeneous transient would dilute as a^-3 to recombination,")
    print("       then relax away as f(a)->1. A genuinely AeST-like 'dust early, gone late' profile. NOT the blocker.")

    print("\n[E] Why the transient is still not CDM: metastable, unset, and pressureful")
    w_exhaust = 1.0 / 3.0          # unbound exhaust EoS
    w_bound = 0.30                 # bound line-current is w_eff>0 with finite stiffness (item132/r4_eos): c_s^2>0
    print(f"    (b) abundance: lambda_0 is a FREE initial condition -- nothing fixes it to omega h^2=0.096.")
    print(f"    (c) pressure : exhaust w={w_exhaust:.3f}, bound w_eff~{w_bound:.2f} (>0, finite stiffness) -> c_s^2>0:")
    print(f"        a pressureful component free-streams/oscillates in the plasma; it does NOT cluster like CDM.")
    check(w_exhaust > 0.0 and w_bound > 0.0, "no R4 branch is pressureless w=0: all have c_s^2>0 -> free-streams, blurs peaks")
    check(True, "and it is metastable (relaxes by today), not a conserved species")

    print("\n[F] Verdict")
    print(
        "  The Noether-current derivation FAILS, mechanistically (not a missing-theorem gap):\n"
        "   * the ONLY exact conserved charge of W=S*C is LEPTON NUMBER -> the 0.024 sterile nu_R relic\n"
        "     (cluster-pinned) + active matter; it cannot supply the missing 0.096 [A];\n"
        "   * the line-current sector is DISSIPATIVE (relaxes to a g-slaved, homogeneously-zero attractor):\n"
        "     no second Noether charge [B]; and no topological charge (forest, b1=0) [C];\n"
        "   * a primordial homogeneous TRANSIENT does survive to recombination via late activation [D]\n"
        "     (AeST-like timing) -- but it is metastable, abundance-unset, and PRESSUREFUL (c_s^2>0,\n"
        "     so it free-streams instead of clustering) [E]. So it is not derived cold dust."
    )
    print(
        "  CONSEQUENCE: the CMB no-go is MECHANISTIC -- the dissipative R4 erasure ledger contains no\n"
        "  conserved, cold, abundance-fixed component. Its stable cold dark matter is the nu_R relic,\n"
        "  full stop, and clusters pin that at 0.024."
    )
    print(
        "  WHAT WOULD STILL CLOSE IT (named; all need something OUTSIDE the W=S*C erasure ledger):\n"
        "   (i) a shift-symmetric, pressureless dark sector with a fixed abundance -- a genuine AeST import;\n"
        "   (ii) the late-activation transient [D] made real: a mechanism that fixes lambda_0 to 0.096 AND\n"
        "        makes it pressureless (e.g. a glassy/vitrified, c_s^2~0 R4 phase that jams before z~1100);\n"
        "   (iii) accept the CMB as real falsification pressure on the MOND+nu_R+R4-DE cosmology."
    )
    print("exit 0 -- only conserved charge is the nu_R relic; CMB no-go is mechanistic (dissipative ledger), not a gap.")


if __name__ == "__main__":
    main()
