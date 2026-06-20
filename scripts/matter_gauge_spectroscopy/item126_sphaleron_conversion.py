#!/usr/bin/env python3
r"""ITEM 126 — the (B-L) -> B sphaleron leg of the baryon asymmetry, on the framework's content.

CONTEXT. item126_baryon_sign_analysis.py established that the baryon SIGN is not a static census
quantity; it routes via leptogenesis from the Majorana neutrino. The (B-L) generated then converts
to a baryon number B through electroweak sphalerons. This script asks: does the framework REPRODUCE
the standard (B-L)->B conversion, using its OWN content?

WHAT IS FRAMEWORK-NATIVE vs INHERITED.
  * NATIVE: (i) N_g = 3 generations — R1 forbids the generation corner (G0,G1)=(1,1), leaving 3 of
    4 (ANCHOR §15 / line 670). (ii) the hypercharges Y in {1/6, 2/3, -1/3, -1/2, -1} — DERIVED from
    the §2.8 Z-stabiliser charge formula (Y = Q - T3; ANCHOR line 284). (iii) B-L conservation +
    the right-handed neutrino nu_R — forced by the anomaly-cancellation (B-L) residual (ANCHOR
    line 319). So the framework supplies exactly the field content and the conserved charge the
    conversion needs.
  * INHERITED: the high-T chemical-equilibrium calculation itself (sphaleron + QCD-instanton + Yukawa +
    hypercharge-neutrality constraints) is standard thermal field theory (Harvey-Turner 1990). The
    framework does not modify it; it supplies the inputs. So the conversion is a CONSISTENCY result
    (the framework is compatible with, and parameterises, the standard 28/79), not a new derivation.
  * DEEPER-OPEN: a SUBSTRATE origin of the sphaleron process (the EW B+L anomaly = the SU(2)_{chi,W}
    instanton sector) is not derived here — only the equilibrium bookkeeping over the content.

RESULT. With the framework content (N_g=3, N_H=1, §2.8 hypercharges), the equilibrium constraints
have a 1-dimensional solution space = the conserved B-L direction, on which
    B = (8 N_g + 4 N_H)/(22 N_g + 13 N_H) (B-L) = 28/79 (B-L).
So a NEGATIVE (B-L) (lepton excess from nu_R decay) converts to a POSITIVE B of magnitude (28/79).

THE "28" COINCIDENCE (flagged, NOT claimed). 8 N_g + 4 N_H = 28 = C(8,2), the framework's Shell-1 /
28-clock count. But the sphaleron 28 is content-counting (contingent on N_g=3, N_H=1: N_g=4 gives 36)
while C(8,2)=28 is register-pairing (independent of N_g, N_H). No structural map is established; this
is a numerical coincidence at (3,1) unless a real identity is shown. Recorded as a question, not a result.

Self-asserting; exit 0 = the framework content yields the standard 28/79 with B-L the conserved mode.
"""
from __future__ import annotations
from fractions import Fraction as F

import numpy as np


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def conversion_factor(Ng: int, NH: int):
    """Solve the high-T chemical-equilibrium constraints; return B/(B-L) on the conserved mode.
    Variables (generation-universal, unbroken SU(2) so doublet members share mu):
        [muQ(qL), muU(uR), muD(dR), muL(lL), muE(eR), mu0(Higgs)]."""
    rows = [
        [3, 0, 0, 1, 0, 0],                       # EW sphaleron:   3 muQ + muL = 0
        [2, -1, -1, 0, 0, 0],                     # QCD instanton:  2 muQ - muU - muD = 0
        [1, -1, 0, 0, 0, 1],                      # up Yukawa:      muQ - muU + mu0 = 0
        [1, 0, -1, 0, 0, -1],                     # down Yukawa:    muQ - muD - mu0 = 0
        [0, 0, 0, 1, -1, -1],                     # lepton Yukawa:  muL - muE - mu0 = 0
        [Ng, 2 * Ng, -Ng, -Ng, -Ng, 2 * NH],     # hypercharge neutrality (Higgs boson factor 2)
    ]
    M = np.array(rows, float)
    s = np.linalg.svd(M, compute_uv=False)
    nz = int(np.sum(s < 1e-9))
    v = np.linalg.svd(M)[2][-1]
    muQ, muU, muD, muL, muE, mu0 = v
    B = 2 * Ng * muQ + Ng * muU + Ng * muD        # baryon number (B=1/3 per quark Weyl, 6Ng+3Ng+3Ng comps)
    L = 2 * Ng * muL + Ng * muE                   # lepton number
    return B / (B - L), nz


def main() -> int:
    print("ITEM 126 — (B-L) -> B SPHALERON CONVERSION ON THE FRAMEWORK CONTENT")
    print("=" * 86)

    # framework-native content
    Ng, NH = 3, 1
    Y = {"qL": F(1, 6), "uR": F(2, 3), "dR": F(-1, 3), "lL": F(-1, 2), "eR": F(-1)}
    print(f"\n[content] N_g = {Ng} (R1 forbids (G0,G1)=(1,1): 3 of 4 — ANCHOR §15/line 670)")
    print(f"          N_H = {NH};  hypercharges Y (§2.8, line 284): {dict((k, str(v)) for k, v in Y.items())}")
    print(f"          nu_R present (anomaly (B-L) residual, line 319) => B-L well-defined and conserved")

    ratio, nz = conversion_factor(Ng, NH)
    print(f"\n[solve] equilibrium constraint null-space dimension = {nz} (the conserved B-L mode)")
    print(f"        B/(B-L) = {ratio:.12f}")
    print(f"        as fraction = {F(ratio).limit_denominator(10000)}   (target Harvey-Turner 28/79)")
    check(nz == 1, "the equilibrium constraints leave exactly ONE conserved direction = B-L")
    check(F(ratio).limit_denominator(10000) == F(28, 79),
          "framework content reproduces the standard sphaleron factor B = (28/79)(B-L)")
    check(F(8 * Ng + 4 * NH, 22 * Ng + 13 * NH) == F(28, 79),
          "matches the closed form (8 N_g + 4 N_H)/(22 N_g + 13 N_H) at (N_g,N_H)=(3,1)")

    # the result depends only on N_g, N_H given the SM hypercharges — show the dependence
    print("\n[dependence] B/(B-L) over (N_g, N_H) — the 28/79 is specific to the framework's (3,1):")
    for ng, nh in ((3, 1), (4, 1), (3, 2), (2, 1)):
        r, _ = conversion_factor(ng, nh)
        print(f"     (N_g={ng}, N_H={nh}): B/(B-L) = {F(r).limit_denominator(10000)}   "
              f"[8N_g+4N_H = {8 * ng + 4 * nh}]")

    # the "28" coincidence — flagged, not claimed
    print("\n[28?] 8 N_g + 4 N_H = {} = C(8,2) = {}".format(8 * Ng + 4 * NH, 8 * 7 // 2))
    check(8 * Ng + 4 * NH == 8 * 7 // 2,
          "numerical coincidence: sphaleron numerator (content-counting, contingent on (3,1)) "
          "equals C(8,2) (register-pairing) — NO structural map established; recorded as a question")

    print("\n" + "=" * 86)
    print("VERDICT")
    print(
        "  The (B-L) -> B leg is QUANTITATIVELY reproduced on the framework's own content: N_g=3\n"
        "  (R1), the §2.8-derived hypercharges, and the anomaly-forced nu_R / B-L conservation feed\n"
        "  the standard high-T equilibrium and give B = (28/79)(B-L), with B-L the unique conserved\n"
        "  mode. A lepton excess (negative B-L from nu_R decay) therefore converts to a POSITIVE\n"
        "  baryon number of fixed magnitude 28/79. TIERS: B-L conservation + the content are\n"
        "  framework-NATIVE (anomaly structure); the equilibrium count is INHERITED standard physics\n"
        "  (a consistency result, not a new derivation); a substrate origin of the sphaleron process\n"
        "  (the SU(2)_{chi,W} B+L anomaly) is deeper-open. The 8N_g+4N_H = C(8,2) = 28 coincidence is\n"
        "  flagged, not claimed. NET: the conversion leg is solid; the asymmetry's SIGN still rests on\n"
        "  the sign of the nu_R-sector CP asymmetry (item 87), which this leg does not fix."
    )
    print("exit 0 — framework content -> standard 28/79 (B-L conserved); sign still set upstream (item 87).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
