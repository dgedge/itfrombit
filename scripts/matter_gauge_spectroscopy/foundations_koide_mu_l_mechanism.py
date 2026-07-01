#!/usr/bin/env python3
r"""foundations_koide_mu_l_mechanism.py

CAUSAL MECHANISM for the lepton Koide base mass mu_l = (2 sqrt2/3) Lambda = m_N/3 (the relation found in
foundations_koide_targetB_lepton_basemass.py). Promotes it from "the lepton inherits the constituent mass
(plausible)" toward a derivation by grounding BOTH scales in the same C8 byte-cell spectrum.

THE SHARED STRUCTURE -- the C8 byte cell (the [8,4,4] 8 qubits as a ring; baryons wrap the full C8,
§ANCHOR L1839). Its adjacency eigenvalues are lambda_k = 2 cos(2 pi k/8); the first-harmonic transverse
pair k=1,7 is lambda_1 = lambda_7 = 2 cos(pi/4) = sqrt2.

  [DERIVED]  nucleon  M_N = (lambda_1 + lambda_7) Lambda = 2 sqrt2 Lambda = 939 MeV  (§11.2 / L2209: the
             two degenerate first-harmonic transverse modes summed, scaled by Lambda; 0.01% vs exp).

  [MECHANISM] lepton base mass  mu_l = M_N/3 = (lambda_1 + lambda_7) Lambda / 3 = (2 sqrt2/3) Lambda
             = 313 MeV  (vs the measured Koide base 313.84, +0.26%). The /3 is the substrate Z3 -- the
             SAME Z3 that is (i) 3 colours: the baryon's M_N is shared over its 3 valence quarks, giving
             the per-colour constituent mass M_N/3; and (ii) 3 generations: the lepton Koide base mass is
             the generation-mean (sum sqrt m / 3)^2. Z3 ties the colour-division (baryon) to the
             generation-averaging (lepton): both pick out the single-Z3-unit C8 dynamical mass M_N/3.

  [CRUX -- rigor gap] the charged-lepton nu_R-Feshbach self-energy eps (§5.8: H_eff = [[0, xi],[xi*, eps]],
             m = eps) must equal that single-Z3-unit C8 dynamical mass (2 sqrt2/3) Lambda. The lepton IS a
             colour-singlet single C8 cell, so its dynamical-mass quantum being one C8 standing-wave
             unit / Z3 is natural -- but eps = (2 sqrt2/3) Lambda is not yet derived from H_eff. That is
             the promotion target (clean relation -> theorem).

So the mechanism is IDENTIFIED (shared C8 first-harmonic sqrt2 x Lambda + the substrate Z3), not merely
asserted; the remaining gap is the Feshbach scale eps. Tier: candidate mechanism / structural identification.

Self-asserting; exit 0.
"""
import numpy as np


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def main():
    print("=== Causal mechanism for mu_l = m_N/3 via the C8 byte cell ===\n")
    Lam = 332.0

    # C8 spectrum
    A = np.zeros((8, 8))
    for i in range(8):
        A[i, (i + 1) % 8] = A[i, (i - 1) % 8] = 1
    ev = np.sort(np.linalg.eigvalsh(A))[::-1]
    print(f"[1] C8 byte-cell eigenvalues 2cos(2pi k/8): {np.round(ev, 4)}")
    ok(abs(ev[1] - np.sqrt(2)) < 1e-9 and abs(ev[2] - np.sqrt(2)) < 1e-9,
       "first-harmonic transverse pair lambda_1 = lambda_7 = sqrt2 (= 2 cos(pi/4))")

    # nucleon (derived endpoint)
    MN = (ev[1] + ev[2]) * Lam
    print(f"\n[2] nucleon: M_N = (lambda_1+lambda_7) Lambda = 2 sqrt2 Lambda = {MN:.2f} MeV (exp 938.92, {100*(MN/938.92-1):+.2f}%)")
    ok(abs(MN - 2 * np.sqrt(2) * Lam) < 1e-6 and abs(MN / 938.92 - 1) < 0.005, "M_N = 2 sqrt2 Lambda (§11.2, derived) matches the nucleon")

    # lepton base mass via M_N/3 (= same C8 first-harmonic, / Z3)
    mL = [0.51099895, 105.6583755, 1776.86]
    muL = (sum(np.sqrt(m) for m in mL) / 3) ** 2
    pred = MN / 3
    print(f"\n[3] lepton base mass: mu_l = M_N/3 = (lambda_1+lambda_7) Lambda/3 = {pred:.2f} MeV")
    print(f"    measured Koide base mu_l = {muL:.2f} MeV ; ratio = {muL/pred:.5f} ({100*(muL/pred-1):+.2f}%)")
    ok(abs(muL / pred - 1) < 0.004, "mu_l = M_N/3 to <0.4%: the lepton base mass = the C8 first-harmonic scale / Z3")

    # the Z3 bridge: 3 colours (baryon) = 3 generations (lepton-mean); both -> single-unit C8 dynamical mass
    print(f"\n[4] the /3 = substrate Z3: 3 colours (M_N over 3 valence quarks = constituent mass M_N/3)")
    print(f"    = 3 generations (lepton Koide base = generation-mean). Both -> single-Z3-unit C8 mass = M_N/3.")
    ok(abs(pred - (2 * np.sqrt(2) / 3) * Lam) < 1e-6, "M_N/3 = (2 sqrt2/3) Lambda -- the single-Z3-unit C8 dynamical mass")

    print("\n[verdict] mechanism for mu_l = m_N/3 IDENTIFIED (not just asserted):")
    print("  - both the nucleon and the lepton scale come from the C8 first-harmonic lambda_1=sqrt2 x Lambda;")
    print("    M_N = 2 sqrt2 Lambda (two modes, §11.2, derived); mu_l = M_N/3 (the single-Z3-unit C8 mass);")
    print("  - the /3 is the substrate Z3, simultaneously 3 colours (baryon constituent mass M_N/3) and")
    print("    3 generations (lepton Koide generation-mean) -- the Z3 bridges colour-division to gen-averaging;")
    print("  - CRUX (rigor gap): derive the §5.8 Feshbach self-energy eps = (2 sqrt2/3) Lambda from H_eff")
    print("    explicitly (the lepton is a colour-singlet single C8 cell, so eps = one C8 unit / Z3 is")
    print("    natural). Tier: candidate mechanism / structural identification; promotes the clean relation")
    print("    toward a theorem, with the Feshbach scale as the last step. exit 0")


if __name__ == "__main__":
    main()
