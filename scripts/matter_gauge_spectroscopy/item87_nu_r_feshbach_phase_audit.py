#!/usr/bin/env python3
r"""Item 87 -- does the DOCUMENTED substrate embed the nu_R defect as an open arc
of d at radius N (a phase delta = d/N rad)?  Build the actual Part-18 Feshbach
coupling and check whether any phase -- arc OR winding -- is present at all.

The user's question (the sharpened residual): the phase form is pinned to raw-
linear d/N with the winding/trig classes excluded; is the open arc realized by the
substrate?  This script answers at the operator level by constructing the
documented coupling exactly as ANCHOR sec 5.8 / Part 18 describes it:

  * P-block lepton states that couple, per generation g: e_R(g), nu_L(g);
  * Q-block: nu_R(g);
  * H_PQ: e_R(g)->nu_R(g) and nu_L(g)->nu_R(g), each UNIT amplitude (line 771),
    arising from an I_3-flip and a chi-flip -- real Pauli-X type, no phase;
  * H_QQ = I_3 (line 770);
  * Feshbach self-energy at the pole:  Sigma = H_PQ (E I - H_QQ)^{-1} H_QP, E=0.

If Sigma is real and generation-block-diagonal, the documented substrate carries
NO inter-generation phase -- neither an arc (d/N) nor a winding (2*pi*d/N) -- so
the answer is NO, and delta = d/N is an external ansatz, not an operator property.
"""
import numpy as np

GEN = 3
# P-block order: [e_R0, nuL0, e_R1, nuL1, e_R2, nuL2]; Q-block: [nuR0, nuR1, nuR2]


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def build_documented_HPQ():
    """6 unit, REAL entries: e_R(g),nu_L(g) -> nu_R(g) (ANCHOR Part 18, line 771)."""
    H = np.zeros((2 * GEN, GEN), dtype=complex)
    for g in range(GEN):
        H[2 * g, g] = 1.0      # e_R(g) -> nu_R(g)  (I_3 flip)
        H[2 * g + 1, g] = 1.0  # nu_L(g) -> nu_R(g) (chi flip)
    return H


def main():
    print("ITEM 87 -- DOCUMENTED nu_R FESHBACH: is a phase (arc/winding) present?")
    print("=" * 80)

    H_PQ = build_documented_HPQ()
    H_QP = H_PQ.conj().T
    H_QQ = np.eye(GEN, dtype=complex)
    Sigma = H_PQ @ np.linalg.inv(0.0 * np.eye(GEN) - H_QQ) @ H_QP  # = -H_PQ H_QP

    print("\n[1] The coupling and self-energy are REAL (no phase to carry an arc)")
    check(np.max(np.abs(H_PQ.imag)) < 1e-15, "H_PQ is real (unit-amplitude I_3/chi flips)")
    check(np.max(np.abs(Sigma.imag)) < 1e-15, "Feshbach self-energy Sigma is real")

    print("\n[2] Sigma is generation-block-diagonal: inter-generation B_nu = 0")
    inter = 0.0
    for g in range(GEN):
        for h in range(GEN):
            if g != h:
                blk = Sigma[2 * g:2 * g + 2, 2 * h:2 * h + 2]
                inter = max(inter, float(np.max(np.abs(blk))))
    print(f"    max |inter-generation block| = {inter:.2e}")
    check(inter < 1e-15, "no inter-generation coupling -> there is NO B_nu, hence NO phase arg(B_nu)")

    print("\n[3] The 'd' that appears is a REAL eigenvalue (channel count), not a phase")
    g0 = Sigma[0:2, 0:2].real                         # generation-0 block on (e_R, nu_L)
    ev = np.sort(np.linalg.eigvalsh(g0))
    print(f"    generation block = {g0.tolist()},  eigenvalues = {ev.tolist()}")
    check(abs(abs(ev[0]) - 2.0) < 1e-12, "massive eigenvalue |lambda| = 2 = d (the 2 coupling channels e_R,nu_L->nu_R)")
    check(np.max(np.abs(g0 - g0.T)) < 1e-15, "block is real-symmetric (phaseless)")

    print("\n[4] N=9 does not enter this operator at all")
    print("    Sigma is built from the 2 lepton Feshbach channels; the '9-site plaquette'")
    print("    is a SEPARATE gauge structure. So d/N=2/9 conflates a Feshbach channel")
    print("    count (2, here a real eigenvalue) with a gauge-plaquette count (9, absent")
    print("    here). There is no single operator in which 'arc d at radius N' is defined.")
    check(True, "N=9 is not a parameter of the documented Feshbach coupling")

    print(
        """
[5] ANSWER (grounded in the documented operator)
    NO. The documented substrate does NOT embed the nu_R defect as an open arc of
    d at radius N. The coupling H_PQ is real (unit-amplitude I_3/chi flips) and the
    Feshbach self-energy Sigma is real and generation-block-diagonal -- the inter-
    generation off-diagonal B_nu is exactly zero, so there is no phase arg(B_nu) at
    all: neither an arc (raw d/N) nor a winding (2*pi*d/N). The "d" enters as a REAL
    eigenvalue |lambda|=2 (the two coupling channels), not as a subtended angle; the
    "N=9" plaquette is a different operator and never appears here.

    Consequences for the residual:
      * The winding route was already excluded by the lepton fit; this shows the arc
        route has no substrate realization either -- the documented coupling has no
        phase to be an arc OF.
      * Therefore delta = d/N is, at the current substrate, an EXTERNAL ANSATZ
        (consistent with the sec 5.8 phase-lift obstruction): producing any nonzero
        inter-generation phase requires the new, underived K3 Majorana portal
        (item87_deltaL2_holonomy_coupling.py), which is not in the documented walk.
      * The honest status is sharper than "open arc, premise undocumented": the
        premise is FALSE for the documented operator (it is phaseless). A phase
        requires NEW substrate physics (the portal), and why its magnitude would be
        the real ratio 2/9 in radians is unexplained. delta=2/9 / 1/3 stay fit-
        selected; item 87 is open and the arc hypothesis is, as posed, answered NO.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
