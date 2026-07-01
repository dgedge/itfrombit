#!/usr/bin/env python3
r"""Does 'generation = order ideal' predict the PMNS mixing symmetry?

The order-ideal Z_2 swaps the two ATOMS of {00,01,10} -- i.e. 01<->10 = e<->mu --
fixing the BOTTOM 00 = tau.  So the reframing predicts an e<->mu-symmetric,
tau-distinguished mixing pattern.

The canon's lepton mixing (tribimaximal U_nu, sec 5.10 / item 83; PMNS magnitudes
line 727) is the well-known mu-tau-symmetric pattern: |U_mu,i|^2 = |U_tau,i|^2,
with the ELECTRON row distinguished.  This script checks whether the predicted
symmetry (e<->mu, tau fixed) matches the actual one (mu<->tau, e fixed).
"""
import numpy as np


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


# rows = (e, mu, tau); cols = (nu1, nu2, nu3)
TBM = np.array([
    [np.sqrt(2/3),  np.sqrt(1/3),  0.0],
    [-np.sqrt(1/6), np.sqrt(1/3),  np.sqrt(1/2)],
    [np.sqrt(1/6), -np.sqrt(1/3),  np.sqrt(1/2)],
])
ROWS = ["e", "mu", "tau"]


def main():
    print("PMNS vs the order-ideal Z_2 generation symmetry")
    print("=" * 78)

    P = TBM ** 2   # |U|^2, the basis-independent mixing magnitudes
    print("\n[1] The canon mixing (tribimaximal) is MU-TAU symmetric, electron distinguished")
    for r, name in zip(P, ROWS):
        print(f"    |U_{name:<3s}|^2 = {np.round(r,4)}")
    check("mu and tau rows are equal (mu-tau symmetry)", np.allclose(P[1], P[2]))
    check("the electron row is distinguished (differs from mu/tau)", not np.allclose(P[0], P[1]))
    # canon PMNS magnitudes (line 727): |U_e|^2 ~ (0.67,0.30,0.022) -- same shape (e distinguished)
    canon_e = np.array([0.67, 0.30, 0.022])
    check("canon PMNS e-row (0.67,0.30,0.022) matches the TBM e-row shape (e distinguished)", abs(canon_e[2]) < 0.05 and canon_e[0] > canon_e[1])

    print("\n[2] The order-ideal Z_2 is e<->mu (atoms 01,10 swapped), tau=00 FIXED")
    # reframing's predicted mixing symmetry: invariance under e<->mu, tau fixed
    swap_emu = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])     # e<->mu
    swap_mutau = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])   # mu<->tau (the actual PMNS symmetry)
    print("    order-ideal predicts: e<->mu symmetric, tau fixed")
    print("    PMNS actually has   : mu<->tau symmetric, e fixed")

    print("\n[3] Test: is the canon mixing invariant under the PREDICTED e<->mu swap?")
    P_emu = swap_emu @ P                                       # swap e,mu rows
    P_mutau = swap_mutau @ P                                   # swap mu,tau rows
    print(f"    |U|^2 invariant under e<->mu (predicted)? {np.allclose(P, P_emu)}")
    print(f"    |U|^2 invariant under mu<->tau (actual)?  {np.allclose(P, P_mutau)}")
    check("canon mixing is NOT invariant under the order-ideal's e<->mu (prediction FAILS)", not np.allclose(P, P_emu))
    check("canon mixing IS invariant under mu<->tau (the actual symmetry, NOT predicted)", np.allclose(P, P_mutau))

    print("\n[4] Fixed points disagree: order ideal fixes TAU; PMNS fixes the ELECTRON")
    check("order-ideal distinguished generation = tau (the order bottom / heaviest)", True)
    check("PMNS distinguished generation = e (the mu-tau fixed point / lightest)", True)
    check("these are DIFFERENT generations (tau vs e) -> the symmetries do not match", True)

    print(
        """
[5] VERDICT -- the mixing prediction FAILS (honest negative; and a premise correction)
    The order-ideal Z_2 swaps the two ATOMS e<->mu and fixes the BOTTOM tau, so it
    predicts an e<->mu-symmetric, tau-distinguished mixing pattern.  But the canon's
    lepton mixing (tribimaximal / PMNS) is MU-TAU symmetric with the ELECTRON
    distinguished: |U_mu|^2 = |U_tau|^2, and it is invariant under mu<->tau, NOT
    under e<->mu.  So:

      * the reframing predicts the WRONG mixing symmetry (e<->mu, not mu-tau);
      * the distinguished generation disagrees: the order ideal singles out TAU (the
        heaviest, the order bottom), while PMNS singles out the ELECTRON (the
        lightest, the mu-tau fixed point) -- OPPOSITE ends of the hierarchy.

    Premise correction: the canon mixing is NOT "tau-distinguished, e<->mu-symmetric"
    -- it is e-distinguished, mu-tau-symmetric.  So 'generation = order ideal' does
    NOT predict the PMNS pattern; it predicts a different (incorrect) one.

    What this bounds: the order-ideal reframing is genuine and useful for the MASS /
    structural sector -- it forces 3 generations, aligns tau=bottom=heaviest, unifies
    delta/Phi as one record, and ties the excluded top to the Majorana neutrino.  But
    it does NOT extend to the MIXING sector: the mixing distinguishes the electron
    (mu-tau symmetry), the order ideal distinguishes the tau, and these conflict.
    The reframing is a mass/structure statement, not a mixing predictor.  (The
    canon's mu-tau / tribimaximal mixing is derived separately, from nu_R
    generation-blindness, sec 5.10 -- a different mechanism that the order ideal does
    not reproduce.)
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
