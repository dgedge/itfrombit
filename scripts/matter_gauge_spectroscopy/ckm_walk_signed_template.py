#!/usr/bin/env python3
r"""CKM WALK-CP ENGINE — a clean SIGNED template (supersedes the ckm_audit.py PART-B negative).

WHY THIS EXISTS. ckm_audit.py PART B reported the bare walk engine as "qualitatively wrong
(suppresses the Cabibbo pair)". That reading was an ARTIFACT: it built only the UP block and used a
generation labeling under which the Hamming-2 pair is called "Cabibbo". Rebuilt correctly — BOTH
quark I3 sectors, the physical generation labels Gen1=(00),Gen2=(01),Gen3=(10), and the physical
(near-identity) CKM ordering — the SAME walk operator reproduces a correct CKM and the canon's bare
Jarlskog, with a DEFINITE sign.

ENGINE. Canonical Part-4 walk amplitudes A_0=sqrt(1-d), A_k=sqrt(d/7) e^{s i k pi/4} (d=2/9, s=+/-1
the phase orientation); coin-projected transfer W = sum_k A_k * CNOT(ctrl=(2-k)%8, targ=(5-k)%8);
mass operator M2=(W^dag W)^2; project the 9 left-handed quark states (3 gen x 3 colour, Boltzmann
colour weight e^{-d*|c|}) to a 3x3 per I3 sector; CKM = U_up^dag U_down.

WHAT THE FIXED ENGINE GIVES (verified below):
  * a correct CKM: |V_us| (Cabibbo) = 0.237 is the LARGEST off-diagonal (~5% over PDG 0.225), the
    Gen-3 couplings are all small; max element deviation from PDG ~ 0.033.
  * |J| = 4.33e-5 — EXACTLY the canon's bare-UV value (ANCHOR line 1032); 1.41x the measured
    3.08e-5 (the bare->IR RG factor, item 88).
  * a DEFINITE Jarlskog SIGN, fixed by the walk-phase orientation: s=+1 gives J=+4.33e-5
    (matching the measured sign J>0); s=-1 gives J=-4.33e-5.

HONEST LIMITS ON "PREDICTING" THE SIGN (this is a signed TEMPLATE, not a parameter-free sign):
  (a) The sign is set by the GLOBAL walk-phase orientation e^{+/- i k pi/4} (the cyclic handedness
      of the 8 bridge directions = the substrate's spatial orientation). It is NOT chirality-locked
      (J keeps its sign under chi: 0<->1, verified). So the sign has the SAME epistemic status as the
      Standard Model's CP sign: fixed by a global orientation convention, not a free per-observable knob.
  (b) The physical near-identity basis is REQUIRED to read the sign. The bare engine in raw
      mass-eigenvalue order gives near-MAXIMAL 1-2 mixing (|V|~0.97 off-diagonal, unphysical) with the
      OPPOSITE J sign; the near-identity (physical) ordering needs the RG flow (item 88) to align
      mass-order with the generation pairing. So the bare engine encodes the sign but does not
      autonomously select the physical basis.
  (c) CP requires the squared mass operator: M=(W^dag W)^1 gives J=0 identically.

NET: the quark sector is now a clean signed template — definite, experiment-consistent CP
(sign + bare magnitude + hierarchy) — refuting the negative result; but the sign rests on the global
orientation convention + the RG-aligned physical basis, not a from-nothing derivation.

Self-asserting; exit 0 = every statement above reproduces. (256-dim; runs in ~1s, no deep box needed.)
"""
from __future__ import annotations
import itertools

import numpy as np

DLT = 2.0 / 9.0
GENS = [(0, 0), (0, 1), (1, 0)]          # Gen1,Gen2,Gen3 (R1 forbids (1,1))
COLS = [(0, 1), (1, 0), (1, 1)]          # the three quark colours
PDG = np.array([[0.97435, 0.22501, 0.003732],
                [0.22487, 0.97349, 0.04183],
                [0.00858, 0.04111, 0.999118]])


def walk(sign: int) -> np.ndarray:
    Ak = np.zeros(8, complex)
    Ak[0] = np.sqrt(1 - DLT)
    for k in range(1, 8):
        Ak[k] = np.sqrt(DLT / 7) * np.exp(sign * 1j * k * np.pi / 4)
    W = np.zeros((256, 256), complex)
    for k in range(8):
        ctrl, targ = (2 - k) % 8, (5 - k) % 8
        for i in range(256):
            if (i >> ctrl) & 1:
                W[i ^ (1 << targ), i] += Ak[k]
            else:
                W[i, i] += Ak[k]
    return W


def sector_eig(W, i3, chi=0, colour_uniform=False, power=2):
    G = W.conj().T @ W
    M = G if power == 1 else G @ G

    def idx(g0, g1, c0, c1):
        return g0 | (g1 << 1) | (1 << 2) | (c0 << 3) | (c1 << 4) | (i3 << 5) | (chi << 6) | (chi << 7)

    if colour_uniform:
        w = np.ones(3) / np.sqrt(3)
    else:
        w = np.array([np.exp(-DLT * (a + b)) for (a, b) in COLS]); w /= np.linalg.norm(w)
    P = np.zeros((9, 3), complex)
    for g in range(3):
        P[g * 3:g * 3 + 3, g] = w
    ix = [idx(a, b, c, d) for (a, b) in GENS for (c, d) in COLS]
    return np.linalg.eigh(P.conj().T @ M[np.ix_(ix, ix)] @ P)


def jarlskog(V):
    return float(np.imag(V[0, 1] * V[1, 2] * np.conj(V[0, 2]) * np.conj(V[1, 1])))


def ckm(W, chi=0, colour_uniform=False, power=2, physical=True):
    _, Uu = sector_eig(W, 0, chi, colour_uniform, power)
    _, Ud = sector_eig(W, 1, chi, colour_uniform, power)
    if not physical:
        return Uu.conj().T @ Ud                     # raw mass-eigenvalue order
    best = None                                      # physical: reorder to near-identity (max diagonal)
    for pu in itertools.permutations(range(3)):
        for pd in itertools.permutations(range(3)):
            V = (Uu[:, pu].conj().T) @ Ud[:, pd]
            d = float(np.sum(np.abs(np.diag(V))))
            if best is None or d > best[0]:
                best = (d, V)
    return best[1]


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("CKM WALK-CP ENGINE — clean signed template")
    print("=" * 84)

    Wp = walk(+1)
    V = ckm(Wp); A = np.abs(V); J = jarlskog(V)
    np.set_printoptions(precision=4, suppress=True)
    print("\n[1] physical CKM (e^{+ikpi/4}, near-identity ordering):")
    print(A)
    print(f"    |V_us| (Cabibbo) = {A[0,1]:.4f} (PDG 0.225);  max|A-PDG| = {np.max(np.abs(A-PDG)):.3f}")
    print(f"    J = {J:+.4e}   |J| (PDG 3.08e-5; ratio {abs(J)/3.08e-5:.2f}, bare->IR)")
    gen3 = [A[0, 2], A[1, 2], A[2, 0], A[2, 1]]
    check(abs(A[0, 1] - 0.237) < 0.005, "|V_us| (Cabibbo) = 0.237, the correct Cabibbo magnitude")
    check(A[0, 1] > 5 * max(gen3), "Cabibbo |V_us| dominates all Gen-3 couplings (correct hierarchy)")
    check(max(gen3) < 0.03, "Gen-3 couplings are all small (<0.03)")
    check(abs(abs(J) - 4.33e-5) < 0.05e-5, "|J| = 4.33e-5 reproduces the canon bare value (ANCHOR line 1032)")
    check(J > 0, "Jarlskog sign is POSITIVE — matches the measured sign (J>0) for the e^{+} orientation")

    print("\n[2] the sign is the walk-phase orientation (s=+/-1), reproduced cleanly:")
    Jm = jarlskog(ckm(walk(-1)))
    print(f"    J(e^+) = {J:+.3e}   J(e^-) = {Jm:+.3e}")
    check(J > 0 and Jm < 0 and abs(J + Jm) < 1e-12, "orientation flips the sign; |J| identical")

    print("\n[3] NOT chirality-locked (J keeps its sign under chi 0<->1):")
    Jr = jarlskog(ckm(Wp, chi=1))
    print(f"    J(chi=0,L) = {J:+.3e}   J(chi=1,R) = {Jr:+.3e}")
    check(np.sign(J) == np.sign(Jr), "same sign for both handedness -> sign is the global orientation, not chirality")

    print("\n[4] the bare mass-order basis is unphysical (near-maximal 1-2 mixing, opposite sign):")
    Vraw = ckm(Wp, physical=False); Araw = np.abs(Vraw); Jraw = jarlskog(Vraw)
    print(f"    raw |V[0,1]| = {Araw[0,1]:.3f} (should be ~0.225);  J_raw = {Jraw:+.3e}")
    check(Araw[0, 1] > 0.9 and np.sign(Jraw) != np.sign(J),
          "raw mass-order gives near-maximal mixing + OPPOSITE sign -> physical basis (RG, item 88) needed")

    print("\n[5] CP requires the squared mass operator (M^1 gives J=0):")
    J1 = jarlskog(ckm(Wp, power=1))
    print(f"    J[(W^dag W)^1] = {J1:+.3e}")
    check(abs(J1) < 1e-12, "M=(W^dag W)^1 gives J=0; CP lives in the squared operator (W^dag W)^2")

    print("\n" + "=" * 84)
    print("VERDICT")
    print(
        "  The walk-CP engine, rebuilt correctly (both I3 sectors, physical labels + near-identity\n"
        "  ordering), is a CLEAN SIGNED TEMPLATE: it reproduces a correct CKM (Cabibbo 0.237 largest,\n"
        "  Gen-3 small, max dev 0.033 from PDG) and the canon's bare |J|=4.33e-5, with a DEFINITE\n"
        "  Jarlskog sign that is POSITIVE — matching the measured J>0 — for the e^{+ikpi/4} phase\n"
        "  orientation. This REFUTES the ckm_audit.py PART-B 'qualitatively wrong' negative (a\n"
        "  labeling + up-sector-only artifact). HONEST LIMITS: the sign is set by the GLOBAL walk\n"
        "  orientation (substrate spatial handedness; NOT chirality-locked) — the same status as the\n"
        "  SM's CP-sign convention, not a free per-observable knob; it is read in the physical\n"
        "  near-identity basis, which the bare mass-order does not autonomously select (RG/item 88\n"
        "  needed); and CP lives in (W^dag W)^2. NET: a definite, experiment-consistent signed CP\n"
        "  template for the quark sector — the leptogenesis-sign argument's template is now clean.\n"
        "  The lepton/nu_R extension (item 87) inherits the SAME global orientation, so the\n"
        "  leptogenesis sign and the CKM sign are CORRELATED, not independent — but item 87 (the\n"
        "  PMNS/nu_R bridge, M12-falsified) remains the open piece."
    )
    print("exit 0 — quark walk-CP is a clean signed template (sign = +J for e^{+}); item-87 lepton bridge open.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
