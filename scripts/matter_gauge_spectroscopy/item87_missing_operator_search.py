#!/usr/bin/env python3
r"""ITEM 87 — search for the missing operator that would carry CP into the lepton sector.

*** PARTIALLY SUPERSEDED (2026-06-16, see `item87_majorana_cp_operator.py`): the conclusion below
that the leptonic CP operator gives "no CP / is external" is CORRECTED. The leptogenesis step here
("Im[(Y^dag Y)^2]=0 -> CP vanishes") used the WRONG (original) weak basis -- it is trivially 0 for a
diagonal Y -- and is not the physical CP measure. Basis-invariant weak-basis invariants show a
complex-symmetric Majorana M_R DOES give physical leptonic CP. STILL VALID below: Route A's
structural no-go, and the root cause that the walk alone gives a real-diagonal Y_nu (so the SUBSTRATE
alone has no leptonic CP). Only the "missing operator gives no CP" framing is wrong: the operator
(a complex M_R) exists and works; its substrate DERIVATION is the open piece. ***

item87_lepton_cp_obstruction.py showed leptonic CP = 0 because the walk's generation-flips are
gated by I3 (G0-flip) and chi (G1-flip). This script tests the two named candidate operators plus
the natural seesaw, to either FIND the missing operator or characterise it precisely.

CANDIDATES TESTED (all self-asserting below):
  A.  Virtual colour bridge (Part 9): a lepton borrows colour paths via virtual coloured states.
  B.  Complement-sourced Majorana M_R: the §2.7 complement couples nu generations (Gen2<->Gen3).
  B'. Type-I seesaw with M_R = the right-handed nu_R walk block (which couples Gen1-Gen2).

RESULTS.
  A  — STRUCTURALLY DEAD. The walk cannot connect ANY two valid neutrino generations at ANY order
       (|W^n[nu_i,nu_j]| = 0 for all n<=12, i!=j). Gating chain: the G0-flip needs I3 (which needs
       LQ=1, absent for leptons); the G1-flip needs chi (which needs C0, which needs G0=1, but the
       G1-flip from Gen3=(1,0) lands on the R1-FORBIDDEN (1,1)). A clean no-go, stronger than the
       M12 angle-falsification.
  B, B' — BOTH give nonzero PMNS ANGLES in the right ballpark (B ~55,36,31; B' ~58,54,27; observed
       33,45,8.5) but BOTH give J = 0. So the type-I seesaw IS a plausible route to the PMNS-angle
       reconciliation (item 87 Q1) -- it produces all three nonzero angles, improving on the
       M12-falsified U_nu^dag U_L cross-talk -- BUT the precise angle values are diagonalisation-
       sensitive (SVD-vs-eigh, ordering), so they are a ballpark direction, not a sharp prediction.
       The ROBUST result across every route is J = 0 (no leptonic CP).

ROOT CAUSE OF THE PERSISTENT J=0 (the real wall). The quark CP comes from the dominant CNOT k=0 =
CNOT(control=LQ, target=I3), which fires for quarks (LQ=1) and makes the down sector's amplitudes
complex relative to the up sector. For LEPTONS (LQ=0) this CNOT is the identity, so the neutrino
Dirac coupling Y_nu is REAL and DIAGONAL -> the PMNS Jarlskog = 0 AND the leptogenesis CP invariant
Im[(Y_nu^dag Y_nu)^2] = 0. The Majorana M_R (Delta L=2, complex-SYMMETRIC) is a different operator
class entirely -- not derivable from the Hermitian / Delta L=0 walk -- so it is not the walk
substrate's own object.

CONCLUSION. The missing CP operator is NOT among the named candidates and NOT in the documented walk
substrate. It must be EITHER (a) a lepton-sector phase source NOT gated by LQ (to make Y_nu complex),
or (b) a complex-symmetric Majorana M_R from new high-scale/boot dynamics (Delta L=2 seesaw input).
The SIGN, if such an operator existed, would still inherit the global walk orientation (the
right-charged CP triangle flips with it), so it would correlate the leptogenesis sign with the CKM.
Partial win: the seesaw (B/B') is the right structural route for the PMNS ANGLES (item 87 Q1); the
CP/sign is the wall.

Self-asserting; exit 0 = the no-go, the angle-ballpark-with-J=0 of B and B', and the real-Y_nu root
cause all verify.
"""
from __future__ import annotations
import numpy as np

GENS = [(0, 0), (0, 1), (1, 0)]


def walk(sign=1):
    dlt = 2 / 9
    Ak = np.zeros(8, complex); Ak[0] = np.sqrt(1 - dlt)
    for k in range(1, 8):
        Ak[k] = np.sqrt(dlt / 7) * np.exp(sign * 1j * k * np.pi / 4)
    W = np.zeros((256, 256), complex)
    for k in range(8):
        c, t = (2 - k) % 8, (5 - k) % 8
        for i in range(256):
            if (i >> c) & 1:
                W[i ^ (1 << t), i] += Ak[k]
            else:
                W[i, i] += Ak[k]
    return W


def blk(M2, chi, i3):
    def idx(a, b):
        return a | (b << 1) | (i3 << 5) | (chi << 6) | (chi << 7)
    ix = [idx(a, b) for a, b in GENS]
    return M2[np.ix_(ix, ix)]


def seesaw_pmns(M2, UL, y, MR):
    mnu = np.diag(y) @ np.linalg.inv(MR) @ np.diag(y)
    _, Un = np.linalg.eigh((mnu + mnu.conj().T) / 2)
    P = Un.conj().T @ UL
    A = np.abs(P)
    ang = (np.degrees(np.arctan2(A[0, 1], A[0, 0])),
           np.degrees(np.arctan2(A[1, 2], A[2, 2])),
           np.degrees(np.arcsin(min(1.0, A[0, 2]))))
    J = float(np.imag(P[0, 1] * P[1, 2] * np.conj(P[0, 2]) * np.conj(P[1, 1])))
    return ang, J


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("ITEM 87 — missing-operator search (carry CP into the lepton sector)")
    print("=" * 88)
    W = walk(+1)
    M2 = (W.conj().T @ W) @ (W.conj().T @ W)

    # ---- A: virtual colour bridge — can the walk connect neutrino generations at any order? ----
    print("\n[A] virtual colour bridge: does any W^n connect different-generation neutrinos?")
    nuc = [0, 2, 1]                                # nu codewords: Gen1=(00)->0, Gen2=(01)->2, Gen3=(10)->1
    Wn = np.eye(256, dtype=complex); maxoff = 0.0
    for _ in range(12):
        Wn = Wn @ W
        for i in range(3):
            for j in range(3):
                if i != j:
                    maxoff = max(maxoff, abs(Wn[nuc[i], nuc[j]]))
    print(f"    max |W^n[nu_i,nu_j]| (n<=12, i!=j) = {maxoff:.2e}")
    check(maxoff < 1e-12, "A: the walk CANNOT connect neutrino generations at any order (structural no-go)")

    _, UL = np.linalg.eigh(blk(M2, 0, 1))          # charged-lepton mixing (walk: 1-3 rotation)
    Yd = blk(M2, 0, 0)                             # light neutrino Dirac block (walk)
    y = np.sqrt(np.abs(np.diag(Yd)))

    # ---- B: complement-sourced Majorana M_R (Gen2<->Gen3) ----
    print("\n[B] complement Majorana M_R (Gen2<->Gen3):")
    ang_b, J_b = seesaw_pmns(M2, UL, y, np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], complex))
    print(f"    PMNS (th12,th23,th13)=({ang_b[0]:.1f},{ang_b[1]:.1f},{ang_b[2]:.1f})  J={J_b:+.2e}")
    check(all(a > 5.0 for a in ang_b) and abs(J_b) < 1e-6,
          "B: nonzero ballpark PMNS angles but J = 0 (mixing yes, CP no)")

    # ---- B': type-I seesaw with M_R = the nu_R (chi=1) walk block ----
    print("\n[B'] type-I seesaw, M_R = nu_R (chi=1) walk block (couples Gen1-Gen2):")
    ang_bp, J_bp = seesaw_pmns(M2, UL, y, blk(M2, 1, 0))
    print(f"    PMNS (th12,th23,th13)=({ang_bp[0]:.1f},{ang_bp[1]:.1f},{ang_bp[2]:.1f})  J={J_bp:+.2e}"
          f"   (observed ~ 33,45,8.5)")
    check(all(a > 5.0 for a in ang_bp) and abs(J_bp) < 1e-6,
          "B': all three PMNS angles nonzero (ballpark) but J = 0 -- seesaw gives angles, not CP")

    # ---- root cause: the CP-source CNOT is LQ-gated; Y_nu is real -> all leptonic CP = 0 ----
    print("\n[root] the CP-source CNOT k=0 is CNOT(LQ -> I3): fires for quarks (LQ=1), identity for leptons")
    c0, t0 = (2 - 0) % 8, (5 - 0) % 8
    check(c0 == 2 and t0 == 5, "the dominant CNOT k=0 has control=LQ(bit2), target=I3(bit5) -> LQ-gated")
    check(np.max(np.abs(Yd.imag)) < 1e-9 and np.max(np.abs(Yd - np.diag(np.diag(Yd)))) < 1e-9,
          "the light neutrino Dirac block Y_nu is REAL and DIAGONAL (no LQ-gated phase reaches leptons)")
    YdY = Yd.conj().T @ Yd
    lepto_eps = float(np.max(np.abs(np.imag(YdY @ YdY))))
    print(f"    leptogenesis CP invariant Im[(Y_nu^dag Y_nu)^2] max = {lepto_eps:.2e}")
    check(lepto_eps < 1e-9, "Im[(Y_nu^dag Y_nu)^2] = 0 -> leptogenesis CP vanishes (real Y_nu)")

    print("\n" + "=" * 88)
    print("VERDICT — the missing CP operator is NOT among the candidates and NOT in the walk substrate")
    print(
        "  A (virtual colour bridge): STRUCTURALLY DEAD -- the walk connects no two neutrino\n"
        "    generations at any order (G0-flip needs LQ=1; the G1-bridge from Gen3 hits R1-forbidden).\n"
        "  B / B' (seesaw, complement-M_R and nu_R-block-M_R): BOTH give all three PMNS angles nonzero\n"
        "    in the right ballpark (~55,36,31 and ~58,54,27 vs 33,45,8.5) -- so the type-I seesaw is a\n"
        "    plausible route to the item-87 Q1 ANGLE reconciliation, improving on the M12-falsified\n"
        "    cross-talk -- but the angle values are diagonalisation-sensitive (a direction, not a sharp\n"
        "    prediction), and BOTH give J = 0.\n"
        "  ROOT (the wall): the quark CP comes from the LQ-gated CNOT(LQ->I3); for leptons (LQ=0) it is\n"
        "    identity, so Y_nu is REAL-diagonal and BOTH the PMNS Jarlskog and the leptogenesis invariant\n"
        "    Im[(Y^dag Y)^2] vanish. The Majorana M_R (Delta L=2, complex-symmetric) is a different\n"
        "    operator class, not derivable from the Hermitian / Delta L=0 walk.\n"
        "  MISSING CP OPERATOR (precisely characterised): EITHER (a) a lepton-sector phase source NOT\n"
        "    gated by LQ (complex Y_nu), OR (b) a complex-symmetric Majorana M_R from new high-scale /\n"
        "    boot dynamics. Both are OUTSIDE the documented walk substrate. IF either existed, the sign\n"
        "    would inherit the global walk orientation (right-charged CP flips with it), correlating the\n"
        "    leptogenesis sign with the CKM. TIER: candidate search -- virtual-colour-bridge falsified\n"
        "    (no-go); the seesaw located as the right ANGLE route; the CP/sign operator located OUTSIDE\n"
        "    the substrate. Item 87 stays open at that external CP operator."
    )
    print("exit 0 — A dead; seesaw gives angles not CP; leptonic CP operator is external (item 87 open).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
