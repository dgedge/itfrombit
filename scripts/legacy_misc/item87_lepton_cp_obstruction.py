#!/usr/bin/env python3
r"""ITEM 87 — why leptonic CP (PMNS / nu_R) vanishes in the walk substrate: a precise root cause.

CONTEXT. Item 87 (the Part-5 <-> Part-9 PMNS reconciliation / UV->IR bridge) is M12-falsified:
U_nu^dag U_L does not give the observed angles, and the documented Part-18 Feshbach is generation-
block-diagonal, so the neutrino inter-generation coupling B_nu = 0 (no phase to lift). This is the
piece that would carry the global walk orientation (the e^{+ik pi/4} that fixes the quark CKM sign,
`ckm_walk_signed_template.py`) into the lepton sector for the leptogenesis sign. This script
root-causes the obstruction at the walk level.

RESULT (verified below). The two generation-flip CNOTs of the canonical walk are controlled by I3
(the G0-flip, k=5) and chi (the G1-flip, k=4). Generation mixing therefore requires I3=1 and/or
chi=1, so the four lepton sectors (chi, I3) mix DIFFERENT generation pairs, and in every NEUTRINO
sector at least one generation decouples:
    light nu      (chi=0, I3=0): neither flip  -> FULLY DIAGONAL (no mixing, no phase)
    light charged (chi=0, I3=1): G0-flip only  -> Gen2 decoupled (effective 2-gen)
    nu_R          (chi=1, I3=0): G1-flip only  -> Gen3 decoupled (effective 2-gen)
    charged_R     (chi=1, I3=1): BOTH flips    -> full 3-gen, nonzero CP triangle
The rephasing-invariant 3-generation CP measure Im(M01 M12 M02*) is ZERO for the light neutrinos,
the light charged leptons, AND the nu_R -- nonzero ONLY for right-handed charged leptons, which are
neither the PMNS source (light, left-handed) nor the leptogenesis source (nu_R). Hence the LIGHT
PMNS Jarlskog = 0 exactly (robust across operator powers), and the nu_R sector carries no
rephasing-invariant CP either (a 2-generation Dirac phase is removable by rephasing).

CONSEQUENCE FOR THE SIGN. The obstruction is about the EXISTENCE of leptonic CP, NOT its sign:
there is no leptonic/nu_R CP invariant in the walk substrate to carry the orientation. The missing
ingredient is precisely a generation-flip for a NEUTRINO (I3=0) sector NOT controlled by I3/chi --
exactly the framework's own "lepton mixing requires the strong force" virtual colour bridge
(Part 9 / ANCHOR §5.10.1), which is M12-falsified on the angles and supplies no CP at the M^2
level. IF such an operator existed, its phase would be the SAME global e^{+ik pi/4} orientation
that fixes the quark CKM sign, so the leptogenesis sign WOULD be correlated with the (now-clean)
CKM sign -- but that is conditional on closing the existence gap. Item 87 stays OPEN; this
sharpens it from "UV->IR bridge open" to a verified structural mechanism.

Self-asserting; exit 0 = the control structure, the per-sector decoupling, J_light=0, and the
nu_R CP-triangle=0 all verify.
"""
from __future__ import annotations
import numpy as np

GENS = [(0, 0), (0, 1), (1, 0)]                 # Gen1, Gen2, Gen3 (R1 forbids (1,1))
NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]


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


def lepton_block(M, chi, i3):
    def idx(a, b):
        return a | (b << 1) | (0 << 2) | (0 << 3) | (0 << 4) | (i3 << 5) | (chi << 6) | (chi << 7)
    ix = [idx(a, b) for (a, b) in GENS]
    return M[np.ix_(ix, ix)]


def triangle(M):                                 # rephasing-invariant 3-generation CP measure
    return float(np.imag(M[0, 1] * M[1, 2] * np.conj(M[0, 2])))


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("ITEM 87 — leptonic-CP obstruction root-caused (I3/chi-controlled generation flips)")
    print("=" * 90)

    # ---- [1] the two generation-flip CNOTs and their controls ----
    print("\n[1] generation-flip CNOTs (target=(5-k)%8, control=(2-k)%8):")
    flips = {}
    for k in range(8):
        c, t = (2 - k) % 8, (5 - k) % 8
        if t in (0, 1):
            flips[NAMES[t]] = NAMES[c]
            print(f"    k={k}: flips {NAMES[t]} controlled by {NAMES[c]}")
    check(flips.get("G0") == "I3" and flips.get("G1") == "chi",
          "G0-flip is controlled by I3, G1-flip by chi -> mixing needs I3=1 and/or chi=1")

    W = walk(+1)
    M2 = (W.conj().T @ W) @ (W.conj().T @ W)

    # ---- [2] per-sector generation-mixing + CP triangle ----
    print("\n[2] per-sector generation mixing (|off-diagonals|) and the CP triangle:")
    print(f"    {'sector':22s} {'g1g2':>6s} {'g2g3':>6s} {'g1g3':>6s} {'CPtri':>10s}  full3gen")
    sectors = [(0, 0, "light nu  (L)"), (0, 1, "light charged (L)"),
               (1, 0, "nu_R      (R)"), (1, 1, "charged_R (R)")]
    tris = {}
    for chi, i3, name in sectors:
        M = lepton_block(M2, chi, i3)
        a, b, c = abs(M[0, 1]), abs(M[1, 2]), abs(M[0, 2])
        tri = triangle(M); tris[name] = (a, b, c, tri)
        print(f"    {name:22s} {a:6.3f} {b:6.3f} {c:6.3f} {tri:+10.2e}  {a>1e-9 and b>1e-9 and c>1e-9}")

    check(all(abs(lepton_block(M2, 0, 0)[i, j]) < 1e-9 for i in range(3) for j in range(3) if i != j),
          "light neutrino block is FULLY DIAGONAL (no mixing, no phase)")
    a, b, c, _ = tris["light charged (L)"]
    check(a < 1e-9 and b < 1e-9 and c > 1e-9, "light charged leptons: only Gen1-Gen3 (Gen2 decoupled)")
    a, b, c, _ = tris["nu_R      (R)"]
    check(a > 1e-9 and b < 1e-9 and c < 1e-9, "nu_R: only Gen1-Gen2 (Gen3 decoupled)")
    check(abs(tris["light nu  (L)"][3]) < 1e-9 and abs(tris["light charged (L)"][3]) < 1e-9
          and abs(tris["nu_R      (R)"][3]) < 1e-9,
          "CP triangle = 0 for light nu, light charged, AND nu_R (each decouples a generation)")
    check(abs(tris["charged_R (R)"][3]) > 1e-6,
          "CP triangle != 0 ONLY for right-handed charged leptons (the lone full-3-gen sector)")

    # ---- [3] the LIGHT PMNS Jarlskog is 0, robust across operator powers ----
    print("\n[3] light (physical) PMNS Jarlskog across operator powers (must be 0):")
    G = W.conj().T @ W
    for p in (1, 2, 3, 4):
        Mp = np.linalg.matrix_power(G, p)
        _, Uc = np.linalg.eigh(lepton_block(Mp, 0, 1))     # light charged
        _, Un = np.linalg.eigh(lepton_block(Mp, 0, 0))     # light neutrino
        P = Uc.conj().T @ Un
        J = float(np.imag(P[0, 1] * P[1, 2] * np.conj(P[0, 2]) * np.conj(P[1, 1])))
        print(f"    power {p}: light PMNS Jarlskog = {J:+.2e}")
        check(abs(J) < 1e-8, f"light PMNS Jarlskog = 0 at power {p} (no light leptonic CP)")

    # ---- [4] the sign question is moot: there is no CP to carry; if it existed it'd be global ----
    print("\n[4] orientation check — the right-charged CP triangle IS the global walk orientation:")
    tri_p = triangle(lepton_block((W.conj().T @ W) @ (W.conj().T @ W), 1, 1))
    Wm = walk(-1); tri_m = triangle(lepton_block((Wm.conj().T @ Wm) @ (Wm.conj().T @ Wm), 1, 1))
    print(f"    charged_R CP triangle: e^+ -> {tri_p:+.2e},  e^- -> {tri_m:+.2e}")
    check(tri_p > 0 and tri_m < 0,
          "the one nonzero lepton CP (right-charged) flips with the walk orientation -> if a "
          "neutrino-sector CP existed it would inherit the SAME global orientation as the CKM")

    print("\n" + "=" * 90)
    print("VERDICT — item 87 sharpened (root cause), not closed")
    print(
        "  Leptonic CP vanishes in every physical lepton sector for a precise structural reason: the\n"
        "  walk's two generation-flip CNOTs are controlled by I3 (G0-flip) and chi (G1-flip), so a\n"
        "  generation decouples in every NEUTRINO sector -- light nu (both flips off -> fully diagonal),\n"
        "  light charged (Gen2 off), nu_R (Gen3 off). The rephasing-invariant CP triangle is therefore\n"
        "  0 for the light neutrinos, the light charged leptons, and the nu_R; it is nonzero ONLY for\n"
        "  right-handed charged leptons, which source neither the PMNS (light) nor leptogenesis (nu_R).\n"
        "  So the LIGHT PMNS Jarlskog = 0 exactly (robust across powers), confirming and root-causing\n"
        "  the item-87 obstruction (= the generation-block-diagonal Feshbach / B_nu=0, now seen to be\n"
        "  the I3/chi control structure, not just a Feshbach accident).\n"
        "  CONSEQUENCE: the global walk orientation CANNOT currently be carried into the lepton sector\n"
        "  -- not because the sign is wrong, but because there is NO leptonic/nu_R CP invariant to carry.\n"
        "  The missing operator is precisely a generation-flip for a NEUTRINO (I3=0) sector NOT gated by\n"
        "  I3/chi -- the framework's 'lepton mixing needs the strong force' virtual colour bridge,\n"
        "  M12-falsified on the angles and CP-free at M^2. The one nonzero lepton CP (right-charged)\n"
        "  DOES flip with the walk orientation, so IF the missing operator supplied a neutrino-sector\n"
        "  CP it would inherit the SAME global orientation as the (now-clean) quark CKM sign -- making\n"
        "  the leptogenesis sign correlated with it. TIER: structural diagnosis / sharpening of item 87;\n"
        "  the existence gap (a neutral-sector inter-generation operator, or a complex Majorana M_R) is\n"
        "  unchanged-open. The leptogenesis SIGN remains open via that existence gap, not via the sign."
    )
    print("exit 0 — leptonic CP=0 root-caused to I3/chi-gated generation flips; item 87 existence gap open.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
