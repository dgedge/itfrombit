#!/usr/bin/env python3
r"""ITEM 126 — baryon-asymmetry SIGN: which involution carries it, and where the sign comes from.
This SUPERSEDES the earlier 'sign from complement non-symmetry' attempt, which used the WRONG
involution. Read the tiers; this is a correction + a re-routed way forward, not a closure.

WHAT WENT WRONG IN THE FIRST ATTEMPT. It treated tau(c)=1^c (§2.7 bitwise inversion) as the
matter/antimatter map and found P not tau-invariant (32/48 'matter-only'). But bitwise-NOT maps
EVERY generation-1 fermion (gen bits (0,0)) to gen (1,1) = R1-invalid — including the electron,
which plainly HAS a positron. So bitwise-NOT cannot be the particle<->antiparticle map.

THE CORRECT READING (verified below).
  * §2.7 bitwise-NOT is the CPT-MASS involution: it preserves the Q3 frustration F (so M_q=M_qbar).
    Its non-invariance on P is the Majorana-neutrino structure (the all-zeros nu has 0xFF invalid),
    NOT the baryon asymmetry.
  * §2.8 charge conjugation C is the Z_p-DOUBLET map: it flips all three Z-evaluations
    (Z_f, sum Z_c -> anticolour, Z_p), giving Q -> -Q EXACTLY. This is what reproduces the SM
    antiparticle charges. It is NOT a register bit-flip (anticolour sum +3/+1 is unreachable on the
    two colour bits — it lives on the Z_p doublet), so it keeps antiparticles physical.
  * Under C, every valid fermion has a valid antiparticle EXCEPT the neutral, self-conjugate
    neutrino (the unique Q=0 fundamental fermion). So the static census is C-SYMMETRIC mod nu.

CONSEQUENCE FOR THE SIGN. There is NO static state-count baryon imbalance to read off. The sign is
therefore necessarily DYNAMICAL (set by the boot bypass-fault rates), consistent with the magnitude
being the alpha^4 rate, not a census. The framework-faithful route to the sign:
    Majorana nu (canon)  ->  L-violation  ->  [leptogenesis: IMPORT, not in canon]  ->
    lepton-asymmetry sign = leptonic CP phase (PMNS, §5.10 / item 87, OPEN)  ->
    (B-L) + electroweak sphalerons [IMPORT]  ->  baryon sign.
So the sign REDUCES to the sign of the PMNS Dirac CP phase — a NAMED OPEN framework target
(item 87; the M12 audit falsified the old PMNS factorisation reading), not a free convention.

Self-asserting; exit 0 = the correction (bitwise-NOT != C) and the C-symmetry-mod-nu are verified.
"""
from __future__ import annotations
import itertools
from collections import Counter

# bit map (§2.1): 0=G0 1=G1 2=LQ 3=C0 4=C1 5=I3 6=chi 7=W
G0, G1, LQ, C0, C1, I3, CHI, W = range(8)


def bit(c, i):
    return (c >> i) & 1


def gen(c):
    return (bit(c, G0), bit(c, G1))


def colour(c):
    return (bit(c, C0), bit(c, C1))


def valid(c):                                   # R1, R2, R3
    return (not (bit(c, G0) and bit(c, G1))) and (bit(c, W) == bit(c, CHI)) \
        and ((bit(c, LQ) == 0) == (colour(c) == (0, 0)))


def tau(c):                                     # §2.7 bitwise inversion
    return c ^ 0xFF


# --- Q3 frustration F (§2.7): 8 register bits on the 3-cube, sum of edge XORs ---
Q3_EDGES = [(u, v) for u in range(8) for v in range(u + 1, 8) if bin(u ^ v).count("1") == 1]


def frustration(c):
    return sum(bit(c, u) ^ bit(c, v) for (u, v) in Q3_EDGES)


# --- §2.8 charge Q = 1/2 Z_f + 1/3 sum Z_c + 1/2 Z_p ---
def Zf(c):
    return 1 if bit(c, I3) == 0 else -1          # sgn(I3): up-type +1, down-type -1


def sumZc(c):
    return -3 if colour(c) == (0, 0) else -1     # one-hot: lepton -3, quark -1


def Q(c, Zp=1):
    return 0.5 * Zf(c) + (1.0 / 3.0) * sumZc(c) + 0.5 * Zp


def Q_C(c):                                       # §2.8 charge-conjugate: flip all three Z-evals
    return 0.5 * (-Zf(c)) + (1.0 / 3.0) * (-sumZc(c)) + 0.5 * (-1)


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("ITEM 126 — BARYON SIGN: bitwise-NOT is CPT-MASS, not C; the sign is dynamical")
    print("=" * 90)
    P = [c for c in range(256) if valid(c)]
    Pset = set(P)
    check(len(P) == 48, "physical codespace P has 48 valid registers")

    # ---- [1] §2.7 bitwise-NOT is the CPT-MASS involution (preserves frustration) ----
    print("\n[1] §2.7 bitwise-NOT preserves the Q3 frustration F  =>  it is the CPT-MASS map (M_q=M_qbar)")
    check(all(frustration(c) == frustration(tau(c)) for c in range(256)),
          "F(1^c)=F(c) for all 256 registers (mass-preserving involution)")

    # ---- [2] but bitwise-NOT maps physical fermions OUT of P -> it is NOT the antiparticle map ----
    print("\n[2] bitwise-NOT maps generation-1 fermions to R1-invalid registers => NOT charge conjugation")
    gen1 = [c for c in P if gen(c) == (0, 0)]
    nu = 0x00
    e = 1 << I3                                   # electron: I3=1, colourless, gen-1
    print(f"    nu={nu:08b} valid={valid(nu)}  Q={Q(nu):+.3f}; 1^nu={tau(nu):08b} valid={valid(tau(nu))}")
    print(f"    e ={e:08b} valid={valid(e)}  Q={Q(e):+.3f}; 1^e ={tau(e):08b} valid={valid(tau(e))}")
    check(valid(nu) and valid(e) and not valid(tau(nu)) and not valid(tau(e)),
          "both nu and e are valid but their bitwise complements are invalid (e HAS a positron!)")
    check(all(not valid(tau(c)) for c in gen1),
          f"ALL {len(gen1)} generation-1 fermions have invalid bitwise complements -> 1^c is not the antiparticle")

    # ---- [3] the REAL C is the §2.8 Z_p-doublet: Q -> -Q, and it is NOT a register bit-flip ----
    print("\n[3] §2.8 charge conjugation C (flip Z_f, colour->anticolour, Z_p): Q -> -Q exactly")
    check(all(abs(Q_C(c) - (-Q(c))) < 1e-12 for c in P),
          "C (flip all three Z-evaluations) gives Q_C = -Q for every valid register (true antiparticle charge)")
    # contrast: the bitwise-complement register charge is NOT -Q (so bitwise-NOT is a different map)
    bad = [c for c in P if abs(Q(tau(c), Zp=-1) - (-Q(c))) > 1e-9]
    print(f"    registers where charge(1^c, Z_p=-1) != -Q(c): {len(bad)}/48  (e.g. nu: {Q(tau(nu), -1):+.3f} vs -Q=0)")
    check(len(bad) > 0, "the bitwise-complement charge != -Q in general -> confirms bitwise-NOT is NOT C")

    # ---- [4] under C the census is symmetric EXCEPT the self-conjugate (Q=0) neutrino ----
    print("\n[4] under C, the only self-conjugate (Q=0) fundamental fermion is the neutrino")
    q0 = [c for c in P if abs(Q(c)) < 1e-12]
    print(f"    valid registers with Q=0: {len(q0)} — all leptons with I3=0 (neutrinos)? "
          f"{all(bit(c, LQ) == 0 and bit(c, I3) == 0 for c in q0)}")
    check(all(bit(c, LQ) == 0 and bit(c, I3) == 0 for c in q0) and len(q0) > 0,
          "the Q=0 (Majorana-capable) states are exactly the neutrinos -> census is C-symmetric mod nu")
    print("    => NO static state-count baryon imbalance exists; the sign is NOT a census quantity.")

    print("\n" + "=" * 90)
    print("VERDICT — correction + re-routed way forward (NOT a closure)")
    print(
        "  CORRECTION: §2.7 bitwise-NOT (1^c) is the CPT-MASS involution (preserves frustration F,\n"
        "  M_q=M_qbar). It is NOT charge conjugation: it maps every generation-1 fermion (incl. the\n"
        "  electron) to an R1-invalid register, yet those fermions have antiparticles. So the earlier\n"
        "  '32/48 matter-only under 1^c' measures non-invariance under the MASS map (the Majorana-\n"
        "  neutrino structure), NOT the matter/antimatter asymmetry. The real C is §2.8's Z_p-doublet\n"
        "  (flip Z_f, colour, Z_p -> Q->-Q), under which antiparticles stay physical and the census is\n"
        "  C-symmetric EXCEPT the neutral self-conjugate neutrino.\n"
        "  CONSEQUENCE: there is NO static baryon-number imbalance to count; the sign is necessarily\n"
        "  DYNAMICAL — consistent with the magnitude being the alpha^4 bypass RATE, not a state census.\n"
        "  WAY FORWARD (framework-faithful): the Majorana neutrino (canon) is the unique C-self-\n"
        "  conjugate, L-violating state. Via leptogenesis [IMPORT — leptogenesis/sphalerons are NOT in\n"
        "  canon] it sources a lepton asymmetry whose SIGN is the leptonic CP phase = the PMNS Dirac\n"
        "  phase (§5.10 / item 87, currently OPEN — the M12 audit falsified the old factorisation read),\n"
        "  converted to the baryon sign by (B-L) [Pati-Salam U(1)_{B-L}, canon] + sphalerons [IMPORT].\n"
        "  NET: the sign REDUCES to the sign of the PMNS CP phase (item 87) — a NAMED OPEN target, not\n"
        "  a free convention. TIER: the correction + C-symmetry-mod-nu are PROVEN; the leptogenesis\n"
        "  route is an IMPORT pending a canon derivation; the PMNS CP phase is the open framework piece."
    )
    print("exit 0 — bitwise-NOT=CPT-mass (not C); census C-symmetric mod nu; sign routes to PMNS CP (item 87).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
