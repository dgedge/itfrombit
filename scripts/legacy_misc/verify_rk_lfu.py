#!/usr/bin/env python3
"""
R_K / R_K* lepton-flavour-universality structural check for TCH.

Question:
  Does the documented TCH neutral-current lepton coupling carry any
  generation (e/mu/tau) dependence -- i.e. can the framework produce a primary
  muon-vs-electron LFU violation in b -> s l l?

Result (asserted):
  No. The electric charge Q and weak isospin T3 are functions of the bits
  {I3, CHI, C0, C1} ONLY; the generation bits {G0, G1} never enter. Hence the
  gamma/Z coupling to a charged lepton is identical for e, mu, tau ->
  C9, C10 are lepton-universal -> R_K = R_K* = 1 at the WET-coefficient level
  (= SM, up to m_l-suppressed scalar/QED effects).

Consequence (the falsifiable structural claim):
  Any TCH contribution to C9/C10 is necessarily LEPTON-UNIVERSAL. A confirmed
  muon-specific LFU violation would require {G0,G1} on the lepton current, which
  the documented charge structure does not have -> would falsify TCH at leading
  order. The signed, q^2-shaped universal-C9 magnitude is NOT derivable here --
  it needs an explicit WET-matching of the Part 04 walk/GIM amplitude (unbuilt).

Documented charge formulas: from check_gauge_dressing.py / ANCHOR §2 charge assignment.
bit order: G0,G1,LQ,C0,C1,I3,CHI,W.
"""
from __future__ import annotations

G0, G1, LQ, C0, C1, I3, CHI, W = range(8)


def Zf(c):    return 1 if c[I3] == 0 else -1
def sumZc(c): return -3 if (c[C0], c[C1]) == (0, 0) else -1
def Q(c):     return 0.5 * Zf(c) + (1 / 3) * sumZc(c) + 0.5
def T3(c):    return (0.5 if c[I3] == 0 else -0.5) if c[CHI] == 0 else 0.0


def state(lq, col, i3, chi, gen=(0, 0)):
    c = [0] * 8
    c[G0], c[G1] = gen
    c[LQ] = lq
    c[C0], c[C1] = col
    c[I3] = i3
    c[CHI] = chi
    return c


def main() -> None:
    # sanity: documented charges reproduce the SM
    e  = Q(state(0, (0, 0), 1, 0))
    nu = Q(state(0, (0, 0), 0, 0))
    u  = Q(state(1, (0, 1), 0, 0))
    d  = Q(state(1, (0, 1), 1, 0))
    print("SM-charge sanity: e=%.3f nu=%.3f u=%.3f d=%.3f" % (e, nu, u, d))
    assert abs(e + 1) < 1e-12 and abs(nu) < 1e-12
    assert abs(u - 2 / 3) < 1e-12 and abs(d + 1 / 3) < 1e-12

    # the test: Q,T3 across the three charged-lepton generations, both chiralities
    gens = [(0, 0), (0, 1), (1, 0)]
    for chi, tag in [(0, "L"), (1, "R")]:
        Qs  = [Q(state(0, (0, 0), 1, chi, gen=g)) for g in gens]
        T3s = [T3(state(0, (0, 0), 1, chi, gen=g)) for g in gens]
        print("charged leptons (%s): Q=%s  T3=%s" % (tag, Qs, T3s))
        assert len(set(Qs)) == 1, "Q is generation-dependent (unexpected)"
        assert len(set(T3s)) == 1, "T3 is generation-dependent (unexpected)"

    # explicit: Q,T3 do not read the generation bits at all
    base = state(0, (0, 0), 1, 0, gen=(0, 0))
    for g in gens:
        c = state(0, (0, 0), 1, 0, gen=g)
        assert Q(c) == Q(base) and T3(c) == T3(base)

    print("\nPASS: neutral-current lepton coupling is generation-blind (LFU).")
    print("  -> C9 = C9^e = C9^mu = C9^tau, C10 likewise -> R_K = R_K* = 1 (= SM).")
    print("  -> TCH cannot produce a primary muon-specific LFU violation; any TCH")
    print("     C9/C10 effect is necessarily lepton-universal. A confirmed muon-")
    print("     specific violation would falsify TCH at leading order.")
    print("  -> the signed universal-C9 magnitude is the open Part 04b WET-matching target.")


if __name__ == "__main__":
    main()
