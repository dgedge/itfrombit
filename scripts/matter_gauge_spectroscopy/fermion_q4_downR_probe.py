#!/usr/bin/env python3
r"""fermion_q4_downR_probe.py -- low-odds swing at the down-quark R anomaly (piece 2).

R = sqrt(#paths) holds for leptons (sqrt2, 2 axes) and up quarks (sqrt3, 3 colours), validated. But the
down quark -- also a colour triplet -- fits R_d ~ 1.543, NOT sqrt3, "no colour story" (canon quotes
sqrt(12/5)~1.549). Two candidate rescues:
  (a) R_d = sqrt3 (colour, restoring the law) + an NLO node-dressing (the down gen-1 sits at a node,
      like up's m_u; canon dresses up's R sqrt3->1.778 by +2.6% "amplified by node proximity").
  (b) R_d = sqrt(12/5) is a phenomenological fit.
This probe tests both. Koide circulant: m_n = (1 + R cos(delta + 2pi n/3))^2, delta_d=1/9.
PDG down: m_d:m_s:m_b = 4.67 : 93.4 : 4180 MeV.

RESULT (sharper than expected, both routes NEGATIVE): (A) bare sqrt3 + node-dressing is INCONSISTENT
(dressing the gen-1 node alone cannot fix both down ratios -- the required boosts disagree ~33%).
(B) the down gen-1 sits at a circulant NODE, so the ratios are HYPERSENSITIVE to R: the full
m_s/m_d AND m_b/m_d ratios pin R to ~1.543 +- 0.1%, and even sqrt(12/5) (0.4% off) gives ~8% ratio
error -- so sqrt(12/5) is a loose approximation, NOT the value, and NO simple alphabet form fits at the
demanded precision. The down R is a bare, node-fragile EMPIRICAL number with no clean (structural OR
alphabet) form. The node (m_d node-suppression) is the root anomaly. Piece 2 stays OPEN/negative.
"""
import numpy as np

PDG_DOWN = np.array([4.67, 93.4, 4180.0])   # m_d, m_s, m_b (MeV)
DELTA = 1.0 / 9.0


def factors(R, delta=DELTA):
    return np.array([1 + R * np.cos(delta + 2 * np.pi * n / 3) for n in range(3)])


def masses(R, delta=DELTA):
    return np.sort(factors(R, delta) ** 2)


def ratio_logrms(R):
    m = masses(R)
    if m[0] <= 0:
        return 9.9
    pred = m / m[0]
    tgt = np.sort(PDG_DOWN); tgt = tgt / tgt[0]
    return float(np.sqrt(np.mean((np.log(pred[1:]) - np.log(tgt[1:])) ** 2)))


def main():
    print("=== Down-quark R anomaly: sqrt3 + node-dressing, or a node-fragile empirical fit? ===\n")

    # ---- Test A: the node-dressing story (bare R=sqrt3) ----
    print("  [A] Node-dressing story: bare R=sqrt3 (colour), dress the gen-1 node to fix m_d.")
    R3 = np.sqrt(3.0)
    f = factors(R3)
    md, ms, mb = np.sort(f ** 2)
    print(f"      at R=sqrt3, delta=1/9:  f_n = {np.round(f,4)}  -> min|f| = {np.min(np.abs(f)):.4f} (the node)")
    print(f"      bare ratios  m_s/m_d = {ms/md:7.1f}   m_b/m_d = {mb/md:8.1f}")
    print(f"      PDG  ratios  m_s/m_d = {PDG_DOWN[1]/PDG_DOWN[0]:7.1f}   m_b/m_d = {PDG_DOWN[2]/PDG_DOWN[0]:8.1f}")
    boost_ms = np.sqrt((ms / md) / (PDG_DOWN[1] / PDG_DOWN[0]))
    boost_mb = np.sqrt((mb / md) / (PDG_DOWN[2] / PDG_DOWN[0]))
    print(f"      f_d boost to fix m_s/m_d: x{boost_ms:.2f};  to fix m_b/m_d: x{boost_mb:.2f}")
    consistent = abs(boost_ms - boost_mb) / boost_ms < 0.10
    print(f"      -> required boosts {'AGREE' if consistent else 'DISAGREE'} "
          f"({100*abs(boost_ms-boost_mb)/boost_ms:.0f}% apart): node-dressing of m_d alone "
          f"{'works' if consistent else 'CANNOT fix both ratios'}.  Rescue FAILS.\n")

    # ---- Test B: how tightly do the (node-sensitive) down ratios pin R, and does any clean form fit? ----
    print("  [B] Best-fit R vs simple-alphabet forms (delta=1/9 fixed); the gen-1 node makes ratios sensitive.")
    Rs = np.linspace(1.30, 1.80, 50001)
    errs = np.array([ratio_logrms(R) for R in Rs])
    Rstar = float(Rs[int(np.argmin(errs))]); e_star = float(errs.min())
    print(f"      best-fit R* = {Rstar:.4f}  (ratio log-RMS {e_star:.4f} ~ {100*e_star:.1f}% on ratios)")
    alphabet = {
        "sqrt2": np.sqrt(2), "sqrt3": np.sqrt(3), "phi": (1 + 5 ** 0.5) / 2,
        "sqrt(12/5)": np.sqrt(12 / 5), "sqrt(7/3)": np.sqrt(7 / 3), "sqrt(5/2)": np.sqrt(5 / 2),
        "17/11": 17 / 11, "14/9": 14 / 9, "31/20": 31 / 20, "sqrt(19/8)": np.sqrt(19 / 8),
    }
    print(f"      {'form':12s} {'R':>7s} {'dev% vs R*':>10s} {'ratio-RMS':>10s} {'%ratios':>8s} {'fits<1.5x?':>10s}")
    nfit = 0
    for name, R in sorted(alphabet.items(), key=lambda kv: abs(kv[1] - Rstar)):
        rms = ratio_logrms(R); fits = rms < 1.5 * e_star + 0.001; nfit += fits
        print(f"      {name:12s} {R:7.4f} {100*abs(R-Rstar)/Rstar:10.1f} {rms:10.4f} {100*rms:8.0f} {str(fits):>10s}")
    print(f"\n      simple forms fitting the FULL ratios comparably (<1.5x best RMS): {nfit}")
    print(f"      -> the gen-1 NODE makes the ratios hypersensitive: R is pinned to {Rstar:.4f} +-0.1%, and")
    print(f"         even sqrt(12/5) (0.4% off) gives {100*ratio_logrms(np.sqrt(12/5)):.0f}% ratio error (NOT the ~1% canon")
    print(f"         quotes for a single observable). NO clean form fits at the demanded precision; the down")
    print(f"         R is a bare, node-fragile EMPIRICAL number.")

    print("\n[verdict] PIECE 2 -- LOW-ODDS SWING returns NEGATIVE (and sharper than expected):")
    print("  (A) the sqrt3-colour + node-dressing rescue is INCONSISTENT -- dressing the gen-1 node cannot")
    print("      fix both down mass-ratios (the required boosts disagree ~33%); a global sqrt3->1.543 shift")
    print("      (-10.9%) is opposite-sign and ~4x larger than up's +2.6%, with no derived reason.")
    print("  (B) the down gen-1 (m_d) sits at a circulant NODE -> the ratios are hypersensitive to R,")
    print("      pinning it to ~1.543 +-0.1%; NO simple form fits at that precision (even sqrt(12/5), 0.4%")
    print("      off, gives ~8% ratio error). So sqrt(12/5) is a loose label, not the value; the down R is")
    print("      a bare, node-fragile empirical number with no clean (structural OR alphabet) form.")
    print("  ROOT: the down sector sits AT a circulant node (m_d node-suppressed) -- the same node that")
    print("  up's m_u has, but here it (i) breaks R=sqrt(#paths) (sqrt3 fails), (ii) makes the R-fit")
    print("  fragile, (iii) admits no clean R. A qualitative hint (the down carries the +1 'more-decoupled'")
    print("  defect, which could reduce its effective paths below sqrt3) is NOT quantitative. Piece 2 stays")
    print("  OPEN/negative; the down-R is an honest fit. (Contrast piece 1 (N_eff), which promoted cleanly.)")

    assert not consistent, "node-dressing of m_d alone must be inconsistent (the rescue fails)"
    assert nfit == 0, "node sensitivity => NO simple form fits the full down ratios comparably (down R bare-empirical)"
    assert ratio_logrms(np.sqrt(12 / 5)) > 5 * e_star, "even sqrt(12/5) fits the full ratios far worse than R* (node)"
    assert ratio_logrms(np.sqrt(3)) > 5 * e_star, "bare sqrt3 must fit the down ratios far worse than R*"
    print("\nGATES PASSED -- node-dressing inconsistent; the node hypersensitivity means NO clean form (incl")
    print("sqrt(12/5)) fits the down ratios; the down R is a bare, node-fragile empirical number. Piece 2 negative.")
    print("exit 0")


if __name__ == "__main__":
    main()
