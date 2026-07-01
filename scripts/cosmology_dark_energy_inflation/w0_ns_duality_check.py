#!/usr/bin/env python3
"""
Check the dark-energy w0 prediction (ANCHOR §15 item 131):
  w(a) = -1 + a/28 ;  w0 = -27/28 = -0.9643 ;  n_s = 1 - 1/28 = 27/28 ;  w0 = -n_s.
28 = 2 (transverse modes) x 14 (weight-4 codewords of [8,4,4]); Delta_1 = 1/28.

Questions:
 1. Reproduce w0 = -27/28, n_s = 27/28.
 2. Is the "w0 = -n_s duality" a coincidence or automatic? (automatic for ANY Delta_1)
 3. Is 28 = 2x14 built from real code invariants? (verify the [8,4,4] weight enumerator)
 4. n_s vs Planck (0.9649+/-0.0042) and the §16.3 competitor (28/29 equally close).
 5. w0 vs DESI DR2 (-0.76+/-0.06), DESI 2024 (-0.73+/-0.07), and LCDM (-1).
 6. The framework's w0 history (three values) and the implied CPL (w0,wa).

Self-asserting. exit 0 == all checks pass. PDG/Planck/DESI cited inline.
"""
import itertools
D1 = 1/28
w0 = -1 + D1
ns = 1 - D1
print(f"Delta_1 = 1/28 = {D1:.5f}")
print(f"w0 = -1 + 1/28 = {w0:.5f}  (= -27/28 = {-27/28:.5f})")
print(f"n_s = 1 - 1/28 = {ns:.5f}  (=  27/28 = {27/28:.5f})")
print(f"w0 = -n_s ?  {abs(w0 + ns) < 1e-12}")

print("\n2) Is w0 = -n_s a coincidence or automatic?")
for d in [1/28, 1/10, 0.2, 1/56]:
    print(f"   Delta_1={d:.4f}:  w0=-1+D1={-1+d:+.4f}, n_s=1-D1={1-d:.4f}, w0+n_s={(-1+d)+(1-d):+.0e}")
print("   => w0 = -(1-Delta_1) = -n_s holds for ANY Delta_1. The 'duality' is a STRUCTURAL")
print("      identity of the two formula forms (DE: -1+a*D1 ; inflation: 1-D1), NOT an")
print("      independent numerical coincidence. Content = whether D1=1/28 + both forms hold.")

print("\n3) Is 28 = 2 x 14 built from real code invariants? Verify [8,4,4] weight enumerator:")
# extended Hamming [8,4,4]: generator (systematic) — count codewords by weight
G = [0b10000111,0b01001011,0b00101101,0b00011110]  # a standard [8,4,4] generator (rows)
def wt(x): return bin(x).count("1")
words = []
for bits in itertools.product([0,1],repeat=4):
    cw = 0
    for b,row in zip(bits,G):
        if b: cw ^= row
    words.append(cw)
from collections import Counter
wd = Counter(wt(w) for w in words)
print(f"   weight distribution of the 16 codewords: {dict(sorted(wd.items()))}")
print(f"   weight-4 codewords = {wd[4]}  (claim: 14);  transverse modes = 2 (photon) ; 2x14 = {2*wd[4]}")
assert wd[0]==1 and wd[4]==14 and wd[8]==1, "weight enumerator != 1+14x^4+x^8"
assert 2*wd[4]==28, "2x14 != 28"
print("   => 28 = 2 x 14 IS built from real invariants (verified). But the IDENTIFICATIONS")
print("      'n_s = 1 - 1/(channel count)' and 'w0 = -1 + 1/(channel count)' are non-standard")
print("      ansaetze mapping code-channel-counts to cosmological observables — unverified.")

print("\n4) n_s vs Planck 2018 (0.9649 +/- 0.0042):")
ns_pl, ns_u = 0.9649, 0.0042
for label,val in [("27/28",27/28),("28/29",28/29),("26/27",26/27)]:
    print(f"   {label} = {val:.5f}  ->  {abs(val-ns_pl)/ns_u:.2f} sigma from Planck")
print("   => 27/28 matches at 0.14 sigma, but 28/29 is equally close (0.15 sigma) on the other")
print("      side — the DATA does not uniquely pick 27/28; the 28=2x14 CONSTRUCTION does.")
print("      Also: 1-n_s=1/28 is the standard slow-roll 2/N with N=56 (generic inflation range).")

print("\n5) w0 = -0.9643 vs measurements:")
for label,val,unc in [("DESI DR2 2025 (CPL, +DES+Planck)",-0.76,0.06),
                       ("DESI 2024 BAO+SN",-0.73,0.07),
                       ("LCDM (cosmological constant)",-1.0,0.06)]:
    print(f"   vs {label:<34} {val:+.2f}+/-{unc}: {abs(w0-val)/unc:.2f} sigma")
print("   => ~3.3-3.4 sigma from DESI's evolving-DE central value, but only ~0.6 sigma from")
print("      LCDM. The prediction is OBSERVATIONALLY INDISTINGUISHABLE from a cosmological")
print("      constant at current precision; its falsifiable content is 'DESI's evolving DE is")
print("      wrong' — a bet on the LCDM-adjacent side, not a framework-unique signature.")

print("\n6) CPL form of w(a) = -1 + a/28  and the framework's w0 history:")
wa = -D1   # w(a)=-1+a/28 = w0 + wa(1-a) with w0=-27/28, wa=-1/28
print(f"   CPL: (w0, wa) = ({w0:.4f}, {wa:.4f}) = (-27/28, -1/28): mild THAW, never phantom.")
print(f"   DESI DR2 prefers (w0,wa) ~ (-0.76, -0.77): phantom-crossing, strong evolution -> opposite.")
print("   w0 HISTORY (DRIFT M10/M11): -3/4 (rule-class) -> -0.761 (item-118 fluid, matched DESI")
print("   at 0.45 sigma!) -> -27/28 (item-131, 28-channel; 3.4 sigma from DESI). The framework")
print("   ABANDONED the DESI-matching value for the substrate-derived one — not robustly pinned.")

assert abs(w0 + 27/28) < 1e-12 and abs(ns - 27/28) < 1e-12
assert abs(w0 - (-0.76))/0.06 > 3, "should be >3 sigma from DESI DR2"
assert abs(w0 - (-1.0))/0.06 < 1, "should be <1 sigma from LCDM"
print("\nexit 0 == arithmetic verified; duality is structural; 28=2x14 real; n_s has a §16.3")
print("competitor; w0 is ~3.4 sigma from DESI but ~LCDM-indistinguishable.")
