#!/usr/bin/env python3
r"""ITEM 126 — does the channel-ledger 3/14 count BARYON-NUMBER VIOLATION, or just colour-silence?

The channel-ledger reposing (item126_channel_ledger.py / item126_weight4_event_algebra.py) dissolved
the old "11 weight-4 states not 14" blocker by reading the 14 as AG(3,2) hyperplane CHANNELS (not
codewords) and the numerator 3 as the channels avoiding the 2-bit colour sector {C0,C1} ("colour-
silent" = unvetoed by confinement). eta = (3/14) alpha^4 then matches Planck.

THE CRACK (this script): the framework's own event-algebra script discriminates 3/14 (matches) against
2/14 (excluded > 30 sigma), where 2/14 is labelled "LQ-flipping colourless only". But LQ-flip IS
baryon-number change: B = LQ/3 per fermion (quark LQ=1 -> B=1/3; lepton LQ=0 -> B=0), so a weight-4
fault changes baryon number IFF it flips the LQ bit. Therefore:
  * 2/14 = colour-silent AND LQ-flipping = the channels that survive confinement AND actually change B
           = the BARYON-NUMBER-VIOLATING surviving channels  -> this is the standard baryogenesis count;
  * 3/14 = colour-silent (survives confinement) REGARDLESS of LQ -> includes ONE channel that does NOT
           flip LQ, i.e. is BARYON-CONSERVING (Delta B = 0).
So the entire 3/14-vs-2/14 difference -- the entire match-vs-exclusion -- is ONE baryon-CONSERVING
channel. The 3/14 match hinges on counting a channel that produces no baryon-number change.

We verify: build the 14 channels, find the 3 colour-silent, show exactly 2 flip LQ (Delta B != 0) and
1 does not (Delta B = 0 on every 48-set transition); then 3/14 matches Planck, 2/14 (the actual
B-violating count) is excluded at ~50 sigma. Self-asserting; exit 0. numpy/itertools only.
"""
import itertools, numpy as np
from collections import Counter

NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
LQ, C0, C1 = 2, 3, 4
alpha0 = 1/137; eta_obs, deta = 6.12e-10, 0.04e-10

# ---------- the 14 AG(3,2) hyperplane channels (same construction as canon scripts) ----------
def bits3(x): return ((x >> 2) & 1, (x >> 1) & 1, x & 1)
hyps = sorted({frozenset(p for p in range(8)
               if (sum(u*v for u, v in zip(bits3(p), bits3(a))) % 2) == b)
               for a in range(1, 8) for b in (0, 1)}, key=sorted)
assert len(hyps) == 14 and all(len(h) == 4 for h in hyps)

# ---------- the 48-set + baryon number ----------
def valid(c):
    g0, g1, lq, c0, c1, i3, chi, w = c
    return not (g0 and g1) and w == chi and ((lq == 0) == (c0 == 0 and c1 == 0))
def tup(i): return tuple((i >> (7 - k)) & 1 for k in range(8))
P48 = [i for i in range(256) if valid(tup(i))]
assert len(P48) == 48
def B_of(i): return tup(i)[LQ] / 3.0            # baryon number of a single fermion codeword = LQ/3

# ---------- classify the colour-silent channels by LQ-flip (= baryon-number change) ----------
def mask(h): return sum(1 << (7 - p) for p in h)
colour_silent = [h for h in hyps if not (h & {C0, C1})]
assert len(colour_silent) == 3
print("the 3 colour-silent channels (avoid C0,C1 -> unvetoed by confinement):")
n_Bviol = 0
for h in colour_silent:
    flips_LQ = LQ in h
    m = mask(h)
    # baryon-number change of the bare weight-4 fault c -> c^h (the deep bypass; LQ-flip lands in Q)
    dBs = [B_of(i ^ m) - B_of(i) for i in P48]
    maxdB = max(abs(x) for x in dBs) if dBs else 0.0
    n_Bviol += flips_LQ
    print(f"   {{{', '.join(NAMES[p] for p in sorted(h))}:>}}  flips_LQ={flips_LQ}  "
          f"max|Delta B| of the fault = {maxdB:.3f}  -> "
          f"{'BARYON-VIOLATING (|Delta B|=1/3)' if flips_LQ else 'BARYON-CONSERVING (Delta B = 0)'}")
print(f"   => {n_Bviol} of 3 colour-silent channels flip LQ (change baryon number); "
      f"{3-n_Bviol} is baryon-CONSERVING.")
check_list = []
def chk(name, cond):
    check_list.append((name, bool(cond)))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")

chk("exactly 2 of the 3 colour-silent channels flip LQ (= change baryon number)", n_Bviol == 2)
# the baryon-conserving one truly has Delta B = 0 on every transition
for h in colour_silent:
    if LQ not in h:
        m = mask(h)
        chk("the non-LQ colour-silent channel changes baryon number on NO state (Delta B = 0 always)",
            all(abs(B_of(i ^ m) - B_of(i)) < 1e-12 for i in P48))

# ---------- the readings vs Planck ----------
print(f"\nreadings (eta = w * alpha_0^4 vs Planck {eta_obs:.2e} +/- {deta:.0e}):")
for nm, w in [("3/14  colour-silent / survives confinement (framework reading)", 3/14),
              ("2/14  colour-silent AND LQ-flip = BARYON-VIOLATING (standard reading)", 2/14)]:
    eta = w * alpha0**4
    print(f"   {nm:<58} eta={eta:.3e}  ({(eta-eta_obs)/deta:+6.1f} sigma)")
chk("3/14 (colour-silence) matches Planck (< 1 sigma)", abs((3/14)*alpha0**4 - eta_obs)/deta < 1)
chk("2/14 (the actual baryon-violating count) is EXCLUDED (> 30 sigma)",
    abs((2/14)*alpha0**4 - eta_obs)/deta > 30)

_ok = all(c for _, c in check_list)
print(f"""
--- VERDICT (item 126 magnitude: the 3/14 counts colour-silence, not baryon violation) ---
The channel-ledger reposing genuinely fixed the OLD mis-attribution (the 14 is the hyperplane-CHANNEL
count, not ideal-code words; the state-count blocker is dissolved). The alpha^4 SCALE (code distance
d=4) and the 14 (= item-131 clock channels) are solid. BUT the NUMERATOR 3 is the colour-SILENT count,
and only 2 of those 3 channels flip LQ -- i.e. actually change baryon number (B = LQ/3). The third
colour-silent channel is BARYON-CONSERVING (Delta B = 0 on every 48-set transition, verified). So:
  * the standard baryogenesis count -- channels that survive confinement AND violate baryon number --
    is 2/14, which the framework's OWN discrimination excludes at ~50 sigma;
  * the matching 3/14 is the colour-SILENCE (survival) count, which equates 'survives confinement' with
    'produces baryon asymmetry' and so includes one Delta B = 0 channel.
The ENTIRE 3/14-vs-2/14 gap -- the entire match-vs-exclusion -- is that one baryon-conserving channel.
NET: the eta magnitude rests on the COLOUR-VETO PREMISE (colour-silence = baryon production), NOT on a
baryon-number-violation count. Under the standard B-violation reading (LQ-flip) the prediction is 2/14
= 4.05e-10, excluded. Honest tier: alpha^4 scale derived; 14 channels derived; the 3 is a colour-silence
count that matches ONLY if colour-silence (not Delta B) is the branching criterion -- an undischarged,
physically non-standard premise. The magnitude match is premise-conditional / data-selected, confirming
'softer ground than the bare number suggests'.
""")
print("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED")
import sys; sys.exit(0 if _ok else 1)
