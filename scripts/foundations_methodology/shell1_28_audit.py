#!/usr/bin/env python3
"""Shell-1 C(8,2)=28 convention audit (the caveat from the item-79 closure, DRIFT K9).

Question: is Shell-1's NO-diagonal pair count C(8,2)=28 vs Shell-2's WITH-diagonal T(16)=136
a per-shell convenience (16.3 formula-freedom), or are both conventions forced by their objects?

Method: build AG(3,2) on the F2^3-labelled octet (the 8 faces/directions) exactly and classify
every 28 in canon. Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

pts = list(range(8))                                   # F2^3 points = octant labels
def xor(a,b): return a ^ b

# --- lines of AG(3,2): cosets of 1-dim subspaces = {x, x+v}, v != 0 ---
lines = set()
for x in pts:
    for v in range(1,8):
        lines.add(frozenset((x, xor(x,v))))
assert len(lines) == 28
# every 2-subset IS a line (over F2 the line through x,y is exactly {x,y}):
assert all(frozenset(p) in lines for p in itertools.combinations(pts,2))
# a "diagonal" {x,x} is not a line: lines require two distinct points BY DEFINITION
print(f"AG(3,2): lines = {len(lines)} = C(8,2) — every 2-subset is a line; self-lines do not exist.")
print( "  => Shell-1's 'no diagonal' is DEFINITIONAL (a line/pair-channel needs 2 distinct points),")
print( "     not a counting convention that could have been chosen otherwise.")

# --- hyperplanes: cosets of 2-dim subspaces ---
subs2 = []
for v1 in range(1,8):
    for v2 in range(v1+1,8):
        S = {0, v1, v2, xor(v1,v2)}
        if S not in subs2: subs2.append(S)
assert len(subs2) == 7
hyps = set()
for S in subs2:
    for x in pts:
        hyps.add(frozenset(xor(x,s) for s in S))
assert len(hyps) == 14
on = {p: sum(1 for h in hyps if p in h) for p in pts}
assert all(v == 7 for v in on.values())                 # each point on 7 of the 14
print(f"  hyperplanes = {len(hyps)} (7 directions x 2 cosets); each point lies on 7.")

# --- the item-131 instrument: flags and channels ---
flags = [(p,h,m) for p in pts for h in hyps if p in h for m in (0,1)]
channels = {(h,m) for (p,h,m) in flags}
pre = {c: sum(1 for f in flags if (f[1],f[2])==c) for c in channels}
assert len(flags) == 112 and len(channels) == 28 and all(v==4 for v in pre.values())
print(f"  item-131 instrument reproduced: flags 8x7x2 = {len(flags)}, channels 14x2 = {len(channels)}, "
      f"preimage 4 each.")

# --- the two 28s are structurally DISTINCT objects ---
# lines factorise as 7 directions x 4 parallel lines; channels as 14 hyperplanes x 2 modes.
by_dir = {}
for L in lines:
    a,b = tuple(L); by_dir.setdefault(xor(a,b), []).append(L)
assert len(by_dir) == 7 and all(len(v)==4 for v in by_dir.values())
# each line lies in exactly 3 hyperplanes (no 1-1 line<->channel bijection is canonical)
in_h = [sum(1 for h in hyps if L <= h) for L in lines]
assert all(n == 3 for n in in_h)
print(f"  TWO distinct 28s in canon: lines (7 dirs x 4) vs service channels (14 hyps x 2 modes);")
print(f"  each line lies in 3 hyperplanes — no canonical bijection. Numerically coincident,")
print(f"  structurally different. CONFLATION HAZARD flagged for future cross-references.")

print("""
VERDICT (closes the K9 Shell-1 caveat):
  Shell-1's 28 counts LINES/CHANNELS on the 8-point register geometry — objects that require
  two distinct points by definition; there is no 'with-diagonal' alternative to have chosen.
  Shell-2's 136 counts FERMION PAIR STATES, where the 3.5 coin legalises the diagonal
  (on-site singlets) and Fermi statistics forces Sym^2(16) = T(16).
  The per-shell convention difference is OBJECT-FORCED, not 16.3 formula-freedom: the
  formula-freedom caveat on the shell-capacity table DISSOLVES.
  One residual (minor): 5.4's gloss 'first-order loop pairs' is loose — as area-loops, 4 of
  the 28 direction-pairs are antipodal/degenerate; the clean statement is 'AG(3,2) lines on
  the F2^3-labelled octet'. Recommend the gloss fix; no number changes.
ALL ASSERTS PASSED""")
