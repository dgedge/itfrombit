#!/usr/bin/env python3
"""
weak_bridge_audit.py
====================
Reproducibility for the weak-bridge transmission audit (§6.6 / §15 item 38).

§6.6 claims the weak parity-flip transmission across the 4-vertex square gauge
bridge is |A_W|^2 = 1/256 = 1/2^8 (the inverse 8-bit state space). Item 38 flags
that this rests on a deterministic localised-sweep postulate, not the canonical
coherent unitary walk W. This script confirms the flag:

  - a COHERENT continuous-time quantum walk on the C4 bridge transmits entry->exit
    with probability 1 (antipodal perfect state transfer; C4 = Q2 has PST), NOT 1/256;
  - 1/256 is the CLASSICAL independent-blocking product (1/4 per vertex)^4, and
    (1/4)^4 = (1/2^2)^4 = 1/2^8 identically -- so the "= inverse 8-bit state space"
    coincidence is trivial arithmetic, not a derived quantum amplitude.

E4 scope: this does NOT prove no walk can give 1/256 (a decohering/measurement-gated
process gives exactly the classical product). It shows the *natural coherent* walk
gives O(1), so 1/256 is the classical-token value and the coherent-W derivation is
open. numpy + scipy.
"""
import numpy as np
import scipy.linalg as la

fails = []
def check(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    if not c: fails.append(m)

A = np.array([[0, 1, 0, 1], [1, 0, 1, 0], [0, 1, 0, 1], [1, 0, 1, 0]], float)  # C4: entry=0, exit=2 (antipode)
ts = np.linspace(0, 2*np.pi, 4000)
P2 = np.array([abs((la.expm(-1j*A*t) @ np.array([1, 0, 0, 0], complex))[2])**2 for t in ts])
print(f"coherent CTQW on C4: max P(entry 0 -> exit 2) = {P2.max():.4f} at t = {ts[P2.argmax()]:.3f}  (pi/2 = {np.pi/2:.3f})")
check(P2.max() > 0.999, "coherent walk gives PERFECT transmission (P=1) -- antipodal perfect state transfer")
check(P2.max() - 1/256 > 0.5, f"coherent transmission (~1) is NOT the classical 1/256 = {1/256:.5f}")

print(f"classical localised token: (1/4 per vertex)^4 = {(1/4)**4:.6f} = 1/256 = 1/2^8 = {1/2**8:.6f}")
check(abs((1/4)**4 - 1/2**8) < 1e-15,
      "(1/4)^4 = 1/2^8 identically -> '1/256 = inverse 8-bit state space' is trivial arithmetic")

print("\n=> 1/256 is the classical independent-blocking product, NOT a coherent-walk amplitude;")
print("   the coherent-W derivation of |A_W|^2 remains open (§15 item 38).")
print("=" * 64)
import sys
if fails:
    print(f"RESULT: {len(fails)} FAILED"); [print("  -", f) for f in fails]; sys.exit(1)
print("RESULT: exit 0 -- coherent C4 walk transmits perfectly; 1/256 is the classical token product.")
print("=" * 64)
