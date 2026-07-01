#!/usr/bin/env python3
r"""Magic budget of the [8,4,4]=RM(1,3) record cell: where the poly(n) wall is.

Context
-------
The record-grammar tractability question is "is the generator-form ledger
poly(n)?".  The edge-loop/transfer results showed the *Wilson* ledger is a 2-D
tensor network with bond dimension 2^width.  But spatial dimension is not the
obstruction: by Gottesman--Knill, stabilizer states are classically poly in ANY
dimension.  The real cost driver is non-Clifford ("magic") content.  This script
measures it directly on the framework's own substrate cell.

Measure
-------
Stabilizer nullity nu (Beverland et al.): a magic monotone with nu(|psi>)=0 iff
|psi> is a stabilizer state.  For an N-qubit pure state,

    #{(a,b) : |<psi| X^a Z^b |psi>| = 1} = 2^{N - nu},

so nu = N - log2(that count).

Result asserted here
--------------------
1. The bare record cell |C> = uniform superposition over RM(1,3)=[8,4,4]
   codewords (the self-dual CSS [[8,0,4]] "complete stabilizer record") has
   nu = 0: it is exactly a stabilizer state.  The substrate's record / syndrome /
   error-correction layer is therefore poly(n) FOR FREE, in any dimension.
2. Magic is opt-in: each non-Clifford (T) injection raises nu by at most one, and
   nu saturates at N/2 = 4 because T commutes with Z so the 4 Z-type stabilizer
   generators always survive.  (Transversal T^{(8)} on this doubly-even self-dual
   code is special -- nu dips back to 3 -- which is exactly why such codes admit a
   transversal T gate.)
3. The stabilizer-rank cost of a state with nullity nu is ~2^{nu}.  So a system
   of cells is poly(n) iff the total magic is sparse: O(log n) injected cells.
   Extensive magic (a finite fraction of cells injected) is exponential -- and
   that is where the strong-sector dynamics (the SU(3) Wilson loop) lives.

Conclusion: the framework does not get universal poly(n) tractability (it cannot,
without collapsing complexity classes).  It gets a stabilizer substrate that is
poly for free, plus an explicit, measurable magic budget that says exactly which
states are tractable.  The open "is it poly?" becomes the concrete "is the magic
sparse?".
"""
from __future__ import annotations

import numpy as np

N = 8
DIM = 1 << N


def rm13_codewords():
    """RM(1,3)=[8,4,4]: span of the all-ones word and the three coordinate words."""
    pos = np.arange(8)
    gens = [np.ones(8, int), (pos >> 0) & 1, (pos >> 1) & 1, (pos >> 2) & 1]
    code = set()
    for m in range(16):
        c = np.zeros(8, int)
        for i in range(4):
            if (m >> i) & 1:
                c ^= gens[i]
        code.add(tuple(int(v) for v in c))
    return sorted(code)


def code_state(code):
    psi = np.zeros(DIM, dtype=complex)
    for c in code:
        psi[sum(b << p for p, b in enumerate(c))] = 1.0
    return psi / np.linalg.norm(psi)


_HAD = np.array([[1 - 2 * (bin(b & x).count("1") & 1) for x in range(DIM)] for b in range(DIM)], dtype=float)
_XOR = np.arange(DIM)[None, :] ^ np.arange(DIM)[:, None]


def stabilizer_nullity(psi):
    conj = np.conj(psi)
    count = 0
    for a in range(DIM):
        c = _HAD @ (conj[_XOR[a]] * psi)        # c[b] = <psi| X^a Z^b |psi>
        count += int(np.sum(np.abs(c) > 1 - 1e-9))
    return N - int(round(np.log2(count))), count


def t_inject(psi, k):
    w = np.exp(1j * np.pi / 4.0)
    x = np.arange(DIM)
    out = psi.copy()
    for q in range(k):
        out = out * np.where((x >> q) & 1, w, 1.0)
    return out


def assert_true(name, value):
    print(f"  {name:<66s} value={value}")
    if not value:
        raise AssertionError(name)


def main():
    print("Magic budget of the [8,4,4]=RM(1,3) record cell")
    print("=" * 92)
    code = rm13_codewords()
    weights = sorted({sum(c) for c in code})
    self_dual = all((sum(a & b for a, b in zip(u, v)) % 2) == 0 for u in code for v in code)
    d_min = min(sum(c) for c in code if sum(c) > 0)
    print(f"  code: |C|={len(code)}  weights={weights}  d_min={d_min}  self-dual={self_dual}")

    psi = code_state(code)
    nus = [stabilizer_nullity(t_inject(psi, k))[0] for k in range(N + 1)]
    print("\n  k (T injections):   " + "  ".join(f"{k}" for k in range(N + 1)))
    print("  nu (magic monotone):" + "  ".join(f"{n}" for n in nus))
    print("  stab-rank cost 2^nu:" + " ".join(f"{2**n:>3d}" for n in nus))

    print("\n[checks]")
    assert_true("code is [8,4,4]: |C|=16, d_min=4", len(code) == 16 and d_min == 4)
    assert_true("code is self-dual and doubly-even (weights in {0,4,8})", self_dual and set(weights) <= {0, 4, 8})
    assert_true("BARE record cell is a stabilizer state (nu=0) -> poly, free", nus[0] == 0)
    assert_true("one T injection is detected as magic (nu=1)", nus[1] == 1)
    assert_true("each injection adds at most one (nu[k] <= k)", all(nus[k] <= k for k in range(N + 1)))
    assert_true("magic saturates at N/2=4 (Z-sector survives T)", max(nus) == N // 2)
    assert_true("transversal T^(8) dips below saturation (doubly-even code)", nus[8] < max(nus))

    print(
        """
VERDICT:
  PASS.  The framework's record cell is, on the nose, a stabilizer state
  (nu=0): the record/syndrome/error-correction layer is poly(n) for free in any
  spatial dimension (Gottesman-Knill).  Magic is opt-in and bounded per cell
  (nu <= N/2), and the classical cost of a state is ~2^{total nu}.  Hence the
  record sector is poly(n) iff the magic is sparse -- O(log n) injected cells --
  and exponential when magic is extensive, which is exactly where the
  strong-sector (SU(3) Wilson-loop) dynamics sits.  The poly(n) question is thus
  resolved into a concrete, measurable condition: sparse magic, not dimension.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
