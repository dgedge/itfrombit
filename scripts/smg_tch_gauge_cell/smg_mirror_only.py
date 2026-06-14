#!/usr/bin/env python3
"""
Frontier step (1): can the X-SMG interaction gap the MIRROR sector ONLY -- gapping the mirror
while leaving the physical fermion chiral, gapless, and correctly coupled?

This is the open 3+1D chiral-lattice-gauge problem; the FULL interacting spatial-x-internal
many-body theory (16 mirror Weyl x spatial sites, strong CSS interaction, dynamical gauge) is
intractable here and is NOT attempted. What IS computable -- and is the NECESSARY enabling
condition -- is the framework's domain-wall spatial separation (ANCHOR item 77 / Part 11 Q5,
bulk_domainwall_overlap.py): physical mode on one wall, mirror (doubler) on the CONJUGATE
wall, exponentially separated. If that separation is real, a mirror-wall-localized interaction
couples to the physical mode only at O(exp(-L_wall)) -> it can be mirror-selective.

This script:
  (A) builds the domain-wall edge-mode structure (SSH, the standard exactly-solvable model for
      two exponentially-separated wall zero modes) and shows the physical<->mirror coupling
      (the near-zero level splitting) decays EXPONENTIALLY with wall separation;
  (B) shows a mirror-wall-localized operator's matrix element on the PHYSICAL mode is
      exp-small -> the X-SMG interaction on the mirror wall preserves the physical mode to
      O(exp(-L));
  (C) combines with smg_construction.py: physical wall = Z-only code (chiral 16-dim
      generation, gapless); mirror wall = full CSS (unique symmetric gap); exp-separated ->
      mirror-only application realised at the FREE + LEADING-interaction level;
  (D) states honestly the parts that remain open (the genuine frontier).

NOTE: SSH is a 1D toy for the SEPARATION mechanism (the framework's wall is the 2D octagonal
domain wall of item 77); and 'free + leading-interaction' is NOT the full non-perturbative
many-body statement. Both flagged. numpy; self-asserting.
"""
import numpy as np
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)

# ----------------------------------------------------------------------------------
Hd("(A) Domain-wall separation: physical<->mirror coupling decays exp with wall length")
# SSH chain: N cells, sites A_n,B_n; intra-cell hop v, inter-cell hop w. Topological for w>v.
# Open chain -> two end zero modes (physical at left, mirror at right), exp-localized.
v, w = 0.5, 1.0                          # v<w => topological => two end modes
def ssh_H(N):
    dim=2*N; H=np.zeros((dim,dim))
    for n in range(N):
        A,B=2*n,2*n+1
        H[A,B]=H[B,A]=v                  # intra-cell
        if n+1<N:
            H[B,2*(n+1)]=H[2*(n+1),B]=w  # inter-cell
    return H
print("  SSH wall model (v=0.5,w=1.0, topological): two end zero modes (physical L / mirror R)")
print(f"  {'N cells':>7} | {'|splitting| = phys<->mirror coupling':>38} | ratio")
prev=None
def splitting(N):                        # gap between the two near-zero (end) modes
    evs=np.sort(np.linalg.eigvalsh(ssh_H(N)))
    return abs(evs[N]-evs[N-1])
for N in [4,6,8,10,12]:
    split=splitting(N)
    r=f"{prev/split:5.1f}x" if prev else ""
    print(f"  {N:>7} | {split:>38.2e} | {r}")
    prev=split
# exponential check: splitting(N) ~ (v/w)^N
N1,N2=6,12
s1,s2=splitting(N1),splitting(N2)
decay=(s2/s1)**(1/(N2-N1))
print(f"  per-cell decay factor of the coupling: {decay:.3f}  (~ v/w = {v/w}); EXPONENTIAL.")
assert decay < 0.8 and s2 < s1*1e-1
print("  => the physical<->mirror coupling is exp-suppressed in the wall separation: a mirror-")
print("     wall interaction reaches the physical mode only at O(exp(-L)).")

# ----------------------------------------------------------------------------------
Hd("(B) Mirror-selectivity: a mirror-wall operator's action on the PHYSICAL mode is exp-small")
# analytic end modes: left (physical) on A-sublattice ~ (-v/w)^n ; right (mirror) on B ~ (-v/w)^(N-1-n)
def end_modes(N):
    r=-v/w
    psiL=np.zeros(2*N); psiR=np.zeros(2*N)
    for n in range(N):
        psiL[2*n]   = r**n           # A sublattice, decays from left
        psiR[2*n+1] = r**(N-1-n)      # B sublattice, decays from right
    return psiL/np.linalg.norm(psiL), psiR/np.linalg.norm(psiR)
for N in [6,10,14]:
    psiL,psiR=end_modes(N)
    # mirror-wall-localized operator: projector onto the last cell (sites 2N-2,2N-1)
    O=np.zeros((2*N,2*N)); O[2*N-1,2*N-1]=1; O[2*N-2,2*N-2]=1
    me_phys=float(psiL@O@psiL)          # how much the mirror operator acts on the physical mode
    me_mirr=float(psiR@O@psiR)          # ... vs on the mirror mode
    ov=abs(float(psiL@psiR))
    print(f"  N={N:>2}: <phys|O_mirror|phys>={me_phys:.2e}   <mirror|O_mirror|mirror>={me_mirr:.3f}   "
          f"<phys|mirror>={ov:.2e}")
assert me_phys<1e-3 and me_mirr>0.3 and ov<1e-3
print("  => the mirror-wall SMG interaction acts O(1) on the mirror, O(exp(-L)) on the physical")
print("     mode and on their overlap -> it gaps the mirror while PRESERVING the physical mode.")

# ----------------------------------------------------------------------------------
Hd("(C) Combined picture: mirror-only application at the free + leading-interaction level")
print("""  Assembling the established pieces:
   * PHYSICAL wall: Z-only code -> 16-dim chiral codespace (the generation), gapless.
                    (rank-16 projector; smg_construction.py / mirror_invalid_subspace.py)
   * MIRROR wall  : full self-dual CSS (Z + X) -> UNIQUE, gapped (gap 2), entangled,
                    zero-bilinear-order state -- symmetric mass generation.
                    (smg_construction.py, exact)
   * SEPARATION   : the X (SMG) interaction is supported on the mirror wall; its coupling to
                    the physical mode is O(exp(-L_wall)) (A,B above).
  => at the free + leading-interaction level the mirror is symmetrically gapped and the
     physical generation survives chiral and gapless, with cross-talk O(exp(-L)). This is the
     mirror-only application -- realised to leading order, using only the framework's own
     domain wall (item 77) + CSS X-half (item 13). No new obstruction appears (and the physical
     generation alone is anomaly-free, §2.11, so gapping the mirror is anomaly-consistent).""")

# ----------------------------------------------------------------------------------
Hd("(D) VERDICT — step (1): advanced, not closed")
print("""ESTABLISHED (necessary enabling conditions, computed):
  * The framework's domain-wall structure (item 77) gives EXPONENTIAL physical<->mirror
    separation: the coupling decays ~ (v/w)^L (A). So a mirror-wall-localized SMG interaction
    is mirror-SELECTIVE -- it acts O(1) on the mirror, O(exp(-L)) on the physical mode (B).
  * Combined with the exact per-cell results (physical = Z-only chiral codespace; mirror =
    full-CSS unique symmetric gap, smg_construction.py), the mirror-only application holds at
    the FREE + LEADING-interaction level: mirror symmetrically gapped, physical preserved
    chiral/gapless, cross-talk O(exp(-L)). The physical generation alone is anomaly-free
    (§2.11), so removing the mirror is anomaly-consistent -- no new obstruction. (C)

OPEN (the genuine 3+1D frontier -- NOT closed here):
  1. Non-perturbative back-reaction: the FULL interacting many-body theory (16 mirror Weyl x
     spatial sites, STRONG CSS coupling) could still flow to the PMS (symmetry-breaking) phase
     rather than the symmetric gapped phase, despite the exact single-cell result. This is the
     field-wide open SMG dynamics question; the spatial-x-internal ED is intractable here.
  2. Dynamical gauge fields: with the SM group GAUGED, the bulk gauge field connects the walls;
     anomaly inflow must cancel correctly (physical anomaly-free helps, but the gauged-mirror
     decoupling is the standard unsolved chiral-lattice-gauge step).
  3. Same-chirality residual (walk_kernel_overlap's 4 even-corner species), untouched by R2 /
     this construction.

NET: step (1)'s ENABLING mechanism -- spatial mirror-selectivity via the item-77 domain wall --
is demonstrated, and combined with the exact per-cell symmetric gap it realises mirror-only
gapping to leading order with no new obstruction. The remaining gap is the NON-PERTURBATIVE
many-body + dynamically-gauged closure (open field-wide). So items 93/143: the SMG route now
has precondition (16), N-N-internal-absence, anti-mirror operator (Z_chi Z_W), constructed
symmetric interaction (CSS X-half / item 13), exact local symmetric gap, AND a demonstrated
spatial mirror-selectivity -- every ingredient EXCEPT the non-perturbative + gauged global
proof. That last step is the genuine frontier; it is not closed, and honestly may not be
closable without the field-wide chiral-lattice-gauge resolution. Continuum-Locked: still NO.""")
print("\nALL ASSERTS PASSED (separation + selectivity computed exactly; the non-perturbative")
print("gauged many-body closure is explicitly NOT claimed).")
