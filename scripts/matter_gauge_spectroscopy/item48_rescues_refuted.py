#!/usr/bin/env python3
"""Item 48 (DRIFT H6): two post-hoc rescues of the rho mass examined and REFUTED.
Self-asserting -- exit 0 means both refutations hold as stated. Logged so they are not re-litigated.

Rescue 1 (bulk-Fano -> sqrt(3)*lambda_res => 874 MeV): claimed a bulk-hybridized Fano resonance at
  lambda_res ~ 1.52 (path-weight ~0.25), rescaled by sqrt(3). REFUTED: at the substrate Delta/t=1.78
  there is no such resonance -- the path-character is a band-wide smear (centre ~0, max single-state
  weight <0.15, spread sigma>2t => 'width' >1 GeV, ~13x the physical 149 MeV). [sqrt(3) also
  double-counts the already-3D-dressed eigenvalue and contradicts sec 9.1's sqrt(2)-plateau; the
  874->775 step flips the sign of sec 9.1's own +15.6 MeV dressing -- see H6 prose.]
Rescue 2 (algorithmic-ledger m_pi + 2^8 * w => 756 MeV; w=alpha*Lambda, sec 5.9): REFUTED on canon's
  own terms -- sec 5.9 (L773) defines w as cost PER ACTIVE BIT, applied as (active-bit count)*w
  (proton 1w, neutron 2w, m_d-m_u=1w). An 8-bit register has <=8 active bits, so the bit-weight
  contributes at most 8w ~ 19 MeV, not 256w=620 MeV. 2^8=256 is the STATE count (its bit/Landauer
  cost is log2(256)=8; cf. canon's own S_vac=-log2(45/256), L801); w is the alpha-suppressed
  fine-splitting currency (2.42 MeV), not a bulk-mass block (bulk masses are O(Lambda)).
"""
import numpy as np

phi = (1 + 5 ** 0.5) / 2; alpha = 1 / 137.036; LQCD = 332.0; mpi = 136.1
w = alpha * LQCD

# --- Rescue 1: no 1.52 resonance at the substrate Delta/t=1.78 (dense SC L=12) ---
L = 12; xs = range(L); C = [(x, y, z) for x in xs for y in xs for z in xs]; idx = {c: i for i, c in enumerate(C)}; N = len(C)
A = np.zeros((N, N))
for c in C:
    for d in range(3):
        for s in (1, -1):
            cc = list(c); cc[d] = (cc[d] + s) % L; A[idx[c], idx[tuple(cc)]] = 1.0
cen = L // 2; path = [idx[(cen + i, cen, cen)] for i in range(4)]
H = A.copy()
for p in path: H[p, p] += 1.78
ew, EV = np.linalg.eigh(H); wt = (EV ** 2)[path, :].sum(0)
centre = float((wt * ew).sum() / wt.sum() - 1.78)               # weighted-mean energy, rel to Delta
sigma = float(np.sqrt((wt * (ew - (centre + 1.78)) ** 2).sum() / wt.sum()))
maxw = float(wt.max())
print(f"Rescue1 (874): centre(rel)={centre:+.3f} [claim 1.52]  max single-state path-weight={maxw:.3f} [claim 0.25]")
print(f"               spread sigma={sigma:.2f}t -> 'width' ~{2.355 * sigma * LQCD:.0f} MeV  [physical rho width 149]")
assert abs(centre) < 0.5,            "no resonance pinned at 1.52: weighted centre is ~0"
assert maxw < 0.15,                  "no sharp 0.25 state: path-weight is smeared across the band"
assert 2.355 * sigma * LQCD > 1000,  "model 'width' is ~GeV, not the narrow 149 MeV"

# --- Rescue 2: the per-active-bit ceiling is 8w, not 256w ---
print(f"Rescue2 (756): w={w:.4f} MeV/active-bit ; 8-bit ceiling 8w={8 * w:.1f} ; 256w={256 * w:.1f} ; m_pi+256w={mpi + 256 * w:.1f}")
assert abs(w - 2.4227) < 1e-3
assert 8 * w < 20,                   "an 8-bit register (<=8 active bits) yields at most 8w ~ 19 MeV"
assert abs(np.log2(256) - 8) < 1e-9, "2^8=256 is the STATE count; its bit/Landauer cost is log2=8"
ncomp = sum(1 for n in range(400) if 0.98 * 760 <= mpi + n * w <= 1.02 * 760)
print(f"               competitors m_pi+n*w within 2% of the refuted 760: {ncomp}")
assert ncomp >= 10,                  "dense competitor grid (s16.3): >=10 integers within 2% of 760"

print("\nBOTH RESCUES REFUTED (exit 0). H6 robust-negative stands; the rho BULK mass is an open 3D problem.")
