#!/usr/bin/env python3
"""
GENUINE discrete 1-loop vacuum polarization on the 16-node 2-octagon graph,
computed as a real spectral integral (Dirac-sea response to gauge flux) rather
than the §5.4 COUNTS (N1 = 2n-1 = 31 ; |C2| = S4/|H1*| = 24/7).

"Restore the integral" done honestly. Per the K3 lesson (a prior integral
approach to dressed-alpha on THIS graph produced a withdrawn false positive),
the discipline is:
  * validate the machinery on a KNOWN/gapped case (two independent routes must
    agree; the exact lattice Ward identity must hold) BEFORE trusting any number;
  * use a physical IR regulator (a fermion mass m) -- the genuine vacuum
    polarization is mass-dependent (continuum QED: Pi ~ (1/3pi) ln(Lambda/m));
  * report normalization-INDEPENDENT structure (ranks, Ward, m-scaling) over a
    single tuned number; default to REFUTE "genuine integral = the count".

Mass regulator: bipartite staggered Dirac mass M = m*Gamma, Gamma=diag((-1)^i)
on sublattices (evens/odds), so {A,Gamma}=0 but {A+mGamma,Gamma}!=0 (a true
chiral-symmetry-breaking Dirac mass). Spectrum -> +-sqrt(lambda^2+m^2): gapped.

numpy. Self-asserting on harness + exact-Hodge; TCH integral values REPORTED.
"""
import numpy as np
np.set_printoptions(suppress=True, linewidth=140)

def build_tch():
    E = []
    for k in range(8): E.append((k, (k+1) % 8))            # octagon 1
    for k in range(8): E.append((8+k, 8+(k+1) % 8))        # octagon 2
    E.append((1, 8)); E.append((9, 0))                     # rungs -> bridge square
    return E, 16

# ----------------------------------------------------------------------
print("="*74); print("PART A — exact discrete-Hodge / cochain counts (normalization-free)")
print("="*74)
EDGES, V = build_tch(); E = len(EDGES)
B = np.zeros((E, V))
for e,(i,j) in enumerate(EDGES): B[e, i] = -1.0; B[e, j] = +1.0   # gradient d0
rank_d0 = np.linalg.matrix_rank(B)
L0 = B.T @ B
dim_ker_L0 = V - np.linalg.matrix_rank(L0)
cycle_space = E - rank_d0
print(f"  V={V}, E={E}, beta_1=E-V+1={E-V+1}")
print(f"  dim ker(L0) = {dim_ker_L0}  (Ward/constant mode = # components)")
print(f"  rank(d0)    = {rank_d0}  (= V-1 independent gauge transformations)")
print(f"  cycle space = {cycle_space}  (= harmonic 1-forms = beta_1)")
print("  Framework sub-item (d): matter 32->31 via 1 Ward; edges 8 null->10 physical.")
print(f"    Ward kernel IS 1-dim (constant): matches the '-1'. But:")
print(f"    edge gauge modes = rank d0 = {rank_d0} (not 8); gauge-invariant edge")
print(f"    modes = E-rank = {cycle_space} = beta_1 (not 10). The '31' is a doubled")
print(f"    count (16x2-1); standard density/Hodge rank is 15 (=16-1).")
assert (V, E, E-V+1) == (16, 18, 3)
assert dim_ker_L0 == 1 and rank_d0 == 15 and cycle_space == 3

# ----------------------------------------------------------------------
print("\n" + "="*74)
print("PART B — genuine 1-loop vacuum polarization (Dirac-sea flux response)")
print("="*74)

def make_H(edges, V, m, flux_edge=None, theta=0.0):
    """Adjacency + staggered Dirac mass m*Gamma; optional Peierls flux on one edge."""
    H = np.zeros((V, V), complex)
    for (i, j) in edges:
        H[i, j] += 1.0; H[j, i] += 1.0
    if flux_edge is not None:
        i, j = flux_edge
        H[i, j] = np.exp(1j*theta); H[j, i] = np.exp(-1j*theta)
    for i in range(V):
        H[i, i] += m * (1.0 if i % 2 == 0 else -1.0)       # Gamma = (-1)^i
    return H

def fill(ev, tol=1e-9):
    occ = np.zeros(len(ev)); occ[ev < -tol] = 1.0
    occ[np.abs(ev) <= tol] = 0.5
    return occ

def pol_perturbative(edges, V, m, flux_edge):
    """d2E0/dth2 = dia + para (genuine current-current bubble), at theta=0."""
    H = make_H(edges, V, m, flux_edge, 0.0)
    ev, U = np.linalg.eigh(H); occ = fill(ev)
    i, j = flux_edge
    J = np.zeros((V, V), complex); J[i, j] = 1j; J[j, i] = -1j   # dH/dth
    K = np.zeros((V, V), complex); K[i, j] = -1.0; K[j, i] = -1.0  # d2H/dth2
    Jm = U.conj().T @ J @ U; Km = U.conj().T @ K @ U
    dia = sum(occ[a]*Km[a, a].real for a in range(V))
    para = 0.0
    for a in range(V):
        if occ[a] <= 0: continue
        for b in range(V):
            if occ[b] >= 1: continue
            if abs(ev[a]-ev[b]) > 1e-9:
                para += (occ[a]-occ[b])*abs(Jm[a, b])**2 * 2.0/(ev[a]-ev[b])
    return dia, para, dia+para

def pol_finitediff(edges, V, m, flux_edge, h=1e-4):
    def E0(th):
        ev = np.linalg.eigvalsh(make_H(edges, V, m, flux_edge, th))
        return np.sum(ev*fill(ev))
    e = [E0(t) for t in (-2*h, -h, 0, h, 2*h)]
    return (-e[0]+16*e[1]-30*e[2]+16*e[3]-e[4])/(12*h*h)

# ---- HARNESS 1: GAPPED ring (staggered mass), two routes must agree ----
def ring(n): return [(k, (k+1) % n) for k in range(n)], n
rE, rV = ring(12)
rp = pol_perturbative(rE, rV, 0.6, (0, 1))[2]
rf = pol_finitediff(rE, rV, 0.6, (0, 1))
print(f"  HARNESS gapped ring C12 (m=0.6): d2E0/dth2 perturbative={rp:.8f} finite-diff={rf:.8f}")
assert abs(rp-rf) < 1e-5, "ring: two routes disagree -> bug"
print("    two independent routes agree -> bubble machinery validated.")

# ---- HARNESS 2: exact lattice Ward identity (static density response Pi.1=0) ----
def density_response(edges, V, m):
    H = make_H(edges, V, m).real
    ev, U = np.linalg.eigh(H); occ = fill(ev)
    Pi = np.zeros((V, V))
    for a in range(V):
        if occ[a] <= 0: continue
        for b in range(V):
            if occ[b] >= 1: continue
            de = ev[a]-ev[b]
            if abs(de) > 1e-9:
                amp = U[:, a].conj()*U[:, b]
                Pi += (occ[a]-occ[b])*2.0/de * np.real(np.outer(amp, amp.conj()))
    return Pi, ev
PiT, evT = density_response(EDGES, V, 0.3)
ward = np.max(np.abs(PiT @ np.ones(V)))
print(f"  HARNESS Ward: max|Pi . 1| = {ward:.2e} (must be 0)")
assert ward < 1e-8, "Ward violated -> bug"
print(f"    rank(density-response Pi) = {np.linalg.matrix_rank(PiT, tol=1e-9)} (<=15, constant in kernel)")

# ---- TCH spectrum (massless) ----
ev0 = np.linalg.eigvalsh(make_H(EDGES, V, 0.0).real)
nz = int(np.sum(np.abs(ev0) < 1e-9))
print(f"\n  TCH massless adjacency spectrum:")
print("   ", np.array2string(np.sort(ev0), precision=4, max_line_width=140))
print(f"    # exact zero modes = {nz}  (=> {'gapless: IR-sensitive, needs mass m' if nz else 'gapped'})")

# ---- THE GENUINE TCH 1-LOOP VACUUM POLARIZATION vs mass m ----
print(f"\n  Genuine 1-loop vacuum polarization (flux through bridge plaquette 0-1-8-9):")
print(f"  {'m':>6} {'diamag':>10} {'para(bubble)':>13} {'total':>10} {'fd-check':>10} {'2pi|total|':>11}")
rows = []
for m in [1.0, 0.5, 0.25, 0.1, 0.05, 0.02]:
    dia, para, tot = pol_perturbative(EDGES, V, m, (1, 8))
    fd = pol_finitediff(EDGES, V, m, (1, 8))
    assert abs(tot-fd) < 1e-2, f"pert vs fd disagree at m={m}"
    rows.append((m, dia, para, tot, fd, 2*np.pi*abs(tot)))
    print(f"  {m:6.2f} {dia:10.5f} {para:13.5f} {tot:10.5f} {fd:10.5f} {2*np.pi*abs(tot):11.4f}")

# scaling diagnostic: does |total| grow like ln(1/m) (continuum QED) or saturate?
ms = np.array([r[0] for r in rows]); tots = np.array([abs(r[3]) for r in rows])
# fit tot vs ln(1/m) over the small-m rows
sl_log = np.polyfit(np.log(1/ms[-4:]), tots[-4:], 1)[0]
print(f"\n  small-m slope d|total|/d ln(1/m) = {sl_log:.4f}  "
      f"({'grows logarithmically (QED-like)' if abs(sl_log)>0.02 else 'flat/saturates'})")

print(f"\n  COMPARISON to the §5.4 COUNTS (1-loop):")
print(f"    framework count  N1 = 2n-1 = 31  (an INTEGER, m-independent)")
print(f"    genuine |total d2E0/dth2| is O(0.1-1), m-DEPENDENT, NOT 31.")
print(f"    To force agreement one sets the photon normalization so 2pi*Pi=31:")
print(f"    that is a CHOICE of coupling, not an output of the integral.")
print(f"    implied N1_eff = 2pi*|total| ranges {min(r[5] for r in rows):.2f}..{max(r[5] for r in rows):.2f}"
      f" across m -> not a fixed 31.")

print("\n" + "="*74); print("VERDICT (reported; harness asserts passed = exit 0)"); print("="*74)
print(" 1. Machinery validated: 2 routes agree (gapped ring + every TCH m); Ward exact.")
print(" 2. The genuine 1-loop vacuum polarization is a mass-dependent O(1) Dirac-sea")
print("    response, NOT the integer 31. N1=2n-1=31 is a (doubled) MODE COUNT; it is")
print("    not the value of the current-current integral. Matching N1/(2pi) requires")
print("    choosing the photon coupling normalization (a free input), so the '1-loop")
print("    = 31' identification is a counting convention, not an integral result.")
print(" 3. Exact Hodge: Ward kernel=1, gauge rank=15, cycle space=3 -- the framework's")
print("    '32->31' and '8 null ->10 physical' (sub-item d) match none of these.")
print(" 4. The genuine integral DOES carry the physical IR (mass) dependence the count")
print("    cannot: whether it is the QED log is the next, sharper question (the m-scan")
print("    above is the start). This is the honest content of 'restore the integral'.")
print("\nexit 0 == harness (2-route agreement at all m + exact Ward) and Hodge asserts pass.")
