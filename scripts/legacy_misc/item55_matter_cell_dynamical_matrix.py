#!/usr/bin/env python3
"""
item55_matter_cell_dynamical_matrix.py

Item 55 (ANCHOR §15): construct the 24x24 mechanical dynamical matrix D(k=0)
of the 8-vertex Q_3 matter cell -- which canon flags as "an open mathematical
problem" -- and extract the A_1g breathing eigenvalue (= M_H^2 per §8.1).

The matter cell is the cube: its 8 vertices are the 8 register bits of §2.1
sitting at corners (+/-1,+/-1,+/-1); the 12 cube edges are the Q_3 frustration
pairs = the §5.2 Z-type stabilizers Z_i Z_j.

Self-asserting (exit 0 iff every check passes). What it establishes:

 1. CONSTRUCTION.  Builds the full O_h (order 48) action on the 24-dim space
    (8 vertices x 3 Cartesian displacements) and verifies, by rigorous
    character projection, the canonical vibrational decomposition
        Gamma_vib = A1g + A2u + Eg + Eu + T1g + 2 T1u + 2 T2g + T2u   (dim 24).
    -> closes the "build D(0)" construction step of item 55.

 2. PARAMETER COUNT.  Confirms a general O_h-invariant D(0) carries exactly
    12 independent parameters (sum over irreps of m(m+1)/2).

 3. SUBSTRATE REDUCTION.  The substrate does NOT supply 12 free constants.
    The cube has exactly THREE O_h orbits of vertex-pairs -- 12 edges
    (the Z-stabilizers), 12 face diagonals, 4 body diagonals (antipodal /
    complement pairs) -- plus an onsite radial term.  A central-force model
    on these orbits has at most 4 parameters; the minimal substrate model
    (edges = Z-stabilizers only) has ONE.  12 -> <=4 is the item-55d reduction.

 4. A1g EIGENVALUE (the "k_1 = M_H^2" of item 55).  Closed form
        lambda_A1g = 2 k_edge + 4 k_face + 2 k_body + k_site,
    verified against numerical diagonalisation AND the A1g projector, and
    shown to be an EXACT eigenvector (A1g is 1-dim, so it must be).

 5. SCALE VERDICT (honest).  lambda_A1g is a pure number times the cell
    stiffness.  At the substrate scale Lambda_QCD the edge-only breathing
    mode is sqrt(2)*Lambda ~ 0.47 GeV -- a factor ~266 below 125 GeV.  So the
    matter-cell breathing mode CANNOT supply the absolute Higgs mass from
    Lambda_QCD; the EW stiffness scale is a mandatory separate input
    (consistent with the §16.2 Single Dimensionful Anchor no-go).  What IS
    derived is the dimensionless A1g eigenvalue and the inter-mode ratios.
"""
import itertools
import sys
import numpy as np

np.set_printoptions(precision=4, suppress=True)
TOL = 1e-9
_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# ---- 1. cube vertices and the O_h group -----------------------------------
verts = np.array(list(itertools.product([-1.0, 1.0], repeat=3)))  # (8,3)
N = 8


def oh_group():
    """All 48 signed 3x3 permutation matrices = the full octahedral group O_h."""
    mats = []
    for perm in itertools.permutations(range(3)):
        P = np.zeros((3, 3))
        for i, p in enumerate(perm):
            P[i, p] = 1.0
        for signs in itertools.product([-1.0, 1.0], repeat=3):
            mats.append(np.diag(signs) @ P)
    return mats


G = oh_group()
check("|O_h| = 48", len(G) == 48)


def vperm(g):
    """Permutation of the 8 cube vertices induced by g."""
    idx = []
    for v in verts:
        w = g @ v
        j = np.where(np.all(np.abs(verts - w) < TOL, axis=1))[0]
        idx.append(int(j[0]))
    return np.array(idx)


def rep24(g):
    """24-dim mechanical rep: (vertex permutation) (x) (3-vector rotation)."""
    perm = vperm(g)
    R = np.zeros((24, 24))
    for a in range(N):
        b = perm[a]
        R[3 * b:3 * b + 3, 3 * a:3 * a + 3] = g
    return R


reps = [rep24(g) for g in G]


# ---- 2. classify elements into O_h classes; project the 24-rep ------------
def order(g):
    M = np.eye(3)
    for n in range(1, 7):
        M = M @ g
        if np.allclose(M, np.eye(3)):
            return n
    return 0


def is_diag(g):
    return np.allclose(g - np.diag(np.diag(g)), 0)


def classify(g):
    d, t, o, dg = round(np.linalg.det(g)), round(np.trace(g)), order(g), is_diag(g)
    if (d, t, o) == (1, 3, 1):
        return 'E'
    if (d, t, o) == (-1, -3, 2):
        return 'i'
    return {
        (1, 0, 3, False): '8C3', (1, 1, 4, False): '6C4',
        (1, -1, 2, True): '3C2', (1, -1, 2, False): '6C2p',
        (-1, 0, 6, False): '8S6', (-1, -1, 4, False): '6S4',
        (-1, 1, 2, True): '3sh', (-1, 1, 2, False): '6sd',
    }[(d, t, o, dg)]


classes = [classify(g) for g in G]
class_order = ['E', '8C3', '6C2p', '6C4', '3C2', 'i', '6S4', '8S6', '3sh', '6sd']
class_size = {c: classes.count(c) for c in class_order}
check("class sizes 1,8,6,6,3,1,6,8,3,6",
      [class_size[c] for c in class_order] == [1, 8, 6, 6, 3, 1, 6, 8, 3, 6])

# standard O_h character table (columns in class_order)
chartab = {
    'A1g': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    'A2g': [1, 1, -1, -1, 1, 1, -1, 1, 1, -1],
    'Eg':  [2, -1, 0, 0, 2, 2, 0, -1, 2, 0],
    'T1g': [3, 0, -1, 1, -1, 3, 1, 0, -1, -1],
    'T2g': [3, 0, 1, -1, -1, 3, -1, 0, -1, 1],
    'A1u': [1, 1, 1, 1, 1, -1, -1, -1, -1, -1],
    'A2u': [1, 1, -1, -1, 1, -1, 1, -1, -1, 1],
    'Eu':  [2, -1, 0, 0, 2, -2, 0, 1, -2, 0],
    'T1u': [3, 0, -1, 1, -1, -3, -1, 0, 1, 1],
    'T2u': [3, 0, 1, -1, -1, -3, 1, 0, 1, -1],
}
dims = {k: v[0] for k, v in chartab.items()}

# character of the 24-rep = (# fixed vertices) * trace(g), a class function
chi24 = {}
for g, cls in zip(G, classes):
    nfix = int(np.sum(vperm(g) == np.arange(N)))
    chi24[cls] = nfix * round(np.trace(g))

mult = {}
for irr, row in chartab.items():
    s = sum(class_size[c] * chi24[c] * row[i] for i, c in enumerate(class_order))
    mult[irr] = round(s / 48)

present = {k: v for k, v in mult.items() if v}
print("Gamma_24 decomposition:", present)
canon = {'A1g': 1, 'Eg': 1, 'T1g': 1, 'T2g': 2, 'A2u': 1, 'Eu': 1, 'T1u': 2, 'T2u': 1}
check("decomposition matches canon Gamma_vib", present == canon)
check("dim sum = 24", sum(mult[k] * dims[k] for k in mult) == 24)

# ---- 3. the 12 symmetry-allowed parameters --------------------------------
nparam = sum(m * (m + 1) // 2 for m in mult.values())
print(f"independent parameters of a general O_h-invariant D(0): {nparam}")
check("12 symmetry-allowed parameters", nparam == 12)

# cube vertex-pair orbits (the substrate's structural bonds)
orbit = {1: 'edge (Z-stabilizer)', 2: 'face diagonal', 3: 'body diagonal (antipodal)'}
orbit_count = {1: 0, 2: 0, 3: 0}
for a in range(N):
    for b in range(a + 1, N):
        orbit_count[int(round(np.sum(np.abs(verts[a] - verts[b]) > TOL)))] += 1
print("vertex-pair orbits:", {orbit[k]: orbit_count[k] for k in orbit})
check("orbit counts 12 edges / 12 face / 4 body", orbit_count == {1: 12, 2: 12, 3: 4})

# ---- 4. substrate dynamical matrix + A1g eigenvalue -----------------------
def spring(D, a, b, k):
    n = verts[b] - verts[a]
    n = n / np.linalg.norm(n)
    P = k * np.outer(n, n)
    for i, j, s in ((a, a, 1), (b, b, 1), (a, b, -1), (b, a, -1)):
        D[3 * i:3 * i + 3, 3 * j:3 * j + 3] += s * P


def build_D(k_edge, k_face, k_body, k_site):
    D = np.zeros((24, 24))
    for a in range(N):
        for b in range(a + 1, N):
            diff = int(round(np.sum(np.abs(verts[a] - verts[b]) > TOL)))
            k = {1: k_edge, 2: k_face, 3: k_body}[diff]
            if k:
                spring(D, a, b, k)
        D[3 * a:3 * a + 3, 3 * a:3 * a + 3] += k_site * np.eye(3)
    return D


# A1g projector (must be rank 1: A1g has multiplicity 1)
PA1g = sum(chartab['A1g'][class_order.index(c)] * R for R, c in zip(reps, classes))
PA1g *= dims['A1g'] / 48
check("A1g projector is rank 1", abs(np.trace(PA1g) - 1) < 1e-9)
w, Vp = np.linalg.eigh(PA1g)
vA = Vp[:, int(np.argmax(w))]              # the breathing coordinate

for ke, kf, kb, ks in [(1.0, 0.0, 0.0, 0.0), (1.0, 0.5, 0.25, 0.3)]:
    D = build_D(ke, kf, kb, ks)
    lam_num = float(vA @ D @ vA / (vA @ vA))
    lam_ana = 2 * ke + 4 * kf + 2 * kb + ks
    check(f"A1g eigenvalue numeric==analytic at (ke,kf,kb,ks)=({ke},{kf},{kb},{ks}): "
          f"{lam_num:.4f}=={lam_ana:.4f}", abs(lam_num - lam_ana) < 1e-6)
    resid = np.linalg.norm(D @ vA - lam_num * vA) / np.linalg.norm(vA)
    check(f"  breathing coordinate is an exact eigenvector (resid {resid:.1e})", resid < 1e-6)

# ---- 5. irrep-resolved spectrum (the ratio calc) --------------------------
# rigorous: diagonalise D inside each irrep isotypic subspace (D commutes with
# every rep, so it block-diagonalises; this avoids cross-irrep degeneracy
# mis-assignment when several irreps share an eigenvalue).
projs = {}
for irr, row in chartab.items():
    P = sum(row[class_order.index(c)] * R for R, c in zip(reps, classes))
    projs[irr] = P * dims[irr] / 48


def irrep_eigs(D, irr):
    w, U = np.linalg.eigh(projs[irr])
    B = U[:, w > 0.5]                       # orthonormal basis of the isotypic space
    if B.shape[1] == 0:
        return []
    e = np.linalg.eigvalsh(B.T @ D @ B)
    return sorted({0.0 if abs(x) < 1e-9 else round(float(x), 4) for x in e})


def show(D, label):
    print(f"  [{label}]")
    table = {}
    for irr in chartab:
        e = irrep_eigs(D, irr)
        if e:
            table[irr] = e
            print(f"    {irr:4s}: {e}")
    return table


print("\n--- edge-only model (k_edge=1 = the Z-stabilizers): A1g NOT isolated ---")
t_edge = show(build_D(1.0, 0.0, 0.0, 0.0), "edge-only")
check("edge-only: the only distinct nonzero eigenvalue is 2 (A1g degenerate)",
      sorted({x for v in t_edge.values() for x in v if x > 0}) == [2.0])

print("\n--- full substrate model (k_edge=k_face=k_body=1, k_site=0) ---")
t_full = show(build_D(1.0, 1.0, 1.0, 0.0), "all bonds equal")
lamA = t_full['A1g'][-1]
allvals = [x for v in t_full.values() for x in v]
check(f"adding face+body bonds splits A1g out: lambda_A1g={lamA} is now non-degenerate",
      sum(abs(lamA - x) < 1e-6 for x in allvals) == 1)
lamEg = [x for x in t_full['Eg'] if x > 1e-9]
if lamEg:
    print(f"  within-cell ratio  M_H^2/M_Eg^2 = lambda_A1g/lambda_Eg "
          f"= {lamA}/{lamEg[-1]} = {lamA / lamEg[-1]:.4f}   (Higgs vs the Eg graviton mode)")

LAMBDA = 0.3317                            # GeV, the single dimensionful anchor (Lambda_p)
mH_edge = np.sqrt(2.0) * LAMBDA
print(f"\nIf cell stiffness = Lambda_QCD: M_H(A1g, edge) = sqrt(2)*Lambda "
      f"= {mH_edge * 1000:.1f} MeV")
print(f"Measured 125 GeV is {125.0 / mH_edge:.0f}x larger -> the EW stiffness "
      f"scale is a mandatory separate input (consistent with §16.2).")
check("125 GeV is NOT reproducible from a Lambda_QCD breathing mode (ratio>100)",
      125.0 / mH_edge > 100)

print("\nDERIVED (dimensionless): Gamma_vib, the 12->{<=4} parameter reduction, "
      "and lambda_A1g = 2 k_edge + 4 k_face + 2 k_body + k_site.")
print("STILL OPEN: the few substrate stiffnesses from Lambda_QCD (55d index "
      "theorem); Delta_1=1/28 dressing (55c); cross-space M_H/M_Z (55d).")
print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
