#!/usr/bin/env python3
"""
item55f_cauchy_born_homogenization.py

Item 55f-ii (ANCHOR §15): the Cauchy-Born dimensional bridge.  ATTEMPT, not a
claimed closure.

Background.  ew_cauchy_born_check.py REFUTES a parameter-free M_Z derivation,
citing free inputs a_0, k_t/k_l, and g.  It never built the homogenisation.
This script builds it and tries the one lever that refutation did not consider:
the framework's OWN Lorentz-invariance requirement (the 4.8.8 reduced-dispersion-
anisotropy theorem, ANCHOR L458) demands the homogenised continuum be ELASTICALLY
ISOTROPIC.  For a square-symmetric lattice that is one equation,
    C11 - C12 - 2 C66 = 0,
which FIXES the bending/stretching ratio k_t/k_l -- removing one of the free
inputs the refutation relied on.

Model (canon's 55e Maxwell-truss framing): the 4.8.8 gauge web as a beam truss
with longitudinal bond stiffness k_l and transverse (bond-bending) stiffness k_t
on every edge.  Cauchy-Born homogenisation WITH internal relaxation gives the 2D
elastic tensor C(k_l, k_t); we solve the isotropy condition for rho = k_t/k_l.

This is honest about scope: Cauchy-Born fixes the ABSOLUTE spring scale from v
(the one anchor, §16.2) and the geometry fixes the DIMENSIONLESS ratios.  Whether
that propagates to M_Z/M_H is addressed in the companion follow-up; here we
establish the homogenisation + the isotropy-fixed k_t/k_l.

Self-asserting (exit 0 iff checks pass).  numpy only.
"""
import sys
import numpy as np

np.set_printoptions(precision=5, suppress=True)
_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# ===========================================================================
# 1. 4.8.8 truncated-square geometry (edge length 1, square cell L = 1 + sqrt2)
# ===========================================================================
S = (1.0 + np.sqrt(2.0)) / 2.0          # octagon half-width
L = 1.0 + np.sqrt(2.0)                   # cell side = 2S
r = 1.0 / np.sqrt(2.0)                   # diamond half-diagonal (side-1 square)

# 4 basis vertices (the small "diamond" square, also octagon vertices)
R = np.array([
    [S - r, S],        # v0  (0.5,   1.2071)
    [S, S - r],        # v1  (1.2071,0.5)
    [S + r, S],        # v2  (1.9142,1.2071)
    [S, S + r],        # v3  (1.2071,1.9142)
])
N = 4
A1, A2 = np.array([L, 0.0]), np.array([0.0, L])

# edges: (i, j, cell-shift n1, n2);  bond vector b = R[j] + n1*A1 + n2*A2 - R[i]
edges = [
    (0, 1, 0, 0), (1, 2, 0, 0), (2, 3, 0, 0), (3, 0, 0, 0),   # 4 E_48 (diamond)
    (0, 2, -1, 0),                                             # E_88 horizontal
    (1, 3, 0, -1),                                             # E_88 vertical
]


def bond_vec(e):
    i, j, n1, n2 = e
    return R[j] + n1 * A1 + n2 * A2 - R[i]


# geometry self-checks
lengths = [np.linalg.norm(bond_vec(e)) for e in edges]
check("all 6 edges have length 1", np.allclose(lengths, 1.0))
coord = np.zeros(N, dtype=int)
for i, j, *_ in edges:
    coord[i] += 1
    coord[j] += 1
check("every vertex has coordination 3", np.all(coord == 3))
check("6 edges per cell (4 E48 + 2 E88)", len(edges) == 6)

# ===========================================================================
# 2. generic Cauchy-Born homogeniser (validated on controls below)
# ===========================================================================
def homogenise(Rv, edge_list, a1, a2, k_l, k_t):
    """2D elastic tensor (C11,C12,C66) by Cauchy-Born with internal relaxation.
    Rv: (n,2) basis positions; edge_list: (i,j,n1,n2); a1,a2 lattice vectors."""
    n = len(Rv)
    Acell = abs(a1[0] * a2[1] - a1[1] * a2[0])

    def bvec(e):
        i, j, n1, n2 = e
        return Rv[j] + n1 * a1 + n2 * a2 - Rv[i]

    Ke, A = [], np.zeros((2 * n, 2 * n))
    for e in edge_list:
        b = bvec(e)
        u = b / np.linalg.norm(b)
        t = np.array([-u[1], u[0]])
        K = k_l * np.outer(u, u) + k_t * np.outer(t, t)
        Ke.append(K)
        i, j = e[0], e[1]
        for (a, c, s) in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
            A[2 * a:2 * a + 2, 2 * c:2 * c + 2] += s * K
    Apinv = np.linalg.pinv(A, rcond=1e-10)

    def W(eps):
        E0, g = 0.0, np.zeros(2 * n)
        for e, K in zip(edge_list, Ke):
            du = eps @ bvec(e)
            E0 += 0.5 * du @ K @ du
            i, j = e[0], e[1]
            g[2 * i:2 * i + 2] += -K @ du
            g[2 * j:2 * j + 2] += K @ du
        w = -Apinv @ g
        return (E0 + g @ w + 0.5 * w @ A @ w) / Acell

    d = 1e-3
    C11 = 2 * W(np.array([[d, 0], [0, 0]])) / d**2
    C12 = W(np.array([[d, 0], [0, d]])) / d**2 - C11
    C66 = 2 * W(np.array([[0, d / 2], [d / 2, 0]])) / d**2
    return C11, C12, C66


def aniso(C):                      # zero iff elastically isotropic (2D)
    C11, C12, C66 = C
    return (C11 - C12 - 2 * C66) / max(abs(C11), 1e-30)


# ---- VALIDATION CONTROLS (known answers) ----------------------------------
# (A) square Bravais, central springs only -> C66=0, strongly ANISOTROPIC
sq = homogenise(np.array([[0., 0.]]), [(0, 0, 1, 0), (0, 0, 0, 1)],
                np.array([1., 0.]), np.array([0., 1.]), 1.0, 0.0)
print(f"control SQUARE (central): C11,C12,C66={tuple(round(x,4) for x in sq)}  "
      f"aniso={aniso(sq):+.3f}")
check("control: square lattice is anisotropic (code CAN detect anisotropy)",
      abs(aniso(sq)) > 0.5 and abs(sq[2]) < 1e-9)

# (B) triangular Bravais, central springs -> known ISOTROPIC
tri = homogenise(np.array([[0., 0.]]),
                 [(0, 0, 1, 0), (0, 0, 0, 1), (0, 0, -1, 1)],
                 np.array([1., 0.]), np.array([0.5, np.sqrt(3) / 2]), 1.0, 0.0)
print(f"control TRIANGULAR (central): C11,C12,C66={tuple(round(x,4) for x in tri)}  "
      f"aniso={aniso(tri):+.2e}")
check("control: triangular lattice is isotropic (code confirms known isotropy)",
      abs(aniso(tri)) < 1e-6)

# ===========================================================================
# 3. the 4.8.8 gauge web: is the isotropy a property at ALL k_t/k_l?
# ===========================================================================
print("\n4.8.8 anisotropy (C11-C12-2C66)/C11 vs rho=k_t/k_l:")
rhos = [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
anis = []
for rho in rhos:
    C = homogenise(R, edges, A1, A2, 1.0, rho)
    anis.append(aniso(C))
    print(f"  rho={rho:.2f}: C11={C[0]:.4f} C12={C[1]:.4f} C66={C[2]:.4f}  "
          f"aniso={aniso(C):+.2e}")

iso_all = all(abs(a) < 1e-8 for a in anis)
check("4.8.8 is elastically ISOTROPIC for ALL k_t/k_l (confirms line-458 theorem)",
      iso_all)

print("\n--- VERDICT ---")
if iso_all:
    print("The 4.8.8 web is elastically isotropic at EVERY k_t/k_l. This independently")
    print("CONFIRMS the framework's Lorentz-isotropy claim (ANCHOR L458) at the elastic-")
    print("homogenisation level -- a genuine positive corroboration. BUT it also means")
    print("isotropy is GEOMETRIC/automatic, NOT a constraint: it does NOT fix k_t/k_l.")
    print("So the lever this script set out to use FAILS: the Cauchy-Born route does not")
    print("become parameter-free via isotropy. k_t/k_l remains a genuine free input")
    print("(as ew_cauchy_born_check.py argued) -- now CONSTRUCTIVELY confirmed, with the")
    print("homogenisation actually built. M_H/M_Z stays blocked on the free k_t/k_l")
    print("(plus the gauge coupling g and the §16.2 absolute scale).")
else:
    print("4.8.8 anisotropy depends on k_t/k_l; an isotropy root would fix it -- revisit.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
