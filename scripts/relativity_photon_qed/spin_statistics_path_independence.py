#!/usr/bin/env python3
"""
spin_statistics_path_independence.py
====================================
Tests ANCHOR §15 item 139 promotion-prerequisite (ii)/(b) [path independence] and
(iii)/(c) [body-diagonal swap paths], described in ANCHOR as "a finite group-
theoretic computation, currently sketched but not exhaustively exhibited."

ANCHOR Part 4 builds the exchange operator's coin Wilson line as the ordered
product of bridge-rotation lifts along the swap path, and exhibits, for ONE
face-diagonal planar path pair, P_12^2 = C4z^4 = exp(-i pi sigma_z) = -1. The
open claim (ii) is that ANY disjoint path-pair in the pi_1=Z2 generator class
gives the same -1.

This script asks, default-to-refute: is the *naive lattice construction*
("ordered product of the per-corner bridge-rotation coin lifts along the path")
actually path-independent -- a topological invariant of the loop -- or is it
geometry/framing-dependent, so that closing (ii)/(iii) needs more than a finite
enumeration?

Method. The per-corner bridge rotation = the rotation aligning consecutive hop
directions (rotation-minimizing frame on the lattice; this is the canonical
reading of "bridge-rotation lift", and it reproduces ANCHOR's C4z^4 for the
planar square). The coin holonomy = ordered product of the SU(2) lifts
exp(-i (theta/2) n.sigma) of those rotations around the closed loop.

numpy only. exit 0 == the MATH HARNESS holds (planar square -> -I reproducing
ANCHOR Part 4; a 2pi coin rotation about any axis = -I). The FINDINGS about
path-(in)dependence are printed and interpreted; the verdict is stated, not
hidden behind a green exit.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)

def Urot(theta, n):
    """SU(2) lift of a spatial rotation by angle theta about unit axis n."""
    n = np.asarray(n, float); n = n / np.linalg.norm(n)
    return np.cos(theta / 2) * I2 - 1j * np.sin(theta / 2) * (n[0]*SX + n[1]*SY + n[2]*SZ)

def align(a, b):
    """Rotation (theta, axis) taking unit vector a to unit vector b; flag if antiparallel (ambiguous)."""
    a = np.asarray(a, float); a /= np.linalg.norm(a)
    b = np.asarray(b, float); b /= np.linalg.norm(b)
    c = float(np.clip(np.dot(a, b), -1, 1))
    ax = np.cross(a, b); nrm = np.linalg.norm(ax)
    if nrm < 1e-12:
        return 0.0, np.array([0.0, 0.0, 1.0]), (c < 0)   # parallel->id ; antiparallel->ambiguous
    return np.arccos(c), ax / nrm, False

def holonomy(verts):
    """Coin Wilson line around the closed polyline `verts` (rotation-minimizing frame)."""
    V = [np.asarray(v, float) for v in verts]; n = len(V)
    dirs = [(V[(i + 1) % n] - V[i]) for i in range(n)]
    dirs = [d / np.linalg.norm(d) for d in dirs]
    U = I2.copy(); total_turn = 0.0; amb = False
    for i in range(n):
        th, ax, a = align(dirs[i - 1], dirs[i]); amb = amb or a; total_turn += th
        U = Urot(th, ax) @ U
    return U, total_turn, amb

def label(U):
    if np.allclose(U, -I2, atol=1e-8): return "-I"
    if np.allclose(U,  I2, atol=1e-8): return "+I"
    return "neither +/-I"

def rot_angle(U):
    """SU(2) rotation angle psi in [0, 2pi]: U = cos(psi/2) I - i sin(psi/2) n.sigma.
    SIGNED via the trace (no abs): +I -> 0, -I -> 2pi, generic -> in between."""
    return 2 * np.arccos(float(np.clip(np.trace(U).real / 2, -1, 1)))

fails = []
def assert_(cond, msg):
    print(f"    [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond: fails.append(msg)

print("=" * 74)
print(" ITEM 139 PREREQ (ii)/(iii): is the naive lattice exchange Wilson line")
print(" path-independent, or geometry/framing-dependent?   (default-to-refute)")
print("=" * 74)

# --------------------------------------------------------------- HARNESS
print("\n[H] Harness -- reproduce ANCHOR Part 4 + the 2pi-rotation fact")
square_xy = [(0,0,0), (1,0,0), (1,1,0), (0,1,0)]
U, turn, _ = holonomy(square_xy)
print(f"    planar unit square (xy): Wilson line = {label(U)}, total turning = {turn/np.pi:.3f} pi")
assert_(label(U) == "-I", "planar square -> -I  (reproduces ANCHOR Part 4 C4z^4 = exp(-i pi sz))")
assert_(abs(turn - 2*np.pi) < 1e-9, "planar simple loop total turning = 2 pi (turning-number theorem)")

for ax, nm in [((1,0,0),"x"), ((0,1,0),"y"), ((0,0,1),"z"), ((1,1,1),"[111]"), ((2,-1,3),"random")]:
    assert_(np.allclose(Urot(2*np.pi, ax), -I2, atol=1e-9),
            f"2pi coin rotation about {nm:7s} = -I  (exp(-i pi n.sigma) = -I for ALL axes n)")

# ---------------------------------------------- (ii) within the planar class
print("\n[1] Prereq (ii), planar class: different shapes / planes / windings")
loops = {
    "square in yz":        [(0,0,0),(0,1,0),(0,1,1),(0,0,1)],
    "square in zx":        [(0,0,0),(0,0,1),(1,0,1),(1,0,0)],
    "2x1 rectangle (xy)":  [(0,0,0),(2,0,0),(2,1,0),(0,1,0)],
    "octagon-ish (xy)":    [(2,0,0),(1.4,1.4,0),(0,2,0),(-1.4,1.4,0),(-2,0,0),(-1.4,-1.4,0),(0,-2,0),(1.4,-1.4,0)],
}
for nm, lp in loops.items():
    U, turn, amb = holonomy(lp)
    print(f"    {nm:22s}: {label(U):12s}  turning = {turn/np.pi:.3f} pi")
all_planar_minusI = all(label(holonomy(lp)[0]) == "-I" for lp in loops.values())
print(f"    -> all planar winding-1 loops give -I (any plane, any shape)?  {all_planar_minusI}")
print("       REASON: a simple planar loop has total turning 2pi, and all corner-axes")
print("       are the SAME (the plane normal), so the lifts commute and sum to a 2pi")
print("       rotation -> -I. This is the abelian case ANCHOR's construction lives in.")

# double winding (figure-tracing the square twice) -> +I ; this is informative:
sq2 = square_xy + square_xy
U2, turn2, _ = holonomy(sq2)
print(f"    square traced TWICE (winding 2): {label(U2)}  turning = {turn2/np.pi:.3f} pi  (4pi -> +I)")

# ----------------------------------------- the encircling-vs-planarity probe
print("\n[2] Does the naive Wilson line detect ENCIRCLING the coincidence locus,")
print("    or merely PLANARITY?  (a loop that does NOT enclose a chosen point)")
far_loop = [(5,5,0),(6,5,0),(6,6,0),(5,6,0)]   # small square far from origin, encloses nothing near 0
Uf, tf, _ = holonomy(far_loop)
print(f"    small square far from origin: {label(Uf)}  turning = {tf/np.pi:.3f} pi")
print("    -> Still -I. The naive frame-self-rotation Wilson line gives -I for EVERY")
print("       simple planar loop, whether or not it encircles the x1=x2 diagonal. So the")
print("       -I in Part 4 tracks the frame's 2pi self-rotation, NOT the loop's homotopy")
print("       class around the coincidence locus. The topology is not being detected.")

# ------------------------------------------- (iii) non-planar / body-diagonal
print("\n[3] Prereq (iii): non-planar loops -- symmetric vs asymmetric")
nonplanar = {
    "cube Petrie hexagon (sym)":  [(1,0,0),(1,1,0),(0,1,0),(0,1,1),(0,0,1),(1,0,1)],
    "cube-edge loop (sym)":       [(0,0,0),(2,0,0),(2,2,0),(2,2,2),(0,2,2),(0,0,2)],
    "asym 3D loop A":             [(0,0,0),(3,0,0),(3,2,0),(3,2,1),(0,1,1),(0,0,1)],
    "asym 3D loop B":             [(0,0,0),(4,0,0),(4,1,0),(4,1,3),(1,1,3),(1,0,1)],
    "asym 3D loop C":             [(0,0,0),(2,0,0),(2,3,0),(2,3,2),(0,2,2),(0,1,1)],
}
labels = {}
for nm, lp in nonplanar.items():
    U, turn, amb = holonomy(lp)
    labels[nm] = label(U)
    print(f"    {nm:26s}: {labels[nm]:14s}  holonomy rot psi = {rot_angle(U)/np.pi:.3f} pi"
          f"   (turning sum = {turn/np.pi:.3f} pi)")
n_minusI = sum(1 for v in labels.values() if v == "-I")
n_neither = sum(1 for v in labels.values() if v == "neither +/-I")
print(f"    -> of {len(labels)} non-planar loops: {n_minusI} give -I, {n_neither} give neither +/-I.")
print("       The closed loop fixes its initial tangent, so the holonomy is a rotation")
print("       ABOUT that tangent by an angle psi = (geometric/solid-angle term). Planar &")
print("       symmetric loops have psi = 2pi -> -I; ASYMMETRIC non-planar loops have")
print("       psi != 2pi -> NEITHER +/-I. So the naive Wilson line is GEOMETRY-DEPENDENT,")
print("       not even valued in {+/-I}, let alone a topological Z2 invariant.")

print("\n" + "=" * 74)
print(" VERDICT")
print("=" * 74)
print(" - HARNESS verified: ANCHOR Part 4's C4z^4 = -I is arithmetically correct, and")
print("   a 2pi coin rotation about ANY axis = -I (so the -I is axis/plane-independent).")
print(" - Prereq (ii) holds ONLY within the planar/coaxial class (abelian: equal corner-")
print("   axes, turning 2pi). That is the class ANCHOR's single exhibited path lives in.")
print(" - REFUTATION of the easy reading: the naive 'product of bridge-rotation lifts'")
print("   is NOT path-independent. The closed loop fixes its initial tangent, so the")
print("   holonomy is a rotation about that tangent by a GEOMETRIC angle psi:")
print("     (a) planar loops: psi = 2pi -> -I, but for EVERY simple planar loop incl.")
print("         ones that do NOT encircle x1=x2 -- so it tracks the frame's 2pi self-")
print("         rotation (tracing parity), not the pi_1(C)=Z2 class (which must be +I")
print("         for a contractible/non-encircling exchange).")
print("     (b) non-planar loops: psi = 2pi (-> -I) only for symmetric ones; asymmetric")
print("         loops give psi != 2pi -> NEITHER +/-I. So it is not even {+/-I}-valued,")
print("         let alone a Z2 invariant.")
print("   Hence (ii)/(iii) are NOT closeable by a finite naive enumeration; they need a")
print("   TOPOLOGICAL FRAMING making the coin holonomy a function of the Z2 class alone")
print("   -- exactly the Finkelstein-Rubinstein ingredient ANCHOR (b)/(c) flags unworked.")
print(" - CONSEQUENCE for item 139: stays 'derivation-with-prerequisites'. The promotion")
print("   criterion (ii)/(iii) is harder than 'a finite group-theoretic computation';")
print("   this computation should sharpen that caveat rather than close it.")
print("\n exit 0 == math harness holds (the refutation is a finding, not a crash).")
import sys
sys.exit(1 if fails else 0)
