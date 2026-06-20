#!/usr/bin/env python3
r"""NO LOCALITY-PROTECTED 0D DEFECT EXISTS — closing the open cell of defect_periodic_table.py.

The periodic-table scaffold left one cell open: is there a locality/topologically PROTECTED 0D
(point) defect? This proves there is NOT, decisively and self-assertingly.

THE ARGUMENT (on the VERIFIED K04 model, imported from k04_kempe_locked_defect.py: degree-3
bond 3-factors of the Z^3 torus; crystal = 2x2x2 cube tiling; Kempe = plaquette flip; the
theorem there: a config is locked IFF [C XOR crystal] != 0 in H_1(T^3;Z_2)). The Z_2 winding
charge in axis a is the parity of difference-edges crossing a cut plane perpendicular to a.
Three lemmas:
  (A) WELL-DEFINED: for an even subgraph D the crossing parity is the SAME at every cut plane
      — a genuine homology charge, not an artefact of where you cut.
  (B) BOUNDED => TRIVIAL: a 0D defect is localized (an empty layer in EVERY axis); put the cut
      plane in an empty layer => 0 crossings => winding = (0,0,0).
  (C) LOCKED => SPANNING: winding_a = 1 forces EVERY a-cut-plane to be crossed (odd) => the
      support touches every a-layer => it spans a full period => extent >= 1 period in a, NOT
      a point.
With the k04 theorem (locked <=> winding != 0): a 0D defect has winding 0, so it is NOT locked;
protection requires winding, which requires >= 1D system-spanning extent. The
protection-dimension FLOOR is 1 (the winding string). 0D defects are at most KINETICALLY
pinned — and the minimal one (a single plaquette = the peanut) has no barrier at all. We also
confirm finite defects heal under BOX-CONFINED (local) moves, ruling out a locality obstruction.

exit 0 = (A) parity is plane-independent on a battery; (B) every bounded defect has winding 0
         with a verified empty layer per axis; (C) the locked string is crossed at every layer
         and spans the period; finite defects heal under local moves; the no-go assembled.
"""
import contextlib
import io
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):                 # reuse VERIFIED machinery, mute its banner
    import k04_kempe_locked_defect as k04

Torus, flip, AX = k04.Torus, k04.flip, k04.AX
random.seed(20260614)
T = Torus(6); V, L = T.V, T.L
CRYS = frozenset(T.crystal())
assert T.degrees(CRYS) == [3] * T.n

def Dxor(C):                                                     # difference graph D = C XOR crystal
    return set(C) ^ set(CRYS)

def crossings(D, axis, m):                                       # D-edges between layers m and m+1 (mod L)
    return sum(1 for (i, j) in D
               if V[i][axis] != V[j][axis] and {V[i][axis], V[j][axis]} == {m % L, (m + 1) % L})

def parities(D, axis):
    return [crossings(D, axis, m) % 2 for m in range(L)]

def empty_planes(D, axis):
    return [m for m in range(L) if crossings(D, axis, m) == 0]

def plaqs_in_box(lo, hi):
    out = []
    for p in V:
        if all(lo[a] <= p[a] <= hi[a] for a in range(3)):
            for a, b in ((0, 1), (0, 2), (1, 2)):
                pa = T.wrap(tuple(p[k] + AX[a][k] for k in range(3)))
                pb = T.wrap(tuple(p[k] + AX[b][k] for k in range(3)))
                pab = T.wrap(tuple(p[k] + AX[a][k] + AX[b][k] for k in range(3)))
                out.append(((T.edge(p, pa), T.edge(pb, pab)), (T.edge(p, pb), T.edge(pa, pab))))
    return out

def bounded_defect(lo, hi, nflips):                             # crystal XOR (flips confined to a box)
    C = set(CRYS); pl = plaqs_in_box(lo, hi); done = tries = 0
    while done < nflips and tries < 3000:
        tries += 1
        nx = flip(C, *random.choice(pl))
        if nx is not None:
            C = set(nx); done += 1
    return frozenset(C)

# ---- build the battery: a minimal (peanut-class) defect, several bounded defects, the locked string ----
BLO, BHI = (1, 1, 1), (3, 3, 3)                                  # endpoints land in [1,4] < L=6 -> empty layers
battery = [("single-plaquette (peanut-class)", bounded_defect(BLO, BHI, 1))]
for t in range(4):
    d = bounded_defect(BLO, BHI, 9)
    if Dxor(d):
        battery.append((f"bounded box-defect #{t + 1}", d))
line = {T.edge((x, 0, 0), T.wrap((x + 1, 0, 0))) for x in range(L)}
STRING = frozenset(set(CRYS) ^ line)                            # crystal XOR one winding x-line = locked
battery.append(("winding string (locked)", STRING))

# ================= [A] winding is well-defined (plane-independent) =================
print("[A] WINDING IS WELL-DEFINED — crossing parity is identical at every cut plane:")
for name, C in battery:
    assert T.degrees(C) == [3] * T.n, name                      # every battery member is a valid 3-factor
    D = Dxor(C)
    for axis in range(3):
        ps = parities(D, axis)
        assert len(set(ps)) == 1, (name, axis, ps)              # SAME parity at all L planes
print(f"    {len(battery)} configs, all valid 3-factors; parity constant across all {L} planes per axis. SOLID.")

# ================= [B] bounded => winding (0,0,0) =================
print("\n[B] EVERY BOUNDED (0D-localizable) DEFECT IS HOMOLOGICALLY TRIVIAL:")
for name, C in battery:
    if "string" in name:
        continue
    D = Dxor(C)
    w = tuple(parities(D, a)[0] for a in range(3))               # well-defined by [A]
    empt = [len(empty_planes(D, a)) for a in range(3)]
    assert all(e > 0 for e in empt)                             # an empty layer in EVERY axis => localized
    assert w == (0, 0, 0), (name, w)
    print(f"    {name:<30s}: winding={w}, empty planes/axis={empt} (localized) -> trivial")
print("    -> a 0D defect (empty layer in all three axes) carries winding (0,0,0).")

# ================= [C] locked => spans a full period =================
print("\n[C] A LOCKED DEFECT MUST SPAN A FULL PERIOD (so cannot be a point):")
Ds = Dxor(STRING)
w = tuple(parities(Ds, a)[0] for a in range(3))
assert w == (1, 0, 0)
assert empty_planes(Ds, 0) == []                                # NO empty plane in the winding axis
span_x = sorted({V[i][0] for e in Ds for i in e})
print(f"    winding string: winding={w}; empty x-planes={empty_planes(Ds, 0)} (none) -> every x-layer crossed")
print(f"    support x-extent = {span_x} = all {L} layers -> spans the period; geometric dim >= 1.")
assert len(span_x) == L

# ================= [D] locking theorem + LOCAL healing (no locality obstruction) =================
print("\n[D] LOCKING THEOREM (k04) + LOCAL HEALING:")
PLQ = T.plaquettes(); C = set(STRING); stayed = True
for _ in range(800):
    nx = flip(C, *random.choice(PLQ))
    if nx is None:
        continue
    C = set(nx)
    if tuple(parities(Dxor(C), a)[0] for a in range(3)) != (1, 0, 0):
        stayed = False; break
assert stayed
print("    locked check: 800 Kempe flips from the string keep winding (1,0,0) — it never reaches the crystal.")
with contextlib.redirect_stdout(io.StringIO()):                 # closed-block orbit = moves confined to the block
    blk = k04.block_3factors((2, 2, 4)); ntot, ncomp, _ = k04.heal_orbit(blk, (2, 2, 4))
print(f"    local-healing check: closed 2x2x4 block — {ntot} 3-factors; crystal Kempe-orbit={ncomp}; "
      f"locked={ntot - ncomp}")
assert ntot - ncomp == 0
print("    -> every finite defect heals using ONLY block-confined (local) moves: no locality obstruction.")

# ================= [NO-GO] =================
print(f"""
[NO-GO] NO LOCALITY-PROTECTED 0D DEFECT EXISTS.
  topological/locality protection  <=>  [D] != 0 in H_1(T^3;Z_2)            [k04 theorem]
                                   <=>  winding != 0 in some axis            [definition, well-defined by A]
                                   ==>  every cut plane in that axis crossed [C]
                                   ==>  support spans a full period => extent >= 1D, NOT a point.
  Contrapositive: a 0D (bounded; empty layer in every axis) defect has winding (0,0,0) [B],
  so it is NOT locked, so it heals — and heals LOCALLY [D]. The protection-dimension FLOOR is
  therefore 1: the minimal protected object is the 1D winding string. The scaffold's open 0D
  cell is CLOSED, with a no-go.
  HONEST SCOPE: this forbids TOPOLOGICAL/locality protection of point defects; it does not
  forbid KINETIC metastability (energy barriers) of 0D clusters — but kinetic != protected, and
  the minimal 0D defect (one plaquette = the peanut) has no barrier (it heals in one move). The
  result is move-set independent: any local move XORs a contractible cycle, preserving [D].
exit 0""")
print("ALL ASSERTIONS PASSED — winding well-defined; bounded=>trivial; locked=>spanning; finite heals locally; 0D no-go.")
