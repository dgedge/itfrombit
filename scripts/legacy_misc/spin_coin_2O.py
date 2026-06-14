#!/usr/bin/env python3
"""
Explicit construction of the C^2 E_{1/2} spin-coin for ANCHOR §15 item 138
(spin-rotation locking / 4pi-periodicity as a lattice operation).

Item 138 states the load-bearing gap: the framework has an internal spinor space
and a lattice rotation group (O, double cover 2O) but "NO demonstrated
homomorphism between them ... the 4pi signature is ASSERTED, not constructed."
This script constructs it (the discrete part) and verifies, exactly:

 (A) 2O subset SU(2): 48 elements, closes as a group, double-covers O (24 SO(3)
     rotations), central element -1 maps to the SO(3) identity = the 2pi rotation.
 (B) E_{1/2} coin = the C^2 defining rep: irreducible (Burnside), and R(2pi) = -1
     on it -> spin lives here.
 (C) The 8 body-diagonal directions carry a SINGLE-VALUED permutation rep:
     R(2pi) = +1, decomposing 8 = A1+A2+T1+T2 -> spin CANNOT live in the
     direction orbit (item 138's load-bearing correction).
 (D) The homomorphism g -> U_g (x) P_g (coin E_{1/2} (x) direction-permutation):
     the 2pi rotation acts as (-1 on coin) (x) (+1 on dirs) -> a one-particle
     spinor state is multiplied by -1 by a 2pi lattice rotation. Spin-rotation
     LOCKED (discretely).
 (E) Canonical bridge rotations C3[111], C4z lift to 2O and GENERATE all 48 ->
     the full double cover (incl. the central 2pi = -1) is realized.

Scope (honest): this closes the *spin-rotation* homomorphism at the DISCRETE
double-group level. It does NOT by itself give (i) spin-STATISTICS (item 139
needs the Finkelstein-Rubinstein exchange<->rotation argument), (ii) continuous
SU(2) spin / 4pi for arbitrary rotations (needs the continuum limit, item 102),
or (iii) the anomalous moment a_e (item 138 wrong-Green's-function wall). g=2 is
the tree-level Dirac value the E_{1/2} coin realizes, not new physics.

numpy only. exit 0 == every assertion verified.
"""
import numpy as np, itertools
np.set_printoptions(suppress=True)

# ---- SU(2) from unit quaternion a + b i + c j + d k ----
def su2(a, b, c, d):
    return np.array([[a + 1j*b, c + 1j*d], [-c + 1j*d, a - 1j*b]])
sx = np.array([[0,1],[1,0]]); sy = np.array([[0,-1j],[1j,0]]); sz = np.array([[1,0],[0,-1]])

# ---- (A) build the 48 quaternions of the binary octahedral group 2O ----
quats = []
for i in range(4):                                   # type 1: (+-1,0,0,0) (8)
    for s in (1,-1):
        q=[0,0,0,0]; q[i]=s; quats.append(tuple(q))
for sg in itertools.product((1,-1),repeat=4):        # type 2: (+-1/2)^4 (16)
    quats.append(tuple(0.5*s for s in sg))
r=1/np.sqrt(2)                                        # type 3: two +-1/sqrt2 (24)
for i,j in itertools.combinations(range(4),2):
    for si,sj in itertools.product((1,-1),repeat=2):
        q=[0,0,0,0]; q[i]=si*r; q[j]=sj*r; quats.append(tuple(q))
assert len(quats)==48
G = [su2(*q) for q in quats]
for U in G:                                          # all genuinely in SU(2)
    assert abs(np.linalg.det(U)-1)<1e-9 and np.allclose(U@U.conj().T, np.eye(2))

def key(U): return tuple(np.round(U.flatten(),6))
def closed_group(gens, cap=200):                     # BFS closure
    seen={key(g):g for g in gens}; seen[key(np.eye(2))]=np.eye(2)
    frontier=list(seen.values())
    while frontier:
        nf=[]
        for a in frontier:
            for g in gens:
                p=a@g; k=key(p)
                if k not in seen: seen[k]=p; nf.append(p)
        frontier=nf
        if len(seen)>cap: break
    return list(seen.values())
grp = closed_group(G)
print(f"(A) |2O| (closure of the 48 quaternions) = {len(grp)}")
assert len(grp)==48, "2O must have 48 elements"

# SO(3) image of each element; count distinct -> double cover
def so3(U):
    R=np.zeros((3,3))
    P=[sx,sy,sz]
    for a in range(3):
        for b in range(3):
            R[a,b]=0.5*np.trace(P[a]@U@P[b]@U.conj().T).real
    return R
so3_imgs={tuple(np.round(so3(U).flatten(),5)) for U in G}
print(f"    distinct SO(3) images = {len(so3_imgs)}  (=|O|; 2:1 double cover: {len(G)//len(so3_imgs)})")
assert len(so3_imgs)==24 and len(G)//len(so3_imgs)==2

minus1 = su2(-1,0,0,0)
print(f"    central element -1: SO(3) image = identity? {np.allclose(so3(minus1),np.eye(3))}"
      f" ; acts on coin as -I? {np.allclose(minus1,-np.eye(2))}")
assert np.allclose(so3(minus1),np.eye(3)) and np.allclose(minus1,-np.eye(2))
print("    => the 2pi rotation (identity in O) lifts to -1 in 2O: the spinor sign.")

# ---- (B) E_{1/2} coin = defining C^2 rep; irreducible by Burnside ----
burnside = sum(abs(np.trace(U))**2 for U in G)/len(G)
print(f"\n(B) E_1/2 coin: (1/|G|) sum |chi|^2 = {burnside:.4f}  (==1 -> irreducible)")
assert abs(burnside-1)<1e-9
print(f"    R(2pi) on the coin = -I  ->  SPIN LIVES HERE (2pi -> -1).")

# ---- (C) 8 body-diagonal directions: single-valued permutation rep ----
dirs = np.array(list(itertools.product((1,-1),repeat=3)),float)/np.sqrt(3)   # 8 dirs
def perm_of(U):                                      # permutation of the 8 dirs under so3(U)
    R=so3(U); out=[]
    for d in dirs:
        Rd=R@d; idx=int(np.argmin(np.linalg.norm(dirs-Rd,axis=1))); out.append(idx)
    return tuple(out)
# character of the permutation rep = # fixed directions; decompose over O via Burnside-per-irrep
perm_minus1 = perm_of(minus1)
print(f"\n(C) 8 body-diagonal directions, permutation rep:")
print(f"    perm of -1 (the 2pi rotation) = identity permutation? {perm_minus1==tuple(range(8))}")
assert perm_minus1==tuple(range(8))
print(f"    => R(2pi) = +1 on the direction orbit: SINGLE-VALUED, spin CANNOT live here.")
# Burnside count of irreducible components (how many irreps, with multiplicity^2 sum):
chi = {}
for U in G:
    p=perm_of(U); chi[key(U)]=sum(1 for i in range(8) if p[i]==i)  # fixed points
n_irr = sum(c*c for c in chi.values())/len(G)
print(f"    sum |chi_perm|^2 /|G| = {n_irr:.2f}  (==4 distinct irreps: 8 = 1+1+3+3 = A1+A2+T1+T2)")
assert abs(n_irr-4)<1e-9

# ---- (D) the homomorphism g -> U_g (x) P_g ; 2pi -> (-I_coin) (x) (I_dir) ----
print(f"\n(D) Total one-particle action g -> U_g (coin E_1/2)  (x)  P_g (8 directions):")
print(f"    under a 2pi lattice rotation the coin gets -1 and the directions are fixed,")
print(f"    so a spinor state |psi> (x) |dir> -> -|psi> (x) |dir>.  SPIN-ROTATION LOCKED.")
# spot-check homomorphism on a pair: U_a U_b 's SO(3) = product of SO(3)'s, and perm composes
import random
for _ in range(200):
    A,B=random.choice(G),random.choice(G)
    assert np.allclose(so3(A@B), so3(A)@so3(B), atol=1e-9)   # so3 is a homomorphism
print("    verified: U -> SO(3)(U) is a homomorphism (200 random pairs); the coin rep and")
print("    the direction permutation are both functorial in g -> consistent total action.")

# ---- (E) bridge rotations generate 2O ----
C4z = su2(np.cos(np.pi/4),0,0,np.sin(np.pi/4))   # 90 deg about z
C3  = su2(0.5,0.5,0.5,0.5)                        # 120 deg about [111]
print(f"\n(E) canonical bridge rotations:")
print(f"    C4z (90 about z) in 2O? {any(np.allclose(C4z,U) for U in G)}")
print(f"    C3 (120 about [111]) in 2O? {any(np.allclose(C3,U) for U in G)}")
assert any(np.allclose(C4z,U) for U in G) and any(np.allclose(C3,U) for U in G)
gen = closed_group([C4z,C3])
print(f"    <C4z, C3> generates {len(gen)} elements")
assert len(gen)==48, "bridge rotations must generate all of 2O"
got_minus1 = any(np.allclose(-np.eye(2),U) for U in gen)
print(f"    the generated group contains -1 (the 2pi rotation)? {got_minus1}")
assert got_minus1
print("    => the canonical bridge rotations REALIZE the full 2O double cover, incl. 2pi=-1.")
# 4pi-periodicity as an EXPLICIT lattice operation: four 90-deg z-turns = 2pi = -1; eight = 4pi = +1
C4z_4 = np.linalg.matrix_power(C4z,4); C4z_8 = np.linalg.matrix_power(C4z,8)
print(f"    C4z^4 (=2pi about z) = -I on coin? {np.allclose(C4z_4,-np.eye(2))}"
      f" ;  C4z^8 (=4pi) = +I? {np.allclose(C4z_8,np.eye(2))}")
assert np.allclose(C4z_4,-np.eye(2)) and np.allclose(C4z_8,np.eye(2))
print("    => 4pi-PERIODICITY IS AN EXPLICIT LATTICE OPERATION: four 90-deg bridge turns send")
print("       the spinor coin to -itself; eight return it (the neutron-interferometry signature).")

print("\n"+"="*72); print("RESULT"); print("="*72)
print(" Constructed (exactly, exit 0): the C^2 E_1/2 coin carrying the double-valued")
print(" irrep of 2O; the 8 body-diagonal directions as a single-valued shift space;")
print(" the homomorphism locking a 2pi lattice rotation to -1 on the spinor coin; and")
print(" the bridge rotations generating the full 2O double cover. This CLOSES item")
print(" 138's spin-rotation-locking homomorphism at the discrete double-group level.")
print(" OPEN (not delivered by this construction):")
print("  - spin-STATISTICS (item 139): ANCHOR §15 139 anchors it as a DERIVATION-WITH-")
print("    PREREQUISITES (NOT closed) -- the Finkelstein-Rubinstein exchange<->rotation")
print("    argument is exhibited for ONE face-diagonal path; promotion needs prereqs (ii)")
print("    path-independence + (iii) body-diagonal, which spin_statistics_path_independence.py")
print("    shows require a TOPOLOGICAL FRAMING (the naive lattice Wilson line is geometry-")
print("    dependent, not even +/-I-valued for asymmetric loops), not a finite enumeration.")
print("    This item-138 coin is 139's prerequisite, not 139 itself.")
print("  - continuous SU(2) spin / 4pi for ARBITRARY rotations: needs the continuum")
print("    limit (the lattice supplies only the finite double group 2O; cf item 102).")
print("  - anomalous moment a_e: item-138 wrong-Green's-function wall stands; g=2 here")
print("    is the tree-level Dirac value the E_1/2 coin realizes, not new physics.")
print("  - freeing chi to pure chirality: this coin makes spin independent of chi, but")
print("    global consistency (R2 W=chi weak coupling, gamma^5 = sigma_y^chi) is unchecked.")
