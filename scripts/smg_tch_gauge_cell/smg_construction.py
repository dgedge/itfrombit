#!/usr/bin/env python3
"""
(e+) ADD and JUSTIFY a symmetric multi-flavour 4-fermion interaction beyond §5's diagonal
stabilizers and the bilinear Feshbach, and PROVE it lands in the symmetric gapped phase.

smg_phase_decisive.py showed §5's CURRENT terms are exhausted by {diagonal commuting Z_iZ_j
-> projection} U {bilinear Feshbach -> breaking}; the symmetric SMG 'third type' was absent.
This script constructs that third type -- and it is NOT ad hoc:

  JUSTIFICATION (ANCHOR §15 item 13, verbatim): "The mass mechanism (§5) uses only Z-type
  stabilizers {Z_iZ_j}. The gauge sector (§7) uses ... X-type stabilizer operations. Under
  the strict self-duality of the [8,4,4] extended Hamming code, a CSS-style combination would
  unify the matter-cell mass code and the gauge-cell electrodynamics code into a single
  quantum-error-correcting structure."  => the missing symmetric off-diagonal interaction is
  the X-type stabilizers of the self-dual [8,4,4] CSS code (the framework's OWN gauge sector).
  Completing the Z-only mass code to the full self-dual CSS code IS item 13's unification.

CLAIM PROVED HERE (exactly, on the cell):
  the full self-dual [8,4,4] CSS Hamiltonian H = -sum(Z-stab) - sum(X-stab) has a UNIQUE,
  GAPPED, ENTANGLED ground state with ZERO single-fermion bilinear order -- symmetric mass
  generation on the cell, exactly. The X-half is exactly what turns §5's Z-only product/
  projection ground space into a unique entangled symmetric state. The construction REQUIRES
  self-duality, which is the structure tied to the anomaly-free 16 (smg_code_projection.py).

HONEST SCOPE (frontier, NOT claimed): applying the X-SMG interaction to the MIRROR sector
only in the full 3+1D spatial lattice (keeping the physical fermion chiral + correctly
gauge-coupled), with DYNAMICAL gauge fields, and the same-chirality residual -- these are the
open 3+1D chiral-lattice-gauge frontier. What is proved is the LOCAL symmetric-gapping
mechanism and that the interaction is the framework's own CSS X-half.

numpy; self-asserting.
"""
import numpy as np, itertools
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]],complex); Z=np.diag([1,-1]).astype(complex)
def pauli(d):
    M=np.array([[1]],complex)
    for k in range(8): M=np.kron(M, d.get(k,I2))
    return M
def Xc(c): return pauli({k:X for k in range(8) if c[k]})
def Zc(c): return pauli({k:Z for k in range(8) if c[k]})
def comm(A,B): return float(np.linalg.norm(A@B-B@A))

# ----------------------------------------------------------------------------------
Hd("(A) The self-dual [8,4,4] code and its CSS stabilizers (Z = mass §5, X = gauge §7)")
G=[[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]   # RM(1,3), self-dual
def span(rows):
    out=[]
    for bits in itertools.product([0,1],repeat=len(rows)):
        c=[0]*8
        for b,r in zip(bits,rows):
            if b: c=[x^y for x,y in zip(c,r)]
        out.append(tuple(c))
    return out
C=span(G)
def dot(a,b): return sum(x&y for x,y in zip(a,b))%2
print(f"  |C| = {len(set(C))}  (= 2^4); min nonzero weight = {min(sum(c) for c in C if any(c))} (distance 4)")
assert len(set(C))==16 and min(sum(c) for c in C if any(c))==4
print(f"  self-dual: all pairwise dot products even = {all(dot(a,b)==0 for a in C for b in C)}")
assert all(dot(a,b)==0 for a in C for b in C)
Zs=[Zc(g) for g in G]; Xs=[Xc(g) for g in G]; stab=Zs+Xs
print(f"  4 Z-stabilizers (§5 mass) + 4 X-stabilizers (§7 gauge); CSS commute = "
      f"{all(comm(a,b)<1e-9 for a in stab for b in stab)}")
assert all(comm(a,b)<1e-9 for a in stab for b in stab)
print("  => Z-half is §5's mass code (item 13); X-half is §7's gauge sector (item 13).")
print("     Their self-dual CSS union is item 13's open unification = the SMG interaction.")

# ----------------------------------------------------------------------------------
Hd("(B) PROOF: full CSS gaps to a UNIQUE state; Z-only is a 16-fold (projection) degeneracy")
def proj(gens):
    P=np.eye(256,dtype=complex)
    for S in gens: P=P@((np.eye(256,dtype=complex)+S)/2)
    return P
rank_full=round(np.real(np.trace(proj(stab))))
rank_Zonly=round(np.real(np.trace(proj(Zs))))
print(f"  rank of full-CSS ground projector  : {rank_full}   ([[8,0,4]] -> UNIQUE state)")
print(f"  rank of Z-only  ground projector   : {rank_Zonly}  (= 2^4 logical codespace)")
assert rank_full==1 and rank_Zonly==16
Hcss=-sum(stab)
ev=np.round(np.real(np.linalg.eigvalsh(Hcss)),6)
gap=sorted(set(ev))[1]-ev.min()
print(f"  full-CSS spectrum gap (E1-E0)      : {gap:.1f}   (commuting stabilizers -> exact gap)")
assert abs(gap-2.0)<1e-9
print("  => §5's Z-only mechanism leaves a 16-fold degeneracy (the chiral logical space, where")
print("     the fermion lives -- 'projection'); ADDING the X-half gaps it to ONE state. The")
print("     X (off-diagonal) stabilizers are exactly what create a dynamical gap, not a projection.")

# ----------------------------------------------------------------------------------
Hd("(C) PROOF: the unique gapped state is ENTANGLED and has ZERO bilinear (mass) order")
w,V=np.linalg.eigh(Hcss); psi=V[:,0]
# entanglement across qubits {0,1,2,3} | {4,5,6,7}
M=psi.reshape(16,16); s=np.linalg.svd(M,compute_uv=False); s=s[s>1e-12]
S_ent=float(-np.sum(s**2*np.log2(s**2)))
print(f"  entanglement entropy (4|4 cut) of the unique state: {S_ent:.2f} bits  (>0 => ENTANGLED,")
print(f"     NOT a product/computational-basis state -- unlike §5's Z-only product ground states)")
assert S_ent>1.0
# Z-only ground space admits product (computational-basis) states -> 'projection' picture
cw_state=np.zeros(256,complex); cw_state[int("".join(map(str,C[1])),2)]=1  # a codeword basis vector
Mc=cw_state.reshape(16,16); sc=np.linalg.svd(Mc,compute_uv=False)
print(f"  (a Z-only codeword basis state has entanglement {float(-np.sum((sc[sc>1e-12]**2)*np.log2(sc[sc>1e-12]**2))):.2f} "
      f"bits -> product, the projection picture)")
# zero bilinear order: any single-fermion mass bilinear -> a Pauli; stabilizer state => 0 unless in group.
# the would-be mirror Dirac mass is the chirality(chi=qubit6)-flip; check X_chi, Y_chi and all 1-qubit Paulis.
ords={}
for k in range(8):
    for nm,Op in (("X",X),("Y",Y),("Z",Z)):
        ords[f"{nm}_{k}"]=abs(np.vdot(psi,pauli({k:Op})@psi))
maxord=max(ords.values())
print(f"  max |<single-qubit Pauli>| over the unique state: {maxord:.1e}")
print(f"     in particular the chirality-flip mass <X_chi>={ords['X_6']:.1e}, <Y_chi>={ords['Y_6']:.1e}  (= 0)")
assert maxord<1e-9
# representative 2-qubit chirality-flip pairing mass (X on chi & W) also zero unless a stabilizer
pair=abs(np.vdot(psi,pauli({6:X,7:X})@psi))
print(f"  representative 2-fermion chirality pairing <X_chi X_W> = {pair:.1e}")
print("  => NO bilinear (and no checked low-order) mass condensate: the gap is SYMMETRIC, not")
print("     spontaneous symmetry breaking. (Stabilizer-state Paulis are 0 off the stabilizer group;")
print("     distance-4 code => no weight-<=2 stabilizers => all 1- and 2-local masses vanish.)")
nonzero_lowweight = [c for c in C if 0<sum(c)<=2]
assert len(nonzero_lowweight)==0

# ----------------------------------------------------------------------------------
Hd("(D) Why this is SMG and not the forbidden projection (reconciles smg_phase_decisive (A))")
print("""  smg_phase_decisive (A): §5's Z_iZ_j are DIAGONAL+commuting -> product ground state ->
  projection -> N-N-forbidden. Correct -- for the Z-only (diagonal) half.
  Here the X-stabilizers are OFF-DIAGONAL (products of X). The full CSS stabilizer group is
  still commuting, but its unique ground state is ENTANGLED (C, {S_ent:.1f} bits), NOT a product
  computational-basis state. An entangled symmetric short-range-entangled gapped state is
  exactly the SMG evasion of Nielsen-Ninomiya (interacting, not free/projective). So:
    diagonal commuting  -> product  -> projection (forbidden)   [§5 as-is]
    full-CSS commuting  -> ENTANGLED -> symmetric gap (SMG)      [+ X-half, item 13]
  The distinction is the entanglement of the ground state, created by the off-diagonal X-half.""".replace("{S_ent:.1f}", f"{S_ent:.1f}"))

# ----------------------------------------------------------------------------------
Hd("VERDICT — the symmetric interaction EXISTS, is the framework's own X-half, gap PROVED locally")
print("""CONSTRUCTED + JUSTIFIED + (locally) PROVED:
  * The 'third type' interaction smg_phase_decisive said was absent EXISTS: the X-type
    stabilizers of the self-dual [8,4,4] CSS code -- the framework's OWN gauge-sector
    stabilizers (§7), whose unification with the §5 Z mass-code is ANCHOR item 13. Not ad hoc.
  * PROVED exactly on the cell: full-CSS H = -sum Z - sum X has a UNIQUE (rank-1), GAPPED
    (gap 2, commuting stabilizers), ENTANGLED (~{E} bits) ground state with ZERO 1- and 2-local
    bilinear order (distance-4 code) -> a SYMMETRIC gapped state, not spontaneous breaking,
    not the forbidden product/projection. The X (off-diagonal) half is precisely what converts
    §5's 16-fold projection degeneracy into a unique symmetric entangled gap.
  * REQUIRES self-duality of [8,4,4] = the structure tied to the anomaly-free 16 -> uniquely-TCH.

So the earlier 'no SMG' was a verdict on §5's Z-ONLY truncation. COMPLETING the mass code to
the full self-dual CSS code (item 13's unification) supplies a symmetric mirror-gapping
interaction whose symmetric gapped phase is LOCALLY EXACT.

FRONTIER (honestly NOT proved -- the 3+1D chiral-lattice-gauge open problem):
  1. apply the X-SMG interaction to the MIRROR sector ONLY in the full spatial lattice, gapping
     the mirror while keeping the physical fermion chiral, gapless, and correctly gauge-coupled;
  2. with DYNAMICAL (gauged) SM fields, show the X-stabilizers remain gauge-invariant (item 13's
     gauge-sector identification) so the gap does not break the gauged symmetry;
  3. the same-chirality residual (walk_kernel_overlap's 4 species), untouched by this.

NET for item 93/143: the SMG route moves from 'dynamically unrealised' to 'symmetric gapping
interaction CONSTRUCTED from the framework's own CSS code (item 13) and the symmetric gapped
phase PROVED locally/exactly; global 3+1D gauge-coupled realisation is the remaining frontier.'
This is a genuine step toward the ultralocal tier -- the mechanism and operator now exist and
work locally -- without yet reaching continuum-Locked (the global step is unproved, and is the
field-wide open problem of a complete chiral lattice gauge theory).""".replace("{E}", f"{S_ent:.1f}"))
print("\nALL ASSERTS PASSED (construction + local symmetric-gap proof are exact; the 3+1D")
print("gauge-coupled mirror-only realisation is explicitly NOT claimed).")
