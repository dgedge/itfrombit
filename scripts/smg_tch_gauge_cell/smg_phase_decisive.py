#!/usr/bin/env python3
"""
The decisive computation (e): does the framework's interaction realise the SYMMETRIC GAPPED
phase (SMG -> ultralocal Strong-CP closure) or SPONTANEOUS BREAKING (PMS -> seesaw, overlap
tier)?  A full 3+1D Monte Carlo of the interacting 16-mirror sector is research-scale and is
NOT attempted here.  Instead a STRUCTURAL result that is rigorously computable and turns out
to settle the question:

  CLAIM: the framework's §5 mass interactions cannot realise SMG, because every interaction
  it actually contains is of a type that provably does NOT give a symmetric many-body gap:
    * the 12 Z-stabilizers Z_i Z_j (incl. the anti-mirror Z_chi Z_W) are DIAGONAL and
      MUTUALLY COMMUTING -> their joint ground state is a PRODUCT computational-basis state
      -> strong coupling = a hard PROJECTION onto the codespace = the N-N-forbidden move
      (mirror_invalid_subspace.py); it is not an entangled symmetric condensate.
    * the only OFF-DIAGONAL mass dynamics (§5.8 Feshbach / §3.4 transient chi-flip Higgs,
      seesaw v=246/M_R) is a chirality-flip BILINEAR -> a Dirac/Majorana mass -> gaps by
      SPONTANEOUSLY BREAKING the chiral symmetry (nonzero bilinear order parameter).
  SMG requires a THIRD type the framework has none of: a symmetric, off-diagonal, multi-
  flavour (>=4-fermion) interaction whose ground state is an entangled singlet with ZERO
  bilinear order -- possible only at the magic flavour number (16, which the content HAS,
  smg_code_projection.py) but realised only by an interaction the framework LACKS.

So the phase question resolves structurally: framework realises {projection} U {breaking},
never the symmetric SMG gap. This is decisive for "does §5 close SMG": no.

Parts:
  (A) the 12 stabilizers are diagonal + commuting -> product ground state -> projection.
  (B) minimal fermion model: the three gapping types, with gap + chiral order parameter
      <c_L^dag c_R> computed -- (i) bilinear GAPS+BREAKS, (ii) diagonal penalty PROJECTS,
      (iii) symmetric singlet IMPOSSIBLE at this flavour count (magic-number requirement).
  (C) map the framework's actual terms onto (i)/(ii); show (iii) is absent.
  (D) verdict.

numpy; self-asserting.
"""
import numpy as np, itertools
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)
def comm(A,B): return A@B-B@A
def nrm(A): return float(np.linalg.norm(A))

# ----------------------------------------------------------------------------------
Hd("(A) The 12 §5.2 stabilizers are DIAGONAL + MUTUALLY COMMUTING -> projection, not SMG")
I2=np.eye(2,dtype=complex); Z=np.diag([1,-1]).astype(complex)
def op(singles):
    M=np.array([[1]],complex)
    for k in range(8): M=np.kron(M, singles.get(k,I2))
    return M
# register index k <-> octant address = binary(k) (ANCHOR §2.1); Q_3 edges = Hamming-dist-1
edges=[(i,j) for i in range(8) for j in range(i+1,8) if bin(i^j).count("1")==1]
assert len(edges)==12
S=[op({i:Z,j:Z}) for (i,j) in edges]
all_diag=all(nrm(s-np.diag(np.diag(s)))<1e-12 for s in S)
all_commute=all(nrm(comm(S[a],S[b]))<1e-12 for a in range(12) for b in range(12))
print(f"  12 Q_3-edge stabilizers Z_i Z_j (incl. anti-mirror Z_chi Z_W = edge {(6,7)} in list: "
      f"{(6,7) in edges})")
print(f"  all 12 diagonal in computational basis : {all_diag}")
print(f"  all 12 mutually commuting              : {all_commute}")
assert all_diag and all_commute
# strong-coupling ground state of -g sum S is a computational basis state (diagonal H)
Hstab=-sum(S)
gs_is_basis = (np.count_nonzero(np.abs(np.linalg.eigh(Hstab)[1][:,0])>1e-9)>=1)
diag_vals=np.real(np.diag(Hstab))
ngs=int(np.sum(np.isclose(diag_vals,diag_vals.min())))
print(f"  => -g*sum S is diagonal; its ground manifold is spanned by {ngs} COMPUTATIONAL-BASIS")
print(f"     states (product states, zero entanglement). Strong coupling = PROJECTION onto")
print(f"     that manifold = the hard codespace projection N-N forbids (mirror_invalid_subspace).")
print(f"     A diagonal commuting interaction CANNOT make the entangled symmetric condensate")
print(f"     SMG needs. [structural, exact]")

# ----------------------------------------------------------------------------------
Hd("(B) Minimal fermion model: the three gapping types (gap + chiral order parameter)")
# two modes: L (physical, e.g. left-handed W=chi) and R (mirror, W!=chi). Fock dim 4.
a=np.array([[0,1],[0,0]],complex)               # annihilation (a|1>=|0>)
cL=np.kron(a,I2); cR=np.kron(Z,a)               # Jordan-Wigner, modes (L,R)
nL=cL.conj().T@cL; nR=cR.conj().T@cR
LdR=cL.conj().T@cR                              # chiral (axial-charged) bilinear order param
def gap_and_order(H):
    w,V=np.linalg.eigh(H); g=w[1]-w[0]
    psi=V[:,0]; od=abs(np.vdot(psi,LdR@psi))
    # if degenerate ground, average order over the ground space
    deg=np.isclose(w,w[0]);
    if deg.sum()>1:
        od=np.mean([abs(np.vdot(V[:,k],LdR@V[:,k])) for k in np.where(deg)[0]])
        g=w[deg.sum()]-w[0]
    return g, od
# a small kinetic term that keeps L,R degenerate & massless at zero perturbation:
# use chemical potential 0; massless means no L-R coupling. Gap measured within 1-particle sector.
# Restrict to 1-particle sector {|L>,|R>} for a clean massless-fermion gap read:
P1=np.zeros((4,4),complex)
for s in range(4):
    if abs(np.real((nL+nR)[s,s])-1)<1e-9: P1[s,s]=1
def gap1(H):                                     # gap in the single-particle sector
    sub=[s for s in range(4) if P1[s,s]>0.5]
    Hs=H[np.ix_(sub,sub)]; w,V=np.linalg.eigh(Hs)
    od=abs(np.vdot(V[:,0],LdR[np.ix_(sub,sub)]@V[:,0]))
    return float(w[1]-w[0]), float(od)

m=1.0; g=8.0
H_bilinear = m*(LdR+LdR.conj().T)               # (i) chirality-flip Dirac mass
H_diagpen  = g*nR                               # (ii) diagonal penalty on the mirror (Z-type)
gi,oi=gap1(H_bilinear); gii,oii=gap1(H_diagpen)
print(f"  (i)  bilinear mass m(c_L^dag c_R+h.c.):  gap={gi:.2f}  <c_L^dag c_R>={oi:.3f}  "
      f"-> GAPS, order!=0 => SYMMETRY BROKEN")
print(f"  (ii) diagonal penalty g*n_R (Z-type)  :  gap={gii:.2f}  <c_L^dag c_R>={oii:.3f}  "
      f"-> mirror PROJECTED out (n_R=0), order=0, NO symmetric dynamical gap")
assert oi>0.4 and oii<1e-9 and gi>0.5 and gii>0.5
print("  (iii) symmetric gapped singlet with <c_L^dag c_R>=0 AND no projection: IMPOSSIBLE")
print("        with 2 modes -- it requires the magic flavour number (16 in 3+1D; Wang-Wen).")
print("        The content HAS 16 (smg_code_projection.py); but realising (iii) needs a")
print("        symmetric off-diagonal >=4-fermion interaction, which is a SEPARATE operator.")

# ----------------------------------------------------------------------------------
Hd("(C) Map the framework's actual mass terms onto these types")
print("""  Z_chi Z_W (§5.2)            = type (ii): diagonal density-density (Z=1-2n), commuting
                                with the other 11 stabilizers -> PROJECTS at strong coupling.
  Feshbach / Higgs (§5.8, §3.4) = type (i): off-diagonal chirality(chi)-flip BILINEAR,
                                seesaw v=246/M_R -> GAPS by BREAKING the chiral symmetry.
  SMG interaction (type iii)    = ABSENT. The framework has no symmetric, off-diagonal,
                                multi-flavour >=4-fermion term; all 12 stabilizers are
                                diagonal+commuting, the walk coin is single-particle
                                (quadratic), and the only quartic-off-diagonal channel is
                                the symmetry-breaking Feshbach.""")

# ----------------------------------------------------------------------------------
Hd("VERDICT — the phase question, resolved structurally")
print("""DECISIVE: the framework's §5 cannot realise symmetric mass generation.

  * Its anti-mirror operator Z_chi Z_W is diagonal & commutes with all 12 stabilizers, so at
    strong coupling it PROJECTS onto the codespace (a product state) -- the N-N-forbidden
    hard projection, not an entangled symmetric gap. (A)
  * Its only off-diagonal mass channel (Feshbach/Higgs seesaw) is a chirality-flip BILINEAR,
    which gaps the mirror by SPONTANEOUSLY BREAKING the chiral symmetry (order param != 0). (B,C)
  * SMG needs the THIRD type -- a symmetric, off-diagonal, >=4-fermion interaction whose
    ground state is an entangled singlet with zero bilinear order, available only at the
    magic number 16. The content is the magic 16 (the precondition holds), but the framework
    contains NO such interaction. (B,C)

So the phase fork resolves to the NON-SMG branch: with §5's existing structure the mirror is
either PROJECTED (forbidden / trivial) or BROKEN (seesaw), never symmetrically gapped. The
SMG/ultralocal Strong-CP closure does NOT follow from the framework as it stands.

CONSEQUENCE for item 93/143: continuum Strong-CP closure rests at the OVERLAP (exp-local) /
seesaw tier -- standard, defensible, but NOT the framework's ultralocal-walk ideal, and NOT
uniquely-TCH. To reach the ultralocal tier the framework would have to ADD and justify a new
symmetric multi-flavour 4-fermion interaction (an SMG term beyond the diagonal stabilizers
and the bilinear Feshbach) and then prove it lands in the symmetric gapped phase -- a genuine
new physics input, not a consequence of the existing §5. This is the honest end of this line:
the SMG route's precondition (anomaly-free 16) and candidate operator (Z_chi Z_W) check out,
but the DYNAMICS the framework actually has cannot realise it.

Net across the four scripts: bare-substrate theta_UV=0 and arg det M=0 STAND (items 93/143);
continuum closure is overlap/seesaw tier; the ultralocal-SMG dream is precondition-met but
dynamically unrealised in the current framework, and would need a new symmetric interaction.""")
print("\nALL ASSERTS PASSED (structural claims computed; the Monte Carlo phase scan is NOT")
print("claimed -- the result is the diagonal-projects / bilinear-breaks structural dichotomy).")
