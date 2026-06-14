#!/usr/bin/env python3
"""
The crux (c): does the framework's §5 mass mechanism supply R2 as a SYMMETRIC INTERACTION
that gaps the W!=chi mirror sector (SMG, branch a) -- or only a symmetry-breaking seesaw
(branch b, not SMG)?  Build the actual operators from canon and test.

Grounded inputs (ANCHOR, read 2026-06-01):
  * §3.1/§3.5 coin:  C = zero-controlled CNOT; U(m)=cos m I - i sin m sigma_x^{(I3)},
    controlled by chi.  §5.3 Yukawa:  p = <psi| i C^dag dC/dlambda |psi> = <psi|H_R2|psi>,
    NAMED 'relaxation of the chirality constraint (W=chi)'.
  * §5.2 mass stabilizers:  S_{(i,j)} = Z_i Z_j on the 12 Q_3 edges; F = sum (1-<Z_iZ_j>)/2;
    M(c)=exp(phi F/2), phi=(sqrt5-1)/2.  chi=octant 110, W=octant 111 -> Hamming-adjacent,
    so Z_chi Z_W = S_{(110,111)} is one of the 12 stabilizers.
  * §2.8/§3.2 gauge action: U(1)_Y diagonal; SU(2)_L = chi-controlled I3 ladder
    (V_weak ~ CNOT); SU(3)_c permutes (C0,C1).

Register bit order (ANCHOR §2.1):  0=G0 1=G1 2=LQ 3=C0 4=C1 5=I3 6=chi 7=W.

Tests:
  (A) compute H_R2 = i C^dag dC/dm explicitly  -> show it is the chi-controlled I3 flip,
      W-BLIND ([H_R2, Z_W]=0), and a weak-ISOSPIN generator ([H_R2, T+]!=0). So §5.3's
      operator does NOT act on the W=chi correlation: the 'W=chi relaxation' label does not
      match the operator. H_R2 is therefore NOT the anti-mirror gapping term.
  (B) the anti-mirror operator is the §5.2 stabilizer Z_chi Z_W: +1 on valid (W=chi),
      -1 on mirror (W!=chi); diagonal/quartic; no VEV.
  (C) gauge invariance: [Z_chi Z_W, G]=0 for G in {Y, SU(2)_L T+, SU(3)_c} -> Z_chi Z_W is
      a GAUGE-INVARIANT SYMMETRIC operator (qualifies as an SMG gapping term, branch a).
  (D) the SCALE fork: the §5.2 per-edge mass factor exp(phi/2)~1.36 is NOT a decoupling;
      the GUT-scale mirror decoupling in canon (E_R2~1e15 GeV) rides the §5.8 Type-I SEESAW
      (branch b, symmetry-breaking, NOT SMG). The off-codespace mirror-gap scale is unpinned.

numpy; self-asserting.
"""
import numpy as np, itertools
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)

I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex); Z=np.array([[1,0],[0,-1]],complex)
P0=np.array([[1,0],[0,0]],complex); P1=np.array([[0,0],[0,1]],complex)
SR=np.array([[0,1],[0,0]],complex)        # |0><1| : raises (|1>-> |0>); I3=1 down -> I3=0 up

# bit order 0..7 = G0,G1,LQ,C0,C1,I3,chi,W
G0,G1,LQ,C0,C1,I3,CHI,W = range(8)
def op(singles):                          # singles: dict pos->2x2, else I2 ; bit0 = leftmost
    M=np.array([[1]],complex)
    for k in range(8): M=np.kron(M, singles.get(k,I2))
    return M
def ctrl(cpos,cval,tpos,top):             # controlled op: when bit cpos==cval apply top on tpos
    Pc=P0 if cval==0 else P1; Pn=P1 if cval==0 else P0
    return op({cpos:Pc,tpos:top}) + op({cpos:Pn,tpos:I2})
def comm(A,B): return A@B-B@A
def n(A): return float(np.linalg.norm(A))

# codeword index (bit0 MSB) and validity (R1-R4 from §2.2)
def idx(c): return int("".join(map(str,c)),2)
def R1(c): return not(c[G0]==1 and c[G1]==1)
def R2(c): return c[W]==c[CHI]
def R3(c): return (c[C0],c[C1])==(0,0) if c[LQ]==0 else (c[C0],c[C1])!=(0,0)
allw=[tuple(b) for b in itertools.product([0,1],repeat=8)]
valid48=[c for c in allw if R1(c) and R2(c) and R3(c)]

# ----------------------------------------------------------------------------------
Hd("(A) §5.3 coin operator H_R2 = i C^dag dC/dm  -- what does it actually act on?")
def Cm(m):                                 # coin: on chi=0, U(m)=cos m I - i sin m X_{I3}; chi=1 -> I
    U=np.cos(m)*I2 - 1j*np.sin(m)*X
    return op({CHI:P0,I3:U}) + op({CHI:P1,I3:I2})
m0=0.37; dm=1e-6
Cder=(Cm(m0+dm)-Cm(m0-dm))/(2*dm)
H_R2=1j*Cm(m0).conj().T@Cder
H_R2_expected=op({CHI:P0,I3:X})           # exact: P0_chi (x) X_I3 ; zero on chi=1 block
print(f"  ||H_R2(numeric i C^dag dC/dm) - P0_chi(x)X_I3|| = {n(H_R2-H_R2_expected):.1e}")
print(f"  m-independence: ||H_R2(m=.37)-H_R2(m=1.1)|| = "
      f"{n(1j*Cm(1.1).conj().T@((Cm(1.1+dm)-Cm(1.1-dm))/(2*dm)) - H_R2):.1e}")
assert n(H_R2-H_R2_expected)<1e-5
ZW=op({W:Z})
print(f"  W-blindness:  ||[H_R2, Z_W]|| = {n(comm(H_R2,ZW)):.1e}   (0 => does NOT involve W)")
assert n(comm(H_R2,ZW))<1e-12
Tplus=ctrl(CHI,0,I3,SR)                    # SU(2)_L raising (left = chi=0), §3.2 V_weak
print(f"  weak-isospin: ||[H_R2, T+]|| = {n(comm(H_R2,Tplus)):.2f}   (!=0 => H_R2 is IN the")
print(f"                weak gauge algebra -- a symmetry GENERATOR, not a gapping term)")
assert n(comm(H_R2,Tplus))>0.1
print("  => §5.3 'H_R2 = relaxation of W=chi' is MISLABELLED: the explicit operator is the")
print("     chi-controlled I3 flip (weak isospin), W-blind. It is NOT the anti-mirror gap.")

# ----------------------------------------------------------------------------------
Hd("(B) The anti-mirror operator is the §5.2 stabilizer Z_chi Z_W")
ZcZw=op({CHI:Z,W:Z})
def mirror(c):
    mm=list(c); mm[CHI]^=1; return tuple(mm)
ev_valid=[np.real(ZcZw[idx(c),idx(c)]) for c in valid48]
ev_mirror=[np.real(ZcZw[idx(mirror(c)),idx(mirror(c))]) for c in valid48]
print(f"  Z_chi Z_W on 48 valid (W=chi)  : eigenvalues all = {set(ev_valid)}  (frustration 0)")
print(f"  Z_chi Z_W on 48 mirror (W!=chi): eigenvalues all = {set(ev_mirror)}  (frustration 1)")
assert set(ev_valid)=={1.0} and set(ev_mirror)=={-1.0}
print("  -> (1 - Z_chi Z_W)/2 is exactly the W-chi edge of §5.2's frustration F: the")
print("     stabilizer that penalises mirrors is already one of the 12 mass-code stabilizers.")
print("     It is DIAGONAL (a quartic Z Z in fermion variables) -- no bilinear, no VEV.")

# ----------------------------------------------------------------------------------
Hd("(C) Gauge invariance of Z_chi Z_W (does it qualify as an SMG gapping term?)")
# U(1)_Y diagonal: Y = Q - T3 ; exact values irrelevant to the commutator (any diagonal commutes)
def Zf(c): return 1 if c[I3]==0 else -1
def sumZc(c): return -3 if (c[C0],c[C1])==(0,0) else -1
def Qf(c): return 0.5*Zf(c)+ (1/3)*sumZc(c) + 0.5   # §2.8, matter
def T3(c):  # weak isospin 3-comp: +-1/2 for left (chi=0) doublet members, 0 else
    if c[CHI]!=0: return 0.0
    return 0.5 if c[I3]==0 else -0.5
Ydiag=np.diag([Qf(c)-T3(c) for c in allw]).astype(complex)
# SU(3)_c colour 3-cycle on (C0,C1): 01->10->11->01, 00 fixed (leptons)
perm=np.zeros((256,256),complex)
cyc={(0,1):(1,0),(1,0):(1,1),(1,1):(0,1),(0,0):(0,0)}
for c in allw:
    t=list(c); t[C0],t[C1]=cyc[(c[C0],c[C1])]; perm[idx(tuple(t)),idx(c)]=1
gens={"U(1)_Y (diag)":Ydiag, "SU(2)_L T+":Tplus, "SU(3)_c 3-cycle":perm}
for name,Gv in gens.items():
    print(f"  ||[Z_chi Z_W, {name:16s}]|| = {n(comm(ZcZw,Gv)):.1e}")
    assert n(comm(ZcZw,Gv))<1e-10
print("  -> Z_chi Z_W commutes with ALL SM gauge generators: it is a GAUGE-INVARIANT,")
print("     diagonal/quartic operator that gaps exactly the W!=chi mirror. It QUALIFIES as")
print("     an SMG mirror-gapping term (branch a). (It only touches chi,W; no gauge gen flips")
print("     chirality or weak charge, so commutation is structural.)")

# ----------------------------------------------------------------------------------
Hd("(D) The SCALE fork -- where branch (a) is not (yet) sufficient")
phi=(np.sqrt(5)-1)/2
print(f"  §5.2 mass factor for ONE frustrated edge: exp(phi/2) = {np.exp(phi/2):.3f}")
print( "    -> a W!=chi mirror, as a §5.2 frustration penalty, is only ~1.36x heavier than")
print( "       its valid partner: a FACTOR, NOT a decoupling. The symmetric operator (a)")
print( "       exists and is gauge-invariant, but at the §5.2 mass scale it does not gap")
print( "       the mirror out of the spectrum.")
print( "  The corpus's only CUTOFF-scale mirror penalty is E_R2 ~ 1e15 GeV (line 341), and it")
print( "  is supplied by the §5.8 Type-I SEESAW (v=246, M_R=E_R2) -- a BILINEAR, symmetry-")
print( "  breaking mass (branch b), which is ordinary mass generation, NOT SMG.")
print( "  => the symmetric operator that COULD do SMG (Z_chi Z_W) is at the wrong (small)")
print( "     scale in §5.2; the mechanism that reaches the cutoff is the non-SMG seesaw.")
print( "     The off-codespace mirror-gap scale (code-ENFORCEMENT coupling) is never pinned;")
print( "     if it were the substrate/cutoff scale, branch (a) would close -- but canon does")
print( "     not state it, and uses the seesaw instead.")

# ----------------------------------------------------------------------------------
Hd("VERDICT — does §5 supply R2 as a symmetric SMG interaction?")
print("""PARTIAL, and the fork now resolves with a clear lean:

  YES, the OPERATOR exists and is correct (branch a is real):
    * §5.2's Z_chi Z_W is exactly the anti-mirror stabilizer (+1 valid / -1 mirror), it is
      diagonal/quartic (no bilinear, no VEV), and it commutes with the entire SM gauge group
      -- a textbook-shape SMG gapping term, already one of the 12 mass-code stabilizers. (B,C)

  NO, §5.3's named operator is the WRONG one:
    * H_R2 = i C^dag dC/dm is the chi-controlled I3 flip -- W-blind and a weak-isospin
      GENERATOR (a symmetry, not a gap). The '(W=chi) relaxation' label does not match the
      explicit operator. So the Yukawa/Higgs operator §5 actually writes down is not the
      mirror-gapping term. (A)  [recordable correction]

  NOT SUFFICIENT, at the scale canon uses:
    * As a §5.2 frustration penalty the mirror is only exp(phi/2)~1.36x heavier -- not gapped
      out. The cutoff-scale decoupling in canon (E_R2~1e15 GeV) is the §5.8 Type-I SEESAW --
      BILINEAR and symmetry-breaking (branch b), i.e. ORDINARY mass, not SMG. The symmetric
      operator is at the wrong scale; the cutoff mechanism is the non-symmetric one. (D)

  STILL OPEN (unchanged, hard):
    * whether the strong-coupling Z_chi Z_W dynamics realise the SYMMETRIC gapped phase vs
      spontaneous breaking (the field-wide PMS-vs-SMG problem); and the same-chirality
      residual (walk_kernel_overlap's 4 species), untouched by R2.

NET: §5 contains a genuine, gauge-invariant, symmetric anti-mirror OPERATOR (Z_chi Z_W) --
so the SMG closure is not missing its central object -- but (i) §5.3's written Yukawa
operator is not it (mislabel), (ii) at the scale §5.2/§5.8 actually use, the mirror is
decoupled by the SEESAW (symmetry-breaking, branch b), not by the symmetric operator, and
(iii) the symmetric-gapped-phase question is unaddressed. So §5 does NOT (yet) close the SMG
route: it supplies the right operator at the wrong scale, and its actual cutoff mechanism is
ordinary symmetry-breaking mass. Continuum chirality therefore still rests at the overlap
tier in practice; the ultralocal-SMG closure would require RE-deriving the mirror gap from
Z_chi Z_W at the code-enforcement (cutoff) scale and proving the symmetric phase -- a
specific, well-posed next computation, not a free consequence of §5.""")
print("\nALL ASSERTS PASSED (operators built from canon; verdict states the lean honestly).")
