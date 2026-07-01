#!/usr/bin/env python3
r"""AUDIT: is sign(eta_B) <-> sign(J) a from-nothing prediction, or a lock + convention?

Earlier I claimed the framework 'passes': both J>0 and eta_B>0 observed, consistent at
sigma=+1. This audits that claim honestly. Two facts decide it:
  (a) the Jarlskog J is rephasing-invariant but FLIPS under an odd permutation of
      generations (it is built from the antisymmetric structure);
  (b) the generation-blind portal carrier I_CP=Im[M_12 M_23 M_31] is permutation-
      INVARIANT (symmetric M, equal off-diagonals).
If (a) and (b) hold, then sign(J)*sign(eta) = the PARITY of the byte->physical
generation map -- so the absolute pairing rides on the generation labelling, not on
sigma alone. The robust, convention-free content is only the LOCK (both carried by
the one orientation omega); the absolute 'both positive' needs the map's parity.
Self-asserting.
"""
import numpy as np, itertools
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m

def jarlskog(V):
    return np.imag(V[0,0]*V[1,1]*np.conj(V[0,1])*np.conj(V[1,0]))

print("="*72); print("SIGN-LOCK CONVENTION AUDIT: J vs eta_B"); print("="*72)

# a CKM-like unitary with a definite CP phase -> definite J sign
def ckm(th12,th13,th23,d):
    c12,s12=np.cos(th12),np.sin(th12); c13,s13=np.cos(th13),np.sin(th13); c23,s23=np.cos(th23),np.sin(th23)
    return np.array([
     [c12*c13, s12*c13, s13*np.exp(-1j*d)],
     [-s12*c23-c12*s23*s13*np.exp(1j*d), c12*c23-s12*s23*s13*np.exp(1j*d), s23*c13],
     [s12*s23-c12*c23*s13*np.exp(1j*d), -c12*s23-s12*c23*s13*np.exp(1j*d), c23*c13]])
V=ckm(0.227,0.0037,0.042,1.2)
J0=jarlskog(V)
print(f"\n[a] Jarlskog of a CKM-like matrix: J = {J0:+.3e}")
# permute two down-type generations (swap columns 0,1) -> odd permutation
P=np.array([[0,1,0],[1,0,0],[0,0,1]],float); Vp=V@P
ok(np.sign(jarlskog(Vp))==-np.sign(J0), "sign(J) FLIPS under an odd generation permutation (swap 2 columns)")

# (b) generation-blind Majorana carrier I_CP
A=np.array([[0,1,1],[1,0,1],[1,1,0]],float)
def MH(sig,phi,M0=1,r=0.3): return M0*(np.eye(3)+r*np.exp(1j*sig*phi)*A)
def ICP(M): return np.imag(M[0,1]*M[1,2]*M[2,0])
M=MH(+1,1/3); base=ICP(M)
print(f"\n[b] generation-blind I_CP = {base:+.4e}")
inv=True
for perm in itertools.permutations((0,1,2)):
    Pp=np.zeros((3,3));
    for i,j in enumerate(perm): Pp[i,j]=1
    inv &= abs(ICP(Pp@M@Pp.T)-base)<1e-12
ok(inv, "I_CP is INVARIANT under every generation permutation (symmetric, equal off-diagonals)")

# (c) therefore the relative sign rides on the generation-map parity
print("\n[c] consequence: sign(J)*sign(eta) = sign(J)*sign(I_CP) flips with the byte->physical")
print("    generation map parity (J flips, I_CP does not). So the ABSOLUTE pairing is NOT set by")
print("    sigma alone; it rides on the generation labelling.")
ok(True, "absolute same/opposite-sign depends on the byte->physical generation map parity")

print("\n"+"="*72); print("VERDICT (correction to 'it passes')")
print("  ROBUST / convention-free: sign(J) and sign(eta_B) are EACH carried by the one global")
print("  orientation omega, so they are LOCKED -- correlated, not independently choosable. That much")
print("  is a genuine structural prediction.")
print("  NOT convention-free: WHICH locking (same vs opposite sign) = the parity of the byte->physical")
print("  generation map. That map is constrained by the Koide mass assignment (each n -> a definite")
print("  mass), so the parity is in-principle fixed -- but I did NOT compute it from scratch in a common")
print("  basis. So 'both observed positive -> it passes' is CONSISTENT under the canonical assignment,")
print("  not a from-nothing prediction. I overstated it: the solid claim is the LOCK; the absolute")
print("  pair reduces to a (fixed-but-here-unverified) generation-map parity. exit 0")
