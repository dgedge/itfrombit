"""Global Z3 colour-symmetry audit.

Checks that the internal walk sectors commute with the colour Z3 action where
expected, and records the remaining bridge-rotation limitation.
"""

import numpy as np
I2=np.eye(2); I4=np.eye(4)
def comm(A,B): return A@B-B@A
def zero(M): return np.allclose(M,0)

# colour qubit space {00,01,10,11} = idx {0,1,2,3}; quark triplet = {01,10,11}={1,2,3}, 00=singlet
# Z3 = GL(2,F2) gen [[0,1],[1,1]] : 01->11->10->01, fixes 00  (perm[in]=out)
perm=[0,3,1,2]
P=np.zeros((4,4));
for i,o in enumerate(perm): P[o,i]=1
# V_strong on colour: A_{K3} on the triplet {1,2,3}, 00 isolated (sec 3.3)
AK3=np.zeros((4,4))
for a in (1,2,3):
    for b in (1,2,3):
        if a!=b: AK3[a,b]=1

print("=== colour-factor (the only nontrivial commutator) ===")
print("[V_strong=A_K3 , Z3] = 0 :", zero(comm(AK3,P)))
print("Z3 is order 3 :", np.allclose(np.linalg.matrix_power(P,3), I4))

# full internal walk on colour(4) x isospin(2) x chirality(2) = 16-dim
def kron(*ms):
    out=ms[0]
    for m in ms[1:]: out=np.kron(out,m)
    return out
Z3f   = kron(P,   I2, I2)                 # Z3 on colour only
Vstr  = kron(AK3, I2, I2)                 # strong: colour
# coin C: zero-controlled CNOT, flip I3 (isospin) when chi=0  (sec 3.1)
X=np.array([[0,1],[1,0]]); P0=np.array([[1,0],[0,0]]); P1=np.array([[0,0],[0,1]])
coin  = kron(I4, X, P0) + kron(I4, I2, P1)        # iso flips iff chi=0
# V_weak: colour-free chirality-mixing (sec 3.2) - any op on iso x chi, NOT colour
Vweak = kron(I4, I2, X) + kron(I4, X, I2)          # representative colour-free
# V_em: diagonal in Q = I3 - 1/2(1-LQ): colour-BLIND diagonal (sec 3.2)
Vem   = kron(I4, np.diag([0.5,-0.5]), I2)          # depends on isospin, colour-blind

print("\n=== full internal walk vs colour Z3 (16-dim) ===")
for name,Op in [("coin C",coin),("V_strong",Vstr),("V_weak",Vweak),("V_em",Vem)]:
    print(f"  [{name:9s}, Z3] = 0 :", zero(comm(Op,Z3f)))
Wint = coin @ (Vstr+Vweak+Vem)             # internal hop content (no spatial R_d/shift)
print("  [coin·(V_strong+V_weak+V_em), Z3] = 0 :", zero(comm(Wint,Z3f)))
print("\n=> internal coin+potential sector is Z3-symmetric.")
print("   Remaining factor in T_d = -(i/√3) R_d (V_em+V_weak+V_strong): the bridge rotation R_d.")
print("   [W,Z3]=0  <=>  [R_d,Z3]=0  <=>  the shift preserves the colour sector {C0,C1}")
print("   = exactly COLOUR CONSERVATION = the sec 2.9 colour-sector F2 closure.")
