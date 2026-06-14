#!/usr/bin/env python3
"""Completion of the user's truncated paste (their EXACT functions) + KS scan & center checks,
to independently reproduce the reported gaps/spreads. Verification, not new construction."""
import numpy as np
def gap(ev): return float(ev[1]-ev[0])

# ---- user's exact functions ----
def u1_rotor_eigs(n_max, beta, g2):
    ns=np.arange(-n_max,n_max+1); H=np.diag((g2/2)*ns**2).astype(complex)
    for i in range(2*n_max):
        H[i,i+1]+=-beta/2; H[i+1,i]+=-beta/2
    return np.linalg.eigvalsh(H)
def su2_character_eigs(n_max, beta, g2, center=1):
    H=np.diag([g2*(n/2)*(n/2+1) for n in range(n_max+1)]).astype(complex)
    for n in range(n_max):
        H[n,n+1]+=-beta*center; H[n+1,n]+=-beta*center
    return np.linalg.eigvalsh(H)
def su3_reps(level): return [(p,q) for p in range(level+1) for q in range(level+1-p)]
def su3_casimir(p,q): return (p*p+q*q+p*q+3*p+3*q)/3
def su3_character_eigs(level, beta, g2, center=1):
    reps=su3_reps(level); idx={r:i for i,r in enumerate(reps)}; d=len(reps)
    M=np.zeros((d,d),complex)
    for col,(p,q) in enumerate(reps):
        for t in [(p+1,q),(p-1,q+1),(p,q-1)]:
            if t in idx: M[idx[t],col]+=1
    E=np.diag([g2*su3_casimir(*r) for r in reps]).astype(complex)
    Mag=-(beta/2)*(center*M+np.conj(center)*M.conj().T)
    H=E+Mag; assert np.linalg.norm(H-H.conj().T)<1e-10
    return np.linalg.eigvalsh(H)

# ---- KS weak-coupling trajectory: beta = 1/g^2  =>  g2 = 1/beta ----
print("KS trajectory beta=1/g^2 (g2=1/beta); gap convergence (two truncations agree => converged):")
for beta in [1.0, 10.0, 100.0]:
    g2=1.0/beta
    u=[gap(u1_rotor_eigs(N,beta,g2)) for N in (400,700)]
    s2=[gap(su2_character_eigs(N,beta,g2)) for N in (400,700)]
    print(f"  beta={beta:6}: U(1) gap={u[1]:.5f} (Nconv {abs(u[0]-u[1]):.1e}) | "
          f"SU(2) gap={s2[1]:.5f} (Nconv {abs(s2[0]-s2[1]):.1e})")
for beta in [1.0, 10.0, 30.0]:
    g2=1.0/beta
    s3=[gap(su3_character_eigs(L,beta,g2)) for L in (70,90)]
    print(f"  beta={beta:6}: SU(3) gap={s3[1]:.5f} (Lconv {abs(s3[0]-s3[1]):.1e})")

print("\nCenter-sector isospectrality (spread over center phases; should be ~0):")
# SU(2): center = +/-1
ev_p=su2_character_eigs(400,10.0,0.1,center=1); ev_m=su2_character_eigs(400,10.0,0.1,center=-1)
print(f"  SU(2) center spread (z=+1 vs -1) = {np.max(np.abs(ev_p-ev_m)):.2e}")
# SU(3): center = cube roots of unity
z3=[1, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)]
evs=[su3_character_eigs(60,10.0,0.1,center=z) for z in z3]
spread=max(np.max(np.abs(evs[0]-evs[k])) for k in (1,2))
print(f"  SU(3) Z3 center spread = {spread:.2e}")
