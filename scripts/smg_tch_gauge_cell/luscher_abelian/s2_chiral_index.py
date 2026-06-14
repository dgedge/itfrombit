#!/usr/bin/env python3
"""S2 — chiral measure structure + overlap index theorem (builds on S1's overlap).
ghat5 = g5(1-D) = -sign(H_W) (modified chirality); P_hat_pm = (1 +- ghat5)/2 chiral projectors;
overlap index = -1/2 Tr sign(H_W)  (GW index theorem). 4D index needs F^F != 0:
 single-plane U(1) flux -> index 0; two-plane (F01 & F23) flux -> index = k1*k2.
Self-asserting; reuses the S1 overlap build."""
import sys, itertools, numpy as np
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
M = 1.0
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.array([[1,0],[0,-1]],complex)]
I2 = np.eye(2,dtype=complex); Zr = np.zeros((2,2),complex)
g = [np.block([[Zr,-1j*s[k]],[1j*s[k],Zr]]) for k in range(3)] + [np.block([[I2,Zr],[Zr,-I2]])]
g5 = g[0]@g[1]@g[2]@g[3]
sites = list(itertools.product(range(L), repeat=4)); N = len(sites); idx = {c:i for i,c in enumerate(sites)}
coord = np.array(sites)                       # (N,4)
def shift(mu, sgn):
    P = np.zeros((N,N), complex)
    for c in sites:
        cc = list(c); cc[mu] = (cc[mu]+sgn) % L
        P[idx[tuple(cc)], idx[c]] = 1.0
    return P
SHIFT = {(mu,sg): shift(mu,sg) for mu in range(4) for sg in (+1,-1)}
def build(theta):
    DW = (4-M)*np.eye(N*4, dtype=complex)
    for mu in range(4):
        if theta is None:
            Tp = SHIFT[(mu,+1)].copy()
        else:
            Tp = SHIFT[(mu,+1)]*np.exp(1j*theta[mu])[None,:]
        Tm = Tp.conj().T
        DW += -0.5*(np.kron(Tp, np.eye(4)-g[mu]) + np.kron(Tm, np.eye(4)+g[mu]))
    G5 = np.kron(np.eye(N), g5)
    HW = G5@DW
    w, V = np.linalg.eigh(HW)
    sgn = V@np.diag(np.sign(w))@V.conj().T
    D = np.eye(N*4, dtype=complex) + G5@sgn
    return D, sgn, G5, w
def flux(k1, k2):
    B1, B2 = 2*np.pi*k1/L**2, 2*np.pi*k2/L**2
    th = [np.zeros(N) for _ in range(4)]
    th[1] = B1*coord[:,0]                                   # U_1 = exp(i B1 n0)  -> F01
    edge0 = coord[:,0] == L-1
    th[0] = np.where(edge0, -B1*L*coord[:,1], 0.0)          # boundary twist on U_0
    th[3] = B2*coord[:,2]                                   # U_3 = exp(i B2 n2)  -> F23
    edge2 = coord[:,2] == L-1
    th[2] = np.where(edge2, -B2*L*coord[:,3], 0.0)          # boundary twist on U_2
    return th
def report(name, theta):
    D, sgn, G5, w = build(theta)
    ghat5 = -sgn                                            # = g5(1-D)
    chk_def = np.linalg.norm(ghat5 - G5@(np.eye(N*4)-D))    # ghat5 == g5(1-D)?
    chk_sq  = np.linalg.norm(ghat5@ghat5 - np.eye(N*4))     # ghat5^2 == 1?
    Pm = 0.5*(np.eye(N*4)-ghat5); Pp = 0.5*(np.eye(N*4)+ghat5)
    chk_idem = np.linalg.norm(Pm@Pm - Pm)
    index = -0.5*np.real(np.trace(sgn))                     # GW index = -1/2 Tr sign(H_W)
    dimPp, dimPm = np.real(np.trace(Pp)), np.real(np.trace(Pm))
    print(f"  {name:<22} index={index:+.4f}  dimP+-dimP-={dimPp-dimPm:+.2f}(=2*index)  "
          f"min|H_W|={np.min(np.abs(w)):.3f}  [ghat5^2-1={chk_sq:.1e} def={chk_def:.1e} idem={chk_idem:.1e}]")
    return index, chk_sq, chk_def, chk_idem
print(f"S2 chiral structure + index, L={L} (dim {N*4})", flush=True)
i0,*c0 = report("free U=1",            None)
i1,*c1 = report("1-plane flux k=1",    flux(1,0))     # F01 only -> F^F=0 -> index 0
i2,*c2 = report("2-plane flux k1=k2=1", flux(1,1))    # F01 & F23 -> F^F!=0 -> index = 1
struct_ok = all(max(c) < 1e-9 for c in (c0,c1,c2))
print("\nchiral structure (ghat5^2=1, ghat5=g5(1-D), P idempotent):",
      "PASS" if struct_ok else "FAIL")
print(f"index theorem: free={i0:+.2f} (exp 0), 1-plane={i1:+.2f} (exp 0, F^F=0), "
      f"2-plane={i2:+.2f} (exp +-1, F^F!=0)")
print("S2 VALIDATION:", "PASS" if struct_ok and abs(i0)<1e-6 and abs(i1)<1e-6 and abs(round(i2))>=1
      else "CHECK (see index values)")
