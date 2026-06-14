#!/usr/bin/env python3
"""S1 — 4D overlap (Neuberger) operator + Ginsparg-Wilson check, free and in a U(1) background.
Foundation for the Luscher abelian chiral-measure build (scope.md). Self-asserting.
D = I + G5 sign(H_W),  H_W = G5 D_W (hermitian),  D_W = (4-M)I - (1/2) sum_mu[(I-g_mu)U_mu T+ + (I+g_mu)U_mu^dag T-].
GW: {G5,D} = D G5 D ;  g5-herm: G5 D G5 = D^dag.  (Both must hold for ANY background U.)"""
import sys, itertools, numpy as np
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
M = 1.0
# Euclidean gammas
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.array([[1,0],[0,-1]],complex)]
I2 = np.eye(2,dtype=complex); Zr = np.zeros((2,2),complex)
g = [np.block([[Zr,-1j*s[k]],[1j*s[k],Zr]]) for k in range(3)] + [np.block([[I2,Zr],[Zr,-I2]])]
g5 = g[0]@g[1]@g[2]@g[3]
for mu in range(4):
    for nu in range(4):
        assert np.allclose(g[mu]@g[nu]+g[nu]@g[mu], 2*(mu==nu)*np.eye(4)), "gamma algebra"
assert np.allclose(g5@g5, np.eye(4)) and np.allclose(g5, g5.conj().T), "g5"
sites = list(itertools.product(range(L), repeat=4)); N = len(sites); idx = {c:i for i,c in enumerate(sites)}
def shift(mu, sgn, phase=None):
    P = np.zeros((N,N), complex)
    for c in sites:
        cc = list(c); cc[mu] = (cc[mu]+sgn) % L
        ph = 1.0 if phase is None else phase[mu][idx[c] if sgn>0 else idx[tuple(cc)]]
        P[idx[tuple(cc)], idx[c]] = ph
    return P
def overlap(theta=None):
    DW = (4-M)*np.eye(N*4, dtype=complex)
    for mu in range(4):
        if theta is None:
            Tp, Tm = shift(mu,+1), shift(mu,-1)
        else:
            U = np.exp(1j*theta[mu])                  # link phase per site (forward link from site)
            Tp = shift(mu,+1); Tp = Tp*U[None,:]      # multiply forward hop by U_mu(x)
            Tm = Tp.conj().T                           # backward = dagger
        DW += -0.5*(np.kron(Tp, np.eye(4)-g[mu]) + np.kron(Tm, np.eye(4)+g[mu]))
    G5 = np.kron(np.eye(N), g5)
    HW = G5@DW
    assert np.allclose(HW, HW.conj().T, atol=1e-10), "H_W hermitian"
    w, V = np.linalg.eigh(HW)
    sgn = V@np.diag(np.sign(w))@V.conj().T
    D = np.eye(N*4, dtype=complex) + G5@sgn
    gw = np.linalg.norm(G5@D + D@G5 - D@G5@D)
    herm = np.linalg.norm(G5@D@G5 - D.conj().T)
    return gw, herm, float(np.min(np.abs(w)))
print(f"S1 overlap GW check, L={L} (dim {N*4})", flush=True)
gw0, h0, gap0 = overlap(None)
print(f"  free U=1     : GW residual={gw0:.2e}  g5-herm={h0:.2e}  min|H_W eig|={gap0:.3f}")
rng = np.random.default_rng(7)
theta = [rng.uniform(-0.4, 0.4, N) for _ in range(4)]      # generic U(1) background
gwU, hU, gapU = overlap(theta)
print(f"  U(1) bgnd    : GW residual={gwU:.2e}  g5-herm={hU:.2e}  min|H_W eig|={gapU:.3f}")
ok = gw0 < 1e-9 and h0 < 1e-9 and gwU < 1e-9 and hU < 1e-9
print("S1 VALIDATION:", "PASS (GW + g5-herm hold, free and gauged)" if ok else "FAIL")
