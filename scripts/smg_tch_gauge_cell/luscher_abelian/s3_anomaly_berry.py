#!/usr/bin/env python3
"""S3 — dynamical abelian anomaly cancellation via the determinant-line-bundle Berry phase.
The measure obstruction is the CONSISTENT cubic anomaly ~Y^3 (the chiral-determinant phase),
not the index (~Y^2). Basis-independent handle: Berry phase around a closed loop {A,B,C} of U(1)
backgrounds, Berry(Y) = -arg[ det(WA^dag WB) det(WB^dag WC) det(WC^dag WA) ], W = ortho basis of
P_hat_-(Y*theta). Family index theorem => curvature ~Y^3 exactly (abelian) => Berry(Y)/Y^3 = const,
and sum_i mult_i Berry(Y_i) = const * sum Y^3 = 0 iff anomaly-free. Tests BOTH the cubic structure
and the cancellation. SM-16 hypercharges below (sum Y^3 = 0)."""
import sys, itertools, numpy as np
np.seterr(divide="ignore", invalid="ignore")            # Y=0 gives 0/0 in the ratio; handled explicitly
L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
M = 1.0
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]],complex), np.array([[1,0],[0,-1]],complex)]
I2 = np.eye(2,dtype=complex); Zr = np.zeros((2,2),complex)
g = [np.block([[Zr,-1j*s[k]],[1j*s[k],Zr]]) for k in range(3)] + [np.block([[I2,Zr],[Zr,-I2]])]
g5 = g[0]@g[1]@g[2]@g[3]
sites = list(itertools.product(range(L), repeat=4)); N = len(sites); idx = {c:i for i,c in enumerate(sites)}
def shift(mu, sgn):
    P = np.zeros((N,N), complex)
    for c in sites:
        cc = list(c); cc[mu] = (cc[mu]+sgn) % L
        P[idx[tuple(cc)], idx[c]] = 1.0
    return P
SHIFTP = {mu: shift(mu,+1) for mu in range(4)}
def subspace(theta):
    """orthonormal basis W of P_hat_-(theta) range, the index, and min|H_W|."""
    DW = (4-M)*np.eye(N*4, dtype=complex)
    for mu in range(4):
        Tp = SHIFTP[mu]*np.exp(1j*theta[mu])[None,:]
        Tm = Tp.conj().T
        DW += -0.5*(np.kron(Tp, np.eye(4)-g[mu]) + np.kron(Tm, np.eye(4)+g[mu]))
    G5 = np.kron(np.eye(N), g5)
    HW = G5@DW
    w, V = np.linalg.eigh(HW)
    sgn = V@np.diag(np.sign(w))@V.conj().T
    ghat5 = -sgn
    ev, U = np.linalg.eigh(ghat5)
    W = U[:, ev < 0]                                   # P_hat_- range
    index = int(round(-0.5*np.real(np.trace(sgn))))
    return W, index, float(np.min(np.abs(w)))
def berry(configs, Y):
    Ws, ixs, mhs = [], [], []
    for cf in configs:
        W, ix, mh = subspace([Y*cf[mu] for mu in range(4)]); Ws.append(W); ixs.append(ix); mhs.append(mh)
    if len(set(ixs)) != 1 or len(set(W.shape[1] for W in Ws)) != 1:
        return None, ixs, min(mhs), 0.0                # rank/index mismatch -> skip
    phase, logabs = 1.0 + 0j, 0.0                       # slogdet: stable phase + |loop|
    for a, b in ((0, 1), (1, 2), (2, 0)):
        sg, la = np.linalg.slogdet(Ws[a].conj().T @ Ws[b])
        phase *= sg; logabs += la
    return float(-np.angle(phase)), ixs[0], min(mhs), float(np.exp(logabs))
rng = np.random.default_rng(11)
base = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
dlt  = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
A0 = [rng.uniform(-base, base, N) for _ in range(4)]                          # generic base background
d1 = [rng.uniform(-dlt, dlt, N) for _ in range(4)]                            # small perturbation 1
d2 = [rng.uniform(-dlt, dlt, N) for _ in range(4)]                            # small perturbation 2
configs = [A0, [A0[mu]+d1[mu] for mu in range(4)], [A0[mu]+d2[mu] for mu in range(4)]]  # small closed loop
charges = [(1/6, 6), (-2/3, 3), (1/3, 3), (-1/2, 2), (1.0, 1), (0.0, 1)]
print(f"S3 anomaly-Berry, L={L} (dim {N*4}), base={base} delta={dlt} (small loop -> |loop|~1)", flush=True)
print(f"  sum Y^3 (SM 16) = {sum(m*Y**3 for Y,m in charges):+.6e}")
tot, tot_abs, ratios = 0.0, 0.0, []
for Y, mult in charges:
    b, ix, mh, mod = berry(configs, Y)
    if b is None:
        print(f"  Y={Y:+.4f} x{mult}: SKIP (index/rank {ix})"); continue
    r = b/Y**3 if abs(Y) > 1e-12 else float('nan')
    if abs(Y) > 1e-12: ratios.append(r)
    tot += mult*b; tot_abs += mult*abs(b)
    print(f"  Y={Y:+.4f} x{mult}: Berry={b:+.4e}  Berry/Y^3={r:+.4e}  index={ix} min|H_W|={mh:.3f} |loop|={mod:.3f}")
rr = np.array(ratios)
print(f"\ncubic structure: Berry/Y^3 = {rr.mean():+.4e} +- {rr.std():.2e}  (rel spread {rr.std()/abs(rr.mean()):.1%})")
print(f"cancellation:    sum mult*Berry = {tot:+.4e}   sum mult*|Berry| = {tot_abs:.4e}   ratio = {abs(tot)/tot_abs:.2e}")
print("S3 VALIDATION:", "PASS (Berry~Y^3 const AND weighted sum cancels)"
      if rr.std()/abs(rr.mean()) < 0.05 and abs(tot)/tot_abs < 0.05 else "CHECK (see numbers)")
