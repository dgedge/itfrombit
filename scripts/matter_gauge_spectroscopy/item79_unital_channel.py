#!/usr/bin/env python3
"""Item 79 (i)-leg, constructive half: EQUIPARTITION DERIVED from the monitored bridge channel.

K9 showed the unitary walk cannot uniformise the 137-state sector (integrable, irregular graph).
This script runs the sharpened promotion criterion: make the bridge an OPEN system under the
canon's own QEC reading (the lattice extracts syndromes = monitors positions every tick), i.e.
site-basis dephasing jumps L_s = sqrt(gamma)|s><s| on the 137-state space (136 confined pair
states + 1 emission mode coupled at the gauge bridge).

Theorem structure (exact, then verified numerically):
  Lemma 1 (unitality): Hermitian jumps => Lindbladian unital => I/137 is stationary.
  Lemma 2 (uniqueness, Evans 1977/Frigerio 1978/Spohn 1977): fixed points of a unital
    Lindbladian = commutant of *-alg{H, L_k}. Commutant of all |s><s| = diagonal algebra;
    a diagonal X commutes with H iff X_ss = X_tt whenever H_st != 0 — i.e. constant on the
    connected components of the H-graph. H-graph connected => commutant = C.I => UNIQUE.
  Corollary: rho_inf = I/137; P(emission mode) = 1/137 = alpha_0. Equipartition is DERIVED
    (microcanonical assumption replaced by an open-channel fixed-point theorem).
Honest scope: this closes the EQUIPARTITION half only. The 137 itself is the state-space
convention (K9 finding 1) — the same theorem on the antisymmetric (fermionic) space gives
P(em)=1/121. Both are run below to pin the residual gap exactly.
Self-asserting; exit 0 = every quoted number verified."""
import itertools, numpy as np

# ---- bridge graph (as in item79_bridge_trace.py) ----
N = 16
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1

def build_space(symmetric):
    """pair basis (with diagonal if symmetric=True), pair Hamiltonian, portal vector."""
    if symmetric: pairs = [(i,j) for i in range(N) for j in range(i,N)]
    else:         pairs = [(i,j) for i in range(N) for j in range(i+1,N)]
    idx = {p:k for k,p in enumerate(pairs)}
    d   = len(pairs)
    Bas = np.zeros((N*N,d))
    for (i,j),k in idx.items():
        v = np.zeros((N,N))
        if i==j: v[i,i]=1.0
        else:
            v[i,j] = 1/np.sqrt(2); v[j,i] = (1 if symmetric else -1)/np.sqrt(2)
        Bas[:,k] = v.reshape(-1)
    Hfull = np.kron(A,np.eye(N)) + np.kron(np.eye(N),A)
    H2 = Bas.T @ Hfull @ Bas
    # portal: the cross-bridge pairs (one walker each side of the gauge bridge interface)
    portal_pairs = [(0,8),(0,9),(1,8),(1,9)]
    port = np.zeros(d)
    for p in portal_pairs: port[idx[(min(p),max(p))]] = 1.0
    port /= np.linalg.norm(port)
    return pairs, H2, port

def run(symmetric, w=0.5, gamma=1.0, label=""):
    pairs, H2, port = build_space(symmetric)
    d = len(pairs); D = d+1                       # + emission mode
    H = np.zeros((D,D)); H[:d,:d] = H2
    H[:d,d] = w*port; H[d,:d] = w*port            # bidirectional gauge coupling (microcanonical)
    # Lemma 2 connectivity: H-graph connected?
    adj = (np.abs(H) > 1e-12); np.fill_diagonal(adj, False)
    seen = {0}; stack=[0]
    while stack:
        u = stack.pop()
        for v in np.nonzero(adj[u])[0]:
            if v not in seen: seen.add(v); stack.append(v)
    connected = (len(seen)==D)
    # diagonal-commutant check (constructive Lemma 2): union-find over H-edges
    parent = list(range(D))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for u in range(D):
        for v in np.nonzero(adj[u])[0]:
            ru,rv = find(u),find(int(v))
            if ru!=rv: parent[ru]=rv
    ncomp = len({find(x) for x in range(D)})
    # unitality: dephasing L_s = sqrt(g)|s><s| -> D(rho) = g(diag(rho)-rho); L(I) = -i[H,I] + 0 = 0
    # evolve from the physical bridge start sym/antisym(0,8) — the state K9 showed CANNOT
    # uniformise unitarily — plus 2 random starts (operational uniqueness)
    i08 = pairs.index((0,8))
    rng = np.random.default_rng(7)
    starts = []
    e0 = np.zeros(D); e0[i08]=1.0; starts.append(np.outer(e0,e0))
    for _ in range(2):
        Vr = rng.standard_normal((D,D)) + 1j*rng.standard_normal((D,D))
        r0 = Vr@Vr.conj().T; starts.append(r0/np.trace(r0).real)
    uni = np.eye(D)/D
    devs, pems = [], []
    for rho in starts:
        rho = rho.astype(complex)
        dt, T = 0.02, 400.0
        for _ in range(int(T/dt)):                # RK4 on rho' = -i[H,rho] + g(diag(rho)-rho)
            def f(r):
                c = -1j*(H@r - r@H)
                return c + gamma*(np.diag(np.diag(r)) - r)
            k1=f(rho); k2=f(rho+dt/2*k1); k3=f(rho+dt/2*k2); k4=f(rho+dt*k3)
            rho = rho + dt/6*(k1+2*k2+2*k3+k4)
        devs.append(np.max(np.abs(rho-uni)))
        pems.append(rho[d,d].real)
    print(f"  {label}: D={D}  H-graph connected={connected}  diagonal-commutant components={ncomp}")
    print(f"     after t=400 (3 starts incl. the K9 bridge state): max|rho - I/{D}| = {max(devs):.2e}")
    print(f"     P(emission) = {np.mean(pems):.7f}  vs 1/{D} = {1/D:.7f}")
    assert connected and ncomp==1
    assert max(devs) < 5e-4
    assert all(abs(p-1/D) < 1e-4 for p in pems)
    return D

print("UNITAL-CHANNEL ROUTE: monitored bridge => unique fixed point I/D => P(em)=1/D derived\n")
D_sym  = run(True,  label="T(16) space (canon convention)")
print()
D_anti = run(False, label="antisymmetric space (Grassmann-consistent)")
print(f"""
ROBUSTNESS (fixed point independent of rates/couplings — only connectivity matters):""")
for w,g in ((0.1,0.3),(2.0,5.0)):
    pairs,H2,port = build_space(True); d=len(pairs); D=d+1
    H = np.zeros((D,D)); H[:d,:d]=H2; H[:d,d]=w*port; H[d,:d]=w*port
    rho = np.zeros((D,D),complex); rho[pairs.index((0,8)),pairs.index((0,8))]=1.0
    dt,T = 0.02, 600.0
    for _ in range(int(T/dt)):
        def f(r):
            return -1j*(H@r-r@H) + g*(np.diag(np.diag(r))-r)
        k1=f(rho);k2=f(rho+dt/2*k1);k3=f(rho+dt/2*k2);k4=f(rho+dt*k3)
        rho=rho+dt/6*(k1+2*k2+2*k3+k4)
    dev = np.max(np.abs(rho-np.eye(D)/D))
    print(f"  w={w}, gamma={g}: max|rho - I/137| = {dev:.2e}")
    assert dev < 5e-3

print(f"""
RESULT:
  EQUIPARTITION IS DERIVED on the monitored bridge: the dephasing (syndrome-extraction)
  channel is unital, and the H-graph is connected, so by Evans-Frigerio the UNIQUE fixed
  point is maximally mixed — the same bridge start that the unitary walk provably could not
  uniformise (K9) now converges to I/137, giving P(emission)=1/137 with no measure assumed.
  The dissipator supplies exactly the thermalisation mechanism the integrable walk lacked.
HONEST RESIDUAL (unchanged by this result):
  the value 137 is the STATE-SPACE CONVENTION, not an output — the identical theorem on the
  Grassmann-consistent antisymmetric space gives P(em)=1/{D_anti}. What remains of item 79
  (i)-leg: (a) derive the unordered-with-diagonal pair space (the 136) from the substrate,
  or accept 1/121; (b) the Sum(gamma_k)=alpha*Lambda normalisation (the (ii)-leg) is untouched.
ALL ASSERTS PASSED""")
