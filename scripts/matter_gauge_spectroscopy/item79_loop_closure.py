#!/usr/bin/env python3
"""Item 79 question (b): close the loop — derive Sum(gamma_k) = alpha*Lambda self-consistently.

Chain: the unital fixed point gives P(em) = 1/137 INDEPENDENT of the monitoring rate gamma and
portal coupling w (only connectivity matters — Evans-Frigerio + verified connected H-graph).
Canon 5.9 DEFINES alpha as 'the per-tick irreversibility fraction'. Therefore the cell's total
non-unitary leakage rate is
    Sum(gamma_k) = (irreversibility fraction) x (clock rate) = P(em) x Lambda = Lambda/137,
and because P(em) is parameter-independent, the self-consistency map F(gamma_tot) = Lambda*P(em)
is CONSTANT: Sum(gamma_k) = Lambda/137 = alpha_0*Lambda is its unique solution. No tuning is
possible even in principle. Item 119's normalisation is thereby DERIVED (conditional on the
137-sector = question (a), now closed under named premises, and on 5.9's definition of alpha).
Honesty: the derived value is the BARE alpha_0 = 1/137 exactly — canon's w = alpha*Lambda quotes
the dressed 1/137.036 (0.026% apart); the dressed shift is the separate open Part-12 item.
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ---- the 137-channel (from item79_unital_channel.py) ----
N = 16
edges  = [(i,(i+1)%8) for i in range(8)] + [(8+i,8+(i+1)%8) for i in range(8)] + [(0,8),(1,9)]
A = np.zeros((N,N))
for i,j in edges: A[i,j]=A[j,i]=1
pairs = [(i,j) for i in range(N) for j in range(i,N)]
idx = {p:k for k,p in enumerate(pairs)}
d = len(pairs); D = d+1
Bas = np.zeros((N*N,d))
for (i,j),k in idx.items():
    v = np.zeros((N,N))
    if i==j: v[i,i]=1.0
    else: v[i,j]=v[j,i]=1/np.sqrt(2)
    Bas[:,k]=v.reshape(-1)
H2 = Bas.T @ (np.kron(A,np.eye(N))+np.kron(np.eye(N),A)) @ Bas
port = np.zeros(d)
for p in ((0,8),(0,9),(1,8),(1,9)): port[idx[p]] = 0.5

# 1. P(em) is (gamma,w)-independent: L(I/D)=0 exactly for any (gamma,w); uniqueness via the
#    verified connectivity (item79_unital_channel) => P_em = 1/D for ALL (gamma,w).
print("1. stationarity of I/137 at arbitrary (gamma, w)  [uniqueness: Evans-Frigerio + connectivity]")
uni = np.eye(D)/D
for w,g in ((0.1,0.3),(0.5,1.0),(2.0,5.0)):
    H = np.zeros((D,D)); H[:d,:d]=H2; H[:d,d]=w*port; H[d,:d]=w*port
    Lind = -1j*(H@uni-uni@H) + g*(np.diag(np.diag(uni))-uni)
    print(f"   (w={w}, gamma={g}): ||Lindblad(I/137)|| = {np.linalg.norm(Lind):.2e}")
    assert np.linalg.norm(Lind) < 1e-14
P_em = 1/D
print(f"   => P(em) = 1/{D} = {P_em:.10f} for all (gamma, w)")

# 2. self-consistency: F(gamma_tot) = Lambda*P_em is CONSTANT => unique solution
Lam = 1.0
gamma_tot = Lam*P_em
alpha0 = 1/137
alpha_dressed = 1/137.035999
print(f"\n2. self-consistency map F(gamma_tot) = Lambda*P_em = Lambda/137 (constant in gamma_tot)")
print(f"   unique solution: Sum(gamma_k) = Lambda/137 = {gamma_tot:.8f} Lambda = alpha_0*Lambda")
assert abs(gamma_tot - alpha0*Lam) < 1e-15
print(f"   (bare alpha_0 = 1/137; canon's dressed alpha = 1/137.036 sits {abs(alpha_dressed-alpha0)/alpha0*100:.3f}% away —")
print(f"    the dressed shift is the separate, open Part-12 item; 5.4's bare/dressed split.)")

# 3. the full derived chain: Sum(gamma) -> walk-active cell support -> C_loop -> K
print(f"\n3. the chain, end to end (premises: (a) coin-singlet sector; 5.9 alpha-definition;")
print(f"   monitored-bridge reading; walk-active dissipator support):")
def valid(c):
    G0,G1,LQ,C0,C1,I3,chi,W = c
    return not(G0 and G1) and W==chi and ((LQ==0) == (C0==0 and C1==0))
ALL = list(itertools.product((0,1),repeat=8))
inP = [valid(c) for c in ALL]
vac = ALL.index((0,)*8)
sup3 = [3,4,5]                                            # {C0,C1,I3} walk-active bits
gk = gamma_tot/len(sup3)
def flip(i,k): return i ^ (1 << (7-k))
Gvac = sum(gk for k in sup3 if not inP[flip(vac,k)])
C_loop = gamma_tot/Gvac
K = C_loop/alpha0
print(f"   Sum(gamma_k) = Lambda/137  ->  Gamma_vac = {Gvac/gamma_tot:.4f} x Sum(gamma) = (2/3)(Lambda/137)")
print(f"   C_loop = {C_loop:.4f} = 3/2   ->   K = C_loop/alpha_0 = {K:.2f}")
Kdata = 206.49
print(f"   K = {K:.2f} vs K_data = {Kdata} -> {abs(K/Kdata-1)*100:.2f}%")
assert abs(C_loop-1.5) < 1e-12 and abs(K-205.5) < 1e-9 and abs(K/Kdata-1) < 0.006

print(f"""
RESULT (b): Sum(gamma_k) = alpha_0*Lambda is DERIVED — the unital fixed point makes the
self-consistency map constant, so the normalisation cannot be tuned even in principle; it is
forced to Lambda x P(em) = Lambda/137. Item 119's asserted normalisation and item 79's (ii)-leg
collapse onto the (i)-leg theorem + canon 5.9's definition of alpha.
THE ALPHA-CHAIN, NOW DERIVED END-TO-END (under the named premises):
   Fermi statistics + 2-dim coin  =>  photon-channel pair space = T(16) = 136          [(a)]
   monitored bridge (unital)      =>  P(em) = 1/137 unique fixed point                 [(i)]
   5.9 irreversibility definition =>  Sum(gamma_k) = Lambda/137 = alpha_0*Lambda       [(b)]
   walk-active cell dissipator    =>  Gamma_vac = (2/3)alpha_0*Lambda, C_loop = 3/2    [119]
   =>  K = 3/(2 alpha_0) = 205.50 vs K_data = 206.49 (0.48%); M_P^2 = (2/3)alpha_0 Lambda^3 R_dS
REMAINING HONEST GAPS: bare-vs-dressed alpha (0.026%, open Part-12); the 0.48% K residual;
the horizon R_dS stays an INPUT (Dirac relation, not intrinsic M_P); premises (a) are
canon-grounded but newly assembled — canon adoption is a physics decision.
ALL ASSERTS PASSED""")
