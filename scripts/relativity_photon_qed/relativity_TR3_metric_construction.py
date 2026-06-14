r"""T-R3 — THE METRIC VARIABLE as coarse strain/clock covariance (constructive).

The matter sector is now an explicit, exactly-isotropic massless Weyl cone
H(k)=sum_d sin(k_d)Gamma_d with H^2=(sum sin^2 k_d) I (this morning). T-R3 asks:
what IS the metric, constructively, and where does the E_g graviton live?

Construction: strain the substrate. A symmetric strain eps_ij distorts bonds,
so the hop momentum maps k -> (I+eps)k and the matter dispersion becomes
E^2(k) = |(I+eps)k|^2 = k^T (I+eps)^2 k.  Reading E^2 = g^{ij} k_i k_j defines
the effective INVERSE METRIC g^{ij}(eps) = [(I+eps)^2]_{ij} = delta + 2 eps + ...
So the metric perturbation IS the strain field, and matter follows geodesics
because its propagation phase is extremised in this distorted clock field.

The symmetric strain (= metric perturbation h_ij) decomposes under O_h as
  A_1g (trace, 1)  +  E_g (2)  +  T_2g (3),
with the spin-2 GRAVITON = E_g (+) T_2g (the full l=2, 5-dim) and the A_1g
trace = the conformal factor / lapse (the clock rate). Matter couples to the
metric through the rank-2 STRESS vertex dH/d eps (not the U(1) current), and
the isotropic Weyl loop dresses all 5 spin-2 modes EQUALLY at leading order
(the rank-2 analogue of the photon's delta_ij lemma) -> spin-2 degeneracy,
the graviton's IR-Lorentz restoration (conditional on the stress-loop scaling,
the T-R5 residual).

exit 0 = metric = delta+2eps verified; irrep split A1g+Eg+T2g confirmed; stress
         vertex rank-2; the 5 spin-2 modes dressed equally at leading order.
"""
import numpy as np

sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.array([[1,0],[0,-1]],complex)
S=[sx,sy,sz]
def Hk_strained(k, eps):
    kp = (np.eye(3)+eps) @ np.array(k)
    return sum(np.sin(kp[d])*S[d] for d in range(3))
def E2(k, eps):
    return float(np.max(np.linalg.eigvalsh(Hk_strained(k, eps)))**2)

print("[1] CONSTRUCTIVE METRIC: strain the Weyl cone, read g^{ij} off the dispersion.")
rng_eps = np.array([[0.03,0.01,-0.02],[0.01,-0.01,0.015],[-0.02,0.015,-0.02]])  # symmetric, ~traceless
rng_eps = (rng_eps+rng_eps.T)/2
h = 1e-4
def hessian_metric(eps):
    g = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            ei = np.zeros(3); ei[i]=h; ej = np.zeros(3); ej[j]=h
            g[i,j] = (E2(ei+ej,eps)-E2(ei-ej,eps)-E2(-ei+ej,eps)+E2(-ei-ej,eps))/(4*h*h)/2
    return g
g_flat = hessian_metric(np.zeros((3,3)))
g_eps  = hessian_metric(rng_eps)
print(f"    flat (eps=0):   g^ij = delta?  ||g-I|| = {np.linalg.norm(g_flat-np.eye(3)):.2e}")
assert np.linalg.norm(g_flat-np.eye(3)) < 1e-6
pred = np.eye(3)+2*rng_eps
print(f"    strained:       g^ij vs delta+2eps  ||g-(I+2eps)|| = {np.linalg.norm(g_eps-pred):.2e}")
assert np.linalg.norm(g_eps-pred) < 5e-3        # linear identification g = delta + 2 eps
print("    -> the metric perturbation h^ij = 2 eps^ij IS the strain field (DERIVED).")

print("\n[2] IRREP CONTENT of the symmetric metric perturbation under O_h:")
# basis: A1g trace; E_g (2); T_2g (3)
A1g = [np.eye(3)/np.sqrt(3)]
Eg  = [np.diag([1,-1,0])/np.sqrt(2), np.diag([1,1,-2])/np.sqrt(6)]
def offsym(i,j):
    M=np.zeros((3,3)); M[i,j]=M[j,i]=1/np.sqrt(2); return M
T2g = [offsym(0,1), offsym(1,2), offsym(0,2)]
basis = A1g+Eg+T2g
# orthonormal & complete on symmetric 3x3 (6-dim)
G = np.array([[np.trace(a.T@b) for b in basis] for a in basis])
print(f"    dims: A1g={len(A1g)}, E_g={len(Eg)}, T_2g={len(T2g)}  (total {len(basis)} = sym 3x3)")
print(f"    orthonormal basis? ||Gram - I|| = {np.linalg.norm(G-np.eye(6)):.2e}")
assert len(basis)==6 and np.linalg.norm(G-np.eye(6))<1e-12
print("    -> metric perturbation = A_1g (conformal/lapse = the CLOCK) (+) E_g (+) T_2g.")
print("       SPIN-2 GRAVITON = E_g (+) T_2g (the full l=2); §10.1's 'E_g graviton' is")
print("       its cubic-visible half — full Lorentz spin-2 needs E_g, T_2g degenerate.")

print("\n[3] CLOCK covariance: the A_1g trace is the local light-cone scale.")
for s in (0.0, 0.05, 0.10):
    eps = s*np.eye(3)                                # pure conformal strain
    v = np.sqrt(E2((0.01,0,0),eps))/0.01            # local speed along x
    print(f"    conformal strain s={s:.2f}:  c_local = {v:.4f}  (= 1+s = {1+s:.4f})")
    assert abs(v-(1+s)) < 1e-3
print("    -> a uniform (A_1g) strain rescales the light cone c_local = 1+s: the metric")
print("       as 'coarse clock covariance' is literal — strain renormalises the tick.")

print("\n[4] STRESS COUPLING + spin-2 dressing (rank-2 analogue of the photon):")
# stress vertex S_ij(k)=cos(k_i)k_j sigma_i; mode coupling O_a=sum_i cos(k_i)(T_a k)_i sigma_i
# leading-order identity: <-,k|O_a|-,k> = -|k| (khat^T T_a khat) -> stiffness density
#   |k|^2 (khat^T T_a khat)^2, whose angular integral = (2/15)Tr(T_a^2): EQUAL for all
#   traceless T_a. Verify the identity, then show E_g/T_2g -> equal as the box shrinks.
def lead_check(Ta, k):
    H=sum(np.sin(k[d])*S[d] for d in range(3)); ev,V=np.linalg.eigh(H)
    lo=V[:,0]; O=sum(np.cos(k[i])*(Ta@k)[i]*S[i] for i in range(3))
    return np.real(lo.conj()@O@lo)
kt=np.array([0.02,0.013,-0.017]); kh=kt/np.linalg.norm(kt)
for Ta in (Eg[0],T2g[0]):
    pred=-np.linalg.norm(kt)*(kh@Ta@kh); got=lead_check(Ta,kt)
    assert abs(pred-got)<1e-4
print("    leading-order identity <-|O_a|-> = -|k| khat^T T_a khat: VERIFIED (1e-4).")
def stiffness(Ta, K):
    N=18; pts=[(-K+2*K*(i+0.5)/N) for i in range(N)]; tot=0.0; cnt=0
    for kx in pts:
        for ky in pts:
            for kz in pts:
                k=np.array([kx,ky,kz])
                if np.linalg.norm(k) > K or np.linalg.norm(k)<1e-6: continue
                O=sum(np.cos(k[i])*(Ta@k)[i]*S[i] for i in range(3))
                H=sum(np.sin(k[d])*S[d] for d in range(3))
                ev,V=np.linalg.eigh(H); Pm=np.outer(V[:,0],V[:,0].conj())
                tot+=np.real(np.trace(O@Pm@O@Pm)); cnt+=1
    return tot/cnt
print("    E_g vs T_2g stiffness as the momentum box K shrinks toward the IR:")
for K in (1.2, 0.6, 0.3, 0.15):
    se=np.mean([stiffness(T,K) for T in Eg]); st=np.mean([stiffness(T,K) for T in T2g])
    print(f"      K={K:<4}:  E_g={se:.6f}  T_2g={st:.6f}  ratio={st/se:.4f}")
rE=np.mean([stiffness(T,0.15) for T in Eg]); rT=np.mean([stiffness(T,0.15) for T in T2g])
print(f"    leading-order (K=0.15) ratio T_2g/E_g = {rT/rE:.4f} -> 1 (isotropic)")
assert abs(rT/rE - 1.0) < 0.05
print("    -> the isotropic Weyl loop dresses all 5 spin-2 modes EQUALLY at leading")
print("       order (SO(3): only invariant is Tr(T^2)); the BZ-scale split (K~1) is the")
print("       cubic higher-order term — the rank-2 analogue of the photon's delta_ij")
print("       lemma, with the same leading-isotropic / cubic-suppressed structure.")

print(f"""
[5] T-R3 VERDICT — constructive metric DELIVERED:
  * METRIC (DERIVED): strain the Weyl cone -> g^{{ij}} = delta^{{ij}} + 2 eps^{{ij}};
    the metric perturbation IS the substrate strain field, and matter geodesics
    are propagation-phase extremals in the distorted clock field. The 'metric as
    coarse strain/clock covariance' is now literal and computed.
  * DECOMPOSITION (DERIVED): h_ij = A_1g (conformal/lapse = the CLOCK) (+) E_g (+)
    T_2g; the SPIN-2 GRAVITON is the full l=2 = E_g (+) T_2g. §10.1's E_g graviton
    is its cubic-visible half; full Lorentz spin-2 requires E_g, T_2g degenerate.
  * CLOCK (DERIVED): a uniform A_1g strain rescales the light cone c_local = 1+s
    — strain renormalises the tick, exactly the SR-clock picture under load.
  * STRESS COUPLING (DERIVED): matter couples to the metric via the rank-2 stress
    vertex dH/d eps (NOT the U(1) current), and the isotropic Weyl loop dresses
    all 5 spin-2 modes EQUALLY at leading order (rank-2 analogue of the delta_ij
    photon lemma) -> E_g/T_2g degeneracy at the induced level.
  * RESIDUAL (honest): whether that equal, growing dressing DOMINATES the bare
    cubic E_g/T_2g splitting in the IR (graviton Lorentz restoration) depends on
    the stress-tensor loop's SCALING — gravity is non-renormalisable, so unlike
    the photon's NN log this is NOT a standard import. That scaling is the
    gravitational-coupling question = T-R5 (Einstein eq / the §10.5 magnitude).
    T-R3 delivers the metric variable and the isotropic-dressing lemma; the
    graviton's full IR closure is handed to T-R5.
exit 0""")
print("ALL ASSERTIONS PASSED — metric=delta+2eps, irreps, clock, equal spin-2 dressing.")
