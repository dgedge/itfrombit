#!/usr/bin/env python3
"""
SU(6) spin-flavour baryon magnetic moments, built on the item-138 E_1/2 spin-coin.

With spin now living in a SEPARATE C^2 E_1/2 coin (spin_coin_2O.py) rather than
tangled in the chi bit, the spin (x) flavour SU(6) wavefunction can be formed
explicitly. The quark magnetic-moment operator (g=2, from the Dirac coin) is
    mu_z / mu_0 = sum_i Q_i sigma_z^(i),   Q_u = +2/3, Q_d = -1/3  (framework charges, §2.8)
sigma_z acting on the E_1/2 coin, Q_i the fractional charge of quark i.

This script (a) builds the totally-symmetric SU(6) proton (uud) and neutron (udd)
spin-1/2, S_z=+1/2 states explicitly, (b) verifies total symmetry / spin-1/2,
(c) computes mu_p, mu_n and the ratio, (d) compares to experiment.

HONEST SCOPE (asserted in the verdict): mu_p/mu_n = -3/2 is the STANDARD SU(6) /
constituent-quark-model result (textbook), reproduced here via the framework's
DERIVED charges + the new spin-coin -- not a framework-unique prediction. It
requires (i) equal u,d masses (mass cancels in the ratio; the 2.7% gap to data is
the m_u!=m_d + relativistic correction, not improvable at CQM level), and (ii) the
totally-symmetric SU(6) wavefunction -- which the framework's spin-statistics
(item 139) would be needed to DERIVE, and item 139 is NOT closed (ANCHOR §15 139
anchors it as a derivation-with-prerequisites; the E_1/2 coin closed only the
spin-ROTATION lock, not the exchange<->rotation FR step). So the SU(6) symmetry is
here an INPUT (standard CQM), not yet a framework OUTPUT.

numpy. exit 0 == mu_p/mu_n = -3/2 exactly, on the explicitly-symmetric wavefunctions.
"""
import numpy as np, itertools

# single-quark basis: 0=u^, 1=u_, 2=d^, 3=d_   (^=up, _=down)
Q  = np.array([ 2/3,  2/3, -1/3, -1/3])     # fractional charge by single-quark index
Sz = np.array([  +1,   -1,   +1,   -1])     # sigma_z eigenvalue (+1 up, -1 down)
def fla(i): return 'u' if i<2 else 'd'

def ket(i,j,k):                              # 3-quark basis vector in C^64
    v=np.zeros(64); v[16*i+4*j+k]=1.0; return v
def idx3(n): return (n//16, (n//4)%4, n%4)

# magnetic-moment operator mu_z/mu_0 = sum_i Q_i sigma_z^(i)  (diagonal in this basis)
mu_diag=np.array([Q[a]*Sz[a]+Q[b]*Sz[b]+Q[c]*Sz[c]
                  for a in range(4) for b in range(4) for c in range(4)])
def mu(psi): return float(np.real(np.vdot(psi, mu_diag*psi))/np.vdot(psi,psi))

# permutation operators on the 3 tensor factors (to test/impose total symmetry)
def Pperm(perm):
    M=np.zeros((64,64))
    for n in range(64):
        a,b,c=idx3(n); t=(a,b,c)
        ap,bp,cp=t[perm[0]],t[perm[1]],t[perm[2]]
        M[16*ap+4*bp+cp, n]=1.0
    return M
perms=list(itertools.permutations(range(3)))
Sym=sum(Pperm(p) for p in perms)/6.0        # total symmetriser

# --- proton seed: uu(slots1,2) coupled to spin-1, + d(slot3), total spin 1/2, Sz=+1/2 ---
# |1/2,1/2> = sqrt(2/3)|1,1>_uu |d_> - sqrt(1/3)|1,0>_uu |d^>
r23, r16 = np.sqrt(2/3), np.sqrt(1/6)
# u^=0,u_=1 ; d^=2,d_=3
phi_p = r23*ket(0,0,3) - r16*(ket(0,1,2)+ket(1,0,2))     # standard textbook proton seed
phi_n = r23*ket(2,2,1) - r16*(ket(2,3,0)+ket(3,2,0))     # neutron: u<->d (udd, odd quark = u)

def build(name, seed):
    seed=seed/np.linalg.norm(seed)
    sym=Sym@seed; sym=sym/np.linalg.norm(sym)            # totally symmetric SU(6) state
    tot_sym=all(np.allclose(Pperm(p)@sym, sym) for p in perms)
    mu_seed, mu_sym = mu(seed), mu(sym)
    print(f"  {name}: totally symmetric (SU(6) 56)? {tot_sym} ;  "
          f"mu(seed)={mu_seed:+.4f}  mu(symmetric)={mu_sym:+.4f}  (must agree)")
    assert tot_sym, "symmetrised state not totally symmetric"
    assert abs(mu_seed-mu_sym)<1e-9, "seed and symmetric state disagree on mu"
    return mu_sym

print("Building totally-symmetric SU(6) spin-flavour wavefunctions (spin on E_1/2 coin):")
mp = build("proton  (uud)", phi_p)
mn = build("neutron (udd)", phi_n)

print(f"\n  mu_p = {mp:+.5f} mu_0      (SU(6): 1/3(4*mu_u - mu_d) = +1)")
print(f"  mu_n = {mn:+.5f} mu_0      (SU(6): 1/3(4*mu_d - mu_u) = -2/3)")
print(f"  mu_p/mu_n = {mp/mn:+.5f}")
print(f"  experiment: mu_p/mu_n = 2.79285/(-1.91304) = {2.79285/-1.91304:+.5f}  "
      f"(dev {abs(mp/mn - 2.79285/-1.91304)/abs(2.79285/-1.91304)*100:.1f}%)")

assert abs(mp - 1.0) < 1e-9 and abs(mn + 2/3) < 1e-9, "SU(6) moments wrong"
assert abs(mp/mn - (-1.5)) < 1e-9, "mu_p/mu_n != -3/2"

# absolute moments need a constituent mass scale m_q (a FITTED Target-B input)
print(f"\n  Absolute values need a constituent mass m_q (fitted): mu[mu_N] = mu * (m_p/m_q).")
mp_mass = 938.272
for mq in (336.0,):   # the standard constituent value that makes mu_p ~ 2.79
    scale = mp_mass/mq
    print(f"    m_q={mq:.0f} MeV -> mu_p={mp*scale:+.3f} mu_N (exp 2.793), "
          f"mu_n={mn*scale:+.3f} mu_N (exp -1.913)  [m_q chosen to match mu_p; ratio is m_q-free]")

print("\n"+"="*72); print("VERDICT"); print("="*72)
print(" CONSTRUCTED & verified (exit 0): the totally-symmetric SU(6) p,n wavefunctions on")
print(" the E_1/2 coin + framework charges (2/3,-1/3) at g=2 give mu_p/mu_n = -3/2 EXACTLY,")
print(" matching experiment (-1.46) to 2.7%. The E_1/2 coin genuinely ENABLED this (a")
print(" separate spin space to tensor with flavour, which the chi-tangled spin could not).")
print(" BUT, honestly:")
print("  * mu_p/mu_n = -3/2 is the STANDARD SU(6)/constituent-quark-model result (textbook),")
print("    reproduced via the framework's derived charges + coin -- NOT a new/unique prediction.")
print("    Same pattern as last turn's Omega_b* ('the baryon sector IS the constituent quark")
print("    model') and the Koide/LFU re-derivations.")
print("  * The ratio is parameter-free ONLY given the standard CQM assumptions: equal u,d")
print("    masses (the 2.7% gap = m_u!=m_d + relativistic, not improvable at this level) and")
print("    the totally-symmetric SU(6) wavefunction.")
print("  * That SU(6) symmetry is an INPUT here, not a framework OUTPUT: deriving it needs the")
print("    spin-STATISTICS theorem (item 139), which ANCHOR §15 139 anchors as a DERIVATION-")
print("    WITH-PREREQUISITES, NOT 'closed' -- the FR exchange argument is exhibited for ONE")
print("    path; promotion needs prereqs (ii)/(iii), which spin_statistics_path_independence.py")
print("    shows need a topological framing (naive lattice Wilson line is geometry-dependent),")
print("    not a finite enumeration. The item-138 coin closed only the spin-ROTATION lock, so a")
print("    'newly closed spin-statistics theorem -> SU(6) forced' framing OVERSTATES the status.")
print("  * Absolute moments need the fitted constituent mass m_q; only the ratio is m_q-free.")
print("\nexit 0 == SU(6) wavefunctions explicitly symmetric; mu_p/mu_n = -3/2 verified.")
