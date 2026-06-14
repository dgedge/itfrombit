#!/usr/bin/env python3
r"""ITEM 7 -- the induced-gravity stiffness K_Eg, a single PRE-PINNED-SCHEME heat-kernel shot.
(DRIFT G7 follow-up: 'build W_QQ(k), BZ-average the E_g-projected stiffness, scheme fixed a priori'.)

§10.3:  K_Eg = d^2 H / d eps^2 |_Eg  ->  M_{P,bare} = sqrt(4 pi K_Eg) * Lambda ~ O(1) GeV   [open comp. target]
§10.4:  the framework ITSELF forbids bare-KINETIC strain ('anti-mass'); the E_g shear must come from the
        t_mix colour-confinement tension on the 3 bipyramid orientations (= 3 QCD colours). That is the
        physically-pinned operator. NO alpha appears in §10.3/10.4 -- the stiffness is a pure spectral number.

================================ PRE-REGISTERED PINS (fixed BEFORE any result) =========================
P1  Operator: Hermitian Bloch Hamiltonian (NOT log of the unitary walk -- avoids the §3.1-unitary vs
    §7.11-Hermitian ambiguity by taking the manifestly-Hermitian tight-binding generator of §3.2).
P2  Internal space: the 3 colour/bipyramid axes (x,y,z) -- the sector §10.4 names as the E_g-shear carrier.
P3  Kinetic part: body-diagonal hopping to the 8 corner-neighbours, amplitude t = 1/3 (Grover-uniqueness
    substrate constant, ANCHOR L1556), energies in units of Lambda (so K_Eg is dimensionless, as §10.3 needs).
P4  Confinement tension: t_mix * A_{K3} on the 3 colours (V_strong, §3.3; A_{K3}=K3 adjacency). Sign pinned
    so the colour SINGLET (1,1,1)/sqrt3 is the ground state (§10.4 'forced into equal superposition'): t_mix<0.
    Magnitude pinned to the §10.4 hierarchy t_mix >> t_hop; canonical g_s=1 (§3.2) -> |t_mix| scanned 1..8.
P5  E_g strain: traceless diagonal eps_i = eps*(1,1,-2)/sqrt6 (one of the two E_g components of O_h), applied
    as on-site energies on the 3 colour axes (§10.4: 'uniaxial squash deforms one bipyramid differently').
P6  Vacuum functional: zero-point energy u(eps) = (1/2) <|E|> over the 3 bands, BZ-averaged (Sakharov induced
    graviton zero-point). K_Eg = d^2 u/d eps^2 by 5-point finite difference at eps=0. (Filled-Dirac-sea variant
    reported too as a robustness check.) BZ = cubic [-pi,pi]^3, grid N pinned at 48.
P7  E_g projection: the strain direction IS the E_g channel (P5); no extra projector needed on the 3-axis space
    (the traceless-diagonal subspace of the 3 colour axes is exactly the E_g doublet).
=======================================================================================================
The shot's PASS criterion (declared in advance): K_Eg is O(1) dimensionless AND M_{P,bare}=sqrt(4pi K_Eg)Lambda
is O(1) GeV (confirms §10.3) -- AND we check whether alpha^2 can emerge (it cannot if alpha is absent from the
operator). The shot RESCUES the parameter-free gravity claim ONLY if alpha^2 falls out of the stiffness.
"""
import sys, itertools as it
import numpy as np

LAM   = 0.332                    # GeV, Lambda_QCD (energy unit; set to 1 inside, restored at the end)
ALPHA = 1/137.035999
t_hop = 1.0/3.0                  # Grover substrate hopping amplitude (P3)
N     = 48                       # BZ grid per axis (P6)

# --- P3: body-diagonal kinetic dispersion eta(k) (8 corner neighbours, energies in units of Lambda) ---
dirs = np.array([[(2*((f>>2)&1)-1), (2*((f>>1)&1)-1), (2*(f&1)-1)] for f in range(8)], float)/np.sqrt(3)
kax  = (np.arange(N)+0.5)/N*2*np.pi - np.pi          # k in [-pi,pi)
KX,KY,KZ = np.meshgrid(kax,kax,kax, indexing='ij')
Kvec = np.stack([KX,KY,KZ], -1).reshape(-1,3)        # (N^3,3)
# eta(k) = -t * sum_f cos(k . d_f)   (real, Hermitian kinetic)
eta = -t_hop * np.cos(Kvec @ dirs.T).sum(1)          # (N^3,)
print(f"[build] body-diagonal kinetic band eta(k): range [{eta.min():.3f},{eta.max():.3f}] * Lambda ; "
      f"grid {N}^3={len(eta)} k-points ; t={t_hop:.4f}")

A_K3 = np.array([[0.,1,1],[1,0,1],[1,1,0]])          # K3 colour adjacency (P4)
EVEC = np.array([1.,1,-2])/np.sqrt(6)                # E_g strain direction (P5), traceless

def colour_mu(eps, t_mix):
    """3 colour eigenvalues of  diag(eps*EVEC) + t_mix*A_K3  (k-independent internal part)."""
    H3 = np.diag(eps*EVEC) + t_mix*A_K3
    return np.linalg.eigvalsh(H3)                     # ascending

def vacuum_u(eps, t_mix, dirac_sea=False):
    """zero-point (P6) [or filled-Dirac-sea] vacuum energy per cell, BZ-averaged, in units of Lambda."""
    mu = colour_mu(eps, t_mix)                        # (3,)
    E  = eta[:,None] + mu[None,:]                     # (N^3, 3) band energies
    if dirac_sea:
        return np.where(E < 0, E, 0.0).sum(1).mean()  # sum of filled (negative) modes
    return 0.5*np.abs(E).sum(1).mean()                # (1/2) sum |E|

def stiffness(t_mix, h=1e-3, dirac_sea=False):
    """K_Eg = d^2 u / d eps^2 at eps=0, 5-point central difference (P6)."""
    f = lambda e: vacuum_u(e, t_mix, dirac_sea)
    return (-f(2*h) + 16*f(h) - 30*f(0) + 16*f(-h) - f(-2*h)) / (12*h*h)

print("\n[shot] K_Eg = d^2 u/d eps^2 (zero-point), and M_{P,bare} = sqrt(4 pi |K_Eg|) * Lambda, vs |t_mix|:")
print(f"  {'|t_mix|':>7} {'K_Eg(zp)':>10} {'K_Eg(sea)':>10} {'M_P,bare[GeV]':>14}")
rows=[]
for tm in [1,2,3,4,6,8]:
    Kzp  = stiffness(-tm, dirac_sea=False)
    Ksea = stiffness(-tm, dirac_sea=True)
    Mbare= np.sqrt(4*np.pi*abs(Kzp))*LAM
    rows.append((tm,Kzp,Ksea,Mbare))
    print(f"  {tm:>7} {Kzp:>10.4f} {Ksea:>10.4f} {Mbare:>14.4f}")

# canonical pin g_s=1 -> |t_mix|=1; report it as THE shot
Kc = stiffness(-1.0); Mc = np.sqrt(4*np.pi*abs(Kc))*LAM
print(f"\n[canonical pin |t_mix|=g_s=1]  K_Eg = {Kc:.4f} (dimensionless)  ->  M_P,bare = {Mc:.4f} GeV")

# --- the decisive test: does alpha^2 fall out? alpha is ABSENT from every pinned ingredient. ---
print(f"""
[alpha test] The operator's ingredients are: t=1/3 (Grover), A_K3 (K3 graph), E_g strain direction,
  body-diagonal geometry, Lambda (unit). The fine-structure constant alpha=1/{1/ALPHA:.1f} appears in NONE of
  them. K_Eg is therefore a pure spectral number O(1); it sets the SCALE M_P,bare ~ Lambda and is
  STRUCTURALLY INCAPABLE of producing the alpha^2 (=  {ALPHA**2:.2e}) that the macroscopic formula needs.
  alpha^2 enters only at the *dilution* step (§10.5 K_eff / Part-20 'double-Landauer'), which is the
  asserted part (DRIFT G7). The heat-kernel does not -- and cannot -- derive it.""")

bare_ok = 0.2 < Mc < 5.0                              # 'O(1) GeV' per §10.3
# the 19-OOM bridge: M_P,macro = M_bare * sqrt( (R_dS/a0) / K ).  The large factor is the HORIZON/lattice ratio.
R_dS_over_a0 = 2.77e41                                # §10.5 node count = R_dS/a0 (the cosmological INPUT)
print(f"\n================================ VERDICT (pre-pinned heat-kernel shot, item 7) ====================")
print(f"  (1) BUILT: the manifestly-Hermitian E_g-shear stiffness operator on the 3-colour/bipyramid sector")
print(f"      (§10.4 confinement-tension carrier) WAS constructible under the pinned scheme -> §10.3's 'open")
print(f"      computational target' is reachable for the BARE STIFFNESS. K_Eg = {Kc:.2f} (O(1) dimensionless),")
print(f"      M_P,bare = {Mc:.2f} GeV = O(1) GeV [{'PASS' if bare_ok else 'FAIL'}] -- this confirms the LAMBDA-END of the bridge.")
print(f"  (2) BUT this is only the LAMBDA-END. M_P,macro = M_bare * sqrt((R_dS/a0)/K): the 19-OOM jump from")
print(f"      {Mc:.2f} GeV to 1.2e19 GeV is the factor sqrt({R_dS_over_a0:.1e}) ~ {np.sqrt(R_dS_over_a0):.1e}, i.e. the")
print(f"      HORIZON/lattice ratio -- supplied as a COSMOLOGICAL INPUT, not from any framework-intrinsic scale.")
print(f"      The heat-kernel computes the one intrinsic piece (~Lambda); the magnitude of M_P is horizon-eaten.")
print(f"  (3) So all three routes are ONE Dirac large-number relation (Lambda_QCD <-> horizon), dressed by three")
print(f"      prefactors (205 / 24pi alpha^2 / alpha^-4 (L_H/a0)^1/4). The CANON ADMITS this: §10.7 L2350 calls")
print(f"      the Part-20 M_P formula a Dirac-LNH 'structural identity' (M_P^2/m_p^2 ~ R_H/r_p); L2647 tiers")
print(f"      item 136 as an 'alpha^-4 large-number coincidence'. One coincidence admits many prefactors -> 3 routes.")
print(f"  (4) alpha is absent from the stiffness, so alpha^2 CANNOT emerge here; the alpha-power dresses the")
print(f"      horizon-consuming DILUTION, not the intrinsic scale. Even a FUTURE derived alpha-power would buy")
print(f"      only a derived Dirac RELATION among {{M_P, Lambda_QCD, horizon}} -- NOT an intrinsic M_P prediction,")
print(f"      because the formula still consumes the horizon. That is the true (weaker) ceiling on W_QQ(k).")
print(f"  NET: the shot does what it can -- the Lambda-end stiffness is now a genuine (pin-dependent, ~3x) O(1)")
print(f"  computation, not an assertion -- and clarifies what it CANNOT do: M_P is the Lambda<->horizon large-")
print(f"  number relation with the horizon as input and the prefactor asserted. The honest residual is thinner")
print(f"  than 'magnitude pinned': the magnitude is pinned by the COSMOLOGICAL INPUT, not the framework's scales.")
print("=" * 100)
assert bare_ok, "bare Planck mass should land O(1) GeV"
assert isinstance(Kc, float)
print("exit 0 -- Lambda-end stiffness built & computed O(1) GeV; M_P magnitude is horizon-supplied; alpha^2 not derivable.")
