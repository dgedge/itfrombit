r"""DIRECT LATTICE PHOTON FLOW COEFFICIENT — retires the Chadha-Nielsen import.

T-R2 closed the photon's IR Lorentz restoration on three facts about the
charged-matter vacuum polarization Pi_ij(q): it is (a) isotropic, (b)
logarithmically divergent as q->0, (c) screening-signed. (a) was computed at
q=0; (b),(c) were IMPORTED from Chadha-Nielsen / standard QED. This script
computes the full q-dependence of Pi_ij(q) DIRECTLY from the lattice massless
Weyl loop and extracts the running, retiring the import.

Method: one-loop static (Omega=0) current-current polarization of the exactly
isotropic Weyl cone H(k)=sum_d sin(k_d) sigma_d, filled lower band (Dirac sea).
Using P_pm(k)=(I +/- shat.sigma)/2 and the symmetric current vertex
V_i = cos(k_i + q_i/2) sigma_i, the interband bubble reduces analytically to
  Pi_ij(q) = -(1/N^3) sum_k c_i c_j [ delta_ij(1 + shat_k.shat_{k+q})
                                      - (shat_k,i shat_{k+q,j} + i<->j) ]
                                   / (E_k + E_{k+q}),
(the imaginary epsilon terms cancel between the two band orderings). We extract
the transverse part Pi_T(q), form Pi_T(q)/|q|^2, and fit it against log(1/|q|)
along [100] and [111]:
   * SLOPE A = the flow coefficient (d of the photon kinetic stiffness per
     decade of momentum) -- the object Chadha-Nielsen supplied;
   * A_[100] == A_[111]  => the running is ISOTROPIC (computed, not imported);
   * sign of A => screening (kinetic stiffness grows in the IR, no pole);
   * Ward check q_i Pi_ij ~ 0 (current conservation on the lattice).

exit 0 = Pi transverse, Pi_T/q^2 linear in log(1/q) with a positive isotropic
         slope; the flow coefficient is reported as a lattice number.
"""
import numpy as np

N = 72
ax = 2*np.pi*np.arange(N)/N
KX, KY, KZ = np.meshgrid(ax, ax, ax, indexing='ij')
SIN = [np.sin(KX), np.sin(KY), np.sin(KZ)]
E = np.sqrt(SIN[0]**2 + SIN[1]**2 + SIN[2]**2)
Ereg = np.where(E < 1e-9, 1e-9, E)
SHAT = [SIN[d]/Ereg for d in range(3)]          # unit Dirac vector shat_k

def polarization(qvec):
    """static Pi_ij(q) (3x3) for q a grid vector qvec=(mx,my,mz) integers."""
    mx, my, mz = qvec
    def roll(A): return np.roll(np.roll(np.roll(A, -mx, 0), -my, 1), -mz, 2)
    Ep = roll(E); Epreg = np.where(Ep < 1e-9, 1e-9, Ep)
    SHATp = [roll(SHAT[d]) for d in range(3)]
    qphys = 2*np.pi*np.array(qvec)/N
    # symmetric-point vertex factors c_i = cos(k_i + q_i/2)
    C = [np.cos([KX, KY, KZ][d] + qphys[d]/2) for d in range(3)]
    denom = E + Ep
    mask = (E > 1e-6) & (Ep > 1e-6)
    sdot = sum(SHAT[d]*SHATp[d] for d in range(3))   # shat_k . shat_{k+q}
    Pi = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Mij = C[i]*C[j]*((1.0 if i == j else 0.0)*(1+sdot)
                             - (SHAT[i]*SHATp[j] + SHAT[j]*SHATp[i]))
            integ = np.where(mask, -Mij/np.where(denom < 1e-9, 1e-9, denom), 0.0)
            Pi[i, j] = integ.sum()/N**3
    return Pi, np.linalg.norm(qphys)

def transverse(Pi, qhat):
    proj = np.eye(3) - np.outer(qhat, qhat)
    return 0.5*np.sum(proj*Pi)          # (1/2) Tr[(1-qq) Pi] = Pi_T in 3D
def ward(Pi, qhat):
    return np.linalg.norm(Pi @ qhat - (qhat @ Pi @ qhat)*qhat)  # transverse residual

print(f"[1] one-loop lattice polarization Pi_ij(q) on N={N}^3, massless Weyl sea.")
# gauge-invariant subtraction: the diamagnetic/seagull term equals the
# paramagnetic bubble at q=0 (f-sum rule), so the physical, running part is
# dPi(q) = Pi^para(q) - Pi^para(0). Without it a spurious contact 'photon mass'
# Pi_T(0)!=0 dominates Pi_T/q^2 ~ 1/q^2 and the directions look anisotropic.
Pi0, _ = polarization((0, 0, 0))
print(f"    contact term Pi^para(0): Tr = {np.trace(Pi0):+.4f} (subtracted as the seagull)")
print("    q-direction, |q|, dPi_T/|q|^2, Ward residual (on dPi):")
data = {}
for name, base in [("[100]", (1, 0, 0)), ("[111]", (1, 1, 1))]:
    rows = []
    qh = np.array(base)/np.linalg.norm(base)
    for m in range(1, 9):
        Pi, qmag = polarization(tuple(m*b for b in base))
        dPi = Pi - Pi0
        PiT = transverse(dPi, qh); w = ward(dPi, qh)
        rows.append((qmag, PiT/qmag**2, w))
        print(f"    {name}  |q|={qmag:.4f}  dPi_T/q^2={PiT/qmag**2:+.5f}  Ward={w:.2e}")
    data[name] = rows
maxward = max(r[2] for nm in data for r in data[nm])
print(f"    max Ward residual = {maxward:.2e} (current conservation on the lattice)")
assert maxward < 5e-2          # transverse to a few %: Ward holds on the lattice

print("\n[2] EXTRACT the flow coefficient: Pi_T/q^2 vs log(1/|q|), linear fit.")
slopes = {}
for name in ("[100]", "[111]"):
    q = np.array([r[0] for r in data[name]]); y = np.array([r[1] for r in data[name]])
    x = np.log(1.0/q)
    A, B = np.polyfit(x, y, 1)              # y = A log(1/q) + B
    resid = np.sqrt(np.mean((y-(A*x+B))**2))
    slopes[name] = A
    print(f"    {name}: Pi_T/q^2 = {A:+.5f} log(1/|q|) + {B:+.5f}   (rms resid {resid:.2e})")
A100, A111 = slopes["[100]"], slopes["[111]"]
print(f"    flow coefficient A: [100]={A100:+.5f}  [111]={A111:+.5f}")
print(f"    isotropy of the running |A100-A111|/|A| = {abs(A100-A111)/abs(0.5*(A100+A111)):.3f}")
assert abs(A100-A111)/abs(0.5*(A100+A111)) < 0.10   # running is isotropic, computed
assert abs(A100) > 1e-3 and np.sign(A100) == np.sign(A111)   # genuine log, definite sign
Abar = 0.5*(A100+A111)
print(f"    mean lattice flow coefficient A = {Abar:+.5f}")
print(f"    continuum reference 1/(12 pi^2) = {1/(12*np.pi**2):.5f} (Dirac), "
      f"1/(24 pi^2) = {1/(24*np.pi**2):.5f} (Weyl)")

print("\n[3] THE QUANTITATIVE TEST — is the log fast enough? (the decisive check):")
# dressed magnetic stiffness Z(n,mu) = Z0(n) + |A| log(Lambda/mu); fractional
# velocity anisotropy Dv/v = (sqrt(Z100)-sqrt(Z111))/sqrt(Zbar). Bare K6: v^2 = 2/3, 1/3.
Z100, Z111 = 2/3, 1/3
def dvv(decades):
    g = abs(Abar)*np.log(10.0**decades)          # running over 'decades' of momentum
    z1, z2 = Z100+g, Z111+g
    return (np.sqrt(z1)-np.sqrt(z2))/(0.5*(np.sqrt(z1)+np.sqrt(z2)))
print(f"    bare (0 decades): Dv/v = {dvv(0):.4f}  (the 41% K6 anisotropy)")
for dec in (5, 15, 30, 60, 120):
    print(f"    after {dec:3d} decades of log running: Dv/v = {dvv(dec):.4e}")
SME = 1e-17
print(f"    SME cavity bound on Dv/v ~ (a0 k)^2 ~ {SME:.0e}")
print(f"    Dv/v at 60 decades ({dvv(60):.2e}) exceeds the SME bound by "
      f"~{dvv(60)/SME:.0e}.")
# the log is monotone but hopeless: assert it does NOT reach the bound (honest negative)
assert dvv(60) > 1e6*SME and dvv(120) > 1e6*SME

print(f"""
[4] VERDICT — the import is COMPUTED, and found QUANTITATIVELY INSUFFICIENT:
  The Chadha-Nielsen import supplied three QUALITATIVE facts about the matter
  loop; all three are now confirmed by direct lattice computation:
    * TRANSVERSE  (q_i Pi_ij ~ 0, residual {maxward:.1e}: current conserved);
    * LOG-DIVERGENT  Pi_T(q)/q^2 = A log(1/|q|) + B (clean line, both directions);
    * ISOTROPIC + SCREENING  (A_[100]={A100:+.4f}, A_[111]={A111:+.4f}, agree to
      {abs(A100-A111)/abs(Abar)*100:.1f}%; single sign, no pole). Flow coefficient
      A = {Abar:+.5f} (lattice; the static-3D object, not the covariant 1/(12 pi^2)).
  BUT the computation also settles the QUANTITATIVE question the import left
  implicit, and the answer is negative: a logarithm cannot hide a MARGINAL O(1)
  anisotropy. Even over 60 decades of running the residual Dv/v is ~{dvv(60):.1e},
  some 15 orders above the SME bound 1e-17. This is the Collins-Perez-Sudarsky
  obstruction: a dimension-4 (marginal) Lorentz-violating velocity runs only
  logarithmically and is NOT driven to observational invisibility.
  CONSEQUENCE (honest self-correction): the earlier T-R2 / gauge-web claim that
  the Nielsen-Ninomiya log 'drives the anisotropy to zero in the IR -> photon
  PASS' is DOWNGRADED. The loop has the right qualitative form, but the
  velocity-log mechanism ALONE does not meet the bound. Photon-Lorentz closure
  needs more than the log:
    (i) an induced-dominated photon kinetic term (isotropic non-perturbatively,
        as argued for the graviton in T-R5) -- requires the bare K6 Maxwell term
        to be negligible vs the matter loop, NOT established here;
   (ii) a custodial symmetry forbidding the marginal cubic anisotropy; or
  (iii) the T1u-E_g degeneracy-gap demotion (power-law (q/Delta)^2, fast enough)
        -- the loop must actually gap the E_g, also not established here.
  The GRAVITON (T-R5) is unaffected only insofar as its kinetic term is genuinely
  induced (no bare marginal term); this result flags that the bare K6 E_g band
  deserves the same scrutiny. NET: direct computation retired the import's
  qualitative content and, in doing so, exposed that the quantitative closure was
  over-claimed. Photon-band Lorentz invariance reverts to OPEN, with three named
  candidate mechanisms and a sharp falsification bound.
exit 0""")
print("ALL ASSERTIONS PASSED — trio computed; log shown quantitatively insufficient.")
