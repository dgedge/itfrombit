#!/usr/bin/env python3
"""
Non-abelian analog of the user's U(1) Gauss-projection test: does the non-abelian
gauge-sector gap-closing SURVIVE Gauss projection (genuine physical degeneracy) or LIFT
(gauge copy), and is the gauge-invariant spectrum gapped?

Background: the U(1) Wilson proxy's inter-sector degeneracies turned out to be gauge copies
(same-flux sectors exactly isospectral) and the flux-fixed gap was positive. The genuine hard
case is NON-abelian (my gap_under_twist found SU(2) T1 / SU(3) background twists close the
CSS gap). This script does the SU(2) version PROPERLY: dynamical gauge + Gauss's law.

Method (standard single-plaquette Kogut-Susskind reduction):
  For one plaquette of pure SU(2) gauge theory, after imposing Gauss's law at every vertex the
  GAUGE-INVARIANT (physical) Hilbert space is spanned by the single Wilson-loop rep
  |j>, j = 0, 1/2, 1, 3/2, ...  The KS Hamiltonian on this physical ladder is:
     H_E |j> = (g^2/2) * c_E * j(j+1) |j>           (electric energy = quadratic Casimir)
     H_B     = -beta * (Tr_F U_plaq)                 hops j -> j +/- 1/2  (SU(2) fusion
               chi_{1/2} chi_j = chi_{j-1/2} + chi_{j+1/2}, coefficient 1)
  => H is a tridiagonal 'discrete Schrodinger operator on the j-ladder'. This IS the
  Gauss-projected physical Hamiltonian (standard result; e.g. the single-plaquette problem).
  The (2j+1)^2 link states behind each |j> are the GAUGE-ORBIT copies removed by projection.

What this decides:
  * un-projected: each physical level has (2j+1)^2-fold multiplicity = the non-abelian gauge
    crossings -> these LIFT under Gauss projection (they are gauge copies, like the U(1) case);
  * physical (Gauss-projected) spectrum: is it gapped, across the coupling?

HONEST SCOPE: this is PURE GAUGE (no matter). It settles the GAUGE-SECTOR question (are the
non-abelian crossings gauge copies / is the physical gauge spectrum gapped). It does NOT test
the MATTER/mirror SMG gap under non-abelian gauging -- that is the model-dependent item-13
construction and is NOT settled here. numpy; self-asserting on construction, reporting the gap.
"""
import numpy as np
from fractions import Fraction as Fr
def Hd(t): print("\n"+"="*78+"\n"+t+"\n"+"="*78)

# physical (Gauss-projected) j-ladder for SU(2) one plaquette
def jladder(beta, cE=1.0, twojmax=40, start_twoj=0):
    # index n = 2j, n = start_twoj, start_twoj+2? no: j steps by 1/2 -> 2j steps by 1 -> n=0,1,2,...
    ns=list(range(start_twoj, twojmax+1))
    d=len(ns); H=np.zeros((d,d))
    for a,n in enumerate(ns):
        j=n/2.0
        H[a,a]=cE*j*(j+1)                      # electric Casimir
        if a+1<d: H[a,a+1]=H[a+1,a]=-beta      # magnetic: j <-> j+1/2 (fusion coeff 1)
    return H, ns
def gap_of(H):
    w=np.linalg.eigvalsh(H); return float(w[1]-w[0]), w[:4]

# ----------------------------------------------------------------------------------
Hd("(A) Gauge-orbit multiplicity: the non-abelian 'crossings' are (2j+1)^2 gauge copies")
print("  In the FULL link Hilbert space each loop rep |J=j> carries (2j+1)^2 states (the two")
print("  magnetic indices m,n of the Wilson loop). Gauss's law (vertex gauge invariance) keeps")
print("  exactly ONE per j (the singlet). So:")
for twoj in range(0,5):
    j=Fr(twoj,2); mult=(twoj+1)**2
    print(f"     j={str(j):>3}: full-link multiplicity = (2j+1)^2 = {mult:>2}  -> Gauss-projected: 1")
print("  => the un-projected degeneracies between these are EXACTLY gauge copies (same physical")
print("     state), and they LIFT under Gauss projection -- the non-abelian analog of the U(1)")
print("     same-P_X isospectrality. (Guaranteed: gauge transformations commute with H.)")

# ----------------------------------------------------------------------------------
Hd("(B) Physical (Gauss-projected) gauge spectrum: is it gapped across the coupling?")
print(f"  {'beta':>6} | {'gap (E1-E0)':>12} | {'low levels (E - E0)':>34}")
betas=[0.0,0.25,0.5,1.0,2.0,5.0,10.0,25.0,100.0]
mingap=1e9
for b in betas:
    H,ns=jladder(b)
    g,levels=gap_of(H); mingap=min(mingap,g)
    rel=np.round(levels-levels[0],3)
    print(f"  {b:>6} | {g:>12.4f} | {str(rel):>34}")
assert mingap>1e-3
print(f"\n  min physical gap over the scan = {mingap:.4f}  -> GAPPED throughout.")
print("  -> the Gauss-projected non-abelian gauge spectrum stays GAPPED; the crossings the raw")
print("     (un-projected) proxy showed were gauge copies that the projection removes.")

# trend note (the 'shrinking gap' the U(1) scan also showed)
g_small,_=gap_of(jladder(0.5)[0]); g_large,_=gap_of(jladder(100.0)[0])
print(f"  trend: gap(beta=0.5)={g_small:.3f} -> gap(beta=100)={g_large:.3f}  -> "
      f"{'shrinks' if g_large<g_small else 'GROWS'} with beta here, UNLIKE the U(1)+matter proxy")
print( "  (which shrank): the beta-trend is model-dependent (likely the matter coupling); pure")
print( "  SU(2) gauge shows no gap-closing across this scan.")

# ----------------------------------------------------------------------------------
Hd("(C) Robustness: a static fundamental charge (half-integer ladder) is also gapped")
# a static charge in the fundamental forces the screening flux -> physical ladder starts at j=1/2
H_q,_=jladder(1.0, start_twoj=1)
gq,_=gap_of(H_q)
print(f"  with a static fundamental charge (ladder j=1/2,3/2,...): gap = {gq:.4f}  -> still GAPPED")
print("  (matter in the fundamental, at the static level, does not close the gauge-invariant gap)")
assert gq>1e-3

# ----------------------------------------------------------------------------------
Hd("VERDICT — does the non-abelian gap-closing survive Gauss projection?")
print(f"""  NO -- it LIFTS, like the U(1) case. Two facts:
   * the non-abelian crossings are (2j+1)^2 gauge-orbit copies (A), removed by Gauss's law;
   * the Gauss-projected physical SU(2) gauge spectrum is GAPPED across the whole coupling
     scan (min gap {mingap:.4f}), and stays gapped with a static fundamental charge (C).

  CONSEQUENCE -- I am downgrading my OWN earlier non-abelian warning:
  gap_under_twist closed the CSS gap by conjugating ONLY the X-half by g=exp(i*theta*T1) --
  an INCONSISTENT (non-gauge-invariant) deformation (matter not co-transformed). The PROPER
  non-abelian gauge treatment (dynamical gauge + Gauss's law) is gapped. So that background
  closing was very likely a PROXY ARTIFACT of the inconsistent twist, NOT a physical gauge
  obstruction -- the same lesson as the U(1) degeneracy. Both gauge-sector fragility signals
  downgrade.

  WHAT REMAINS OPEN (now sharper, and it is the real SMG question):
  This is PURE GAUGE. It shows the non-abelian GAUGE SECTOR is fine (gapped, crossings are
  gauge copies). It does NOT test the MATTER/MIRROR SMG gap under non-abelian gauging -- the
  thing that actually has to survive to gap the mirror. That coupling (matter in the right
  reps + the actual CSS X-half as the gauge-cell plaquette operator + Gauss's law, then the
  gap scan) is the genuine item-13 construction and is NOT settled here.
  On the gap trend: pure SU(2) gauge gap GROWS with beta here (no gap-closing in scan), whereas
  the U(1)+matter proxy SHRANK -- so the 'shrinking' is model-dependent (likely the matter
  coupling), not a universal gauge feature. A proper continuum-limit statement still needs a
  physical-units analysis, but the pure non-abelian gauge sector shows no closing.

  NET: the gauge-sector warnings (U(1) dynamical degeneracy AND my non-abelian background
  closing) were gauge/proxy artifacts and downgrade. The honest open problem narrows to:
  does the MATTER/mirror gap survive when the ACTUAL TCH X-plaquette operator is coupled to
  matter under Gauss's law? That is item 13, still unbuilt.""")
print("\nConstruction asserts passed; gauge-sector gap is GAPPED (reported); matter gap OPEN.")
