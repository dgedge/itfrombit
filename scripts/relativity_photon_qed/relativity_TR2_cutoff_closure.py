r"""T-R2 PHOTON CUTOFF — closure by consistency: route (A) costs the CC, route (B) is canonical.

The cutoff residual (relativity_TR2_cutoff_residual.py) posed two routes to bridge the
photon UV:
  (A) derive a refined gauge-web spacing a_gamma << a0, or
  (B) demote the SC Wilson/Maxwell photon to an IR effective action whose microscopic
      cutoff is NOT the matter-cell spacing.

This script CLOSES the residual at the photon-sector level by showing the two routes are
not free choices and route (B) is already canon:

  1. The SC Wilson/Maxwell dispersion is, by the photon paper's own construction
     (photon_paper2_final.tex sec:sc_web), the |k| a0 << 1 long-wavelength ENVELOPE of the
     microscopic line-graph gauge theory -- it was never the UV photon. So demanding it
     support TeV modes at its Brillouin zone is a category error (the same one the
     object-identity theorem fixed at the IR/identity level, now at the UV/cutoff level).

  2. Route (A) is NOT a free choice. Any a_gamma fine enough to even SUPPORT a TeV photon
     as a mode forces a gauge UV scale 1/a_gamma >> Lambda_QCD, which (a) contradicts
     item-118's load-bearing premise "substrate UV cutoff = Lambda_QCD, no Planck-scale
     physics", and (b) re-inflates the photon vacuum energy by (Lambda_gamma/Lambda_QCD)^4,
     reopening the cosmological-constant problem item 118 resolves. So route (A) and the CC
     resolution are MUTUALLY EXCLUSIVE.

  3. Therefore the only item-118-consistent closure is route (B), and it is already canon:
     photon paper sec:sc_web (SC Maxwell = long-wavelength limit) + the item-118
     stroboscopic-envelope validity bound (ALL continuum/EFT descriptions hold only for
     omega << Lambda_QCD). The photon obeys the framework's universal IR-envelope scope,
     not a photon-specific patch.

CONSEQUENCE: the cutoff residual is CLOSED as a T-R2 photon target (dissolved as a category
error). The genuine remainder -- a Lorentz-invariant representation of trans-Lambda_QCD
(omega >> Lambda_QCD) quanta -- is NOT a photon-sector problem; it is a FOUNDATIONAL
single-scale item, shared with and constrained by the CC resolution's Lambda_QCD UV cutoff.

exit 0 = SC EFT-envelope validity verified; route-(A) UV scale computed and shown
         >> Lambda_QCD; CC reopening factor computed; mutual-exclusivity / re-filing holds.
"""
import numpy as np

HBARC = 0.1973269804      # GeV fm
LAMBDA = 0.332            # GeV, Lambda_QCD (sec 1.4 / item 118)
A0 = HBARC / LAMBDA       # canonical substrate spacing (sec 1.4, item 49)

print(f"[0] canon scales: Lambda_QCD={LAMBDA} GeV, a0 = hbar c/Lambda = {A0:.4f} fm (sec 1.4, item 49).")


def omega_SC(K, n):
    """dimensionless SC Wilson/Maxwell lattice dispersion, K = a0|k|, direction n."""
    n = np.array(n, float); n = n / np.linalg.norm(n)
    return 2.0 * np.sqrt(sum(np.sin(K * ni / 2.0) ** 2 for ni in n))


print("\n[1] SC Wilson/Maxwell dispersion is the long-wavelength ENVELOPE (photon paper sec:sc_web):")
print("    K=a0|k|    omega_SC/K (->1 = continuum)   |deviation from linear|")
for K in (1e-2, 1e-1, 1.0, np.pi):
    r = omega_SC(K, (1, 0, 0)) / K
    print(f"    {K:7.4f}    {r:.6f}                  {abs(1 - r):.3e}")
assert abs(omega_SC(1e-2, (1, 0, 0)) / 1e-2 - 1) < 1e-3     # clean continuum only for K << 1
assert abs(omega_SC(np.pi, (1, 0, 0)) / np.pi - 1) > 0.1    # O(1) failure at the BZ edge
E_BZ = np.pi * HBARC / A0
print(f"    => SC EFT validity domain: omega << Lambda_QCD. Brillouin ceiling E_BZ = pi hbar c/a0")
print(f"       = {E_BZ:.3f} GeV. The microscopic photon is the line-graph 12-band theory at the SAME")
print(f"       ~a0 scale (ceiling O(Lambda_QCD)); no substrate mode -- SC web OR line graph -- exists")
print(f"       above ~1 GeV. The SC lattice is the IR regulator, never the photon's UV support.")

print("\n[2] route (A) cost: a_gamma fine enough to SUPPORT high-E photons forces 1/a_gamma >> Lambda_QCD")
print("    (Nyquist E < pi hbar c/a_gamma is the WEAKEST necessary condition -- the mode must merely exist):")
for E, label in ((1e3, "1 TeV"), (1e5, "100 TeV")):
    a_gamma = np.pi * HBARC / E                 # spacing needed for the mode to exist at all
    Lam_gamma = HBARC / a_gamma                 # characteristic gauge UV energy scale
    cc_factor = (Lam_gamma / LAMBDA) ** 4       # photon vacuum-energy inflation vs item-118 Lambda^4
    print(f"    {label:>7}: a_gamma <= {a_gamma:.2e} fm  =>  1/a_gamma ~ {Lam_gamma:7.1f} GeV "
          f"= {Lam_gamma / LAMBDA:.0f}x Lambda_QCD ;  CC inflation (Lam_g/Lam)^4 = {cc_factor:.1e}")
a_TeV = np.pi * HBARC / 1e3
Lam_TeV = HBARC / a_TeV
assert Lam_TeV > 100 * LAMBDA                    # even bare TeV support needs a UV scale >> Lambda_QCD
assert (Lam_TeV / LAMBDA) ** 4 > 1e6            # ... and reopens many OOM of the cosmological constant
print("    => route (A) introduces a UV scale >> Lambda_QCD (>=~960x for mere TeV support; far more under")
print("       real high-energy Lorentz bounds), CONTRADICTING item-118 'substrate UV = Lambda_QCD', and")
print("       re-inflating the photon vacuum energy by >=~1e12 -- REOPENING the cosmological constant.")
print("       Route (A) and the item-118 CC resolution are MUTUALLY EXCLUSIVE.")

print("\n[3] CLOSURE — residual dissolved at the photon-sector level (route B, forced and canonical):")
print("    * SC Wilson/Maxwell photon = |k|a0<<1 IR envelope (photon paper sec:sc_web). The Nyquist/TeV")
print("      'falsification' is a category error: the SC lattice was never the photon's UV support.")
print("    * item-118 stroboscopic validity bound already scopes ALL continuum/EFT descriptions to")
print("      omega << Lambda_QCD. The photon obeys the framework's universal IR-envelope scope.")
print("    * GENUINE REMAINDER (re-filed, NOT a T-R2 photon target): a Lorentz-invariant representation")
print("      of trans-Lambda_QCD (omega >> Lambda_QCD) quanta. FOUNDATIONAL single-scale item -- the very")
print("      a0 = hbar c/Lambda_QCD commitment that RESOLVES the CC (item 118) is what bars a finer photon")
print("      UV. Observed Lorentz invariance of >> Lambda_QCD photons is unexplained as the framework")
print("      stands; closing it needs a second scale, which costs the CC resolution.")
print("\nVERDICT: T-R2 photon cutoff residual CLOSED at the photon-sector level (category error dissolved);")
print("         remainder re-filed as the foundational trans-Lambda_QCD representation problem.")
print("exit 0")
print("ALL ASSERTIONS PASSED — SC = IR envelope; route (A) reopens the CC; route (B) forced & canonical.")
