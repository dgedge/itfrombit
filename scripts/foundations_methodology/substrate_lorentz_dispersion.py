#!/usr/bin/env python3
r"""LORENTZ INVARIANCE ON THE LATTICE — verifying the photon-dispersion claims.

A discrete substrate has a preferred (lattice/CMB) frame, so Lorentz invariance can only be
EMERGENT. The relativity paper (papers/relativity/) argues it is restored for the photon at
accessible energies. This script checks the load-bearing numbers and states the UV scope sharply.

Claims checked (photon = transverse mode of the simple-cubic gauge web, |g(k)|^2 = 6 - 2 sum cos k_i a):
  (1) NO odd (linear) Lorentz violation -- inversion symmetry omega(k)=omega(-k) forbids it (this is
      the term GRB time-of-flight bounds most tightly);
  (2) the dispersion is isotropic at leading order; the FIRST anisotropy is O(k^4) in the action
      (a dimension-six, irrelevant velocity correction Dv/v ~ (a0 k)^2/12), not a marginal O(k^2);
  (3) magnitude for optical photons is ~1e-18, below the SME cavity bound ~1e-17;
  (4) the scope limit: the gauge-web Brillouin ceiling is ~Lambda_QCD ~ 1 GeV, so trans-cutoff
      (e.g. TeV) quanta are not photon Bloch modes and must be represented by the separate framed
      causal-set/null-chain sector (item 150).

exit 0 = inversion holds (no linear LV); the leading anisotropy power is 4 (no O(k^2) anisotropy);
         the optical Dv/v is below the SME bound; the cutoff is ~1 GeV, so TeV quanta are outside
         the lattice-photon support and delegated to the item-150 null-chain sector.
"""
import numpy as np

HBARC = 197.327          # MeV*fm
A0 = 0.595               # fm (gauge-web spacing ~ hbar c / Lambda_QCD)

def lap(kvec):           # simple-cubic gauge-web dispersion omega^2(k) = (2/a^2) sum (1 - cos(k_i a))
    return (2.0 / A0 ** 2) * sum(1 - np.cos(ki * A0) for ki in kvec)

# ================= [1] no linear (odd) Lorentz violation -- inversion symmetry =================
print("[1] NO LINEAR LORENTZ VIOLATION (inversion omega(k)=omega(-k) forbids odd terms):")
rng = np.random.default_rng(0)
maxasym = max(abs(lap(k) - lap(-k)) for k in rng.standard_normal((200, 3)))
print(f"    max|omega^2(k) - omega^2(-k)| over 200 random k = {maxasym:.2e}  -> ZERO odd content")
assert maxasym < 1e-12
print("    -> the dangerous LINEAR-in-energy LV (the GRB time-of-flight term) is structurally absent.")

# ================= [2] leading anisotropy is O(k^4): no marginal O(k^2) anisotropy =================
print("\n[2] LEADING ANISOTROPY IS O(k^4) (axis vs body-diagonal direction):")
axis = np.array([1.0, 0, 0])
diag = np.array([1, 1, 1.0]) / np.sqrt(3)
kappas = np.array([0.02, 0.04, 0.08, 0.16]) / A0           # small k a0
dD = np.array([abs(lap(kp * axis) - lap(kp * diag)) for kp in kappas])
slope = np.polyfit(np.log(kappas), np.log(dD), 1)[0]
print(f"    |omega^2(axis) - omega^2(diag)| ~ k^{slope:.2f}  (=> first anisotropy at O(k^4) in the action)")
assert abs(slope - 4.0) < 0.1                               # power 4: NO O(k^2) anisotropy (it cancels)
print("    -> isotropic at leading order; the residual anisotropy is a dimension-six irrelevant term,")
print("       Dv/v ~ (a0 k)^2/12, not the marginal O(1) effect a naive lattice would give.")

# ================= [3] magnitude for optical photons vs the SME bound =================
print("\n[3] MAGNITUDE for optical photons (Dv/v ~ (a0 k)^2 / 12) vs the SME cavity bound:")
E_opt_eV = 2.0
k_opt = (E_opt_eV * 1e-6) / HBARC                          # 1/fm  (E in MeV / hbarc)
dv_over_v = (A0 * k_opt) ** 2 / 12
SME = 1e-17
print(f"    optical 2 eV: a0 k = {A0*k_opt:.2e};  Dv/v = {dv_over_v:.1e}  vs SME bound ~ {SME:.0e}")
assert dv_over_v < SME
print("    -> photon Lorentz invariance is restored at accessible energies (well below experiment).")

# ================= [4] the Lambda_QCD lattice-photon scope limit =================
print("\n[4] LATTICE-PHOTON SCOPE LIMIT -- the gauge-web UV cutoff:")
E_cut_GeV = np.pi * HBARC / A0 / 1000                      # Brillouin ceiling pi hbarc / a0, in GeV
print(f"    Brillouin ceiling = pi hbar c / a0 = {E_cut_GeV:.2f} GeV ~ Lambda_QCD")
for label, E in [("optical 2 eV", 2e-9), ("GeV photon", 1.0), ("TeV blazar photon", 1e3)]:
    print(f"      {label:<20s} E={E:>7.0e} GeV:  {'WITHIN support' if E < E_cut_GeV else 'TRANS-CUTOFF (outside photon support)'}")
assert E_cut_GeV < 2.0                                      # the cutoff is ~1 GeV, not Planck
print("    -> trans-cutoff (TeV+) quanta -- where LV is tested hardest (GRB, blazars, GZK) -- are NOT")
print("       represented by the gauge-web photon.  They are delegated to the framed causal-set/null-chain")
print("       sector audited under item 150; this script only checks the lattice-photon IR envelope.")

print(f"""
[verdict] LORENTZ INVARIANCE: EMERGENT AND RESTORED FOR THE LATTICE-PHOTON IR ENVELOPE.
  The substrate has a preferred frame, so Lorentz invariance is emergent -- but the photon sector
  passes the tests that kill naive lattices: inversion symmetry forbids the dangerous LINEAR LV [1];
  the leading anisotropy is the irrelevant O(k^4) dimension-six term, not a marginal O(k^2) [2]; and
  the magnitude (Dv/v ~ 1e-18 optical) sits below SME bounds [3]. The Collins-Perez-Sudarsky
  radiative-anisotropy obstruction does not apply to the photon (it was the auxiliary K6 band -- see
  papers/relativity/). The gauge web's cutoff is ~Lambda_QCD (~1 GeV) [4], so trans-cutoff TeV+
  quanta are outside the lattice-photon representation; the current canon handles them as framed
  causal-set/null-chain external events under item 150, with endpoint LSZ residue fixed by the
  service-GNS event isometry.
  TIER: [1]-[3] verify the relativity paper's accessible-energy claims (rigorous); [4] is the
  lattice-envelope scope boundary, shared with the cosmological-constant cutoff (same Lambda_QCD).
exit 0""")
print("ALL ASSERTIONS PASSED -- no linear LV; anisotropy O(k^4); optical Dv/v below SME; cutoff ~1 GeV (trans-cutoff handled by item 150).")
