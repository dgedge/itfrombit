#!/usr/bin/env python3
"""(i) Is the 0.48% K residual a signal? Propagate the framework's own input uncertainties.
(ii) Horizon rank theorem: are the live large-scale relations sufficient to determine M_P
     intrinsically, or is the horizon structurally an input?
Self-asserting; exit 0 = every number verified."""
import numpy as np

# ---------- inputs (sources quoted) ----------
MP   = 1.220890e19        # GeV, CODATA Planck mass
hbar = 6.582120e-25       # GeV s
H0_planck = 67.36         # km/s/Mpc, Planck 2018
H0_shoes  = 73.04         # km/s/Mpc, SH0ES 2022 (Hubble tension)
OmL  = 0.6847             # Planck 2018
Mpc  = 3.085678e19        # km
alpha0 = 1/137
Lam_anchor   = 0.332      # GeV (canonical 1.4 chiral anchor)
Lam_nucleon  = 0.939565/(2*np.sqrt(2))     # M_N/(2 sqrt2) anchor
Lam_sectors  = (0.3271, 0.3320, 0.3365)    # H7 per-sector best fits (pion, nucleon, glueball)

def K_data(Lam, H0_kms):
    H = H0_kms/Mpc * hbar                  # GeV
    RdS = 1/(H*np.sqrt(OmL))               # GeV^-1
    return Lam**3 * RdS / MP**2            # = C/alpha with C the loop coefficient

K0 = K_data(Lam_anchor, H0_planck)
K_pred = 1.5/alpha0                        # 3/(2 alpha_0) = 205.50
res = K0/K_pred - 1
print(f"(i) K residual vs input floor")
print(f"   K_data(Lam=332 MeV, Planck H0) = {K0:.2f};  K_pred = 3/(2 alpha_0) = {K_pred:.2f};  "
      f"residual = {res*100:.2f}%")
assert abs(K0-206.5) < 0.7 and abs(res) < 0.006
# sensitivity: K ~ Lam^3 / H0
dLam_for_res = abs(res)/3                  # fractional Lam shift absorbing the residual
print(f"   K ~ Lam^3:   absorbing the residual needs dLam/Lam = {dLam_for_res*100:.2f}% "
      f"= {dLam_for_res*332:.1f} MeV")
spread = (max(Lam_sectors)-min(Lam_sectors))/Lam_anchor
print(f"   framework's own Lam softness: anchor variants 332 vs {Lam_nucleon*1000:.1f} (M_N/2sqrt2); "
      f"H7 sector spread = {spread*100:.1f}%  ->  K swings {3*spread*100:.1f}%")
K_shoes = K_data(Lam_anchor, H0_shoes)
print(f"   Hubble tension: K_data(SH0ES H0) = {K_shoes:.2f}  ({(K_shoes/K0-1)*100:+.1f}% vs Planck)")
assert 3*spread > 5*abs(res)               # Lam spread alone dwarfs the residual
assert abs(K_shoes/K0-1) > 10*abs(res)     # Hubble tension dwarfs it further
print(f"   VERDICT: the 0.48% residual is ~{3*spread/abs(res):.0f}x below the Lam-anchor floor and "
      f"~{abs(K_shoes/K0-1)/abs(res):.0f}x below the Hubble-tension floor —")
print(f"   it is SUB-FLOOR: not a signal; chasing a prefactor for it would be 16.3 malpractice.")

# ---------- (ii) horizon rank theorem ----------
print(f"\n(ii) horizon rank theorem: live canon relations as exponent rows over (ln M_P, ln H)")
print(f"     (Lam, alpha treated as substrate-derived knowns; O(1) prefactors and Omega_L dropped)")
rows = {
  "gravity (2/3)alpha route":  (2, 1),     # M_P^2 H = (2/3) alpha Lam^3 / sqrt(OmL)
  "rho_Lambda (item 123) + Friedmann": (2, 1),   # 3 OmL H^2 M_P^2 = c alpha Lam^3 H  ->  M_P^2 H = ...
  "MOND a0 = cH/2pi (item 132)":       (0, 0),   # pins a0, adds no (M_P,H) constraint
}
M_live = np.array([v for v in rows.values()])
r_live = np.linalg.matrix_rank(M_live)
for k,v in rows.items(): print(f"   {k:36s} -> {v}")
print(f"   rank(live) = {r_live} over 2 unknowns (M_P, H)  =>  UNDERDETERMINED: one-parameter family")
assert r_live == 1
# the retired item-136 row would have closed it — and was retired for cause (G7)
M_with136 = np.vstack([M_live, (1, 0.25)]) # M_P = Lam (Lam/H)^{1/4} alpha^-4
r_136 = np.linalg.matrix_rank(M_with136)
print(f"   adding RETIRED item-136 (1, 1/4): rank = {r_136} -> would close — but 136 was retired for")
print(f"   cause (G7: violates the Sakharov codim-1 s=1/2; alpha^4 a category error). The only")
print(f"   closure canon ever had was the relation it correctly retired.")
assert r_136 == 2
# numerical check of the single live constraint
H_pl = H0_planck/Mpc*hbar
lhs = MP**2 * H_pl * np.sqrt(OmL); rhs = (2/3)*alpha0*Lam_anchor**3
print(f"   live constraint check: M_P^2 H sqrt(OmL) = {lhs:.4e} vs (2/3)alpha_0 Lam^3 = {rhs:.4e} "
      f"(ratio {lhs/rhs:.4f} — the same sub-floor residual)")
assert abs(lhs/rhs - 1) < 0.006

print(f"""
THEOREM (within the live canon relation set): all large-scale relations are projections of the
SINGLE constraint  M_P^2 H = O(alpha) Lam^3  — rank 1 over the two unknowns (M_P, H). An
intrinsic M_P is therefore IMPOSSIBLE without one new INDEPENDENT (M_P,H) relation. What would
suffice (each named, none currently exists): (a) a substrate-derived de Sitter entropy S_dS
with a non-circular count (not horizon-area-in-lattice-units, which re-imports R_dS); (b) a
derived rho_Lambda coefficient with a DIFFERENT (M_P,H) scaling; (c) a boot-sequence derivation
of H itself from substrate dynamics. The horizon-as-input status (G7 Dirac relation) is hereby
upgraded from an observation to a RANK statement with the closure requirements named.
ALL ASSERTS PASSED""")
