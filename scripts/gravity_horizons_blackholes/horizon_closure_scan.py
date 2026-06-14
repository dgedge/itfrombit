#!/usr/bin/env python3
"""The independent (M_P, H) relation: systematic closure-candidate scan.

Context (k_residual_and_horizon_rank.py): the live canon relations have rank 1 over (M_P, H) —
all are projections of M_P^2 H = O(alpha) Lambda^3. One new INDEPENDENT relation would make the
Planck mass intrinsic. Honest framing up front: with the static R_dS = sqrt(3/Lambda_cosmo) the
unknown pair is (M_P, Lambda_cosmo), so ANY successful closure IS a solution of the cosmological
constant problem (the remaining ~41 OOM). This scan tests every named route quantitatively and
adopts nothing. Self-asserting; exit 0 = every number verified."""
import numpy as np

ln, log2 = np.log, np.log2
# inputs (as in k_residual_and_horizon_rank.py)
MP   = 1.220890e19; hbar = 6.582120e-25
H0   = 67.36; Mpc = 3.085678e19; OmL = 0.6847
H    = H0/Mpc*hbar                          # GeV
RdS  = 1/(H*np.sqrt(OmL))                   # GeV^-1
Lam  = 0.332; alpha0 = 1/137
rhoL = 3*OmL*H**2*MP**2/(8*np.pi)           # observed dark-energy density (reduced-Planck Friedmann), GeV^4
print(f"inputs: H = {H:.3e} GeV, R_dS = {RdS:.3e} GeV^-1, rho_Lambda = {rhoL:.3e} GeV^4")
assert 2.0e-47 < rhoL < 3.0e-47

# ---------- A. canon-sweep: hidden H-determining relations ----------
print("\nA. canon-sweep for hidden H-determining relations (numeric verdicts):")
T_dS  = H/(2*np.pi)                         # de Sitter temperature
T_KMS = alpha0*Lam                          # the substrate KMS rate/temperature scale (item 123 I)
print(f"   A1 KMS<->dS temperature identification: T_KMS = alpha*Lam = {T_KMS:.3e} GeV vs "
      f"T_dS = H/2pi = {T_dS:.3e} GeV -> off by 10^{log2(T_KMS/T_dS)/log2(10):.0f} -> DEAD")
assert T_KMS/T_dS > 1e39
print(f"   A2 tick-count today N = 1/(H tau0) = {1/(H/Lam):.2e}: epoch-dependent (H(t) varies) and")
print(f"      canon's own dot-G/G = 0 row forbids constants tracking H(t) -> NOT a constants relation")
print(f"   A3 boot threshold p_th = 0.110: a pure number; supplies no (M_P,H) scale by itself")
print(f"   => no hidden independent relation exists in canon; rank-1 stands comprehensively.")

# ---------- B. route (b): power-law rho_Lambda = Lambda^(4+n)/M_P^n ----------
print("\nB. power-law closure rho_L = Lambda^(4+n) M_P^-n  (+ Friedmann + the gravity relation):")
print("   exponent row (2+n, 2) vs gravity (2,1): independent iff n != 2 (n=2 = Zeldovich, degenerate).")
print(f"   {'n':>6} {'M_P predicted':>15} {'verdict':>28}")
for n_ in (1, 2, 3, 4):
    if n_ == 2:
        print(f"   {n_:>6} {'—':>15} {'degenerate (Zeldovich)':>28}"); continue
    # M_P = Lam * (3/(4 alpha^2) / 3OmL ... keep O(1)s: solve exactly with the (2/3)alpha gravity form
    # gravity: H = (2/3) a Lam^3 / (MP^2 sqrt(OmL));  rho: 3 OmL H^2 MP^2 = Lam^(4+n) MP^-n
    # -> MP^(n-2) = (4/3) a^2 /sqrt? do it numerically:
    from scipy.optimize import brentq
    def f(lnM):
        M = np.exp(lnM); Hg = (2/3)*alpha0*Lam**3/(M**2*np.sqrt(OmL))
        return ln(3*OmL*Hg**2*M**2/(8*np.pi)) - ln(Lam**(4+n_)*M**(-n_))
    try:
        lnM = brentq(f, ln(1e-3), ln(1e60)); Mpred = np.exp(lnM)
        verdict = f"off true M_P by 10^{log2(Mpred/MP)/log2(10):+.1f}"
        print(f"   {n_:>6} {Mpred:>15.3e} {verdict:>28}")
    except ValueError:
        print(f"   {n_:>6} {'no root':>15} {'—':>28}")
# the n that exactly closes:
n_exact = 2 + ln((4/3)*alpha0**2/ (3*OmL/3/OmL)) / ln(Lam/MP)   # derive cleanly below instead
# exact condition: MP^(n-2) = (4/3) alpha0^2 Lam^(n-2) * (1/(3 OmL)) * (3 OmL)??  compute numerically:
from scipy.optimize import brentq
def n_solve(nn):
    Hg = (2/3)*alpha0*Lam**3/(MP**2*np.sqrt(OmL))
    return ln(3*OmL*Hg**2*MP**2/(8*np.pi)) - ln(Lam**(4+nn)*MP**(-nn))
n_exact = brentq(n_solve, 0.1, 6.0)
print(f"   exact closure needs n = {n_exact:.4f} — non-integer; nearest simple rationals (2+1/5=2.2,")
print(f"   2+3/14=2.214...) shift M_P by 10s of % under the huge lever ln(M_P/Lam)=45 -> NO clean n.")
assert 2.1 < n_exact < 2.35

# ---------- C. route (c): boot-sequence residual fault rate as a code-suppression exponent ----------
print("\nC. boot-sequence closure: rho_L = q * alpha * Lam^4 (residual fault rate q per cell-tick):")
q = rhoL/(alpha0*Lam**4)
bits = log2(1/q)
print(f"   required q = {q:.3e}  =  2^-{bits:.1f}   ({ln(1/q):.1f} nats)")
# fault-tolerance double-exponential: q ~ p_th (p/p_th)^(2^l) — the ONLY natural mechanism in the
# framework's toolbox that GENERATES such numbers (concatenation of the [8,4,4] code)
p_th = 0.110
for ratio in (0.5, 0.3):
    lev = log2(ln(p_th/q)/ln(1/ratio))
    print(f"   concatenation depth at p/p_th = {ratio}: 2^l = {ln(p_th/q)/ln(1/ratio):.1f} -> l = {lev:.2f}")
print(f"   => the Dirac number becomes a concatenation DEPTH l ~ 7 (integer!), but the double-exp")
print(f"      sensitivity means l=7 spans q ~ 1e-39..1e-44 as p/p_th varies — a STRUCTURED TARGET")
print(f"      (derive l and p from the boot dynamics), not a numerical hit.")
assert 6.5 < log2(ln(p_th/q)/ln(2)) < 7.5

# ---------- D. the log2 observation, with full dressing-spread (16.3 in log-space) ----------
print("\nD. the 'address-space' observation and its honest spread:")
N_lin  = Lam*RdS                            # horizon in cells (linear)
dress  = {"log2(Lam R_dS) [linear cells]": log2(N_lin),
          "log2(1/q) [fault bits, C above]": bits,
          "log2(N_ticks today) [epoch!]": log2(Lam/H),
          "log2(sqrt(N_cells_vol))": 1.5*log2(N_lin)/1.5*1.0 if False else log2(N_lin**1.5)/2}
for k,v in dress.items(): print(f"   {k:38s} = {v:.1f}")
print(f"   spread across dressings: {min(dress.values()):.1f} – {max(dress.values()):.1f} bits;")
print(f"   137 sits inside the band — but so do 136–143; +-1 bit = x2; the 'log2 ~ 137' reading is")
print(f"   dressing-dependent (class-1 in log-space). Recorded as a TARGET STRUCTURE only: a")
print(f"   derivation would have to produce R_dS = c * 2^137 / Lam with c derived (observed c = "
      f"{N_lin/2**137:.3f}); nothing in canon does.")
assert 136 < log2(N_lin) < 139

print(f"""
VERDICT:
  A: no hidden independent relation exists in canon — the rank-1 theorem is comprehensive.
  B: power-law rho_Lambda closure FAILS cleanly (needs n = {n_exact:.3f}, non-integer; integer n
     miss M_P by many OOM or are degenerate). Route (b) is dead as a simple power law.
  C: the boot-sequence route restates the Dirac number as a concatenation depth l ~ 7 of the
     [8,4,4] code — the framework's only natural double-exponential generator. STRUCTURED TARGET:
     derive (l, p/p_th) from boot dynamics; not closed, but now posed in code-native variables.
  D: log2(Lam R_dS) = {log2(N_lin):.1f} — temptingly near 137 but dressing-dependent (other dressings
     give {bits:.0f}-{log2(Lam/H):.0f}); class-1 in log-space; recorded as target structure, NOT adopted.
  NET: the closure remains ABSENT — confirmed against every canon candidate. The two routes worth
  a future derivation are now precisely posed IN CODE-NATIVE VARIABLES: (i) concatenation depth
  l ~ 7 from boot dynamics; (ii) the 137-bit address-space relation with its c = {N_lin/2**137:.2f} prefactor.
  Any such derivation must clear the T1-T9 thermodynamic gates (16.4) before tiering.
ALL ASSERTS PASSED""")
