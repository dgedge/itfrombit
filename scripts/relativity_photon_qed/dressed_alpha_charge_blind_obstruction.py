#!/usr/bin/env python3
r"""DRESSED-ALPHA: the charge-blind vs charge^2-weighted obstruction (why it is open).

This consolidates the dressed-alpha no-go into one structural statement, and
explains *why* the bare integer is derivable but the dressing is not.

  bare:     alpha_0^{-1} = Sym^2(16)+1 = 137   -- a CHARGE-BLIND count
  dressing: alpha_phys^{-1}-137 = 0.035999...  -- a CHARGE^2-WEIGHTED (Ward) kernel

The observed shift, written delta = N_eff * alpha_0/(2 pi), needs N_eff ~ 31.
N1 = 2*16-1 = 31 reproduces it, but 31 is a mode count (charge-blind: every
fermion slot, including Q=0, counts once). The physical low-energy coupling is
Ward/Kubo billed by the photon self-energy, i.e. charge^2-weighted: sum Q_f^2 N_c.
Those two functionals of the same spectrum cannot be equal -- a charge-blind
count over-weights neutral modes and ignores Q^2 -- so no affine sector-billing
map can carry one into the other (the no_go result). Self-asserting.
"""
import math

ALPHA0 = 1.0/137.0
ALPHA_PHYS_INV = 137.035999084
def ok(c,m): print(("  PASS " if c else "  FAIL ")+m); assert c,m

print("="*72); print("DRESSED ALPHA -- charge-blind vs charge^2 obstruction"); print("="*72)

# 1) bare is a charge-blind count
sym2_16_plus1 = 16*17//2 + 1
ok(sym2_16_plus1 == 137, "bare alpha_0^{-1} = Sym^2(16)+1 = 137 is a charge-BLIND count")

# 2) the dressing written as a one-loop-style kernel needs N_eff ~ 31
delta_inv = ALPHA_PHYS_INV - 137.0
N_req = delta_inv * 2.0*math.pi / ALPHA0
print(f"\n  observed shift  delta(1/alpha) = {delta_inv:.6f}")
print(f"  required kernel N_eff = delta*2pi/alpha0 = {N_req:.3f}")
ok(30.8 < N_req < 31.1, "observed dressing needs N_eff ~ 31")

# 3) N1 = 2*16-1 = 31 reproduces it -- but it is a MODE count (charge-blind)
N1 = 2*16 - 1
ok(N1 == 31, "N1 = 2*16-1 = 31 (mode count) reproduces the shift")
ok(abs(N1 - N_req) < 0.2, "the mode-count 31 matches the required kernel numerically")

# 4) the Ward-legitimate charge^2-weighted kernels (sum Q^2 N_c)
lep  = 3 * (1.0**2) * 1          # e,mu,tau
up   = 3 * ((2/3)**2) * 3        # u,c,t
down = 3 * ((1/3)**2) * 3        # d,s,b
dirac = lep + up + down          # = 8
weyl_doubled = 2*dirac           # = 16  (chiral doubling, generous)
generous = 17.62                 # + QCD/two-loop headroom (from consolidation)
print(f"\n  charge^2-weighted Dirac sum Q^2 N_c = {dirac:.2f}")
print(f"  Weyl-doubled (generous)             = {weyl_doubled:.2f}")
print(f"  generous max (QCD/two-loop)         = {generous:.2f}")
ok(abs(dirac-8.0) < 1e-9, "honest charge^2 Dirac kernel = 8")
ok(generous < 0.60*N_req, "every charge^2-weighted kernel UNDERSHOOTS 31 by >1.7x")

# 5) the obstruction: charge-blind count != charge^2-weighted kernel
print("\n  OBSTRUCTION:")
print("   * 31 counts modes charge-BLIND: a Q=0 mode counts the same as the electron,")
print("     and fractional charges are NOT down-weighted by Q^2.")
print("   * the physical coupling is Ward/Kubo charge^2-weighted: sum Q_f^2 N_c <= 17.6.")
print("   * no affine map sends a charge-blind count to a charge^2 kernel (no_go).")
ok(N1 > generous, "the charge-blind count exceeds the largest charge^2 kernel -> not identifiable")

print("\n"+"="*72); print("VERDICT")
print("  Bare 137 derives BECAUSE it is a charge-blind count, which is exactly the")
print("  kind of object the record-monitor produces. The dressing is charge^2-weighted")
print("  (Ward), a different functional of the spectrum; the honest kernel undershoots")
print("  by ~1.7-3.9x, and N1=31 only matches by being the WRONG (charge-blind) weight.")
print("  => dressed alpha is open by a STRUCTURAL obstruction, not a missing integral.")
print("  Recommendation: carry as a bounded EM residual / external QED renormalization;")
print("  do NOT promote N1=31. Closing it needs a Ward-compatible charge^2 mechanism")
print("  that reaches 31 -- which the charged spectrum cannot supply. exit 0")
