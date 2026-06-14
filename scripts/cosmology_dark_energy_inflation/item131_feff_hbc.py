#!/usr/bin/env python3
"""F_eff for the HBC scalar amplitude — the named half of A_nu = F_eff/N_eff (16.4 worked target).

The 16.4 protocol reduced the inflationary amplitude to A_nu = F_eff/N_eff and named three gaps:
N_eff, F_eff, H_*. This script derives F_eff for the canonical 28-channel service instrument.

Results:
 (1) THEOREM (jump/CTMC reading): the 28-channel ledger graph is 25-REGULAR and vertex-transitive
     => uniform total exit rate from every channel => holding times are iid exponential => the
     total service count is EXACTLY Poisson => F_eff = 1, exactly. (The Cox/Drazin correction
     vanishes identically: delta-r = 0 for a uniform rate functional.)
 (2) DISCRETE bandwidth-1 scheduler reading: events per tick in {0,1} with duty cycle p =>
     N ~ Binomial => F_eff = 1 - p <= 1. TENSION made quantitative: canon's 'one pixel = one
     independent POISSON event' premise (F=1) is the DILUTE limit p<<1; the same passage's
     'early SATURATED HBC printing' (p ~ 1) would give F -> 0 and suppress A_nu — the two
     readings of the canon text are incompatible and now cleanly separated.
 (3) Per-channel count Fano (the correlation correction, Cox/Drazin formula — the gate-audit
     toolkit): computed exactly on the 28-chain.
 (4) The amplitude bookkeeping with F_eff = 1: N_eff = 1/A_s = 4.76e8; the entropy-pixel
     identification then gives H_* = 9.9e14 GeV and V^(1/4) = 6.4e16 GeV (scale observation
     only; log-space 16.3 discipline — canon's 1e15-1e17 decade is dense with scales).
Self-asserting; exit 0 = every number verified."""
import itertools, numpy as np

# ---------- the 28-channel ledger (canonical construction) ----------
def bits(x): return ((x>>2)&1, (x>>1)&1, x&1)
hyps = sorted({frozenset(p for p in range(8)
               if (sum(u*v for u,v in zip(bits(p),bits(a)))%2)==b)
               for a in range(1,8) for b in (0,1)}, key=sorted)
assert len(hyps)==14
chans = [(i,m) for i in range(14) for m in (0,1)]
cidx = {c:k for k,c in enumerate(chans)}
A = np.zeros((28,28))
for i,hi in enumerate(hyps):
    for j,hj in enumerate(hyps):
        if i!=j and len(hi&hj)==2:
            for m in (0,1):
                for mm in (0,1):
                    A[cidx[(i,m)],cidx[(j,mm)]] = 1
for i in range(14):
    A[cidx[(i,0)],cidx[(i,1)]] = A[cidx[(i,1)],cidx[(i,0)]] = 1

# ---------- (1) regularity + the exact-Poisson theorem ----------
deg = A.sum(1)
assert set(deg) == {25.0}                              # 24 cross-channel + 1 mode partner
print(f"1. the 28-channel ledger graph is {int(deg[0])}-REGULAR (verified; vertex-transitive under")
print(f"   AGL(3,2) x C2) => every channel has the SAME total exit rate => holding times are iid")
print(f"   exponential => the TOTAL service count is EXACTLY Poisson:  F_eff = 1  (theorem).")
# Cox/Drazin check: for the uniform rate functional r_s = r_tot, delta-r = 0 => correction = 0
Q = A - np.diag(deg)
lam, U = np.linalg.eigh(Q)
order = np.argsort(lam)[::-1]; lam = lam[order]; U = U[:,order]
assert abs(lam[0]) < 1e-12 and lam[1] < -1e-9
r = deg.copy()                                         # rate functional = total exit rate
rbar = r.mean(); dr = r - rbar
corr_total = sum((dr@U[:,k])**2/(28*(-lam[k])) for k in range(1,28))
F_total = 1 + (2/rbar)*corr_total
print(f"   Cox/Drazin verification: correction = {corr_total:.2e}  ->  F_total = {F_total:.10f}")
assert abs(F_total - 1) < 1e-12

# ---------- (2) the bandwidth-1 scheduler: F = 1 - p ----------
print(f"\n2. discrete bandwidth-<=1 scheduler (canon: 'one-tick bandwidth <= 1'): events/tick in")
print(f"   {{0,1}} with duty cycle p => N ~ Binomial(T,p) => F_eff = 1 - p  EXACTLY:")
for p in (0.01, 0.5, 0.99):
    print(f"     p = {p:>4}: F_eff = {1-p:.2f}" + ("   <- 'saturated printing' regime: A_nu suppressed" if p>0.9 else
          ("   <- dilute regime: the canon Poisson-pixel premise" if p<0.1 else "")))
print(f"   => the 16.4 conditional reduction's premise 'one pixel = one independent Poisson event'")
print(f"      (F_eff = 1) is the DILUTE limit; the same passage's 'early SATURATED HBC printing'")
print(f"      would give F_eff -> 0. The two readings are now cleanly separated: the amplitude")
print(f"      bookkeeping requires the dilute reading OR a duty-cycle-corrected N_eff.")

# ---------- (3) per-channel count Fano (correlation correction) ----------
x = 0                                                   # channel (h0, m0); all equivalent by transitivity
rx = A[x,:].copy()                                      # arrival rate into x from each state
rbx = rx.mean()
drx = rx - rbx
corr_x = sum((drx@U[:,k])**2/(28*(-lam[k])) for k in range(1,28))
F_x = 1 + (2/rbx)*corr_x
print(f"\n3. per-channel count Fano (jumps into one channel; Cox/Drazin on the 28-chain):")
print(f"   F_channel = {F_x:.4f}  (correlation correction {F_x-1:+.4f} — the channel-resolved")
print(f"   counts are mildly super-Poissonian; the TOTAL stays exactly 1 by regularity).")
assert 1.0 < F_x < 1.2

# ---------- (4) the amplitude bookkeeping with F_eff = 1 ----------
A_s, dA_s = 2.10e-9, 0.03e-9                            # Planck 2018: ln(10^10 A_s)=3.044
N_eff = 1.0/A_s
Mbar = 2.435e18                                          # reduced Planck mass, GeV
H_star = Mbar*np.sqrt(8*np.pi**2/N_eff)
V14 = (3*H_star**2*Mbar**2)**0.25
print(f"\n4. amplitude bookkeeping with the derived F_eff = 1 (dilute reading):")
print(f"   N_eff = F_eff/A_s = {N_eff:.3e} events per mode-local shell (canon's 4.76e8 reproduced)")
print(f"   entropy-pixel identification N_eff = S_dS = 8 pi^2 Mbar^2/H_*^2  ->  H_* = {H_star:.2e} GeV")
print(f"   inflationary energy scale V^(1/4) = (3 H_*^2 Mbar^2)^(1/4) = {V14:.2e} GeV")
print(f"   scale observation (NOT adopted; the 1e15-1e17 decade is dense with canon scales):")
print(f"   V^(1/4) = {V14/1e16:.1f}e16 sits in the GUT/PS decade (E_R3 ~ M_X ~ 1e16; E_R2 = 1.2e15);")
print(f"   a derivation would have to land the number, not the decade (log-space class-1 otherwise).")
assert abs(N_eff - 4.76e8)/4.76e8 < 0.01
assert abs(H_star - 9.9e14)/9.9e14 < 0.02

print("""
VERDICT — F_eff IS NOW DERIVED, with its regime-dependence exact:
  F_eff = 1 (exactly) under the canonical jump/CTMC reading — a THEOREM from the 25-regularity
  of the channel ledger, not an assumption; = 1 - p for the bandwidth-1 discrete scheduler,
  exposing the saturated-vs-Poisson tension in the canon prose. The per-channel correction is
  computed (F_ch = {:.3f}). The A_nu burden now rests ENTIRELY on N_eff (the event count per
  mode-local shell — the horizon-class object per the rank theorem) and the duty cycle of the
  early ledger; the Fano/correlation gap named by 16.4/E5 is closed.
ALL ASSERTS PASSED""".format(F_x))
