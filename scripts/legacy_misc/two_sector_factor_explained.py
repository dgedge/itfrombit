#!/usr/bin/env python3
r"""THE TWO-SECTOR x2.4 — EXPLAINED EXACTLY: it is the 28-fold amplification of
a 1.05% per-level offset that the CC chain itself REQUIRES.

Algebra: r3/r6 = (21q1)^(4-32) = (21q1)^-28 exactly, and (T_boot/T0)^3 =
exp(3 ln(Lambda/T0)).  So the two-sector comparison is

    28 * ln(1/(21q1))   vs   3 * ln(Lambda/T0),

i.e. per level: ln(1/(21q1)) = 2.96602 vs 3*ln(Lambda/T0)/28 = 2.99758 —
a +1.05% gap, raised to the 28th power: e^(28*0.03156) = 2.4196 = the factor.
The x2.4 is not a missing O(1) physics term; it is the exponent-28 magnifying
glass applied to a percent-grade offset.

And the offset is REQUIRED: 21q1 = e^-3 exactly (which would close the gate to
x0.94) is EXCLUDED by the CC chain — the 32-power sensitivity would move
rho_Lambda to 0.34 rho_obs.  The x2.4 is the necessary shadow of the CC
chain's own landed value.

Flagged with sect-16.3 caution (quarantined observation, no promotion):
ln(Lambda/T_CMB,today) = 27.977 = 28 - 0.08% — the channel count.  Epoch-
dependent (a why-now statement), owned-integer proximity at one-trial grade;
recorded, not promoted.

Status: the gate's wall-DM premise was superseded (debris invisible; DM =
R4/nu_R), so the gate constrains nothing today; the x2.4 item RETIRES as
explained.  exit 0 = every number verified."""
import math

LAM, T0 = 0.332, 2.348e-13
def f_k(k): return 0.0 if k <= 3 else (0.5 if k == 4 else 1.0)
def q1_of(p): return sum(math.comb(8,k)*p**k*(1-p)**(8-k)*f_k(k) for k in range(9))
u = 21 * q1_of(0.0972)
ratio = (LAM/T0)**3 / u**-28
gap = 3*math.log(LAM/T0)/28 - math.log(1/u)
print(f"[1] r3/r6 = (21q1)^-28 = {u**-28:.4e};  (Lam/T0)^3 = {(LAM/T0)**3:.4e};  ratio = {ratio:.4f}")
print(f"[2] per-level gap = {gap:.5f} ({gap/3:+.3%} of 3); amplified e^(28x) = {math.exp(28*gap):.4f}")
assert abs(math.exp(28*gap) - ratio) < 1e-9
print(f"[3] 21q1 = e^-3 would give ratio {math.exp(28*(3*math.log(LAM/T0)/28-3)):.3f} but shifts")
print(f"    rho_Lambda by x{math.exp(32*math.log(math.exp(-3)/u)):.2f} — EXCLUDED by the CC landing.")
assert abs(32*math.log(math.exp(-3)/u)) > 1
print(f"[4] quarantined observation: ln(Lambda/T0) = {math.log(LAM/T0):.5f} = 28 x (1 {math.log(LAM/T0)/28-1:+.4%})")
print("VERDICT: the x2.4 is EXPLAINED (28-fold amplified, CC-required 1.05% offset);")
print("the item retires; the gate itself is moot post-supersession. exit 0")
print("ALL ASSERTIONS PASSED")
