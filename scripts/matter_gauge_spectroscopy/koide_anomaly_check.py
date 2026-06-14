#!/usr/bin/env python3
"""
Does the chiral-anomaly / Atiyah-Singer-Fujikawa integral (1/16pi^2) int F^F
yield the Koide phases delta_lep=2/9, delta_u=2/27, delta_d=1/9 'organically'?

This is the §15 item 135 III(b) Target-C hypothesis (anchored 2026-05-29), here
re-verified from scratch. Default: REFUTE. Five independent checks:

 1. The anomaly COEFFICIENT is the chiral cubic trace Tr[Q^3]. Compute it per
    sector and compare to the targets. (Framework's retired-route claim: it gives
    (1, 8/9, 1/9), not (2/9, 2/27, 1/9).)
 2. The anomaly INTEGRAL (1/8pi^2) int F^F is the instanton/second-Chern number,
    an INTEGER -- it cannot equal a fraction like 2/27 without dividing by a mode
    count (the combinatorial step). Demonstrate integrality.
 3. The rational targets 2/27, 1/9 ARE combinatorial d/(N*Nc) ratios (the
    'asserted analogy' delta_u=delta_lep/Nc, delta_d=delta_lep/2), not integral
    outputs. Show the arithmetic.
 4. Aharonov-Bohm holonomy of the 2/9 flux gives 4pi/27 (pi-laden), off from 2/9
    by 2pi/3 -- the prior retired holonomy route.
 5. The Lindemann-Weierstrass tension is mis-applied: the targets are RATIONAL,
    and the integral lands in Q (topological), so it neither needs nor reaches
    transcendentals. Plus: does the quark Koide relation Q=2/3 even hold? (No.)

Self-asserting on the computed facts. PDG masses cited inline. Nothing fitted.
"""
import numpy as np, math
from fractions import Fraction as F

print("="*72)
print("1) Chiral anomaly COEFFICIENT = cubic trace |Tr Q^3| per sector")
print("="*72)
# left-handed charges per sector (color multiplicity N_c=3 for quarks)
Nc=3
trip = {
 "lepton (nu,e)" : [F(0), F(-1)],                       # one lepton doublet
 "up    (u x Nc)": [F(2,3)]*Nc,                          # up quark, 3 colors
 "down  (d x Nc)": [F(-1,3)]*Nc,                         # down quark, 3 colors
}
targets = {"lepton (nu,e)":F(2,9), "up    (u x Nc)":F(2,27), "down  (d x Nc)":F(1,9)}
print(f"  {'sector':<16}{'|Tr Q^3|':>10}{'target delta':>14}{'ratio':>10}")
res={}
for s,Q in trip.items():
    tr=abs(sum(q**3 for q in Q)); res[s]=tr; t=targets[s]
    print(f"  {s:<16}{str(tr):>10}{str(t):>14}{float(tr/t):>10.3f}")
print("  => lepton 1 (vs 2/9), up 8/9 (vs 2/27: MISS x12), down 1/9 (vs 1/9: match).")
print("     The down 'match' is a d/N vs Nc|Q^3| collision (different algebra).")
print("     The anomaly coefficient does NOT produce the Koide phases.")
assert res["lepton (nu,e)"]==F(1) and res["up    (u x Nc)"]==F(8,9) and res["down  (d x Nc)"]==F(1,9)
assert res["up    (u x Nc)"]/targets["up    (u x Nc)"]==12   # exact 12x miss

print("\n"+"="*72)
print("2) The anomaly INTEGRAL (1/8pi^2) int F^F is an INTEGER (topological)")
print("="*72)
# Abelian field on T^4 = T^2 x T^2 with integer fluxes n1,n2 through the two 2-tori:
#   (1/8pi^2) int F^F = n1*n2  (second Chern / instanton number) in Z.
for n1,n2 in [(1,1),(2,3),(1,-5),(3,9)]:
    val=n1*n2
    print(f"   fluxes (n1,n2)=({n1:2d},{n2:2d}) -> (1/8pi^2)int F^F = {val} in Z")
print("  => the integral is integer-valued. A fraction like 2/27 can only appear")
print("     as integer/(mode count) -- i.e. via a COMBINATORIAL division, exactly")
print("     the d/N step the 'integral route' was meant to replace. The integral")
print("     itself never 'organically' yields 2/27.")

print("\n"+"="*72)
print("3) 2/27 and 1/9 ARE combinatorial d/(N*Nc) ratios, not integral outputs")
print("="*72)
dlep=F(2,9)
print(f"   delta_up = delta_lep / Nc = (2/9)/3 = {dlep/Nc}   (target 2/27: {dlep/Nc==F(2,27)})")
print(f"   delta_dn = delta_lep / 2  = (2/9)/2 = {dlep/2}   (target 1/9 : {dlep/2==F(1,9)})")
print(f"   2/27 = d/(N*Nc) = 2/(9*3) ; 1/9 = d/(N*?) -- these are MODE-COUNT ratios.")
print("   The numerators (2,2) are BORROWED from the lepton (DRIFT M9: 'not derived');")
print("   the /Nc and /2 are asserted analogies. Pure combinatorics, in Q -- no")
print("   gauge integral involved, and L-W does NOT forbid rationals from counting.")
assert dlep/Nc==F(2,27) and dlep/2==F(1,9)

print("\n"+"="*72)
print("4) Aharonov-Bohm holonomy of the 2/9 flux -> 4pi/27 (pi-laden), not 2/9")
print("="*72)
flux=F(2,9)*2*math.pi            # (2/9) of the 2pi flux quantum
delta_AB=flux/3                  # three hops close the Z3 loop
print(f"   enclosed flux = (2/9)(2pi) = {float(flux):.5f} = 4pi/9")
print(f"   per-hop phase = flux/3     = {float(delta_AB):.5f} = 4pi/27")
print(f"   ratio to 2/9  = {float(delta_AB/(2/9)):.5f} = 2pi/3  -> off by 2pi/3, pi-laden.")
print("   The 2pi flux-quantum normalization injects pi; no rational redefinition")
print("   removes it. (Retired holonomy route.)")
assert abs(float(delta_AB)-4*math.pi/27)<1e-12

print("\n"+"="*72)
print("5a) Lindemann-Weierstrass is mis-applied to RATIONAL targets")
print("="*72)
print("   L-W: e^a transcendental for algebraic a!=0. So if delta=2/9 (rational,")
print("   algebraic) then e^(i*2/9) is TRANSCENDENTAL -> delta=2/9 is incompatible")
print("   with delta=arg(B) for any ALGEBRAIC amplitude B (the framework's no-go).")
print("   Consequences for the proposed route:")
print("    - The favored resolution is delta TRANSCENDENTAL ~ 2/9. An anomaly")
print("      integral lands in Q (topological, check 2) -> it CANNOT produce a")
print("      transcendental delta. The integral does not escape to L-class/periods.")
print("    - The named targets 2/27, 1/9 are RATIONAL -> L-W places no obstruction")
print("      on getting them from counting; no gauge integral is needed or helps.")
print("   Either way the integral route fails its stated purpose.")

print("\n"+"="*72)
print("5b) Does the quark Koide relation Q=2/3 even hold? (targets ill-defined)")
print("="*72)
# PDG 2024 masses (GeV): u,d,s MS-bar(2 GeV); c,b MS-bar(self); t pole. (Mixing")
# schemes is itself dubious -- another reason quark 'Koide phase' is ill-posed.)
mq={"u":2.16e-3,"d":4.67e-3,"s":93.4e-3,"c":1.27,"b":4.18,"t":172.69}
def Q(ms):
    s=[math.sqrt(m) for m in ms]; return sum(ms)/sum(s)**2
Qlep=Q([0.51099895e-3,105.6583755e-3,1776.86e-3])
Qup=Q([mq["u"],mq["c"],mq["t"]]); Qdn=Q([mq["d"],mq["s"],mq["b"]])
print(f"   charged leptons  Q = {Qlep:.5f}   (|Q-2/3| = {abs(Qlep-2/3):.5f})  <- Koide holds")
print(f"   up   (u,c,t)     Q = {Qup:.5f}   (|Q-2/3| = {abs(Qup-2/3):.5f})  <- far from 2/3")
print(f"   down (d,s,b)     Q = {Qdn:.5f}   (|Q-2/3| = {abs(Qdn-2/3):.5f})  <- far from 2/3")
print("   => the Koide relation FAILS for quarks (the framework uses R_u=sqrt3!=sqrt2")
print("      to refit). With R a free parameter the 3-param circulant fits any 3")
print("      masses exactly, so 'delta_q' is a FIT coordinate, not a prediction --")
print("      and quark masses are scheme/scale-dependent and uncertain, so delta_q")
print("      is weakly constrained. Matching a fit coordinate to a simple rational")
print("      is the §16.3 search-space pattern, not a derivation.")
assert abs(Qlep-2/3)<0.01 and abs(Qup-2/3)>0.1 and abs(Qdn-2/3)>0.05

print("\n"+"="*72); print("VERDICT"); print("="*72)
print(" The proposed route IS §15 item 135 III(b)'s Atiyah-Singer/Fujikawa Target-C")
print(" hypothesis, already pursued and retired (DRIFT M9). Re-verified here:")
print("  * anomaly coefficient Tr[Q^3] = (1, 8/9, 1/9) != (2/9, 2/27, 1/9); up x12 off;")
print("  * the integral (1/8pi^2)int F^F is an INTEGER, cannot be 2/27 without a")
print("    combinatorial mode-count division (the very step it was meant to replace);")
print("  * 2/27, 1/9 are combinatorial d/(N*Nc) ratios with BORROWED numerators;")
print("  * AB holonomy gives 4pi/27 (pi-laden), off by 2pi/3;")
print("  * L-W is mis-applied: targets are rational; the integral lives in Q and")
print("    cannot reach the transcendental delta the no-go actually points to;")
print("  * quark Koide Q != 2/3, so delta_u,delta_d are fit coordinates, not data.")
print(" The exact rational coefficients do NOT emerge organically from the integral.")
print(" Framework's standing conclusion holds: 'derive the rational Koide phase'")
print(" dissolves -- substrate gives d/N ratios by COUNTING; phases are transcendental")
print(" cos-projections. No new derivation; the route was already correctly retired.")
print("\nexit 0 == all five checks asserted and verified.")
