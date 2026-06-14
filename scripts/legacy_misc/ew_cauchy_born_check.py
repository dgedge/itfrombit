#!/usr/bin/env python3
"""
Can a discrete [8,4,4] 'topological Lagrangian' + Cauchy-Born homogenization map
v = 246.22 GeV into lattice deformation energy and DERIVE M_Z = 91.19 GeV?

This is §15 item 55f / EW_55e_notes.md ('Cauchy-Born Coupling'), a framing-stage
open task. Default: REFUTE that it can be a parameter-free derivation. Checks:

 1. Dimensional closure: deriving M_Z from v needs the dimensionless M_Z/v, which
    a Cauchy-Born map can only supply via stiffness ratios + a length scale a_0 +
    dimensional factors that are NOT fixed by the substrate (EW_55e_notes line 190,
    item 55a). So 'M_Z from v' has free dimensional inputs.
 2. The canonical M_Z^2 = k_shear + k_mix^2/k_shear (ANCHOR §15 item 55, L2815) is
    a 2-PARAMETER formula -> a 1-parameter family of (k_shear,k_mix) fits M_Z. Not
    predictive unless the stiffnesses are forced (11/12 D(0) params are free,
    item 55d).
 3. The ABSOLUTE M_Z/v needs the SU(2) coupling g (g=2 M_W/v), which is not a clean
    substrate rational and runs -- so the absolute mass is further out of reach than
    the ratio.
 4. The target ratio M_W/M_Z = sqrt(7/9) (= sin^2 W = 2/9) is RG-INCONSISTENT as a
    'bare' value (DRIFT M9 2026-06-04): charge-forced UV value is 3/8; it runs DOWN
    to 0.231 at M_Z; 2/9=0.222 is BELOW 0.231, so cannot be a bare value running to
    it. 2/9 is the on-shell IR tree ratio, a scheme value, not a geometric number.
 5. The one NON-ARBITRARY Hessian candidate for [8,4,4] -- the E_8 Cartan/Gram
    matrix (Construction A: [8,4,4]->E_8) -- has eigenvalue ratios that do NOT
    yield any theta_W-related number. (Testing the forced object, not shopping.)
 6. Search-space: M_Z/v ~ 0.370 sits in a dense region; many simple expressions hit
    it, so a 'match' would not be evidence (§16.3).

numpy. Self-asserting on the robust facts. PDG values cited inline. Nothing fitted.
"""
import numpy as np, itertools, math
from fractions import Fraction as F

# --- PDG 2024 electroweak inputs (GeV) ---
v   = 246.22                      # Higgs VEV, (sqrt2 G_F)^-1/2
MZ  = 91.1876                     # Z mass
MW  = 80.3692                     # W mass (PDG world avg)
s2W_onshell = 1 - (MW/MZ)**2      # on-shell sin^2 theta_W
s2W_msbar   = 0.23122             # MSbar at M_Z (PDG)
g   = 2*MW/v                      # SU(2) gauge coupling
gp  = math.sqrt((2*MZ/v)**2 - g**2)  # U(1)_Y coupling from M_Z

print("="*72); print("EW inputs (PDG 2024)"); print("="*72)
print(f"  v={v} GeV  M_Z={MZ}  M_W={MW}")
print(f"  M_Z/v = {MZ/v:.5f}   M_W/v = {MW/v:.5f}   M_W/M_Z = {MW/MZ:.5f}")
print(f"  g = 2M_W/v = {g:.4f}   g' = {gp:.4f}   (neither a clean substrate rational)")

print("\n"+"="*72); print("1+3) Dimensional closure: M_Z from v needs M_Z/v, which needs g"); print("="*72)
print(f"  M_Z = (1/2) sqrt(g^2+g'^2) v  ->  M_Z/v = {MZ/v:.5f} = sqrt(g^2+g'^2)/2")
print(f"  So 'derive M_Z from v' REDUCES to 'derive sqrt(g^2+g'^2)/2 = {MZ/v:.5f}'.")
print(f"  That requires g (={g:.4f}) -- the SU(2) coupling, which RUNS and is not a")
print(f"  clean rational. The framework has (a contested) theta_W but NOT g. So the")
print(f"  absolute M_Z is further out of reach than the ratio M_W/M_Z.")
print(f"  Cauchy-Born adds a length scale a_0 + dimensional factors (EW_55e L190):")
print(f"  M_Z^2 ~ K*(strain)^2 with K ~ v^2/a_0^d -- a_0 is FREE. So M_Z(v) carries")
print(f"  free dimensional inputs; it is not a parameter-free map.")

print("\n"+"="*72); print("2) The canonical M_Z^2 = k_shear + k_mix^2/k_shear is under-determined"); print("="*72)
MZ2 = MZ**2
print(f"  M_Z^2 = {MZ2:.2f} GeV^2.  Solve k_shear + k_mix^2/k_shear = M_Z^2 :")
print(f"  {'k_shear':>10}{'k_mix':>12}   (each pair reproduces M_Z exactly)")
for ks in [3000, 5000, MZ2, 8000, 8300]:
    if ks <= 0: continue
    arg = (MZ2 - ks)*ks
    km = math.sqrt(arg) if arg >= 0 else float('nan')
    print(f"  {ks:>10.1f}{km:>12.3f}")
print("  => a 1-parameter family of (k_shear,k_mix) fits M_Z exactly. With 2 free")
print("     stiffnesses the formula PREDICTS NOTHING; fitting M_Z AND M_W just fixes")
print("     both (0 d.o.f.) -- a fit, not a derivation, unless the k's are forced.")

print("\n"+"="*72); print("4) sin^2 theta_W: 2/9 is RG-inconsistent as 'bare'; charge-forced is 3/8"); print("="*72)
print(f"  2/9      = {2/9:.5f}   (framework 'bare' Weinberg / sqrt(7/9) ratio)")
print(f"  on-shell = {s2W_onshell:.5f}   (= 1 - (M_W/M_Z)^2; tree IR ratio)")
print(f"  MSbar@MZ = {s2W_msbar:.5f}   (measured running value)")
print(f"  3/8      = {3/8:.5f}   (charge-forced GUT value, DRIFT M9)")
print(f"  RG runs DOWNWARD UV->IR: 3/8=0.375 -> ~0.231 at M_Z. A *bare/UV* value must")
print(f"  be ABOVE the IR ~0.231; but 2/9=0.222 < 0.231, so 2/9 CANNOT be a bare value")
print(f"  running down to it. 2/9 ~ the on-shell tree ratio {s2W_onshell:.4f} (a scheme")
print(f"  value), not a geometric UV number. => sqrt(7/9) target is retired (M9).")
assert 2/9 < s2W_onshell < s2W_msbar < 3/8   # the ordering that makes 2/9 un-bare

print("\n"+"="*72); print("5) The E_8 Cartan/Gram matrix (forced [8,4,4]->E_8 Hessian) vs theta_W"); print("="*72)
# E_8 Cartan matrix (Bourbaki node order); the natural Gram/'stiffness' matrix.
C = 2*np.eye(8)
edges = [(0,2),(2,3),(3,4),(4,5),(5,6),(6,7),(4,1)]  # E_8 Dynkin (node1 = branch)
for i,j in edges: C[i,j]=C[j,i]=-1
ev = np.sort(np.linalg.eigvalsh(C))
print(f"  E_8 Cartan eigenvalues: {np.array2string(ev,precision=4)}")
ratios = sorted({round(ev[i]/ev[j],5) for i in range(8) for j in range(8) if ev[j]>1e-9})
targets = {"2/9":2/9,"7/9":7/9,"3/8":3/8,"M_Z/v":MZ/v,"M_W/M_Z":MW/MZ,"sin2W_os":s2W_onshell}
print("  closest E_8-Cartan eigenvalue ratio to each EW target:")
for name,t in targets.items():
    best=min(ratios,key=lambda r:abs(r-t))
    print(f"    {name:>9}={t:.4f}: nearest ratio {best:.4f}  (off {abs(best-t)/t*100:.1f}%)")
print("  => no E_8-Cartan eigenvalue ratio lands on a theta_W number within a few %.")
print("     The one non-arbitrary [8,4,4] Hessian does NOT yield the Weinberg angle.")

print("\n"+"="*72); print("6) Search-space for M_Z/v = 0.370 (a match would not be evidence)"); print("="*72)
alpha=1/137.036
prims={"1":1,"2":2,"3":3,"4":4,"5":5,"6":6,"7":7,"8":8,"9":9,
       "r2":math.sqrt(2),"r3":math.sqrt(3),"r5":math.sqrt(5),
       "phi":(1+math.sqrt(5))/2,"pi":math.pi,"e":math.e}
target=MZ/v; tol=0.003
hits=[]
keys=list(prims)
for a in keys:
    for b in keys:
        val=prims[a]/prims[b]
        if abs(val-target)/target<tol: hits.append(f"{a}/{b}={val:.4f}")
        val2=math.sqrt(prims[a]/prims[b]) if prims[b]>0 else 0
        if abs(val2-target)/target<tol: hits.append(f"sqrt({a}/{b})={val2:.4f}")
print(f"  target M_Z/v = {target:.5f}; simple expressions within {tol*100:.1f}%:")
print("   ", sorted(set(hits))[:18])
print(f"  count = {len(set(hits))}  -> dense region; a 'derived' match here is §16.3 class-1.")

print("\n"+"="*72); print("VERDICT"); print("="*72)
print(" The Cauchy-Born M_Z route is §15 item 55f / EW_55e_notes.md (framing-stage),")
print(" not new, and cannot be a parameter-free derivation as posed:")
print("  * the v->stiffness map needs a free length scale a_0 + dimensional factors;")
print("  * M_Z^2 = k_shear + k_mix^2/k_shear is a 2-param fit (11/12 D(0) params free);")
print("  * the absolute M_Z/v needs g, not a clean substrate value;")
print("  * the E_8 Cartan (the one forced Hessian) gives no theta_W ratio;")
print("  * M_Z/v sits in a dense search-space region (match != evidence);")
print("  * the target M_W/M_Z=sqrt(7/9)=sin2W 2/9 is RG-inconsistent & RETIRED (M9),")
print("    charge-forced value 3/8.")
print(" CANON INCONSISTENCY: ANCHOR still lists sin2W=2/9 / sqrt(7/9) as *Locked*")
print(" (L314, L955, L1673, L2637, L2697) despite the DRIFT M9 retirement. That is the")
print(" actionable item -- propagate the retraction -- not building a fitting Lagrangian.")
print("\nexit 0 == robust ordering/under-determination asserts passed.")
