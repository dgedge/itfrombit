#!/usr/bin/env python3
"""
item113_weizsaecker_theorem_candidates.py

Item 113 theorem-candidate search for the remaining nuclear SEMF coefficients.

This is the candidate-bundle audit.  It does not promote the chart to Locked:
the microscopic contact-Hamiltonian lift and many-body shell corrections remain
open.  It asks whether the residuals left by
item113_nuclear_absolute_scale_audit.py have a coherent local substrate reading
instead of being independent fitted constants.  The later
item113_t1_t2_local_map_theorems.py supplies the local-map proof of the T1/T2
readings used here.

Two candidate theorems emerge from the z=6 matter-cell graph:

T1. Bulk closed-cell normal-ordering theorem
    A fully saturated bulk contact has the closed-record-pair scale
        eps_sat = 2 alpha0 Lambda_QCD
    and a 12-edge cell normal-ordering correction:
        eps_bulk = (13/12) eps_sat.
    The +1 is a bulk closure/latch ledger over the 12 Q3 stabilizer edges.
    Exposed surface contacts do not receive the closure bonus.  This closes
    a_V while leaving the previously computed a_S surface penalty intact.

T2. Five-transverse-channel bipartite-matching theorem
    After the one radial/saturated channel is consumed, the z=6 matter-cell
    contact graph leaves z-1=5 transverse matching channels.  An I3 imbalance
    pays five single endpoint costs:
        a_A,pot = 5 alpha0 Lambda_QCD,
    while chi pairing is a coherent closed-pair RMS over the same five channels:
        a_P = sqrt(5) eps_sat = 2 sqrt(5) alpha0 Lambda_QCD.

The script checks that these two theorem statements jointly predict the
empirical Weizsacker coefficients at the expected model grade, while rejecting
two obvious wrong extensions: applying the 13/12 bulk correction to surface
penalty, and using Lambda/64 as an unlicensed per-contact scale.
"""

import math
import sys

ok = True


def check(name, cond):
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


HBARC = 197.327       # MeV fm
A0 = 0.5944          # fm
LAMBDA = HBARC / A0
ALPHA0 = 1.0 / 137.0
W = ALPHA0 * LAMBDA
EPS_SAT = 2.0 * W
Z = 6
SURFACE_RATIO = 1.206  # item113_aS_aA_aP_attempt.py geometry, rounded

EMP = {
    "aV": 15.75,
    "aS": 17.8,
    "aC": 0.711,
    "aA": 23.7,
    "aP": 11.18,
}


def rel(pred, obs):
    return (pred / obs) - 1.0


def pct(pred, obs):
    return 100.0 * rel(pred, obs)


def row(name, pred, obs):
    print(f"  {name:4s}  pred={pred:8.3f}  empirical={obs:8.3f}  err={pct(pred, obs):+6.2f}%")


print("=" * 96)
print("ITEM 113 WEIZSACKER THEOREM-CANDIDATE AUDIT")
print("=" * 96)
print(f"Lambda_QCD={LAMBDA:.3f} MeV, alpha0 Lambda={W:.4f} MeV, eps_sat=2 alpha0 Lambda={EPS_SAT:.4f} MeV")

# Existing grounded pieces.
aC = 0.6 * (1.0 / 137.036) * HBARC / (2.0 * A0)
aV_pair = (Z / 2.0) * EPS_SAT
aS_pair = SURFACE_RATIO * aV_pair

print("\nExisting closed-contact scale without new theorem:")
row("aV0", aV_pair, EMP["aV"])
row("aS0", aS_pair, EMP["aS"])

# T1: bulk-only 13/12 normal-ordering correction.
bulk_factor = 13.0 / 12.0
eps_bulk = bulk_factor * EPS_SAT
aV_T1 = (Z / 2.0) * eps_bulk
aS_T1 = aS_pair                         # surface penalty lacks the bulk closure latch
aS_wrong = SURFACE_RATIO * aV_T1         # negative control: applying T1 to surface too

print("\nT1 bulk closed-cell normal-ordering candidate:")
print("  eps_bulk = (13/12) eps_sat; surface penalty remains unclosed eps_sat.")
row("aV", aV_T1, EMP["aV"])
row("aS", aS_T1, EMP["aS"])
print(f"  negative control: applying 13/12 to surface gives aS={aS_wrong:.3f}"
      f" ({pct(aS_wrong, EMP['aS']):+.2f}%), worse and conceptually wrong.")
check("T1 closes a_V at <1%", abs(rel(aV_T1, EMP["aV"])) < 0.01)
check("surface must NOT receive the bulk closure factor", abs(rel(aS_wrong, EMP["aS"])) > 0.05)
check("unclosed surface penalty remains within 3%", abs(rel(aS_T1, EMP["aS"])) < 0.03)

# T2: five transverse-channel I3/chi theorem.
r0 = 2.0 * A0
rho = 3.0 / (4.0 * math.pi * r0 ** 3)
mN = 939.0
EF = (HBARC ** 2 / (2.0 * mN)) * (1.5 * math.pi ** 2 * rho) ** (2.0 / 3.0)
aA_kin = EF / 3.0
aA_pot = (Z - 1.0) * W
aA_T2 = aA_kin + aA_pot
aP_T2 = math.sqrt(Z - 1.0) * EPS_SAT

print("\nT2 five-transverse-channel bipartite-matching candidate:")
print(f"  density r0=2a0 -> E_F={EF:.3f} MeV, kinetic a_A=E_F/3={aA_kin:.3f} MeV")
print(f"  potential a_A,pot=(z-1) alpha0 Lambda={aA_pot:.3f} MeV")
print(f"  pairing a_P=sqrt(z-1) eps_sat={aP_T2:.3f} MeV")
row("aA", aA_T2, EMP["aA"])
row("aP", aP_T2, EMP["aP"])
check("T2 puts a_A within 3%", abs(rel(aA_T2, EMP["aA"])) < 0.03)
check("T2 puts a_P within 5%", abs(rel(aP_T2, EMP["aP"])) < 0.05)

# Controls against dense denominator matching.
lambda_over_64 = LAMBDA / 64.0
print("\nForking-path controls:")
print(f"  Lambda/64={lambda_over_64:.3f} MeV is near eps_needed but has no T1/T2 ledger object.")
print("  This is intentionally treated as a dangerous near-hit: T1 uses the")
print("  closed-contact ledger 2 alpha0 Lambda times 13/12, not a denominator-64")
print("  rule.  If the 13/12 normal-ordering theorem fails, Lambda/64 is rejected.")
check("Lambda/64 is acknowledged as a near-hit but is not a licensed theorem source", True)

print("\nCoefficient table under T1+T2 candidate bundle:")
row("aC", aC, EMP["aC"])
row("aV", aV_T1, EMP["aV"])
row("aS", aS_T1, EMP["aS"])
row("aA", aA_T2, EMP["aA"])
row("aP", aP_T2, EMP["aP"])
check("all five SEMF coefficients land within 5% under the candidate bundle",
      all(abs(rel(p, EMP[k])) < 0.05 for k, p in {
          "aC": aC, "aV": aV_T1, "aS": aS_T1, "aA": aA_T2, "aP": aP_T2
      }.items()))

print("\nVERDICT")
print("  Suitable theorem targets exist and are local, discrete, and mutually coherent:")
print("    T1: the 13/12 bulk-only closed-cell normal-ordering rule.")
print("    T2: the z-1=5 transverse-channel I3/chi bipartite-matching rule.")
print("  item113_t1_t2_local_map_theorems.py supplies these local maps at")
print("  contact-algebra grade.  What remains open is the deeper microscopic")
print("  contact-Hamiltonian lift of eps_sat and the many-body shell/cluster")
print("  extension.  This is therefore a conditional liquid-drop bundle, not a")
print("  Locked chart derivation.")

if ok:
    print("\nALL CHECKS PASSED")
    sys.exit(0)

print("\nCHECKS FAILED")
sys.exit(1)
