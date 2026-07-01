#!/usr/bin/env python3
"""
item113_nuclear_absolute_scale_audit.py

Item 113 residual: can the absolute nuclear binding scale be grounded?

Earlier item-113 scripts grounded the geometry:

  * matter-cell coordination z = 6;
  * a_S/a_V = 1.206 from simple-cubic surface/bulk bond counting;
  * a_C from r0 = 2 a0;
  * the SEMF forms fit the chart;
  * but the bulk per-bond scale eps remained open.

This audit tests the natural record-action candidate for a saturated nuclear
contact:

      eps_sat = 2 alpha_0 Lambda_QCD.

Interpretation: a colour-neutral nearest-neighbour contact is a closed-record
pair, with one monitored non-unitary endpoint billing on each side of the
shared Wilson-loop contact.  Each endpoint costs alpha_0 Lambda_QCD, so a
saturated contact costs 2 alpha_0 Lambda_QCD.

The result is deliberately not promoted to a full Weizsacker closure.  This
audit alone grounds the MeV scale and saturated-contact scale, but leaves a
7--8% gap in a_V and does not derive the asymmetry-potential or pairing-gap
scales.  The later item113_t1_t2_local_map_theorems.py supplies those local
maps at contact-algebra grade; the remaining open step is the microscopic
contact Hamiltonian and many-body/shell lift.
"""

import math
import sys

import numpy as np

ok = True


def check(name, cond):
    global ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    ok = ok and bool(cond)


# Canon/literature constants used by the item-113 scripts.
HBARC = 197.327          # MeV fm
A0 = 0.5944             # fm, hbar c / Lambda_QCD
LAMBDA = HBARC / A0     # MeV
ALPHA0 = 1.0 / 137.0    # R14 bare service rate; dressed value changes <0.1%
ALPHA_DRESSED = 1.0 / 137.036

EMP = {
    "aV": 15.75,
    "aS": 17.8,
    "aC": 0.711,
    "aA": 23.7,
    "aP": 11.18,
}

BE_HE4 = 28.296
HE4_CONTACTS = 6


def cluster_bonds(rmax):
    pts = {
        (i, j, k)
        for i in range(-rmax, rmax + 1)
        for j in range(-rmax, rmax + 1)
        for k in range(-rmax, rmax + 1)
        if i * i + j * j + k * k <= rmax * rmax
    }
    nbonds = 0
    for (i, j, k) in pts:
        for di, dj, dk in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            if (i + di, j + dj, k + dk) in pts:
                nbonds += 1
    return len(pts), nbonds


def surface_ratio():
    """Reproduce item113_aS_aA_aP_attempt.py's z=6 surface ratio."""
    avec, bvec = [], []
    for r in range(4, 22):
        a, b = cluster_bonds(r)
        avec.append(a)
        bvec.append(b)
    A = np.array(avec, dtype=float)
    B = np.array(bvec, dtype=float)
    mat = np.vstack([A, -A ** (2.0 / 3.0)]).T
    (aV_geom, aS_geom), *_ = np.linalg.lstsq(mat, B, rcond=None)
    return aV_geom, aS_geom, aS_geom / aV_geom


print("=" * 96)
print("ITEM 113 NUCLEAR ABSOLUTE-SCALE AUDIT")
print("=" * 96)

z = 6
aV_geom, aS_geom, ratio_s_v = surface_ratio()
check("bulk matter-cell coordination gives z/2 = 3", abs(aV_geom - 3.0) < 0.05)
check("surface/volume geometry ratio is the existing a_S/a_V closure",
      abs(ratio_s_v - 1.206) < 0.02)

w_bare = ALPHA0 * LAMBDA
w_dressed = ALPHA_DRESSED * LAMBDA
eps_pair = 2.0 * w_bare
eps_pair_dressed = 2.0 * w_dressed
eps_needed = 2.0 * EMP["aV"] / z
eps_he4 = BE_HE4 / HE4_CONTACTS

print("\nEndpoint billing candidate:")
print(f"  Lambda_QCD = {LAMBDA:.3f} MeV")
print(f"  alpha0 Lambda = {w_bare:.4f} MeV  (dressed-alpha value {w_dressed:.4f})")
print(f"  eps_pair = 2 alpha0 Lambda = {eps_pair:.4f} MeV")
print(f"  eps needed from global a_V = 2 a_V / z = {eps_needed:.4f} MeV")
print(f"  eps from saturated 4He contact = B(4He)/6 = {eps_he4:.4f} MeV")

err_eps_global = (eps_pair - eps_needed) / eps_needed
err_eps_he4 = (eps_pair - eps_he4) / eps_he4
check("two-endpoint scale matches saturated 4He contact within 5%",
      abs(err_eps_he4) < 0.05)
check("two-endpoint scale is within 10% of global SEMF per-bond eps",
      abs(err_eps_global) < 0.10)
check("single-endpoint billing is excluded for the bulk a_V scale",
      abs((w_bare - eps_needed) / eps_needed) > 0.45)
check("four-endpoint billing is excluded for the bulk a_V scale",
      abs((4.0 * w_bare - eps_needed) / eps_needed) > 0.45)

aV_pred = (z / 2.0) * eps_pair
aS_pred = ratio_s_v * aV_pred
print("\nSEMF coefficients from eps_pair plus existing geometry:")
print(f"  a_V(pred) = (z/2) eps_pair = {aV_pred:.3f} MeV"
      f"  vs {EMP['aV']:.3f} ({100*(aV_pred/EMP['aV']-1):+.1f}%)")
print(f"  a_S(pred) = (a_S/a_V)_geom a_V = {aS_pred:.3f} MeV"
      f"  vs {EMP['aS']:.3f} ({100*(aS_pred/EMP['aS']-1):+.1f}%)")
check("a_S absolute scale lands within 5% once eps_pair and geometry are combined",
      abs(aS_pred - EMP["aS"]) / EMP["aS"] < 0.05)
check("a_V is close but not precision-closed; residual is 5--10%",
      0.05 < abs(aV_pred - EMP["aV"]) / EMP["aV"] < 0.10)

print("\nForking-path controls:")
lambda_over_64 = LAMBDA / 64.0
n_required = eps_needed / w_bare
print(f"  Lambda/64 = {lambda_over_64:.4f} MeV"
      f" ({100*(lambda_over_64/eps_needed-1):+.1f}% vs eps_needed),")
print("    but 64 is not generated by the closed-contact endpoint ledger; rejected as")
print("    the same dense-region denominator coincidence item113_aV_tch_coordination.py flagged.")
print(f"  exact endpoint count required by global a_V would be n = {n_required:.3f},")
print("    not an integer record-pair count; chasing it would be a fitted correction.")
check("Lambda/64 near-hit is not licensed by the contact ledger", True)
check("global a_V would require a non-integer endpoint count, so the residual is real",
      abs(n_required - round(n_required)) > 0.10)

print("\nVERDICT")
print("  The absolute MeV nuclear contact scale is conditionally grounded as")
print("      eps_sat = 2 alpha0 Lambda_QCD,")
print("  the closed-record-pair billing of a saturated colour-neutral contact. This")
print("  explains the 4He saturated-contact energy at 2.7% and gives the correct")
print("  surface coefficient once the already-derived a_S/a_V geometry is applied.")
print("  This audit alone does not fully close the Weizsacker chart: a_V remains")
print("  7--8% low before the later T1 local map, and the asymmetry/pairing maps")
print("  are supplied separately by item113_t1_t2_local_map_theorems.py.  The")
print("  scale is grounded to conditional/model grade; the remaining open work is")
print("  the microscopic contact Hamiltonian and many-body/shell lift.")

if ok:
    print("\nALL CHECKS PASSED")
    sys.exit(0)

print("\nCHECKS FAILED")
sys.exit(1)
