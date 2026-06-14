#!/usr/bin/env python3
r"""T-R2 photon: Ward/Gauss projection gate for the K6 anisotropy.

The open photon problem is not just "does the K6 T1u/Eg anisotropy flow away?"
There is a prior object-identity question:

    Is the anisotropic K6 T1u/Eg branch actually the physical photon?

A physical U(1) photon is a Gauss-projected transverse gauge cohomology. Its
inverse propagator must obey the Ward identity

    Gamma_ij(k) k_j = 0.

The K6 item-102 object is a first-order degenerate k.p matrix on
T1u (+) Eg. This script checks whether that object can be the photon, and then
compares the two possible readings:

  A. K6 branch is photon: gapping Eg leaves a dimension-5 k/Delta residual
     and fails the SME bound at the natural glueball gap.
  B. Photon is Gauss-projected SC Wilson/Maxwell web: the K6 branch is an
     auxiliary pre-Gauss/gauge-web mode; Ward transversality makes the O(k^2)
     kinetic term unique and pushes cubic anisotropy to O(k^4) in Gamma, i.e.
     Delta v/v = O((q/Lambda)^2), which passes at optical q.

exit 0 = K6 coupling reproduced; Ward failure of K6-as-photon verified;
         first-order demotion fails; Gauss-projected Maxwell reading passes.
"""
from __future__ import annotations

import math
import numpy as np


def unit(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)


def c_block(n):
    """Eg <- T1u coupling block fixed by O_h selection rules at unit |k|."""
    nx, ny, nz = unit(n)
    return np.array(
        [
            [-nx / math.sqrt(6), -ny / math.sqrt(6), 2 * nz / math.sqrt(6)],
            [nx / math.sqrt(2), -ny / math.sqrt(2), 0.0],
        ]
    )


def i4(n):
    x, y, z = unit(n)
    return x * x * y * y + y * y * z * z + z * z * x * x


def transverse_basis(n):
    n = unit(n)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(trial, n)) > 0.9:
        trial = np.array([0.0, 1.0, 0.0])
    e1 = trial - np.dot(trial, n) * n
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    e2 /= np.linalg.norm(e2)
    return e1, e2


def helicity_plus(n):
    e1, e2 = transverse_basis(n)
    return (e1 + 1j * e2) / math.sqrt(2)


def k6_velocities(n):
    """Singular values of the K6 T1u/Eg block: item-102 velocity branches."""
    return np.linalg.svd(c_block(n), compute_uv=False)


def ward_residual(n):
    """Longitudinal response of the K6 mixing. A Ward-valid photon has zero."""
    n = unit(n)
    return np.linalg.norm(c_block(n) @ n)


def first_order_shift_coeff(n):
    """Coefficient of the Eg-Schur frequency shift for a Maxwell helicity.

    In a first-order wave Hamiltonian H_eff = curl - C^dag C / Delta, the
    positive-helicity photon gets delta omega = -k^2 coeff / Delta.
    """
    e = helicity_plus(n)
    ctc = c_block(n).T @ c_block(n)
    return float(np.real(np.vdot(e, ctc @ e)))


def sc_group_velocity(q_over_lam, n):
    """Group velocity from the SC Wilson/Maxwell structure factor.

    K(k)=4 sum_i sin^2(k_i/2), omega=sqrt(K), with |k|=q/Lambda.
    """
    n = unit(n)

    def omega(kmag):
        k = kmag * n
        K = 4.0 * sum(math.sin(ki / 2.0) ** 2 for ki in k)
        return math.sqrt(K)

    h = q_over_lam * 1e-4
    return (omega(q_over_lam + h) - omega(q_over_lam - h)) / (2 * h)


print("[1] K6 item-102 velocity harness:")
for name, n, target in [
    ("[100]", (1, 0, 0), math.sqrt(2 / 3)),
    ("[110]", (1, 1, 0), 1 / math.sqrt(2)),
    ("[111]", (1, 1, 1), 1 / math.sqrt(3)),
]:
    s = k6_velocities(n)
    print(f"    {name}: upper singular velocity {s[0]:.9f}  target {target:.9f}")
    assert abs(s[0] - target) < 1e-12
print("    -> the script is testing the canonical item-102 K6 object.")

print("\n[2] Ward/Gauss test on the K6 T1u/Eg mixing:")
for name, n in [("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1))]:
    res = ward_residual(n)
    print(f"    || C(n) n || for {name} = {res:.9f}")
assert ward_residual((1, 0, 0)) > 0.8
assert ward_residual((1, 1, 0)) > 0.2
assert ward_residual((1, 1, 1)) < 1e-12
print("    -> K6 mixing is not Ward-transverse in general. As a photon inverse")
print("       propagator it would couple to longitudinal A, so it is pre-Gauss/")
print("       auxiliary gauge-web data, not automatically the physical photon.")

print("\n[3] If K6 is nevertheless kept as a first-order photon branch, it fails:")
for name, n in [("[100]", (1, 0, 0)), ("[110]", (1, 1, 0)), ("[111]", (1, 1, 1))]:
    coeff = first_order_shift_coeff(n)
    print(f"    helicity Schur coefficient <C^dag C> for {name}: {coeff:.9f}  (I4={i4(n):.9f})")
    assert abs(coeff - i4(n)) < 1e-12
Lam = 0.33
q_opt = 1.0e-9
SME = 1.0e-17
linear_residual = (2.0 / 3.0) * q_opt / Lam
print(f"    [100]-[111] velocity residual from the first-order Schur term:")
print(f"      Delta v/v ~= (2/3) q/Delta = {linear_residual:.2e}")
print(f"      SME bound                 = {SME:.1e}")
assert linear_residual > 1e7 * SME
print("    -> degeneracy-gap demotion is dimension-5 (q/Delta), so the K6-as-")
print("       photon reading is experimentally dead at the natural glueball gap.")

print("\n[4] Ward-projected Maxwell reading:")
print("    At O(k^2), a symmetric gauge inverse propagator has only")
print("      a k^2 delta_ij + b k_i k_j.")
print("    Ward transversality forces a+b=0, hence the unique Maxwell tensor")
print("      k^2 delta_ij - k_i k_j.")
print("    Cubic anisotropy first enters the gauge action at O(k^4), so")
print("      Delta v/v = O((q/Lambda)^2), not O(q/Lambda).")
q_test = 1.0e-3
v100 = sc_group_velocity(q_test, (1, 0, 0))
v111 = sc_group_velocity(q_test, (1, 1, 1))
dvv_test = abs(v100 - v111) / (0.5 * (v100 + v111))
analytic_test = q_test * q_test / 12.0
assert abs(dvv_test - analytic_test) / analytic_test < 1e-3
q_lam = q_opt / Lam
analytic = q_lam * q_lam / 12.0
print(f"    SC Wilson/Maxwell group-velocity anisotropy:")
print(f"      series check at q/Lambda={q_test:.1e}: numeric {dvv_test:.2e}, q^2/12 {analytic_test:.2e}")
print(f"      optical q/Lambda={q_lam:.2e}: small-k q^2/12 = {analytic:.2e}")
print(f"      SME bound          = {SME:.1e}")
assert analytic < SME

print(
    """
[5] VERDICT:
  The hard problem separates cleanly.

  If item-102's K6 T1u/Eg branch is identified with the photon, the route fails:
  the coupling is not Ward-transverse and, after the Eg glueball gap is inserted,
  leaves a dimension-5 q/Delta velocity residual, about eight orders above the
  SME bound for optical photons.

  If the photon is instead the Gauss-projected compact-U(1) Wilson/Maxwell mode
  on the SC gauge web, the T1u/Eg K6 branch is an auxiliary pre-Gauss/gauge-web
  multiplet. Then Ward transversality makes the O(k^2) kinetic term uniquely
  isotropic and the first cubic correction is O(k^4), giving Delta v/v ~ q^2/12
  Lambda^2, below the optical SME bound.

  So T-R2 is not closed by a hidden cancellation in the K6 branch. It is closed
  only by the object-identity theorem: physical photon = Gauss-projected
  transverse Maxwell cohomology, not the raw K6 T1u/Eg velocity branch. That
  theorem is the remaining canon statement to prove or reject.
exit 0"""
)
print("ALL ASSERTIONS PASSED -- Ward fork sharpened; K6-as-photon fails; Gauss photon passes.")
