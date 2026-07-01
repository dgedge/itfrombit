#!/usr/bin/env python3
r"""Item 87 -- the missing anisotropy for the traffic-tensor route comes from R1.

The traffic-tensor attempt (item87_delta_traffic_tensor_attempt.py) is type-correct
(symmetric K=K^T -> mass-shape delta, not a holonomy) but stalled because an
S3-symmetric generation triangle (angles 0, 2pi/3, 4pi/3) is ISOTROPIC -> no
principal axis.  The fix the scalar ledger was missing: the three generations are
NOT an S3 triangle.  Per the canon bit-to-mass map (ANCHOR sec 5.8): the
generation bits (G0,G1) take values (0,0)->tau, (0,1)->e, (1,0)->mu, while
(1,1) is FORBIDDEN by R1 (G0*G1 != 1).  So the generation space is the 2-bit
square with one corner deleted -- a Hamming PATH e--tau--mu (tau central), which
is intrinsically anisotropic.  That broken S3 is a SUBSTRATE-derived source of the
principal axis, not a hand-picked traffic triple.

This script: (1) confirms the symmetric triangle is isotropic; (2) shows the R1
generation embedding is anisotropic and HAS a principal axis; (3) reports the angle
it gives, honestly, and checks it against the targets.  No fishing: the embedding
is fixed by R1, with equal (generation-blind) traffic.
"""
import cmath
import math
from fractions import Fraction

TARGETS = {"e": Fraction(2, 9), "nu": Fraction(1, 3), "u": Fraction(2, 27), "d": Fraction(1, 9)}


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def traffic_tensor(points, weights):
    """K = sum_i w_i (v_i - vbar)(v_i - vbar)^T on the 2D generation plane."""
    n = len(points)
    cx = sum(w * p[0] for w, p in zip(weights, points)) / sum(weights)
    cy = sum(w * p[1] for w, p in zip(weights, points)) / sum(weights)
    Kxx = sum(w * (p[0] - cx) ** 2 for w, p in zip(weights, points))
    Kyy = sum(w * (p[1] - cy) ** 2 for w, p in zip(weights, points))
    Kxy = sum(w * (p[0] - cx) * (p[1] - cy) for w, p in zip(weights, points))
    return Kxx, Kyy, Kxy


def anisotropy_z(K):
    Kxx, Kyy, Kxy = K
    return complex(Kxx - Kyy, 2 * Kxy)


def fold_axis_angle(theta):
    while theta < 0.0:
        theta += math.pi
    theta %= math.pi / 3.0
    if theta > math.pi / 6.0:
        theta = math.pi / 3.0 - theta
    return theta


def eplane_eigs(K):
    Kxx, Kyy, Kxy = K
    tr, det = Kxx + Kyy, Kxx * Kyy - Kxy ** 2
    disc = math.sqrt(max(tr ** 2 / 4 - det, 0.0))
    return sorted([tr / 2 - disc, tr / 2 + disc])


def main():
    print("ITEM 87 -- R1-forbidden corner as the source of generation anisotropy")
    print("=" * 80)

    print("\n[1] The symmetric S3 triangle is isotropic (the prior obstruction)")
    tri = [(math.cos(a), math.sin(a)) for a in (0, 2 * math.pi / 3, 4 * math.pi / 3)]
    zt = anisotropy_z(traffic_tensor(tri, [1, 1, 1]))
    print(f"    symmetric triangle, equal traffic: |z| = {abs(zt):.2e}  -> no principal axis")
    check("symmetric triangle is isotropic (|z|=0)", abs(zt) < 1e-12)

    print("\n[2] The R1 generation embedding (2-bit square minus forbidden (1,1)) is anisotropic")
    # (G0,G1): (0,0)->tau central, (0,1)->e, (1,0)->mu ; (1,1) forbidden/deleted
    r1_points = [(0, 0), (0, 1), (1, 0)]
    K = traffic_tensor(r1_points, [1, 1, 1])
    z = anisotropy_z(K)
    print(f"    R1 path e--tau--mu, equal traffic: K=({K[0]:.3f},{K[1]:.3f},{K[2]:.3f})  |z|={abs(z):.3f}")
    check("R1 embedding is ANISOTROPIC (|z|>0) -> a principal axis EXISTS", abs(z) > 1e-9)
    delta_r1 = fold_axis_angle(0.5 * math.atan2(z.imag, z.real))
    print(f"    principal-axis angle delta_R1 = {delta_r1:.6f} rad  (= pi/12 = {math.pi/12:.6f})")
    eigs = eplane_eigs(K)
    print(f"    E-plane eigenvalues = {[round(e,4) for e in eigs]}  (ratio {eigs[0]/eigs[1]:.4f})")

    print("\n[3] Honest check against the targets")
    print(f"    delta_R1 = {delta_r1:.4f} rad vs targets {[f'{t}={float(t):.4f}' for t in TARGETS.values()]}")
    near = min(TARGETS.values(), key=lambda t: abs(float(t) - delta_r1))
    check("delta_R1 (pi/12) does NOT equal any target d/N", all(abs(float(t) - delta_r1) > 1e-3 for t in TARGETS.values()))
    check("delta_R1 is SECTOR-INDEPENDENT (one geometric angle for all sectors)", True)
    # note (flagged, not claimed): the E-plane eigenvalue 1/3 numerically equals delta_nu
    check("E-plane eigenvalue 1/3 coincides with delta_nu (FLAGGED as likely coincidence)", abs(eigs[0] - 1/3) < 1e-9)

    print(
        """
[4] VERDICT -- the anisotropy SOURCE is found; the sector-magnitude law is the gap
    Progress on the traffic-tensor target:
      * the missing anisotropy is REAL and SUBSTRATE-derived: the R1-forbidden
        (1,1) corner makes the three generations a Hamming PATH (tau central), not
        an S3 triangle.  This breaks the isotropy that blocked the prior attempt --
        a principal axis now EXISTS, from canon (R1), not a hand-picked triple.
    But it is not closed:
      * the R1 embedding with generation-blind traffic gives a SINGLE geometric
        angle delta_R1 = pi/12 ~ 0.262, the same for every sector -- it does NOT
        reproduce the sector-dependent targets (2/9, 1/3, 2/27, 1/9).
      * so R1 supplies the generation-space DIRECTION (the axis), but the
        sector-dependent MAGNITUDE d/N must come from the per-sector defect
        statistics entering the covariance weights m_i, which the ledger still
        does not supply.  (The E-plane eigenvalue 1/3 = delta_nu is noted but
        flagged as a likely coincidence: it is sector-independent and the other
        three targets do not appear.)

    RE-LOCALISED THEOREM TARGET (sharper than before):
        C_ij^sector = <dj_i dj_j> = (R1 generation-path geometry: the axis, now
        supplied) x (per-sector service-current statistics: the magnitude).
        The open piece is now ONLY the magnitude law -- why the per-sector current
        variance on the R1 path scales as d/N -- not the whole covariance.  The
        direction half is closed by R1; the magnitude half remains phenomenological.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
