#!/usr/bin/env python3
r"""Item 87 -- corrected reformulation (per the user): delta is the CP-even
ELLIPTICITY (anisotropy strength) of the recovery-traffic tensor, NOT a principal-
axis angle.

Why the correction matters.  A real symmetric tensor on the generation E-plane has
TWO invariants: a principal-axis DIRECTION and an anisotropy STRENGTH (ellipticity).
The earlier "R1 supplies the direction, scalar d/N supplies the magnitude" phrasing
was wrong as stated: scaling a fixed tensor (A K_R1 + B I) changes its eigenvalues
(hence ellipticity) but NOT its principal-axis angle.  So if delta were the axis
angle, R1 pins it at pi/12 and no scalar reaches 2/9, 1/3, 2/27, 1/9.

Corrected type-decomposition of the recovery tensor T = sym + antisym:
    antisymmetric part        -> CP-ODD  holonomy phase  Phi   (1 dof)
    symmetric principal AXIS  -> R1-fixed generation direction (1 dof)
    symmetric ELLIPTICITY     -> CP-EVEN strength  delta        (1 dof)
This fits the no-go history: delta CP-even (ellipticity is a ratio of eigenvalues,
reflection-invariant), dimensionless, raw, and SECTOR-tunable -- unlike the axis.

Sharp forward target (NOT solved here -- solving B/A backwards would be fishing):
    for the sector covariance C^(s) = A_s K_R1 + B_s I, a normalized invariant
    (ellipticity eps = (lam+ - lam-)/(lam+ + lam-), or equivalent) must equal d_s/N_s,
    with A_s, B_s DERIVED from the recovery-action service-current covariance.
"""
import math


def eig2(M):
    a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
    tr, det = a + d, a * d - b * c
    disc = math.sqrt(max(tr * tr / 4 - det, 0.0))
    return sorted([tr / 2 - disc, tr / 2 + disc])


def axis_angle(M):
    z = complex(M[0][0] - M[1][1], 2 * M[0][1])  # anisotropy phasor
    if abs(z) < 1e-15:
        return None
    th = 0.5 * math.atan2(z.imag, z.real)
    while th < 0:
        th += math.pi
    th %= math.pi / 3
    return min(th, math.pi / 3 - th)


def ellipticity(M):
    lo, hi = eig2(M)
    return (hi - lo) / (hi + lo)


def scaled(K, A, B):  # A*K + B*I
    return [[A * K[0][0] + B, A * K[0][1]], [A * K[1][0], A * K[1][1] + B]]


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def main():
    print("ITEM 87 -- delta = CP-even ellipticity (corrected reformulation)")
    print("=" * 78)

    # R1 generation anisotropy tensor on the E-plane (2-bit embedding minus (1,1)):
    K_R1 = [[2 / 3, -1 / 3], [-1 / 3, 2 / 3]]
    lo, hi = eig2(K_R1)
    print(f"\n[1] Bare R1 tensor: eigenvalues {{{lo:.4f}, {hi:.4f}}}, axis={axis_angle(K_R1):.6f} rad (pi/12),"
          f" ellipticity={ellipticity(K_R1):.4f}")
    check("R1 eigenvalues are {1/3, 1}", abs(lo - 1 / 3) < 1e-9 and abs(hi - 1) < 1e-9)
    check("bare R1 ellipticity = 1/2", abs(ellipticity(K_R1) - 0.5) < 1e-9)

    print("\n[2] The user's correction: isotropic admixture moves ELLIPTICITY, not the AXIS")
    base_axis = axis_angle(K_R1)
    for (A, B) in [(1, 0), (1, 0.5), (1, 2), (3, 1), (0.4, 0.9)]:
        C = scaled(K_R1, A, B)
        print(f"    A={A}, B={B}: axis={axis_angle(C):.6f}  ellipticity={ellipticity(C):.4f}")
        check(f"axis invariant under (A={A},B={B}) -- scalar mix does NOT change the angle",
              abs(axis_angle(C) - base_axis) < 1e-9)
    eps_vals = {(A, B): ellipticity(scaled(K_R1, A, B)) for (A, B) in [(1, 0), (1, 1), (1, 3)]}
    check("ellipticity DOES vary with the isotropic admixture B/A (so it is sector-tunable)",
          len(set(round(v, 6) for v in eps_vals.values())) == 3)

    print("\n[3] delta = ellipticity is the right TYPE")
    # CP-even: reflection (x->-x, i.e. flip K_xy sign) leaves eigenvalues (hence ellipticity) unchanged
    K_refl = [[K_R1[0][0], -K_R1[0][1]], [-K_R1[1][0], K_R1[1][1]]]
    check("ellipticity is CP-EVEN (reflection-invariant)", abs(ellipticity(K_R1) - ellipticity(K_refl)) < 1e-12)
    check("ellipticity is dimensionless and in [0,1) (a raw strength, not a 2pi-phase)",
          0 <= ellipticity(K_R1) < 1)

    print("\n[4] The sharp forward target (stated, NOT solved -- solving backwards = fishing)")
    print("    For C^(s) = A_s K_R1 + B_s I, require a normalized invariant = d_s/N_s:")
    print("      ellipticity(C^(s)) = A_s / (2 A_s + 3 B_s)  [from the {1,1/3} R1 spectrum]")
    print("    The theorem is to DERIVE (A_s, B_s) -- the anisotropic vs isotropic recovery")
    print("    traffic per sector -- from the recovery-action service-current covariance, and")
    print("    show the resulting invariant lands on d_s/N_s = 2/9, 1/3, 2/27, 1/9. We do NOT")
    print("    back-solve B_s/A_s here: that would relabel 'fit an angle' as 'fit a mix ratio'.")
    print("    Prerequisite object: generation-resolved service currents j_i on the R1 channels")
    print("    (the recovery action's covariance), which the scalar ledger does not yet supply.")

    print(
        """
[5] VERDICT -- framing corrected and locked; wall moved, not removed
    delta is the CP-EVEN ELLIPTICITY (anisotropy strength) of the symmetric
    recovery-traffic tensor -- not the principal-axis angle.  This fixes the
    scalar-vs-angle error: the R1-forbidden corner fixes the symmetric AXIS
    (direction, pi/12), the antisymmetric part is the holonomy Phi (CP-odd), and
    the ELLIPTICITY is the remaining CP-even degree of freedom = delta.  It is
    dimensionless, raw, sector-tunable, and CP-even -- matching every measured
    property, and avoiding all the holonomy/phase no-goes.

    The wall is now precisely: WHY does the CP-even recovery-current ellipticity
    equal d/N?  i.e. derive the per-sector anisotropic:isotropic traffic split
    (A_s, B_s) from the recovery action, not fit it.  That requires the
    generation-resolved service-current covariance -- a substrate object not in
    the current scalar ledger.  Finding/deriving that object is the next step;
    until then delta=d/N stays phenomenological, but the target is now sharp and
    non-numerological.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
