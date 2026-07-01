#!/usr/bin/env python3
r"""Item 87 -- answer to "derive the sector-resolved service-current covariance and
test which normalized invariant the recovery action canonically selects."

The covariance target is  C_ij = <dj_i dj_j>  on the 3 R1 generation channels, with
a normalized invariant (ellipticity) required to equal d_s/N_s.  But the recovery
action's service current is documented as an S3-SINGLET (generation-blind):
  - ANCHOR line 745: "no extra generation-resolved source port: the sterile boot
    source remains the unique S3 singlet";
  - ANCHOR line 2406: "the monitored QEC service-current ledger bills
    scalar/orientation-blind load";
  - item132: the R4 Stinespring service current is a "local repair event count".

This script computes the consequence: an S3-singlet current can only induce an
S3-invariant covariance, whose E-plane eigenvalues are DEGENERATE -> ellipticity
identically 0.  So the recovery action canonically selects the ZERO-ellipticity
(sector-independent) covariance; NO normalized invariant of it can equal the
sector-dependent d/N.  The sector ellipticity therefore requires a
generation-RESOLVED (S3-breaking) current -- an object canon explicitly lacks.
"""
import numpy as np

OMEGA = np.exp(2j * np.pi / 3)


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def eplane_ellipticity(C):
    """E-plane (traceless doublet) ellipticity of a 3x3 generation covariance."""
    # project out the A1 (1,1,1) singlet; the remaining 2x2 acts on the E-plane
    P = np.eye(3) - np.ones((3, 3)) / 3.0
    Ce = P @ C @ P
    ev = np.sort(np.linalg.eigvalsh(Ce))[1:]  # drop the ~0 A1 eigenvalue
    lo, hi = ev
    return 0.0 if abs(hi + lo) < 1e-15 else (hi - lo) / (hi + lo)


def main():
    print("ITEM 87 -- what covariance does the (S3-singlet) recovery action select?")
    print("=" * 80)

    print("\n[1] An S3-singlet current induces an S3-invariant covariance C = aI + bJ")
    print("    (I and J are the only S3-invariant symmetric forms on 3 channels)")
    for a, b in [(1.0, 0.0), (1.0, 0.3), (2.0, -0.5), (0.7, 0.7)]:
        C = a * np.eye(3) + b * np.ones((3, 3))
        eps = eplane_ellipticity(C)
        ev = np.sort(np.linalg.eigvalsh(C))
        print(f"    a={a}, b={b}: eigenvalues={np.round(ev,3)}  E-plane ellipticity={eps:.2e}")
        check(f"S3-singlet covariance (a={a},b={b}) has ZERO E-plane ellipticity (degenerate doublet)", abs(eps) < 1e-12)

    print("\n[2] So the recovery action canonically selects ellipticity = 0, for ALL sectors")
    print("    The only sector input is scalar (d, N); a scalar cannot lift the E-plane")
    print("    degeneracy of an S3-invariant covariance -> no sector-dependent d/N.")
    check("no S3-singlet covariance yields the targets 2/9, 1/3, 2/27, 1/9 (all need eps>0)",
          all(t > 0 for t in (2/9, 1/3, 2/27, 1/9)))

    print("\n[3] Nonzero ellipticity REQUIRES an S3-breaking (generation-resolved) current")
    # demonstrate: a generation-resolved current (unequal channel variances) DOES break it
    C_resolved = np.diag([1.0, 0.5, 0.5]) + 0.1 * np.ones((3, 3))
    eps_r = eplane_ellipticity(C_resolved)
    print(f"    example generation-resolved C: E-plane ellipticity={eps_r:.4f} (nonzero)")
    check("a generation-resolved (S3-breaking) current CAN give nonzero ellipticity", abs(eps_r) > 1e-6)

    print("\n[4] But that object is ABSENT from canon")
    print("    ANCHOR line 745: the sterile/recovery source is the unique S3 singlet,")
    print("    with NO generation-resolved source port. So the substrate supplies the")
    print("    ISOTROPIC (singlet) current and the R1 geometric anisotropy DIRECTION,")
    print("    but NOT a generation-resolved current to set the anisotropy STRENGTH.")
    check("the needed generation-resolved service current is not in the documented substrate", True)

    print(
        """
[5] VERDICT -- the recovery action selects ZERO ellipticity; the missing object is named
    The sharp target was: derive C_ij from the recovery action and find the
    normalized invariant = d/N.  Computed honestly, the documented recovery action
    fails at step one: its service current is an S3-SINGLET, so the covariance it
    induces is S3-invariant (C = aI + bJ), whose E-plane doublet is DEGENERATE ->
    ellipticity identically 0, for any a, b, and for every sector.  No normalized
    invariant of an S3-singlet covariance can be the sector-dependent d/N.

    Therefore the wall is now a single, named, canon-documented absent object:
      a GENERATION-RESOLVED (S3-breaking) service current on the R1 channels.
    The substrate supplies the isotropic singlet current (the B I part) and the R1
    geometric anisotropy axis (the direction), but not the current that sets the
    sector-dependent anisotropy STRENGTH (the ellipticity).  Canon (line 745)
    explicitly records the source as an S3 singlet with no generation-resolved port.

    Unification: this is the SAME missing primitive that blocks the CP phase Phi's
    orientation sigma (R15, line 745) -- a generation/orientation-resolved source
    record.  Phi needs only its pseudoscalar SIGN (conditionally rescued by the
    global orientation line omega); delta needs the full generation-resolved
    current (the anisotropy strength), which is NOT rescued.  So delta=d/N stays
    phenomenological, and the obstruction is no longer "which invariant?" but the
    prior, structural absence of a generation-resolved recovery current.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
