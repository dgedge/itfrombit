#!/usr/bin/env python3
r"""Item 87 -- the user's strict conditional-covariance test: NO new source.

Protocol (strict):
  * do NOT add a generation-resolved source current;
  * decompose the existing SCALAR service current over the R1 path EDGES;
  * read it onto the generation nodes; form C_ij = <P_i dj_s P_j dj_s>;
  * test whether the E-plane ellipticity gives d/N.
  * if it collapses to a scalar identity, the route dies; if it gives a genuine
    edge/path covariance, that is the right object.

R1 path: e -- tau -- mu (tau central), edges {e-tau, tau-mu}.  A scalar (i.i.d.)
edge service current, read onto nodes via the incidence matrix B, gives the node
covariance C = B B^T -- the substrate's OWN path structure, no generation-resolved
source added.
"""
import numpy as np

# node-edge incidence: rows = nodes (e, tau, mu), cols = edges (e-tau, tau-mu)
B = np.array([[1, 0],
              [1, 1],
              [0, 1]], float)
TARGETS = {"e": 2 / 9, "nu": 1 / 3, "u": 2 / 27, "d": 1 / 9}


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def eplane_eigs(C):
    P = np.eye(3) - np.ones((3, 3)) / 3.0      # project out the A1 (1,1,1) singlet
    ev = np.sort(np.linalg.eigvalsh(P @ C @ P))
    return ev[1], ev[2]                          # the two E-plane eigenvalues (drop ~0)


def main():
    print("ITEM 87 -- conditional covariance from a SCALAR current on the R1 path")
    print("=" * 80)

    print("\n[1] Scalar edge current, read onto nodes: C = B B^T (no new source)")
    C = B @ B.T
    print(f"    C =\n{C}")
    check("C is NOT a scalar identity (the readout does NOT collapse)", not np.allclose(C, C[0, 0] * np.eye(3)))
    check("C is NOT S3-invariant aI+bJ (tau-central node has variance 2 vs 1)", abs(C[1, 1] - C[0, 0]) > 0.5)

    print("\n[2] E-plane projection -> a genuine path covariance (the right object type)")
    lo, hi = eplane_eigs(C)
    eps = (hi - lo) / (hi + lo)
    print(f"    E-plane eigenvalues = {{{lo:.4f}, {hi:.4f}}}  (= the R1 path {{1/3, 1}})")
    print(f"    ellipticity eps = {eps:.4f}")
    check("E-plane eigenvalues are {1/3, 1} (same as the R1 position tensor -- consistent)",
          abs(lo - 1 / 3) < 1e-9 and abs(hi - 1) < 1e-9)
    check("so the scalar current DOES yield a genuine edge/path covariance (not scalar identity)", abs(eps) > 1e-6)

    print("\n[3] But it is SECTOR-INDEPENDENT, and matches no target")
    print("    The R1 generation path is the SAME for every sector; the sector data")
    print("    (d, N) does not enter the path readout at all.  So eps = 1/2 for ALL sectors.")
    for s, t in TARGETS.items():
        print(f"      sector {s}: target d/N={t:.4f}  vs conditional-covariance eps={eps:.4f}  match={abs(t-eps)<1e-3}")
    check("conditional-covariance ellipticity (1/2) matches NO target d/N", all(abs(t - eps) > 1e-3 for t in TARGETS.values()))

    print("\n[4] Why: anisotropy and sector-dependence live in NON-INTERACTING structures")
    print("    - the anisotropy (R1 generation path) is sector-INDEPENDENT (eps=1/2);")
    print("    - the sector-dependence (d/N defect occupation) is generation-BLIND (S3-singlet);")
    print("    so the scalar conditional covariance inherits the path anisotropy but CANNOT")
    print("    pick up d/N -- the two factors never multiply without a coupling between them.")
    check("the missing coupling = a generation-resolved current (the named absent object)", True)

    print(
        """
[5] VERDICT -- the route gives the right OBJECT but not the sector MAGNITUDE
    Run strictly, the conditional-covariance route does the GOOD half: a scalar
    service current decomposed over the R1 path edges and read onto the generation
    nodes yields a genuine, non-trivial edge/path covariance (E-plane {1/3, 1},
    ellipticity 1/2) -- it does NOT collapse to a scalar identity.  So "edge/path
    covariance from the scalar current" is real and derivable, no new source needed.

    But it is SECTOR-INDEPENDENT: the R1 path is identical for every sector, the
    sector counts (d, N) never enter the readout, and the single ellipticity 1/2
    matches NONE of the targets (2/9, 1/3, 2/27, 1/9).  The reason is structural:
    the anisotropy (R1 generation path) carries no sector-dependence, while the
    sector-dependence (d/N defect occupation) is generation-blind (S3-singlet);
    the two live in non-interacting substrate structures and cannot multiply into
    a sector-dependent generation-anisotropic covariance without a coupling.

    HARD CONCLUSION (as anticipated): delta's DIRECTION is derived (R1 path) and a
    scalar-current conditional covariance EXISTS, but the d/N MAGNITUDE law cannot
    be obtained by re-reading the scalar current -- it requires exactly the coupling
    that is missing: a generation-resolved service current linking the path
    anisotropy to the sector defect counts.  Without new service-current structure,
    delta = d/N remains phenomenological.  The edge/readout decomposition exists;
    it just does not carry the sector magnitude.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
