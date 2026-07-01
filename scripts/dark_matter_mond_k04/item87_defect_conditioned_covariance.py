#!/usr/bin/env python3
r"""Item 87 -- the user's Defect-Conditioned R1 Covariance proposal.

Compact form (NO new generation-resolved source; generation structure appears only
in the 2nd moment of the readout, after conditioning on defect occupation):

    C_s = d_s * K_R1 + ((N_s - 2 d_s)/3) * I
          \_____/        \____________/
        defect-polarised   isotropic completion of the
        path traffic       scalar service ledger

with K_R1 = B B^T the R1 path covariance (E-plane spectrum {1, 1/3}).  In the
ellipticity normalization eps = A/(2A+3B), setting A_s=d_s, B_s=(N_s-2d_s)/3 gives
eps = d_s/N_s for all four sectors.

This script (1) verifies that; (2) checks the construction is physically consistent
(B_s>=0); (3) separates EARNED structure from ASSUMED premises; (4) flags honestly
that the four-way agreement is TAUTOLOGICAL (A,B were set to produce it), so only a
ledger derivation of the premises -- not the agreement -- can close delta.
"""
from fractions import Fraction as F
import numpy as np

# R1 path incidence (nodes e,tau,mu x edges e-tau, tau-mu) and path covariance
B = np.array([[1, 0], [1, 1], [0, 1]], float)
K_R1 = B @ B.T  # = [[1,1,0],[1,2,1],[0,1,1]]
SECTORS = {"e": (2, 9), "nu": (3, 9), "down": (1, 9), "up": (2, 27)}  # (d, N)


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def eplane_eigs(C):
    P = np.eye(3) - np.ones((3, 3)) / 3.0
    return np.sort(np.linalg.eigvalsh(P @ C @ P))[1:]


def ellipticity(C):
    lo, hi = eplane_eigs(C)
    return (hi - lo) / (hi + lo)


def main():
    print("ITEM 87 -- defect-conditioned R1 covariance: C_s = d K_R1 + (N-2d)/3 I")
    print("=" * 80)

    print("\n[1] The construction reproduces eps = d/N for all four sectors")
    for s, (d, N) in SECTORS.items():
        C = d * K_R1 + F(N - 2 * d, 3).__float__() * np.eye(3)
        eps = ellipticity(C)
        tgt = d / N
        print(f"    {s:<5s} d={d} N={N:2d}: eps={eps:.5f}  target d/N={tgt:.5f}  match={abs(eps-tgt)<1e-9}")
        check(f"{s}: eps == d/N", abs(eps - tgt) < 1e-9)

    print("\n[2] Physical consistency: the isotropic completion B_s=(N-2d)/3 must be >= 0")
    for s, (d, N) in SECTORS.items():
        Bs = F(N - 2 * d, 3)
        print(f"    {s:<5s}: N-2d = {N-2*d} >= 0 -> B_s = {Bs} >= 0  ({'OK' if Bs >= 0 else 'UNPHYSICAL'})")
        check(f"{s}: N_s >= 2 d_s (defects don't over-consume the service budget)", N >= 2 * d)

    print("\n[3] EARNED: the '2' in (N-2d) is the R1 path structure, not a fudge")
    lo, hi = eplane_eigs(K_R1)
    ratio = (hi + lo) / (hi - lo)  # E-plane trace/gap of K_R1
    deg_tau = int(K_R1[1, 1])      # central-node degree = incident edges
    print(f"    K_R1 E-plane spectrum = {{{lo:.4f}, {hi:.4f}}};  trace/gap = {ratio:.4f}")
    print(f"    R1 path central-node (tau) degree = {deg_tau} incident edges")
    check("the '2' = K_R1 E-plane trace/gap ratio (forced by the R1 path spectrum)", abs(ratio - 2) < 1e-9)
    check("...and equals the path central-node degree (2 edges) -- same 2-edge structure", deg_tau == 2)

    print("\n[4] HONEST: eps depends ONLY on the ratio B/A -> the content is a RATIO law")
    # ellipticity is scale-invariant: scaling (A,B)->(cA,cB) leaves eps unchanged.
    for s, (d, N) in SECTORS.items():
        A, Bs = d, F(N - 2 * d, 3).__float__()
        eps1 = ellipticity(A * K_R1 + Bs * np.eye(3))
        eps2 = ellipticity((5 * A) * K_R1 + (5 * Bs) * np.eye(3))  # same ratio, 5x scale
        check(f"{s}: eps is scale-invariant (depends only on B/A = (N-2d)/(3d))", abs(eps1 - eps2) < 1e-12)
    print("    So for delta=eps the content is exactly  B_s/A_s = (N_s - 2 d_s)/(3 d_s).")
    print("    A_s = d_s as an ABSOLUTE amplitude is extra (it would set the mass SCALE,")
    print("    testable via masses, not via eps).")

    print(
        """
[5] VERDICT -- promising, R1-grounded, NOT closed (and the agreement is not evidence)
    The construction is internally consistent (B_s>=0 for all sectors; N_s>=2 d_s
    holds non-trivially) and its single non-trivial coefficient -- the '2' in
    (N-2d) -- is EARNED by the R1 path: it is the path E-plane trace/gap ratio,
    equivalently the central-node (tau) degree of 2 edges.  So the user's
    "two record contacts per path-polarised defect" is the same 2-edge path
    structure, not a free parameter.

    But for delta=eps the predictive content is only the RATIO
        B_s/A_s = (N_s - 2 d_s)/(3 d_s),
    and the four-way agreement eps=d/N is TAUTOLOGICAL (A,B were chosen to give it),
    so it is NOT evidence.  What WOULD close delta is a ledger derivation of the
    remaining premises, none of which is in canon yet:
      (a) defect traffic is R1-path-polarised with variance proportional to d_s
          (independent-defect variance additivity -> the A_s ~ d_s term);
      (b) the total scalar service budget per sector is N_s;
      (c) the non-defect remainder (N_s - 2 d_s slots) completes ISOTROPICALLY.
    The crux is (a)+(c) together: that DEFECT occupation conditions the scalar
    service current onto the R1 path while the rest averages isotropically -- a
    defect-conditioned readout, not a new source port (consistent with canon's
    S3-singlet source).  This is a sharp, finite, falsifiable theorem target, with
    the '2' already earned.  The delta thread is NOT exhausted: it is reduced to the
    Defect-Conditioned R1 Covariance Theorem on the record-action/QEC ledger.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
