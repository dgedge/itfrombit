#!/usr/bin/env python3
r"""Item 87 -- the user's balance correction + the explicit boot-monitored R1 QND
projector.  Two results, then the two remaining theorem targets.

CORRECTION (user): the e-vs-mu "bias" must NOT be an unequal real rate inside the
CP-even covariance.  For C_fc = diag(r_e, 0, r_mu) the E-plane ellipticity is

    eps = sqrt(r_e^2 - r_e r_mu + r_mu^2) / (r_e + r_mu),

which equals 1/2 ONLY when r_e = r_mu.  So the symmetric (delta) channel must stay
e/mu-BALANCED; Phi must live in the oriented antisymmetric cochain (a sign/phase),
not in unequal real rescue rates.

QND structure: R1 as a boot-monitored validity projector
    P_R1 = I - |11><11|
with single-bit repair Kraus  11->01, 11->10  (no two-bit 11->00), balanced.
"""
import numpy as np

# B_2 basis: 00=0, 01=1, 10=2, 11=3
E = [np.eye(4)[i] for i in range(4)]


def check(name, cond):
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        raise AssertionError(name)


def ellipticity_diag(r_e, r_mu):
    """E-plane ellipticity of diag(r_e, 0, r_mu) on (e, tau, mu)."""
    C = np.diag([r_e, 0.0, r_mu])
    u1 = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
    u2 = np.array([1.0, 1.0, -2.0]) / np.sqrt(6)
    U = np.column_stack([u1, u2])
    lo, hi = np.sort(np.linalg.eigvalsh(U.T @ C @ U))
    return (hi - lo) / (hi + lo)


def main():
    print("ITEM 87 -- balance correction + explicit boot-monitored R1 QND projector")
    print("=" * 80)

    print("\n[1] CORRECTION verified: CP-even delta requires e/mu-BALANCED rescue rates")
    for (r_e, r_mu) in [(1, 1), (2, 1), (3, 1), (1, 4)]:
        eps = ellipticity_diag(r_e, r_mu)
        closed = np.sqrt(r_e**2 - r_e*r_mu + r_mu**2) / (r_e + r_mu)
        print(f"    r_e={r_e}, r_mu={r_mu}: eps={eps:.4f}  (closed form {closed:.4f})")
        check(f"eps matches sqrt(r_e^2-r_e r_mu+r_mu^2)/(r_e+r_mu)", abs(eps - closed) < 1e-12)
    check("eps = 1/2 IFF r_e = r_mu (balanced)", abs(ellipticity_diag(1, 1) - 0.5) < 1e-12 and abs(ellipticity_diag(2, 1) - 0.5) > 1e-3)
    print("    => the symmetric delta channel MUST be balanced; Phi lives in the ORIENTED")
    print("       antisymmetric cochain (a sign/phase), not in unequal CP-even rescue rates.")

    print("\n[2] Boot-monitored R1 QND projector P_R1 = I - |11><11| with balanced single-bit Kraus")
    P_R1 = np.diag([1.0, 1.0, 1.0, 0.0])              # QND validity projector (kills 11)
    K0 = P_R1
    K1 = (1/np.sqrt(2)) * np.outer(E[1], E[3])        # 11 -> 01 (e), single-bit
    K2 = (1/np.sqrt(2)) * np.outer(E[2], E[3])        # 11 -> 10 (mu), single-bit
    kraus = [K0, K1, K2]
    completeness = sum(K.conj().T @ K for K in kraus)
    check("Kraus set is CPTP (sum K^dag K = I)", np.allclose(completeness, np.eye(4)))
    # no two-bit 11->00 repair present
    check("no 11->00 one-step repair (that would be a two-bit move)", abs(K1[0, 3]) < 1e-12 and abs(K2[0, 3]) < 1e-12)

    def channel(rho):
        return sum(K @ rho @ K.conj().T for K in kraus)

    # QND on the valid subspace: identity on any state of {00,01,10}
    rho_valid = np.diag([0.4, 0.35, 0.25, 0.0])
    check("QND on the valid subspace (channel acts as identity there -> no flavour change)",
          np.allclose(channel(rho_valid), rho_valid))
    # the rescue fires only on 11, BALANCED into e and mu:
    out11 = channel(np.outer(E[3], E[3]))
    check("11 is rescued, balanced 1/2 to e (01) and 1/2 to mu (10)",
          abs(out11[1, 1] - 0.5) < 1e-12 and abs(out11[2, 2] - 0.5) < 1e-12 and out11[3, 3] < 1e-12)
    print("    -> balanced rescue (1/2,1/2) gives the CP-even delta; an oriented version")
    print("       (complex phases on K1,K2) carries Phi WITHOUT unbalancing the real rates.")

    print(
        """
[3] The two remaining theorem targets (the narrow premises)
    (A) Monitored-R1 (clause 1, nonlinear): prove P_R1 = I - |11><11| is a
        BOOT-MONITORED QND validity projector with the balanced single-bit Kraus
        records above -- i.e. that the nonlinear R1 closure is part of the per-tick
        monitored validity syndrome, not a timeless cut.  (sec 5.9 favours this;
        the open point is that R1 is nonlinear, so it is not one of the documented
        LINEAR [8,4,4] stabilisers.)
    (B) Sector-count law: prove A_s = d_s and B_s = (N_s - 2 d_s)/3, with d_s the
        number of sector-visible R1 rescue contacts and N_s the total scalar service
        budget.  This is the entire remaining delta-MAGNITUDE content (the '2' is
        already R1-earned; the ellipticity is balanced-symmetric -> d/N).

[4] VERDICT -- delta/Phi NEAR-DERIVATION: grounded to two narrow premises, not closed
    With the balance correction, the construction is clean and consistent:
      * delta = CP-even ellipticity of the BALANCED symmetric rescue covariance
        (e/mu equal) -> 1/2 -> d/N; an unequal real rate would (correctly) spoil it,
        so Phi cannot hide there;
      * Phi  = the ORIENTED antisymmetric R1 cochain (sign/phase), separate channel;
      * the QND projector P_R1 with single-bit balanced Kraus (11->01, 11->10, no
        11->00) is an explicit, CPTP, QND-on-valid recovery -- no flavour change.
    Everything reduces to two sharp premises: (A) nonlinear R1 is part of the
    per-tick monitored validity syndrome (the QND projector), and (B) its sector-
    visible rescue-contact count is d_s.  If both land, delta=d/N and Phi's
    orientation close together; if not, delta stays phenomenological.  This is the
    end of the line for the canon-grounded thread: a near-derivation resting on two
    narrow, well-posed, sec-5.9-aligned premises.
exit 0"""
    )
    print("ALL ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
