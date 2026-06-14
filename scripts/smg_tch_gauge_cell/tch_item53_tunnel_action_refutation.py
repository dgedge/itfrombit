#!/usr/bin/env python3
"""
Refute/prove candidate item-53 lemma:

    S_tunnel = (# nonzero H_PQ channels) x d_code = 6 x 4 = 24.

Verdict encoded here:
  The strong lemma is false for the actual Part 18 Feshbach operator.

Reason:
  The six H_PQ entries are not six barriers in series. They are parallel
  matrix elements. Feshbach elimination uses

      Sigma(E) = H_PQ (E I - H_QQ)^-1 H_QP,

  so the channel count contributes to ranks/eigenvalue prefactors, not as an
  additive exponent 6*d_code for a single mass eigenvalue.

The only honest place where "6*d_code" appears is the determinant/product of
the three nonzero self-energy eigenvalues under an added per-channel
attenuation. That is not the individual neutrino mass scale.
"""

from __future__ import annotations

import math


D_CODE = 4
N_Q = 3
CHANNELS_PER_Q = 2
N_CHANNELS = N_Q * CHANNELS_PER_Q


def coupling_matrix(g: float = 1.0) -> list[list[float]]:
    """
    Minimal Part 18 coupling matrix A = H_PQ restricted to the six active
    lepton states that couple to the three nu_R pseudocodewords.

    Rows:
      e_R(0), nu_L(0), e_R(1), nu_L(1), e_R(2), nu_L(2)
    Columns:
      nu_R(0), nu_R(1), nu_R(2)

    Each nu_R has two parallel unit couplings.
    """

    A = [[0.0 for _ in range(N_Q)] for _ in range(N_CHANNELS)]
    for q in range(N_Q):
        A[2 * q][q] = g
        A[2 * q + 1][q] = g
    return A


def spectra(g: float = 1.0) -> tuple[list[float], list[float], list[float]]:
    """
    Return singular values of H_PQ, Q-space H_QP H_PQ eigenvalues, and
    P-space H_PQ H_QP eigenvalues.
    """

    # Each Q column has two entries equal to g and is orthogonal to the other
    # columns. Therefore H_QP H_PQ = 2 g^2 I_3 exactly.
    singular = [math.sqrt(2.0) * abs(g)] * N_Q
    q_eigs = [2.0 * g * g] * N_Q
    p_eigs = [0.0] * (N_CHANNELS - N_Q) + [2.0 * g * g] * N_Q
    return singular, q_eigs, p_eigs


def action(x: float) -> float:
    return -math.log(abs(x))


def fmt(values: list[float]) -> str:
    return "[" + " ".join("%.12g" % v for v in values) + "]"


def main() -> None:
    print("=" * 78)
    print("Item-53 tunnel-action lemma audit")
    print("=" * 78)
    print("Candidate: S_tunnel = (# nonzero H_PQ channels) x d_code = 6 x 4 = 24")
    print()

    singular, q_eigs, p_eigs = spectra(1.0)
    nonzero_p = [x for x in p_eigs if abs(x) > 1e-12]
    print("A. Actual Part 18 coupling algebra")
    print("   H_PQ shape in coupled subspace       = 6 x 3")
    print("   nonzero entries                      = %d" % N_CHANNELS)
    print("   singular values                      = %s" % fmt(singular))
    print("   eigenvalues of H_QP H_PQ             = %s" % fmt(q_eigs))
    print("   nonzero eigenvalues of H_PQ H_QP     = %s" % fmt(nonzero_p))
    print("   interpretation: three rank-one columns, each with two parallel")
    print("   couplings. The six entries do not multiply in series.")
    print()

    g = math.exp(-D_CODE)
    singular_g, q_eigs_g, p_eigs_g = spectra(g)
    nonzero = [x for x in q_eigs_g if abs(x) > 1e-300]
    eig_action = action(nonzero[0])
    trace_action = action(sum(nonzero))
    det_value = 1.0
    for value in nonzero:
        det_value *= value
    det_action = action(det_value)
    geom_action = det_action / len(nonzero)
    print("B. If each one-bit H_PQ amplitude is assigned action d_code = 4")
    print("   coupling amplitude g = exp(-4)")
    print("   self-energy eigenvalue per generation = 2 g^2")
    print("   action per nonzero eigenvalue         = -ln(2 exp(-8))")
    print("                                        = 2*d_code - ln 2")
    print("                                        = %.6f" % eig_action)
    print("   trace action over all channels        = -ln(6 exp(-8)) = %.6f" % trace_action)
    print("   determinant/product action            = -ln((2 exp(-8))^3)")
    print("                                        = 6*d_code - ln 8")
    print("                                        = %.6f" % det_action)
    print("   geometric-mean eigenvalue action      = %.6f" % geom_action)
    print()

    print("C. Comparison to the proposed S=24")
    print("   proposed S=24 is not:")
    print("     eigenvalue action       %.6f" % eig_action)
    print("     trace action            %.6f" % trace_action)
    print("     determinant action      %.6f  (even this has the ln 8 prefactor)" % det_action)
    print("     geometric mean action   %.6f" % geom_action)
    print()

    print("D. Algebraic proof of refutation")
    print("   With per-entry attenuation g, the actual matrix is:")
    print("       H_QP H_PQ = 2 g^2 I_3")
    print("   Therefore each physical nonzero Feshbach eigenvalue scales as 2 g^2,")
    print("   not as g^6. In exponent language:")
    print("       S_eigen = 2 S_entry - ln 2")
    print("   Setting S_entry = d_code gives S_eigen = 2*d_code - ln 2, not 6*d_code.")
    print("   The factor 6 appears only in counts/traces or in the determinant product,")
    print("   neither of which is the light-neutrino mass eigenvalue scale.")
    print()

    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    print("  Refuted as stated for the actual Part 18 Feshbach operator.")
    print("  The six H_PQ entries are parallel channels. Feshbach projection squares")
    print("  couplings and sums parallel paths; it does not add six code-distance")
    print("  barriers into one tunnelling exponent.")
    print()
    print("  Weaker surviving possibility:")
    print("    A separate collective instanton/stabilizer-action rule could multiply")
    print("    the whole neutrino block by exp(-24). But that is an additional")
    print("    dimensionalisation postulate, not a consequence of H_PQ as written.")


if __name__ == "__main__":
    main()
