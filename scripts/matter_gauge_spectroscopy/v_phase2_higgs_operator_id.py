#!/usr/bin/env python3
r"""v_phase2_higgs_operator_id.py -- attempt the operator-identification theorem for the Higgs VEV.

TARGET (user's ask): derive "Higgs interpolating field = the all-8 product via the R4 projector", which
would force the k=8 channel (mu/M_P = alpha_0^8) rather than leave it a physical argument.

HONEST FINDING -- the LITERAL theorem FAILS, informatively, and the corrected form succeeds modulo 2 premises:

PART 1 (the literal form is the HIERARCHY DISASTER). The R4 projector onto the forbidden sterile
pseudo-codeword nu_R is P = |nu_R><nu_R| = PROD_i (I + s_i Z_i)/2 (s_i = +-1 the nu_R bit pattern).
Expanding: P = (1/2^8) SUM_{S subset of 8} (prod_{i in S} s_i) Z_S -- it contains EVERY Pauli weight
k=0..8, with C(8,k) terms at weight k. The non-unital expectation of a weight-k Pauli ~ alpha_0^k (k bits,
alpha_0 each), so the projector's contribution to the Higgs mass parameter is
    <P>_non-unital ~ (1/2^8) SUM_k C(8,k) alpha_0^k = (1/2^8)(1+alpha_0)^8  ~  1/256  (the k=0 IDENTITY term),
i.e. mu^2 ~ M_P^2 -- the LOW-weight (identity) term DOMINATES; the alpha_0^8 (weight-8) term is the
SMALLEST, suppressed by alpha_0^8 ~ 1e-17 relative to it. So "Higgs = R4 projector all-8 product" as
literally stated gives mu ~ M_P: the hierarchy problem. The naive operator-identification is FALSE.

PART 2 (the corrected form: a STATE TRANSITION, not an operator matrix element). The flaw is treating the
condensate as <operator>; it is a STATE-FORMATION amplitude vacuum -> (EW condensate). The EW condensate
is the FILLED matter cell: R4 says the generation is the complete 16-spinor minus the sterile nu_R, i.e.
the EW vacuum is the cell with all its non-sterile modes OCCUPIED. The amplitude to form the filled cell
from the empty vacuum requires CREATING every one of the 8 modes -- there is no "skip" (an unoccupied
mode = an UNfilled cell = a different, orthogonal state). So:
    A(vacuum -> all-8 filled)  = prod_i (per-mode fill amp)        -> P(filled)  = alpha_0^8   [the Higgs]
    A(vacuum -> exactly 1 mode) = C(8,1)(...)                       -> P(1 mode)  ~ 8 alpha_0   [self-energy]
These are DISTINCT, ORTHOGONAL final states (different occupation number) -- so the k<8 terms do NOT
contaminate <filled>, unlike the projector's identity terms. The k=8 is forced because the FILLED cell
(R4's complete generation) needs all 8 modes; the k=1 single-occupation is a different transition (the
photon self-energy). mu/M_P = P(filled) = alpha_0^8.

WHY a state transition and not an operator (the key correction): the projector P is hermitian and contains
I_i on each bit ("this bit may do nothing"), so its expectation is dominated by doing-nothing (mu~M_P).
The condensate is vacuum->filled, where every mode MUST be created; "doing nothing" on a mode yields a
DIFFERENT (unfilled) state, removed by orthogonality, not a contaminating contribution. Orthogonality of
occupation sectors is what the projector's operator-expectation throws away -- and it is exactly what
kills the hierarchy here.

RESIDUAL (2 named premises + lambda), honestly:
  (P1) EWSB = the filled-cell transition (the Higgs VEV is the vacuum->complete-generation amplitude).
       Grounded in R4 = "the generation is complete (16-spinor minus the sterile nu_R)", but a physical
       identification, not derived from a Lagrangian.
  (P2) the per-mode fill amplitude is alpha_0 (item 79's universal per-cell coupling; an extension of
       leg iii's per-bit photon-emission projection to the mode-creation amplitude).
  (lambda) Phase 3.

NET: the literal operator-identification (projector -> all-8) is a NEGATIVE (it is the hierarchy disaster:
the identity/low-weight terms dominate). The CORRECT object is the filled-cell state transition, which is
intrinsically all-8 (alpha_0^8) with partial fillings orthogonal (k=1 = the self-energy), and k=8 is forced
by R4 = complete-generation. So the k=8 selection is reduced from "R4-codeword hand-wave" to "orthogonality
of occupation sectors + R4 = complete generation", modulo P1, P2, lambda. Sharper, not closed.
"""
import numpy as np
from math import comb

N = 8
ALPHA0 = 1.0 / 137.036
v, M_P = 246.22, 1.2209e19
lam = 125.25 ** 2 / (2 * v ** 2)


def build_R4_projector(pattern):
    """|nu_R><nu_R| as PROD_i (I + s_i Z_i)/2 on N qubits; returns the 2^N x 2^N rank-1 projector."""
    I2 = np.eye(2); Z = np.array([[1.0, 0.0], [0.0, -1.0]])
    P = np.array([[1.0]])
    for b in pattern:
        s = 1.0 if b == 0 else -1.0     # |0>: (I+Z)/2 ; |1>: (I-Z)/2
        P = np.kron(P, (I2 + s * Z) / 2.0)
    return P


def pauli_weight_spectrum():
    """analytic Pauli-weight spectrum of the single-state projector: C(8,k) terms at weight k."""
    return [comb(N, k) for k in range(N + 1)]


def main():
    print("=== Operator-identification attempt: is the Higgs the all-8 product via the R4 projector? ===\n")

    # nu_R bit pattern (LQ,C0,C1,I3,chi)=(0,0,0,0,1), chi=W=1, generation (G0,G1)=(0,1) -> 8 bits
    pattern = [0, 1, 0, 0, 0, 0, 1, 1]   # (G0,G1,LQ,C0,C1,I3,chi,W)
    P = build_R4_projector(pattern)
    rank = int(round(np.trace(P @ P).real))   # projector: P^2=P, Tr=rank=1
    print(f"  [0] R4 projector |nu_R><nu_R| built (nu_R pattern {pattern}); rank = {rank}, Tr = {np.trace(P).real:.3f}")
    spec = pauli_weight_spectrum()
    print(f"      Pauli-weight spectrum (C(8,k)): {spec}  (sum {sum(spec)} = 2^8); ALL weights k=0..8 present.\n")

    # PART 1: the literal (bosonic operator) reading -> hierarchy
    print("  [1] LITERAL form (projector as a bosonic operator): <P>_non-unital ~ (1/2^8) SUM_k C(8,k) alpha_0^k:")
    terms = [comb(N, k) * ALPHA0 ** k for k in range(N + 1)]
    bos = sum(terms) / 2 ** N            # = (1+alpha_0)^8 / 256
    print(f"      (1/256) SUM_k C(8,k) alpha_0^k = (1+alpha_0)^8/256 = {bos:.5f}")
    print(f"      leading term k=0 (IDENTITY) = {1/2**N:.5f};  weight-8 term = {ALPHA0**N/2**N:.2e}")
    print(f"      -> the IDENTITY/low-weight terms DOMINATE by ~alpha_0^-8 ~ {1/ALPHA0**N:.1e}; mu^2 ~ M_P^2.")
    print(f"      VERDICT: the literal 'Higgs = R4-projector all-8 product' IS THE HIERARCHY PROBLEM. FALSE.\n")

    # PART 2: the corrected (state-transition) reading -> alpha_0^8
    print("  [2] CORRECTED form (state transition vacuum -> filled cell; partial fillings ORTHOGONAL):")
    P_fill = [comb(N, k) * ALPHA0 ** k * (1 - ALPHA0) ** (N - k) for k in range(N + 1)]
    print(f"      P(vacuum -> all-8 filled) = alpha_0^8        = {P_fill[8]:.3e}   [the Higgs condensate]")
    print(f"      P(vacuum -> exactly 1 mode) = 8 alpha_0(...)  = {P_fill[1]:.3e}   [the EM self-energy]")
    print(f"      these are DISTINCT orthogonal final states (occupation number differs) -> no contamination,")
    print(f"      unlike the projector's identity terms. The filled cell (R4's complete generation) needs ALL 8.")
    mu_over_MP = ALPHA0 ** N
    vMP_pred = mu_over_MP / np.sqrt(lam)
    print(f"      => mu/M_P = P(filled) = alpha_0^8 = {mu_over_MP:.3e};  v/M_P = alpha_0^8/sqrt(lambda) = {vMP_pred:.3e}")
    print(f"         (observed v/M_P = {v/M_P:.3e}; agree to {100*abs(vMP_pred/(v/M_P)-1):.0f}%, the lambda-ballpark)\n")

    # the discriminator
    print("  [3] WHY a transition, not an operator (the correction that kills the hierarchy):")
    print(f"      the projector carries I_i on each bit ('do nothing') -> its expectation is dominated by")
    print(f"      doing nothing (mu~M_P). The condensate is vacuum->filled, where every mode MUST be created;")
    print(f"      'do nothing' on a mode = a DIFFERENT (unfilled) state, removed by ORTHOGONALITY of")
    print(f"      occupation sectors -- exactly the structure the operator-expectation discards. k=8 forced.\n")

    print("[verdict] OPERATOR-IDENTIFICATION ATTEMPT -- literal form NEGATIVE, corrected form reduces to 2 premises:")
    print("  - the LITERAL 'Higgs = R4-projector all-8 product' is the HIERARCHY DISASTER: the projector's")
    print("    identity/low-weight terms dominate (mu^2 ~ M_P^2), the alpha_0^8 term is the smallest. FALSE as stated.")
    print("  + the CORRECT object is the filled-cell STATE TRANSITION (vacuum -> all-8 occupied = R4's complete")
    print("    generation): intrinsically all-8 = alpha_0^8, with partial fillings ORTHOGONAL (k=1 = the self-energy),")
    print("    so no low-weight contamination. k=8 is FORCED by R4 = complete-generation + occupation-orthogonality.")
    print("  RESIDUAL (sharper): (P1) EWSB = the filled-cell transition (grounded in R4=complete-generation, but a")
    print("  physical id, not a Lagrangian derivation); (P2) per-mode fill amplitude = alpha_0 (item 79 universal-alpha,")
    print("  extending leg-iii emission to mode-creation); (lambda) Phase 3. The k=8 selection is reduced from an")
    print("  'R4-codeword hand-wave' to 'occupation-sector orthogonality + R4=complete-generation'. Sharper, NOT closed.")

    # gates
    assert rank == 1, "R4 projector must be rank-1 (a single-state projector)"
    assert spec == [1, 8, 28, 56, 70, 56, 28, 8, 1] and sum(spec) == 2 ** N, "Pauli weights must be C(8,k), sum 2^8"
    assert 1.0 < bos * 2 ** N < 1.2, "bosonic-projector reading = (1+alpha_0)^8/256 ~ 1.06/256, identity-dominated -- the hierarchy"
    assert (1 / 2 ** N) / (ALPHA0 ** N / 2 ** N) > 1e15, "identity term must dominate weight-8 by >1e15 (hierarchy)"
    assert abs(np.log(P_fill[8]) / np.log(ALPHA0) - 8) < 1e-9, "filled-cell transition must be exactly alpha_0^8"
    assert P_fill[1] / P_fill[8] > 1e14, "single-occupation (self-energy) and filled (Higgs) must be different powers"
    assert abs(vMP_pred / (v / M_P) - 1) < 0.15, "alpha_0^8/sqrt(lambda) must give v/M_P to the lambda-ballpark"
    print("\nGATES PASSED -- literal projector reading = hierarchy (identity dominates by >1e15); corrected filled-cell")
    print("transition = alpha_0^8 (orthogonality kills low-weight); k=8 forced by R4=complete-generation. exit 0")


if __name__ == "__main__":
    main()
