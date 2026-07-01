#!/usr/bin/env python3
r"""foundations_g7_gksl_reduction_attempt.py

G7 — the LAST IDENTIFICATION attempt: derive the gravity erasure dissipator as the GKSL reduction of the
canonical walk W, and try to force gamma_k uniform (=> C_loop=3/2) from W directly.

Setup (foundations_g7_cloop_support_from_walk.py + DRIFT G7): C_loop = (sum gamma)/Gamma_vac is the gravity
O(1) coefficient; =3/2 IFF the erasure rates gamma_k are UNIFORM on the walk-active support {C0,C1,I3}.
Last turn: support DERIVED from §3.2; golden-rule billing (gamma ∝ coherent amplitude^2) gives 10/9 and was
"data-excluded at 26%". The open step the user asked to attempt: build the GKSL reduction of W and show its
dissipative part is the unital (uniform) channel on the support, NOT the golden-rule one.

RESULT (decisive, and it RESOLVES the route rather than closing the coefficient):

  The canonical coherent walk H_int = V_em + V_weak + V_strong is CODE-PRESERVING: Pi_Q H_int Pi_P = 0
  EXACTLY (verified). V_em is diagonal; V_weak (zero-controlled I3 flip) flips the FREE bit I3 -> stays in P;
  V_strong (colour permutation among {01,10,11}) never changes (C==00) -> never changes validity -> stays in P.
  So the coherent walk causes NO P->Q transitions: H_PQ = 0.

  CONSEQUENCE 1 -- the GKSL-from-W route is EMPTY. The second-order self-energy Sigma = H_PQ G_QQ H_QP = 0
  because H_PQ = 0. There is nothing to reduce: W is the protected LOGICAL dynamics; it does not generate
  errors. (This is the correct QEC picture: logical ops preserve the code; NOISE causes errors.)

  CONSEQUENCE 2 -- golden-rule is a CATEGORY ERROR, not merely data-excluded. Golden-rule bills the coherent
  amplitudes (weak sqrt(2/9), strong 1); but those amplitudes drive IN-CODE logical transitions, not errors.
  On the vacuum the golden-rule P->Q rate is identically 0 -- nonsensical (the code DOES correct colour
  errors). So the §3.2 coherent amplitudes are simply IRRELEVANT to gamma_k. The 10/9 worry is removed.

  CONSEQUENCE 3 -- the erasure dissipator is therefore the QEC NOISE/ERROR MODEL (item 119's
  L_k = Pi_Q X_k Pi_P with X_k the Pauli error), supported on the walk-active bits. Pauli-X jumps are
  Hermitian, so uniform-rate noise on {C0,C1,I3} is UNITAL (sum gamma X^dag X = sum gamma X X^dag ∝ I) --
  the SAME unital class item 79 proves gives the alpha_0=1/137 equipartition. Under that symmetric/unital
  noise, C_loop = 3/2 (support from W; uniformity from the unital error model).

VERDICT: the GKSL-from-W identification does NOT yield a parameter-free C_loop -- it shows the route is empty
(W is code-preserving), which (a) KILLS the golden-rule objection as a category error and (b) pins the one
remaining premise to the NOISE model: that the gravity-loop error model is the item-79 unital monitor
restricted to the walk-active support. That is an identification of two existing canon dissipators, not a
derivation from W; W cannot supply it (it makes no errors). So C_loop=3/2 stands on the item-79 unital noise
premise; the magnitude stays the T7 wall. Honest bottom line: further progress IS the noise-model
identification + new physics -- the last internal route (reduce-from-W) is now shown to be empty.

Self-asserting; exit 0. Run under ~/bin/py13_7 (numpy).
"""
import itertools
import numpy as np

ALPHA = 1 / 137.035999
NAMES = ["G0", "G1", "LQ", "C0", "C1", "I3", "chi", "W"]
B = {n: i for i, n in enumerate(NAMES)}
ALL = list(itertools.product((0, 1), repeat=8))
IDX = {c: i for i, c in enumerate(ALL)}


def ok(c, m):
    print(f"  [{'PASS' if c else 'FAIL'}] {m}")
    assert c, m


def valid(c):
    G0, G1, LQ, C0, C1, I3, chi, W = c
    return (not (G0 and G1)) and W == chi and ((LQ == 0) == (C0 == 0 and C1 == 0))


inP = np.array([valid(c) for c in ALL], bool)
PiP = np.diag(inP.astype(float)); PiQ = np.diag((~inP).astype(float))
assert inP.sum() == 48


def fl(c, *ns):
    d = list(c)
    for n in ns:
        d[B[n]] ^= 1
    return tuple(d)


def main():
    print("=== build the canonical COHERENT walk H_int = V_em + V_weak + V_strong (256-dim) ===")
    Vem = np.diag([c[B["I3"]] - 0.5 * (1 - c[B["LQ"]]) for c in ALL]).astype(complex)   # diagonal
    Vw = np.zeros((256, 256), complex); cw = np.sqrt(2 / 9)
    for i, c in enumerate(ALL):
        if c[B["chi"]] == 0:                                   # zero-controlled CNOT (§3.1) -> flip I3 when chi=0
            Vw[IDX[fl(c, "I3")], i] = cw
    Vs = np.zeros((256, 256), complex); colours = [(0, 1), (1, 0), (1, 1)]
    for i, c in enumerate(ALL):
        cc = (c[B["C0"]], c[B["C1"]])
        if cc in colours:
            for t in colours:
                if t != cc:
                    d = list(c); d[B["C0"]], d[B["C1"]] = t; Vs[IDX[tuple(d)], i] = 1.0   # K3 colour perm, g_s=1
    H = Vem + Vw + Vs

    print("[1] is the coherent walk CODE-PRESERVING?  (the GKSL route needs H_PQ != 0)")
    Hpq = PiQ @ H @ PiP
    ok(np.linalg.norm(Hpq) < 1e-12, f"||Pi_Q H Pi_P|| = {np.linalg.norm(Hpq):.2e}  ->  H_PQ = 0: the coherent walk PRESERVES the code space P")
    # also each piece individually
    ok(np.linalg.norm(PiQ @ Vem @ PiP) < 1e-12, "  V_em is diagonal (no P->Q)")
    ok(np.linalg.norm(PiQ @ Vw @ PiP) < 1e-12, "  V_weak flips the FREE bit I3 -> stays in P (no P->Q)")
    ok(np.linalg.norm(PiQ @ Vs @ PiP) < 1e-12, "  V_strong permutes colour (never changes C==00) -> stays in P (no P->Q)")

    print("\n[2] CONSEQUENCE: the GKSL-from-W self-energy Sigma = H_PQ G_QQ H_QP is identically ZERO")
    # mock resolvent G_QQ (any bounded op on Q); Sigma must vanish because H_PQ=0
    rng = np.random.default_rng(0); Gqq = PiQ @ rng.standard_normal((256, 256)) @ PiQ
    Sigma = (PiP @ H @ PiQ) @ Gqq @ (PiQ @ H @ PiP)
    ok(np.linalg.norm(Sigma) < 1e-10, f"||Sigma_GKSL|| = {np.linalg.norm(Sigma):.2e} -> the reduce-from-W route is EMPTY (W makes no errors)")

    print("\n[3] CONSEQUENCE: golden-rule rate from the COHERENT walk on the vacuum = 0 (category error)")
    vac = IDX[(0,) * 8]
    gr_vac = sum(abs(H[q, vac]) ** 2 for q in range(256) if not inP[q])   # Fermi-golden-rule P->Q weight
    ok(gr_vac < 1e-12, f"golden-rule Gamma_vac(coherent) = {gr_vac:.2e} = 0 -> golden-rule bills the WRONG object (the code DOES correct colour errors)")
    print("    => the §3.2 coherent amplitudes (weak sqrt(2/9), strong 1) are IRRELEVANT to the erasure rates gamma_k.")

    print("\n[4] the TRUE dissipator is the QEC NOISE model (Pauli errors X_k); uniform noise is UNITAL")
    def Xk(name):
        M = np.zeros((256, 256));
        for i, c in enumerate(ALL):
            M[IDX[fl(c, name)], i] = 1.0
        return M
    sup = ["C0", "C1", "I3"]; g = ALPHA / 3
    lhs = sum(g * (Xk(n).conj().T @ Xk(n)) for n in sup)
    rhs = sum(g * (Xk(n) @ Xk(n).conj().T) for n in sup)
    ok(np.allclose(lhs, rhs) and np.allclose(lhs, ALPHA * np.eye(256)),
       "uniform Pauli-X noise on {C0,C1,I3} is UNITAL (sum g X^dag X = sum g X X^dag = alpha I) -- item-79's unital class")

    # under that unital noise, C_loop = 3/2 (support from W; uniformity from the unital error model)
    exit_k = {n: (not inP[IDX[fl(ALL[vac], n)]]) for n in sup}
    Gamma_vac = sum(g for n in sup if exit_k[n]); C_loop = ALPHA / Gamma_vac
    ok(exit_k == {"C0": True, "C1": True, "I3": False}, f"vacuum exits: {exit_k} (colour prints, I3 free)")
    ok(abs(C_loop - 1.5) < 1e-12, f"=> C_loop = alpha/Gamma_vac = {C_loop:.4f} = 3/2 (K = {C_loop/ALPHA:.1f}) under unital noise")

    print("\n[verdict] G7 GKSL-from-W identification:")
    print("  - The route is EMPTY: the coherent walk W is CODE-PRESERVING (H_PQ=0), so Sigma_GKSL=0 -- W")
    print("    generates no errors to reduce. (Correct QEC picture: W is the protected logical dynamics.)")
    print("  - This KILLS golden-rule as a CATEGORY ERROR (its coherent P->Q rate is identically 0), removing")
    print("    the 10/9 worry: the §3.2 amplitudes do not bill gamma_k at all.")
    print("  - The erasure dissipator is the QEC NOISE model; uniform Pauli noise on {C0,C1,I3} is UNITAL")
    print("    (item-79's class) -> C_loop = 3/2. But uniformity is a NOISE-MODEL premise (= the item-79")
    print("    monitor restricted to the support), an identification of two canon dissipators -- W cannot")
    print("    supply it, because W makes no errors.")
    print("  - HONEST: C_loop=3/2 is NOT made parameter-free. The last internal route (reduce-from-W) is shown")
    print("    EMPTY; the residual is the noise-model identification + the T7 magnitude wall -> further progress")
    print("    is the item-79-monitor identification or new physics. Confirms the prior assessment.")
    print("  exit 0")


if __name__ == "__main__":
    main()
