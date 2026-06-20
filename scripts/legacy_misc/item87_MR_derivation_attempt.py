#!/usr/bin/env python3
r"""ITEM 87 — can the complex Majorana M_R (the delta_nu phase) be DERIVED from boot/walk dynamics?
Answer: NO, and the obstruction is EXACT. The delta_nu phase is a geometric primitive, not a
dynamical output of the documented walk.

SETUP. item87_majorana_cp_operator.py showed the missing CP operator IS a complex-symmetric Majorana
M_R = the framework's delta_nu=1/3 neutrino structure. The remaining item-87 gap is to DERIVE that
complex M_R (its phase) from the boot, not posit it (item 53 supplies only the scale). The natural
boot origin is the Delta L=2 Feshbach coupling of nu_R to its complement (the conjugate 0xFF^c, which
lives in the invalid/Q subspace), walk-dressed -- the phase coming from the walk's e^{ik pi/4}.

EXACT NO-GO (verified below). The all-zeros electron-neutrino nu_e = 0x00 is an EXACT EIGENSTATE of
the walk: W|0x00> = (sum_k A_k)|0x00> with sum_k A_k = sqrt(7/9) - sqrt(2/63) = 0.70374... REAL, and
the off-diagonal column W[1:,0] is EXACTLY zero -- because every CNOT control bit is 0 in the
all-zeros register, so no flip ever fires on it. It remains an eigenstate of (W^dag W) and of W^n for
all n. Therefore:
  * nu_e (Gen1) NEVER mixes with the other neutrino generations under the boot/walk dynamics;
  * the walk-built neutrino sector is forever effectively <=2-generation with Gen1 decoupled, so it
    carries NO 3-generation CP and NO inter-generation phase on Gen1;
  * the Delta L=2 walk-through-Q amplitude <0xFF^reg_j|W^n|reg_i> is moreover NON-SYMMETRIC (invalid
    as a Majorana mass) and has its Gen1 row identically zero (nu_e inert), giving CP = 0.

CONSEQUENCE. The complex symmetric M_R that carries delta_nu CANNOT be a dynamical output of the
documented walk: (i) the walk is Hermitian / Delta L=0 (cannot make a Delta L=2 symmetric mass), and
(ii) even setting that aside, nu_e is an inert walk eigenstate so no walk-dressed process gives Gen1
an inter-generation phase. The delta_nu phase is therefore a GEOMETRIC PRIMITIVE -- the defect
fraction d/N = 3/9 (and the nu_R Berry phase 2/9, item 86) -- NOT a dynamical product of the walk.
The only way to DERIVE it dynamically would be a genuinely new UNCONDITIONAL (non-CNOT-gated) Delta L=2
coupling that moves the all-zeros nu_e -- which is, by construction, outside the documented substrate.

NET for item 87: the "ansatz -> derived" gap is RESOLVED in the negative for the walk-dynamics route
(exact no-go: nu_e is a walk eigenstate). delta_nu stands as a geometric/combinatorial primitive
(d/N=3/9), at the same tier as the framework's other defect-fraction inputs; a dynamical derivation
requires new Delta L=2 physics, not more cleverness with the existing walk.

Self-asserting; exit 0 = nu_e is an exact (real-eigenvalue) walk eigenstate, robust across powers and
to W^dag W; and the Delta L=2 walk-through-Q M_R attempt is non-symmetric + Gen1-inert + CP=0.
"""
from __future__ import annotations
import numpy as np


def walk(sign=1):
    dlt = 2 / 9
    A = np.zeros(8, complex); A[0] = np.sqrt(1 - dlt)
    for k in range(1, 8):
        A[k] = np.sqrt(dlt / 7) * np.exp(sign * 1j * k * np.pi / 4)
    W = np.zeros((256, 256), complex)
    for k in range(8):
        c, t = (2 - k) % 8, (5 - k) % 8
        for i in range(256):
            W[(i ^ (1 << t)) if (i >> c) & 1 else i, i] += A[k]
    return W, A


def check(cond, msg):
    print(f"  [{'PASS' if cond else 'FAIL'}] {msg}")
    if not cond:
        raise AssertionError(msg)


def main() -> int:
    print("ITEM 87 — is the delta_nu phase derivable from boot/walk dynamics? (exact no-go)")
    print("=" * 88)
    W, A = walk()

    # ---- [1] nu_e = 0x00 is an EXACT walk eigenstate (all CNOT controls absent in all-zeros) ----
    print("\n[1] nu_e = 0x00: exact eigenstate of the walk?")
    col0 = W[:, 0]
    lam = col0[0]
    print(f"    W|0x00> nonzero entries = {int(np.count_nonzero(np.abs(col0) > 1e-14))};  "
          f"eigenvalue = {lam.real:.6f}{lam.imag:+.1e}i;  sum_k A_k = {np.sum(A).real:.6f}")
    check(np.max(np.abs(col0[1:])) < 1e-12, "W|0x00> = lambda|0x00> exactly (off-diagonal column = 0)")
    check(abs(lam.imag) < 1e-12 and abs(lam - np.sum(A)) < 1e-12,
          "the eigenvalue is REAL = sum_k A_k = sqrt(7/9)-sqrt(2/63) (no phase on nu_e)")
    G = W.conj().T @ W
    check(np.max(np.abs(G[1:, 0])) < 1e-12, "0x00 is also an exact eigenstate of the mass operator W^dag W")
    W12 = np.linalg.matrix_power(W, 12)
    check(int(np.count_nonzero(np.abs(W12[:, 0]) > 1e-12)) == 1,
          "W^n|0x00> stays purely |0x00> for all n -> nu_e NEVER mixes with other generations")

    # ---- [2] the Delta L=2 walk-through-Q Majorana attempt: non-symmetric + Gen1-inert + CP=0 ----
    print("\n[2] Delta L=2 boot Majorana M_R[i,j] = <0xFF^reg_j | W^n | reg_i> (walk-dressed complement):")
    regs = [0, 2, 1]                                   # nu_R Gen1,Gen2,Gen3 codewords
    Wn = np.linalg.matrix_power(W, 10)
    M = np.array([[Wn[0xFF ^ regs[j], regs[i]] for j in range(3)] for i in range(3)])
    H = np.diag([0.5, 0.8, 1.3]) ** 2
    Md = M.conj().T
    cp = float(np.imag(np.trace(H @ Md @ M @ Md @ H @ M)))
    asym = float(np.max(np.abs(M - M.T)))
    gen1row = float(np.max(np.abs(M[0, :])))
    print(f"    max|M| = {np.max(np.abs(M)):.2e};  non-symmetry max|M-M^T| = {asym:.2e};  "
          f"Gen1 row max = {gen1row:.2e};  CP invariant = {cp:+.2e}")
    check(asym > 0.4 * np.max(np.abs(M)), "the walk-through-Q amplitude is NON-symmetric -> invalid Majorana mass")
    check(gen1row < 1e-12, "Gen1 (nu_e) row is identically zero -> nu_e inert (consistent with [1])")
    check(abs(cp) < 1e-9, "CP invariant = 0 (Gen1 decoupled -> effectively 2-generation -> no CP)")

    print("\n" + "=" * 88)
    print("VERDICT — delta_nu is NOT derivable from the walk; it is a geometric primitive (exact no-go)")
    print(
        "  The all-zeros nu_e (0x00) is an EXACT walk eigenstate with a REAL eigenvalue (sum_k A_k =\n"
        "  0.70374), because every CNOT control is absent in the all-zeros register -- so no flip ever\n"
        "  fires on it. It stays an eigenstate of W^dag W and of W^n for all n: nu_e NEVER mixes with\n"
        "  the other neutrino generations under the boot/walk. Hence the walk-built neutrino sector is\n"
        "  permanently <=2-generation (Gen1 decoupled) and carries no 3-generation CP; and the natural\n"
        "  Delta L=2 walk-through-Q Majorana attempt is non-symmetric (invalid) with a zero Gen1 row\n"
        "  and CP=0. So the complex symmetric M_R that carries delta_nu CANNOT be a dynamical output of\n"
        "  the documented walk (which is also Hermitian / Delta L=0 and so cannot make a Delta L=2 mass\n"
        "  at all). delta_nu therefore stands as a GEOMETRIC PRIMITIVE -- the defect fraction d/N=3/9\n"
        "  (and the nu_R Berry phase 2/9, item 86) -- not a derived dynamical quantity. A dynamical\n"
        "  derivation would require a NEW unconditional (non-CNOT-gated) Delta L=2 coupling that moves\n"
        "  the all-zeros nu_e -- by construction outside the documented substrate.\n"
        "  ITEM 87 RESOLUTION: the 'ansatz -> derived' gap is closed in the NEGATIVE for the walk route\n"
        "  (exact no-go). The CP operator (complex M_R) exists and works (prior result); its phase\n"
        "  delta_nu is a geometric primitive (d/N=3/9), at the tier of the framework's other defect-\n"
        "  fraction inputs -- not walk-derivable. TIER: exact no-go for the walk-derivation; delta_nu\n"
        "  reclassified as geometric primitive; a dynamical origin needs new Delta L=2 physics."
    )
    print("exit 0 — nu_e is an exact walk eigenstate; delta_nu not walk-derivable, it is a geometric primitive.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
