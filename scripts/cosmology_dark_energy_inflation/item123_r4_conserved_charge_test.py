#!/usr/bin/env python3
"""
item123_r4_conserved_charge_test.py

The CMB third-peak closure (DRIFT M15) needs a CONSERVED, pressureless dust
charge of abundance omega ~ 0.096 (the "R4 zero-mode reservoir").  M15 contains
a tension: the "Noether-charge attempt" found the R4 repair dynamics DISSIPATIVE
(only lepton number conserved -> the 0.024 nu_R relic), while the later
"reservoir-lift" entries POSIT a conserved N_tot.  This script adjudicates it by
computing the exact conserved-charge structure of the documented R4 Lindbladian.

Documented R4 dynamics (M15): an immigration-death process per active cell,
    d<N>/dt = Gamma (x - <N>),   x = |g|/a0   (the MOND source),
with "finite observable nullspace consisting only of constants".  At homogeneous
recombination g=0 -> x=0.  We build this as a Lindbladian and compute:
  (a) its conserved charges (kernel of the Heisenberg generator), and
  (b) its steady state at x=0,
then contrast with a NON-dissipative mode (death channel removed) for which a
number charge IS conserved.

Result determines whether the dust charge is DERIVED from W=S.C or IMPORTED.
Self-asserting; numpy only.
"""
import sys
import numpy as np

D = 8                                   # Fock cutoff
_ok = True


def check(name, cond):
    global _ok
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    _ok = _ok and bool(cond)


# Fock operators
a = np.diag(np.sqrt(np.arange(1, D)), 1)          # annihilation
ad = a.conj().T                                    # creation
Nop = ad @ a                                        # number
I = np.eye(D)


def lindblad_super(H, jumps):
    """Schrodinger generator L (acts on vec(rho))."""
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(I, LdL) - 0.5 * np.kron(LdL.T, I))
    return L


def heisenberg_super(H, jumps):
    """Adjoint (Heisenberg) generator L^dagger; kernel = conserved charges."""
    Ld = 1j * (np.kron(I, H) - np.kron(H.T, I))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        Ld += (np.kron(Lk.T, Lk.conj().T)
               - 0.5 * np.kron(I, LdL) - 0.5 * np.kron(LdL.T, I))
    return Ld


def kernel_dim(M, tol=1e-9):
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s < tol * max(s.max(), 1.0)))


# ===========================================================================
# 1. documented R4 (immigration-death) -- conserved charges
# ===========================================================================
Gamma = 1.0
for x in (0.0, 2.0):                               # x=0 is homogeneous recombination
    g_in, g_out = Gamma * x, Gamma
    jumps = [np.sqrt(g_in) * ad, np.sqrt(g_out) * a] if g_in > 0 else [np.sqrt(g_out) * a]
    Ld = heisenberg_super(0 * Nop, jumps)
    kd = kernel_dim(Ld)
    print(f"R4 immigration-death, x={x}: #conserved charges = {kd}")
    check(f"R4 (x={x}) has NO nontrivial conserved charge (only identity, dim 1)",
          kd == 1)

# steady state at x=0 (homogeneous recombination): vacuum -> <N>=0, no dust
Ls = lindblad_super(0 * Nop, [a])                  # death-only (x=0)
s = np.linalg.svd(Ls, compute_uv=False)
# steady state = right null vector of Ls
_, _, Vh = np.linalg.svd(Ls)
rho_ss = Vh[-1].reshape(D, D)
rho_ss = rho_ss / np.trace(rho_ss)
Nbar_ss = float(np.real(np.trace(rho_ss @ Nop)))
print(f"R4 steady state at x=0: <N> = {Nbar_ss:.3e}  (relaxes to vacuum -> no dust)")
check("R4 at homogeneous recombination relaxes to <N>=0 (no surviving density)",
      abs(Nbar_ss) < 1e-6)

# ===========================================================================
# 2. contrast: a NON-dissipative mode conserves a number charge (dust possible)
# ===========================================================================
Ld_free = heisenberg_super(1.0 * Nop, [])          # H=eps N, NO jumps
kd_free = kernel_dim(Ld_free)
# is N itself conserved?  L^dagger(N) = 0 ?
Nvec = Nop.reshape(-1)
resid = np.linalg.norm(Ld_free @ Nvec)
print(f"\nnon-dissipative mode (no death channel): #conserved charges = {kd_free} "
      f"(all diagonal); ||L^dag(N)|| = {resid:.2e}")
check("a number charge IS conserved ONLY when the death channel is removed",
      kd_free >= D and resid < 1e-9)

# the decisive difference: add ANY death channel back -> N stops being conserved
Ld_leak = heisenberg_super(1.0 * Nop, [0.3 * a])
resid_leak = np.linalg.norm(Ld_leak @ Nvec)
print(f"add back a death/erasure channel: ||L^dag(N)|| = {resid_leak:.3e} "
      f"(N no longer conserved)")
check("any erasure (Landauer-emit) channel destroys the number charge",
      resid_leak > 1e-6)

# ===========================================================================
# verdict
# ===========================================================================
print("\n--- VERDICT ---")
print("The DOCUMENTED R4 dynamics (immigration-death, dissipative) has exactly ONE")
print("conserved charge: the identity. It carries NO conserved pressureless number")
print("charge, and at homogeneous recombination (x=0) it relaxes to the vacuum")
print("(<N>=0) -- so it supplies NO CMB dust. (Lepton number, conserved elsewhere,")
print("gives only the nu_R relic omega=0.024.)")
print()
print("A conserved number charge -> dust exists ONLY if the death/erasure channel")
print("is REMOVED (records made permanent). But R4 erasure is Landauer-DISSIPATIVE")
print("(it emits exhaust); removing the death channel is an IMPORT, not a property")
print("of W=S.C. So the omega=0.096 'zero-mode reservoir' conserved charge is")
print("POSITED, not derived: the reservoir-lift assumes the erased records are kept")
print("as a permanent, local, massive dust species rather than emitted.")
print()
print("NET: the M15 conserved-charge derivation is NOT closed. The precise missing")
print("ingredient is sharp -- a non-dissipative conservation law (permanent local")
print("records that gravitate as a^-3 dust), which CONTRADICTS the dissipative,")
print("radiation-emitting character of R4 repair. The dust-Hamiltonian form and the")
print("alpha0/208 abundance are both downstream of this unmet assumption. The CMB")
print("third peak therefore remains a genuine, live falsification risk.")

print("\n" + ("ALL CHECKS PASSED" if _ok else "SOME CHECKS FAILED"))
sys.exit(0 if _ok else 1)
