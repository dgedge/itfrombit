"""
TCH lattice hydrogen: convergence with lattice spacing a.

Re-run the calculation at several values of a to see how:
  - the 1s binding approaches -13.6057 eV
  - the 2s-2p "Lamb-like" splitting decreases
as a → 0.

We hold N·a (box size) approximately constant so that the same physical
volume is sampled. Larger N is more expensive, so we keep a modest box.
"""

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import time

HARTREE_eV = 27.211386245988
WATSON     = 0.252731009858

def tch_green_function(Nk):
    k = np.fft.fftfreq(Nk) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    K = 6 - 2*(np.cos(kx) + np.cos(ky) + np.cos(kz))
    inv_K = np.zeros_like(K)
    nz = K > 1e-9
    inv_K[nz] = 1.0/K[nz]
    G = np.fft.ifftn(inv_K).real
    return G + (WATSON - G[0,0,0])

def solve_tch_hydrogen(a, N, Nk=256, n_states=6):
    """Solve and return the lowest n_states eigenvalues in Hartree."""
    L = N//2
    G = tch_green_function(Nk)
    xs = np.arange(N) - L
    IX, IY, IZ = np.meshgrid(xs, xs, xs, indexing='ij')
    alpha = 4*np.pi/a
    V = -alpha * G[IX % Nk, IY % Nk, IZ % Nk]
    V_flat = V.ravel()

    t = 1.0/(2*a**2)
    I_N  = eye(N, format='csr')
    T_1d = diags([np.ones(N-1), np.ones(N-1)], [-1,1], format='csr')
    Tx = kron(kron(T_1d, I_N), I_N, format='csr')
    Ty = kron(kron(I_N, T_1d), I_N, format='csr')
    Tz = kron(kron(I_N, I_N), T_1d, format='csr')
    H = -t*(Tx+Ty+Tz) + diags([6*t + V_flat],[0], format='csr')

    evals, _ = eigsh(H, k=n_states, which='SA', tol=1e-8)
    return np.sort(evals)

# Sweep lattice spacing while keeping box size ~ constant
# (smaller a → larger N, but we cap N by computational cost)
configs = [
    (0.60, 41),   # box 24 a0
    (0.50, 51),   # box 25 a0
    (0.40, 61),   # box 24 a0
    (0.30, 71),   # box 21 a0
    (0.25, 81),   # box 20 a0
]

results = []
print(f"{'a':>5} {'N':>4} {'box':>5} | {'E_1s':>8}  {'E_2s':>8}  {'E_2p':>8} | "
      f"{'1s shift':>9}  {'2s-2p':>8}")
print(f"{'(a₀)':>5} {'':>4} {'(a₀)':>5} | {'(eV)':>8}  {'(eV)':>8}  {'(eV)':>8} | "
      f"{'(eV)':>9}  {'(eV)':>8}")
print("-"*78)

for a, N in configs:
    print(f"{a:>5.2f} {N:>4d} {(N-1)*a:>5.1f} ", end="", flush=True)
    t0 = time.time()
    ev = solve_tch_hydrogen(a, N, n_states=6)
    elapsed = time.time() - t0
    E1s, E2s = ev[0], ev[1]
    E2p = np.mean(ev[2:5])  # triplet mean
    E1s_eV = E1s * HARTREE_eV
    E2s_eV = E2s * HARTREE_eV
    E2p_eV = E2p * HARTREE_eV
    shift_1s = E1s_eV - (-13.6057)
    split_2s2p = E2s_eV - E2p_eV
    print(f"| {E1s_eV:>+8.3f}  {E2s_eV:>+8.3f}  {E2p_eV:>+8.3f} | "
          f"{shift_1s:>+9.4f}  {split_2s2p:>+8.4f}    [{elapsed:.0f}s]")
    results.append((a, N, E1s_eV, E2s_eV, E2p_eV))

results = np.array(results)
a_arr   = results[:, 0]
E1s_arr = results[:, 2]
E2s_arr = results[:, 3]
E2p_arr = results[:, 4]
shift_1s_arr = E1s_arr - (-13.6057)
split_arr    = E2s_arr - E2p_arr

# Plot convergence
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

ax = axes[0]
ax.plot(a_arr**2, shift_1s_arr, 'o-', color='C0', label='E_1s − (−13.61) eV')
ax.plot(a_arr**2, split_arr,   's-', color='C1', label='E_2s − E_2p')
ax.axhline(0, color='k', lw=0.6, ls=':')
ax.set_xlabel("a² (Bohr²)"); ax.set_ylabel("shift (eV)")
ax.set_title("TCH lattice corrections vs spacing a²")
ax.legend(); ax.grid(True, alpha=0.3)

ax = axes[1]
# log-log to see scaling exponent
ax.loglog(a_arr, -shift_1s_arr, 'o-', color='C0', label='|E_1s − ideal|')
ax.loglog(a_arr, -split_arr,    's-', color='C1', label='|E_2s − E_2p|')
# reference power laws
a_ref = np.array([0.15, 0.6])
ax.loglog(a_ref, 12*a_ref**2, 'k:', alpha=0.6, label='∝ a²')
ax.loglog(a_ref, 1.7*a_ref**3, 'k--', alpha=0.6, label='∝ a³')
ax.set_xlabel("a (Bohr)"); ax.set_ylabel("|shift| (eV)")
ax.set_title("Power-law scaling")
ax.legend(); ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig('/home/claude/tch_convergence.png', dpi=130, bbox_inches='tight')
print("\nSaved convergence plot to tch_convergence.png")
