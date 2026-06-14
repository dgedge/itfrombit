"""
TCH lattice hydrogen.

Geometry (per the TCH Gauge Web document):
  - matter sector  = octahedra centres (simple cubic)
  - gauge sector   = truncated-cube centres (simple cubic, offset by (½,½,½)·a)
  - matter-matter separations are integer multiples of a
  - Coulomb interaction between matter sites is mediated by the gauge web,
    whose static propagator is the SC lattice Green's function
        G_TCH(ρ) = (1/(2π)³) ∫_BZ d³k  e^{ik·ρ}  /  K(k)
    with K(k) = 6 − 2(cos kx + cos ky + cos kz).

We solve the single-electron Schrödinger problem on the matter sublattice
with a proton fixed at the centre and Coulomb potential V(ρ) = −α·G_TCH(ρ).
Coupling α is fixed by matching the long-range tail to atomic-unit −1/r:
    long range:   G_TCH(ρ) → 1/(4π|ρ|)
    continuum:    V(r) = −1/r,  r = ρ·a
    ⇒ α = 4π / a.

Output is in Hartree; multiply by 27.2114 for eV.
Reference Rydberg series: E_n = −13.6057 / n²  eV.
"""

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt
import time

HARTREE_eV = 27.211386245988

# ====================== TCH GAUGE-WEB GREEN'S FUNCTION ======================
def tch_green_function(Nk):
    """SC lattice Green's function on a periodic Nk³ grid, computed via FFT.

    Returns G with G[i,j,k] = G_TCH(i,j,k) for integer (i,j,k), with the
    k=0 zero-mode discarded (this only shifts G by an unimportant constant).
    """
    k = np.fft.fftfreq(Nk) * 2 * np.pi
    kx, ky, kz = np.meshgrid(k, k, k, indexing='ij')
    K = 6 - 2*(np.cos(kx) + np.cos(ky) + np.cos(kz))
    inv_K = np.zeros_like(K)
    nz = K > 1e-9
    inv_K[nz] = 1.0 / K[nz]
    return np.fft.ifftn(inv_K).real

Nk = 256
print(f"Computing TCH Green's function on {Nk}³ grid (FFT)...")
t0 = time.time()
G = tch_green_function(Nk)
print(f"  done in {time.time()-t0:.1f} s")

# Sanity check against known values
print("\nGreen's-function check (FFT result vs analytic large-ρ asymptote):")
print(f"  G(0,0,0) = {G[0,0,0]:.6f}    (Watson integral 0.252731 — ours is offset by")
print(f"                                 the discarded k=0 mode; constant shift only)")
for d in [1, 2, 3, 5, 10]:
    asymptote = 1/(4*np.pi*d)
    print(f"  G({d},0,0) = {G[d,0,0]:.6f}   asymptote 1/(4π·{d}) = {asymptote:.6f}  "
          f"diff = {G[d,0,0] - asymptote:+.5f}")

# Pin G(0) to the Watson value by shifting the constant
# (this absorbs the discarded zero-mode and makes physical predictions absolute)
WATSON = 0.252731009858
shift = WATSON - G[0, 0, 0]
G += shift
print(f"\nApplied uniform shift {shift:+.6f} to set G(0)=Watson value.")
print(f"  After shift: G(1,0,0) = {G[1,0,0]:.6f}, asymptote {1/(4*np.pi):.6f}")

# ====================== MATTER LATTICE & POTENTIAL ======================
N = 51                      # matter sites per axis
a = 0.5                     # matter-lattice spacing in Bohr radii
L = N // 2                  # central index (proton sits here)

xs = np.arange(N) - L
IX, IY, IZ = np.meshgrid(xs, xs, xs, indexing='ij')
rho = np.sqrt(IX**2 + IY**2 + IZ**2)            # |ρ| in lattice units

# V_i = -α · G_TCH(displacement from proton site at origin)
alpha = 4*np.pi / a
V = -alpha * G[IX % Nk, IY % Nk, IZ % Nk]
V_flat = V.ravel()

print(f"\nMatter lattice: N={N}, a={a} a₀, box across = {(N-1)*a:.1f} a₀")
print(f"Coulomb coupling α = 4π/a = {alpha:.4f} (atomic units)")
print(f"V at proton site (ρ=0): {V_flat[L*N*N + L*N + L]:.4f} Ha")
print(f"V at first neighbour (ρ=1): {V[L+1, L, L]:.4f} Ha   (continuum −1/a = {-1/a:.4f})")
print(f"V at ρ=10: {V[L+10, L, L]:.6f} Ha   (continuum = {-1/(10*a):.6f})")

# ====================== HAMILTONIAN ======================
# H = (1/(2a²)) (6 I − Σ neighbour-shifts) + diag(V)
#   = -t (Tx + Ty + Tz) + diag(6t + V)
t = 1.0 / (2 * a**2)

I_N  = eye(N, format='csr')
T_1d = diags([np.ones(N-1), np.ones(N-1)], [-1, 1], format='csr')
Tx = kron(kron(T_1d, I_N), I_N, format='csr')
Ty = kron(kron(I_N, T_1d), I_N, format='csr')
Tz = kron(kron(I_N, I_N), T_1d, format='csr')

H = -t*(Tx + Ty + Tz) + diags([6*t + V_flat], [0], format='csr')
print(f"\nHamiltonian: {N**3} × {N**3} sparse matrix, t = {t:.3f} Ha")

# ====================== DIAGONALISE ======================
n_states = 16
print(f"\nDiagonalising for lowest {n_states} eigenstates...")
t0 = time.time()
evals, evecs = eigsh(H, k=n_states, which='SA', tol=1e-7)
order = np.argsort(evals); evals = evals[order]; evecs = evecs[:, order]
print(f"  done in {time.time()-t0:.1f} s")

print(f"\n  state    E (Ha)      E (eV)     guessed n     hydrogen exact (eV)")
print("  " + "-"*68)
for i, E in enumerate(evals):
    Eev = E * HARTREE_eV
    if E < -5e-4:
        n_guess = max(1, int(round(1.0/np.sqrt(-2*E))))
        ideal = -HARTREE_eV/(2*n_guess**2)
        ideal_str = f"n={n_guess}: {ideal:+.3f}"
    else:
        ideal_str = "(unbound or near continuum)"
    print(f"  {i+1:>5d}   {E:>9.5f}   {Eev:>+8.3f}     n ≈ {n_guess}      {ideal_str}")

# ====================== ANALYSE DEGENERACIES & SPLITTINGS ======================
print("\nShell groupings (states grouped by E within 0.5 mHa = 0.014 eV):")
groups, current = [], [0]
for i in range(1, len(evals)):
    if abs(evals[i] - evals[current[-1]]) < 5e-4:
        current.append(i)
    else:
        groups.append(current); current = [i]
groups.append(current)
shell_n = 1
for g in groups:
    Emean = np.mean(evals[g]) * HARTREE_eV
    n_guess = max(1, int(round(1.0/np.sqrt(-2*np.mean(evals[g]))))) if np.mean(evals[g]) < -5e-4 else None
    deg = len(g)
    label = (f"degeneracy {deg}, ⟨E⟩ = {Emean:+.4f} eV"
             + (f"  (≈ n={n_guess}, ideal {-HARTREE_eV/(2*n_guess**2):+.3f} eV)" if n_guess else ""))
    print(f"  states {g[0]+1}–{g[-1]+1}: {label}")

# ====================== PLOTS ======================
psi  = evecs.T.reshape(n_states, N, N, N)
dens = np.abs(psi)**2

fig = plt.figure(figsize=(13, 10))

# --- spectrum (lattice vs continuum) ---
ax = fig.add_subplot(2, 3, 1)
for E in evals:
    ax.hlines(E*HARTREE_eV, 0, 1, colors='C0', linewidth=2)
for n in (1, 2, 3, 4):
    Eid = -HARTREE_eV/(2*n**2)
    ax.hlines(Eid, 1.4, 2.4, colors='crimson', linestyles='--', linewidth=1.3)
    ax.text(2.5, Eid, f"  n={n}: {Eid:+.2f}", va='center', fontsize=9)
ax.set_xticks([0.5, 1.9]); ax.set_xticklabels(['TCH lattice', 'ideal H'])
ax.set_xlim(0, 4.6); ax.set_ylabel("E (eV)"); ax.set_title("Spectrum")
ax.grid(True, alpha=0.3)

# --- 1s density (xy slice through proton) ---
ax = fig.add_subplot(2, 3, 2)
im = ax.pcolormesh(xs*a, xs*a, dens[0][:, :, L], shading='auto', cmap='viridis')
ax.set_aspect('equal')
ax.set_title(f"1s   E = {evals[0]*HARTREE_eV:.3f} eV")
ax.set_xlabel("x (a₀)"); ax.set_ylabel("y (a₀)")
plt.colorbar(im, ax=ax)

# --- 2p density (one of the three, xy slice) ---
ax = fig.add_subplot(2, 3, 3)
vmax = np.max(np.abs(psi[1]))
im = ax.pcolormesh(xs*a, xs*a, psi[1][:, :, L], shading='auto',
                   cmap='RdBu_r', vmin=-vmax, vmax=vmax)
ax.set_aspect('equal')
ax.set_title(f"2p   E = {evals[1]*HARTREE_eV:.3f} eV")
ax.set_xlabel("x (a₀)"); ax.set_ylabel("y (a₀)")
plt.colorbar(im, ax=ax)

# Find candidate 3d states (states whose ⟨E⟩ is near n=3 = -1.51 eV)
n3_indices = [i for i, E in enumerate(evals) if -0.07 < E < -0.045]
if len(n3_indices) >= 5:
    # 3d shell expected at indices (after 3s, 3p): try the last 5 of n=3 group
    n3 = n3_indices[-5:]
    print(f"\nCandidate 3d-like states (indices, energies in eV):")
    for i in n3:
        print(f"  state {i+1}: {evals[i]*HARTREE_eV:.4f} eV")

    # --- 3d cubic-field splitting plot ---
    ax = fig.add_subplot(2, 3, 4)
    for i in n3:
        ax.hlines(evals[i]*HARTREE_eV, 0, 1, colors='C0', linewidth=2)
    Eid_3 = -HARTREE_eV/(2*9)
    ax.hlines(Eid_3, 1.4, 2.4, colors='crimson', linestyles='--', linewidth=1.3)
    ax.text(2.5, Eid_3, f"  ideal n=3", va='center', fontsize=9)
    ax.set_xticks([0.5, 1.9]); ax.set_xticklabels(['TCH 3d', 'ideal'])
    ax.set_xlim(0, 4.6); ax.set_ylabel("E (eV)")
    ax.set_title(f"3d shell — Oₕ splitting\nΔ = {(evals[n3[-1]]-evals[n3[0]])*HARTREE_eV*1000:.2f} meV")
    ax.grid(True, alpha=0.3)

    # --- 3d eg-like density (one with lobes along axes) ---
    ax = fig.add_subplot(2, 3, 5)
    state = n3[-1]   # likely eg
    vmax = np.max(np.abs(psi[state]))
    im = ax.pcolormesh(xs*a, xs*a, psi[state][:, :, L], shading='auto',
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_aspect('equal')
    ax.set_title(f"state {state+1} (3d-like)\nE = {evals[state]*HARTREE_eV:.4f} eV")
    ax.set_xlabel("x (a₀)"); ax.set_ylabel("y (a₀)")
    plt.colorbar(im, ax=ax)

    # --- 3d t2g-like density (one with lobes between axes) ---
    ax = fig.add_subplot(2, 3, 6)
    state = n3[0]   # likely t2g
    vmax = np.max(np.abs(psi[state]))
    im = ax.pcolormesh(xs*a, xs*a, psi[state][:, :, L], shading='auto',
                       cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    ax.set_aspect('equal')
    ax.set_title(f"state {state+1} (3d-like)\nE = {evals[state]*HARTREE_eV:.4f} eV")
    ax.set_xlabel("x (a₀)"); ax.set_ylabel("y (a₀)")
    plt.colorbar(im, ax=ax)

plt.suptitle("TCH lattice hydrogen — gauge-web Coulomb potential", fontsize=13)
plt.tight_layout()
plt.savefig('/home/claude/tch_hydrogen.png', dpi=130, bbox_inches='tight')
print("\nSaved plot to tch_hydrogen.png")
