import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------
# PANEL A: Analytical Band Structure (High-Symmetry Path)
# ---------------------------------------------------------
def H_exact(kx, ky):
    """Exact Tight-Binding Bloch Hamiltonian for 4.8.8 Unit Cell"""
    H = np.zeros((4,4), dtype=complex)
    H[0,1] = H[1,0] = 1; H[1,2] = H[2,1] = 1
    H[2,3] = H[3,2] = 1; H[3,0] = H[0,3] = 1
    H[1,3] = np.exp(-1j*kx); H[3,1] = np.exp(1j*kx)
    H[0,2] = np.exp(1j*ky);  H[2,0] = np.exp(-1j*ky)
    return np.sort(np.linalg.eigvalsh(H))

# BZ Path: Gamma(0,0) -> X(pi,0) -> M(pi,pi) -> Gamma(0,0)
pts = 100
path_GX = [(k, 0) for k in np.linspace(0, np.pi, pts)]
path_XM = [(np.pi, k) for k in np.linspace(0, np.pi, pts)]
path_MG = [(k, k) for k in np.linspace(np.pi, 0, pts)]

# Stitch together without duplicating the high-symmetry points
path = path_GX + path_XM[1:] + path_MG[1:]
bands = np.zeros((4, len(path)))

for i, (kx, ky) in enumerate(path):
    bands[:, i] = H_exact(kx, ky)

fig, ax1 = plt.subplots(figsize=(8, 6))
x_ticks = [0, pts-1, 2*pts-2, len(path)-1]
labels = ['$\Gamma\ (0,0)$', '$X\ (\pi,0)$', '$M\ (\pi,\pi)$', '$\Gamma\ (0,0)$']

# Plot the 4 bands (Notice the perfectly flat stationary segments!)
ax1.plot(bands[0, :], lw=2.5, color='#1f77b4', label='Dispersive (Slow)')
ax1.plot(bands[1, :], lw=2.5, color='#ff7f0e', label='Partially Flat (Stationary $\Gamma \\rightarrow X$)')
ax1.plot(bands[2, :], lw=2.5, color='#2ca02c', label='Partially Flat (Stationary $X \\rightarrow M$)')
ax1.plot(bands[3, :], lw=2.5, color='#d62728', label='Dispersive (Fast)')

ax1.set_title("Walk Spectrum of the 4.8.8 Topology\n(High-Symmetry Path $\Gamma \\rightarrow X \\rightarrow M \\rightarrow \Gamma$)", fontsize=13)
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(labels, fontsize=12)
ax1.set_ylabel("Energy $E(\mathbf{k})$", fontsize=12)
ax1.grid(True, alpha=0.3)
for xc in x_ticks:
    ax1.axvline(x=xc, color='k', linestyle='--', alpha=0.3)
ax1.legend(loc='lower left', fontsize=10)
plt.tight_layout()
plt.show()